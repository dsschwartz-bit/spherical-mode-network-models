%% COMPARE_OOK_CST_ADS_SIMULINK
% OOK direct antenna modulation (DAM, rectangular patch): radiated far-field time signals
% r*E_theta(t), r*E_phi(t) from CST, ADS and the Simulink RF Blockset Tx-only model with the
% switch at the internal port (as run_ook_simulink_rf.m), compared per observation direction
% (time signal, difference, spectrum) and over radial (spoke) cuts.
%
% Only results shared by ALL the requested sources (show_sources) are compared and plotted:
%   directions : those present in every requested data set (CST probes of ff_files, ADS grid of
%                ads_file); Simulink is simulated at exactly these directions (one model build,
%                Fast Restart; the direction only changes gains: link_update_gains). With Simulink
%                and only one of CST / ADS, all that source's directions are used.
%   components : zero fields (peak < 1 mV, e.g. E_phi in the phi = 0 / 180 planes) are skipped.
%   On the z axis (theta = 0) all phi are the same point: the ADS components are rotated to
%   the CST phi there.
% Simulink: RF Blockset model with the switch (as run_ook_simulink_rf.m), the S-parameter block's own
% fit, driven by the CST port V/I exports or the user excitation (excitation_source). Results are
% kept in <out_dir>/simulink_out/<...>.mat and reused while the inputs (files, dates, settings) are
% unchanged; simulink_mode = "rerun" simulates again.
%
% Errors over the overlap of the records, relative to the reference (the first of CST, ADS,
% Simulink requested): relative MSE = integral (x - ref)^2 dt / integral ref^2 dt, the mean |dB|
% spectrum magnitude error over the band where the reference is within bw_level_dB of its peak,
% and the total rel. MSE over all compared signals = sum of error energies / sum of reference energies.
%
% Outputs (in out_dir): plots/ (PNG) and simulink_out/ (Simulink result cache); the RF model and
% its _vars.mat in model_dir.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(fullfile(here, '..'));                      % plot_rE_radial_simulink

%% ===== USER SETTINGS =====
cst_root  = "C:\path\to\CST\";                                      % folder containing the CST projects
ads_file  = "C:\path\to\ADS\<workspace>_wrk\rectangularpatch_rE_filteredsine_2.45GHz.csv";   % ADS Data Display CSV export
snp_dir   = fullfile(repo_root, 'output', 'SnP');                       % Touchstone files
out_dir   = fullfile(repo_root, 'output', 'simulink', 'ook');           % plots, Simulink cache
model_dir = out_dir;                                                    % RF model + <model>_vars.mat
% =========================

set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
try, set(groot, 'defaultFigureWindowStyle', 'docked'); catch, end

%% ---- what to present -------------------------------------------------------------------
show_sources  = ["CST" "ADS" "Simulink"];   % any of "CST", "ADS", "Simulink"
simulink_mode = "auto";      % "auto": saved results, simulate missing directions; "rerun"; "cache"
rect_dirs     = "all";       % per-direction figures: "all" (every shared direction) or [theta phi; ...] (shared ones only)
excitation_source = "user";   % Simulink input: "cst" = the CST external-port V/I exports below; "user" = user_exc_file
fast_restart  = true;        % Simulink Fast Restart over the directions (compile once)
radial_cuts   = struct('cut',     {"theta", "phi"}, ...
                       'angle',   {[0],       [90]}, ...
                       'spokes',  {[],      []});
% one cut per column (commas, not semicolons: every {...} must be 1 x N with the same N):
%   'cut'    : the angle that varies along the cut ("theta": spokes in theta; "phi": spokes in phi)
%   'angle'  : the angle held fixed: phi (deg) of a "theta" cut, theta (deg) of a "phi" cut;
%              a vector gives one cut per value, e.g. {[0 90], 45} = theta cuts at phi 0 and 90, phi cut at theta 45
%   'spokes' : [] = every shared direction on the cut, or a list of spoke angles to keep
rc = radial_cuts([]);                                        % a vector 'angle' -> one cut per value
for q = radial_cuts(:).', for a = q.angle(:).', q1 = q; q1.angle = a; rc(end+1) = q1; end, end %#ok<AGROW>
radial_cuts = rc;

%% ---- files and settings ---------------------------------------------------------------------
antenna  = "rectangularpatch_symmetry_DAM";
sig_dir  = fullfile(cst_root, antenna, "Export", "time_signals");
voltage_file = "excitation_filteredsine_voltage_2.45GHz.txt";
current_file = "excitation_filteredsine_current_2.45GHz.txt";
ff_files     = fullfile(sig_dir, ["ff_Etheta_filteredsine_port1_2.45GHz.txt" "ff_Ephi_filteredsine_port1_2.45GHz.txt"]);
ads_phi  = ads_phi_default();     % deg, phi of the ADS sweep, used only if the CSV has no phi column
user_exc_file = fullfile(repo_root, 'output', 'cst_ads', 'user_excitations', ...
    'excitation_sine2.45GHz_30ns_step66.7MHz_30ns_Vinc1_50Ohms_Vs2.mat');   % from generate_user_excitations.m
snp      = fullfile(snp_dir, 'S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p');
mdl      = 'Txrectangularpatch_internal_ports_OOK_RF';
rf_port = 1;  ctrl_port = 2;
step_s     = 2e-12;               % RF Blockset max step
rf_fitopt  = 'Share all poles';
% Simulink solver (RF platform), user choice: "ode45" = non-stiff (linear network) | "ode15s" = stiff (switches, time-varying R/L/C).
% The RF circuit itself is stepped at MaxStep (default Nyquist/6) whichever solver is chosen.
rf_solver  = 'ode15s';            % ode15s here: the OOK model has a switch at the internal port
max_poles  = 99;                  % RF block limit
fit_tol_dB = -40;
components = ["Theta" "Phi"];
angle_tol  = 0.1;                 % deg
bw_level_dB = 20;
f_plot_GHz  = [0 5];
save_plots  = true;
plot_dir    = fullfile(out_dir, 'plots');
sim_cache   = fullfile(out_dir, 'simulink_out', "ook_rf_" + antenna + "_filteredsine_2.45GHz.mat");

show_sources = string(show_sources);
order = ["CST" "ADS" "Simulink"];  shown = order(ismember(order, show_sources));
useS = any(shown == "Simulink");  useA = any(shown == "ADS");  useC = any(shown == "CST");
sty = struct('CST',      struct('col', [0 0 0],       'ls', '-', 'lw', 1.2), ...
             'ADS',      struct('col', [0.85 0 0],    'ls', ':', 'lw', 1.6), ...
             'Simulink', struct('col', [0.55 0.75 1], 'ls', '-', 'lw', 3.5));

if save_plots && isfolder(plot_dir)       % the previous run's comparison plots (other direction sets) are replaced
    delete(fullfile(plot_dir, 'ook_compare_E*.png'));  delete(fullfile(plot_dir, 'ook_radial_*'));
end

%% ---- excitation: CST exports (as ook_cst_to_ads_mdf.m) or the user excitation ------------------
% the Simulink input
switch lower(string(excitation_source))       % only the selected source is read
    case "cst",  E = ook_load_excitation("cst", struct('sig_dir', sig_dir, 'voltage_file', voltage_file, ...
                     'current_file', current_file, 'rf_port', rf_port, 'ctrl_port', ctrl_port));  exc_file = fullfile(sig_dir, voltage_file);
    case "user", E = ook_load_excitation("user", struct('user_file', user_exc_file));  exc_file = user_exc_file;
    otherwise,   error('excitation_source must be "cst" or "user".');
end
t = E.t;  Vinc = E.Vinc;  Vs = E.Vs;  Z = E.Z;  t0 = E.t0;

%% ---- CST and ADS ----------------------------------------------------------------------------
cst = readCstFarfieldProbes(ff_files, components, 1, 1e-3);
r_probe = cst(1).r;  a_tx = readCstA(ff_files(1), 1e-3);
ads = [];  if useA, ads = readAdsFarfieldSweeps(ads_file, components, ads_phi); end
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
same = @(th1, ph1, th2, ph2) abs(th1 - th2) < angle_tol && (th1 < angle_tol || abs(th1 - 180) < angle_tol || angd(ph1, ph2) < angle_tol);
exact  = @(S, th, ph) find(arrayfun(@(s) abs(s.theta - th) < angle_tol && angd(s.phi, ph) < angle_tol, S), 1);
findIn = @(S, th, ph) [exact(S, th, ph), find(arrayfun(@(s) same(s.theta, s.phi, th, ph), S), 1)];   % exact phi first, else z axis at any phi
findIn = @(S, th, ph) firstOf(findIn(S, th, ph));

% shared directions: present in every requested data set (CST / ADS); Simulink can be simulated anywhere.
% In the CST frame where CST is requested (z axis: ADS matched at any phi and rotated), else in the ADS frame.
cstDirs = [[cst.theta].' [cst.phi].'];
if useA, adsDirs = [[ads.theta].' [ads.phi].']; else, adsDirs = zeros(0, 2); end
if useC, allDirs = cstDirs; elseif useA, allDirs = adsDirs; else, allDirs = cstDirs; end   % Simulink only: the CST probes
if useC && useA
    allDirs = allDirs(arrayfun(@(d) ~isempty(findIn(ads, allDirs(d, 1), allDirs(d, 2))), 1:size(allDirs, 1)), :);
end
if isnumeric(rect_dirs)
    rdirs = rect_dirs(arrayfun(@(d) any(abs(allDirs(:, 1) - rect_dirs(d, 1)) < angle_tol & angd(allDirs(:, 2), rect_dirs(d, 2)) < angle_tol), ...
        1:size(rect_dirs, 1)), :);
else
    rdirs = allDirs;
end
fprintf('%s, OOK DAM, r = %g m: %d CST probe(s), %d ADS direction(s); %d direction(s) shared by %s\n', ...
    antenna, r_probe, numel(cst), size(adsDirs, 1), size(allDirs, 1), strjoin(shown, ', '));
if isempty(allDirs), error('No direction is shared by %s.', strjoin(shown, ', ')); end

%% ---- Simulink: one build, every direction --------------------------------------------------
simAll = struct('theta', {}, 'phi', {}, 't', {}, 'e', {});
if useS
    fd = @(f) getfield(dir(f), 'datenum'); %#ok<GFLD>
    key = struct('snp', string(snp), 'snp_date', fd(snp), 'exc', string(excitation_source), 'exc_file', string(exc_file), ...
        'exc_date', fd(exc_file), 'r', r_probe, 'step', step_s, 'fit', fit_tol_dB, 'cst_stop', useC, 'ff_date', fd(ff_files(1)), ...
        'fitopt', string(rf_fitopt), 'solver', string(rf_solver), 'poles', max_poles, 'builder', fd(which('build_tx_rx_link_simulink')));
    if isfile(sim_cache) && simulink_mode ~= "rerun"
        S = load(sim_cache, 'key', 'cache');
        if isequaln(S.key, key), simAll = S.cache; else, fprintf('Simulink cache made with other inputs: simulating again.\n'); end
    end
    miss = arrayfun(@(d) ~any(arrayfun(@(s) abs(s.theta - allDirs(d, 1)) < angle_tol && angd(s.phi, allDirs(d, 2)) < angle_tol, simAll)), 1:size(allDirs, 1));
    if any(miss)
        if simulink_mode == "cache", error('Simulink results missing in %s.', sim_cache); end
        assignin('base', 'Vtx_wave', timeseries(Vinc, t - t0, 'Name', 'Vtx_wave'));
        assignin('base', 'Vs_ctrl',  timeseries(Vs,   t - t0, 'Name', 'Vs_ctrl'));
        assignin('base', 'Z0', Z);
        todo = allDirs(miss, :);
        build_tx_rx_link_simulink(snp, '', 'TxOnly', true, 'Platform', 'rf', 'LoadSwitch', true, 'ControlVar', 'Vs_ctrl', ...
            'Waveform', 'workspace', 'WaveformVar', 'Vtx_wave', 'SourceType', 'generator', 'Z0', Z, 'FeedZ', Z, ...
            'r', r_probe, 'theta_deg', todo(1, 1), 'phi_deg', todo(1, 2), 'ModelName', mdl, 'MaxStep', step_s, ...
            'MaxPoles', max_poles, 'FitTolerance', fit_tol_dB, 'RFSolver', rf_solver, 'FitOptions', rf_fitopt, ...
            'OutputDir', model_dir, 'OpenModel', false);
        tStop = [];                                            % with CST shown: run to the end of the CST records
        if useC, tStop = max(arrayfun(@(c) max(cellfun(@(x) max([x(:); -Inf]), c.t)), cst)) - t0; end
        rs = ook_observation_loop(mdl, todo, r_probe, [], t0, 'RF Blockset', fast_restart, false, tStop);   % simulate only (no figures)
        for d = 1:numel(rs)
            simAll(end+1) = struct('theta', rs(d).theta, 'phi', rs(d).phi, 't', rs(d).t, 'e', rs(d).e); %#ok<SAGROW>
        end
        if ~isfolder(fileparts(sim_cache)), mkdir(fileparts(sim_cache)); end
        cache = simAll;  save(sim_cache, 'key', 'cache');
    else
        fprintf('Simulink: all %d direction(s) from %s\n', size(allDirs, 1), sim_cache);
    end
end
simExact = @(th, ph) find(arrayfun(@(s) abs(s.theta - th) < angle_tol && angd(s.phi, ph) < angle_tol, simAll), 1);

%% ---- per direction ---------------------------------------------------------------------------
results = struct('theta', {}, 'phi', {}, 'comp', {}, 'src', {}, 'peak', {}, 'relmse', {}, 'magerr', {}, 'eerr', {}, 'eref', {});
for d = 1:size(rdirs, 1)
    th = rdirs(d, 1);  ph = rdirs(d, 2);
    for c = 1:numel(components)
        names = strings(0);  T = {};  E = {};
        if useC, j = findIn(cst, th, ph);
            if ~isempty(j) && angd(cst(j).phi, ph) < angle_tol, names(end+1) = "CST"; T{end+1} = cst(j).t{c}(:); E{end+1} = cst(j).e{c}(:); end %#ok<SAGROW>
        end
        if useA, k = findIn(ads, th, ph);
            if ~isempty(k)
                e = rotPhi(ads(k).e, ads(k).phi, ph);            % on the z axis: components at this phi
                names(end+1) = "ADS"; T{end+1} = ads(k).t(:); E{end+1} = e(:, c); %#ok<SAGROW>
            end
        end
        if useS, m = simExact(th, ph);
            if ~isempty(m), names(end+1) = "Simulink"; T{end+1} = simAll(m).t(:); E{end+1} = simAll(m).e(:, c); end %#ok<SAGROW>
        end
        nS = numel(names);  if nS < 2 || ~isequal(sort(names), sort(shown)), continue, end   % every selected source must have data
        for s = 1:nS, [T{s}, iu] = unique(T{s}, 'last'); E{s} = E{s}(iu); end
        pk = max(cellfun(@(x) max(abs(x)), E));
        if pk < 1e-3, continue, end                                  % null component (symmetry plane)
        tB = max(cellfun(@(x) x(1), T));  tE = min(cellfun(@(x) x(end), T));
        dt = min(cellfun(@(x) median(diff(x)), T)) / 4;  tg = (tB:dt:tE).';
        X = zeros(numel(tg), nS);
        for s = 1:nS, X(:, s) = interp1(T{s}, E{s}, tg, 'pchip', 0); end
        xr = X(:, 1);  rm = nan(1, nS);  eM = nan(1, nS);  eE = nan(1, nS);  eR = trapz(tg, xr.^2);
        for s = 2:nS, eE(s) = trapz(tg, (X(:, s) - xr).^2);  rm(s) = eE(s) / eR;  [eM(s), fl, fh] = bandMagError(tg, xr, X(:, s), bw_level_dB); end
        results(end+1) = struct('theta', th, 'phi', ph, 'comp', components(c), 'src', names, 'peak', max(abs(X)), 'relmse', rm, 'magerr', eM, ...
            'eerr', eE, 'eref', eR); %#ok<SAGROW>

        fig = figure('Color', 'w', 'Name', sprintf('E_%s th %g ph %g', lower(components(c)), th, ph), 'NumberTitle', 'off');
        tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
        title(tl, sprintf('$rE_{\\%s}$, $\\theta = %g^\\circ$, $\\phi = %g^\\circ$, $r = %g$ m, OOK DAM: %s', lower(components(c)), th, ph, r_probe, strjoin(names, ' vs ')), 'Interpreter', 'latex');
        dOrd = [find(names == "Simulink"), find(names ~= "Simulink")];
        nexttile(tl); hold on;  h = gobjects(1, nS);
        for s = dOrd, st = sty.(names(s)); h(s) = plot(tg*1e9, X(:, s), st.ls, 'Color', st.col, 'LineWidth', st.lw); end
        grid on; box on; xlabel('Time (ns)'); ylabel('$rE$ (V)'); legend(h, cellstr(names), 'Location', 'best'); xlim([tB tE]*1e9);
        title(sprintf('rel. MSE compared to %s: %s', names(1), strjoin(compose('%s %.2g\\%%', names(2:end).', 100*rm(2:end).'), ', ')));
        nexttile(tl); hold on
        for s = 2:nS, st = sty.(names(s)); plot(tg*1e9, X(:, s) - xr, '-', 'Color', st.col, 'LineWidth', min(st.lw, 2)); end
        grid on; box on; xlabel('Time (ns)'); ylabel(sprintf('difference to %s (V)', names(1))); xlim([tB tE]*1e9);
        legend(cellstr(names(2:end) + " $-$ " + names(1)), 'Location', 'best');  title('differences');
        nexttile(tl); hold on
        Nf = 2^nextpow2(8*numel(tg));  f = (0:Nf/2).'/(Nf*dt);  Xr = max(abs(fft(xr, Nf)));  h = gobjects(1, nS);
        for s = dOrd, Xs = abs(fft(X(:, s), Nf)); Xs = Xs(1:Nf/2+1); st = sty.(names(s));
            h(s) = plot(f/1e9, 20*log10(max(Xs, eps)/Xr), st.ls, 'Color', st.col, 'LineWidth', st.lw); end
        xlim(f_plot_GHz); ylim([-80 5]); grid on; box on; xlabel('Frequency (GHz)'); ylabel('dB normalized');
        legend(h, cellstr(names), 'Location', 'south');
        title(sprintf('mean $|X|$ error (%.3g--%.3g GHz): %s', fl/1e9, fh/1e9, strjoin(compose('%s %.2f dB', names(2:end).', eM(2:end).'), ', ')));
        if save_plots
            if ~isfolder(plot_dir), mkdir(plot_dir); end
            stem = sprintf('ook_compare_E%s_th%g_ph%g', lower(components(c)), th, ph);
            exportgraphics(fig, fullfile(plot_dir, stem + ".png"), 'Resolution', 120);
        end
    end
end

%% ---- radial (spoke) cuts ---------------------------------------------------------------------
for q = 1:numel(radial_cuts)
    src = shown;  ang = radial_cuts(q).angle;  isTh = radial_cuts(q).cut == "theta";
    % directions of the cut (exact phi; the z axis only at phi = ang, so no frame rotation is needed)
    inCut = @(th, ph) (isTh && (angd(ph, ang) < angle_tol || (th > angle_tol && angd(ph, ang + 180) < angle_tol))) || ...
                      (~isTh && abs(th - ang) < angle_tol);
    has = @(S, th, ph) any(arrayfun(@(s) abs(s.theta - th) < angle_tol && angd(s.phi, ph) < angle_tol, S));
    cand = allDirs(arrayfun(@(d) inCut(allDirs(d, 1), allDirs(d, 2)), 1:size(allDirs, 1)), :);
    keep = true(size(cand, 1), 1);
    for d = 1:size(cand, 1)
        if useC, keep(d) = keep(d) && has(cst, cand(d, 1), cand(d, 2)); end
        if useA, keep(d) = keep(d) && has(ads, cand(d, 1), cand(d, 2)); end
        if useS, keep(d) = keep(d) && has(simAll, cand(d, 1), cand(d, 2)); end
    end
    cand = cand(keep, :);
    if isTh, sp = cand(:, 1) .* (1 - 2*(angd(cand(:, 2), ang) >= angle_tol)); else, sp = cand(:, 2); end   % signed theta / phi
    if ~isempty(radial_cuts(q).spokes), sp = sp(arrayfun(@(x) any(abs(x - radial_cuts(q).spokes) < angle_tol), sp)); end
    if numel(sp) < 2
        fprintf('%s cut %g: %d direction(s) with data from all of %s: not plotted\n', radial_cuts(q).cut, ang, numel(sp), strjoin(src, '/'));
        continue
    end
    sf = '';  if save_plots, sf = fullfile(plot_dir, sprintf('ook_radial_%s_%s%g.png', strjoin(src, '_'), radial_cuts(q).cut, ang)); end
    inC = arrayfun(@(s) inCut(s.theta, s.phi), simAll);
    [~, res] = plot_rE_radial_simulink(simAll(inC), cst, ads, 'Sources', src, 'Cut', radial_cuts(q).cut, 'Angle', ang, ...
        'Spokes', sp.', 'r', r_probe, 'a', a_tx, 'AngleTol', angle_tol, 'SaveFile', sf, 'Verbose', false, ...
        'Title', sprintf('Rectangular patch, OOK DAM: %s', strjoin(src, ' vs ')));
    if numel(src) > 1 && ~isempty(res) && isfield(res, 'eref')
        er = vertcat(res.eref);  ed = vertcat(res.ediff);  ok = ~isnan(ed(:, end));
        ln = sprintf('%s cut %g, whole-cut rel. MSE compared to %s:', radial_cuts(q).cut, radial_cuts(q).angle, src(1));
        col = struct('CST', 1, 'ADS', 2, 'Simulink', 3);
        for s = src(2:end), j = col.(s); ok = ~isnan(ed(:, j)); if any(ok), ln = [ln sprintf(' %s %.3g %%', s, 100*sum(ed(ok, j))/sum(er(ok)))]; end, end %#ok<AGROW>
        fprintf('%s\n', ln);
    end
end

%% ---- summary ----------------------------------------------------------------------------------
fprintf('\n theta   phi  comp   peak (V)                  rel. MSE / mean |X| error compared to the first source\n');
for j = 1:numel(results)
    q = results(j);
    pk = strjoin(compose('%s %.3g', q.src.', q.peak.'), ', ');
    er = strjoin(compose('%s %.2f %% / %.2f dB', q.src(2:end).', 100*q.relmse(2:end).', q.magerr(2:end).'), ', ');
    fprintf(' %5g %5g  %-5s  %-26s %s\n', q.theta, q.phi, q.comp, pk, er);
end
% total over every compared (shared, non-zero) signal: sum of error energies / sum of reference energies
if ~isempty(results) && numel(shown) > 1
    eRefAll = [results.eref];  ln = sprintf('Total rel. MSE compared to %s over %d shared non-zero signal(s):', shown(1), numel(results));
    for s = 2:numel(shown)
        eErrAll = arrayfun(@(q) q.eerr(q.src == shown(s)), results);
        ln = [ln sprintf(' %s %.3f %%', shown(s), 100*sum(eErrAll)/sum(eRefAll))]; %#ok<AGROW>
    end
    fprintf('%s\n', ln);
end
if save_plots, fprintf('Plots in %s\n', plot_dir); end

%% ---- local functions ------------------------------------------------------------------------------
function k = firstOf(v)
if isempty(v), k = []; else, k = v(1); end
end

function p = ads_phi_default()
p = 0;            % deg: the ADS sweep has only E_theta, i.e. the E plane phi = 0/180; its sign matches the CST/Simulink E_theta at phi = 0
end

function e = rotPhi(e, phFrom, phTo)
% field components [E_theta E_phi] at the z axis expressed at another phi (same point, rotated unit vectors)
dp = (phTo - phFrom)*pi/180;
e = [cos(dp)*e(:, 1) + sin(dp)*e(:, 2), -sin(dp)*e(:, 1) + cos(dp)*e(:, 2)];
end

function a = readCstA(file, unit)
a = 0;  L = readlines(file);  L = L(startsWith(strtrim(L), "#Parameters"));
if isempty(L), return, end
tok = regexp(L(1), '[{;\s]a\s*=\s*([-+\deE.]+)', 'tokens', 'once');
if ~isempty(tok), a = str2double(tok{1}) * unit; end
end

function [eM, fLo, fHi] = bandMagError(tg, xRef, x, lev)
dt = tg(2) - tg(1);  Nf = 2^nextpow2(8*numel(tg));  f = (0:Nf/2).'/(Nf*dt);
Xr = abs(fft(xRef, Nf));  Xr = Xr(1:Nf/2 + 1);  Xx = abs(fft(x, Nf));  Xx = Xx(1:Nf/2 + 1);
[pk, ip] = max(Xr);  thr = pk*10^(-lev/20);
il = find(Xr(1:ip) < thr, 1, 'last');  if isempty(il), il = 1; end
iu = ip - 1 + find(Xr(ip:end) < thr, 1, 'first');  if isempty(iu), iu = numel(f); end
m = il:iu;  eM = mean(abs(20*log10(max(Xx(m), eps) ./ max(Xr(m), eps))));  fLo = f(il);  fHi = f(iu);
end
