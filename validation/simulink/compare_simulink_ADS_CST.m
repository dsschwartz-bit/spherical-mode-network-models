%% COMPARE_SIMULINK_ADS_CST
% Compare the radiated far-field time signals r*E_theta(t) and r*E_phi(t) of the
% Simulink Tx spherical-mode model (build_tx_rx_link_simulink.m), the ADS result and the
% CST far-field export, at observation directions (theta, phi) at distance r (the CST
% probe distance).
%
% Inputs : CST time-signal exports (cst_root), ADS Data Display CSV export (ads_dir) and the
%          Touchstone file (snp_dir); see USER SETTINGS.
% Outputs: figures (PNG + .fig) in plot_dir, the model in model_dir, struct array "results".
%
% Sources: show_sources picks which of "CST", "ADS", "Simulink" are presented (any
% combination).  Errors are relative to the reference source: the first of CST, ADS,
% Simulink that is shown.  Simulink is only built and simulated when it is shown.
%
% Antenna: antenna_case selects a block of settings below ("rectangularpatch",
% "circularloop"); add a case for another antenna.
%
% Simulink : the model is simulated with the CST port signal (iowaves export, signal i1)
%            as Tx waveform; for an antenna with an internal (load) port that port is open
%            ('OC') or shorted ('SC'); Vtheta / Vphi are r*E_theta / r*E_phi in V.
% ADS      : CSV Data Display export <antenna>[_<OC|SC>]_rE_<excitation>_<band>.csv
%            (readAdsFarfieldSweeps.m), columns Vtheta_prime / Vphi_prime.
% CST      : far-field time-signal exports ff_Etheta_* / ff_Ephi_* in
%            <CST project>\Export\time_signals\ (with a load port: ...\time_signals\port<N>_<OC|SC>\)
%            (readCstFarfieldProbes.m), r*E in V.  Always read (probe distance r and the
%            available directions), shown only if selected.
%
% Model: the radiated fields come from the TX-ONLY copy of the Simulink model
% (build_tx_rx_link_simulink(..., 'TxOnly', true), named Tx<antenna>[_internal_ports] like the ADS cell: the full model
% name without the Rx antenna).  Its delay is (r - a_Tx)/c0, exactly the probe delay of the CST/ADS
% far-field data, so no Rx de-embedding shift is needed.
%
% Time axis: every plot (per-direction figures and radial plots) ends at the minimum of the last
% times of the records shown.
%
% Figures are docked as tabs of the standard MATLAB figure window (WindowStyle 'docked').
%
% For every direction and component, one figure row shows: the time signals, the differences
% to the reference, and the spectrum magnitudes.  Printed errors are relative to the reference
% over the overlap of the records: relative MSE = integral (x - ref)^2 dt / integral ref^2 dt,
% the peaks, the best-fit delay and the in-band (-bw_level_dB) mean magnitude error.
%
% Radial (spoke) plots of the same sources over a cut: plot_rE_radial_simulink.m (called below;
% see radial_cuts).
%
% EXAMPLES (copy the lines between %{ and %} into the MATLAB console):
%{
compare_simulink_ADS_CST                       % settings below
% edit  show_sources = ["CST" "ADS" "Simulink"];  and run again for all three
results(1).rel_mse                             % relative MSE of each shown source vs the reference
%}

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(here);                                % plot_rE_radial_simulink

%% ===== USER SETTINGS =====
cst_root  = "C:\path\to\CST\";                 % folder containing the CST projects (<antenna>_<variant>\Export\time_signals)
ads_dir   = "C:\path\to\ADS\<workspace>_wrk\"; % folder with the ADS Data Display CSV exports
snp_dir   = fullfile(repo_root, 'output', 'SnP');           % Touchstone files
model_dir = fullfile(repo_root, 'output', 'simulink');      % Simulink model + <model>_vars.mat
plot_dir  = fullfile(repo_root, 'output', 'simulink', 'plots');
% =========================

set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');   % LaTeX fonts
set(groot, 'defaultFigureWindowStyle', 'docked');     % every new figure becomes a tab of the one figure window

%% ---- sources to present -----------------------------------------------------------
show_sources = ["CST" "ADS" "Simulink"];          % any of "CST", "ADS", "Simulink"; reference = first of CST, ADS, Simulink shown

%% ---- antenna, excitation, load case ---------------------------------------------
antenna_case = "rectangularpatch";          % "rectangularpatch" | "circularloop"
components   = ["Theta" "Phi"];
switch antenna_case
    case "rectangularpatch"
        antenna_name = "rectangularpatch";
        cst_variant  = "symmetry_DAM";
        cst_fc       = "";                % CST project folder suffix ("" = none)
        excitation   = "filteredsine";
        band         = "2GHz_3GHz";             % band tag in the file names
        load_case    = "OC";                    % internal (load) port: "OC" open, "SC" short
        port         = 2;                       % the load port that is open / shorted in the export folder name
        directions   = [30 90; 45 135; 90 90];  % [theta phi] in deg; must exist in the CST (and ADS) files
        snp_name     = "S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p";
        radial_cuts  = struct('cut', {"theta", "phi"}, 'angle', {90, 45});   % theta cut at phi = 90/270, phi cut at theta = 45
    case "circularloop"
        antenna_name = "circularloop";
        cst_variant  = "FDTD";
        cst_fc       = "";
        excitation   = "gaussiansine";
        band         = "0.2GHz_2GHz";
        load_case    = "";                      % no internal (load) port: exports directly in Export\time_signals
        port         = [];
        directions   = [45 0; 90 0; 45 90];     % CST probes: theta 0, 45, 90, 135, 180 at phi 0, 90 (and 180, 270 for theta 45..135)
        snp_name     = "S_dom_psi_circularloop_FDTD_CST.s8p";
        radial_cuts  = struct('cut', {"theta", "phi"}, 'angle', {0, 90});    % theta cut at phi = 0/180, phi cut at theta = 90
    otherwise
        error('Unknown antenna_case "%s".', antenna_case);
end
% radial_cuts: the directions of the cut present in the CST / ADS data shown (Simulink: those of CST, and of ADS if shown); [] = none

%% ---- files ------------------------------------------------------------------------
tx_snp   = fullfile(snp_dir, snp_name);
cst_proj = antenna_name + "_" + cst_variant; if strlength(cst_fc) > 0, cst_proj = cst_proj + "_" + cst_fc; end
cst_dir  = fullfile(cst_root, cst_proj, "Export", "time_signals");
if ~isempty(port), cst_dir = fullfile(cst_dir, "port" + port + "_" + load_case); end
case_tag = "";  if strlength(load_case) > 0, case_tag = "_" + load_case; end
if strlength(load_case) > 0, load_txt = "load " + load_case; else, load_txt = "no load port"; end
ads_file = fullfile(ads_dir, antenna_name + case_tag + "_rE_" + excitation + "_" + band + ".csv");
cst_ff   = fullfile(cst_dir, ["ff_Etheta_" "ff_Ephi_"] + excitation + "_" + band + ".txt");
cst_io   = fullfile(cst_dir, "iowaves_" + excitation + "_" + band + ".txt");   % Tx waveform (signal i1)

show_sources = string(show_sources);
bad = setdiff(show_sources, ["CST" "ADS" "Simulink"]);
if ~isempty(bad) || isempty(show_sources), error('show_sources: choose from "CST", "ADS", "Simulink" (got %s).', strjoin(show_sources, ', ')); end
useC = any(show_sources == "CST");  useA = any(show_sources == "ADS");  useS = any(show_sources == "Simulink");
order = ["CST" "ADS" "Simulink"];  shown = order(ismember(order, show_sources));   % fixed order, reference first
ref = shown(1);
need = cst_ff;
if useS, need = [need, cst_io, string(tx_snp)]; end
if useA, need = [need, ads_file]; end
for f = need
    if ~isfile(f), error('File not found: %s', f); end
end

%% ---- settings --------------------------------------------------------------------------
cst_signal      = "i1";     % port incident wave used as Tx waveform
cst_excitation  = 1;        % excitation [n] of the CST far-field export
cst_length_unit = 1e-3;     % m
angle_tol       = 0.1;      % deg
save_plots      = true;     % also write every figure as PNG (and .fig) into plot_dir
rebuild         = true;     % false: reuse model_name, only update workspace variables
platform        = "rf";         % "rf": RF Blockset circuit + S-parameters block (name gets _RF; used for the paper);  "simscape": Simscape Electrical + rational-fit block
rf_fitopt       = "Share all poles";   % platform "rf" only: "Share all poles" | "Share poles by columns" | "Fit individually". If the simulation is inaccurate, try by columns or individually (slower), or tighten fit_tol_dB
% Simulink solver (RF platform), user choice: "ode45" = non-stiff (linear network) | "ode15s" = stiff (switches, time-varying R/L/C).
% The RF circuit itself is stepped at MaxStep (default Nyquist/6) whichever solver is chosen.
rf_solver       = "ode45";
model_name      = '';          % '' = Tx<antenna>[_internal_ports] (the ADS cell name without the Rx antenna), chosen by the builder
rf_tag = "";
if platform == "rf" && useS, rf_tag = "_RF"; end
fit_tol_dB      = -40;
bw_level_dB     = 20;       % in-band magnitude error: reference spectrum within this many dB of its peak
f_plot_GHz      = [0 5];    % spectrum plot range
Rload_OC = 1e9;  Rload_SC = 1e-6;     % Ohm used for "open" and "short" (internal ports only)
if load_case == "SC", Rload = Rload_SC; else, Rload = Rload_OC; end
tag    = antenna_name + case_tag + rf_tag + "_" + strjoin(shown, "_");   % file-name tag of the plots
srcTxt = strjoin(shown, " vs ");
% one colour code for every plot (as in the radial plots): CST black, ADS red dotted, Simulink light blue, wide, underneath
sty = struct('CST',      struct('col', [0 0 0],       'ls', '-', 'lw', 1.2), ...
             'ADS',      struct('col', [0.85 0 0],    'ls', ':', 'lw', 1.6), ...
             'Simulink', struct('col', [0.55 0.75 1], 'ls', '-', 'lw', 3.5));
drawOrder = [find(shown == "Simulink"), find(shown ~= "Simulink")];   % Simulink underneath

%% ---- read CST and ADS ---------------------------------------------------------------------
cst = readCstFarfieldProbes(cst_ff, components, cst_excitation, cst_length_unit);
if useA, ads = readAdsFarfieldSweeps(ads_file, components, 0); else, ads = []; end
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
nD = size(directions, 1);
iC = zeros(nD, 1);  iA = zeros(nD, 1);
for d = 1:nD
    th = directions(d, 1);  ph = directions(d, 2);
    j = find(arrayfun(@(c) abs(c.theta - th) < angle_tol && angd(c.phi, ph) < angle_tol, cst), 1);
    k = 0;
    if useA, k = find(arrayfun(@(a) abs(a.theta - th) < angle_tol && angd(a.phi, ph) < angle_tol, ads), 1); end
    if isempty(j) || isempty(k)
        error('Direction (theta %g, phi %g) missing in %s.', th, ph, ternary(isempty(j), 'the CST export', 'the ADS file'));
    end
    iC(d) = j;  iA(d) = k;
end
r_probe = cst(iC(1)).r;
fprintf('Antenna %s, %s, r = %g m, directions: %s\nShown: %s (reference %s)\n', antenna_name, load_txt, r_probe, ...
    mat2str(directions), srcTxt, ref);

%% ---- build / update the Simulink model and simulate (only if Simulink is shown) ------------------------
sim_sig = cell(nD, 1);
if useS
    if rebuild
        bargs = {'TxOnly', true, 'Waveform', 'cst', 'CSTFile', char(cst_io), 'CSTSignal', char(cst_signal), 'r', r_probe, ...
            'theta_deg', directions(1, 1), 'phi_deg', directions(1, 2), 'TxLoadR', Rload, 'FitTolerance', fit_tol_dB, ...
            'OpenModel', false, 'Platform', char(platform), 'RFSolver', char(rf_solver), 'FitOptions', char(rf_fitopt), ...
            'OutputDir', model_dir};
        if ~isempty(model_name), bargs = [bargs, {'ModelName', model_name}]; end
        mdl = build_tx_rx_link_simulink(tx_snp, '', bargs{:});
    else
        mdl = model_name;
        if isempty(mdl), error('Set model_name to reuse a model (rebuild = false).'); end
        if ~bdIsLoaded(mdl), load_system(fullfile(model_dir, mdl)); end
    end
    assignin('base', 'Rload_Tx', Rload);
    wave = evalin('base', 'Vtx_wave');
    tend = wave.Time(end);
    for d = 1:nD
        sim_sig{d} = simulateDirection(mdl, directions(d, 1), directions(d, 2), r_probe, tend);
    end
end

%% ---- compare ----------------------------------------------------------------------------------------
results = struct('theta', {}, 'phi', {}, 'comp', {}, 'ref', {}, 'sources', {}, 'peak', {}, 'rel_mse', {}, ...
    'lag_ps', {}, 'magerr_dB', {});
nS = numel(shown);
for d = 1:nD
    cp = cst(iC(d));
    for c = 1:numel(components)
        % records of the shown sources for this direction / component
        T = cell(1, nS);  E = cell(1, nS);
        for s = 1:nS
            switch shown(s)
                case "CST",      T{s} = cp.t{c}(:);  E{s} = cp.e{c}(:);
                case "ADS",      T{s} = ads(iA(d)).t(:);  E{s} = ads(iA(d)).e(:, c);
                case "Simulink", T{s} = sim_sig{d}.t(:);  E{s} = sim_sig{d}.e(:, c);
            end
            if isempty(T{s}), continue, end
            [T{s}, iu] = unique(T{s}, 'last');  E{s} = E{s}(iu);
        end
        if any(cellfun(@isempty, T)), continue, end
        pk = max(cellfun(@(x) max(abs(x)), E));
        pkDir = max(cellfun(@(x) max([abs(x(:)); 0]), cp.e));       % strongest CST component in this direction
        if pk < 1e-4*pkDir || pk < 1e-9                             % null component (numerical noise only)
            fprintf('theta %g phi %g: E_%s is zero (peak %.2g V, %.0f dB below E of this direction), skipped\n', ...
                directions(d, :), lower(components(c)), pk, -20*log10(max(pk, realmin)/pkDir));
            continue
        end
        % common uniform grid over the overlap of the records
        tB = max(cellfun(@(x) x(1), T));  tE = min(cellfun(@(x) x(end), T));
        dt = min(cellfun(@(x) median(diff(x)), T)) / 4;
        tg = (tB:dt:tE).';
        X = zeros(numel(tg), nS);
        for s = 1:nS, X(:, s) = interp1(T{s}, E{s}, tg, 'pchip', 0); end
        nrm = @(x) trapz(tg, x.^2);
        relmse = @(x, y) nrm(x - y) / nrm(y);
        xr = X(:, 1);
        rm = nan(1, nS);  lagps = nan(1, nS);  eM = nan(1, nS);  fl = NaN;  fh = NaN;
        for s = 2:nS
            rm(s) = relmse(X(:, s), xr);
            [rr, lags] = xcorr(xr, X(:, s), round(0.5e-9/dt));  [~, im] = max(rr);  lagps(s) = lags(im)*dt*1e12;
            [eM(s), fl, fh] = bandMagError(tg, xr, X(:, s), bw_level_dB);
        end
        results(end+1) = struct('theta', directions(d, 1), 'phi', directions(d, 2), 'comp', components(c), 'ref', ref, ...
            'sources', shown, 'peak', max(abs(X), [], 1), 'rel_mse', rm, 'lag_ps', lagps, 'magerr_dB', eM); %#ok<SAGROW>

        % ---- figure row: time, difference, spectrum ----
        fig = figure('Color', 'w', 'WindowStyle', 'docked', 'Name', sprintf('E_%s theta %g phi %g (%s)', lower(components(c)), ...
            directions(d, :), strjoin(shown, '/')), 'NumberTitle', 'off');
        tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
        title(tl, sprintf('$rE_{\\%s}$, $\\theta = %g^\\circ$, $\\phi = %g^\\circ$, $r = %g$ m, %s: %s', ...
            lower(components(c)), directions(d, :), r_probe, strrep(antenna_name, '_', '\_') + ", " + load_txt, srcTxt), 'Interpreter', 'latex');
        tmaxns = tE*1e9;        % plots end at the minimum of the last times of the records
        nexttile(tl); hold on
        h = gobjects(1, nS);
        for s = drawOrder, st = sty.(shown(s));  h(s) = plot(tg*1e9, X(:, s), st.ls, 'Color', st.col, 'LineWidth', st.lw); end
        grid on; box on; xlabel('Time (ns)'); ylabel('$rE$ (V)');  legend(h, cellstr(shown), 'Location', 'best');
        xlim([tg(1)*1e9 tmaxns]);
        if nS > 1
            title(sprintf('rel. MSE vs %s: %s', ref, strjoin(compose('%s %.2g\\%%', shown(2:end).', 100*rm(2:end).'), ', ')));
        else
            title(ref);
        end
        nexttile(tl); hold on
        for s = 2:nS, st = sty.(shown(s));  plot(tg*1e9, X(:, s) - xr, '-', 'Color', st.col, 'LineWidth', min(st.lw, 2)); end
        grid on; box on; xlabel('Time (ns)'); ylabel(sprintf('difference to %s (V)', ref));
        if nS > 1, legend(cellstr(shown(2:end) + " - " + ref), 'Location', 'best'); end
        xlim([tg(1)*1e9 tmaxns]);
        title('differences');
        nexttile(tl); hold on
        Nf = 2^nextpow2(8*numel(tg));  f = (0:Nf/2).'/(Nf*dt);
        Xr = max(abs(fft(xr, Nf)));
        h = gobjects(1, nS);
        for s = drawOrder
            Xs = abs(fft(X(:, s), Nf));  Xs = Xs(1:Nf/2 + 1);  st = sty.(shown(s));
            h(s) = plot(f/1e9, 20*log10(max(Xs, eps)/Xr), st.ls, 'Color', st.col, 'LineWidth', st.lw);
        end
        xlim(f_plot_GHz);  ylim([-80 5]);  grid on; box on; xlabel('Frequency (GHz)'); ylabel(sprintf('dB re %s peak', ref));
        legend(h, cellstr(shown), 'Location', 'south');
        if nS > 1
            title(sprintf('in-band $|X|$ error vs %s: %s (%.3g--%.3g GHz)', ref, ...
                strjoin(compose('%s %.2f dB', shown(2:end).', eM(2:end).'), ', '), fl/1e9, fh/1e9));
        else
            title('spectrum');
        end
        % (no Position change: setting it would undock the figure from the tabbed window)
        if save_plots
            if ~isfolder(plot_dir), mkdir(plot_dir); end
            stem = sprintf('compare_%s_E%s_th%g_ph%g', tag, lower(components(c)), directions(d, :));
            ax_ = findall(fig, 'Type', 'axes');  for a_ = ax_.', if ~isempty(a_.Toolbar), a_.Toolbar.Visible = 'off'; end, end   % no axes toolbar in the PNG
            exportgraphics(fig, fullfile(plot_dir, stem + ".png"), 'Resolution', 120);
            savefig(fig, fullfile(plot_dir, stem + ".fig"));
        end
    end
end

%% ---- radial (spoke) plots over a cut ------------------------------------------------
if ~isempty(radial_cuts)
    simAll = struct('theta', {}, 'phi', {}, 't', {}, 'e', {});
    if useS, for d = 1:nD, simAll(end+1) = sim_sig{d}; end, end %#ok<SAGROW>
    for q = 1:numel(radial_cuts)
        cutq = radial_cuts(q).cut;  angq = radial_cuts(q).angle;
        dirs = cutDirections(cst, ads, cutq, angq, angle_tol);
        if useS
            for i = 1:size(dirs, 1)
                have = arrayfun(@(s) abs(s.theta - dirs(i, 1)) < angle_tol && angd(s.phi, dirs(i, 2)) < angle_tol, simAll);
                if ~any(have)
                    simAll(end+1) = simulateDirection(mdl, dirs(i, 1), dirs(i, 2), r_probe, tend); %#ok<SAGROW>
                end
            end
            a_tx = evalin('base', 'a_Tx');
        else
            a_tx = readCstA(cst_ff(1), cst_length_unit);
        end
        sf = '';
        if save_plots, sf = fullfile(plot_dir, sprintf('radial_%s_%s%g.png', tag, cutq, angq)); end
        plot_rE_radial_simulink(simAll, cst, ads, 'Sources', shown, 'Cut', cutq, 'Angle', angq, 'r', r_probe, 'a', a_tx, ...
            'AngleTol', angle_tol, 'SaveFile', sf, ...
            'Title', sprintf('%s, %s, %s %s: %s', strrep(antenna_name, '_', '\_'), load_txt, excitation, strrep(band, '_', '--'), srcTxt));
    end
end

%% ---- summary -------------------------------------------------------------------------------------------------
if save_plots, fprintf('\nPlots saved in %s\n', plot_dir); end
fprintf('\nsummary (errors relative to %s over the overlap of the records; %s)\n', ref, load_txt);
hdr = sprintf('  theta   phi  comp ');
hdr = [hdr, sprintf('  peak %-8s', shown{:})];
for s = 2:nS, hdr = [hdr, sprintf('  relMSE %-8s  delay %-8s  in-band %-8s', shown(s), shown(s), shown(s))]; end %#ok<AGROW>
fprintf('%s\n', hdr);
for j = 1:numel(results)
    q = results(j);
    ln = sprintf('  %5g %5g  %-5s', q.theta, q.phi, q.comp);
    ln = [ln, sprintf('  %11.4g V', q.peak)]; %#ok<AGROW>
    for s = 2:nS, ln = [ln, sprintf('  %13.3f %%  %11.1f ps  %13.2f dB', 100*q.rel_mse(s), q.lag_ps(s), q.magerr_dB(s))]; end %#ok<AGROW>
    fprintf('%s\n', ln);
end

% make sure every figure of this run sits in the tabbed figure window
set(findall(groot, 'Type', 'figure'), 'WindowStyle', 'docked');
set(groot, 'defaultTextInterpreter', 'factory', 'defaultAxesTickLabelInterpreter', 'factory', 'defaultLegendInterpreter', 'factory');

%% ---- local functions -------------------------------------------------------------------------------------------
function s = simulateDirection(mdl, th, ph, r_probe, tend)
% simulate the Tx-only model for the observation direction (theta, phi) at distance r_probe and
% return r*E_theta, r*E_phi (Vtheta, Vphi) in V
assignin('base', 'r', r_probe);
assignin('base', 'theta', th*pi/180);
assignin('base', 'phi', ph*pi/180);
link_update_gains;
pt = evalin('base', 'prop_t');
assignin('base', 't_exc', tend);  set_param(mdl, 'StopTime', 'prop_t + t_exc');     % excitation duration + propagation delay
fprintf('theta = %g, phi = %g: prop_t = %.4g ns, simulating to %.4g ns...\n', th, ph, pt*1e9, (pt + tend)*1e9);
out = sim(mdl);
vth = out.yout.getElement('Vtheta').Values;
vph = out.yout.getElement('Vphi').Values;
s = struct('theta', th, 'phi', ph, 't', vth.Time, 'e', [vth.Data(:), vph.Data(:)]);
end

function dirs = cutDirections(cst, ads, cut, ang, tol)
% directions [theta phi] that lie in the cut and are present in the CST data (and in the ADS data, if any)
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
dirs = zeros(0, 2);
for j = 1:numel(cst)
    th = cst(j).theta;  ph = cst(j).phi;
    if string(cut) == "theta"
        inCut = angd(ph, ang) < tol || angd(ph, ang + 180) < tol;
    else
        inCut = abs(th - ang) < tol;
    end
    if inCut && (isempty(ads) || any(arrayfun(@(a) abs(a.theta - th) < tol && angd(a.phi, ph) < tol, ads)))
        dirs(end+1, :) = [th ph]; %#ok<AGROW>
    end
end
end

function a = readCstA(file, unit)
% enclosing-sphere radius a (m) from the '#Parameters = {... a=40; ...}' line of a CST export (0 if absent)
a = 0;
L = readlines(file);  L = L(startsWith(strtrim(L), "#Parameters"));
if isempty(L), return, end
tok = regexp(L(1), '[{;\s]a\s*=\s*([-+\deE.]+)', 'tokens', 'once');
if ~isempty(tok), a = str2double(tok{1}) * unit; end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function [eM, fLo, fHi] = bandMagError(tg, xRef, x, lev)
% mean |20 log10(|X|/|Xref|)| over the band where the reference spectrum is within lev dB of its peak
dt = tg(2) - tg(1);
Nf = 2^nextpow2(8*numel(tg));  f = (0:Nf/2).'/(Nf*dt);
Xr = abs(fft(xRef, Nf));  Xr = Xr(1:Nf/2 + 1);
Xx = abs(fft(x, Nf));     Xx = Xx(1:Nf/2 + 1);
[pk, ip] = max(Xr);
thr = pk*10^(-lev/20);
il = find(Xr(1:ip) < thr, 1, 'last');  if isempty(il), il = 1; end
iu = ip - 1 + find(Xr(ip:end) < thr, 1, 'first');  if isempty(iu), iu = numel(f); end
m = il:iu;
eM = mean(abs(20*log10(max(Xx(m), eps) ./ max(Xr(m), eps))));
fLo = f(il);  fHi = f(iu);
end
