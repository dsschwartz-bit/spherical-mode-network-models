%COMPARE_RADIATED_FIELDS  Compare radiated far-field time signals from CST, ADS and Simulink.
% Radiated far-field time signals r*E_theta(t), r*E_phi(t) of an antenna from
% any of CST, ADS and Simulink, compared at observation directions (theta, phi) at the CST
% probe distance r, and over radial (spoke) cuts.
%
%   show_sources : which data to present, any of "CST", "ADS", "Simulink".  Errors are
%                  relative to the reference: the first of CST, ADS, Simulink shown.
%   antenna_case : a block of settings per antenna (files, directions, cuts); add a case
%                  for another antenna.
%
% Data:
%   CST      : far-field time-signal exports ff_Etheta_* / ff_Ephi_* in ...\Export\time_signals\
%              (with a load port: ...\time_signals\port<N>_<OC|SC>\), readCstFarfieldProbes.m.
%              Always read (probe distance r and the available directions).
%   ADS      : CSV Data Display export <antenna>[_<OC|SC>]_rE_<excitation>_<band>.csv
%              (readAdsFarfieldSweeps.m).
%   Simulink : Tx-only spherical-mode model (matlab/simulink/build_tx_rx_link_simulink.m) driven
%              by the CST port signal i1, results kept in output/cst_ads/simulink_out/<...>.mat.  Directions
%              not in that file are SIMULATED AUTOMATICALLY (the model is built first) and
%              added to it (simulinkRadiatedField.m); simulink_mode = "rerun" simulates all again.
%
% Output: per direction and component one figure (time signals, differences to the reference,
% spectra), radial plots per cut (plot_rE_radial_simulink.m), a summary table in the console,
% and PNG + .fig files in output/cst_ads/plots (save_plots).  Figures are docked in the figure window.
%
% Errors over the overlap of the records: relative MSE = integral (x - ref)^2 dt / integral ref^2 dt,
% best-fit delay, and the in-band (-bw_level_dB) mean spectrum magnitude error.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));     % adds matlab/ and third_party/
addpath(fullfile(repo_root, 'validation', 'simulink'));  % plot_rE_radial_simulink
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');   % LaTeX fonts
set(groot, 'defaultFigureWindowStyle', 'docked');     % every new figure becomes a tab of the one figure window

%% ===== USER SETTINGS =====
cst_root        = "C:\path\to\CST\";                    % folder containing the CST project folders
ads_dir         = "C:\path\to\ADS\workspace\data\";   % ADS workspace data folder (CSV exports)
plot_dir        = fullfile(repo_root, 'output', 'cst_ads', 'plots');   % saved figures (save_plots)
% =========================

%% ---- what to present -------------------------------------------------------------------
show_sources  = ["CST" "ADS"];   % any of "CST", "ADS", "Simulink"
simulink_mode = "auto";      % "auto": saved Simulink results, simulate missing directions;
                             % "rerun": simulate everything again;  "cache": saved results only

%% ---- antenna, excitation, load case ---------------------------------------------------
antenna_case = "circularloop";          % "rectangularpatch" | "circularloop"
components   = ["Theta" "Phi"];
switch antenna_case
    case "rectangularpatch"
        antenna_name = "rectangularpatch";
        cst_variant  = "symmetry_FDTD";
        cst_fc       = "2.5GHz";                % CST project folder suffix ("" = none)
        excitation   = "gaussian";
        band         = "2GHz_3GHz";             % band tag in the file names
        load_case    = "OC";                    % internal (load) port: "OC" open, "SC" short
        port         = 2;                       % the load port that is open / shorted in the export folder name
        directions   = [0 0; 45 135; 90 90];  % [theta phi] in deg; must exist in the CST (and ADS) files
        snp_name     = "S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p";
        radial_cuts  = struct('cut', {"theta", "phi"}, 'angle', {90, 45}, 'spokes', {[], []});   % theta cut at phi = 90/270, phi cut at theta = 45
    case "circularloop"
        antenna_name = "circularloop";
        cst_variant  = "FDTD";
        cst_fc       = "";
        excitation   = "gaussiansine";
        band         = "0.2GHz_2GHz";
        load_case    = "";                      % no internal (load) port: exports directly in Export\time_signals
        port         = [];
        directions   = [0 45; 0 90 ; 90 90 ];     % CST probes: theta 0, 45, 90, 135, 180 at phi 0, 90 (and 180, 270 for theta 45..135)
        snp_name     = "S_dom_psi_circularloop_FDTD_CST.s8p";
        radial_cuts  = struct('cut', {"theta", "phi"}, 'angle', {0, 90}, 'spokes', {[0 45 90], []});   % theta cut at phi = 0, spokes theta = 0, 45, 90 (no phi = 180 side); phi cut at theta = 90
        % radial_cuts = struct('cut', {"theta", "phi"}, 'angle', {0, 90}, 'spokes', {[], []});   % alternative: full theta cut and phi cut at theta = 90
    otherwise
        error('Unknown antenna_case "%s".', antenna_case);
end
% radial_cuts: spokes = directions of the cut in the CST data (and in the ADS data, if shown); [] = none.
%   'spokes' limits them: signed theta values for a theta cut (deg; negative = the half-plane phi + 180,
%   e.g. [-90 -45 0 45 90]), phi values for a phi cut; [] = all

%% ---- settings ------------------------------------------------------------------------------
cst_signal      = "i1";     % port incident wave: the Simulink Tx waveform
cst_excitation  = 1;        % excitation [n] of the CST far-field export
cst_length_unit = 1e-3;     % m
angle_tol       = 0.1;      % deg
save_plots      = true;     % write every figure as PNG and .fig into plot_dir
bw_level_dB     = 20;       % in-band magnitude error: reference spectrum within this many dB of its peak
f_plot_GHz      = [0 5];    % spectrum plot range
% Simulink model
platform        = "rf";         % "rf": RF Blockset circuit + S-parameters block (used for the paper);  "simscape": Simscape Electrical + rational fit
rf_fitopt       = "Share all poles";   % platform "rf": "Share all poles" | "Share poles by columns" | "Fit individually"
% Simulink solver (RF platform), user choice: "ode45" = non-stiff (linear network) | "ode15s" = stiff (switches, time-varying R/L/C).
% The RF circuit itself is stepped at MaxStep (default Nyquist/6) whichever solver is chosen.
rf_solver       = "ode45";
fit_tol_dB      = -40;          % rational-fit tolerance
Rload_OC = 1e9;  Rload_SC = 1e-6;     % Ohm for "open" / "short" internal ports

%% ---- files ---------------------------------------------------------------------------------
cst_proj = antenna_name + "_" + cst_variant; if strlength(cst_fc) > 0, cst_proj = cst_proj + "_" + cst_fc; end
cst_dir  = string(fullfile(cst_root, cst_proj, "Export", "time_signals")) + filesep;
if ~isempty(port), cst_dir = cst_dir + "port" + port + "_" + load_case + filesep; end
case_tag = "";  if strlength(load_case) > 0, case_tag = "_" + load_case; end
if strlength(load_case) > 0, load_txt = "load " + load_case; else, load_txt = "no load port"; end
ads_file = string(fullfile(ads_dir, antenna_name + case_tag + "_rE_" + excitation + "_" + band + ".csv"));
cst_ff   = cst_dir + ["ff_Etheta_" "ff_Ephi_"] + excitation + "_" + band + ".txt";
cst_io   = cst_dir + "iowaves_" + excitation + "_" + band + ".txt";
tx_snp   = fullfile(repo_root, 'output', 'SnP', snp_name);   % spherical-mode Touchstone model
if load_case == "SC", Rload = Rload_SC; else, Rload = Rload_OC; end
rf_tag = "";  if platform == "rf", rf_tag = "_RF"; end
sim_cache = fullfile(repo_root, 'output', 'cst_ads', 'simulink_out', antenna_name + case_tag + "_rE_" + excitation + "_" + band + rf_tag + ".mat");

show_sources = string(show_sources);
bad = setdiff(show_sources, ["CST" "ADS" "Simulink"]);
if ~isempty(bad) || isempty(show_sources), error('show_sources: choose from "CST", "ADS", "Simulink" (got %s).', strjoin(show_sources, ', ')); end
order = ["CST" "ADS" "Simulink"];  shown = order(ismember(order, show_sources));   % fixed order, reference first
useA = any(shown == "ADS");  useS = any(shown == "Simulink");
ref = shown(1);  nS = numel(shown);
need = cst_ff;
if useS, need = [need, cst_io, string(tx_snp)]; end
if useA, need = [need, ads_file]; end
for f = need
    if ~isfile(f), error('File not found: %s', f); end
end
tag    = antenna_name + case_tag + "_" + strjoin(shown, "_");     % file-name tag of the plots
srcTxt = strjoin(shown, " vs ");
% one colour code for every plot (as in the radial plots): CST black, ADS red dotted, Simulink light blue, wide,
% drawn first (underneath) so that all three stay visible where they overlap
sty = struct('CST',      struct('col', [0 0 0],       'ls', '-', 'lw', 1.2), ...
             'ADS',      struct('col', [0.85 0 0],    'ls', ':', 'lw', 1.6), ...
             'Simulink', struct('col', [0.55 0.75 1], 'ls', '-', 'lw', 3.5));
drawOrder = [find(shown == "Simulink"), find(shown ~= "Simulink")];   % Simulink underneath

%% ---- read CST and ADS --------------------------------------------------------------------
cst = readCstFarfieldProbes(cst_ff, components, cst_excitation, cst_length_unit);
if useA, ads = readAdsFarfieldSweeps(ads_file, components, 0); else, ads = []; end
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
findDir = @(S, th, ph) find(arrayfun(@(c) abs(c.theta - th) < angle_tol && angd(c.phi, ph) < angle_tol, S), 1);
nD = size(directions, 1);
iC = zeros(nD, 1);  iA = zeros(nD, 1);
for d = 1:nD
    j = findDir(cst, directions(d, 1), directions(d, 2));
    k = 0;  if useA, k = findDir(ads, directions(d, 1), directions(d, 2)); end
    if isempty(j) || isempty(k)
        warning('Direction (theta %g, phi %g) is not in %s: skipped.', directions(d, 1), directions(d, 2), ternary(isempty(j), 'the CST export', 'the ADS file'));
        continue
    end
    iC(d) = j;  iA(d) = k;
end
ok = iC > 0;  directions = directions(ok, :);  iC = iC(ok);  iA = iA(ok);  nD = size(directions, 1);
r_probe = cst(1).r;
a_tx    = readCstA(cst_ff(1), cst_length_unit);
fprintf('Antenna %s, %s, r = %g m\nShown: %s (reference %s)\n', antenna_name, load_txt, r_probe, srcTxt, ref);

% every direction needed: the listed ones and those of the radial cuts
cutDirs = cell(1, numel(radial_cuts));
allDirs = directions;
for q = 1:numel(radial_cuts)
    cutDirs{q} = cutDirections(cst, ads, radial_cuts(q).cut, radial_cuts(q).angle, angle_tol, radial_cuts(q).spokes);
    allDirs = [allDirs; cutDirs{q}]; %#ok<AGROW>
end
[~, iu] = uniquetol([allDirs(:, 1), mod(allDirs(:, 2), 360)], angle_tol, 'ByRows', true, 'DataScale', 1);
allDirs = allDirs(sort(iu), :);

%% ---- Simulink: saved results, missing directions simulated -------------------------------
simAll = struct('theta', {}, 'phi', {}, 't', {}, 'e', {});
if useS
    sopt = struct('tx_snp', tx_snp, 'cst_io', cst_io, 'cst_signal', cst_signal, 'r', r_probe, 'Rload', Rload, ...
        'fit_tol_dB', fit_tol_dB, 'platform', platform, 'rf_solver', rf_solver, 'rf_fitopt', rf_fitopt, ...
        'cache_file', sim_cache, 'mode', simulink_mode, 'angle_tol', angle_tol);
    simAll = simulinkRadiatedField(allDirs, sopt);
end

%% ---- compare per direction ------------------------------------------------------------------
results = struct('theta', {}, 'phi', {}, 'comp', {}, 'ref', {}, 'sources', {}, 'peak', {}, 'rel_mse', {}, ...
    'lag_ps', {}, 'magerr_dB', {});
for d = 1:nD
    cp = cst(iC(d));
    for c = 1:numel(components)
        T = cell(1, nS);  E = cell(1, nS);
        for s = 1:nS
            switch shown(s)
                case "CST",      T{s} = cp.t{c}(:);  E{s} = cp.e{c}(:);
                case "ADS",      T{s} = ads(iA(d)).t(:);  E{s} = ads(iA(d)).e(:, c);
                case "Simulink", m = findDir(simAll, directions(d, 1), directions(d, 2));
                                 T{s} = simAll(m).t(:);  E{s} = simAll(m).e(:, c);
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
        xr = X(:, 1);
        rm = nan(1, nS);  lagps = nan(1, nS);  eM = nan(1, nS);  fl = NaN;  fh = NaN;
        for s = 2:nS
            rm(s) = nrm(X(:, s) - xr) / nrm(xr);
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
        nexttile(tl); hold on
        h = gobjects(1, nS);
        for s = drawOrder, st = sty.(shown(s));  h(s) = plot(tg*1e9, X(:, s), st.ls, 'Color', st.col, 'LineWidth', st.lw); end
        grid on; box on; xlabel('Time (ns)'); ylabel('$rE$ (V)');  legend(h, cellstr(shown), 'Location', 'best');
        xlim([tg(1) tE]*1e9);
        if nS > 1, title(sprintf('rel. MSE vs %s: %s', ref, strjoin(compose('%s %.2g\\%%', shown(2:end).', 100*rm(2:end).'), ', ')));
        else,      title(ref); end
        nexttile(tl); hold on
        for s = 2:nS, st = sty.(shown(s));  plot(tg*1e9, X(:, s) - xr, '-', 'Color', st.col, 'LineWidth', min(st.lw, 2)); end
        grid on; box on; xlabel('Time (ns)'); ylabel(sprintf('difference to %s (V)', ref));
        if nS > 1, legend(cellstr(shown(2:end) + " - " + ref), 'Location', 'best'); end
        xlim([tg(1) tE]*1e9);  title('differences');
        nexttile(tl); hold on
        Nf = 2^nextpow2(8*numel(tg));  f = (0:Nf/2).'/(Nf*dt);
        Xr = max(abs(fft(xr, Nf)));
        h = gobjects(1, nS);
        for s = drawOrder
            Xs = abs(fft(X(:, s), Nf));  Xs = Xs(1:Nf/2 + 1);  st = sty.(shown(s));
            h(s) = plot(f/1e9, 20*log10(max(Xs, eps)/Xr), st.ls, 'Color', st.col, 'LineWidth', st.lw);
        end
        xlim(f_plot_GHz);  ylim([-80 5]);  grid on; box on; xlabel('Frequency (GHz)'); ylabel(sprintf('dB rel. %s peak', ref));
        legend(h, cellstr(shown), 'Location', 'south');
        if nS > 1
            title(sprintf('in-band $|X|$ error vs %s: %s (%.3g--%.3g GHz)', ref, ...
                strjoin(compose('%s %.2f dB', shown(2:end).', eM(2:end).'), ', '), fl/1e9, fh/1e9));
        else
            title('spectrum');
        end
        if save_plots
            if ~isfolder(plot_dir), mkdir(plot_dir); end
            stem = sprintf('compare_%s_E%s_th%g_ph%g', tag, lower(components(c)), directions(d, :));
            ax_ = findall(fig, 'Type', 'axes');  for a_ = ax_.', if ~isempty(a_.Toolbar), a_.Toolbar.Visible = 'off'; end, end
            exportgraphics(fig, fullfile(plot_dir, stem + ".png"), 'Resolution', 120);
            savefig(fig, fullfile(plot_dir, stem + ".fig"));
        end
    end
end

%% ---- radial (spoke) plots over the cuts -------------------------------------------------------
for q = 1:numel(radial_cuts)
    cutq = radial_cuts(q).cut;  angq = radial_cuts(q).angle;
    sf = '';
    if save_plots, sf = fullfile(plot_dir, sprintf('radial_%s_%s%g.png', tag, cutq, angq)); end
    plot_rE_radial_simulink(simAll, cst, ads, 'Sources', shown, 'Cut', cutq, 'Angle', angq, 'Spokes', radial_cuts(q).spokes, 'r', r_probe, 'a', a_tx, ...
        'AngleTol', angle_tol, 'SaveFile', sf, ...
        'Title', sprintf('%s, %s, %s %s: %s', strrep(antenna_name, '_', '\_'), load_txt, excitation, strrep(band, '_', '--'), srcTxt));
end

%% ---- summary ------------------------------------------------------------------------------------
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

set(findall(groot, 'Type', 'figure'), 'WindowStyle', 'docked');
set(groot, 'defaultTextInterpreter', 'factory', 'defaultAxesTickLabelInterpreter', 'factory', 'defaultLegendInterpreter', 'factory');

%% ---- local functions ----------------------------------------------------------------------------
function dirs = cutDirections(cst, ads, cut, ang, tol, spokes)
% directions [theta phi] that lie in the cut and are present in the CST data (and in the ADS data, if any);
% spokes (optional): only these signed theta values (theta cut; negative = half-plane phi + 180) or phi values (phi cut)
if nargin < 6, spokes = []; end
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
dirs = zeros(0, 2);
for j = 1:numel(cst)
    th = cst(j).theta;  ph = cst(j).phi;
    if string(cut) == "theta"
        inCut = angd(ph, ang) < tol || angd(ph, ang + 180) < tol;
        if angd(ph, ang) < tol, us = th; else, us = -th; end      % signed theta along the cut
        if ~isempty(spokes), inCut = inCut && any(abs(us - spokes) < tol); end
    else
        inCut = abs(th - ang) < tol;
        if ~isempty(spokes), inCut = inCut && any(angd(ph, spokes) < tol); end
    end
    onAxis = th < tol || abs(th - 180) < tol;            % z axis: the same point for every phi
    if inCut && (isempty(ads) || any(arrayfun(@(a) abs(a.theta - th) < tol && (onAxis || angd(a.phi, ph) < tol), ads)))
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
