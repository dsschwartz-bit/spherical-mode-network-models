%COMPARE_S11_CST_ADS  Compare S11 with port 2 shorted or open: CST, ADS and the spherical-mode model.
% Input reflection coefficient S11 of a two-port antenna (port 1 feed, port 2
% load) with port 2 terminated in a short (SC) or an open (OC), from up to
% four sources:
%
%   CST co-simulation  <cst_sp_dir>\<antenna_name>_<termination>.s1p   (Touchstone)
%   ADS                <ads_dir>\<antenna_name>_<termination>_S11.csv  (freq, mag, phase in deg)
%   CST two-port       each of s2p_files, and
%   spherical-mode     <snp_file>: the feed/load block (snp_ports) of the combined
%   multiport          port + spherical-mode matrix, mode ports terminated in
%                      mode_load (radiating into free space); each with port 2
%                      terminated here:
%                        Z   = Z-parameters of the s2p (removes its own reference)
%                        Zin = Z11 - Z12*Z21/(Z22 + ZL)   (ZL = 0 short; open: Zin = Z11)
%                        S11 = (Zin - z_ref)/(Zin + z_ref)
%
% Every S11 is referred to z_ref. The phase is that of S11 at each source's
% own port-1 reference plane (no de-embedding).
%
% Sources are compared over the overlap of their frequency ranges (limited
% to the plot range), on the coarser of the two grids (the finer one is
% interpolated, complex, onto it), by the relative RMS error
%   rel. RMS error = sqrt( sum |Sa - Sb|^2 / sum |Sb|^2 )     (b = reference)
% printed for every pair and shown in the legends: figure 1 (S11 with port 2
% terminated) compares ADS and the CST results (re fig1_ref), figure 2 (two-port
% S-parameters) re diff_ref. Output: figures and console tables only.
% Requires RF Toolbox.
clear; close all;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
cst_root = "C:\path\to\CST\";                % folder containing the CST project folders
ads_dir  = "C:\path\to\ADS\workspace\data\";   % ADS workspace data folder (CSV exports)
% =========================

termination  = "SC";       % port 2 load: "SC" (short circuit) or "OC" (open circuit)
antenna_name = "rectangularpatch";
cst_sp_dir   = fullfile(cst_root, "rectangularpatch_symmetry_FDTD_2.5GHz", "SP");
s2p_files    = [fullfile(cst_sp_dir, "rectangularpatch_symmetry_FDTD.s2p") ...         % CST two-port results (port 1 feed, port 2 load); [] = none
                fullfile(cst_root, "rectangularpatch_symmetry_FDTD", "SP", "rectangularpatch_symmetry_FDTD.s2p")];
s2p_names    = ["CST 2.5 GHz project" "CST FDTD project"];    % legend label for each s2p file
snp_file     = fullfile(repo_root, "output", "SnP", "S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p");  % spherical-mode model; "" = skip
snp_ports    = [1 2];      % [feed load] port numbers in snp_file
mode_load    = 376.73;     % ohm, load on the other (spherical-mode) ports of snp_file: eta0, as in the ADS model
z_ref        = 50;         % ohm, reference impedance of every S11 shown
ads_z0       = 50;         % ohm, reference impedance of the ADS S11
plot_frange_GHz = [];      % [fmin fmax] GHz; [] = range of the CST co-simulation, [0 Inf] = all data
diff_ref     = "CST FDTD project";   % reference for the RMS error in the figure-2 legends: the source whose name starts with this; "" = none
fig1_ref     = "CST FDTD project";   % figure 1: reference for the RMS error (the source whose name starts with this)
fig1_exclude = ["ADS s2p" "CST co-sim"];   % figure 1: sources whose name starts with any of these are not drawn (still in the printed table)

cst_s1p_file = fullfile(cst_sp_dir, antenna_name + "_" + termination + ".s1p");
ads_file     = fullfile(ads_dir, antenna_name + "_" + termination + "_S11.csv");
termName     = struct('SC', "shorted", 'OC', "open");
if ~isfield(termName, termination), error('termination must be "SC" or "OC".'); end

%% ---- read
src = struct('name', {}, 'f', {}, 'S', {}, 'col', {}, 'ls', {}, 'file', {});
co  = orderedcolors("gem");          % default MATLAB line colours (readable in light and dark themes)

s1 = sparameters(char(cst_s1p_file));
if s1.Impedance ~= z_ref, s1 = sparameters(s1, z_ref); end
src(end+1) = struct('name', "CST co-sim", 'f', s1.Frequencies(:), 'S', squeeze(s1.Parameters(1,1,:)), ...
    'col', co(1,:), 'ls', '-', 'file', cst_s1p_file);

[fA, SA] = readAdsS11(ads_file);
if ads_z0 ~= z_ref                       % renormalise the ADS S11 through Z
    Z = ads_z0 * (1 + SA) ./ (1 - SA);  SA = (Z - z_ref) ./ (Z + z_ref);
end
src(end+1) = struct('name', "ADS", 'f', fA, 'S', SA, 'col', co(2,:), 'ls', ':', 'file', ads_file);

s2pCol = co([3 4 6], :);
for k = 1:numel(s2p_files)
    s2 = sparameters(char(s2p_files(k)));
    src(end+1) = struct('name', s2p_names(k) + " s2p, port 2 " + termName.(termination), 'f', s2.Frequencies(:), ...
        'S', terminatePort2(s2, termination, z_ref), 'col', s2pCol(mod(k-1, 3)+1, :), 'ls', '--', 'file', s2p_files(k)); %#ok<SAGROW>
    fprintf('%s s2p reference %.5g ohm (renormalised to %g ohm through Z)\n', s2p_names(k), s2.Impedance, z_ref);
end

if snp_file ~= ""
    sn = sparameters(char(snp_file));
    s2snp = reduceSnp(sn, snp_ports, mode_load);    % "ADS s2p": mode ports terminated in mode_load
    src(end+1) = struct('name', "ADS s2p, port 2 " + termName.(termination), ...
        'f', s2snp.Frequencies(:), 'S', terminatePort2(s2snp, termination, z_ref), 'col', co(5,:), 'ls', '-.', 'file', snp_file);
    fprintf('s%dp reference %.5g ohm; ports %s terminated in %g ohm (ADS s2p), renormalised to %g ohm\n', ...
        size(sn.Parameters, 1), sn.Impedance, mat2str(setdiff(1:size(sn.Parameters, 1), snp_ports)), mode_load, z_ref);
end

for k = 1:numel(src)
    [~, n, e] = fileparts(src(k).file);
    fprintf('%-42s %s%s: %d points, %.4g to %.4g GHz\n', src(k).name, n, e, numel(src(k).f), src(k).f(1)/1e9, src(k).f(end)/1e9);
end

%% ---- plot range
if isempty(plot_frange_GHz), fr = [src(1).f(1) src(1).f(end)];
else,                        fr = plot_frange_GHz * 1e9;
end
allF = vertcat(src.f);  fr = [max(fr(1), min(allF)) min(fr(2), max(allF))];

%% ---- resonance (|S11| minimum) in the plot range
fprintf('\n|S11| minimum between %.4g and %.4g GHz\n', fr/1e9);
for k = 1:numel(src)
    m = src(k).f >= fr(1) & src(k).f <= fr(2);
    if ~any(m), fprintf('  %-42s no data in range\n', src(k).name); continue, end
    [v, i] = min(abs(src(k).S(m)));  fk = src(k).f(m);
    fprintf('  %-42s %.4f GHz, |S11| = %.4f (%.2f dB)\n', src(k).name, fk(i)/1e9, v, 20*log10(v));
end

%% ---- relative RMS error of each pair (second source = reference)
fprintf('\nrelative RMS error, sqrt(sum|Sa - Sb|^2 / sum|Sb|^2), over the overlap (within the plot range)\n');
fprintf('  %-80s %-17s  %9s\n', 'a / b', 'band (GHz)', 'rel. RMS');
for p = nchoosek(1:numel(src), 2)'
    a = src(p(1));  b = src(p(2));
    [g, A, B] = commonGrid(a.f, a.S, b.f, b.S, fr);
    if isempty(g), fprintf('  %-80s no overlap\n', a.name + " / " + b.name); continue, end
    fprintf('  %-80s %6.4g to %-7.4g  %8.3g %%\n', a.name + " / " + b.name, g(1)/1e9, g(end)/1e9, relRMS(A, B));
end

%% ---- figure 1: S11 magnitude and phase, ADS vs the CST simulations (legend: rel. RMS error re fig1_ref)
in1  = find(~startsWith([src.name], fig1_exclude));
iRef = find(startsWith([src.name], fig1_ref), 1);
if isempty(iRef), fprintf('fig1_ref "%s" matches no source: no RMS error in the legend\n', fig1_ref); end
fig = figure('Name', "S11, port 2 " + termination + ", ADS vs CST");   % default style: docked as a tab of the figure window
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');   % magnitude | phase
tOpt = {'Interpreter', 'latex', 'FontSize', 11};
title(tl, {sprintf('Rectangular patch, port 2 %s', termName.(termination)), ...
    sprintf('$S_{11}$ re %g $\\Omega$', z_ref)}, 'Interpreter', 'latex', 'FontSize', 12);

ax1 = nexttile(tl); hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
for k = in1
    m = src(k).f >= fr(1) & src(k).f <= fr(2);
    nm = src(k).name + rmsTag(src(k).f, src(k).S, src(max([iRef 1])).f, src(max([iRef 1])).S, fr, k, iRef);
    sty = lineStyle(src(k).col, src(k).ls, k == iRef);
    plot(ax1, src(k).f(m)/1e9, 20*log10(abs(src(k).S(m))), sty{:}, 'DisplayName', nm);
    plot(ax2, src(k).f(m)/1e9, rad2deg(angle(src(k).S(m))), sty{:}, 'DisplayName', nm);
end
ylabel(ax1, '$|S_{11}|$ (dB)', tOpt{:});  title(ax1, 'Magnitude', tOpt{:});
ylabel(ax2, '$\angle S_{11}$ (deg)', tOpt{:});  title(ax2, 'Phase', tOpt{:});
xlabel(ax1, 'Frequency (GHz)', tOpt{:});  xlabel(ax2, 'Frequency (GHz)', tOpt{:});
lg = legend(ax1, 'Location', 'southwest', 'Interpreter', 'none', 'FontSize', 8);
if ~isempty(iRef), title(lg, 'rel. RMS error re the reference (ref)', 'FontWeight', 'normal'); end
linkaxes([ax1 ax2], 'x');  xlim(ax1, fr/1e9);

%% ---- figure 2: two-port S-parameters of the CST simulations and the ADS s2p
% (ports 1 feed, 2 load; no termination applied), all referred to z_ref
two = struct('name', {}, 's', {}, 'col', {}, 'ls', {});
for k = 1:numel(s2p_files)
    two(end+1) = struct('name', s2p_names(k) + " s2p", 's', sparameters(sparameters(char(s2p_files(k))), z_ref), ...
        'col', s2pCol(mod(k-1, 3)+1, :), 'ls', '--'); %#ok<SAGROW>
end
if snp_file ~= ""
    two(end+1) = struct('name', string(sprintf('ADS s2p (s%dp, mode ports %g ohm)', size(sn.Parameters, 1), mode_load)), ...
        's', sparameters(s2snp, z_ref), 'col', co(5,:), 'ls', '-.');
end
prm   = {[1 1], [2 1], [2 2]};  prmNm = ["S_{11}" "S_{21}" "S_{22}"];
iRef2 = find(startsWith([two.name], diff_ref), 1);

fig2 = figure('Name', 'Two-port S-parameters');
tl2  = tiledlayout(fig2, 2, numel(prm), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, {'Rectangular patch, two-port (port 1 feed, port 2 load)', ...
    sprintf('re %g $\\Omega$; ADS s2p: spherical-mode ports terminated in %g $\\Omega$', z_ref, mode_load)}, ...
    'Interpreter', 'latex', 'FontSize', 12);
axM = gobjects(1, numel(prm)); axP = axM;
for c = 1:numel(prm), axM(c) = nexttile(tl2, c);              hold(axM(c), 'on'); grid(axM(c), 'on'); box(axM(c), 'on'); end
for c = 1:numel(prm), axP(c) = nexttile(tl2, numel(prm) + c); hold(axP(c), 'on'); grid(axP(c), 'on'); box(axP(c), 'on'); end

fprintf('\ntwo-port S-parameters: relative RMS error re %s (%.4g to %.4g GHz)\n', two(max([iRef2 1])).name, fr/1e9);
for c = 1:numel(prm)
    ij = prm{c};
    fR = two(max([iRef2 1])).s.Frequencies(:);  SR = squeeze(two(max([iRef2 1])).s.Parameters(ij(1), ij(2), :));
    for k = 1:numel(two)
        f = two(k).s.Frequencies(:);  S = squeeze(two(k).s.Parameters(ij(1), ij(2), :));
        m = f >= fr(1) & f <= fr(2);
        tag = rmsTag(f, S, fR, SR, fr, k, iRef2);
        sty = lineStyle(two(k).col, two(k).ls, k == iRef2);
        plot(axM(c), f(m)/1e9, 20*log10(abs(S(m))), sty{:}, 'DisplayName', two(k).name + tag);
        plot(axP(c), f(m)/1e9, rad2deg(angle(S(m))),  sty{:}, 'DisplayName', two(k).name + tag);
        if ~isempty(iRef2) && k ~= iRef2, fprintf('  %-6s %-42s %s\n', erase(prmNm(c), ["_" "{" "}"]), two(k).name, strtrim(tag)); end
    end
    title(axM(c), sprintf('$|%s|$', prmNm(c)), tOpt{:});        ylabel(axM(c), 'dB', tOpt{:});
    title(axP(c), sprintf('$\\angle %s$', prmNm(c)), tOpt{:});  ylabel(axP(c), 'deg', tOpt{:});
    xlabel(axP(c), 'Frequency (GHz)', tOpt{:});
    legend(axM(c), 'Location', 'southwest', 'Interpreter', 'none', 'FontSize', 8);   % each panel: RMS error of that parameter
    linkaxes([axM(c) axP(c)], 'x');  xlim(axM(c), fr/1e9);
end

%% ---- local functions

function s = lineStyle(col, ls, isRef)
% plot options: the RMS error reference is a solid black line 0.5 pt thicker than the others
lw = 1.4;
if isRef, s = {'Color', 'k', 'LineStyle', '-', 'LineWidth', lw + 0.5};
else,     s = {'Color', col, 'LineStyle', ls,  'LineWidth', lw};
end
end

function s2 = reduceSnp(sn, keep, ZL)
% Multiport sn (sparameters object) reduced to the ports in keep, all other
% ports terminated in ZL ohm:  S = Saa + Sab*G*(I - Sbb*G)^-1*Sba, with
% G = (ZL - Z0)/(ZL + Z0) the reflection of the load re the file's Z0.
S = sn.Parameters;  t = setdiff(1:size(S, 1), keep);
G = (ZL - sn.Impedance) / (ZL + sn.Impedance);
P = S(keep, keep, :);
if abs(G) > 1e-12
    for n = 1:size(S, 3)
        P(:, :, n) = P(:, :, n) + S(keep, t, n) * ((eye(numel(t))/G - S(t, t, n)) \ S(t, keep, n));
    end
end
s2 = sparameters(P, sn.Frequencies, sn.Impedance);
end

function e = relRMS(A, B)
% relative root-mean-square error of A re B, in percent
e = 100 * sqrt(sum(abs(A - B).^2) / sum(abs(B).^2));
end

function t = rmsTag(f, S, fR, SR, fr, k, iRef)
% legend suffix: " (ref)" for the reference, " (RMS x %)" re the reference
% otherwise, "" when there is no reference or no overlap
t = "";
if isempty(iRef), return, end
if k == iRef, t = " (ref)"; return, end
[g, A, B] = commonGrid(f, S, fR, SR, fr);
if ~isempty(g), t = string(sprintf(' (RMS %.3g %%)', relRMS(A, B))); end
end

function [g, A, B] = commonGrid(fa, a, fb, b, fr)
% a(fa) and b(fb) on the coarser of their two grids over the common range
% (limited to fr); the finer one is interpolated (complex, linear)
lo = max([fa(1) fb(1) fr(1)]);  hi = min([fa(end) fb(end) fr(2)]);
ma = fa >= lo & fa <= hi;  mb = fb >= lo & fb <= hi;
if nnz(ma) < 2 || nnz(mb) < 2, g = []; A = []; B = []; return, end
if nnz(ma) <= nnz(mb), g = fa(ma); A = a(ma); B = interp1(fb, b, g);
else,                  g = fb(mb); B = b(mb); A = interp1(fa, a, g);
end
end

function S11 = terminatePort2(s2, termination, z_ref)
% S11 at port 1 of the two-port s2 (sparameters object) with port 2 shorted
% ("SC", ZL = 0) or open ("OC", ZL = Inf), referred to z_ref. Going through Z
% removes the two-port's own reference impedance.
Z  = zparameters(s2).Parameters;
Z11 = squeeze(Z(1,1,:)); Z12 = squeeze(Z(1,2,:)); Z21 = squeeze(Z(2,1,:)); Z22 = squeeze(Z(2,2,:));
if termination == "SC", Zin = Z11 - Z12 .* Z21 ./ Z22;    % ZL = 0
else,                   Zin = Z11;                        % ZL = Inf
end
S11 = (Zin - z_ref) ./ (Zin + z_ref);
end

function [f, S] = readAdsS11(file)
% ADS data-display CSV export: rows "freq [unit],mag,phase(deg)"; header and
% blank lines are skipped.
txt = splitlines(fileread(file));
tok = regexp(txt, '^\s*([-+\d.eE]+)\s*(GHz|MHz|kHz|Hz)?\s*,\s*([-+\d.eE]+)\s*,\s*([-+\d.eE]+)\s*$', 'tokens', 'once');
tok = tok(~cellfun(@isempty, tok));
if isempty(tok), error('No "freq,mag,phase" rows found in %s', file); end
tok = vertcat(tok{:});
scale = struct('GHz', 1e9, 'MHz', 1e6, 'kHz', 1e3, 'Hz', 1);
u = string(tok(:, 2));  u(u == "") = "Hz";
f = str2double(tok(:, 1)) .* arrayfun(@(x) scale.(char(x)), u);
S = str2double(tok(:, 3)) .* exp(1j * deg2rad(str2double(tok(:, 4))));
end
