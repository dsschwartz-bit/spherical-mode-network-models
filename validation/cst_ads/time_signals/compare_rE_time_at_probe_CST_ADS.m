%COMPARE_RE_TIME_AT_PROBE_CST_ADS  Compare CST and ADS radiated far-field time signals at probe directions.
% Compare the radiated far-field time signals r*E_theta(t) and r*E_phi(t)
% from CST and ADS at one or more observation points (far-field probe
% directions theta, phi at distance r): time signals, difference, and
% spectra (magnitude and phase). One figure (2 x 2 tiled layout) per
% observation point and field component.
% (For far-field patterns over angle at fixed frequencies, see
% ../farfield/compare_cst_ads_farfields.m.) Output: figures and console tables.
%
% ADS : CSV export from a Data Display. Column names are read from the
%       header: time [s], Vtheta_prime, Vphi_prime (r*E in V) and the sweep
%       angles theta_deg / phi_deg [deg] (or theta / phi [rad]). Without a
%       phi column the sweep is taken at phi = ads_phi. Repeated time
%       points (ADS transient breakpoints) are removed
%       (readAdsFarfieldSweeps.m).
% CST : farfield time-signal ASCII export with blocks (Abs), (Phi),
%       (Theta) of E in V/m for each probe and excitation [n], header
%       "(Spherical) (theta phi r)". Only excitation cst_excitation is used.
%       r*E = E * r_m, with r taken from the probe header (in the CST
%       length unit, normally mm) and converted to metres
%       (readCstFarfieldProbes.m).
%
% Every direction (theta, phi) present in both files is compared (within
% angle_tol; phi must match on the z axis too, because the theta/phi unit
% vectors there still depend on phi); directions in only one file
% are listed and skipped. Signals are compared as given (no time shift).
% Errors are relative to CST, over the overlap of the two records (from the
% later start to the earlier end), with ||x||^2 = integral x(t)^2 dt:
%   relative MSE = ||ADS - CST||^2 / ||CST component||^2
%                = integral (ADS - CST)^2 dt / integral CST^2 dt
%   MSE/total    = ||ADS - CST||^2 / ||CST total field||^2, ||total||^2 = ||E_theta||^2 + ||E_phi||^2
% The last stays meaningful when a component is close to zero. The
% best-fit delay and amplitude ratio are printed for information.
% In the frequency domain, the magnitude error |20 log10(|X_ADS|/|X_CST|)|
% averaged over the -bw_level_dB band of the CST spectrum is reported (also
% for the port waves); the band edges are marked in the magnitude panel.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
cst_root = "C:\path\to\CST\";                    % folder containing the CST project folders
ads_dir  = "C:\path\to\ADS\workspace\data\";   % ADS workspace data folder (CSV exports)
% =========================

%% ---- antenna, excitation and folders ------------------------------------
antenna_name = "rectangularpatch";
cst_variant  = "symmetry_DAM";     % CST model / solver variant
excitation   = "filteredsine";      % excitation tag in the file names, e.g. "gaussiansine", "sinestep", "filteredsine"
cst_fc       = "";                  % centre frequency in the CST project folder name, e.g. "2.5GHz"
                                     % ("" = project folder without an excitation suffix)
ads_band     = "2.45GHz";         % frequency tag in the ADS file names (-20 dB band)
cst_band     = "2.45GHz";         % frequency tag in the CST file names
                                     % (e.g. gaussiansine "2.4GHz_2.6GHz"; sinestep "2.5GHz")
port         = 1;                    % excited port
cst_project  = antenna_name + "_" + cst_variant;
if strlength(cst_fc) > 0, cst_project = cst_project + "_" + excitation + "_" + cst_fc; end

cst_dir  = string(fullfile(cst_root, cst_project, "Export", "time_signals")) + filesep;

%% ---- files and settings ---------------------------------------------------
ads_file  = string(fullfile(ads_dir, antenna_name + "_rE_" + excitation + "_" + ads_band + ".csv"));
cst_file  = cst_dir + ["ff_Etheta_" + excitation + "_port" + port + "_" + cst_band + ".txt" ...
                       "ff_Ephi_"   + excitation + "_port" + port + "_" + cst_band + ".txt"];   % [1] Theta / Phi blocks
cst_file  = cst_file(isfile(cst_file));        % use whichever components were exported
if isempty(cst_file)                           % else any far-field export ff_E*_<band>.txt in cst_dir
    L = dir(cst_dir + "ff_E*_" + cst_band + ".txt");
    cst_file = string(fullfile({L.folder}, {L.name}));
    if isempty(cst_file), error('No CST far-field export ff_E*_%s.txt in %s', cst_band, cst_dir); end
    fprintf('CST far-field file(s): %s\n', strjoin(string({L.name}), ', '));
end
% CST ASCII exports hold at most 25 curves: split many probes over several
% exports and list them all, e.g. cst_file = ["...part1.txt" "...part2.txt"]

% Port incident/outgoing waves (optional): compared in an extra figure if
% BOTH files exist. CST signals (i1, o1,1, o2,1, ...) in W^(1/2) are
% converted to volts, V = value * sqrt(Z0) with Z0 from the file's
% #Parameters line (importCSTSignal.m); CST signals in V are used as-is.
% CST names map to ADS columns as i1 -> Vi_1, o1,1 -> Vo_1_1, o2,1 -> Vo_2_1.
% CST file: iowaves_<exc>_<band>.txt or iowaves_<exc>_port<N>_<band>.txt, or else the excitation export
% excitation_<exc>_port_<band>.txt (it also holds the o signals)
cst_io_file = cst_dir + ["iowaves_" + excitation + "_" + cst_band + ".txt", ...
                         "iowaves_" + excitation + "_port" + port + "_" + cst_band + ".txt", ...
                         "excitation_" + excitation + "_port_" + cst_band + ".txt", ...
                         "excitation_" + excitation + "_" + cst_band + ".txt"];
cst_io_file = cst_io_file(find(isfile(cst_io_file), 1));
if isempty(cst_io_file), cst_io_file = cst_dir + "iowaves_" + excitation + "_" + cst_band + ".txt"; end
ads_io_file = string(fullfile(ads_dir, antenna_name + ["_SC_iowaves_" "_iowaves_"] + excitation + "_" + ads_band + ".csv"));
if any(isfile(ads_io_file)), ads_io_file = ads_io_file(find(isfile(ads_io_file), 1)); else, ads_io_file = ads_io_file(1); end
components = ["Theta" "Phi"];   % CST block (Theta)/(Phi), ADS column V<comp>_prime
cst_excitation = 1;      % CST excitation [n] to use (ADS excites port 1 only)
ads_theta = 90;          % deg, theta of the ADS sweep if the file has no theta column
ads_phi   = 180;         % deg, phi of the ADS sweep if the file has no phi column
cst_length_unit = 1e-3;  % CST project length unit in metres (mm)
a_enclosing = [];        % m, radius a of the sphere enclosing the antenna; [] = parameter 'a' of the
                         % CST export's #Parameters line (CST length unit). The field can first reach
                         % a probe at r at t0 = (r - a)/c0 (marked in the time plots)
angle_tol = 0.1;         % deg; ADS and CST directions match within this
arrival_tol = 1e-4;      % for information only: the time |rE| first exceeds
                         % arrival_tol x peak is printed (compare with t0 = (r - a)/c0)
bw_level_dB = 20;        % band for the in-band magnitude error: frequencies where the CST
                         % spectrum is within bw_level_dB of its peak (first crossings either
                         % side of the peak); the average of |20 log10(|X_ADS| / |X_CST|)|
                         % over that band is reported
Threshold = 60;          % dB below peak that sets the plotted frequency range
                         % (the lower of the two signals' cutoffs is used, so
                         % neither simulator's numerical noise floor is shown)

nC = numel(components);

% enclosing-sphere radius from the CST export
if isempty(a_enclosing)
    hdr = readlines(cst_file(1));
    pl  = hdr(find(startsWith(strtrim(hdr), "#Parameters"), 1));
    tk  = regexp(pl, '[{;]\s*a\s*=\s*([-+0-9.eE]+)', 'tokens', 'once');
    if isempty(tk)
        a_enclosing = 0;
        warning('No parameter a in %s; t0 = r/c0 is marked.', cst_file(1));
    else
        a_enclosing = str2double(tk{1}) * cst_length_unit;
        fprintf('Enclosing-sphere radius a = %g m (CST parameter a)\n', a_enclosing);
    end
end

%% ---- read both files --------------------------------------------------------
ads = readAdsFarfieldSweeps(ads_file, components, ads_phi);
for k = 1:numel(ads)                                   % no theta column: use ads_theta
    if isnan(ads(k).theta), ads(k).theta = ads_theta; end
end
if any(~[ads.phiInFile]) || any(isnan(ads(1).theta))
    fprintf('ADS file has no angle column(s): taken as theta = %g deg, phi = %g deg\n', ads_theta, ads_phi);
end
cst = readCstFarfieldProbes(cst_file, components, cst_excitation, cst_length_unit);
for n = 1:numel(cst)
    has = components(~cellfun(@isempty, cst(n).e));
    fprintf('CST [%d]: theta = %g deg, phi = %g deg, r = %g m, E_%s\n', cst_excitation, ...
        cst(n).theta, cst(n).phi, cst(n).r, strjoin(lower(has), ', E_'));
end

% pair directions: same theta and phi (phi ignored on the z axis). An ADS
% sweep without a theta column (theta NaN) matches any theta, i.e. CST
% probes are matched on phi only.
angd   = @(a, b) abs(mod(a - b + 180, 360) - 180);
sameDir = @(p, q) (isnan(p.theta) || abs(p.theta - q.theta) < angle_tol) && ...
    angd(p.phi, q.phi) < angle_tol;   % phi matters on the z axis too (components rotate with phi)
if any(isnan([ads.theta]) & [ads.phiInFile])
    fprintf('ADS file has no theta column: CST probes are matched on phi only\n');
end
pairs = zeros(0, 2);                           % [cst index, ads index]
for j = 1:numel(cst)
    if isscalar(ads) && isnan(ads(1).theta) && ~ads(1).phiInFile   % ADS file is one direction
        if j == 1, pairs(end+1, :) = [1 1]; end %#ok<AGROW>
        continue
    end
    k = find(arrayfun(@(a) sameDir(a, cst(j)), ads), 1);
    if isempty(k)
        fprintf('CST (%g, %g) deg: no ADS sweep; skipped\n', cst(j).theta, cst(j).phi);
    else
        pairs(end+1, :) = [j k]; %#ok<AGROW>
    end
end
for k = 1:numel(ads)
    if (~isnan(ads(k).theta) || ads(k).phiInFile) && ~any(pairs(:, 2) == k)
        fprintf('ADS (%g, %g) deg: no CST probe; skipped\n', ads(k).theta, ads(k).phi);
    end
end
if isempty(pairs), error('No direction is present in both files.'); end

%% ---- compare each observation point -----------------------------------------
o = struct('components', components, 'arrival_tol', arrival_tol, ...
    'Threshold', Threshold, 'ads_file', ads_file, 'cst_file', cst_file, 'bw_level_dB', bw_level_dB, ...
    'a', a_enclosing);
res = struct('theta', {}, 'phi', {}, 'r', {}, 'comp', {}, 'peakADS', {}, 'peakCST', {}, ...
    'errL2', {}, 'errTot', {}, 'lag', {}, 'g', {}, 'errL2fit', {}, 'magErrdB', {}, ...
    'fLo', {}, 'fHi', {});
for p = 1:size(pairs, 1)
    j = pairs(p, 1);  k = pairs(p, 2);
    res = [res, compareLocation(ads(k).t, ads(k).e, cst(j), o)]; %#ok<AGROW>
end

%% ---- summary ------------------------------------------------------------------
if numel(res) > 1
    fprintf('\nsummary\n');
    fprintf(['  (errors relative to CST over the overlap of the records)\n' ...
        '  theta   phi     r    comp    peak ADS    peak CST   rel. MSE  MSE/total  best-fit scale  ' ...
        'rel. MSE after fit   mean in-band |X| error (-%g dB band)\n'], bw_level_dB);
    for j = 1:numel(res)
        fprintf('  %5g  %4g  %5g m  %-5s  %9.4g V  %9.4g V  %s  %6.2f %%      %.4f        %s       %5.2f dB (%.3g-%.3g GHz)\n', ...
            res(j).theta, res(j).phi, res(j).r, res(j).comp, res(j).peakADS, res(j).peakCST, ...
            pct(res(j).errL2^2), 100*res(j).errTot^2, res(j).g, pct(res(j).errL2fit^2), ...
            res(j).magErrdB, res(j).fLo/1e9, res(j).fHi/1e9);
    end
end

%% ---- port incident / outgoing waves ---------------------------------------
if isfile(cst_io_file) && isfile(ads_io_file)
    compareIOWaves(cst_io_file, ads_io_file, bw_level_dB, port);
else
    fprintf('\nport-wave comparison skipped (file not found: %s)\n', ...
        strjoin([cst_io_file(~isfile(cst_io_file)), ads_io_file(~isfile(ads_io_file))], ', '));
end

%% ---- local functions ------------------------------------------------------
function s = pct(x)
% percentage for the tables; n/a when the reference norm was zero
if isfinite(x), s = sprintf('%6.2f %%', 100*x); else, s = '   n/a  '; end
end

function compareIOWaves(cst_io_file, ads_io_file, lev, port)
% Compare CST and ADS port incident/outgoing wave voltages, one tile per
% CST signal that has a matching ADS column. Both records are cut to their
% overlap (later start to earlier end) for the error, spectra and plot.
[tc, vc, info] = importCSTSignal(cst_io_file, []);   % seconds, volts
[names, A] = readAdsCsv(ads_io_file, false);
iT = find(strcmpi(names, "time"), 1);
if isempty(iT), error('No "time" column in %s', ads_io_file); end

% pair each CST signal with its ADS column: i1 -> Vi_1, oN,M -> Vo_N_M
pairs = zeros(0, 2);
for k = 1:numel(info.names)
    an = "V" + regexprep(string(info.names{k}), '^([io])(\d)', '$1_$2');
    an = strrep(an, ",", "_");
    if ~isempty(regexp(info.names{k}, '^o\d+$', 'once'))   % o1 (one excitation) -> Vo_1_<port>
        an = an + "_" + port;
    end
    j = find(strcmpi(names, an), 1);
    if isempty(j)
        fprintf('  CST %s: no ADS column %s, skipped\n', info.names{k}, an);
    else
        pairs(end+1, :) = [k j]; %#ok<AGROW>
    end
end
if isempty(pairs), warning('No matching port-wave signals between the two files.'); return; end

[ta, iu] = unique(A(:, iT), 'last');            % drop repeated breakpoints
n  = size(pairs, 1);
[~, cn, ce] = fileparts(cst_io_file);  [~, an0, ae] = fileparts(ads_io_file);
esc = @(s) strrep(s, '_', '\_');
fig = figure('Color', 'w', 'NumberTitle', 'off', 'Name', 'Port waves');
tl  = tiledlayout(fig, ceil(n/2), 2, 'TileSpacing', 'compact');   % two columns
title(tl, {'Port incident ($V^+$) / outgoing ($V^-$) waves; CST in W$^{1/2}$ converted to V by $\times\sqrt{Z_0}$', ...
    sprintf('CST (%s%s) vs ADS (%s%s)', esc(cn), ce, esc(an0), ae)}, ...
    'Interpreter', 'latex', 'FontSize', 11);
fprintf('\nport waves (CST converted to V, Z0 = %g ohm):\n', info.Z0(1));
ax = gobjects(n, 1);
for p = 1:n
    k = pairs(p, 1);  t1 = tc{k}(:);  v1 = vc{k}(:);
    v2 = A(iu, pairs(p, 2));  t2 = ta;
    % both records cut to their overlap: later start to the earlier last time
    tBeg = max(t1(1), t2(1));  tEnd = min(t1(end), t2(end));
    m1 = t1 >= tBeg & t1 <= tEnd;  m2 = t2 >= tBeg & t2 <= tEnd;
    t1 = t1(m1);  v1 = v1(m1);  t2 = t2(m2);  v2 = v2(m2);
    % relative MSE on a common uniform grid over the overlap:
    % integral (ADS - CST)^2 dt / integral CST^2 dt
    dt   = min(median(diff(t1)), median(diff(t2))) / 4;
    tg   = (tBeg:dt:tEnd).';
    x1 = interp1(t1, v1, tg, 'pchip', 0);  x2 = interp1(t2, v2, tg, 'pchip', 0);
    mse = trapz(tg, (x2 - x1).^2) / trapz(tg, x1.^2);
    [eM, fl, fh] = bandMagError(t1, v1, t2, v2, lev);
    fprintf(['  %-6s / %-7s  peak CST %.4g V, ADS %.4g V, relative MSE %.3g %%, ' ...
        'mean in-band |X| error %.2f dB (-%g dB band %.3g-%.3g GHz), %.4g to %.4g ns\n'], ...
        info.names{k}, names(pairs(p, 2)), max(abs(v1)), max(abs(v2)), 100*mse, eM, lev, fl/1e9, fh/1e9, ...
        tBeg*1e9, tEnd*1e9);

    ax(p) = nexttile(tl);
    plot(t1*1e9, v1, 'LineWidth', 1.3); hold on
    plot(t2*1e9, v2, '--', 'LineWidth', 1.3);
    grid on; box on; ylabel('V'); xlabel('Time (ns)');
    xlim([tBeg, tEnd] * 1e9);                    % overlap of the two records
    % short title: i1 -> V^+_{1}, o1,1 -> V^-_{1,1}
    sym = regexprep(info.names{k}, '^i(.*)$', 'V^+_{$1}');
    sym = regexprep(sym, '^o(.*)$', 'V^-_{$1}');
    title(sprintf('$%s$, rel. MSE %.2g\\%%, mean in-band $|X|$ error %.2f dB', ...
        sym, 100*mse, eM), 'Interpreter', 'latex');
    legend('CST', 'ADS', 'Location', 'southeast');
end
linkaxes(ax, 'x');
end

function [names, A] = readAdsCsv(file, allowOld)
% Read an ADS Data Display CSV export: numeric rows into A, column names
% from the header ("name, then "," name for each further column). With
% allowOld, a 3-column file without readable names is taken as
% (time, Vtheta_prime, Vphi_prime).
L = strtrim(readlines(file));
L = L(strlength(L) > 0);
isnum = ~cellfun(@isempty, regexp(cellstr(L), '^"?\s*[-+0-9.]', 'once'));
rows = cellfun(@(s) sscanf(strrep(strrep(s, '"', ''), ',', ' '), '%f').', ...
    cellstr(L(isnum)), 'UniformOutput', false);
A = vertcat(rows{:});
tok = regexp(strjoin(L(~isnum), newline), '(?:^"|",")(\w+)', 'tokens');
names = string(cellfun(@(c) c{1}, tok, 'UniformOutput', false));
if numel(names) ~= size(A, 2)
    if allowOld && size(A, 2) == 3
        names = ["time" "Vtheta_prime" "Vphi_prime"];
    else
        error('Could not read %d column names from %s.', size(A, 2), file);
    end
end
end

function res = compareLocation(ta, Ea, cp, o)
% Compare ADS (time ta, r*E columns Ea per component, NaN = missing) with
% one CST probe cp (fields theta, phi, r, t{c}, e{c}) for every component
% present in both: print the error figures and draw one figure with a row
% of four panels per component.
c0    = 299792458;
comps = o.components;  nC = numel(comps);
% cut both records to their overlap (later start to the earlier last time):
% errors, spectra and plots all use the same time span
hasC = find(~cellfun(@isempty, cp.e));
tB = max([ta(1), cellfun(@(x) x(1), cp.t(hasC))]);
tE = min([ta(end), cellfun(@(x) x(end), cp.t(hasC))]);
m  = ta >= tB & ta <= tE;  ta = ta(m);  Ea = Ea(m, :);
for c = hasC
    m = cp.t{c} >= tB & cp.t{c} <= tE;  cp.t{c} = cp.t{c}(m);  cp.e{c} = cp.e{c}(m);
end
r_m   = cp.r;  loc = [cp.theta cp.phi];
t0    = (r_m - o.a) / c0;                     % earliest arrival from the enclosing sphere
have  = false(1, nC);
pkc   = zeros(1, nC);
for c = 1:nC
    have(c) = ~all(isnan(Ea(:, c))) && ~isempty(cp.e{c});
    if have(c), pkc(c) = max([max(abs(Ea(:, c))), max(abs(cp.e{c}))]); end
end
% a component that is numerically zero in both files (e.g. E_phi in a
% symmetry plane) is skipped: its relative errors are meaningless
zeroC = have & pkc < 1e-9 * max(pkc);
for c = find(zeroC)
    fprintf('theta = %g, phi = %g: E_%s is zero in both files (peak %.2g V), skipped\n', ...
        loc(1), loc(2), lower(comps(c)), pkc(c));
end
have = have & ~zeroC;
res = struct('theta', {}, 'phi', {}, 'r', {}, 'comp', {}, 'peakADS', {}, 'peakCST', {}, ...
    'errL2', {}, 'errTot', {}, 'lag', {}, 'g', {}, 'errL2fit', {}, 'magErrdB', {}, ...
    'fLo', {}, 'fHi', {});
if ~any(have)
    warning('theta = %g, phi = %g: no component present in both files; skipped.', loc(1), loc(2));
    return
end

% ---- time-domain comparison per component ----
S = repmat(struct('tg', [], 'xa', [], 'xc', [], 'd', [], 'nA', 0, 'nC', 0, 'tBeg', NaN, 'nD', 0, ...
    'tOn', NaN, 'tEnd', NaN, 'lag', 0, 'g', NaN, 'nFit', NaN, 'tc', [], 'ec', [], 'ea', []), 1, nC);
for c = find(have)
    tc = cp.t{c};  ec = cp.e{c};  ea = Ea(:, c);
    tBeg = max(ta(1), tc(1));  tEnd = min(ta(end), tc(end));   % overlap of the two records
    dt   = min(median(diff(ta)), median(diff(tc))) / 4;
    tg   = (tBeg:dt:tEnd).';
    xa   = interp1(ta, ea, tg, 'pchip', 0);
    xc   = interp1(tc, ec, tg, 'pchip', 0);
    d    = xa - xc;
    % error window: the overlap of the two records; trapezoid integration
    % in time (uniform grid)
    win = true(size(tg));
    L2  = @(x) sqrt(trapz(tg(win), x(win).^2));
    pk  = max(max(abs(xa)), max(abs(xc)));
    iOn = find(abs(xa) > o.arrival_tol*pk | abs(xc) > o.arrival_tol*pk, 1);
    % best-fit delay and scale (information only)
    [rr, lags] = xcorr(xc, xa, round(0.5e-9/dt));
    [~, im] = max(rr);  lag = lags(im) * dt;   % positive correlation (abs would allow a half-cycle, sign-flipped fit)
    xs  = interp1(tg - lag, xc, tg, 'linear', 0);
    g   = (xa' * xs) / max(xs' * xs, realmin);
    S(c) = struct('tg', tg, 'xa', xa, 'xc', xc, 'd', d, 'nA', L2(xa), 'nC', L2(xc), 'tBeg', tBeg, 'nD', L2(d), ...
        'tOn', tg(max([iOn 1])), 'tEnd', tEnd, 'lag', lag, 'g', g, 'nFit', L2(xa - g*xs), ...
        'tc', tc, 'ec', ec, 'ea', ea);
end
nTot = sqrt(sum([S(have).nC].^2));             % CST total-field norm

where = sprintf('theta = %g deg, phi = %g deg, r = %g m', loc(1), loc(2), r_m);
fprintf('\n==== %s ====\n', where);
fprintf('t0 = (r - a)/c0 = %.4g ns (a = %g m)\n', t0*1e9, o.a);
for c = find(have)
    s = S(c);  cn = lower(comps(c));
    fprintf('-- E_%s (error window %.4g to %.4g ns, overlap of the records)\n', cn, s.tBeg*1e9, s.tEnd*1e9);
    fprintf('  |rE| first exceeds %g x peak at %.4g ns\n', o.arrival_tol, s.tOn*1e9);
    fprintf('  peak |rE|        ADS %.4g V, CST %.4g V (ADS/CST = %.4f)\n', ...
        max(abs(s.xa)), max(abs(s.xc)), max(abs(s.xa))/max(abs(s.xc)));
    fprintf('  error            relative MSE ||ADS-CST||^2/||CST E_%s||^2 = %s, /||CST total field||^2 = %.3g %%\n', ...
        cn, strtrim(pct((s.nD/s.nC)^2)), 100*(s.nD/nTot)^2);
    fprintf('  best-fit delay   CST lags ADS by %.4g ps; best-fit ADS = %.4f x CST,\n', s.lag*1e12, s.g);
    fprintf('                   leaving %s relative MSE\n', strtrim(pct((s.nFit/s.nC)^2)));
end

% ---- spectra (non-uniform DFT, trapezoid weights) ----
% plotted frequency range set by the strongest component
X = @(t, v, f) nufft(v .* ([diff(t(1:2)); diff(t(1:end-1)) + diff(t(2:end)); diff(t(end-1:end))]/2), t, f);
nA = zeros(1, nC);  nA(have) = [S(have).nA];
[~, cD] = max(nA);
fAll = linspace(0, 1/(2*max(median(diff(ta)), median(diff(S(cD).tc)))), 4001).';
Xa0 = abs(X(ta, S(cD).ea, fAll));  Xc0 = abs(X(S(cD).tc, S(cD).ec, fAll));
last = @(m) fAll(find(m >= max(m) * 10^(-o.Threshold/20), 1, 'last'));
fEff = min(last(Xa0), last(Xc0));
fprintf('spectrum of E_%s within %g dB of peak up to %.3g GHz (ADS) / %.3g GHz (CST)\n', ...
    lower(comps(cD)), o.Threshold, last(Xa0)/1e9, last(Xc0)/1e9);
dec  = 10^floor(log10(fEff));
fMax = min(fAll(end), ceil(2 * fEff / dec) / 2 * dec);   % round up to 0.5 x decade
f    = linspace(0, fMax, 2001).';
% phase de-embedded by the propagation to the probe: multiply by exp(+jkr),
% k = 2*pi*f/c0, which removes the r/c0 delay (spectra use exp(-j2pi f t))
deemb = exp(1j * 2*pi*f/c0 * r_m);

% ---- figures: one 2 x 2 tiled layout per component ----
shown = find(have);  nR = numel(shown);
[~, cn0, ce] = fileparts(o.cst_file(1));
if numel(o.cst_file) > 1, ce = ce + " + " + (numel(o.cst_file) - 1) + " more"; end
[~, an, ae] = fileparts(o.ads_file);
esc = @(s) strrep(s, '_', '\_');              % file names in LaTeX text
dB  = @(z) 20*log10(abs(z) + realmin);
% common magnitude scale: from Threshold dB below the strongest component's
% peak (the same level that sets the frequency range) to 5 dB above it, so
% a near-zero component shows as near-zero rather than rescaling its panel
pkdB = max([dB(X(ta, S(cD).ea, f)); dB(X(S(cD).tc, S(cD).ec, f))]);
for ir = 1:nR
    c = shown(ir);  s = S(c);  sym = sprintf('rE_\\%s', lower(comps(c)));
    Xa = X(ta, s.ea, f);  Xc = X(s.tc, s.ec, f);
    [eM, fl, fh] = bandMagError(s.tc, s.ec, ta, s.ea, o.bw_level_dB);
    fprintf('  E_%s mean in-band |X| error (ADS vs CST) %.2f dB (-%g dB band of CST: %.3g-%.3g GHz)\n', ...
        lower(comps(c)), eM, o.bw_level_dB, fl/1e9, fh/1e9);

    % short tab label, e.g. "rEθ θ=30 φ=0" (Greek letters by code point)
    cg  = struct('theta', char(952), 'phi', char(966));
    fig = figure('Color', 'w', 'NumberTitle', 'off', 'Name', sprintf('rE%s %s=%g %s=%g', ...
        cg.(lower(char(comps(c)))), char(952), loc(1), char(966), loc(2)));
    tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
    title(tl, {sprintf('$%s$ at $\\theta = %g^\\circ$, $\\phi = %g^\\circ$, $r = %g$ m', ...
        sym, loc(1), loc(2), r_m), ...
        sprintf('CST (%s%s) vs ADS (%s%s)', esc(cn0), ce, esc(an), ae)}, ...
        'Interpreter', 'latex', 'FontSize', 11);
    axT = gobjects(1, 2);

    axT(1) = nexttile(tl);
    plot(s.tc*1e9, s.ec, 'LineWidth', 1.3); hold on
    plot(ta*1e9, s.ea, '--', 'LineWidth', 1.3);
    grid on; box on; xlabel('Time (ns)'); ylabel(sprintf('$%s$ (V)', sym), 'Interpreter', 'latex');
    xline(t0*1e9, ':', sprintf('$t_0 = (r-a)/c_0 = %.3g$ ns', t0*1e9), 'Interpreter', 'latex', ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    legend('CST', 'ADS', 'Location', 'best');
    title('Time signal');

    axT(2) = nexttile(tl);
    plot(s.tg*1e9, s.d, 'LineWidth', 1.1);
    grid on; box on; xlabel('Time (ns)'); ylabel('ADS - CST (V)');
    xline(t0*1e9, ':', '$t_0$', 'Interpreter', 'latex', 'LabelVerticalAlignment', 'bottom');
    if s.nC > 0
        eTxt = sprintf('rel. MSE %.2g\\%%', 100*(s.nD/s.nC)^2);
    else
        eTxt = sprintf('MSE %.2g\\%% of total', 100*(s.nD/nTot)^2);
    end
    title(sprintf('Difference: %s', eTxt), 'Interpreter', 'latex');
    linkaxes(axT, 'x');
    xlim(axT(1), [min(s.tc(1), ta(1)), max(s.tc(end), ta(end))] * 1e9);   % span of the time data

    nexttile(tl);
    plot(f/1e9, dB(Xc), 'LineWidth', 1.3); hold on
    plot(f/1e9, dB(Xa), '--', 'LineWidth', 1.3);
    grid on; box on; xlabel('Frequency (GHz)');
    ylabel(sprintf('$|%s(f)|$ (dB re 1 V/Hz)', sym), 'Interpreter', 'latex');
    xlim([f(1) f(end)]/1e9);                                % span of the frequency data
    ylim([pkdB - o.Threshold, pkdB + 5]);
    xline([fl fh]/1e9, ':', 'HandleVisibility', 'off');      % -bw_level_dB band
    legend('CST', 'ADS', 'Location', 'best');
    title(sprintf('Spectrum magnitude; mean in-band error %.2f dB', eM));

    nexttile(tl);
    plot(f/1e9, wrapBreak(angle(Xc.*deemb)*180/pi), 'LineWidth', 1.3); hold on
    plot(f/1e9, wrapBreak(angle(Xa.*deemb)*180/pi), '--', 'LineWidth', 1.3);
    grid on; box on; xlabel('Frequency (GHz)');
    ylabel(sprintf('$\\angle\\, %s(f)\\, e^{jkr}$ (deg)', sym), 'Interpreter', 'latex');
    ylim([-180 180]); yticks(-180:90:180); xlim([f(1) f(end)]/1e9);
    legend('CST', 'ADS', 'Location', 'best');
    title('Phase, de-embedded by $e^{-jkr}$', 'Interpreter', 'latex');

    res(end+1) = struct('theta', loc(1), 'phi', loc(2), 'r', r_m, 'comp', comps(c), ...
        'peakADS', max(abs(s.xa)), 'peakCST', max(abs(s.xc)), 'errL2', s.nD/s.nC, ...
        'errTot', s.nD/nTot, 'lag', s.lag, 'g', s.g, 'errL2fit', s.nFit/s.nC, ...
        'magErrdB', eM, 'fLo', fl, 'fHi', fh); %#ok<AGROW>
end
end

function [eAvg, fl, fh] = bandMagError(tc, vc, ta, va, lev)
% Magnitude error |20 log10(|X_ADS| / |X_CST|)| (dB) averaged over the band
% (uniform in frequency) where the CST spectrum is within lev dB of its peak: the first
% frequencies on each side of the peak where |X_CST| drops lev dB below it.
% Spectra by non-uniform DFT with trapezoid weights.
X  = @(t, v, f) nufft(v .* ([diff(t(1:2)); diff(t(1:end-1)) + diff(t(2:end)); diff(t(end-1:end))]/2), t, f);
fN = 1 / (2 * median(diff(tc)));
T  = tc(end) - tc(1);
fk = linspace(0, fN, max(4001, ceil(fN * 4 * T))).';
M  = abs(X(tc, vc, fk));
[pk, ip] = max(M);  thr = pk * 10^(-lev/20);
il = find(M(1:ip) < thr, 1, 'last');  if isempty(il), il = 1; end
iu = ip - 1 + find(M(ip:end) < thr, 1, 'first');  if isempty(iu), iu = numel(fk); end
f  = linspace(fk(il), fk(iu), 4001).';                 % refine
Mc = abs(X(tc, vc, f));  Ma = abs(X(ta, va, f));
[~, jp] = max(Mc);
jl = find(Mc(1:jp) < thr, 1, 'last');  if isempty(jl), jl = 0; end
ju = jp - 1 + find(Mc(jp:end) < thr, 1, 'first');  if isempty(ju), ju = numel(f) + 1; end
in = (jl + 1):(ju - 1);
e  = 20 * log10(Ma(in) ./ Mc(in));
eAvg = mean(abs(e));  fl = f(in(1));  fh = f(in(end));
end

function p = wrapBreak(p)
% Insert NaN at +-180 deg wraps so the plot line does not draw vertical
% jumps across the axes.
p(find(abs(diff(p)) > 180) + 1) = NaN;
end
