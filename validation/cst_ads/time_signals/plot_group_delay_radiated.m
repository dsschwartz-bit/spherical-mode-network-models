%PLOT_GROUP_DELAY_RADIATED  Group delay of the radiated far field (ADS AC, CST and ADS transient).
% Group delay of the radiated far field r*E_comp(f) at one or more far-field
% directions, after removing the propagation delay r/c0 (multiplying by
% exp(+jkr)), so what remains is the antenna's own delay relative to the
% coordinate origin:
%   tau_g(f) = -d/domega angle( H(f) * exp(+j*2*pi*f*r/c0) )
%
% Sources (each optional except the ADS AC export):
%   ADS AC     : far-field export, columns freq, mag/phase(V<comp>_prime)
%                [deg], and optionally theta_deg / phi_deg (one curve per
%                direction). Without angle columns the export is taken as
%                the single direction ads_ac_dir. H = r*E per source phasor.
%   CST trans. : far-field probe time signals divided by the port incident
%                wave i<port> in the frequency domain: H = rE(f) / V_inc(f).
%   ADS trans. : same, with the ADS rE CSV and incident wave Vi_<port>.
% The transient transfer functions are only meaningful where the
% excitation has energy: they are drawn only where |rE(f)| is within
% band_level_dB of its peak. They also need COMPLETE records: rE(t) and
% V_inc(t) must both have died away before the record ends (rE = H*V_inc
% only holds for the whole signals). A sine step whose simulation stops
% while the source is still on gives a wrong group delay; use a pulse
% such as the Gaussian sine. The script warns when a record ends above
% -40 dB of its peak.
%
% Phase: unwrapped along frequency before differentiating. Keep the AC
% frequency step small enough that the phase changes by well under 180 deg
% per step (the r/c0 term alone is 360*df*r/c0 deg per step: 12 deg at
% df = 10 MHz, r = 1 m).
% Output: figures and a console table.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ (read_ads_farfield_csv in matlab/io/cst) and third_party/

%% ===== USER SETTINGS =====
cst_root = "C:\path\to\CST\";                    % folder containing the CST project folders
ads_dir  = "C:\path\to\ADS\workspace\data\";   % ADS workspace data folder (CSV exports)
% =========================

%% ---- antenna, excitation and folders ------------------------------------
antenna_name = "rectangularpatch";
cst_variant  = "symmetry_FDTD";     % CST model / solver variant
excitation   = "gaussiansine";      % "gaussiansine" or "sinestep" (see note on complete records)
cst_fc       = "2.5GHz";            % centre frequency in the CST project folder name
ads_band     = "2.4GHz_2.6GHz";     % frequency tag in the ADS file names
cst_band     = "2.4GHz_2.6GHz";     % frequency tag in the CST file names
                                     % (gaussiansine: "2.4GHz_2.6GHz"; sinestep: "2.5GHz")
port         = 1;                    % excited port
cst_project  = antenna_name + "_" + cst_variant + "_" + excitation + "_" + cst_fc;

cst_dir  = string(fullfile(cst_root, cst_project, "Export", "time_signals")) + filesep;

%% ---- files ------------------------------------------------------------------
ads_ac_file  = string(fullfile(ads_dir, antenna_name + "_FF_export.csv"));
ads_rE_file  = string(fullfile(ads_dir, antenna_name + "_rE_" + excitation + "_" + ads_band + ".csv"));
ads_io_file  = string(fullfile(ads_dir, antenna_name + "_iowaves_" + excitation + "_" + ads_band + ".csv"));
cst_inc_file = cst_dir + ["iowaves_" + excitation + "_" + cst_band + ".txt", ...
                          "excitation_" + excitation + "_port_" + cst_band + ".txt"];
cst_inc_file = cst_inc_file(find(isfile(cst_inc_file), 1));

%% ---- settings ---------------------------------------------------------------
comp          = "Theta";       % "Theta" or "Phi"
cst_rE_file   = cst_dir + "ff_E" + lower(comp) + "_" + excitation + "_port" + port + "_" + cst_band + ".txt";   % CST export of comp
ads_ac_dir    = [90 0];        % deg (theta, phi) of the ADS AC export if it has no angle columns
plot_dirs     = [90 0; 90 180];% deg (theta, phi) rows: directions plotted (nearest on each grid)
vinc_source   = "auto";        % ADS transient incident wave: "ads", "cst" or "auto"
                               % (auto: CST i<port> if the ADS export looks corrupted)
r_m           = 1;             % m, observation distance of r*E (ADS AC and both transients)
cst_length_unit = 1e-3;        % CST project length unit in metres (mm)
use_transient = true;          % also show CST / ADS transient transfer functions
f_plot        = [1.5 3.5]*1e9; % Hz, plotted frequency range
f_report      = 2.5e9;         % Hz, group delay printed at this frequency
band_level_dB = 20;            % transient curves shown within this of their peak
c0 = 299792458;
angd  = @(a, b) abs(mod(a - b + 180, 360) - 180);          % angular distance, deg
wantD = @(th, ph) any(abs(plot_dirs(:, 1) - th) < 0.5 & ...
    (angd(plot_dirs(:, 2), ph) < 0.5 | th < 0.5 | th > 179.5));   % direction requested?

%% ---- ADS AC -------------------------------------------------------------------
F  = read_ads_farfield_csv(ads_ac_file);
mk = "mag" + "V" + lower(comp) + "prime";   pk = "phase" + "V" + lower(comp) + "prime";
if ~isfield(F, mk), error('ADS AC file has no mag(V%s_prime) column', lower(comp)); end
if isfield(F, 'thetadeg') && isfield(F, 'phideg')
    allD = unique([F.thetadeg(:) F.phideg(:)], 'rows');
    dirs = zeros(0, 2);
    for q = 1:size(plot_dirs, 1)                 % nearest grid direction to each request
        [~, j] = min(hypot(allD(:, 1) - plot_dirs(q, 1), angd(allD(:, 2), plot_dirs(q, 2))));
        dirs(end+1, :) = allD(j, :); %#ok<AGROW>
    end
    dirs = unique(dirs, 'rows', 'stable');
else
    dirs = ads_ac_dir;  F.thetadeg = repmat(dirs(1), size(F.freq));  F.phideg = repmat(dirs(2), size(F.freq));
end
curves = struct('src', {}, 'theta', {}, 'phi', {}, 'f', {}, 'H', {}, 'show', {});
for d = 1:size(dirs, 1)
    k = find(F.thetadeg == dirs(d, 1) & F.phideg == dirs(d, 2));
    [f, o] = sort(F.freq(k));  k = k(o);
    H = F.(mk)(k) .* exp(1j * F.(pk)(k) * pi/180);
    curves(end+1) = struct('src', "ADS AC", 'theta', dirs(d, 1), 'phi', dirs(d, 2), ...
        'f', f(:), 'H', H(:), 'show', true(numel(f), 1)); %#ok<SAGROW>
end
fprintf('ADS AC: %d direction(s), %.4g-%.4g GHz, df = %.3g MHz\n', size(dirs, 1), ...
    curves(1).f(1)/1e9, curves(1).f(end)/1e9, median(diff(curves(1).f))/1e6);

%% ---- transient transfer functions rE(f) / V_inc(f) -------------------------
if use_transient
    fT = linspace(f_plot(1), f_plot(2), 801).';
    % CST
    if isfile(cst_rE_file) && ~isempty(cst_inc_file)
        [tc, vc, info] = importCSTSignal(cst_inc_file, []);
        ki = find(strcmp(info.names, "i" + port), 1);
        tCi = tc{ki}(:);  vCi = vc{ki}(:);
        checkComplete("CST i" + port, vCi);
        Vi = spec(tCi, vCi, fT);
        cst = readCstFarfieldProbes(cst_rE_file, comp, 1, cst_length_unit);
        for n = 1:numel(cst)
            if isempty(cst(n).e{1}) || ~wantD(cst(n).theta, cst(n).phi), continue, end
            checkComplete(sprintf('CST rE (%g,%g)', cst(n).theta, cst(n).phi), cst(n).e{1});
            curves(end+1) = transferCurve("CST trans.", cst(n).theta, cst(n).phi, fT, ...
                spec(cst(n).t{1}(:), cst(n).e{1}(:), fT), Vi, band_level_dB); %#ok<SAGROW>
        end
    else
        fprintf('CST transient skipped (file not found)\n');
    end
    % ADS
    if isfile(ads_rE_file) && isfile(ads_io_file)
        [names, A] = readAdsTable(ads_io_file);
        iT = find(strcmpi(names, "time"), 1);  iV = find(strcmpi(names, "Vi_" + port), 1);
        iP = find(strcmpi(names, "phi_deg"), 1);
        if ~isempty(iP), A = A(A(:, iP) == A(1, iP), :); end       % one sweep is enough
        [ta, iu] = unique(A(:, iT), 'last');  va = A(iu, iV);
        useCst = vinc_source == "cst";
        if vinc_source == "auto" && exist('vCi', 'var') && max(abs(va)) > 1.5*max(abs(vCi))
            fprintf(['ADS Vi_%d peaks at %.4g V, CST i%d at %.4g V: the ADS export looks ' ...
                'corrupted; using CST i%d as the ADS incident wave\n'], port, max(abs(va)), ...
                port, max(abs(vCi)), port);
            useCst = true;
        end
        if useCst, Vi = spec(tCi, vCi, fT); else, Vi = spec(ta, va, fT); checkComplete("ADS Vi_" + port, va); end
        ads = readAdsFarfieldSweeps(ads_rE_file, comp, ads_ac_dir(2));
        for n = 1:numel(ads)
            th = ads(n).theta;  if isnan(th), th = ads_ac_dir(1); end   % no theta column
            if ~wantD(th, ads(n).phi), continue, end
            checkComplete(sprintf('ADS rE (%g,%g)', th, ads(n).phi), ads(n).e(:, 1));
            curves(end+1) = transferCurve("ADS trans.", th, ads(n).phi, fT, ...
                spec(ads(n).t, ads(n).e(:, 1), fT), Vi, band_level_dB); %#ok<SAGROW>
        end
    else
        fprintf('ADS transient skipped (file not found)\n');
    end
end

%% ---- de-embed r/c0 and differentiate ------------------------------------------
for n = 1:numel(curves)
    f = curves(n).f;
    curves(n).Hd  = curves(n).H .* exp(1j*2*pi*f*r_m/c0);      % remove exp(-jkr)
    curves(n).tau = groupDelay(f, curves(n).Hd);
end

%% ---- figure ---------------------------------------------------------------------
dirKey = unique(string(compose('%g,%g', [[curves.theta]; [curves.phi]].')), 'stable');
cols   = lines(max(numel(dirKey), 1));
sty    = containers.Map(["ADS AC" "CST trans." "ADS trans."], {'-', '--', ':'});
esc    = @(s) strrep(s, '_', '\_');
fig = figure('Color', 'w', 'NumberTitle', 'off', 'Name', 'Group delay');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');   % magnitude | phase; group delay below, full width
title(tl, {sprintf('Radiated $rE_\\%s$: group delay with $r/c_0 = %.4g$ ns removed ($r = %g$ m)', ...
    lower(comp), r_m/c0*1e9, r_m), sprintf('%s; %s excitation', esc(cst_project), esc(excitation))}, ...
    'Interpreter', 'latex', 'FontSize', 11);
ax = gobjects(3, 1);
ax(1) = nexttile(tl, 1);  ax(2) = nexttile(tl, 2);  ax(3) = nexttile(tl, 3, [1 2]);
for p = 1:3, hold(ax(p), 'on'); grid(ax(p), 'on'); box(ax(p), 'on'); end
lg = strings(0);  y1min = 0;
for n = 1:numel(curves)
    c = curves(n);  m = c.show & c.f >= f_plot(1) & c.f <= f_plot(2);
    col = cols(dirKey == string(sprintf('%g,%g', c.theta, c.phi)), :);
    fG = c.f/1e9;  y1 = 20*log10(abs(c.H)/max(abs(c.H(m))));  y2 = angle(c.Hd)*180/pi;  y3 = c.tau*1e9;
    y1(~m) = NaN;  y2(~m) = NaN;  y3(~m) = NaN;
    y1min = min([y1min; y1(isfinite(y1))]);
    y2(find(abs(diff(y2)) > 180) + 1) = NaN;                   % no vertical wrap lines
    plot(ax(1), fG, y1, sty(c.src), 'Color', col, 'LineWidth', 1.3);
    plot(ax(2), fG, y2, sty(c.src), 'Color', col, 'LineWidth', 1.3);
    plot(ax(3), fG, y3, sty(c.src), 'Color', col, 'LineWidth', 1.3);
    lg(end+1) = sprintf('%s  \\theta=%g\\circ \\phi=%g\\circ', c.src, c.theta, c.phi); %#ok<SAGROW>
end
ylabel(ax(1), '|H| (dB re max)');
ylim(ax(1), [max(floor(y1min/5)*5, -60), 2]);            % fit the plotted data (floor -60 dB)
title(ax(1), 'Magnitude (AC: $rE$; transient: $rE/V_{inc}$), each normalised to its maximum', 'Interpreter', 'latex');
ylabel(ax(2), 'deg');  ylim(ax(2), [-180 180]);  yticks(ax(2), -180:90:180);
title(ax(2), 'Phase with propagation removed, $\angle\, H\, e^{jkr}$', 'Interpreter', 'latex');
ylabel(ax(3), '\tau_g (ns)');  xlabel(ax(3), 'Frequency (GHz)');
yline(ax(3), 0, 'k:', 'HandleVisibility', 'off');
title(ax(3), 'Group delay relative to the origin, $\tau_g = -d\angle(H e^{jkr})/d\omega$', 'Interpreter', 'latex');
linkaxes(ax, 'x');  xlim(ax(1), f_plot/1e9);
lgd = legend(ax(1), lg, 'NumColumns', 3, 'Interpreter', 'tex');
lgd.Layout.Tile = 'south';                                % at the bottom, below the whole layout
xlabel(ax(1), 'Frequency (GHz)');  xlabel(ax(2), 'Frequency (GHz)');

%% ---- report at f_report ----------------------------------------------------------
fprintf('\ngroup delay at %.4g GHz (r/c0 = %.4f ns removed):\n', f_report/1e9, r_m/c0*1e9);
for n = 1:numel(curves)
    c = curves(n);  [~, k] = min(abs(c.f - f_report));
    if c.show(k), s = sprintf('%8.1f ps', c.tau(k)*1e12); else, s = '   (outside excitation band)'; end
    fprintf('  %-10s theta %3g phi %3g : %s\n', c.src, c.theta, c.phi, s);
end

%% ========================================================================
function X = spec(t, v, f)
% spectrum by non-uniform DFT, trapezoid weights
t = t(:);  v = v(:);  dt = diff(t);
w = [dt(1); dt(1:end-1) + dt(2:end); dt(end)] / 2;
X = nufft(v .* w, t, f(:));
end

function checkComplete(name, v)
% warn when a record has not died away by its end (rE = H*V_inc needs
% complete signals)
v = abs(v(:));  nEnd = max(5, round(0.02*numel(v)));
lev = 20*log10(max(v(end-nEnd+1:end)) / max(v) + realmin);
if lev > -40
    warning('%s ends at %.0f dB re its peak: record not complete, transient group delay unreliable', name, lev);
end
end

function c = transferCurve(src, th, ph, f, E, Vi, levelDB)
% transfer function rE/V_inc, flagged where the radiated spectrum is
% within levelDB of its peak (inside the excitation band)
show = abs(E) >= max(abs(E)) * 10^(-levelDB/20);
c = struct('src', src, 'theta', th, 'phi', ph, 'f', f(:), 'H', E(:) ./ Vi(:), 'show', show(:));
end

function tau = groupDelay(f, H)
% -d(phase)/d(omega), central differences on the unwrapped phase
tau = -gradient(unwrap(angle(H(:))), 2*pi*f(:));
end

function [names, A] = readAdsTable(file)
% numeric rows into A, column names from the ADS CSV header
L = strtrim(readlines(file));  L = L(strlength(L) > 0);
isnum = ~cellfun(@isempty, regexp(cellstr(L), '^"?\s*[-+0-9.]', 'once'));
rows = cellfun(@(s) sscanf(strrep(strrep(s, '"', ''), ',', ' '), '%f').', ...
    cellstr(L(isnum)), 'UniformOutput', false);
A = vertcat(rows{:});
tok = regexp(strjoin(L(~isnum), newline), '(?:^"|",")(\w+)', 'tokens');
names = string(cellfun(@(c) c{1}, tok, 'UniformOutput', false));
if numel(names) ~= size(A, 2), error('Could not read %d column names from %s.', size(A, 2), file); end
end
