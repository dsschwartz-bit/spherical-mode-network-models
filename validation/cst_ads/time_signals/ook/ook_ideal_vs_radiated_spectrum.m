%OOK_IDEAL_VS_RADIATED_SPECTRUM  Ideal OOK spectrum vs the radiated r*E spectrum from CST and ADS.
% On-off keying by direct antenna modulation: compare the ideal OOK
% spectrum with the realized radiated spectrum r*E(f) at every observation
% point (far-field probe direction) of the CST and ADS far-field exports
% (the same files as ../compare_rE_time_at_probe_CST_ADS.m).
%
% Ideal OOK signal: the carrier is switched off instantly when the control
% (switch) voltage at the external port ctrl_port first passes
% switch_fraction (1/2) of its final (high) value:
%     x_ideal(t) = carrier(t) * [t < t_sw]
%   carrier = "incident" : the port rf_port incident wave from the CST
%                          i/o-wave export (the excitation actually applied)
%           = "sine"     : sin(2 pi f_carrier t) from t = 0
% The radiated field of an ideal OOK transmitter is this signal delayed and
% scaled, so at each observation point and field component the ideal signal
% is fitted to the realized r*E: one delay tau for all observation points
% (same r), found on the strongest CST signal (largest energy over all points
% and components) as the delay, from (r - a)/c0 to (r - a)/c0 + max_delay,
% whose scaled ideal signal fits it best over the ON interval; tau is used for
% every component and both simulators, and
% a least-squares scale g per component and simulator over the ON interval
% t < t_sw + tau only (for the time plot and the ON/OFF figures); the spectra
% are compared each normalized to its own peak.
% Plots: one figure per observation point, docked as tabs of the MATLAB figure
% window; one row per field component: time signal | normalized spectrum. The
% time panels of E_theta and E_phi share the same y limits (same scale).
%
% Reported per observation point, component and simulator:
%   fit       : delay tau, scale g, relative RMS difference over the ON interval
%   OFF energy: energy of r*E after t_sw + tau, % of its total energy (ideal: 0)
%   bandwidths: -bw_level_dB bandwidth and the occ_pct % occupied bandwidth of
%               the ideal and the realized spectra
%   spectrum  : mean difference of the peak-normalized spectra (dB) over the
%               ideal's -bw_level_dB band
% All records are cut to the common span: t = 0 to the earliest last time of
% the CST, ADS, carrier and switch records.
%
% Helpers: importCSTSignal.m, readCstFarfieldProbes.m, readAdsFarfieldSweeps.m
% (matlab/io). Output: figures and a console table.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
cst_root = "C:\path\to\CST\";                    % folder containing the CST project folders
ads_dir  = "C:\path\to\ADS\workspace\data\";   % ADS workspace data folder (CSV exports)
% =========================

%% ---- antenna, excitation and folders ------------------------------------
antenna_name = "rectangularpatch";
cst_variant  = "symmetry_DAM";
excitation   = "filteredsine";
band         = "2.45GHz";        % frequency tag in the CST and ADS file names
rf_port      = 1;
ctrl_port    = 2;
cst_dir = string(fullfile(cst_root, antenna_name + "_" + cst_variant, "Export", "time_signals")) + filesep;

cst_ff_file = cst_dir + ["ff_Etheta_" "ff_Ephi_"] + excitation + "_port" + rf_port + "_" + band + ".txt";
ads_ff_file = string(fullfile(ads_dir, antenna_name + "_rE_" + excitation + "_" + band + ".csv"));   % "" = CST only
cst_io_file = cst_dir + "iowaves_" + excitation + "_port" + rf_port + "_" + band + ".txt";   % carrier
cst_v_file  = cst_dir + "excitation_" + excitation + "_voltage_" + band + ".txt";         % switch voltage

%% ---- settings -------------------------------------------------------------
carrier         = "incident";   % "incident" (port incident wave) or "sine"
f_carrier       = 2.45e9;       % Hz, for carrier = "sine"
switch_fraction = 0.5;          % carrier off when the switch voltage first passes this x final value
components      = ["Theta" "Phi"];
cst_excitation  = 1;
cst_length_unit = 1e-3;         % m
ads_phi         = 180;          % deg, ADS sweep phi if the file has no phi column
angle_tol       = 0.1;          % deg
dt              = 1e-12;        % s, common uniform grid (pchip interpolation)
nfft_min        = 2^18;         % zero padding: df = 1/(nfft*dt)
bw_level_dB     = 20;           % bandwidth level below the peak
occ_pct         = 99;           % occupied bandwidth (% of the energy)
plot_range_dB   = 60;           % spectrum plots: magnitude range below the peak
n_sidelobes     = 3;            % spectrum plots: frequency range = ideal main lobe + this many
                                %   sidelobes on each side (to the null after the last one)
max_delay       = 2e-9;         % s, delay search window beyond (r - a)/c0 (antenna group delay)
c0 = 299792458;

%% ---- carrier and switch ---------------------------------------------------
[tV, vV, infoV] = importCSTSignal(cst_v_file, []);
k2 = find(infoV.port == ctrl_port, 1);
if isempty(k2), error('No port %d voltage in %s', ctrl_port, cst_v_file); end
tS = tV{k2};  vS = vV{k2};
vFinal = max(vS);                                   % high level of the switch voltage
iS = find(vS >= switch_fraction * vFinal, 1);
if isempty(iS) || iS == 1, error('Switch voltage never passes %g x its final value', switch_fraction); end
t_sw = interp1(vS(iS-1:iS), tS(iS-1:iS), switch_fraction * vFinal);   % linear between samples
fprintf('Switch: port %d voltage passes %g x %.4g V at t_sw = %.5g ns\n', ctrl_port, switch_fraction, vFinal, t_sw*1e9);

[tW, vW, infoW] = importCSTSignal(cst_io_file, []);
ki = find(strcmp(infoW.names, "i" + rf_port), 1);
if isempty(ki), error('No i%d in %s', rf_port, cst_io_file); end
tCar = tW{ki};  vCar = vW{ki};

%% ---- far-field data -------------------------------------------------------
cst_ff_file = cst_ff_file(isfile(cst_ff_file));
if isempty(cst_ff_file), error('No CST far-field export in %s', cst_dir); end
cst = readCstFarfieldProbes(cst_ff_file, components, cst_excitation, cst_length_unit);
useADS = strlength(ads_ff_file) > 0 && isfile(ads_ff_file);
if useADS
    ads = readAdsFarfieldSweeps(ads_ff_file, components, ads_phi);
else
    ads = struct('theta', {}, 'phi', {}, 't', {}, 'e', {});
    fprintf('No ADS file: CST only\n');
end
hdr = readlines(cst_ff_file(1));
tk  = regexp(hdr(find(startsWith(strtrim(hdr), "#Parameters"), 1)), '[{;]\s*a\s*=\s*([-+0-9.eE]+)', 'tokens', 'once');
if isempty(tk), a_enc = 0; else, a_enc = str2double(tk{1}) * cst_length_unit; end

% common span: 0 to the earliest last time of every record used
tEnd = min([tCar(end), tS(end), cellfun(@(x) x(end), [cst.t])]);
if useADS, tEnd = min(tEnd, min(arrayfun(@(s) s.t(end), ads))); end
tg = (0:dt:tEnd).';
fprintf('Common span 0 to %.4g ns, dt = %g ps\n', tEnd*1e9, dt*1e12);

switch carrier
    case "incident", car = interp1(tCar, vCar, tg, 'pchip', 0);
    case "sine",     car = sin(2*pi*f_carrier*tg);
    otherwise, error('carrier must be "incident" or "sine"');
end
gate   = double(tg < t_sw);
xIdeal = car .* gate;

% spectra on a zero-padded FFT grid (one-sided, V/Hz)
nfft = max(nfft_min, 2^nextpow2(numel(tg)));
f    = (0:nfft/2).' / (nfft*dt);
spec = @(x) oneSided(dt * fft(x, nfft));
dB   = @(z) 20*log10(abs(z) + realmin);

Xid0 = spec(xIdeal);
[bwI, fLoI, fHiI] = levelBW(f, abs(Xid0), bw_level_dB);
[obI, oLoI, oHiI] = occBW(f, abs(Xid0), occ_pct);
fprintf('Ideal OOK (%s carrier): -%g dB band %.4g-%.4g GHz (%.4g MHz), %g %% occupied %.4g-%.4g GHz (%.4g MHz)\n', ...
    carrier, bw_level_dB, fLoI/1e9, fHiI/1e9, bwI/1e6, occ_pct, oLoI/1e9, oHiI/1e9, obI/1e6);
% plotted frequency range: the ideal spectrum's main lobe plus n_sidelobes
% sidelobes either side, i.e. out to the (n_sidelobes + 1)-th null from the peak
Mid = abs(Xid0);  [~, ip] = max(Mid);
null = find(islocalmin(Mid));
nLo = flip(null(null < ip));  nHi = null(null > ip);
if numel(nLo) > n_sidelobes, fPlot(1) = f(nLo(n_sidelobes + 1)); else, fPlot(1) = 0; end
if numel(nHi) > n_sidelobes, fPlot(2) = f(nHi(n_sidelobes + 1)); else, fPlot(2) = f(end); end
fprintf('Plotted range: main lobe + %d sidelobes, %.4g-%.4g GHz\n', n_sidelobes, fPlot/1e9);
inBand = f >= fLoI & f <= fHiI;

%% ---- one delay for every observation point ---------------------------------
% strongest CST signal over all probes and components
best = [0 0 0];
for j = 1:numel(cst)
    for c = 1:numel(components)
        if isempty(cst(j).e{c}), continue, end
        E = trapz(cst(j).t{c}, cst(j).e{c}.^2);
        if E > best(1), best = [E j c]; end
    end
end
pD = cst(best(2));  cD = best(3);
[tr, iu] = unique(pD.t{cD}(:), 'last');  er = pD.e{cD}(:);
xD   = interp1(tr, er(iu), tg, 'pchip', 0);
tMin = (pD.r - a_enc) / c0;
lagMin = max(0, floor(tMin/dt));  lagMax = min(numel(tg) - 1, ceil((tMin + max_delay)/dt));
% tau: the delay whose least-squares fit (scale free) leaves the smallest
% relative residual over the ON interval t < t_sw + tau
res = inf(1, lagMax + 1);
for L = lagMin:lagMax
    xs = [zeros(L, 1); xIdeal(1:end-L)];  on = tg < t_sw + L*dt;
    g  = (xD(on)' * xs(on)) / max(xs(on)' * xs(on), realmin);
    res(L+1) = norm(xD(on) - g*xs(on)) / norm(xD(on));
end
[~, im] = min(res);  tau = (im - 1) * dt;
fprintf(['Delay tau = %.4g ns for every point (best ON fit to CST E_%s at theta = %g, phi = %g deg;\n' ...
    '  (r - a)/c0 = %.4g ns)\n'], tau*1e9, lower(components(cD)), pD.theta, pD.phi, tMin*1e9);
xs = interp1(tg + tau, xIdeal, tg, 'linear', 0);          % ideal delayed by tau
on = tg < t_sw + tau;                                    % ON interval of the radiated field

%% ---- each observation point ------------------------------------------------
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
cg   = struct('theta', char(952), 'phi', char(966));
co   = orderedcolors("gem");
R = struct('theta', {}, 'phi', {}, 'comp', {}, 'src', {}, 'tau', {}, 'g', {}, 'rmsOn', {}, ...
    'offPct', {}, 'bw', {}, 'ob', {}, 'magErr', {});
for j = 1:numel(cst)
    p = cst(j);
    src = {"CST", p.t, p.e};
    if useADS
        k = find(arrayfun(@(a) abs(a.theta - p.theta) < angle_tol && angd(a.phi, p.phi) < angle_tol, ads), 1);
        if isempty(k)
            fprintf('(%g, %g) deg: no ADS sweep, CST only\n', p.theta, p.phi);
        else
            src(end+1, :) = {"ADS", repmat({ads(k).t}, 1, numel(components)), ...
                num2cell(ads(k).e, 1)}; %#ok<SAGROW>
        end
    end
    % components present and not zero in every source
    pkc = zeros(size(src, 1), numel(components));
    for s = 1:size(src, 1)
        for c = 1:numel(components)
            e = src{s, 3}{c};
            if ~isempty(e) && ~all(isnan(e)), pkc(s, c) = max(abs(e)); else, pkc(s, c) = NaN; end
        end
    end
    show = find(all(~isnan(pkc), 1) & max(pkc, [], 1) > 1e-9 * max(pkc(:)));
    for c = setdiff(find(all(~isnan(pkc), 1)), show)
        fprintf('(%g, %g) deg: E_%s is zero, skipped\n', p.theta, p.phi, lower(components(c)));
    end
    if isempty(show), continue, end

    fprintf('\n==== theta = %g deg, phi = %g deg, r = %g m ====\n', p.theta, p.phi, p.r);
    % realized signals on the common grid
    xs_all = cell(size(src, 1), numel(components));
    for s = 1:size(src, 1)
        for c = show
            tr = src{s, 2}{c}(:);  er = src{s, 3}{c}(:);
            [tr, iu] = unique(tr, 'last');
            xs_all{s, c} = interp1(tr, er(iu), tg, 'pchip', 0);
        end
    end
    fig = figure('Name', sprintf('OOK %s=%g %s=%g', char(952), p.theta, char(966), p.phi), ...
        'WindowStyle', 'docked');
    tl  = tiledlayout(fig, numel(show), 2, 'TileSpacing', 'compact');
    title(tl, sprintf('Ideal OOK vs radiated $rE$ at $\\theta = %g^\\circ$, $\\phi = %g^\\circ$, $r = %g$ m; switch at $t_{sw} = %.4g$ ns', ...
        p.theta, p.phi, p.r, t_sw*1e9), 'Interpreter', 'latex');
    axTall = gobjects(0);  yMax = 0;                % time panels share one y scale
    for c = show
        cn = lower(components(c));  sym = "rE_\" + cn;
        axT = nexttile(tl); hold on; box on; grid on;  axTall(end+1) = axT; %#ok<SAGROW>
        axN = nexttile(tl); hold on; box on; grid on;
        for s = 1:size(src, 1)
            x  = xs_all{s, c};
            g  = (x(on)' * xs(on)) / max(xs(on)' * xs(on), realmin);  % least squares over ON
            xi = g * xs;
            rmsOn  = 100 * norm(x(on) - xi(on)) / norm(x(on));
            offPct = 100 * sum(x(~on).^2) / sum(x.^2);
            X  = spec(x);  Xi = g * Xid0;                             % |X| independent of delay
            [bw, ~, ~] = levelBW(f, abs(X), bw_level_dB);
            [ob, ~, ~] = occBW(f, abs(X), occ_pct);
            nX = dB(X) - max(dB(X));  nXi = dB(Xid0) - max(dB(Xid0));     % peak-normalized
            mErr = mean(abs(nX(inBand) - nXi(inBand)));
            fprintf(['  E_%-5s %s: scale %.4g, rel. RMS difference (ON) %.2f %%, OFF energy %.3g %%; ' ...
                '-%g dB BW %.4g MHz (ideal %.4g), %g %% occupied BW %.4g MHz (ideal %.4g), in-band |X| error %.2f dB\n'], ...
                cn, src{s, 1}, g, rmsOn, offPct, bw_level_dB, bw/1e6, bwI/1e6, occ_pct, ob/1e6, obI/1e6, mErr);
            R(end+1) = struct('theta', p.theta, 'phi', p.phi, 'comp', components(c), 'src', src{s, 1}, ...
                'tau', tau, 'g', g, 'rmsOn', rmsOn, 'offPct', offPct, 'bw', bw, 'ob', ob, 'magErr', mErr); %#ok<SAGROW>

            ls = "-";  if s > 1, ls = "--"; end
            if s == 1
                plot(axT, tg*1e9, xi, 'k', 'LineWidth', 1.9, 'DisplayName', 'ideal OOK');
                plot(axN, f/1e9, nXi, 'k', 'LineWidth', 1.9, 'DisplayName', 'ideal OOK');
                xline(axT, (t_sw + tau)*1e9, ':', '$t_{sw} + \tau$', 'Interpreter', 'latex', ...
                    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
            end
            plot(axT, tg*1e9, x, ls, 'Color', co(s,:), 'LineWidth', 1.4, 'DisplayName', src{s, 1});
            yMax = max([yMax; abs(x); abs(xi)]);
            plot(axN, f/1e9, nX, ls, 'Color', co(s,:), 'LineWidth', 1.4, 'DisplayName', src{s, 1});
        end
        xlabel(axT, 'Time (ns)');  ylabel(axT, sprintf('$%s$ (V)', sym), 'Interpreter', 'latex');
        xlim(axT, [0 tEnd*1e9]);  title(axT, 'Time signal (ideal fitted over ON)');
        legend(axT, 'Location', 'best');
        xlabel(axN, 'Frequency (GHz)');  ylabel(axN, sprintf('$|%s(f)|$ (dB re peak)', sym), 'Interpreter', 'latex');
        xlim(axN, fPlot/1e9);  ylim(axN, [-plot_range_dB 5]);
        xline(axN, [fLoI fHiI]/1e9, ':', 'HandleVisibility', 'off');
        title(axN, 'Spectrum (each normalized to its peak)');  legend(axN, 'Location', 'best');
    end
    set(axTall, 'YLim', 1.1 * yMax * [-1 1]);     % E_theta and E_phi on the same scale
end

%% ---- summary ----------------------------------------------------------------
fprintf(['\nsummary (ideal OOK: -%g dB BW %.4g MHz, %g %% occupied BW %.4g MHz)\n' ...
    '  theta   phi  comp   sim   delay (ns)   scale   rel. RMS ON   OFF energy   -%g dB BW (MHz)   occ. BW (MHz)   in-band |X| error\n'], ...
    bw_level_dB, bwI/1e6, occ_pct, obI/1e6, bw_level_dB);
for r = R
    fprintf('  %5g  %4g  %-5s  %-4s  %9.4g   %7.4g   %8.2f %%   %8.3g %%   %12.4g   %13.4g   %10.2f dB\n', ...
        r.theta, r.phi, r.comp, r.src, r.tau*1e9, r.g, r.rmsOn, r.offPct, r.bw/1e6, r.ob/1e6, r.magErr);
end

%% ---- local functions -----------------------------------------------------
function X = oneSided(X)
% 0 to Nyquist of a full FFT (column)
X = X(1:numel(X)/2 + 1);  X = X(:);
end

function [bw, fl, fh] = levelBW(f, M, lev)
% band around the peak where M is within lev dB of its peak (first crossings
% either side, linear interpolation)
[pk, ip] = max(M);  thr = pk * 10^(-lev/20);
il = find(M(1:ip) < thr, 1, 'last');
iu = ip - 1 + find(M(ip:end) < thr, 1, 'first');
if isempty(il), fl = f(1); else, fl = interp1(M(il:il+1), f(il:il+1), thr); end
if isempty(iu), fh = f(end); else, fh = interp1(M(iu-1:iu), f(iu-1:iu), thr); end
bw = fh - fl;
end

function [bw, fl, fh] = occBW(f, M, pct)
% occupied bandwidth: band holding pct % of the energy, (100 - pct)/2 % cut
% off on each side
E  = cumtrapz(f, M.^2);  E = E / E(end);
q  = (1 - pct/100) / 2;
[Eu, iu] = unique(E);
fl = interp1(Eu, f(iu), q);  fh = interp1(Eu, f(iu), 1 - q);
bw = fh - fl;
end
