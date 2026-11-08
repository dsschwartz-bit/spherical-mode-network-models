%% plot_transmission_comparison.m
% Tx-loop / Rx-bowtie transmission coefficient, ADS vs FEKO, for several (theta, phi) directions
% of the Rx antenna.
%
% Inputs: one 4-port Touchstone file per direction, angles in the name as "th<deg>" / "ph<deg>",
%   stored as <dataFolder>/<antennaType>/<antennaSubtype>/<file> (see buildDataFilePath.m):
%     ports 1-2: ADS Tx/Rx pair        -> transmission S21
%     ports 3-4: FEKO full-wave pair   -> transmission S43
% Outputs:
%   |S21| error ADS vs FEKO in errBandHz per direction, printed.
%   Figures:
%     transmission_magnitude            |S21| (dB), all directions overlaid
%     transmission_phase_th<..>ph<..>   phase (deg), one figure per direction
%     transmission_legend               shared legend
%   ADS: solid lines; FEKO: thinned filled markers in the same colour (one colour per direction).
% Requires RF Toolbox (sparameters, rfparam).

clear; clc; close all;

scriptFolder = fileparts(mfilename("fullpath"));
run(fullfile(scriptFolder, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(scriptFolder, fullfile(scriptFolder, '..', 'common'));

%% ===== USER SETTINGS =====
dataFolder = scriptFolder;   % root of <antennaType>/<antennaSubtype>/<file>.s4p
s4pFiles = {
    "loop_bowtie_ph30th60r5.s4p"
    "loop_bowtie_ph60th60r5.s4p"
    "loop_bowtie_ph90th60r5.s4p"
};

freqUnit              = "MHz";       % "Hz", "MHz", or "GHz"
freqLim               = [0 2500];    % x-axis limits in freqUnit
fekoMarker            = "d";         % marker style used for all FEKO points
fekoMagDownsampleHz   = 50e6;        % FEKO magnitude markers every this frequency step
fekoPhaseDownsampleHz = 25e6;        % FEKO phase markers every this frequency step
namePrefix            = "transmission";
errBandHz             = [0.1e9 2e9]; % band of the printed |S21| error, ADS vs FEKO (Hz)

S = plotStyle();
% =========================

%% ---- IMPORT -----------------------------------------------------
nCurves = numel(s4pFiles);
for k = nCurves:-1:1
    dataset(k) = buildTransmissionDataset(buildDataFilePath(dataFolder, s4pFiles{k}), ...
        fekoMagDownsampleHz, fekoPhaseDownsampleHz);
end

switch freqUnit
    case "GHz", freqScale = 1e-9;
    case "MHz", freqScale = 1e-6;
    otherwise,  freqScale = 1;
end
colors = S.colors(mod(0:nCurves-1, size(S.colors, 1)) + 1, :);

%% ---- MAGNITUDE (all angles) ---------------------------------------
[fig, ax] = newFigure(namePrefix + "_magnitude");
[legH, legLabels] = plotTransmissionTile(ax, dataset, "mag", colors, freqScale, ...
    fekoMarker, S.markerSize, S.fontSize);
hold(ax, 'on');
styleAxes(ax, S, "Frequency (" + freqUnit + ")", '$|S_{21}|$ (dB)');
maxMagAll = max(arrayfun(@(d) max([d.ads.mag(:); d.fekoDS.mag.data(:)]), dataset));
xlim(ax, freqLim);  ylim(ax, [-80, 5 + maxMagAll]);

%% ---- PHASE (one figure per angle) ---------------------------------
for k = 1:nCurves
    d = dataset(k);
    name = sprintf('%s_phase_th%sph%s', namePrefix, d.theta, d.phi);
    [fig, ax] = newFigure(name);
    plotSingleAnglePhaseTile(ax, d, colors(k,:), freqScale, fekoMarker, S.markerSize, S.fontSize);
    styleAxes(ax, S, "Frequency (" + freqUnit + ")", 'Phase ($^\circ$)');
    xlim(ax, freqLim);  ylim(ax, [-180 180]);  yticks(ax, -180:90:180);
    obsLabel(ax, S, d.theta, d.phi);    % observation point, top-right corner
end

%% ---- SEPARATE LEGEND: one row per (theta,phi), ADS | FEKO side by side (nCurves x 2) ----------
legLabels = string(legLabels);
iAds  = find(startsWith(legLabels, "ADS"));
iFeko = find(startsWith(legLabels, "FEKO"));
ord = zeros(1, 2*numel(iAds));
for k = 1:numel(iAds)                      % pair each ADS entry with the FEKO entry of the same angle
    angleTxt = extractAfter(legLabels(iAds(k)), "ADS");
    ord(2*k-1:2*k) = [iAds(k), iFeko(extractAfter(legLabels(iFeko), "FEKO") == angleTxt)];
end
legendFigure(namePrefix + "_legend", legH(ord), legLabels(ord), 2);   % fills row by row

%% ---- |S21| ERROR, ADS (s4p S21) vs FEKO (s4p S43), full-resolution s4p data in errBandHz ----------
fprintf('\n|S21| error ADS vs FEKO, %.3g-%.3g GHz:\n', errBandHz/1e9);
worst = struct('err', -Inf, 'name', "", 'f', NaN);  eAll = [];
for k = 1:nCurves
    sp = sparameters(buildDataFilePath(dataFolder, s4pFiles{k}));
    f  = sp.Frequencies;  b = f >= errBandHz(1) - 1 & f <= errBandHz(2) + 1;
    fe = abs(squeeze(sp.Parameters(4,3,b)));  ad = abs(squeeze(sp.Parameters(2,1,b)));  fb = f(b);
    eDB = abs(20*log10(ad./fe));  [mx, i] = max(eDB);
    relMSE = sum((ad - fe).^2) / sum(fe.^2);
    lbl = sprintf('theta %s, phi %s', dataset(k).theta, dataset(k).phi);
    eAll = [eAll; eDB]; %#ok<AGROW>
    fprintf('  %-18s mean %.2f dB, max %.2f dB at %6.1f MHz (FEKO %.1f dB, ADS %.1f dB), max rel. err %.1f %%, rel. MSE %.3f %%\n', ...
        lbl, mean(eDB), mx, fb(i)/1e6, 20*log10(fe(i)), 20*log10(ad(i)), 100*max(abs(ad - fe)./fe), 100*relMSE);
    if mx > worst.err, worst = struct('err', mx, 'name', string(lbl), 'f', fb(i)); end
end
fprintf('  All angles: mean %.2f dB, max %.2f dB (%s, %.1f MHz)\n', mean(eAll), worst.err, worst.name, worst.f/1e6);
