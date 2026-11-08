%% plot_current_comparison.m
% Load current of the bowtie (45 deg edge port) under plane-wave excitation, ADS circuit model vs
% FEKO full wave, for several incidence directions (phi, theta).
%
% Inputs (in dataFolder):
%   ADS  - one CSV per incidence direction: freq [Hz], mag(I_Probe1.i) [A], phase(I_Probe1.i) [deg]
%   FEKO - two .dat load-current exports (magnitude / phase), one column per direction,
%          headed "Load_<tag> [A]" / "Load_<tag> [deg]"
%   The tag (e.g. ph30th60 = incidence phi / theta) appears in the CSV names and the FEKO
%   headers; the two sources are matched by tag.
% Outputs:
%   Mean / largest complex-current error |I_ADS - I_FEKO| per direction, printed.
%   Figures:
%     load_current_magnitude   |I| (uA) vs frequency
%     load_current_phase       phase (deg) vs frequency
%     load_current_phasor      polar plot of the phasor I (uA) over frequency
%     load_current_legend      shared legend
%   ADS: solid lines; FEKO: thinned filled markers in the same colour (one colour per direction).

clear; clc; close all;

scriptFolder = fileparts(mfilename("fullpath"));
run(fullfile(scriptFolder, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(scriptFolder, fullfile(scriptFolder, '..', 'common'));

%% ===== USER SETTINGS =====
dataFolder = scriptFolder;   % folder containing the ADS .csv and FEKO .dat exports listed below

csvFiles = {
    "bowtie_edgeport_45degrees_ph30th60_Eph_p707_Eth_p707.csv";
    "bowtie_edgeport_45degrees_ph60th60_Eph_p707_Eth_p707.csv";
    "bowtie_edgeport_45degrees_ph90th60_Eph_p707_Eth_p707.csv"
};
fekoMagFile   = "bowtie_edgeport_45degrees_FEKO_loadcurrents_magnitude.dat";
fekoPhaseFile = "bowtie_edgeport_45degrees_FEKO_loadcurrents_phase.dat";

freqUnit         = "GHz";  % "Hz", "MHz" or "GHz"
freqLim          = [0 2.5];% x-axis limits in freqUnit ([] = data range)
fekoMarker       = "d";    % FEKO marker
fekoMarkerStepHz = 125e6;  % FEKO markers every this frequency step (staggered per direction)
fekoPolarMarkers = 22;     % FEKO markers per trace on the phasor plot (evenly along the curve)
namePrefix       = "load_current";

S = plotStyle();            % line widths, markers, colours (edit there or override here)
% =========================

%% ---- IMPORT -----------------------------------------------------
circuit   = importCircuitModelSet(dataFolder, csvFiles);
fekoMag   = importFekoLoadCurrents(fullfile(dataFolder, fekoMagFile));     % full resolution; markers are
fekoPhase = importFekoLoadCurrents(fullfile(dataFolder, fekoPhaseFile));   % thinned while plotting

switch freqUnit
    case "GHz", freqScale = 1e-9;
    case "MHz", freqScale = 1e-6;
    otherwise,  freqScale = 1;
end
nLoads = numel(circuit);
colors = S.colors(mod(0:nLoads-1, size(S.colors, 1)) + 1, :);

%% ---- LOAD-CURRENT ERROR (ADS vs FEKO), per incidence direction ----
% Complex current I = |I| exp(j*phase), so magnitude and phase errors both count. ADS (10 MHz grid)
% is interpolated onto the FEKO frequencies (25 MHz grid) in real and imaginary parts. Bands: from
% the first common frequency up to each value of errUpperHz (the last one is the whole plotted range).
% "all" = over every frequency of all the directions.
errUpperHz = [1.8e9, freqLim(end)/freqScale];
errScale = 1e6;  errUnit = "uA";
nB = numel(errUpperHz);  f0 = max(fekoMag.freq(1), min(arrayfun(@(c) c.freq(1), circuit)));
bandTxt = arrayfun(@(u) sprintf('%.1f-%.1f GHz', f0/1e9, u/1e9), errUpperHz, 'UniformOutput', false);
E = struct('tag', {}, 'f', {}, 'err', {}, 'ref', {});
for k = 1:nLoads
    ci = find(fekoMag.tags == circuit(k).tag, 1);  if isempty(ci), continue; end
    fF = fekoMag.freq;  keep = fF >= f0 & fF <= errUpperHz(end) + 1;  fF = fF(keep);
    eF = fekoMag.data(keep, ci) .* exp(1j*deg2rad(fekoPhase.data(keep, ci)));
    eA = circuit(k).mag .* exp(1j*deg2rad(circuit(k).phase));
    eA = interp1(circuit(k).freq, real(eA), fF, 'spline') + 1j*interp1(circuit(k).freq, imag(eA), fF, 'spline');
    E(end+1) = struct('tag', circuit(k).tag, 'f', fF, 'err', abs(eA - eF)*errScale, 'ref', abs(eF)*errScale); %#ok<SAGROW>
end
lbl = @(t) sprintf('I, %s', regexprep(t, 'ph(\d+)th(\d+)', 'phi = $1, theta = $2 deg'));

fprintf('\nMean current error mean|I_ADS - I_FEKO| (%s; FEKO frequencies, ADS interpolated):\n', errUnit);
fprintf('%-32s%s\n', 'direction', sprintf('  %13s', bandTxt{:}));
allErr = cell(1, nB);
for k = 1:numel(E)
    row = sprintf('%-32s', lbl(E(k).tag));
    for b = 1:nB
        in = E(k).f <= errUpperHz(b) + 1;  allErr{b} = [allErr{b}; E(k).err(in)];
        row = [row sprintf('  %#13.3g', mean(E(k).err(in)))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end
fprintf('%-32s%s\n', 'I, all directions', sprintf('  %#13.3g', cellfun(@mean, allErr)));

fprintf('\nLargest current error max|I_ADS - I_FEKO| (%s, at frequency)  [max|I_FEKO| in the band]:\n', errUnit);
fprintf('%-32s%s\n', 'direction', sprintf('  %22s', bandTxt{:}));
worst = repmat([-Inf NaN], nB, 1);
for k = 1:numel(E)
    row = sprintf('%-32s', lbl(E(k).tag));
    for b = 1:nB
        in = find(E(k).f <= errUpperHz(b) + 1);  [m, i] = max(E(k).err(in));  fm_ = E(k).f(in(i));
        if m > worst(b, 1), worst(b, :) = [m fm_]; end
        row = [row sprintf('  %#5.3g (%4.2f GHz) [%.3g]', m, fm_/1e9, max(E(k).ref(in)))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end
row = sprintf('%-32s', 'I, all directions');
for b = 1:nB, row = [row sprintf('  %#5.3g (%4.2f GHz)      ', worst(b, 1), worst(b, 2)/1e9)]; end %#ok<AGROW>
fprintf('%s\n\n', row);

%% ---- MAGNITUDE --------------------------------------------------
[fig, ax] = newFigure(namePrefix + "_magnitude");
[legH, legLabels] = plotCurrentTile(ax, circuit, fekoMag, "mag", colors, freqScale, ...
    fekoMarker, S.markerSize, 1e6, S.fontSize, fekoMarkerStepHz);  % A -> uA
styleAxes(ax, S, "Frequency (" + freqUnit + ")", '$|I|$ ($\mu$A)');
if ~isempty(freqLim), xlim(ax, freqLim); end

%% ---- PHASE ------------------------------------------------------
[fig, ax] = newFigure(namePrefix + "_phase");
plotCurrentTile(ax, circuit, fekoPhase, "phase", colors, freqScale, fekoMarker, S.markerSize, 1, S.fontSize, ...
    fekoMarkerStepHz);
styleAxes(ax, S, "Frequency (" + freqUnit + ")", 'Phase ($^\circ$)');
ylim(ax, [-180 180]);  yticks(ax, -180:90:180);
if ~isempty(freqLim), xlim(ax, freqLim); end

%% ---- PHASOR (polar) ---------------------------------------------
[fig, pax] = newFigure(namePrefix + "_phasor", "polar");
[cP, fmP, fpP] = clipBand(circuit, fekoMag, fekoPhase, freqLim/freqScale);   % same band as the other panels
plotCurrentPolarTile(pax, cP, fmP, fpP, colors, fekoMarker, S.markerSize, S.fontSize, 1e6, ...
    fekoPolarMarkers);
styleAxes(pax, S);
maxI = max(cellfun(@(m) max(m(:)), {cP.mag}));  maxI = max(maxI, max(fmP.data(:)))*1e6;
pax.RLim = [0, 1.1*maxI];
pax.ThetaTickLabel = compose('$%d^\\circ$', pax.ThetaTick);
polarLabel(fig, S, '$I$ ($\mu$A)');     % phasor plot: the complex current I itself

%% ---- SEPARATE LEGEND: one column per direction, ADS on the top row, FEKO underneath --------------
isAds = contains(string(legLabels), "ADS");
legendFigure(namePrefix + "_legend", [legH(isAds), legH(~isAds)], ...
    [string(legLabels(isAds)), string(legLabels(~isAds))], nnz(isAds));

function [c, fm, fp] = clipBand(c, fm, fp, band)
% Keep only the samples inside band = [fmin fmax] (Hz); [] keeps everything.
if isempty(band), return; end
for k = 1:numel(c)
    in = c(k).freq >= band(1) & c(k).freq <= band(2);
    c(k).freq = c(k).freq(in);  c(k).mag = c(k).mag(in);  c(k).phase = c(k).phase(in);
end
in = fm.freq >= band(1) & fm.freq <= band(2);  fm.freq = fm.freq(in);  fm.data = fm.data(in, :);
in = fp.freq >= band(1) & fp.freq <= band(2);  fp.freq = fp.freq(in);  fp.data = fp.data(in, :);
end
