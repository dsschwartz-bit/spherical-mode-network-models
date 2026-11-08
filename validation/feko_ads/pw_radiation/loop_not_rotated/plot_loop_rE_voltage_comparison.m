%% plot_loop_rE_voltage_comparison.m
% rE far-field voltage (Vtheta, Vphi) of the loop antenna (axis along z, no loads), ADS circuit
% model vs FEKO, for three phi cuts. Standalone variant of ../plot_rE_voltage_comparison.m
% (loopCase = "not_rotated") with its own marker settings and figure names; the ADS sweep
% variable (phi_deg) is in degrees.
%
% Inputs (in dataFolder), one file per (quantity q = mag|phase, component c = Vtheta|Vphi):
%   ADS  - loop_noloads_ADS_rE_<q>_<c>.csv, long format with a quoted multi-line header:
%          phi_deg, freq [Hz], value -- the phi cuts are stacked, not in columns.
%   FEKO - loop_noloads_FEKO_rE_<q>_<c>.dat, tab-delimited, one column per cut headed
%          "Phi = 60 deg [V]" ...; curves are matched by phi value, never by position.
% Outputs:
%   Mean / largest complex-voltage error |E_ADS - E_FEKO| per component and cut, printed.
%   Figures:
%     loop_rE_voltage_<comp>_magnitude | _phase | _phasor, and loop_rE_voltage_legend.

clear; clc; close all;

scriptFolder = fileparts(mfilename("fullpath"));
if isempty(scriptFolder), scriptFolder = pwd; end
run(fullfile(scriptFolder, '..', '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
% Helper functions (importFieldData, plotComparisonTile, ...) are in the parent folder pw_radiation.
addpath(scriptFolder, fullfile(scriptFolder, '..'), fullfile(scriptFolder, '..', '..', 'common'));

%% ===== USER SETTINGS =====
dataFolder = scriptFolder;          % folder containing the loop_noloads_* exports

components = ["Vtheta", "Vphi"];    % rE far-field voltage components to plot
phiList    = [30, 60, 90];          % phi cuts [deg] -- 3 curves per axes

% Filename builders: quantity q is "mag" or "phase", component c is e.g. "Vtheta"
adsFile  = @(q,c) "loop_noloads_ADS_rE_"  + q + "_" + c + ".csv";
fekoFile = @(q,c) "loop_noloads_FEKO_rE_" + q + "_" + c + ".dat";
outPrefix = "loop_rE_voltage";   % figure names: <outPrefix>_<comp>_magnitude|phase|phasor, <outPrefix>_legend

freqUnit = "GHz";   % "Hz", "MHz" or "GHz" -- x-axis unit
magUnit  = "mV";    % "mV" (x1e3) or "V" (x1)  -- rE voltage is exported in V

fekoMarker       = "^";      % marker style for all FEKO points
fekoMarkerStepHz = 25e6;    % [] = every FEKO sample as a marker (no thinning); else marker step [Hz]
fekoPolarMarkers = 49;      % [] = every FEKO sample on the phasor plots; else markers per trace (49 = 50 MHz over 0.1-2.5 GHz)

% If the two solvers were swept over different bands, true clips both to the
% overlap so the comparison is like-for-like (matters most on the polar plot,
% where the extra FEKO points would look like disagreement).
restrictToCommonBand = true;
freqMaxHz = 2.5e9;        % upper plot/comparison limit [Hz] (the FEKO export runs to 6 GHz; ADS stops at 2.5 GHz)

% Phase is plotted as exported: no sign flip, offset, wrapping or unwrapping.

% Angular unit of the phi sweep variable in each source. phiList is in degrees; the phase
% data are in degrees in both sources and are not converted.
adsPhiUnit  = "deg";                % phi_deg sweep in the ADS schematic
fekoPhiUnit = "deg";

% Fallback for files whose unit does not match the declaration above (e.g. an
% export that had already converted to degrees, which then reads as
% 1719/3438/5157). Applied after the unit conversion:
%   "auto" - rescale only when a known factor reproduces phiList exactly
%   1      - off; use the values as they stand
%   <num>  - apply that factor to every file
% "auto" prints a line naming any file it rescales.
phiScaleMode = "auto";

S = plotStyle();             % line widths, markers, colours (validation/feko_ads/common/plotStyle.m)
% =========================

%% ---- SCALING -------------------------------------------------------
switch freqUnit
    case "GHz", freqScale = 1e-9;
    case "MHz", freqScale = 1e-6;
    otherwise,  freqScale = 1;
end

switch magUnit
    case "mV", magScale = 1e3;   % V -> mV
    otherwise, magScale = 1;     % V, as exported
end

magUnitLatex = ['\mathrm{' char(magUnit) '}'];

colors = S.colors(mod((0:numel(phiList)-1), size(S.colors,1)) + 1, :);   % Wong (2011)

%% ---- IMPORT + PREPARE EVERY COMPONENT --------------------------------
nComp = numel(components);
D = struct("comp",{}, "compSym",{}, "compLatex",{}, "adsMag",{}, ...
    "adsPhase",{}, "fekoMag",{}, "fekoPhase",{}, "fLim",{});

for ci = 1:nComp
    comp = components(ci);
    [compSym, compLatex] = componentLatex(comp);   % 'rE_{\theta}' and '$rE_{\theta}$'

    % importFieldData picks the parser from the file content and warns if it
    % does not match the expected source.
    adsMag    = importFieldData(fullfile(dataFolder, adsFile("mag",   comp)), "ADS");
    adsPhase  = importFieldData(fullfile(dataFolder, adsFile("phase", comp)), "ADS");
    fekoMag   = importFieldData(fullfile(dataFolder, fekoFile("mag",   comp)), "FEKO");
    fekoPhase = importFieldData(fullfile(dataFolder, fekoFile("phase", comp)), "FEKO");

    % Sweep variable to degrees, then report any requested cut still missing.
    adsMag    = normalizePhi(adsMag,    phiList, phiScaleMode, adsPhiUnit);
    adsPhase  = normalizePhi(adsPhase,  phiList, phiScaleMode, adsPhiUnit);
    fekoMag   = normalizePhi(fekoMag,   phiList, phiScaleMode, fekoPhiUnit);
    fekoPhase = normalizePhi(fekoPhase, phiList, phiScaleMode, fekoPhiUnit);

    checkPhiCoverage(adsMag,    phiList);
    checkPhiCoverage(adsPhase,  phiList);
    checkPhiCoverage(fekoMag,   phiList);
    checkPhiCoverage(fekoPhase, phiList);

    fprintf("%s: ADS phi = [%s] deg, %g-%g GHz (%d pts) | FEKO phi = [%s] deg, %g-%g GHz (%d pts)\n", ...
        comp, join(string(adsMag.phi), ", "), ...
        min(adsMag.freq)/1e9, max(adsMag.freq)/1e9, numel(adsMag.freq), ...
        join(string(fekoMag.phi), ", "), ...
        min(fekoMag.freq)/1e9, max(fekoMag.freq)/1e9, numel(fekoMag.freq));

    % ---- Common frequency band ------------------------------------------
    if restrictToCommonBand
        fLim = [max(min(adsMag.freq), min(fekoMag.freq)), ...
                min([max(adsMag.freq), max(fekoMag.freq), freqMaxHz])];
        adsMag    = restrictFreqRange(adsMag,    fLim);
        adsPhase  = restrictFreqRange(adsPhase,  fLim);
        fekoMag   = restrictFreqRange(fekoMag,   fLim);
        fekoPhase = restrictFreqRange(fekoPhase, fLim);
    else
        fLim = [min([adsMag.freq; fekoMag.freq]), min(max([adsMag.freq; fekoMag.freq]), freqMaxHz)];
        adsMag    = restrictFreqRange(adsMag,    fLim);
        adsPhase  = restrictFreqRange(adsPhase,  fLim);
        fekoMag   = restrictFreqRange(fekoMag,   fLim);
        fekoPhase = restrictFreqRange(fekoPhase, fLim);
    end

    D(ci).comp      = comp;
    D(ci).compSym   = compSym;
    D(ci).compLatex = compLatex;
    D(ci).adsMag    = adsMag;
    D(ci).adsPhase  = adsPhase;
    D(ci).fekoMag   = fekoMag;
    D(ci).fekoPhase = fekoPhase;
    D(ci).fLim      = fLim;
end

%% ---- MEAN VOLTAGE ERROR (ADS vs FEKO), per observation point ----
% mean_f |E_ADS(f) - E_FEKO(f)| of the complex rE (E = |rE| exp(j*phase), so magnitude and phase
% errors both count), in magUnit, over the FEKO frequencies in each band. ADS (10 MHz grid) is
% interpolated onto the FEKO frequencies (25 MHz grid) in real and imaginary parts. Bands: start of
% the common band up to each value of mseUpperHz (the last one is the whole range, i.e. freqMaxHz).
% "all phi" = mean over every frequency of all the cuts.
mseUpperHz = [1.8e9, freqMaxHz];
fprintf('\nMean voltage error mean|E_ADS - E_FEKO| (%s; FEKO frequencies, ADS interpolated):\n', magUnit);
hdr = sprintf('%-22s', 'component, phi');
for b = 1:numel(mseUpperHz), hdr = [hdr sprintf('  %13s', sprintf('%.1f-%.1f GHz', D(1).fLim(1)/1e9, min(mseUpperHz(b), D(1).fLim(2))/1e9))]; end %#ok<AGROW>
fprintf('%s\n', hdr);
for ci = 1:nComp
    d = D(ci);
    errAll = cell(1, numel(mseUpperHz));                % errors of all cuts, per band
    for p = phiList
        cplx = @(M, P) M.data(:, findPhiColumn(M.phi, p)) .* exp(1j*deg2rad(P.data(:, findPhiColumn(P.phi, p))));   % same phi tolerance as the plots
        eF = cplx(d.fekoMag, d.fekoPhase);  fF = d.fekoMag.freq;
        eA = cplx(d.adsMag,  d.adsPhase);   fA = d.adsMag.freq;
        eA = interp1(fA, real(eA), fF, 'spline') + 1j*interp1(fA, imag(eA), fF, 'spline');
        err = abs(eA - eF) * magScale;
        row = sprintf('%-22s', sprintf('E_%s, phi = %g deg', erase(d.comp, "V"), p));
        for b = 1:numel(mseUpperHz)
            in = fF <= mseUpperHz(b) + 1;                  % +1 Hz: include the band edge
            errAll{b} = [errAll{b}; err(in)];
            row = [row sprintf('  %#13.3g', mean(err(in)))]; %#ok<AGROW>
        end
        fprintf('%s\n', row);
    end
    fprintf('%-22s%s\n', sprintf('E_%s, all phi', erase(d.comp, "V")), sprintf('  %#13.3g', cellfun(@mean, errAll)));
end
fprintf('\n');

%% ---- LARGEST VOLTAGE ERROR (ADS vs FEKO), per observation point ----
% max_f |E_ADS(f) - E_FEKO(f)| of the complex rE (magnitude and phase together), in magUnit,
% over the FEKO frequencies in each band (ADS interpolated as above), with the frequency where
% it occurs; max|E_FEKO| of the cut in the band is printed alongside for scale.
% "all phi" = worst over the cuts.
fprintf('Largest voltage error max|E_ADS - E_FEKO| (%s, at frequency)  [max|E_FEKO| in the band]:\n', magUnit);
hdrM = sprintf('%-22s', 'component, phi');
for b = 1:numel(mseUpperHz), hdrM = [hdrM sprintf('  %22s', sprintf('%.1f-%.1f GHz', D(1).fLim(1)/1e9, min(mseUpperHz(b), D(1).fLim(2))/1e9))]; end %#ok<AGROW>
fprintf('%s\n', hdrM);
for ci = 1:nComp
    d = D(ci);
    worst = -inf(1, numel(mseUpperHz));  fWorst = nan(1, numel(mseUpperHz));
    for p = phiList
        cplx = @(M, P) M.data(:, findPhiColumn(M.phi, p)) .* exp(1j*deg2rad(P.data(:, findPhiColumn(P.phi, p))));   % same phi tolerance as the plots
        eF = cplx(d.fekoMag, d.fekoPhase);  fF = d.fekoMag.freq;
        eA = cplx(d.adsMag,  d.adsPhase);   fA = d.adsMag.freq;
        eA = interp1(fA, real(eA), fF, 'spline') + 1j*interp1(fA, imag(eA), fF, 'spline');
        err = abs(eA - eF) * magScale;
        row = sprintf('%-22s', sprintf('E_%s, phi = %g deg', erase(d.comp, "V"), p));
        for b = 1:numel(mseUpperHz)
            in = find(fF <= mseUpperHz(b) + 1);
            [mx, k] = max(err(in));
            row = [row sprintf('  %#5.3g (%4.2f GHz) [%.3g]', mx, fF(in(k))/1e9, max(abs(eF(in)))*magScale)]; %#ok<AGROW>
            if mx > worst(b), worst(b) = mx;  fWorst(b) = fF(in(k)); end
        end
        fprintf('%s\n', row);
    end
    row = sprintf('%-22s', sprintf('E_%s, all phi', erase(d.comp, "V")));
    for b = 1:numel(mseUpperHz), row = [row sprintf('  %#5.3g (%4.2f GHz)      ', worst(b), fWorst(b)/1e9)]; end %#ok<AGROW>
    fprintf('%s\n', row);
end
fprintf('\n');

%% ---- FIGURES: one panel per figure, component named in the y-label, shared legend ----
legH = gobjects(0);  legLabels = strings(0);
for ci = 1:nComp
    d = D(ci);
    base = outPrefix + "_" + d.comp;

    % ---- Magnitude ------------------------------------------------------
    [fig, ax] = newFigure(base + "_magnitude");
    [h, lbl] = plotComparisonTile(ax, d.adsMag, d.fekoMag, phiList, colors, freqScale, ...
        fekoMarker, S.markerSize, magScale, S.fontSize, fekoMarkerStepHz);
    if numel(h) > numel(legH), legH = h; legLabels = lbl; end
    styleAxes(ax, S, "Frequency (" + freqUnit + ")", ['$|' d.compSym '|$ ($' magUnitLatex '$)']);
    xlim(ax, d.fLim * freqScale);

    % ---- Phase ----------------------------------------------------------
    [fig, ax] = newFigure(base + "_phase");
    plotComparisonTile(ax, d.adsPhase, d.fekoPhase, phiList, colors, freqScale, ...
        fekoMarker, S.markerSize, 1, S.fontSize, fekoMarkerStepHz);
    styleAxes(ax, S, "Frequency (" + freqUnit + ")", ['$\angle ' d.compSym '$ ($^\circ$)']);
    xlim(ax, d.fLim * freqScale);
    phaseAll = [d.adsPhase.data(:); d.fekoPhase.data(:)];  phaseAll = phaseAll(~isnan(phaseAll));
    if ~isempty(phaseAll) && min(phaseAll) >= -180 && max(phaseAll) <= 180   % never clip unwrapped data
        ylim(ax, [-180 180]);  yticks(ax, -180:90:180);
    end

    % ---- Phasor (polar) -------------------------------------------------
    [fig, pax] = newFigure(base + "_phasor", "polar");
    plotComparisonPolarTile(pax, d.adsMag, d.adsPhase, d.fekoMag, d.fekoPhase, phiList, ...
        colors, fekoMarker, S.markerSize, S.fontSize, magScale, fekoPolarMarkers);
    styleAxes(pax, S);
    allMag = [d.adsMag.data(:); d.fekoMag.data(:)];
    maxMag = max(allMag(~isnan(allMag))) * magScale;
    pax.RLim = [0, 1.1*maxMag];
    pax.ThetaTickLabel = compose('$%d^\\circ$', pax.ThetaTick);
    polarLabel(fig, S, ['$|' d.compSym '|$ ($' magUnitLatex '$)']);
end

% ---- Shared legend: one column per phi, ADS on the top row, FEKO underneath ----
isAds = startsWith(legLabels, "ADS");
legendFigure(outPrefix + "_legend", [legH(isAds), legH(~isAds)], ...
    [legLabels(isAds), legLabels(~isAds)], nnz(isAds));
