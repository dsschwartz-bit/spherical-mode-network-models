function [legH, legLabels] = plotComparisonPolarTile(pax, adsMag, adsPhase, ...
    fekoMag, fekoPhase, phiList, colors, marker, markerSize, ...
    axisLabelFontSize, magScale, nMarkers)
%PLOTCOMPARISONPOLARTILE  Plot the rE far-field voltage phasor on polar axes, ADS vs FEKO.
%   [legH, legLabels] = PLOTCOMPARISONPOLARTILE(pax, adsMag, adsPhase,
%   fekoMag, fekoPhase, phiList, colors, marker, markerSize,
%   axisLabelFontSize, magScale, nMarkers) draws onto the polar axes pax one
%   curve per entry of phiList: angle = rE voltage phase, radius =
%   |rE| * magScale (1e3 for V -> mV). Over frequency each curve traces the
%   voltage phasor.
%
%   ADS is a solid line; FEKO is markers only in the same colour. Columns are
%   matched to phiList by phi value. nMarkers: [] plots every FEKO sample,
%   otherwise nMarkers markers spaced evenly along the curve (fekoMarkerIdx).
%
%   If the magnitude and phase structs of a source are on different
%   frequency grids, the phase is interpolated onto the magnitude grid over
%   their overlap.

    if nargin < 12, nMarkers = []; end
    legH = gobjects(0);
    legLabels = strings(0);

    hold(pax, "on");

    for k = 1:numel(phiList)
        p   = phiList(k);
        col = colors(mod(k-1, size(colors,1)) + 1, :);

        % --- ADS ----------------------------------------------------
        iam = findPhiColumn(adsMag.phi,   p);
        iap = findPhiColumn(adsPhase.phi, p);
        if ~isempty(iam) && ~isempty(iap)
            [rho, theta] = alignMagPhase(adsMag, adsPhase, iam, iap, magScale);
            h = polarplot(pax, theta, rho, "LineStyle","-", "LineWidth",1.8, ...
                "Color",col);
            legH(end+1)      = h;                                      %#ok<AGROW>
            legLabels(end+1) = "ADS, " + phiLatexLabel(p);              %#ok<AGROW>
        end

        % --- FEKO ---------------------------------------------------
        ifm = findPhiColumn(fekoMag.phi,   p);
        ifp = findPhiColumn(fekoPhase.phi, p);
        if ~isempty(ifm) && ~isempty(ifp)
            [rho, theta] = alignMagPhase(fekoMag, fekoPhase, ifm, ifp, magScale);
            m = fekoMarkerIdx((1:numel(rho)).', nMarkers, k, numel(phiList), ...
                rho.*cos(theta), rho.*sin(theta));
            h = polarplot(pax, theta(m), rho(m), "LineStyle","none", ...
                "Marker",marker, "MarkerSize",markerSize, ...
                "MarkerEdgeColor",col, "MarkerFaceColor",col);
            legH(end+1)      = h;                                      %#ok<AGROW>
            legLabels(end+1) = "FEKO, " + phiLatexLabel(p);             %#ok<AGROW>
        end
    end

    hold(pax, "off");

    grid(pax, "on");
    set(pax, "ThetaGrid","on", "RGrid","on", ...
        "FontSize", axisLabelFontSize, "TickLabelInterpreter","latex");
end

% ----------------------------------------------------------------------
function [rho, theta] = alignMagPhase(magS, phaseS, iMag, iPhase, magScale)
%ALIGNMAGPHASE  Put a magnitude column and a phase column on one grid.

    fMag   = magS.freq;
    fPhase = phaseS.freq;

    if numel(fMag) == numel(fPhase) && ...
            max(abs(fMag - fPhase)) <= 1e-6 * max(abs(fMag))
        rho   = magS.data(:, iMag) * magScale;
        theta = deg2rad(phaseS.data(:, iPhase));
        return;
    end

    % Grids differ: interpolate the phase onto the magnitude grid over the
    % overlap. The unwrapped phase is interpolated so no sample is taken
    % across a +-180 deg jump; on polar axes the plotted angle is the same.
    ph = unwrap(deg2rad(phaseS.data(:, iPhase)));
    keep = fMag >= min(fPhase) & fMag <= max(fPhase);

    rho   = magS.data(keep, iMag) * magScale;
    theta = interp1(fPhase, ph, fMag(keep), "linear");
end