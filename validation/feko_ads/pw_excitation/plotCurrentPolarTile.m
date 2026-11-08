function plotCurrentPolarTile(pax, circuit, fekoMag, fekoPhase, colors, ...
    marker, markerSize, axisLabelFontSize, magScale, nMarkers)
%PLOTCURRENTPOLARTILE  Plot the load-current phasor on polar axes, ADS vs FEKO.
%   PLOTCURRENTPOLARTILE(pax, circuit, fekoMag, fekoPhase, colors, marker,
%   markerSize, axisLabelFontSize, magScale, nMarkers) draws onto the polar
%   axes pax one curve per element of the circuit struct array (from
%   importCircuitModelSet): angle = current phase (deg, converted to rad),
%   radius = |I| * magScale (e.g. 1e6 for A -> uA). Over frequency each
%   curve traces the current phasor.
%
%   ADS data is a solid line; the FEKO column with the same tag (from
%   importFekoLoadCurrents) is drawn as markers only in the same colour.
%   nMarkers: [] plots every FEKO sample, otherwise nMarkers markers spaced
%   evenly along the curve (see fekoMarkerIdx).
%
%   fekoMag and fekoPhase must share a frequency grid (rows line up).

    if nargin < 10, nMarkers = []; end
    hold(pax,"on");

    for k = 1:numel(circuit)
        c = circuit(k);
        col = colors(k,:);

        thetaCirc = deg2rad(c.phase);
        rhoCirc = c.mag * magScale;
        polarplot(pax, thetaCirc, rhoCirc, "LineWidth", 1.8, "Color", col, ...
            "LineStyle","-");

        colIdxMag   = find(fekoMag.tags   == c.tag, 1);
        colIdxPhase = find(fekoPhase.tags == c.tag, 1);
        if ~isempty(colIdxMag) && ~isempty(colIdxPhase)
            thetaFeko = deg2rad(fekoPhase.data(:,colIdxPhase));
            rhoFeko = fekoMag.data(:,colIdxMag) * magScale;
            m = fekoMarkerIdx(fekoMag.freq, nMarkers, k, numel(circuit), ...
                rhoFeko.*cos(thetaFeko), rhoFeko.*sin(thetaFeko));
            polarplot(pax, thetaFeko(m), rhoFeko(m), "LineStyle","none", ...
                "Marker", marker, "MarkerSize", markerSize, ...
                "MarkerEdgeColor", col, "MarkerFaceColor", col);
        end
    end

    hold(pax,"off");
    set(pax, "FontSize", axisLabelFontSize);
end