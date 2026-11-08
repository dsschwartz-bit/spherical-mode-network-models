function [legH, legLabels] = plotCurrentTile(ax, circuit, feko, field, ...
    colors, freqScale, fekoMarker, fekoMarkerSize, valueScale, axisLabelFontSize, markerStepHz)
%PLOTCURRENTTILE  Plot ADS load-current lines and matching FEKO markers vs frequency.
%   [legH, legLabels] = PLOTCURRENTTILE(ax, circuit, feko, field, colors,
%   freqScale, fekoMarker, fekoMarkerSize, valueScale, axisLabelFontSize,
%   markerStepHz) draws onto axes ax one solid line per element of the
%   circuit struct array (from importCircuitModelSet) and, for the FEKO
%   column with the same tag (from importFekoLoadCurrents), a marker-only
%   curve in the same colour.
%
%   field             - "mag" or "phase" (field of circuit to plot)
%   freqScale         - multiplies the frequency axis (e.g. 1e-9 for GHz)
%   valueScale        - multiplies ADS and FEKO values (e.g. 1e6 for A -> uA;
%                       1 for phase in degrees). Default 1.
%   axisLabelFontSize - tick-label font size
%   markerStepHz      - FEKO marker step [Hz], staggered per load
%                       (see fekoMarkerIdx); [] plots every sample
%
%   No legend is drawn: the handles and LaTeX labels (tagToLatexLabel) are
%   returned so the caller can build one shared legend.

    if nargin < 9
        valueScale = 1;
    end
    if nargin < 10
        axisLabelFontSize = get(groot,"defaultAxesFontSize");
    end
    if nargin < 11
        markerStepHz = [];
    end

    hold(ax,"on");
    legH = gobjects(0); legLabels = {};

    for k = 1:numel(circuit)
        c = circuit(k);
        col = colors(k,:);

        h1 = plot(ax, c.freq*freqScale, c.(field)*valueScale, ...
            "LineWidth", 1.8, "Color", col, "LineStyle","-");
        legH(end+1) = h1; %#ok<AGROW>
        legLabels{end+1} = tagToLatexLabel(c.tag, "ADS"); %#ok<AGROW>

        colIdx = find(feko.tags == c.tag, 1);
        if ~isempty(colIdx)
            m = fekoMarkerIdx(feko.freq, markerStepHz, k, numel(circuit));  % staggered per load
            h2 = plot(ax, feko.freq(m)*freqScale, feko.data(m,colIdx)*valueScale, ...
                "LineStyle","none", "Marker", fekoMarker, ...
                "MarkerSize", fekoMarkerSize, ...
                "MarkerEdgeColor", col, "MarkerFaceColor", col);
            legH(end+1) = h2; %#ok<AGROW>
            legLabels{end+1} = tagToLatexLabel(c.tag, "FEKO"); %#ok<AGROW>
        end
    end

    hold(ax,"off");
    grid(ax,"on"); box(ax,"on");
    set(ax, "FontSize", axisLabelFontSize);
end