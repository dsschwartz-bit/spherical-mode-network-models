function [h1, h2] = plotSingleAnglePhaseTile(ax, d, col, freqScale, marker, ...
    markerSize, axisLabelFontSize)
%PLOTSINGLEANGLEPHASETILE  Plot ADS and FEKO phase vs frequency for a
%   single (theta, phi) observation angle onto axes ax.
%   [h1, h2] = PLOTSINGLEANGLEPHASETILE(ax, d, col, freqScale, marker, markerSize,
%   axisLabelFontSize) where d is a single dataset entry (as returned by
%   buildTransmissionDataset), and col is the [R G B] color shared by
%   both the ADS line (h1) and the FEKO markers (h2).
%
%   No legend is drawn here: the caller builds one shared legend from the returned
%   handles, and names the (theta, phi) pair in the panel (e.g. a corner label).

    hold(ax,"on");
    h1 = plot(ax, d.freq*freqScale, d.ads.phase, ...
        "LineWidth", 1.8, "Color", col, "LineStyle","-");
    h2 = plot(ax, d.fekoDS.phase.freq*freqScale, d.fekoDS.phase.data, ...
        "LineStyle","none", "Marker", marker, "MarkerSize", markerSize, ...
        "MarkerEdgeColor", col, "MarkerFaceColor", col);

    grid(ax,"on"); box(ax,"on");
    set(ax, "FontSize", axisLabelFontSize);
    ylim(ax, [-180 180]);
end
