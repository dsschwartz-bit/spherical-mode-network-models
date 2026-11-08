function plotSingleThetaPhaseTile(ax, d, col, freqScale, marker, ...
    markerSize, axisLabelFontSize)
%PLOTSINGLETHETAPHASETILE  Plot circuit-model and FEKO phase vs frequency for one angle, with legend.
%   PLOTSINGLETHETAPHASETILE(ax, d, col, freqScale, marker, markerSize,
%   axisLabelFontSize) draws d.circuit.phase as a line and the thinned FEKO
%   phase d.fekoDS.phase as markers, both in colour col, and adds an
%   "ADS"/"FEKO" legend to the axes. Not used by the scripts in this folder;
%   see plotSingleAnglePhaseTile for the dataset layout of
%   buildTransmissionDataset (which names the ADS field d.ads, not d.circuit).

    hold(ax,"on");
    h1 = plot(ax, d.freq*freqScale, d.circuit.phase, ...
        "LineWidth", 1.8, "Color", col, "LineStyle","-");
    h2 = plot(ax, d.fekoDS.phase.freq*freqScale, d.fekoDS.phase.data, ...
        "LineStyle","none", "Marker", marker, "MarkerSize", markerSize, ...
        "MarkerEdgeColor", col, "MarkerFaceColor", col);
    hold(ax,"off");

    grid(ax,"on"); box(ax,"on");
    set(ax, "FontSize", axisLabelFontSize);
    ylim(ax, [-180 180]);
    legend(ax, [h1 h2], ["ADS","FEKO"], "Location","best", ...
        "FontSize", max(axisLabelFontSize-2, 8));
end