function [legH, legLabels] = plotTransmissionTile(ax, dataset, field, ...
    colors, freqScale, marker, markerSize, axisLabelFontSize)
%PLOTTRANSMISSIONTILE  Plot ADS transmission lines and matching
%   downsampled FEKO markers.
%   [legH, legLabels] = PLOTTRANSMISSIONTILE(ax, dataset, field, colors,
%   freqScale, marker, markerSize, axisLabelFontSize) draws, onto axes
%   ax, one solid line per element of the dataset struct array (dense
%   ADS transmission data) plus a marker-only curve (the
%   downsampled FEKO transmission data), using the same color per
%   dataset element.
%
%   dataset(k) fields (see buildTransmissionDataset):
%     theta                    - theta angle (string/number) for the legend
%     phi                      - phi angle (string/number) for the legend
%     freq                     - dense frequency vector [Hz] (ADS grid)
%     ads.(mag|phase)          - dense ADS transmission data
%     fekoDS.(mag|phase).freq  - downsampled FEKO frequency vector [Hz]
%     fekoDS.(mag|phase).data  - downsampled FEKO transmission data
%
%   field - "mag" or "phase", selecting which quantity to plot
%
%   No legend is created inside this function -- it returns the plotted
%   handles (legH) and LaTeX labels (legLabels) so a single shared
%   legend can be built once by the caller.

    hold(ax,"on");
    legH = gobjects(0); legLabels = {};

    for k = 1:numel(dataset)
        d = dataset(k);
        col = colors(k,:);

        h1 = plot(ax, d.freq*freqScale, d.ads.(field), ...
            "LineWidth", 1.8, "Color", col, "LineStyle","-");
        legH(end+1) = h1; %#ok<AGROW>
        legLabels{end+1} = angleTagToLatexLabel(d.theta, d.phi, "ADS"); %#ok<AGROW>

        h2 = plot(ax, d.fekoDS.(field).freq*freqScale, d.fekoDS.(field).data, ...
            "LineStyle","none", "Marker", marker, "MarkerSize", markerSize, ...
            "MarkerEdgeColor", col, "MarkerFaceColor", col);
        legH(end+1) = h2; %#ok<AGROW>
        legLabels{end+1} = angleTagToLatexLabel(d.theta, d.phi, "FEKO"); %#ok<AGROW>
    end

    hold(ax,"off");
    grid(ax,"on"); box(ax,"on");
    set(ax, "FontSize", axisLabelFontSize);
end