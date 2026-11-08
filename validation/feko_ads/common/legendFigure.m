function fig = legendFigure(name, h, labels, nCols)
%LEGENDFIGURE  One shared legend, in a figure of its own, for traces drawn on several figures.
%   h: line handles whose style the legend shows; labels: one label per handle; nCols: legend columns.

if nargin < 4 || isempty(nCols), nCols = numel(h); end
[fig, ~] = newFigure(name, "none");
ax = axes(fig, 'Visible', 'off');  hold(ax, 'on');
hL = gobjects(1, numel(h));
for k = 1:numel(h)
    hL(k) = plot(ax, NaN, NaN, 'LineStyle', h(k).LineStyle, 'LineWidth', h(k).LineWidth, 'Color', h(k).Color, ...
        'Marker', h(k).Marker, 'MarkerSize', h(k).MarkerSize, 'MarkerFaceColor', h(k).MarkerFaceColor, ...
        'MarkerEdgeColor', h(k).MarkerEdgeColor);
end
legend(ax, hL, cellstr(string(labels)), 'Interpreter', 'latex', 'NumColumns', nCols, 'Location', 'north');
end
