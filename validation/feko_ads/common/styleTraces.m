function styleTraces(h, S)
%STYLETRACES  Line width for lines; filled markers (darker edge) for marker traces.
%   Lines with Tag 'keepWidth' keep their own width.

for k = 1:numel(h)
    if ~isgraphics(h(k)) || ~isprop(h(k), 'LineStyle') || strcmp(h(k).Tag, 'keepWidth'), continue, end
    if ~strcmp(h(k).LineStyle, 'none'), h(k).LineWidth = S.lineWidth; end
    if ~strcmp(h(k).Marker, 'none')
        fc = h(k).MarkerFaceColor;
        if ischar(fc) || isstring(fc), fc = h(k).Color; end   % 'none' / 'auto': fill with the line colour
        h(k).MarkerSize = S.markerSize;
        h(k).MarkerFaceColor = fc;
        h(k).MarkerEdgeColor = 0.6*fc;
        if strcmp(h(k).LineStyle, 'none'), h(k).LineWidth = S.markerEdge; end
    end
end
end
