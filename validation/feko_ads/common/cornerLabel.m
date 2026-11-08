function h = cornerLabel(ax, S, txt, where)
%CORNERLABEL  LaTeX text in a corner of a rectangular axes: where = "ne" (default) | "nw" | "se" | "sw".

if ~S.cornerLabel, h = gobjects(0); return, end
if nargin < 4, where = "ne"; end
switch where
    case "nw", p = [0.02 0.96];  ha = 'left';   va = 'top';
    case "se", p = [0.98 0.04];  ha = 'right';  va = 'bottom';
    case "sw", p = [0.02 0.04];  ha = 'left';   va = 'bottom';
    otherwise, p = [0.98 0.96];  ha = 'right';  va = 'top';
end
h = text(ax, p(1), p(2), txt, 'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', S.labelSize, ...
    'HorizontalAlignment', ha, 'VerticalAlignment', va, 'BackgroundColor', 'w', 'Margin', 1);
end
