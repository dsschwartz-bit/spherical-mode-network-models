function styleAxes(ax, S, xlab, ylab)
%STYLEAXES  Grid, font size and LaTeX axis labels; applies styleTraces to the lines on the axes.

for a = ax(:).'
    set(a, 'FontSize', S.fontSize, 'TickLabelInterpreter', 'latex');
    grid(a, 'on');
    if ~isa(a, 'matlab.graphics.axis.PolarAxes'), box(a, 'on'); end
    styleTraces(findobj(a, 'Type', 'line'), S);
end
if nargin >= 3 && ~isempty(xlab), xlabel(ax(1), xlab, 'Interpreter', 'latex', 'FontSize', S.labelSize); end
if nargin >= 4 && ~isempty(ylab), ylabel(ax(1), ylab, 'Interpreter', 'latex', 'FontSize', S.labelSize); end
end
