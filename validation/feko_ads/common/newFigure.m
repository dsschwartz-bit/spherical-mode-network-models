function [fig, ax] = newFigure(name, type)
%NEWFIGURE  Docked figure with one axes (type "rect", default), polar axes ("polar") or none ("none").

if nargin < 2 || isempty(type), type = "rect"; end
fig = figure('Name', name, 'NumberTitle', 'off', 'WindowStyle', 'docked');
switch type
    case "polar", ax = polaraxes(fig);
    case "none",  ax = [];
    otherwise,    ax = axes(fig);
end
if ~isempty(ax), hold(ax, 'on'); end
end
