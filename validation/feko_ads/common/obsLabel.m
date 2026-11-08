function h = obsLabel(ax, S, theta, phi, where)
%OBSLABEL  Observation point "theta = .., phi = .." (degrees) in a corner of a rectangular axes.

if nargin < 5, where = "ne"; end
txt = sprintf('$\\theta = %s^\\circ$, $\\phi = %s^\\circ$', string(theta), string(phi));
h = cornerLabel(ax, S, txt, where);
end
