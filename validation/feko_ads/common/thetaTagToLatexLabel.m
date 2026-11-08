function label = thetaTagToLatexLabel(theta, source)
%THETATAGTOLATEXLABEL  LaTeX legend label from a theta angle and a source name.
%   label = THETATAGTOLATEXLABEL(theta, source) returns e.g.
%   '$\theta=30^\circ \mathrm{Circuit}$' for theta = "30", source = "Circuit".
%   The format string is single-quoted so sprintf resolves the backslashes.

    label = sprintf('$\\theta=%s^\\circ \\mathrm{%s}$', theta, source);
end
