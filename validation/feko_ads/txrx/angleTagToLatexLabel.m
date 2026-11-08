function label = angleTagToLatexLabel(theta, phi, source)
%ANGLETAGTOLATEXLABEL  Legend label from a (theta, phi) observation-angle pair and a source name
%   (e.g. "ADS" or "FEKO"):  "ADS, $(\theta,\phi) = (60^\circ,30^\circ)$".
%   The source is plain text outside the math (MATLAB's LaTeX interpreter drops spaces in math mode).

    label = string(source) + ", $(\theta,\phi) = (" + theta + "^\circ," + phi + "^\circ)$";
end
