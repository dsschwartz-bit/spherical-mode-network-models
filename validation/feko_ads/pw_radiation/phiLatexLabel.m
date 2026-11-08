function lbl = phiLatexLabel(phiDeg)
%PHILATEXLABEL  LaTeX legend label for a phi cut, e.g. "$\phi = 30^{\circ}$".
%   The ADS sweep variable is called phi_prime; the prime is dropped here
%   so ADS and FEKO labels read identically.

    lbl = string(sprintf('$\\phi = %g^{\\circ}$', phiDeg));
end
