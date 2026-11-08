function js = sphericalBesselj(nu,x)
%SPHERICALBESSELJ  Spherical Bessel function of the first kind j_nu(x).
%   js = sphericalBesselj(nu,x) = sqrt(pi/(2x)) * J_{nu+1/2}(x)
%   (1.2533141373155 = sqrt(pi/2)).
js = 1.253314137315500*besselj(nu+0.5,x)./sqrt(x);
end