function [h] = sphericalHankel(nu,k,z)
%SPHERICALHANKEL  Spherical Hankel function h_nu^(k)(z), k = 1 or 2.
%   h = sphericalHankel(nu,k,z) = sqrt(pi/(2z)) * H_{nu+1/2}^(k)(z)
%   (1.2533141373155 = sqrt(pi/2)). With the exp(+jwt) time convention,
%   k = 2 is an outgoing wave and k = 1 an incoming wave.
h=1.253314137315500.*besselh(nu+0.5,k,z) ./ sqrt(z);
end