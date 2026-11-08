function d_H = besselH_prime(nu,x)
%BESSELH_PRIME  Derivative of the Hankel function H_nu^(1)(x) with respect to x.
%   d_H = besselH_prime(nu,x) uses H'_nu = (H_{nu-1} - H_{nu+1})/2 with
%   MATLAB's besselh (first kind).
d_H = 0.5*( besselh(nu-1,x) - besselh(nu+1,x) ) ;
end