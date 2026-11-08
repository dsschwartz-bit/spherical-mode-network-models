function d_J = besselJ_prime(nu,x)
%BESSELJ_PRIME  Derivative of the Bessel function J_nu(x) with respect to x.
%   d_J = besselJ_prime(nu,x) uses J'_nu = (J_{nu-1} - J_{nu+1})/2.
d_J = 0.5*( besselj(nu-1,x) - besselj(nu+1,x) ) ;
end