function d_Y = besselY_prime(nu,x)
%BESSELY_PRIME  Derivative of the Bessel function Y_nu(x) with respect to x.
%   d_Y = besselY_prime(nu,x) uses Y'_nu = (Y_{nu-1} - Y_{nu+1})/2.
d_Y = 0.5*( bessely(nu-1,x) - bessely(nu+1,x) ) ;
end