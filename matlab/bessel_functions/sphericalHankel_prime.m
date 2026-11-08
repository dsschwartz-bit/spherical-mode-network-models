function dh = sphericalHankel_prime(nu,k,z)
%SPHERICALHANKEL_PRIME  Derivative of h_nu^(k)(z) with respect to z.
%   dh = sphericalHankel_prime(nu,k,z) uses h'_nu = h_{nu-1} - (nu+1)/z * h_nu.
dh = sphericalHankel(nu-1,k,z) - ((nu+1)./z).*sphericalHankel(nu,k,z);
end