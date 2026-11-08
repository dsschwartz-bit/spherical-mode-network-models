function dh_hat = sphericalRicatti_h_prime(nu,k,z)
%SPHERICALRICATTI_H_PRIME  Derivative of the Riccati-Hankel function z*h_nu^(k)(z).
%   dh_hat = sphericalRicatti_h_prime(nu,k,z) uses
%   d/dz[z h_nu] = z h_{nu-1} - nu h_nu.
dh_hat = sphericalRicatti_h(nu-1,k,z) - nu.*sphericalHankel(nu,k,z) ;
end