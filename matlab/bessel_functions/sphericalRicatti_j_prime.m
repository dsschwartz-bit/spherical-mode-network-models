function dj_hat = sphericalRicatti_j_prime(nu,z)
%SPHERICALRICATTI_J_PRIME  Derivative of the Riccati-Bessel function z*j_nu(z).
%   dj_hat = sphericalRicatti_j_prime(nu,z) uses
%   d/dz[z j_nu] = z j_{nu-1} - nu j_nu.
dj_hat = sphericalRicatti_j(nu-1,z) - nu.*sphericalBesselj(nu,z) ;
end