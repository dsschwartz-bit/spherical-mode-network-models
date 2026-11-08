function dy_hat = sphericalRicatti_y_prime(nu,z)
%SPHERICALRICATTI_Y_PRIME  Derivative of the Riccati-Bessel function z*y_nu(z).
%   dy_hat = sphericalRicatti_y_prime(nu,z) uses
%   d/dz[z y_nu] = z y_{nu-1} - nu y_nu.
dy_hat = sphericalRicatti_y(nu-1,z) - nu.*sphericalBessely(nu,z) ;
end