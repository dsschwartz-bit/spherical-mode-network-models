function j_hat = sphericalRicatti_j(nu,z)
%SPHERICALRICATTI_J  Riccati-Bessel function z*j_nu(z).
%   j_hat = sphericalRicatti_j(nu,z) = sqrt(pi*z/2) * J_{nu+1/2}(z).
j_hat = 1.253314137315500*sqrt(z).*besselj(nu+0.5,z);
end