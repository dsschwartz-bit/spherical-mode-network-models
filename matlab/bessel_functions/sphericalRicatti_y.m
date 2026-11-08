function y_hat = sphericalRicatti_y(nu,z)
%SPHERICALRICATTI_Y  Riccati-Bessel function z*y_nu(z) (second kind).
%   y_hat = sphericalRicatti_y(nu,z) = sqrt(pi*z/2) * Y_{nu+1/2}(z).
y_hat = 1.253314137315500*sqrt(z).*bessely(nu+0.5,z);
end