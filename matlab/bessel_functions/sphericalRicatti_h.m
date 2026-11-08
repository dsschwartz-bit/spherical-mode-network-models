function h_hat = sphericalRicatti_h(nu,k,z)
%SPHERICALRICATTI_H  Riccati-Hankel function z*h_nu^(k)(z), k = 1 or 2.
%   h_hat = sphericalRicatti_h(nu,k,z) = sqrt(pi*z/2) * H_{nu+1/2}^(k)(z).
h_hat = 1.253314137315500*sqrt(z).*besselh(nu+0.5,k,z);
end