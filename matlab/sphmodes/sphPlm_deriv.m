function [y, d] = sphPlm_deriv(lmax, x)
%SPHPLM_DERIV  Normalised associated Legendre functions and derivatives.
%   [y, d] = sphPlm_deriv(lmax, x)
% Computes the normalized associated Legendre polynomial 
% \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) 
% and its derivative, suitable for use in spherical harmonics,
% for l=0,...,lmax and m=0,...,l.
%
% Vectorized in x, returns matrix with scaled P_l^m(x) at y(:, l+1,m+1) and
% scaled d/dx P_l^m(x) at d(:, l+1,m+1).
%
% Based on recursion formulas found in:
%   W.H. Press, Numerical Recipes 3rd edition, p. 292
%   T. Limpanuparb, J. Milthorpe, arXiv:1410.1748v1 [physics.chem-ph]
% using the implementation gsl_sf_legendre_sphPlm_deriv_array found in the
% Gnu Scientific Library as a template.
%
% Ludvig af Klinteberg, November 2014 

y = zeros(numel(x), lmax+1, lmax+1);
d = zeros(numel(x), lmax+1, lmax+1);
y(:, 1, 1) = 1./sqrt(4*pi);
d(:, 1, 1) = 0;
x = x(:);
u = sqrt(1-x.*x);
ud = -x./u;
for m=0:lmax
    if m>0
        c = -sqrt(1+1/(2*m));
        y(:,m+1,m+1) = c*u.*y(:,m,m);
        d(:,m+1,m+1) = c*(ud.*y(:,m,m) + u.*d(:,m,m));
    end
    if m<lmax
        c = sqrt(2*m+3);
        y(:,m+2,m+1) = c*x.*y(:,m+1,m+1);
        d(:,m+2,m+1) = c*(y(:,m+1,m+1) + x.*d(:,m+1,m+1));
    end
    for l=m+2:lmax
        a = sqrt((4*l^2-1)/(l^2-m^2));
        b = -sqrt(((l-1)^2-m^2)/(4*(l-1)^2-1));
        y(:,l+1,m+1) = a*(x.*y(:,l,m+1)+b*y(:,l-1,m+1));
        d(:,l+1,m+1) = a*(y(:,l,m+1)+x.*d(:,l,m+1)+b*d(:,l-1,m+1));
    end
end

if numel(x)==1
    y = squeeze(y);
    d = squeeze(d);
end
