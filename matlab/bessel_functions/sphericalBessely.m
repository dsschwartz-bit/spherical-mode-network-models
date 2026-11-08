function ys = sphericalBessely(nu,x)
%SPHERICALBESSELY  Spherical Bessel function of the second kind y_nu(x).
%   ys = sphericalBessely(nu,x) = (-1)^(nu+1) * j_{-nu-1}(x), for integer nu.
ys = ((-1).^(nu+1)).*sphericalBesselj(-(nu+1),x);
end
