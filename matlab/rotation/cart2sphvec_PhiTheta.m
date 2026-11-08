function F_ptr = cart2sphvec_PhiTheta(F_xyz,phi,theta)
%CART2SPHVEC_PHITHETA  Cartesian vector -> [phi; theta; r] components.
%   F_ptr = cart2sphvec_PhiTheta(F_xyz,phi,theta) converts the 3-element
%   vector F_xyz = [Fx Fy Fz] at direction (theta,phi) (radians, theta
%   measured from +z) to F_ptr = [F_phi; F_theta; F_r]. Wraps MATLAB's
%   cart2sphvec (azimuth/elevation in degrees); the theta component is
%   negated because elevation increases in the -theta direction.
F_xyz = reshape(F_xyz,3,1);

% Elevation = 90 deg - theta, azimuth = phi wrapped to (-180,180] deg.
% theta in [0,pi] is used as is (theta = pi is the -z pole). Outside it,
% e.g. theta < 0 -> (-theta, phi+pi) and theta in (pi,2pi) -> (2pi-theta,
% phi+pi): the same direction, where the theta and phi unit vectors point
% the other way, so those two components change sign.
theta=mod(theta,2*pi);
flip=theta>pi;
if flip
    theta=2*pi-theta;
    phi=phi+pi;
end
phi=mod(phi,2*pi);
if phi>pi
    phi=phi-(2*pi) ;
end

F_ptr = cart2sphvec(F_xyz,phi*180/pi,90-(theta*180/pi));
F_ptr(2)=-F_ptr(2);
if flip, F_ptr(1:2)=-F_ptr(1:2); end
end