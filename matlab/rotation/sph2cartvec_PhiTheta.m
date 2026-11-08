function F_xyz = sph2cartvec_PhiTheta(F_ptr,phi,theta)
%SPH2CARTVEC_PHITHETA  [phi; theta; r] components -> Cartesian vector.
%   F_xyz = sph2cartvec_PhiTheta(F_ptr,phi,theta) converts
%   F_ptr = [F_phi F_theta F_r] at direction (theta,phi) (radians, theta
%   measured from +z) to Cartesian [Fx; Fy; Fz]. Inverse of
%   cart2sphvec_PhiTheta; wraps MATLAB's sph2cartvec.
F_ptr = reshape(F_ptr,3,1);

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
if flip, F_ptr(1:2)=-F_ptr(1:2); end

F_xyz = sph2cartvec([F_ptr(1);-F_ptr(2);F_ptr(3)],phi*180/pi,90-(theta*180/pi));

end