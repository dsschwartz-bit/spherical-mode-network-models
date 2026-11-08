%ROTATE_PW  Polarisation angle of rotated plane waves.
%   Rotates a set of plane waves (directions r_pt = [phi; theta], fields
%   E_ptr in [phi; theta; r] components) by 45 deg about the y-axis and,
%   for each rotated wave, computes the angle nu between its polarisation
%   and the default FEKO plane-wave polarisation E0 = -theta_hat. nu_sph
%   (spherical components) and nu_cart (Cartesian) should agree.

clear
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

r_pt_deg = [ 90,10 ; 90,45 ; 90,90 ].' % columns are [phi; theta] (deg)
r_pt=[deg2rad(r_pt_deg.')].'; %radians
E0=1;
E_ptr = [ 0 ,-1 ,0 ;...
          0 ,-1 ,0 ;...
          0 ,-1 ,0 ;...
    ].'*E0;
num_PWs=size(E_ptr,2);

[x,y,z] = sph2cart(r_pt(1,:),pi/2-r_pt(2,:),1);r_xyz=[x;y;z];
E_xyz=zeros(size(E_ptr));
for n=1:num_PWs
    E_xyz(:,n) = sph2cartvec_PhiTheta(E_ptr(:,n),r_pt(1,n),r_pt(2,n));
end
E_xyz
Axis_xyz = [0;1;0];
angle = 45*(pi/180);
Rot=rotationmat3D(angle,Axis_xyz);
E_xyz_rot=Rot*E_xyz;
r_xyz_rot=Rot*r_xyz;

[az,el,radius]=cart2sph(r_xyz_rot(1,:),r_xyz_rot(2,:),r_xyz_rot(3,:));
r_pt_rot=[az;pi/2-el];
r_pt_rot_deg=rad2deg(r_pt_rot)
E_ptr_rot=zeros(size(E_ptr));
for n=1:num_PWs
    E_ptr_rot(:,n) = cart2sphvec_PhiTheta(E_xyz_rot(:,n),r_pt_rot(1,n),r_pt_rot(2,n)); 
end
E_ptr_rot
% Default FEKO plane-wave polarisation (-theta_hat); nu is the polarisation
% angle about k needed to obtain E_ptr_rot from it.
E0_ptr_rot = [0;-1;0];
E0_xyz_rot = zeros(3,num_PWs);
nu_sph=zeros(1,num_PWs);nu_cart=nu_sph;
s=zeros(1,num_PWs);

for n=1:num_PWs
    E0_xyz_rot(:,n) = sph2cartvec([E0_ptr_rot(1);-E0_ptr_rot(2);E0_ptr_rot(3)],r_pt_rot_deg(1,n),90-r_pt_rot_deg(2,n));
    s(n)=sign( dot( cross( E_xyz_rot(:,n),E0_xyz_rot(:,n)) , r_xyz(:,n)  ) ); 
    nu_sph(n) = atan2d( s(n)*sqrt(sum(cross(E0_ptr_rot,E_ptr_rot(:,n)).^2)) , dot(E0_ptr_rot,E_ptr_rot(:,n)) ) ;
    nu_cart(n) = atan2d( s(n)*sqrt(sum(cross(E0_xyz_rot(:,n),E_xyz_rot(:,n)).^2)) , dot(E0_xyz_rot(:,n),E_xyz_rot(:,n)) ) ;
end
% nu_sph and nu_cart should be identical
nu_sph
nu_cart