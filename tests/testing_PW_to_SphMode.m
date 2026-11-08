%TESTING_PW_TO_SPHMODE  Plane-wave field in local spherical components.
%   Converts a plane-wave field E0 given in (phi, theta, r) components at
%   direction r_pt to Cartesian, then projects that uniform field onto the
%   local (phi, theta, r) unit vectors over a spherical grid and plots the
%   position and field vectors. Uses sph2cartvec (Phased Array System
%   Toolbox).

clear
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

r_pt_deg = [ -30;60 ] % direction [phi; theta] (deg)
r_pt=[deg2rad(r_pt_deg.')].'; %radians
E0=1;
E0_ptr = [ 1 ;-1 ;1 ]*E0;   % [E_phi; E_theta; E_r]

[x,y,z] = sph2cart(r_pt(1,:),pi/2-r_pt(2,:),1);r_xyz=[x;y;z];

E0_xyz = sph2cartvec([E0_ptr(1);-E0_ptr(2);E0_ptr(3)],r_pt_deg(1),90-r_pt_deg(2));

E0_xyz
%%
numth=10;numph=10;numr=1;
theta=linspace(0,pi,numth);
phi=linspace(0,2*pi,numph);
r=linspace(1,2,numr);

[PHI,THETA,R]=meshgrid(phi,theta,r);

X=R.*sin(THETA).*cos(PHI);
Y=R.*sin(THETA).*sin(PHI);
Z=R.*cos(THETA);
E_xyz = zeros(numth,numph,numr,3);E_ptr=E_xyz; r_xyz=E_xyz;
for l=1:numth
    for m=1:numph
        for n=1:numr
            r_xyz(l,m,n,:) = R(n)*[sin(THETA(l,m,n))*cos(PHI(l,m,n));...
                                   sin(THETA(l,m,n))*sin(PHI(l,m,n));...
                                   cos(THETA(l,m,n))]; 
            E_xyz(l,m,n,:)=E0_xyz;
            E_ptr(l,m,n,:)=cart2sphvec_PhiTheta(E0_xyz,PHI(l,m,n),THETA(l,m,n));
        end
    end
end
figure
quiver3(X,Y,Z,r_xyz(:,:,:,1),r_xyz(:,:,:,2),r_xyz(:,:,:,3))
xlabel("x")
ylabel("y")

figure
quiver3(X,Y,Z,E_xyz(:,:,:,1),E_xyz(:,:,:,2),E_xyz(:,:,:,3))
xlabel("x")
ylabel("y")

surf(X,Y,Z,E_ptr(:,:,:,3))
xlabel("x")
ylabel("y")
colorbar
