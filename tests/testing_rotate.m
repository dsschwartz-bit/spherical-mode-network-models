%TESTING_ROTATE  Visual check of RotateVectorField.
%   Rotates the tangential field F = theta_hat/2, sampled on a sphere, by
%   alpha about axis k and overlays the original and rotated fields with
%   quiver3.

clear
close all
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

k=[1;1;0];             % rotation axis
alpha=60*(pi/180);     % rotation angle (rad)

numth=10;numph=10;numr=1;
theta=linspace(0,pi,numth);
phi=linspace(0,2*pi,numph);
r=linspace(1,2,numr);

[PHI,THETA,R]=meshgrid(phi,theta,r);

X=R.*sin(THETA).*cos(PHI);
Y=R.*sin(THETA).*sin(PHI);
Z=R.*cos(THETA);

F=zeros(numth,numph,numr,3);
F(:,:,:,1)=cos(THETA).*cos(PHI)/2;
F(:,:,:,2)=cos(THETA).*sin(PHI)/2;
F(:,:,:,3)=-sin(THETA)/2;
figure
quiver3(X,Y,Z,F(:,:,:,1),F(:,:,:,2),F(:,:,:,3))
xlabel("x")
ylabel("y")

[Fprime,Xprime,Yprime,Zprime] = RotateVectorField(k,alpha,X,Y,Z,F);

hold on
quiver3(Xprime,Yprime,Zprime,Fprime(:,:,:,1),Fprime(:,:,:,2),Fprime(:,:,:,3))
hold off
