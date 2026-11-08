%TESTING_ADDITION_THEOREM  Plane-wave expansion in spherical harmonics.
%   Evaluates the plane wave exp(j*k.r) for direction inc_pt on a
%   Cartesian grid of kr, and compares it with the truncated expansion
%     sum_lm 4*pi * j^l * j_l(kr) * Y_lm(rhat') * Y_lm(rhat)
%   using real (tesseral) harmonics up to lmax. Plots phase and magnitude
%   in the z = 0 plane.

clear
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

%%
inc_pt_deg = [ 120;30 ] % plane-wave direction [phi; theta] (deg)
inc_pt=[deg2rad(inc_pt_deg.')].'; %radians

numr=1;
kr=linspace(10,10,numr);   % grid extent in kr

numx=50;numy=50;numz=3;

x=linspace(-max(kr),max(kr),numx);
y=linspace(-max(kr),max(kr),numy);
z=linspace(-max(kr),max(kr),numz);
[X,Y,Z]=meshgrid(x,y,z);

KR=sqrt(X.^2+Y.^2+Z.^2);
THETA = acos(Z./KR);
PHI = mod(atan2(Y,X),2*pi);

kr_xyz=zeros(numx,numy,numz,3);rhatprime_xyz=kr_xyz;
for l=1:numx
    for m=1:numy
        for n=1:numz
            kr_xyz(l,m,n,:) = [X(l,m,n);Y(l,m,n);Z(l,m,n)];
            rhatprime_xyz(l,m,n,:) = sph2cartvec_PhiTheta( [0;0;1], inc_pt(1), inc_pt(2) );
        end
    end
end
phase=dot(kr_xyz,rhatprime_xyz,4);

[~,indexZ]=min(abs(z-0));

PW=exp(1j*phase);
surf(X(:,:,indexZ),Y(:,:,indexZ),angle(PW(:,:,indexZ)))
xlabel("x")
ylabel("y")
colorbar
view([0 0 1])
%%
lmax = 12;
lval = (0:1:lmax).';   % mode orders l = 0..lmax
mval = (-lmax:lmax);   % m ranges from -lmax to +lmax
L=length(lval);M=length(mval);

Clm = zeros(L,M);
PW_modes=zeros(numx,numy,numz,L,M);
for lindx=1:L
    for mindx=1:M
        Clm(lindx,mindx) = 4*pi*TesseralHarmonic(lval(lindx),mval(mindx),inc_pt(2),inc_pt(1));
        PW_modes(:,:,:,lindx,mindx) = Clm(lindx,mindx)*(1j^lval(lindx))*sphericalBesselj(lval(lindx),KR).*TesseralHarmonic(lval(lindx),mval(mindx),THETA,PHI);
    end
end
PW_total=sum(PW_modes,[4,5],"omitnan");
surf(X(:,:,indexZ),Y(:,:,indexZ),angle(PW_total(:,:,indexZ)))
xlabel("x")
ylabel("y")
colorbar
view([0 0 1])

PW_total=sum(PW_modes,[4,5],"omitnan");
surf(X(:,:,indexZ),Y(:,:,indexZ),abs(PW_total(:,:,indexZ)))
xlabel("x")
ylabel("y")
colorbar
view([0 0 1])