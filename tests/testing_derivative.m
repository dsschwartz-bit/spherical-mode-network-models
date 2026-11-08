%TESTING_DERIVATIVE  Compare associated Legendre values and derivatives.
%   For one (l,m), plots MATLAB's normalised legendre() against
%   FastLegendre_Plm, and legendre_derivative (third_party) against
%   legendre_derivative_modified, versus x = cos(theta).

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

numth=101;
th=linspace(pi/10,pi-pi/10,numth).';
l=4;lmax=max(l);
m = -1;             % azimuthal order, -lmax <= m <= lmax
L    = length(l);
M=length(m);
% Normalisation factors (tesseral harmonics depend on |m|)
l_factor=l*(l+1);
LgndNorm = sqrt( 2* factorial(l+abs(m)) ./ factorial(l-abs(m)) ./ (2*l+1) );
PhiNorm2 = (2.^(m==0))*pi; 
use_cs_phase = false;
if use_cs_phase
    CSphase=(-1).^m;   % Condon-Shortley phase
else
    CSphase=ones(size(m));  % no Condon-Shortley phase; (-1)^m folded into LgndNorm
    LgndNorm = ((-1).^m).*LgndNorm ;
end
LgndPowNorm = LgndNorm.*sqrt(l_factor);

Plm=CSphase.*legendre(l,cos(th),'norm');
Plm2=FastLegendre_Plm(l,abs(m),cos(th))/LgndNorm;
dPlm=zeros(size(Plm));
dPlm2=legendre_derivative_modified(l,abs(m),cos(th),Plm2,'norm');
for t=1:length(th)
    dPlm(:,t)=legendre_derivative(l,cos(th(t)),'norm');
end

plot(cos(th),Plm(1+abs(m),:),"linewidth",1)
hold on
plot(cos(th),Plm2,"--")
plot(cos(th),dPlm(1+abs(m),:))
plot(cos(th),dPlm2,"--")
hold off
grid