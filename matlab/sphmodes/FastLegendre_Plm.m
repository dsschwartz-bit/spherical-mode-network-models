function [p, err]=FastLegendre_Plm(l,m,x)
%FASTLEGENDRE_PLM  Associated Legendre function P_l^m(x) for one (l,m).
%   [p, err] = FastLegendre_Plm(l,m,x) evaluates the unnormalised
%   associated Legendre function, including the Condon-Shortley phase
%   (-1)^m (same convention as MATLAB's legendre(l,x)), at every element of
%   x; p has the shape of x. p = FastLegendre_Plm(l,x) gives the Legendre
%   polynomial P_l(x).
%
%   The polynomial coefficients are computed directly rather than by
%   recursion, which is fast but loses accuracy as l grows because the
%   coefficients become large. err = (largest coefficient)*eps is a
%   worst-case error estimate; a message is printed if it exceeds 0.01.
%   Returns zeros if |m| > |l|, uses P_{-l-1} for l < 0, and (for m ~= 0)
%   NaN where |x| > 1.
if(nargin==2)
    x=m;
    clear m;
    m=0;
end
if(abs(m)>abs(l)), p=zeros(size(x)); err=0; return; end
if(l<0), l=-l-1; end
% For m<0, P_l^m is proportional to P_l^|m|; cfnm is the factor.
cfnm=1;
if(m<0)
    m=-m;
    cfnm=(-1)^m*factorial(l-m)/factorial(l+m);
end
% Coefficient of the highest power of x (explicit formula).
cl=(-1)^m*cfnm*factorial(2*l)/((2^l)*factorial(l)*factorial(l-m));
maxcf=abs(cl);
px=l-m;
% Polynomial part P_l^m(x)/(1-x^2)^(m/2), by Horner's rule in x^2
% (exponents of x step down by 2; px is the current exponent).
x2=x.*x;
p=cl*ones(size(x));
for j=l-1:-1:0
    if(px<2), break; end
    cl=-(j+j+2-l-m)*(j+j+1-l-m)/(2*(j+j+1)*(l-j))*cl;
    if(maxcf<abs(cl)), maxcf=abs(cl); end
    p=p.*x2 + cl;
    px=px-2;
end
err=maxcf*eps;
if err>0.01, fprintf('Coef is %.16f, err %.16f\n',maxcf, err); end
% Odd polynomial: remaining factor x.
if(px==1), p=p.*x; end
% Multiply by (1-x^2)^(m/2).
if(m==0), return; end
x2=1-x2;
for j=1:floor(m/2), p=p.*x2; end
if(m~=2*floor(m/2)), p=p.*sqrt(x2); end
% Not defined for |x|>1.
p(abs(x)>1)=NaN;
