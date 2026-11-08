function Pnu_m = Legendre_Pnu_m(nu,m,x)
%LEGENDRE_PNU_M  Associated Legendre function P_nu^m(x) of non-integer degree nu.
%   Pnu_m = Legendre_Pnu_m(nu,m,x) evaluates the Ferrers function of the
%   first kind for real (non-integer) degree nu and integer order m >= 0,
%   on -1 < x <= 1, from the hypergeometric representation
%       P_nu^0(x) = 2F1(nu+1, -nu; 1; (1-x)/2)
%       P_nu^m(x) = (-1)^m Gamma(nu+m+1)/(Gamma(nu-m+1) m!)
%                   * ((1-x)/(1+x))^(m/2) * 2F1(nu+1, -nu; m+1; (1-x)/2)
%   The (-1)^m is the Condon-Shortley phase, so for integer nu this equals
%   MATLAB's legendre(nu,x) row m+1. P_nu(x) is singular at x = -1 for
%   non-integer nu.
%
%   nu may be symbolic (hypergeom is used), which is how coneHarmonic is
%   solved with vpasolve for the biconical eigenvalues nu_n. x is
%   elementwise; m must be an integer (error otherwise).
%
%   See also coneHarmonic.
if ~mod(m,1)  
    if m==0
        Pnu_m = hypergeom([nu+1,-nu],m+1,(1-x)/2) ;
    else
        Pnu_m = (-1)^m * (gamma(nu+m+1)/gamma(nu-m+1))* ((1-x)./(1+x)).^(m/2) ...
            .* hypergeom([nu+1,-nu],m+1,(1-x)/2) / gamma(m+1) ;
    end
else
    error("this function is for integer m values")
end
end
