function Plm = Tesseral_Plm(l,m,THETA)
%TESSERAL_PLM  Normalised associated Legendre function of a tesseral harmonic.
%   Plm = Tesseral_Plm(l,m,THETA) returns
%       sqrt((l+1/2)(l-|m|)!/(l+|m|)!) * P_l^|m|(cos THETA)
%   without the Condon-Shortley phase (the (-1)^m in LgndNorm cancels the
%   one included by FastLegendre_Plm), so that int_0^pi Plm^2 sin(theta)
%   dtheta = 1. This is the theta factor of TesseralHarmonic; Tesseral_Phim
%   is the phi factor. THETA (radians) may be any array; returns 0 if |m| > l.
%
%   See also TesseralHarmonic, Tesseral_Phim, gradYlm.
if abs(m)>l
    Plm=0;
else
    LgndNorm = (-1)^m * sqrt( factorial(l + abs(m)) / factorial(l - abs(m)) / (l + (1/2)) );

    Plm = FastLegendre_Plm(l, abs(m), cos(THETA)) / LgndNorm;
    
end

end