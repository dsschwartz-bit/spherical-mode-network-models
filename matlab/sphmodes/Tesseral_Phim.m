function Phim = Tesseral_Phim(m,PHI)
%TESSERAL_PHIM  Orthonormal azimuthal factor of a tesseral harmonic.
%   Phim = Tesseral_Phim(m,PHI) returns, for PHI in radians,
%       1/sqrt(2*pi)          m = 0
%       cos(m*PHI)/sqrt(pi)   m > 0
%       sin(|m|*PHI)/sqrt(pi) m < 0
%   so that int_0^2pi Phim^2 dphi = 1.
%
%   See also Tesseral_Plm, TesseralHarmonic.

if m == 0
    Phim = 1/sqrt(2*pi);
elseif m < 0
    Phim = sin((-m)*PHI)/sqrt(pi);
else  % m > 0
    Phim = cos(m*PHI)/sqrt(pi);
end

end