function Ylm = TesseralHarmonic(l,m,THETA,PHI)
%TESSERALHARMONIC  Real-valued (tesseral) spherical harmonic Y_lm(theta,phi).
%   Ylm = TesseralHarmonic(l,m,THETA,PHI) returns
%       Ylm = Phi_m(phi) * Pbar_l^|m|(cos theta)
%   with Phi_m = 1/sqrt(2pi) (m=0), cos(m phi)/sqrt(pi) (m>0),
%   sin(|m| phi)/sqrt(pi) (m<0) and Pbar as in Tesseral_Plm (no
%   Condon-Shortley phase). Ylm is orthonormal over the unit sphere.
%
%   THETA, PHI (radians) are arrays of the same size, e.g. from
%   [PHI,THETA] = meshgrid(phi,theta). Two vectors of different
%   size/orientation (theta Nx1, phi 1xM) are expanded to an NxM grid.
%   Returns zeros if |m| > l.
%
%   See also gradYlm, Tesseral_Plm, Tesseral_Phim.

if isvector(THETA) && isvector(PHI) && ~isequal(size(THETA), size(PHI))
    [PHI, THETA] = meshgrid(PHI(:).', THETA(:));
end

assert(isequal(size(THETA), size(PHI)), 'gradYlm:ShapeMismatch', ...
    ['THETA and PHI must be the same size, or be two vectors that ' ...
     'can be expanded via meshgrid (e.g. theta as an Nx1 column and ' ...
     'phi as a 1xM row).']);
if l<0
    error("expected l to be nonegative")
elseif abs(m)>l
    Ylm=zeros(size(THETA));
else
    % (-1)^m cancels the Condon-Shortley phase included by FastLegendre_Plm.
    LgndNorm = (-1)^m * sqrt( factorial(l + abs(m)) / factorial(l - abs(m)) / (l + (1/2)) );
    Plm = FastLegendre_Plm(l, abs(m), cos(THETA)) / LgndNorm;

    if m == 0
        Ylm = (1/sqrt(2*pi)) .* Plm;
    elseif m < 0
        Ylm = (sin((-m)*PHI)/sqrt(pi)) .* Plm;
    else  % m > 0
        Ylm = (cos(m*PHI)/sqrt(pi)) .* Plm;
    end
end

end