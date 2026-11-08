function [gradYlm_th,gradYlm_ph] = gradYlm(l,m,THETA,PHI,include_CSphase)
%GRADYLM  Angular gradient of a real (tesseral) spherical harmonic.
%   [gradYlm_th,gradYlm_ph] = gradYlm(l,m,THETA,PHI,include_CSphase)
%   returns the theta and phi components of
%       r * grad( Y_lm(theta,phi) ) / sqrt(l(l+1))
%   i.e.  dY/dtheta / sqrt(l(l+1))  and  (1/sin theta) dY/dphi / sqrt(l(l+1)),
%   where Y_lm = Phi_m(phi) * Pbar_l^|m|(cos theta) is the orthonormal
%   tesseral harmonic of TesseralHarmonic (Phi_m = 1/sqrt(2pi), cos(m phi)/
%   sqrt(pi) for m > 0, sin(|m| phi)/sqrt(pi) for m < 0). With the
%   1/sqrt(l(l+1)) factor the gradients are orthonormal over the sphere.
%
%   THETA, PHI (radians) are arrays of the same size, e.g. from
%   [PHI,THETA] = meshgrid(phi,theta). Two vectors of different
%   size/orientation (theta Nx1, phi 1xM) are expanded to an NxM grid.
%   Exact pole points (THETA == 0 or pi) use the closed-form limits.
%
%   Condon-Shortley phase: FastLegendre_Plm returns P_l^m with the
%   (-1)^m phase included. That phase is always cancelled, so the result
%   is the phase-free orthonormal harmonic whatever the setting:
%   include_CSphase = false (default) folds (-1)^m into the normalisation
%   LgndNorm, and true applies it as a separate factor CSphase. Either way
%   Plm is the same Pbar_l^|m|(cos theta), and the pole limits use
%   CSphase/LgndNorm = (-1)^m/N in both cases. All callers in this toolbox
%   use false.
%
%   See also TesseralHarmonic, legendre_derivative_modified, Spherical_Mode_Field.

if isvector(THETA) && isvector(PHI) && ~isequal(size(THETA), size(PHI))
    [PHI, THETA] = meshgrid(PHI(:).', THETA(:));
end

assert(isequal(size(THETA), size(PHI)), 'gradYlm:ShapeMismatch', ...
    ['THETA and PHI must be the same size, or be two vectors that ' ...
     'can be expanded via meshgrid (e.g. theta as an Nx1 column and ' ...
     'phi as a 1xM row).']);

l_factor = l * (l + 1);
LgndNorm = sqrt( 2 * factorial(l + abs(m)) / factorial(l - abs(m)) / (2*l + 1) );

% Cancel the Condon-Shortley phase of FastLegendre_Plm: through LgndNorm
% (false) or through CSphase (true). Plm is the same phase-free
% orthonormal Pbar_l^|m|(cos theta) for both settings.
if nargin < 5 || ~include_CSphase
    CSphase  = 1;
    LgndNorm = ((-1)^m) * LgndNorm;
else
    CSphase = (-1)^m;
end

Plm = CSphase * FastLegendre_Plm(l, abs(m), cos(THETA))  / LgndNorm;

% dPbar/dx at x = cos(theta). legendre_derivative_modified's 'norm' Pbar
% (and the adjacent order it builds) is phase-free like Plm, so Plm is
% passed as is for both settings.
dPlm = legendre_derivative_modified(l, abs(m), cos(THETA), Plm, 'norm');

% d/dtheta = -sin(theta) d/dx; the phi component carries 1/sin(theta)
% from the spherical gradient. Phi_m normalisation (1/sqrt(2pi) or
% 1/sqrt(pi)) is grouped with the azimuthal term, 1/sqrt(l(l+1)) applied
% separately. Pole points are overwritten below.

if m == 0
    gradYlm_th = -sin(THETA) .* dPlm / sqrt(2*pi) / sqrt(l_factor);
    gradYlm_ph = zeros(size(THETA));
elseif m < 0
    gradYlm_th = -sin(THETA) .* dPlm .* (sin((-m)*PHI)/sqrt(pi)) / sqrt(l_factor);
    gradYlm_ph = (-m) * ((cos((-m)*PHI)/sqrt(pi)) ./ sin(THETA)) .* Plm / sqrt(l_factor);
else  % m > 0
    gradYlm_th = -sin(THETA) .* dPlm .* (cos(m*PHI)/sqrt(pi)) / sqrt(l_factor);
    gradYlm_ph = -m * ((sin(m*PHI)/sqrt(pi)) ./ sin(THETA)) .* Plm / sqrt(l_factor);
end

% Closed-form limits at theta = 0 and pi, where the expressions above
% are 0/0. Only |m| = 1 is non-zero there.

theta0_mask  = (THETA == 0);
thetapi_mask = (THETA == pi);

if m == -1
    gradYlm_th(theta0_mask)  = -(1/2) * l_factor           * (CSphase/LgndNorm) * (sin((-m)*PHI(theta0_mask))/sqrt(pi))  / sqrt(l_factor);
    gradYlm_th(thetapi_mask) = -(1/2) * l_factor * (-1)^l   * (CSphase/LgndNorm) * (sin((-m)*PHI(thetapi_mask))/sqrt(pi)) / sqrt(l_factor);

    gradYlm_ph(theta0_mask)  = -(1/2) * l_factor           * (CSphase/LgndNorm) * (cos(PHI(theta0_mask))/sqrt(pi))  / sqrt(l_factor);
    gradYlm_ph(thetapi_mask) =  (1/2) * l_factor * (-1)^l   * (CSphase/LgndNorm) * (cos(PHI(thetapi_mask))/sqrt(pi)) / sqrt(l_factor);

elseif m == 1
    gradYlm_th(theta0_mask)  = -(1/2) * l_factor           * (CSphase/LgndNorm) * (cos(m*PHI(theta0_mask))/sqrt(pi))  / sqrt(l_factor);
    gradYlm_th(thetapi_mask) = -(1/2) * l_factor * (-1)^l   * (CSphase/LgndNorm) * (cos(m*PHI(thetapi_mask))/sqrt(pi)) / sqrt(l_factor);

    gradYlm_ph(theta0_mask)  =  (1/2) * l_factor           * (CSphase/LgndNorm) * (sin(PHI(theta0_mask))/sqrt(pi))  / sqrt(l_factor);
    gradYlm_ph(thetapi_mask) = -(1/2) * l_factor * (-1)^l   * (CSphase/LgndNorm) * (sin(PHI(thetapi_mask))/sqrt(pi)) / sqrt(l_factor);

else
    gradYlm_th(theta0_mask | thetapi_mask) = 0;
    gradYlm_ph(theta0_mask | thetapi_mask) = 0;
end

end