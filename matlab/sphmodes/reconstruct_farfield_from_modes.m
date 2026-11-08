function [Eth, Eph] = reconstruct_farfield_from_modes(psiTM, psiTE, theta, phi, lval, mval, eta0)
%RECONSTRUCT_FARFIELD_FROM_MODES  Far field rE_theta, rE_phi from TM/TE mode coefficients.
%   [Eth, Eph] = reconstruct_farfield_from_modes(psiTM, psiTE, theta, phi, lval, mval, eta0)
%   rebuilds the far field at one frequency by inverting the projection in
%   compute_mode_coeffs_farfield / compute_port_mode_coeffs_farfield:
%       psiTM = 1/sqrt(eta0) * INT [E_ph(-gradY_ph) + E_th(-gradY_th)] dOmega
%       psiTE = sqrt(eta0)   * INT [H_ph(-gradY_ph) + H_th(-gradY_th)] dOmega
%   with H = (rhat x E)/eta0. Since the gradYlm basis is orthonormal and
%   the gradient (TM) and rhat x gradient (TE) parts are orthogonal,
%       E_th = sqrt(eta0) * SUM [ -psiTM gradY_th - psiTE gradY_ph ]
%       E_ph = sqrt(eta0) * SUM [ -psiTM gradY_ph + psiTE gradY_th ]
%   (rhat x grad Y = gradY_th phihat - gradY_ph thetahat). The result has
%   the same units and phase reference as the rE data the coefficients
%   were computed from.
%
%   Inputs:
%     psiTM, psiTE  (L x M) coefficients at one frequency (NaN where |m| > l)
%     theta         (numth x 1) radians;  phi (1 x numph) radians
%     lval, mval    mode index vectors matching psiTM/psiTE (mval = -lmax:lmax)
%     eta0          free-space impedance (ohm)
%   Outputs:
%     Eth, Eph      (numth x numph) complex
%
%   See also gradYlm, compute_mode_coeffs_farfield.

theta = theta(:); phi = phi(:).';
[PHI, THETA] = meshgrid(phi, theta);
have_CSphase = false;   % same convention as the projection

Eth = zeros(numel(theta), numel(phi));
Eph = zeros(numel(theta), numel(phi));

for lindx = 1:numel(lval)
    for mindx = 1:numel(mval)
        if abs(mval(mindx)) > lval(lindx)
            continue
        end
        a = psiTM(lindx, mindx);
        b = psiTE(lindx, mindx);
        if isnan(a) && isnan(b)
            continue
        end
        a(isnan(a)) = 0; b(isnan(b)) = 0;
        [gth, gph] = gradYlm(lval(lindx), mval(mindx), THETA, PHI, have_CSphase);
        Eth = Eth + (-a*gth - b*gph);
        Eph = Eph + (-a*gph + b*gth);
    end
end

Eth = sqrt(eta0) * Eth;
Eph = sqrt(eta0) * Eph;
end
