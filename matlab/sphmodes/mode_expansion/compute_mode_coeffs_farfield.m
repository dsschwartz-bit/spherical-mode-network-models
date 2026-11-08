function ff = compute_mode_coeffs_farfield(fileprefix, far_request, Freq, Zindx, lval, mval, L, M, use_FF_method)
%COMPUTE_MODE_COEFFS_FARFIELD  TM/TE mode coefficients from a FEKO far-field (.ffe) export.
%   ff = compute_mode_coeffs_farfield(fileprefix, far_request, Freq, Zindx,
%        lval, mval, L, M, use_FF_method)
%   reads <fileprefix>_<far_request>.ffe and, if the grid covers the full
%   sphere, projects the far field onto the gradYlm basis.
%
%   Builds its own (theta,phi) grid and gradYlm basis from the .ffe file,
%   so the far-field grid need not match the near-field grid. The pattern
%   rE has the 1/r dependence removed, so no radius or wavenumber enters:
%       psiTM = 1/sqrt(eta0) * INT [E_ph(-gradY_ph) + E_th(-gradY_th)] dOmega
%       psiTE = sqrt(eta0)   * INT [H_ph(-gradY_ph) + H_th(-gradY_th)] dOmega
%   with H = (rhat x E)/eta0, evaluated by trapezoidal integration.
%
%   fileprefix   : string, path prefix of the FEKO output files
%   far_request  : far-field request name; if empty, the file is missing or
%                  the grid is partial, ff.available is false
%   Freq, Zindx  : (1 x numf) frequencies and the indices into the file's
%                  frequency list selected by the caller
%   lval, mval, L, M : mode-index grids built by the caller from lmax
%   use_FF_method : if false, return immediately with ff.available = false
%
% Returns struct ff:
%   .available     - true if the file was read, covers the full sphere
%                    (theta 0..pi, phi 0..2*pi within 1 deg) and
%                    psiTM/psiTE were computed
%   .psiTM, .psiTE - (L x M x numf) coefficients; NaN if not available
%
% See also gradYlm, reconstruct_farfield_from_modes, ReadFF_spherical.

have_CSphase = false;   % same gradYlm convention as the workflow scripts

numf = numel(Freq);
ff = struct('available', false, 'psiTM', nan(L,M,numf), 'psiTE', nan(L,M,numf));

if ~use_FF_method || strlength(far_request) == 0
    return
end

c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0);

try
    ffe_file = fileprefix + "_" + far_request + ".ffe";
    if ~isfile(ffe_file)
        error('compute_mode_coeffs_farfield:ffe_not_found', 'no .ffe file found: %s', ffe_file);
    end
    [rEpt_FF, ~, theta, phi, ~, ~] = ReadFF_spherical(ffe_file);
    theta = theta(:);    % (numth x 1)
    phi   = phi(:).';    % (1 x numph)
    rEpt_FF = rEpt_FF(:,:,:,Zindx);

    if ~local_full_angular_range(theta, phi)
        fprintf('compute_mode_coeffs_farfield: far-field theta/phi grid does not cover the full angular range — skipping far-field method.\n');
        return
    end
catch ME
    fprintf('compute_mode_coeffs_farfield: could not read far-field data (%s) — skipping far-field method.\n', ME.message);
    return
end

numth = length(theta); numph = length(phi);
[PHI, THETA] = meshgrid(phi, theta);

% rHpt_FF: H_far = (r_hat x E_far)/eta0, from the far-field E-pattern
% (H_phi = E_theta/eta0, H_theta = -E_phi/eta0)
rHpt_FF = zeros(size(rEpt_FF));
rHpt_FF(:,:,1,:) =  rEpt_FF(:,:,2,:) / eta0;   % H_phi   =  E_theta/eta0
rHpt_FF(:,:,2,:) = -rEpt_FF(:,:,1,:) / eta0;   % H_theta = -E_phi/eta0

psiTM = nan(L,M,numf); psiTE = psiTM;
gradYlm_th = zeros(numth,numph,L,M); gradYlm_ph = gradYlm_th;

for lindx = 1:L
    for mindx = 1:M
        if abs(mval(mindx)) > lval(lindx)
            continue
        end

        [gradYlm_th(:,:,lindx,mindx), gradYlm_ph(:,:,lindx,mindx)] = ...
            gradYlm(lval(lindx), mval(mindx), THETA, PHI, have_CSphase);

        for f = 1:numf
            psiTM(lindx,mindx,f) = 1/sqrt(eta0) * ...
                trapz(phi,...
                        trapz( theta, ( rEpt_FF(:,:,1,f).*(-gradYlm_ph(:,:,lindx,mindx)) + ...
                                        rEpt_FF(:,:,2,f).*(-gradYlm_th(:,:,lindx,mindx)) ...
                                       ).*sin(THETA)...
                        ,1) ...
                ,2) ;

            psiTE(lindx,mindx,f) = sqrt(eta0) * ...
                trapz(phi,...
                        trapz( theta, ( rHpt_FF(:,:,1,f).*(-gradYlm_ph(:,:,lindx,mindx)) + ...
                                        rHpt_FF(:,:,2,f).*(-gradYlm_th(:,:,lindx,mindx)) ...
                                       ).*sin(THETA)...
                        ,1) ...
                ,2) ;
        end
    end
end

ff.available = true;
ff.psiTM = psiTM;
ff.psiTE = psiTE;

end

function tf = local_full_angular_range(theta_chk, phi_chk)
% True if theta spans 0..pi and phi spans 0..2*pi (within 1 degree).
tol = 1*pi/180;   % 1 degree tolerance
theta_chk = theta_chk(:);
phi_chk   = phi_chk(:);
tf = numel(theta_chk) > 1 && numel(phi_chk) > 1 && ...
     min(theta_chk) < tol           && max(theta_chk) > (pi   - tol) && ...
     min(phi_chk)   < tol           && max(phi_chk)   > (2*pi - tol);
end
