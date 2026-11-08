function result = compute_spherical_mode_coefficients_farfield(fileprefix, lmax, fmin, fmax)
%COMPUTE_SPHERICAL_MODE_COEFFICIENTS_FARFIELD  TM/TE mode coefficients
%from FEKO far-field data.
%   Projects the FEKO far-field pattern onto the angular gradient basis
%   gradYlm (same far-field projection as the far-field cross-check in
%   matlab/workflows/general_vector_spherical_modes_FEKO.m) to obtain psiTM
%   and psiTE.
%
%   result = COMPUTE_SPHERICAL_MODE_COEFFICIENTS_FARFIELD(fileprefix, lmax)
%   result = COMPUTE_SPHERICAL_MODE_COEFFICIENTS_FARFIELD(fileprefix, lmax, fmin, fmax)
%
%   fileprefix : string, path prefix of the FEKO output files, e.g.
%                fullfile(feko_dir, "cylinder_rotated_yaxis_45deg");
%                <prefix>_FarField1.ffe and <prefix>.out must exist
%   lmax       : maximum spherical mode order to compute (l = 1..lmax;
%                there is no l=0 TE/TM mode)
%   fmin, fmax : optional frequency range in Hz to restrict computation
%                to (applied after the negative-real-impedance Zindx
%                filtering). If omitted or empty, the full frequency
%                range found in the imported files is used.
%
%   Returns a struct with fields:
%     psiTM, psiTE : (L x M x numf) complex coefficient arrays.
%                    Row lindx <-> l = lval(lindx) = lindx (l = 1..lmax).
%                    Column mindx <-> m = mval(mindx), mval = -lmax:lmax,
%                    FIXED across all rows (centered convention: m=0 is
%                    always at column lmax+1). For a given row l, only
%                    columns with |m|<=l are valid; others are NaN.
%     lval, mval   : as above
%     L, M         : size(psiTM,1), size(psiTM,2)
%     Freq         : (1 x numf) frequency vector in Hz, after removing
%                    points with negative real port impedance (Zindx) and
%                    applying fmin/fmax
%     fileprefix   : echoed back
%
%   The far-field pattern has the 1/r dependence removed, so no extraction
%   radius or wavenumber normalisation is needed; gradYlm is built on the
%   far-field file's own (theta,phi) grid. The projection integral needs
%   the full sphere, so a single-direction or partial angular sweep raises
%   an error.

if nargin < 3 || isempty(fmin)
    fmin = -Inf;
end
if nargin < 4 || isempty(fmax)
    fmax = Inf;
end

fileprefix = string(fileprefix);

%% --- Physical Constants ---
c0   = physconst('lightspeed'); %#ok<NASGU>
mu0  = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0);

%% --- FEKO Import: Read Port Parameters from .out File ---
[~, Z0, ~, ~, ~, ~, Zant, Freq, ~] = ReadFEKOout(fileprefix+".out");

Zant = transpose(Zant(:,1)); Freq = transpose(Freq);
Zindx = find(~(real(Zant)<0));
Freq  = Freq(Zindx);

if abs(Z0-eta0) > 0.01
    error("antenna port impedance should be eta0 for now")
end

%% --- Import FEKO Far Field Data ---
[rEpt_FF, ~, theta, phi, ~, ~] = ReadFF_spherical(fileprefix+"_FarField1.ffe");
theta = theta(:);    % (numth x 1)
phi   = phi(:).';    % (1 x numph)
rEpt_FF = rEpt_FF(:,:,:,Zindx);

if ~local_full_angular_range(theta, phi)
    error('compute_spherical_mode_coefficients_farfield:partialSweep', ...
          ['Far-field theta/phi grid does not cover the full angular range ' ...
           '(theta: 0 to pi, phi: 0 to 2*pi) -- a single-direction gain query ' ...
           'or a truncated sweep cannot be used for mode-coefficient extraction.'])
end

%% --- Restrict to requested frequency range ---
frange_idx = find(Freq >= fmin & Freq <= fmax);
if isempty(frange_idx)
    error('compute_spherical_mode_coefficients_farfield:noFreqInRange', ...
          'No frequency points found in [%.4g, %.4g] Hz (available range: [%.4g, %.4g] Hz)', ...
          fmin, fmax, min(Freq), max(Freq));
end
Freq    = Freq(frange_idx);
rEpt_FF = rEpt_FF(:,:,:,frange_idx);

numf  = length(Freq);
numth = length(theta); numph = length(phi); %#ok<NASGU>

%% --- Mode index setup ---
lval = (1:1:lmax).';   % mode orders l = 1..lmax
mval = (-lmax:lmax);   % m ranges from -lmax to +lmax
L    = length(lval);
M    = length(mval);

psiTM = nan(L,M,numf); psiTE = psiTM;

have_CSphase = false;   % must match the convention used by the workflow scripts

[PHI, THETA] = meshgrid(phi, theta);

gradYlm_th = zeros(numth,numph,L,M); gradYlm_ph = gradYlm_th;

% rHpt_FF: H_far = (r_hat x E_far)/eta0, from the far-field E-pattern
% (H_phi = E_theta/eta0, H_theta = -E_phi/eta0)
rHpt_FF = zeros(size(rEpt_FF));
rHpt_FF(:,:,1,:) =  rEpt_FF(:,:,2,:) / eta0;   % H_phi   =  E_theta/eta0
rHpt_FF(:,:,2,:) = -rEpt_FF(:,:,1,:) / eta0;   % H_theta = -E_phi/eta0

%% --- Determine psiTM, psiTE via far-field projection ---
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

%% --- Package result ---
result.psiTM = psiTM;
result.psiTE = psiTE;
result.lval  = lval;
result.mval  = mval;
result.L     = L;
result.M     = M;
result.Freq  = Freq;
result.fileprefix = fileprefix;

end

function tf = local_full_angular_range(theta_chk, phi_chk)
% LOCAL_FULL_ANGULAR_RANGE  True if theta_chk spans [0, pi] and phi_chk
% spans [0, 2*pi] (full angular sweep) to within 1 degree.
tol = 1*pi/180;   % 1 degree tolerance
theta_chk = theta_chk(:);
phi_chk   = phi_chk(:);
tf = numel(theta_chk) > 1 && numel(phi_chk) > 1 && ...
     min(theta_chk) < tol           && max(theta_chk) > (pi   - tol) && ...
     min(phi_chk)   < tol           && max(phi_chk)   > (2*pi - tol);
end
