function result = compute_spherical_mode_coefficients_nearfield(fileprefix, lmax, fmin, fmax)
%COMPUTE_SPHERICAL_MODE_COEFFICIENTS_NEARFIELD  TM/TE mode coefficients
%from FEKO near-field data.
%   Projects the FEKO near-field E and H onto the spherical vector mode
%   fields to obtain psiTM and psiTE (same vector-projection method as
%   matlab/workflows/general_vector_spherical_modes_FEKO.m).
%
%   result = COMPUTE_SPHERICAL_MODE_COEFFICIENTS_NEARFIELD(fileprefix, lmax)
%   result = COMPUTE_SPHERICAL_MODE_COEFFICIENTS_NEARFIELD(fileprefix, lmax, fmin, fmax)
%
%   fileprefix : string, path prefix of the FEKO output files, e.g.
%                fullfile(feko_dir, "cylinder_rotated_yaxis_45deg");
%                <prefix>_NearField1.efe, <prefix>_NearField1.hfe and
%                <prefix>.out must exist
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
%   Used by compare_rotated_antenna_coefficients.m.

if nargin < 3 || isempty(fmin)
    fmin = -Inf;
end
if nargin < 4 || isempty(fmax)
    fmax = Inf;
end

fileprefix = string(fileprefix);

%% --- Physical Constants ---
c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0);

%% --- FEKO Import: Read Port Parameters from .out File ---
[V0, Z0, R0, symmetry_planes, even_odd_excitation, a, Zant, Freq, port_pol] = ...
    ReadFEKOout(fileprefix+".out"); %#ok<ASGLU>

Zant = transpose(Zant(:,1)); Freq = transpose(Freq);
Zindx = find(~(real(Zant)<0));
Freq  = Freq(Zindx);
Zant(find(real(Zant)<1e-4)) = Zant(find(real(Zant)<1e-4)) + 1e-4; %#ok<FNDSB>
Zant  = Zant(Zindx);

if abs(Z0-eta0) > 0.01
    error("antenna port impedance should be eta0 for now")
else
    Z0 = eta0; %#ok<NASGU>
end

%% --- Import FEKO Near Field Data ---
[Eptr, Rs, theta, phi, ~] = ReadFE_spherical(fileprefix+"_NearField1.efe"); Eptr = Eptr(:,:,:,:,Zindx);
[Hptr, ~,  ~,     ~,   ~] = ReadFE_spherical(fileprefix+"_NearField1.hfe"); Hptr = Hptr(:,:,:,:,Zindx);
theta = theta(:);    % (numth x 1)
phi   = phi(:).';    % (1 x numph)

%% --- Restrict to requested frequency range ---
frange_idx = find(Freq >= fmin & Freq <= fmax);
if isempty(frange_idx)
    error('compute_spherical_mode_coefficients:noFreqInRange', ...
          'No frequency points found in [%.4g, %.4g] Hz (available range: [%.4g, %.4g] Hz)', ...
          fmin, fmax, min(Freq), max(Freq));
end
Freq = Freq(frange_idx);
Eptr = Eptr(:,:,:,:,frange_idx);
Hptr = Hptr(:,:,:,:,frange_idx);

numf = length(Freq);

[~, index_r] = min(abs(a+c0/2/pi/min(Freq)-Rs), [], 'all');

% Symmetry factors
if max(theta) > (pi/2)
    sym_th = 1;
else
    sym_th = 2;
end
if max(phi) > (pi/2)
    sym_phi = 1;
elseif max(phi) > (pi/4)
    sym_phi = 2;
else
    sym_phi = 4;
end

k = 2*pi*Freq/c0;
numth = length(theta); numph = length(phi); %#ok<NASGU>

%% --- Mode index setup ---
lval = (1:1:lmax).';   % mode orders l = 1..lmax
mval = (-lmax:lmax);   % m ranges from -lmax to +lmax
L    = length(lval);
M    = length(mval);

psiTM = nan(L,M,numf); psiTE = psiTM;

l_factor = lval .* (lval + 1);

use_jl1_convention = true;
if use_jl1_convention
    j_l_1 = (1j) .^ (lval + 1);
else
    j_l_1 = ones(L, 1);
end

kr = k * Rs(index_r);
[PHI, THETA] = meshgrid(phi, theta);

prop_factor_norm = ones(L,numf);
jrFlm_phi = zeros(numth,numph,L,M,numf); jrFlm_theta = jrFlm_phi; jrFlm_r = jrFlm_phi;
denom_vec = ones(L,numf);

%% --- Determine psiTM, psiTE via vector-method projection ---
for lindx = 1:L
    ell = lval(lindx);
    for f = 1:numf
        prop_factor_norm(lindx,f) = sphericalHankel(lval(lindx),2,kr(f)) / kr(f) / j_l_1(lindx);
        denom_vec(lindx,f) = l_factor(lindx)*sphericalHankel(ell,2,kr(f))*sphericalHankel(ell,1,kr(f)) ...
                    + sphericalRicatti_h_prime(ell,2,kr(f))*sphericalRicatti_h_prime(ell,1,kr(f));
    end

    for mindx = 1:M
        if abs(mval(mindx)) > lval(lindx)
            continue
        end

        for f = 1:numf
            [jrFlm_phi(:,:,lindx,mindx,f), jrFlm_theta(:,:,lindx,mindx,f), jrFlm_r(:,:,lindx,mindx,f)] ...
                = Spherical_Mode_Field(lval(lindx),mval(mindx),1,THETA,PHI,kr(f));

            psiTM(lindx,mindx,f) = sym_th * sym_phi * 1j*Rs(index_r)/sqrt(eta0) * ...
                trapz(phi, ...
                        trapz( theta, ( Eptr(:,:,index_r,1,f).*jrFlm_phi(:,:,lindx,mindx,f) + ...
                                        Eptr(:,:,index_r,2,f).*jrFlm_theta(:,:,lindx,mindx,f) + ...
                                        Eptr(:,:,index_r,3,f).*jrFlm_r(:,:,lindx,mindx,f) ...
                                       ).*sin(THETA) ...
                        ,1) ...
                ,2) / denom_vec(lindx,f);

            psiTE(lindx,mindx,f) = sym_th * sym_phi * 1j*Rs(index_r)*sqrt(eta0) * ...
                trapz(phi, ...
                        trapz( theta, ( Hptr(:,:,index_r,1,f).*jrFlm_phi(:,:,lindx,mindx,f) + ...
                                        Hptr(:,:,index_r,2,f).*jrFlm_theta(:,:,lindx,mindx,f) + ...
                                        Hptr(:,:,index_r,3,f).*jrFlm_r(:,:,lindx,mindx,f) ...
                                       ).*sin(THETA) ...
                        ,1) ...
                ,2) / denom_vec(lindx,f);
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