function nf = compute_port_mode_coeffs_nearfield(fileprefix, near_name, near_occ, near_type, Freq, Zindx, a, ...
        lval, mval, L, M, use_vector_method, l_factor, port_label)
%COMPUTE_PORT_MODE_COEFFS_NEARFIELD  Per-port TM/TE mode coefficients from FEKO .efe/.hfe exports.
%   nf = compute_port_mode_coeffs_nearfield(fileprefix, near_name, near_occ,
%        near_type, Freq, Zindx, a, lval, mval, L, M, use_vector_method,
%        l_factor, port_label)
%   Multi-port form of compute_mode_coeffs_nearfield: the .efe file of one
%   internal port is located with resolve_field_file (the file name can
%   differ from <fileprefix>_<near_name>.efe) and the matching .hfe is
%   <fileprefix>_<resolved name>.hfe.
%
%   The near-field files are always read when found (the caller also
%   uses Eptr/Hptr and the Poynting-vector power Pr_NF); use_vector_method
%   only controls the projection. The sphere of radius Rs(index_r) closest
%   to a + lambda_max/(2*pi) is used, and the fields are projected onto the
%   outgoing-wave basis Spherical_Mode_Field(l,m,1,...):
%       psiTM = j Rs/sqrt(eta0) * INT E . jrF_lm dOmega / D_l
%       psiTE = j Rs sqrt(eta0) * INT H . jrF_lm dOmega / D_l
%       D_l   = l(l+1) h_l^(2) h_l^(1) + (kr h_l^(2))' (kr h_l^(1))'
%   at kr = k Rs, with trapezoidal integration.
%
%   fileprefix    : string, path prefix of the FEKO output files
%   near_name, near_occ, near_type : the port's near-field request name,
%                   occurrence count and reuse type (see
%                   resolve_field_file); if near_name is empty or no file
%                   is found, nf.available is false
%   Freq, Zindx   : (1 x numf) frequencies and the indices into the
%                   file's frequency list selected by the caller
%   a             : antenna apothem (m); selects the extraction radius
%   lval, mval, L, M, l_factor : mode-index grids built by the caller from
%                   lmax (l_factor = l(l+1))
%   use_vector_method : compute psiTM/psiTE by projection; forced off
%                   (with a message) if the grid does not cover the full
%                   sphere (theta 0..pi, phi 0..2*pi within 1 deg)
%   port_label    : port name used in printed messages
%
% Returns struct nf (fields are []/NaN/false when .available is false):
%   .available     - near-field files found and read
%   .used_vector   - psiTM/psiTE were computed by projection
%   .psiTM, .psiTE - (L x M x numf) coefficients; NaN if not computed
%   .Pr_NF         - (numf x 1) radiated power from the radial Poynting vector
%   .theta, .phi, .Rs, .index_r, .k, .kr, .THETA, .PHI - grid used
%                    (THETA/PHI are numth x numph)
%   .jrFlm_r       - radial component of the mode basis on that grid
%   .Eptr, .Hptr   - near-field data, (numth x numph x numr x 3 x numf),
%                    dim 4 = [phi, theta, r]
%
% See also Spherical_Mode_Field, ReadFE_spherical.

c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0);

numf = numel(Freq);

nf = struct('available', false, 'used_vector', false, ...
    'psiTM', nan(L,M,numf), 'psiTE', nan(L,M,numf), ...
    'Pr_NF', nan(numf,1), ...
    'theta', [], 'phi', [], 'Rs', [], 'index_r', NaN, 'k', [], 'kr', [], ...
    'THETA', [], 'PHI', [], 'jrFlm_r', [], 'Eptr', [], 'Hptr', []);

if strlength(near_name) == 0
    fprintf('compute_port_mode_coeffs_nearfield (%s): no near-field request name given -- this port has no spherical-mode coupling.\n', port_label);
    return
end

[efe_file, resolved_near_name] = resolve_field_file(fileprefix, near_name, near_occ, near_type, "efe");
if strlength(efe_file) == 0
    fprintf('compute_port_mode_coeffs_nearfield (%s): no near-field .efe file found for request "%s" (occurrence %d; tried the plain name and disambiguated candidates -- see resolve_field_file) -- this port has no spherical-mode coupling.\n', port_label, near_name, near_occ);
    return
end
hfe_file = fileprefix + "_" + resolved_near_name + ".hfe";
if ~isfile(hfe_file)
    fprintf('compute_port_mode_coeffs_nearfield (%s): found near-field .efe file (%s) but the matching .hfe file was not found (%s) -- this port has no spherical-mode coupling.\n', port_label, efe_file, hfe_file);
    return
end

use_vector_method_local = use_vector_method;

% Eptr/Hptr: (numth x numph x numr x 3 x numf); dim 4 = [E_phi, E_theta, E_r]
[Eptr, Rs, theta, phi, ~] = ReadFE_spherical(efe_file); Eptr = Eptr(:,:,:,:,Zindx);
[Hptr, ~,  ~,     ~,   ~] = ReadFE_spherical(hfe_file); Hptr = Hptr(:,:,:,:,Zindx);
theta = theta(:);    % (numth x 1)
phi   = phi(:).';    % (1 x numph)

if use_vector_method_local && ~local_full_angular_range(theta, phi)
    use_vector_method_local = false;
    fprintf('compute_port_mode_coeffs_nearfield (%s): near-field theta/phi grid does not cover the full angular range (theta: 0 to pi, phi: 0 to 2*pi) — forcing vector method off for this port.\n', port_label);
end

[~,index_r] = min(abs(a+c0/2/pi/min(Freq)-Rs),[],'all');
k = 2*pi*Freq/c0;
numth = length(theta); numph = length(phi); numr = length(Rs);

Poynting = NaN(numth,numph,numr,3,numf);
Poynting(:,:,:,1,:) = ( Eptr(:,:,:,3,:).*conj(Hptr(:,:,:,2,:)) - Eptr(:,:,:,2,:).*conj(Hptr(:,:,:,3,:)) )/2;
Poynting(:,:,:,2,:) = ( Eptr(:,:,:,1,:).*conj(Hptr(:,:,:,3,:)) - Eptr(:,:,:,3,:).*conj(Hptr(:,:,:,1,:)) )/2;
Poynting(:,:,:,3,:) = ( Eptr(:,:,:,2,:).*conj(Hptr(:,:,:,1,:)) - Eptr(:,:,:,1,:).*conj(Hptr(:,:,:,2,:)) )/2;

psiTM_vec  = nan(L,M,numf);  psiTE_vec  = psiTM_vec;
jrFlm_phi = zeros(numth,numph,L,M,numf); jrFlm_theta=jrFlm_phi; jrFlm_r=jrFlm_phi;
denom_vec = ones(L,numf);

kr  = k * Rs(index_r);
[PHI, THETA] = meshgrid(phi, theta);

for lindx = 1:L
    ell = lval(lindx);
    for f = 1:numf
        denom_vec(lindx,f) = l_factor(lindx)*sphericalHankel(ell,2,kr(f))*sphericalHankel(ell,1,kr(f))...
                    + sphericalRicatti_h_prime(ell,2,kr(f))*sphericalRicatti_h_prime(ell,1,kr(f)) ;
    end

    for mindx = 1:M
        if abs(mval(mindx)) > lval(lindx)
            continue
        end

        for f = 1:numf
            [jrFlm_phi(:,:,lindx,mindx,f),...
            jrFlm_theta(:,:,lindx,mindx,f),...
            jrFlm_r(:,:,lindx,mindx,f)]...
            = Spherical_Mode_Field(lval(lindx),mval(mindx),1,THETA,PHI,kr(f)) ;

            if ~use_vector_method_local
                continue
            end

            psiTM_vec(lindx,mindx,f) = 1j*Rs(index_r)/sqrt(eta0) * ...
                trapz(phi,...
                        trapz( theta, ( Eptr(:,:,index_r,1,f).*jrFlm_phi(:,:,lindx,mindx,f) + ...
                                        Eptr(:,:,index_r,2,f).*jrFlm_theta(:,:,lindx,mindx,f) + ...
                                        Eptr(:,:,index_r,3,f).*jrFlm_r(:,:,lindx,mindx,f)...
                                       ).*sin(THETA)...
                        ,1) ...
                ,2) / denom_vec(lindx,f) ;

            psiTE_vec(lindx,mindx,f) = 1j*Rs(index_r)*sqrt(eta0) * ...
                trapz(phi,...
                        trapz( theta, ( Hptr(:,:,index_r,1,f).*jrFlm_phi(:,:,lindx,mindx,f) + ...
                                        Hptr(:,:,index_r,2,f).*jrFlm_theta(:,:,lindx,mindx,f) + ...
                                        Hptr(:,:,index_r,3,f).*jrFlm_r(:,:,lindx,mindx,f)...
                                       ).*sin(THETA)...
                        ,1) ...
                ,2) / denom_vec(lindx,f) ;
        end
    end
end
clearvars jrFlm_phi jrFlm_theta

Pr_NF = zeros(numf,1);
for f = 1:numf
    Pr_NF(f) = (Rs(index_r))^2 * ...
               trapz(phi, trapz(theta, Poynting(:,:,index_r,3,f).*sin(THETA), 1), 2);
end

nf.available   = true;
nf.used_vector = use_vector_method_local;
nf.psiTM = psiTM_vec; nf.psiTE = psiTE_vec;
nf.Pr_NF = Pr_NF;
nf.theta = theta; nf.phi = phi; nf.Rs = Rs; nf.index_r = index_r;
nf.k = k; nf.kr = kr; nf.THETA = THETA; nf.PHI = PHI;
nf.jrFlm_r = jrFlm_r;
nf.Eptr = Eptr; nf.Hptr = Hptr;

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
