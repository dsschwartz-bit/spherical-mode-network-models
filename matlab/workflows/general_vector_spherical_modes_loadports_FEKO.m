%GENERAL_VECTOR_SPHERICAL_MODES_LOADPORTS_FEKO  Multi-port spherical-mode network model from FEKO.
%
% Builds an S-parameter network model of an antenna with one driven "Feed"
% port and any number of internal load ports ("load1", "load2", ...) from a
% FEKO solution, and exports it as a Touchstone file.
%
% Required FEKO exports (in the FEKO project folder, prefix = fileprefix):
%   <fileprefix>.out      multi-port S-parameter configuration plus one
%                         fields configuration per excited port
%   *.efe / *.hfe         near fields (vector method)
%   *.ffe                 far fields (far-field method, optional)
%   *.sph                 FEKO spherical-mode Q coefficients (optional)
% Field files are located from the request names reported in the .out file
% (resolve_field_file). A port without field requests gets zero
% spherical-mode coupling.
%
% Processing:
%   1. Read ports, S-matrix and geometry (ReadFEKOout_loadports).
%   2. Per port, compute psi coefficients (near-field and/or far-field
%      method, averaged) and the converted FEKO Q coefficients.
%   3. Rank the n_dom dominant (l,m,TE/TM) modes by accepted-power coupling
%      |T_lm| of the primary port (or the first port that has field data).
%   4. Build the (1+Nsec+n_dom)-port S-matrix: internal-port block from the
%      FEKO S-parameter configuration, port-to-mode coupling from psi (or Q
%      when Q is the only source). Mode ports are mutually uncoupled.
%   5. Plot impedance, field reconstruction, T_lm and radiated power per port.
%   6. Write S_dom_<psi|Q>_<antenna_name>.s<N>p (with a DC point) to snp_dir.
%
% Output port order: 1 = Feed, 2..1+Nsec = load ports, then the n_dom
% spherical-mode ports (ordered by snp_order_method).

clear
close all
clc
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
% fileprefix: FEKO project path without extension (the .out file is
% fileprefix + ".out"). Uncomment one antenna_name/fileprefix pair.
antenna_name = "loop_loadport";
fileprefix = fullfile("C:\path\to\FEKO\loop", antenna_name);
% antenna_name = "bowtie_edgeport_45degrees";
% fileprefix = fullfile("C:\path\to\FEKO\bowtie", antenna_name);
% antenna_name = "dipole_cylinder_edgeport_symmetry_full";
% fileprefix = fullfile("C:\path\to\FEKO\dipole", antenna_name);

% Output folder for the Touchstone file (created if missing).
snp_dir = fullfile(repo_root, 'output', 'SnP');
% =========================

%% --- Physical Constants ---
c0  = physconst('lightspeed');
mu0 = 4*pi*1e-7;          % Permeability of free space (H/m)
eps0 = 1/c0^2/mu0;        % Permittivity of free space (F/m)
eta0 = sqrt(mu0/eps0);    % Characteristic impedance of free space (Ohm)

%% --- Simulation Parameters ---
f0     = 1e9;         % Reference frequency for single-frequency plots (Hz)
fmin   = 0.1e9;       % Minimum frequency for plot x-axis (Hz)
fmax   = 1e9;         % Maximum frequency for plot x-axis (Hz)
lmax   = 5;           % Maximum spherical mode order to compute
n_dom  = 6;           % Number of dominant modes to identify and plot

% psi- or Q-derived coefficients for ranking/S-parameters/SnP export:
% set automatically below (Q only when it is the sole data source).
snp_use_Q = false;   % placeholder; overwritten below

% Coefficient sources (at least one must be true):
%   use_vector_method - project near-field E/H (.efe/.hfe) onto vector
%                       spherical waves. Near fields are read regardless,
%                       for the Poynting power check and field plots.
%   use_FF_method     - fit far-field (.ffe) data; averaged with the
%                       vector-method result.
%   use_Q_method      - read FEKO .sph Q coefficients. Q enters the psi
%                       average only when no near/far-field data is
%                       available; otherwise it is used for psi-vs-Q
%                       comparison plots. false disables all Q outputs.
use_vector_method = true;
use_FF_method = true;
use_Q_method = true;

% Order of the spherical-mode ports in the SnP file:
%   'peak_freq'  - ascending frequency of peak |T_lm|
%   'mode_index' - ascending FEKO mode index j = 2*(l*(l+1)+m-1)+s,
%                  s=1 for TE, s=2 for TM
%snp_order_method = 'peak_freq';
snp_order_method = 'mode_index';

% true: the j^(l+1) far-field phase is absorbed into psi/Q (via j_l_1), so
% the coupling formula uses jl1 = 1. false: the coupling formula applies
% jl1 = (1j)^(l+1) explicitly. Passed to the coefficient functions so both
% use the same convention.
use_jl1_convention = true;

% true: multiply the port-to-mode coupling by exp(-j*k*a), moving the mode
% reference plane to r = a (antenna apothem). false: no de-embedding phase.
apply_ka_deembedding = false;

% DC (0 Hz) point of the SnP export, per internal port (validated below):
%   'parallel_loads' - DC resistance = parallel combination of the other
%                      internal ports' load impedances (e.g. a loop with a
%                      feed gap and load gaps)
%   'series_loads'   - DC resistance = series sum of the other ports' loads
%   'Rmin'           - each port's input resistance at the lowest simulated
%                      frequency (radiation resistance vanishes as f->0)
%   'user'           - prompt for each value with input()
% Off-diagonal DC terms use the N-port divider formula for
% 'parallel_loads'/'series_loads', otherwise the real-valued magnitude of
% the lowest-frequency S-parameter (or a prompted value for 'user').
dc_resistance_method = 'parallel_loads';
%% --- Plot Settings ---
LW = 2;    % LineWidth
FS = 14;   % FontSize

set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% --- FEKO Import: Read Port Parameters from .out File ---
% V0: complex excitation voltage of the primary (Feed) port; Z0: port
% reference impedance; R0: primary port location;
% symmetry_planes.physical_gnd: PEC ground-plane flag; even_odd_excitation:
% log2(number of symmetric ports) if mode symmetry is used; a: antenna apothem.
% secondary_ports: one struct per internal load port (V0_own, geometry,
%   field-request names, sparam_port_index, load impedances).
% S: N x N x numf S-matrix of the internal ports from the S-parameter
%   configuration; the only impedance source used here:
%   Zin_n = Z0*(1+S(n,n,:))./(1-S(n,n,:)).
% near_field_name/far_field_name (+ occurrence/type): the primary port's
%   field requests, used to locate its .efe/.hfe/.ffe/.sph files.
[V0, Z0, R0, symmetry_planes, even_odd_excitation, a, Freq, port_pol, ...
    secondary_ports, S, fields_config_name, near_field_name, far_field_name, primary_sparam_index, ...
    near_field_occurrence, far_field_occurrence, near_field_type, far_field_type, sparam_diagnostic] ...
    = ReadFEKOout_loadports(fileprefix+".out");

if isempty(S) || isempty(primary_sparam_index)
    % sparam_diagnostic says why: no S-parameter configuration, a 1-port
    % one, or one that does not include the primary port.
    if isempty(sparam_diagnostic)
        sparam_diagnostic = 'No S-parameter configuration was found in the .out file, for an unspecified reason (see the ReadFEKOout_loadports console output above for parsing details).';
    end
    error("general_vector_spherical_modes_loadports: %s", sparam_diagnostic)
end

if abs(Z0-eta0) > 0.01
    error("antenna port impedance should be eta0 for now")
else
    Z0 = eta0;Y0 = 1/Z0;
end

% Feed input impedance from S; frequencies with negative input resistance
% are discarded (Zindx) and S is trimmed to match.
Freq=transpose(Freq);
Zant = Z0 * (1 + reshape(S(primary_sparam_index, primary_sparam_index, :), 1, [])) ./ ...
            (1 - reshape(S(primary_sparam_index, primary_sparam_index, :), 1, []));
indicies_negR=find(real(Zant)<0);
Zindx=find(~(real(Zant)<0));Freq=Freq(Zindx);
Zant(find(real(Zant)<1e-4))=Zant(find(real(Zant)<1e-4))+1e-4;
Zant=Zant(Zindx);
S = S(:,:,Zindx);

numf = length(Freq);

if fmax>Freq(numf)
    fmax=Freq(numf);
end
if fmin<Freq(1)
    fmin=Freq(1);
end

%% --- Build the internal-port list: Feed (primary) + load1, load2, ... ---
% Primary first, then secondary_ports in the order returned by
% ReadFEKOout_loadports. The n_dom mode ports follow in the SnP export.
n_sec  = numel(secondary_ports);
num_internal = 1 + n_sec;

int_name       = strings(1, num_internal);
int_near_name  = strings(1, num_internal);
int_far_name   = strings(1, num_internal);
int_near_occ   = zeros(1, num_internal);
int_far_occ    = zeros(1, num_internal);
int_near_type  = strings(1, num_internal);
int_far_type   = strings(1, num_internal);
int_sparam_idx = nan(1, num_internal);

int_name(1)       = "Feed";
int_near_name(1)  = string(near_field_name);
int_far_name(1)   = string(far_field_name);
int_near_occ(1)   = near_field_occurrence;
int_far_occ(1)    = far_field_occurrence;
int_near_type(1)  = string(near_field_type);
int_far_type(1)   = string(far_field_type);
int_sparam_idx(1) = primary_sparam_index;

for k_sec = 1:n_sec
    int_name(k_sec+1)       = "load" + num2str(k_sec);
    int_near_name(k_sec+1)  = string(secondary_ports(k_sec).near_field_name);
    int_far_name(k_sec+1)   = string(secondary_ports(k_sec).far_field_name);
    int_near_occ(k_sec+1)   = secondary_ports(k_sec).near_field_occurrence;
    int_far_occ(k_sec+1)    = secondary_ports(k_sec).far_field_occurrence;
    int_near_type(k_sec+1)  = string(secondary_ports(k_sec).near_field_type);
    int_far_type(k_sec+1)   = string(secondary_ports(k_sec).far_field_type);
    int_sparam_idx(k_sec+1) = secondary_ports(k_sec).sparam_port_index;
end

fprintf('general_vector_spherical_modes_loadports: %d internal port(s): %s (FEKO port: %s)\n', ...
    num_internal, strjoin(int_name, ', '), strjoin(["Feed=" + int_name(1), ...
    arrayfun(@(k) int_name(k+1)+"="+secondary_ports(k).name, 1:n_sec)], ', '));

if any(isnan(int_sparam_idx))
    fprintf('general_vector_spherical_modes_loadports: WARNING - not every internal port was found in the S-parameter configuration; the internal-port block of the exported S-matrix will be incomplete for: %s.\n', ...
        strjoin(int_name(isnan(int_sparam_idx)), ', '));
end

%% --- Compute Spherical-Mode Coefficients for Every Internal Port ---
% Each port's near-field, far-field and Q coefficients come from
% compute_port_mode_coeffs_nearfield/_farfield/_Q; the available psi
% estimates are averaged here. Field files are located by
% resolve_field_file from the request name, occurrence count and reuse
% type; a port without a usable file gets zero mode coupling.
%
% Mode-index grids and normalisation depend only on lmax and the
% conventions below, so they are computed once and shared by all ports.
numf = numel(Freq);
lval = (1:1:lmax).'; mval = (-lmax:lmax); L = length(lval); M = length(mval);
l_factor = lval .* (lval + 1);
if use_jl1_convention
    j_l_1 = (1j) .^ (lval + 1);
else
    j_l_1 = ones(L, 1);
end
Norm = zeros(L,M);
for lindx=1:L
    for mindx=1:M
        if abs(mval(mindx)) > lval(lindx)
            Norm(lindx,mindx)=NaN;
        else
            Norm(lindx,mindx) = sqrt( (2^(mval(mindx)==0))*pi * 2*...
                    factorial(lval(lindx)+abs(mval(mindx))) / factorial(lval(lindx)-abs(mval(mindx)))...
                        / (2*lval(lindx)+1) );
        end
    end
end
have_CSphase = false;
if have_CSphase
    CSphase=(-1).^mval;
else
    CSphase=ones(size(mval));
    Norm = ((-1).^mval).*Norm ;
end
Normp = Norm.*sqrt(l_factor);
midx_0 = find(mval == 0);

% port_data{p}: per-port coefficients and near-field data, used by the
% S-parameter build and the per-port plots below.
port_data = cell(1, num_internal);
for p_idx = 1:num_internal
    nf_p = compute_port_mode_coeffs_nearfield(fileprefix, int_near_name(p_idx), int_near_occ(p_idx), int_near_type(p_idx), ...
        Freq, Zindx, a, lval, mval, L, M, use_vector_method, l_factor, int_name(p_idx));
    ff_p = compute_port_mode_coeffs_farfield(fileprefix, int_far_name(p_idx), int_far_occ(p_idx), int_far_type(p_idx), ...
        Freq, Zindx, lval, mval, L, M, use_FF_method, int_name(p_idx));
    qd_p = compute_port_mode_coeffs_Q(fileprefix, int_far_name(p_idx), int_far_occ(p_idx), int_far_type(p_idx), ...
        Freq, Zindx, lmax, lval, mval, L, M, symmetry_planes, j_l_1, CSphase, use_Q_method, int_name(p_idx));

    used_vector_p = nf_p.used_vector;
    used_FF_p     = ff_p.available;
    used_Q_p      = qd_p.available;
    has_field_data_p = used_vector_p || used_FF_p;

    psi_stack_TM = []; psi_stack_TE = [];
    if used_vector_p
        psi_stack_TM = cat(4, psi_stack_TM, nf_p.psiTM);
        psi_stack_TE = cat(4, psi_stack_TE, nf_p.psiTE);
    end
    if used_FF_p
        psi_stack_TM = cat(4, psi_stack_TM, ff_p.psiTM);
        psi_stack_TE = cat(4, psi_stack_TE, ff_p.psiTE);
    end
    if used_Q_p && ~has_field_data_p
        psi_stack_TM = cat(4, psi_stack_TM, qd_p.QTMconv);
        psi_stack_TE = cat(4, psi_stack_TE, qd_p.QTEconv);
    end

    pd = struct('available', ~isempty(psi_stack_TM), ...
        'psiTM', nan(L,M,numf), 'psiTE', nan(L,M,numf), ...
        'QTMconv', qd_p.QTMconv, 'QTEconv', qd_p.QTEconv, ...
        'Q', qd_p.Q, 'freqs_Q', qd_p.freqs_Q, 'mmax', qd_p.mmax, 'num_modes', qd_p.num_modes, ...
        'Pr_NF', nf_p.Pr_NF, ...
        'theta', nf_p.theta, 'phi', nf_p.phi, 'Rs', nf_p.Rs, 'index_r', nf_p.index_r, ...
        'k', nf_p.k, 'kr', nf_p.kr, ...
        'jrFlm_r', nf_p.jrFlm_r, 'Eptr', nf_p.Eptr, 'Hptr', nf_p.Hptr, ...
        'used_vector', used_vector_p, 'used_FF', used_FF_p, 'used_Q', used_Q_p);
    if pd.available
        pd.psiTM = mean(psi_stack_TM, 4, 'omitnan');
        pd.psiTE = mean(psi_stack_TE, 4, 'omitnan');
    else
        fprintf('general_vector_spherical_modes_loadports (%s): no usable method (vector/FF/Q) produced data for this port -- it has no spherical-mode coupling.\n', int_name(p_idx));
    end
    port_data{p_idx} = pd;
end
clearvars nf_p ff_p qd_p pd psi_stack_TM psi_stack_TE used_vector_p used_FF_p used_Q_p has_field_data_p

% --- Port used for dominant-mode ranking: the primary port if it has
% field data, otherwise the first port that does (the primary fields
% configuration may only compute Zin/S-parameters). ---
rank_port = 1;
if ~port_data{1}.available
    rank_port = find(cellfun(@(pd) pd.available, port_data), 1);
    if isempty(rank_port)
        error('general_vector_spherical_modes_loadports: no internal port has usable spherical-mode data (near-field/far-field/Q all unavailable for every port -- check that at least one fields configuration requests near-field or far-field output).')
    end
    fprintf('general_vector_spherical_modes_loadports: WARNING - primary port "%s" has no field data; ranking dominant modes from port "%s" instead.\n', ...
        int_name(1), int_name(rank_port));
end
pd_rank = port_data{rank_port};

% --- Unpack the ranking port's data into top-level variables. ---
psiTM = pd_rank.psiTM;       psiTE   = pd_rank.psiTE;
QTMconv = pd_rank.QTMconv;   QTEconv = pd_rank.QTEconv;
theta = pd_rank.theta; phi = pd_rank.phi; Rs = pd_rank.Rs;
index_r = pd_rank.index_r;
[~,index_f] = min(abs(Freq-f0));
k = pd_rank.k; kr = pd_rank.kr;
Pr_NF = pd_rank.Pr_NF;
Q = pd_rank.Q; freqs_Q = pd_rank.freqs_Q; mmax = pd_rank.mmax; num_modes = pd_rank.num_modes;
use_vector_method = pd_rank.used_vector;
use_FF_method     = pd_rank.used_FF;
use_Q_method      = pd_rank.used_Q;
has_field_data    = use_vector_method || use_FF_method;

% Use Q-derived coefficients only when Q is the ranking port's sole source.
snp_use_Q = use_Q_method && ~has_field_data;
if snp_use_Q
    fprintf('general_vector_spherical_modes_loadports: only Q data is available for port "%s" — using Q-derived coefficients for S-parameters/ranking/SnP export.\n', int_name(rank_port));
end

%% --- Build S-parameters for EVERY (l,m,type), for Ranking (ranking port only) ---
% Modes are ranked by accepted-power coupling
% T_lm = S_{mode,port}/sqrt(1-|S_pp|^2) rather than raw |psi|, so a mode is
% not counted as dominant where the port is poorly matched. Same coupling
% formula as the SnP build below, applied to every mode of the ranking port.
rank_p_idx = int_sparam_idx(rank_port);
S_pp_rank = reshape(S(rank_p_idx, rank_p_idx, :), 1, numf);   % ranking port's own S_pp
if rank_port == 1
    rank_V0 = V0(1);
else
    rank_V0 = secondary_ports(rank_port-1).V0_own;
end
port_factor = sqrt(Z0 / 2^abs(even_odd_excitation));
if apply_ka_deembedding
    ka_factor = exp(-1j * k(:).' * a);                   % (1 x numf)
else
    ka_factor = ones(1, numf);                           % no de-embedding phase
end
ka_factor_r = reshape(ka_factor, 1, 1, numf);
if use_jl1_convention
    jl1_l = ones(L,1);                                   % phase already absorbed into psiTM/psiTE/QTMconv/QTEconv via j_l_1
else
    jl1_l = (1j).^(lval(:) + 1);                         % (L x 1)
end

one_plus_Spp_r = reshape(1 + S_pp_rank, 1, 1, numf);
jl1_r          = reshape(jl1_l, L, 1, 1);                % broadcasts over M and numf

psiTM_safe = psiTM; psiTM_safe(isnan(psiTM_safe)) = 0;
psiTE_safe = psiTE; psiTE_safe(isnan(psiTE_safe)) = 0;

S_TM_all   = port_factor .* one_plus_Spp_r .* (psiTM_safe ./ rank_V0) .* ka_factor_r ./ jl1_r;   % (L x M x numf)
S_TE_all   = port_factor .* one_plus_Spp_r .* (psiTE_safe ./ rank_V0) .* ka_factor_r ./ jl1_r;

T_denom_rank = sqrt(1 - abs(S_pp_rank).^2);              % (1 x numf)
T_TM_all   = S_TM_all   ./ reshape(T_denom_rank, 1, 1, numf);   % (L x M x numf)
T_TE_all   = S_TE_all   ./ reshape(T_denom_rank, 1, 1, numf);

% Q-derived coupling is only needed when Q was read.
if use_Q_method
    S_TM_Q_all = port_factor .* one_plus_Spp_r .* (QTMconv ./ rank_V0) .* ka_factor_r ./ jl1_r;
    S_TE_Q_all = port_factor .* one_plus_Spp_r .* (QTEconv ./ rank_V0) .* ka_factor_r ./ jl1_r;
    T_TM_Q_all = S_TM_Q_all ./ reshape(T_denom_rank, 1, 1, numf);
    T_TE_Q_all = S_TE_Q_all ./ reshape(T_denom_rank, 1, 1, numf);
end

%% --- Identify Dominant Modal Coefficients across TE and TM jointly ---
% TE and TM peak |T_lm| over [fmin, fmax] are pooled into one list and
% the top n_dom entries (either type) are selected.

fband_idx = find(Freq >= fmin & Freq <= fmax);

% Max |T_lm| in band per (l,m) and type, and its index into Freq.
if use_Q_method
    [QTM_mag, QTM_fmax_band] = max(abs(reshape(T_TM_Q_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
    [QTE_mag, QTE_fmax_band] = max(abs(reshape(T_TE_Q_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
    QTM_fmax = fband_idx(QTM_fmax_band(:)).';
    QTE_fmax = fband_idx(QTE_fmax_band(:)).';
end
[psiTM_mag, psiTM_fmax_band] = max(abs(reshape(T_TM_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
[psiTE_mag, psiTE_fmax_band] = max(abs(reshape(T_TE_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
psiTM_fmax = fband_idx(psiTM_fmax_band(:)).';
psiTE_fmax = fband_idx(psiTE_fmax_band(:)).';

if snp_use_Q
    rank_TM_mag = QTM_mag;   rank_TM_fmax = QTM_fmax;
    rank_TE_mag = QTE_mag;   rank_TE_fmax = QTE_fmax;
    rank_src = 'Q';
else
    rank_TM_mag = psiTM_mag; rank_TM_fmax = psiTM_fmax;
    rank_TE_mag = psiTE_mag; rank_TE_fmax = psiTE_fmax;
    rank_src = 'psi';
end

% Replace NaN with 0 so invalid (|m|>l) modes sort to the bottom
rank_TM_mag(isnan(rank_TM_mag)) = 0;
rank_TE_mag(isnan(rank_TE_mag)) = 0;

% Pool [TM; TE] entries as columns (type 0 = TM, 1 = TE); (:) is needed
% because the *_fmax vectors are rows.
pool_mag     = [rank_TM_mag(:);  rank_TE_mag(:)  ];
pool_fmax    = [rank_TM_fmax(:); rank_TE_fmax(:) ];
pool_type    = [zeros(L*M,1); ones(L*M,1)];
pool_lmidx   = [(1:L*M).';   (1:L*M).'  ];

% Top n_dom by magnitude, then sorted by ascending peak frequency.
[~, sort_pool] = sort(pool_mag, 'descend');
dom_pool_mag   = sort_pool(1:n_dom);
dom_peak_freqs = Freq(pool_fmax(dom_pool_mag));
[~, freq_order] = sort(dom_peak_freqs(:), 'ascend');
dom_pool = dom_pool_mag(freq_order);

% Decode each dominant entry into type, (lindx, mindx)
dom_type = pool_type(dom_pool);      % 0=TM, 1=TE, length n_dom
dom_lm   = pool_lmidx(dom_pool);     % linear index into (L x M), length n_dom
[dom_lidx, dom_midx] = ind2sub([L, M], dom_lm);

% Separate TM and TE index lists (lengths sum to n_dom)
dom_TM_lidx = dom_lidx(dom_type == 0);
dom_TM_midx = dom_midx(dom_type == 0);
dom_TE_lidx = dom_lidx(dom_type == 1);
dom_TE_midx = dom_midx(dom_type == 1);
n_dom_TM = length(dom_TM_lidx);
n_dom_TE = length(dom_TE_lidx);

% Report dominant modes, in dB relative to the strongest mode at each peak
fprintf('\n--- Top %d dominant modes (TE and TM ranked jointly by max |T_%s|, accepted-power coupling) ---\n', n_dom, rank_src)
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0
        if snp_use_Q
            mode_vs_freq = abs(reshape(T_TM_Q_all(li,mi,:), 1, numf));
        else
            mode_vs_freq = abs(reshape(T_TM_all(li,mi,:), 1, numf));
        end
        type_str = 'TM';
    else
        if snp_use_Q
            mode_vs_freq = abs(reshape(T_TE_Q_all(li,mi,:), 1, numf));
        else
            mode_vs_freq = abs(reshape(T_TE_all(li,mi,:), 1, numf));
        end
        type_str = 'TE';
    end
    [~, f_max_idx_band] = max(mode_vs_freq(fband_idx));
    f_max_idx = fband_idx(f_max_idx_band);
    f_max = Freq(f_max_idx);
    if snp_use_Q
        all_at_f = [abs(reshape(T_TM_Q_all(:,:,f_max_idx), [], 1)); abs(reshape(T_TE_Q_all(:,:,f_max_idx), [], 1))];
    else
        all_at_f = [abs(reshape(T_TM_all(:,:,f_max_idx), [], 1)); abs(reshape(T_TE_all(:,:,f_max_idx), [], 1))];
    end
    global_max_at_f = max(all_at_f);
    rel_dB = 20*log10(mode_vs_freq(f_max_idx) / global_max_at_f);
    if rel_dB == 0
        fprintf('  Rank %d: %s l=%d, m=%d  — peak @ f=%.3g GHz, DOMINANT\n', ...
            kk, type_str, lval(li), mval(mi), f_max/1e9)
    else
        fprintf('  Rank %d: %s l=%d, m=%d  — peak @ f=%.3g GHz, %.1f dB rel. to dominant\n', ...
            kk, type_str, lval(li), mval(mi), f_max/1e9, rel_dB)
    end
end
fprintf('\n')

%% --- Build Spherical Wave S-parameters for Dominant Modes ---
% S_psi and S_Q are (num_internal+n_dom) x (num_internal+n_dom) x numf.
%   1..num_internal: internal ports (Feed, load1, ...). This block,
%     including port-to-port coupling, is taken from the FEKO S-parameter
%     configuration (reordered by int_sparam_idx).
%   num_internal+1..end: dominant modes (dom_lidx/dom_midx/dom_type). Each
%     port's coupling uses its own psi/Q coefficients (zero without field
%     data) and its own S_pp:
%       S_{mode,p} = sqrt(Z0/2^|even_odd_excitation|) * (1 + S_pp)
%                    * (psi_p/V0_p) * ka_factor / jl1
%     Mode ports are not coupled to each other.
n_tot = num_internal + n_dom;
S_psi = zeros(n_tot, n_tot, numf);
S_Q   = zeros(n_tot, n_tot, numf);

S_psi(1:num_internal, 1:num_internal, :) = S(int_sparam_idx, int_sparam_idx, :);
S_Q(1:num_internal, 1:num_internal, :)   = S_psi(1:num_internal, 1:num_internal, :);

% Excitation voltage per internal port (NaN/0 -> zero mode coupling).
int_V0 = nan(1, num_internal);
int_V0(1) = V0(1);
for k_sec = 1:n_sec
    int_V0(k_sec+1) = secondary_ports(k_sec).V0_own;
end

port_factor = sqrt(Z0 / 2^abs(even_odd_excitation));

for f = 1:numf
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        exp_ka = ka_factor(f);
        if use_jl1_convention
            jl1 = 1;                        % phase already absorbed into psiTM/psiTE via j_l_1
        else
            jl1 = (1j)^(lval(li) + 1);
        end

        for p = 1:num_internal
            pd = port_data{p};
            if dom_type(kk) == 0   % TM
                psi_val = NaN;
                if pd.available, psi_val = pd.psiTM(li,mi,f); end
                Q_val = 0;
                if pd.used_Q, Q_val = pd.QTMconv(li,mi,f); end
            else                   % TE
                psi_val = NaN;
                if pd.available, psi_val = pd.psiTE(li,mi,f); end
                Q_val = 0;
                if pd.used_Q, Q_val = pd.QTEconv(li,mi,f); end
            end
            if isnan(psi_val), psi_val = 0; end

            if isnan(int_V0(p)) || int_V0(p) == 0
                coupling_psi = 0; coupling_Q = 0;
            else
                coupling_psi = port_factor * (1 + S_psi(p,p,f)) * (psi_val / int_V0(p)) * exp_ka / jl1;
                coupling_Q   = port_factor * (1 + S_Q(p,p,f))   * (Q_val   / int_V0(p)) * exp_ka / jl1;
            end

            S_psi(num_internal+kk, p, f) = coupling_psi;
            S_psi(p, num_internal+kk, f) = coupling_psi;
            S_Q(num_internal+kk, p, f)   = coupling_Q;
            S_Q(p, num_internal+kk, f)   = coupling_Q;
        end
    end
end
% Sanitise: replace any remaining NaN with 0 before export
S_psi(isnan(S_psi)) = 0;
S_Q(isnan(S_Q))     = 0;

%% --- Per-Port Plots: Impedance, Reconstruction, T_lm, Power ---
% One set of figures per internal port that has field data
% (port_data{p}.available); the port name is added to each title. fignum
% is a single running figure counter (plot_field_comparison_polar returns
% the next value).

mode_labels = cell(1, n_dom);
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0
        mode_labels{kk} = sprintf('TM $l=%d,m=%d$', lval(li), mval(mi));
    else
        mode_labels{kk} = sprintf('TE $l=%d,m=%d$', lval(li), mval(mi));
    end
end

clr_S = get(groot,'defaultAxesColorOrder');
n_clr_S = size(clr_S, 1);

% snp_order: order of the n_dom mode ports in the SnP file and in every
% T_lm plot legend; based on the Feed port's T_lm so it is the same for
% all ports.
T_psi_feed = zeros(n_dom, numf);
T_Q_feed   = zeros(n_dom, numf);
denom_psi_feed = sqrt(1 - abs(reshape(S_psi(1,1,:), 1, numf)).^2);
denom_Q_feed   = sqrt(1 - abs(reshape(S_Q(1,1,:),   1, numf)).^2);
for kk = 1:n_dom
    T_psi_feed(kk,:) = reshape(S_psi(num_internal+kk,1,:), 1, numf) ./ denom_psi_feed;
    T_Q_feed(kk,:)   = reshape(S_Q(num_internal+kk,1,:),   1, numf) ./ denom_Q_feed;
end
if snp_use_Q
    T_for_order = T_Q_feed;
else
    T_for_order = T_psi_feed;
end
[~, T_fmax_band] = max(abs(T_for_order(:,fband_idx)), [], 2);
T_fmax = fband_idx(T_fmax_band);

if strcmpi(snp_order_method, 'mode_index')
    % FEKO mode index: j_ = 2*(l*(l+1) + m - 1) + s, s=1 for TE, s=2 for TM
    dom_l = lval(dom_lidx);
    dom_m = mval(dom_midx);
    s_val = ones(n_dom,1); s_val(dom_type == 0) = 2;   % TM=2, TE=1
    dom_j = 2*(dom_l.*(dom_l+1) + dom_m(:) - 1) + s_val;
    [~, snp_order] = sort(dom_j, 'ascend');
else   % 'peak_freq' (default)
    [~, snp_order] = sort(T_fmax, 'ascend');
end

% Dominant-mode masks for the radiated-power plot (Pr_PSIdom).
mask_TM_dom = false(L, M);
mask_TE_dom = false(L, M);
for kk = 1:n_dom_TM
    mask_TM_dom(dom_TM_lidx(kk), dom_TM_midx(kk)) = true;
end
for kk = 1:n_dom_TE
    mask_TE_dom(dom_TE_lidx(kk), dom_TE_midx(kk)) = true;
end

fignum = 0;   % running figure counter across all ports
for p = 1:num_internal
    pd = port_data{p};
    has_sparam = ~isnan(int_sparam_idx(p));

    % Plot only ports with their own fields configuration (field data).
    if ~pd.available
        fprintf('general_vector_spherical_modes_loadports: skipping plots for port "%s" -- no dedicated fields configuration (only S-parameter data).\n', int_name(p));
        continue
    end

    % --- Input Impedance (from S; skipped if the port is not in S) ---
    if has_sparam
        Zin_p = Z0 * (1 + reshape(S(int_sparam_idx(p), int_sparam_idx(p), :), 1, numf)) ./ ...
                     (1 - reshape(S(int_sparam_idx(p), int_sparam_idx(p), :), 1, numf));
        fignum = fignum+1; figure(fignum)
        plot(Freq, real(Zin_p), "Color","b","LineWidth",LW)
        hold on
        plot(Freq, imag(Zin_p), "Color","r","LineWidth",LW)
        hold off
        title("Impedance ("+int_name(p)+")","FontSize",FS)
        legend("$R$","$X$","FontSize",FS)
        xlim([fmin,fmax])
        ylim([-eta0,eta0])
        grid on
    else
        fprintf('general_vector_spherical_modes_loadports: port "%s" has no S-parameter data -- skipping its Input Impedance plot.\n', int_name(p));
    end

    % --- Hr / Er Reconstruction from the dominant modes ---
    % Each dominant entry contributes to Er (TM) or Hr (TE) only. Compared
    % with this port's FEKO near field at r = Rs(index_r).
    numth_p = numel(pd.theta); numph_p = numel(pd.phi);
    [~,index_phi_p] = min(abs(pd.phi - 0));
    [~,index_th_p]  = min(abs(pd.theta - pi/2));

    field_title_base = "@r="+num2str(pd.Rs(pd.index_r)*1000,2)+"mm, "+ ...
                       num2str(Freq(index_f)/1e9,2)+"GHz, \#dom="+num2str(n_dom)+", "+int_name(p);
    phi_label   = "$\phi="+num2str(pd.phi(index_phi_p)*180/pi,3)+"^\circ$";
    theta_label = "$\theta="+num2str(pd.theta(index_th_p)*180/pi,3)+"^\circ$";

    Er_modes = zeros(numth_p,numph_p,n_dom,numf); Hr_modes = Er_modes;
    for f = 1:numf
        for kk = 1:n_dom
            li = dom_lidx(kk); mi = dom_midx(kk);
            if dom_type(kk) == 0   % TM -- contributes to Er only
                if isnan(pd.psiTM(li,mi,1)), continue; end
                Er_modes(:,:,kk,f) = sqrt(eta0)*pd.psiTM(li,mi,f)*conj(pd.jrFlm_r(:,:,li,mi,f))/(1j*pd.Rs(pd.index_r));
            else                    % TE -- contributes to Hr only
                if isnan(pd.psiTE(li,mi,1)), continue; end
                Hr_modes(:,:,kk,f) = (1/sqrt(eta0))*pd.psiTE(li,mi,f)*conj(pd.jrFlm_r(:,:,li,mi,f))/(1j*pd.Rs(pd.index_r));
            end
        end
    end
    Hr_total = sum(Hr_modes, 3, 'omitnan');
    Er_total = sum(Er_modes, 3, 'omitnan');
    clearvars Er_modes Hr_modes

    fignum = plot_field_comparison_polar(fignum, ...
        pd.theta/pi, Hr_total(:,index_phi_p,index_f), [], pd.Hptr(:,index_phi_p,pd.index_r,3,index_f), ...
        pd.phi/pi,   Hr_total(index_th_p,:,index_f),  [], pd.Hptr(index_th_p,:,pd.index_r,3,index_f), ...
        eta0, '$\eta_0|H_r|$ (V/m)', {"(dominant-modes) $H_r$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
    clearvars Hr_total

    fignum = plot_field_comparison_polar(fignum, ...
        pd.theta/pi, Er_total(:,index_phi_p,index_f), [], pd.Eptr(:,index_phi_p,pd.index_r,3,index_f), ...
        pd.phi/pi,   Er_total(index_th_p,:,index_f),  [], pd.Eptr(index_th_p,:,pd.index_r,3,index_f), ...
        1, '$|E_r|$ (V/m)', {"(dominant-modes) $E_r$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
    clearvars Er_total

    % --- T_lm Magnitude / Phase ---
    T_psi_p = zeros(n_dom, numf);
    denom_psi_p = sqrt(1 - abs(reshape(S_psi(p,p,:), 1, numf)).^2);
    for kk = 1:n_dom
        T_psi_p(kk,:) = reshape(S_psi(num_internal+kk,p,:), 1, numf) ./ denom_psi_p;
    end
    if pd.used_Q
        T_Q_p = zeros(n_dom, numf);
        denom_Q_p = sqrt(1 - abs(reshape(S_Q(p,p,:), 1, numf)).^2);
        for kk = 1:n_dom
            T_Q_p(kk,:) = reshape(S_Q(num_internal+kk,p,:), 1, numf) ./ denom_Q_p;
        end
    end

    fignum = fignum+1; figure(fignum);
    ax_S = gca;
    set(ax_S, 'LineStyleOrder', {'-','--',':','-.'}')
    if pd.used_Q
        yyaxis left
    end
    set(ax_S,'YColor',get(groot,'defaultAxesXColor'))
    hold on
    plot(Freq, db(reshape(S_psi(p,p,:),1,numf)), "LineStyle","--", "LineWidth",LW, ...
        "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{"+int_name(p)+int_name(p)+"}$")
    for snp_k = 1:n_dom
        kk = snp_order(snp_k);
        c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
        plot(Freq, db(T_psi_p(kk,:)), "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
    end
    hold off
    if pd.used_Q
        ylabel("$|T_{lm}|$ dB ($\psi$)","FontSize",FS)
    else
        ylabel("$|T_{lm}|$ dB","FontSize",FS)
    end

    % Q-derived (right) axis only when this port's Q was read; otherwise
    % T_Q_p is zero and db() would be -Inf.
    if pd.used_Q
        yyaxis right
        set(ax_S,'YColor',get(groot,'defaultAxesXColor'))
        hold on
        plot(Freq, db(reshape(S_Q(p,p,:),1,numf)), "LineStyle","--", "LineWidth",LW+0.5, ...
            "Color",get(groot,'defaultAxesXColor'), "HandleVisibility","off")
        for snp_k = 1:n_dom
            kk = snp_order(snp_k);
            c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
            plot(Freq, db(T_Q_p(kk,:)), "LineStyle",":", "LineWidth",LW+0.5, "Color",c, "HandleVisibility","off")
        end
        hold off
        ylabel("$|T_{lm}|$ dB ($Q$)","FontSize",FS)
    end

    if pd.used_Q
        yyaxis left;  ylim([round(min(max(db(T_psi_p(:,fband_idx)), [], 2)))-5 0])
        yyaxis right; ylim([round(min(max(db(T_Q_p(:,fband_idx)),   [], 2)))-5 0])
    else
        ylim([round(min(max(db(T_psi_p(:,fband_idx)), [], 2)))-5 0])
    end

    if pd.used_Q
        title("Accepted Power Coupling $T_{lm} = S_{lm}/\sqrt{1-|S_{"+int_name(p)+int_name(p)+"}|^2}$ ("+int_name(p)+"): $\psi$ (solid) vs $Q$ (dotted)","FontSize",FS)
    else
        title("Accepted Power Coupling $T_{lm} = S_{lm}/\sqrt{1-|S_{"+int_name(p)+int_name(p)+"}|^2}$ ("+int_name(p)+")","FontSize",FS)
    end
    legend("FontSize",FS)
    xlabel("Frequency (Hz)","FontSize",FS)
    grid on; xlim([fmin,fmax])

    fignum = fignum+1; figure(fignum);
    ax_Sph = gca;
    set(ax_Sph, 'LineStyleOrder', {'-','--',':','-.'}')
    if pd.used_Q
        yyaxis left
    end
    set(ax_Sph,'YColor',get(groot,'defaultAxesXColor'))
    hold on
    plot(Freq, angle(reshape(S_psi(p,p,:),1,numf))*180/pi, "LineStyle","--", "LineWidth",LW, ...
        "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{"+int_name(p)+int_name(p)+"}$")
    for snp_k = 1:n_dom
        kk = snp_order(snp_k);
        c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
        plot(Freq, angle(T_psi_p(kk,:))*180/pi, "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
    end
    hold off
    if pd.used_Q
        ylabel("$\angle T_{lm}$ deg ($\psi$)","FontSize",FS)
    else
        ylabel("$\angle T_{lm}$ deg","FontSize",FS)
    end

    if pd.used_Q
        yyaxis right
        set(ax_Sph,'YColor',get(groot,'defaultAxesXColor'))
        hold on
        plot(Freq, angle(reshape(S_Q(p,p,:),1,numf))*180/pi, "LineStyle","--", "LineWidth",LW+0.5, ...
            "Color",get(groot,'defaultAxesXColor'), "HandleVisibility","off")
        for snp_k = 1:n_dom
            kk = snp_order(snp_k);
            c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
            plot(Freq, angle(T_Q_p(kk,:))*180/pi, "LineStyle",":", "LineWidth",LW+0.5, "Color",c, "HandleVisibility","off")
        end
        hold off
        ylabel("$\angle Q_{lm}$ deg ($Q$)","FontSize",FS)
        yyaxis left;  ylim([-180 180])
        yyaxis right; ylim([-180 180])
    else
        ylim([-180 180])
    end

    if pd.used_Q
        title("Accepted Power Coupling Phase $\angle T_{lm}$ ("+int_name(p)+"): $\psi$ (solid) vs $Q$ (dotted)","FontSize",FS)
    else
        title("Accepted Power Coupling Phase $\angle T_{lm}$ ("+int_name(p)+")","FontSize",FS)
    end
    legend("FontSize",FS)
    xlabel("Frequency (Hz)","FontSize",FS)
    grid on; xlim([fmin,fmax])

    % --- Radiated Power Comparison (this port's own psi/Q/Pr_NF and S) ---
    Pr_PSI_p = zeros(numf,1); Pr_PSIdom_p = zeros(numf,1);
    for f = 1:numf
        Pr_PSI_p(f) = (1/2) * sum(abs(pd.psiTM(:,:,f)).^2 + abs(pd.psiTE(:,:,f)).^2, 'all', 'omitnan');
        psiTM_f = pd.psiTM(:,:,f); psiTE_f = pd.psiTE(:,:,f);
        Pr_PSIdom_p(f) = (1/2) * (sum(abs(psiTM_f(mask_TM_dom)).^2, 'omitnan') + sum(abs(psiTE_f(mask_TE_dom)).^2, 'omitnan'));
    end

    if pd.used_Q
        Pr_Q_p = (1/2) * 8*pi * sum(abs(pd.Q).^2, 1).';

        Pr_Qdom_p = zeros(numf,1);
        for f = 1:numf
            Pr_Qdom_f = 0;
            for kk = 1:n_dom
                li = dom_lidx(kk); mi = dom_midx(kk);
                if dom_type(kk) == 0
                    j_pos = SphModeIndex(lval(li),  mval(mi), 'TM');
                    Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(pd.Q(j_pos,f))^2;
                    if mval(mi) ~= 0
                        j_neg = SphModeIndex(lval(li), -mval(mi), 'TM');
                        Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(pd.Q(j_neg,f))^2;
                    end
                else
                    j_pos = SphModeIndex(lval(li),  mval(mi), 'TE');
                    Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(pd.Q(j_pos,f))^2;
                    if mval(mi) ~= 0
                        j_neg = SphModeIndex(lval(li), -mval(mi), 'TE');
                        Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(pd.Q(j_neg,f))^2;
                    end
                end
            end
            Pr_Qdom_p(f) = Pr_Qdom_f;
        end
    else
        % No Q for this port: NaN traces
        Pr_Q_p    = nan(numf,1);
        Pr_Qdom_p = nan(numf,1);
    end

    Pr_Spsi_p = zeros(numf,1);
    for f = 1:numf
        % Radiated power = power into the mode ports only (power into other
        % internal ports is delivered to loads). The (1+S_pp) factor
        % cancels the one in S_psi(mode,p,f).
        port_voltage_factor_p = abs(int_V0(p) / (1 + S_psi(p,p,f)))^2;
        Pr_Spsi_p(f) = (1/2) * Y0 * port_voltage_factor_p * 2^abs(even_odd_excitation) * ...
                     sum(abs(S_psi(num_internal+1:end,p,f)).^2);
    end
    if pd.used_Q
        Pr_SQ_p = zeros(numf,1);
        for f = 1:numf
            port_voltage_factor_p = abs(int_V0(p) / (1 + S_Q(p,p,f)))^2;
            Pr_SQ_p(f) = (1/2) * Y0 * port_voltage_factor_p * 2^abs(even_odd_excitation) * ...
                         sum(abs(S_Q(num_internal+1:end,p,f)).^2);
        end
    else
        Pr_SQ_p = nan(numf,1);   % no Q for this port
    end

    fignum = fignum+1; figure_P = figure(fignum);
    hold on
    plot(Freq, real(pd.Pr_NF), "LineWidth",LW)
    plot(Freq, Pr_PSI_p,       "--", "LineWidth",LW+1)
    plot(Freq, Pr_PSIdom_p,    "s",  "LineWidth",LW+1, "MarkerSize",6)
    plot_labels_p = {"NF Poynting","$|\psi|^2$ all","$|\psi|^2$ dom"};
    if pd.used_Q
        plot(Freq, Pr_Q_p,    "--", "LineWidth",LW)
        plot(Freq, Pr_Qdom_p, "^",  "LineWidth",LW,   "MarkerSize",3)
        plot_labels_p = [plot_labels_p, {"$|Q|^2$ all","$|Q|^2$ dom"}];
    end
    plot(Freq, Pr_Spsi_p,      "o",  "LineWidth",LW,   "MarkerSize",4)
    plot_labels_p = [plot_labels_p, {"$|S|^2$ $\psi$"}];
    if pd.used_Q
        plot(Freq, Pr_SQ_p, "d",  "LineWidth",LW,   "MarkerSize",4)
        plot_labels_p = [plot_labels_p, {"$|S|^2$ $Q$"}];
    end
    hold off
    legend(plot_labels_p{:},"FontSize",FS)
    title("Radiated Power ("+int_name(p)+"): Different Methods","FontSize",FS)
    ylabel("Watts","FontSize",FS); xlabel("Freq (Hz)","FontSize",FS)
    grid on; xlim([fmin,fmax])

    % Relative L2 error over [fmin, fmax] of each method against the
    % near-field Poynting power (reference). Q methods only if Q was read.
    Pr_NF_band_p  = real(pd.Pr_NF(fband_idx));
    Pr_ref_norm_p = norm(Pr_NF_band_p);

    L2_err_PSI_p    = norm(Pr_PSI_p(fband_idx)    - Pr_NF_band_p) / Pr_ref_norm_p;
    L2_err_PSIdom_p = norm(Pr_PSIdom_p(fband_idx) - Pr_NF_band_p) / Pr_ref_norm_p;
    L2_err_Spsi_p   = norm(Pr_Spsi_p(fband_idx)   - Pr_NF_band_p) / Pr_ref_norm_p;

    fprintf('\nRadiated power L2 error vs NF Poynting (%s, %.3g-%.3g GHz band):\n', int_name(p), fmin/1e9, fmax/1e9);
    fprintf('  |psi|^2 all : %6.2f%%\n', 100*L2_err_PSI_p);
    fprintf('  |psi|^2 dom : %6.2f%%\n', 100*L2_err_PSIdom_p);
    if pd.used_Q
        L2_err_Q_p    = norm(Pr_Q_p(fband_idx)    - Pr_NF_band_p) / Pr_ref_norm_p;
        L2_err_Qdom_p = norm(Pr_Qdom_p(fband_idx) - Pr_NF_band_p) / Pr_ref_norm_p;
        fprintf('  |Q|^2 all   : %6.2f%%\n', 100*L2_err_Q_p);
        fprintf('  |Q|^2 dom   : %6.2f%%\n', 100*L2_err_Qdom_p);
    end
    fprintf('  |S|^2 psi   : %6.2f%%\n', 100*L2_err_Spsi_p);
    if pd.used_Q
        L2_err_SQ_p = norm(Pr_SQ_p(fband_idx) - Pr_NF_band_p) / Pr_ref_norm_p;
        fprintf('  |S|^2 Q     : %6.2f%%\n', 100*L2_err_SQ_p);
    end
end

%% --- Export Spherical Wave S-parameters to Touchstone File ---
% Ports 1..num_internal: internal FEKO ports (Feed, load1, ...);
% ports num_internal+1..n_tot: dominant modes in snp_order.
% snp_dir is set in USER SETTINGS.
if ~exist(snp_dir, 'dir'), mkdir(snp_dir); end

if snp_use_Q
    S_snp_base = S_Q;
    snp_label  = 'Q';
else
    S_snp_base = S_psi;
    snp_label  = 'psi';
end

snp_port_order = [1:num_internal, snp_order(:).' + num_internal];
S_snp = S_snp_base(snp_port_order, snp_port_order, :);

% Reject unknown dc_resistance_method values (they would otherwise fall
% through to the 'Rmin' behaviour silently).
valid_dc_methods = {'parallel_loads', 'series_loads', 'Rmin', 'user'};
if ~ismember(dc_resistance_method, valid_dc_methods)
    error('general_vector_spherical_modes_loadports: dc_resistance_method = ''%s'' is not recognized. Valid values are: %s.', ...
        dc_resistance_method, strjoin(valid_dc_methods, ', '));
end

% DC (0 Hz) point, per dc_resistance_method (see Simulation Parameters).
% Mode ports stay at S_DC = 0 (matched at DC).
S_DC = zeros(n_tot, n_tot);

% Load resistance of each internal port (DC is purely resistive): Z0 for
% the Feed; for load ports the first available of Z_field_load,
% Z_sparam_ref, Z0_ref_own as reported by FEKO.
Z_load_all = nan(1, num_internal);
Z_load_all(1) = Z0;
for k_sec = 1:n_sec
    if ~isnan(secondary_ports(k_sec).Z_field_load)
        Z_load_all(k_sec+1) = real(secondary_ports(k_sec).Z_field_load);
    elseif ~isnan(secondary_ports(k_sec).Z_sparam_ref)
        Z_load_all(k_sec+1) = secondary_ports(k_sec).Z_sparam_ref;
    elseif ~isnan(secondary_ports(k_sec).Z0_ref_own)
        Z_load_all(k_sec+1) = secondary_ports(k_sec).Z0_ref_own;
    end
end

% --- Feed ---
Zin_feed = Zant(1);
if ismember(dc_resistance_method, {'parallel_loads', 'series_loads'}) && isempty(Z_load_all(2:end))
    % Feed only, no other DC path: use the reactive limit, S = -1 if
    % inductive (short as f->0), +1 if capacitive (open as f->0).
    if imag(Zin_feed) > 0
        dc_reactance_type = 'inductive';
        S11_DC = -1;
    else
        dc_reactance_type = 'capacitive';
        S11_DC = +1;
    end
    fprintf('Spherical modes: DC point (Feed) — %s, no other internal port to combine with (isolated at DC), S_DC = %.4f\n', ...
        dc_reactance_type, S11_DC);
else
    Rin_DC_feed = determine_dc_resistance(dc_resistance_method, Zin_feed, Z_load_all(2:end));
    S11_DC = (Rin_DC_feed - Z0) / (Rin_DC_feed + Z0);
end
S_DC(1,1) = S11_DC;

% --- Secondary ports (diagonal DC entries) ---
for k_sec = 1:n_sec
    p_idx = int_sparam_idx(k_sec+1);
    if isnan(p_idx)
        continue   % port not in the S-parameter configuration: S_DC stays 0
    end
    Sk1 = S(p_idx, p_idx, 1);
    Zin_k1 = Z0 * (1 + Sk1) / (1 - Sk1);   % Zin at the lowest retained frequency

    other_idx = setdiff(1:num_internal, k_sec+1);   % always contains the Feed
    Rin_DC = determine_dc_resistance(dc_resistance_method, Zin_k1, Z_load_all(other_idx));
    Sk_DC = (Rin_DC - Z0) / (Rin_DC + Z0);
    S_DC(k_sec+1, k_sec+1) = Sk_DC;
end

% --- DC coupling between every pair of internal ports ---
% Resistive N-port with port resistances Ri (R1 = Z0 for the Feed) all in
% series, or conductances Gi = 1/Ri all in parallel, at a common node:
%   series:   S_nm = 2*sqrt(Rn*Rm) / sum(Ri)
%   parallel: S_nm = 2*sqrt(Gn*Gm) / sum(Gi)
% For 'Rmin'/'user' (or missing load data) each pair (n<m) is set from
% the lowest-frequency S_nm: |S_nm| with sign -1 if imag(S_nm) > 0, else
% +1 (DC S-parameters are real), or entered by the user; then mirrored
% (reciprocal network).
if strcmp(dc_resistance_method, 'series_loads') && ~any(isnan(Z_load_all))
    Total_R = sum(Z_load_all);
    for n = 1:num_internal
        for m = 1:num_internal
            if n == m, continue; end
            S_DC(n,m) = 2 * sqrt(Z_load_all(n) * Z_load_all(m)) / Total_R;
        end
    end
elseif strcmp(dc_resistance_method, 'parallel_loads') && ~any(isnan(Z_load_all))
    G_load_all = 1 ./ Z_load_all;
    Total_G = sum(G_load_all);
    for n = 1:num_internal
        for m = 1:num_internal
            if n == m, continue; end
            S_DC(n,m) = 2 * sqrt(G_load_all(n) * G_load_all(m)) / Total_G;
        end
    end
else
    for n = 1:num_internal
        p_n = int_sparam_idx(n);
        if isnan(p_n), continue; end
        for m = (n+1):num_internal
            p_m = int_sparam_idx(m);
            if isnan(p_m), continue; end
            Snm_lowfreq = S(p_n, p_m, 1);
            if strcmp(dc_resistance_method, 'user')
                Snm_DC = input(sprintf("enter value for DC coupling S(%s,%s): ", int_name(n), int_name(m)));
            elseif imag(Snm_lowfreq) > 0
                Snm_DC = -abs(Snm_lowfreq);
            else
                Snm_DC = +abs(Snm_lowfreq);
            end
            S_DC(n,m) = Snm_DC;
            S_DC(m,n) = Snm_DC;
        end
    end
end

% Export only frequencies up to fmax.
export_freq_idx = find(Freq <= fmax);
Freq_export = Freq(export_freq_idx);
S_snp_export = S_snp(:,:,export_freq_idx);

snp_path = fullfile(snp_dir, "S_dom_"+snp_label+"_"+antenna_name+".s"+num2str(n_tot)+"p");
rfwrite(cat(3, S_DC, S_snp_export), [0, Freq_export], snp_path, "ReferenceResistance", Z0)

% ADS-compatible port names (no spaces/special characters): Feed, load1,
% ..., then one per mode, e.g. TM_L1_M0 or TE_L2_Mneg1 (m = -1).
port_names = strings(1, n_tot);
port_names(1) = "Feed";
for k_sec = 1:n_sec
    port_names(k_sec+1) = int_name(k_sec+1);   % "load1", "load2", ...
end
for snp_k = 1:n_dom
    kk = snp_order(snp_k);
    li = dom_lidx(kk); mi = dom_midx(kk);
    type_str = 'TM'; if dom_type(kk) == 1, type_str = 'TE'; end
    if mval(mi) < 0
        m_str = sprintf('Mneg%d', abs(mval(mi)));
    else
        m_str = sprintf('M%d', mval(mi));
    end
    port_names(num_internal+snp_k) = sprintf('%s_L%d_%s', type_str, lval(li), m_str);
end

if apply_ka_deembedding
    ka_factor_str = "exp(-jka) applied";
else
    ka_factor_str = "1 (no de-embedding)";
end
comment_lines = "! Touchstone file version 1.1" + newline;
comment_lines = comment_lines + "! Written by " + string(mfilename) + ".m (" + snp_label + "-derived)," +...
                " a="+num2str(a,3)+"m, Z0="+num2str(Z0,5)+"ohms, ka_factor="+ka_factor_str + newline;
comment_lines = comment_lines + "! Date : " + string(datetime('now')) + newline;

% Port names as "! Port[n] name" comment lines (the form ADS writes in
% Touchstone 1.1). ADS shows them only if no other comment lines (besides
% the three above) precede the '#' line; a [Port Names] keyword is not
% Touchstone 1.1 and ADS rejects it.
for pn = 1:length(port_names)
    comment_lines = comment_lines + sprintf('! Port[%d] %s\n', pn, port_names(pn));
end

fid_in  = fopen(snp_path, 'r');
existing = fread(fid_in, '*char')';
fclose(fid_in);

fid_out = fopen(snp_path, 'w');
fprintf(fid_out, '%s', comment_lines);
fprintf(fid_out, '%s', existing);
fclose(fid_out);
clearvars existing

fprintf('Spherical modes: SnP port names written: %s\n', strjoin(port_names, ', '));

%% ================= Local functions ================= %%

function [Rin_DC] = determine_dc_resistance(method, Zin_lowfreq, Z_load_others)
%DETERMINE_DC_RESISTANCE  DC (0 Hz) input resistance of one internal port.
%   method        - 'parallel_loads' | 'series_loads' | 'Rmin' | 'user'
%   Zin_lowfreq   - port's complex Zin at the lowest simulated frequency
%                   ('Rmin' uses its real part)
%   Z_load_others - real load resistances of the other internal ports,
%                   combined as 1/sum(1./Ri) ('parallel_loads') or
%                   sum(Ri) ('series_loads')
%   Rin_DC        - DC input resistance (Ohm); prompted for 'user'.
% Unknown methods fall back to 'Rmin'.
Rmin = real(Zin_lowfreq);
switch method
    case 'user'
        Rin_DC = input("enter value DC input resistance: ");
        return
    case 'Rmin'
        Rin_DC=Rmin;
    case 'parallel_loads'
        Rin_DC = 1 / sum(1 ./ Z_load_others);
    case 'series_loads'
        Rin_DC = sum(Z_load_others);
   otherwise
            Rin_DC=Rmin;
end

end
