%GENERAL_VECTOR_SPHERICAL_MODES_FEKO  Single-port FEKO antenna -> spherical-mode network model.
%
% Extracts the spherical-mode coefficients psiTM/psiTE of a FEKO antenna
% with one driven port, ranks the (l,m,TE/TM) modes by accepted-power
% coupling |T_lm| and exports an (n_dom+1)-port Touchstone file: port 1 is
% the antenna feed, ports 2..n_dom+1 are the n_dom dominant spherical modes.
%
% Inputs (FEKO exports, all named <fileprefix>...):
%   .out                              port parameters, input impedance Zant,
%                                     internal loads (read by ReadFEKOout)
%   _<near_field_request>.efe/.hfe    near-field E/H on a sphere (required)
%   _<far_field_request>.ffe          far field (optional, use_FF_method)
%   _<far_field_request>.sph          spherical-mode Q coefficients
%                                     (optional, use_Q_method)
%
% Outputs:
%   <snp_dir>/S_dom_<psi|Q>_<antenna_name>.s<n_dom+1>p with ADS-compatible
%   port names; figures of impedance, coefficients, field reconstruction,
%   T_lm and radiated power.

clear
close all
clc

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
% FEKO project file prefix (folder + antenna name). Uncomment one antenna.
antenna_name="loop_loadport";
fileprefix="C:\path\to\FEKO\loop\"+antenna_name;
% antenna_name="bowtie_edgeport_45degrees";
% fileprefix="C:\path\to\FEKO\bowtie\"+antenna_name;
% antenna_name="dipole_cylinder_edgeport_symmetry_full";
% fileprefix="C:\path\to\FEKO\dipole\"+antenna_name;

% FEKO field-request names: files are read from <fileprefix>_<name>.<ext>.
% far_field_request = "" skips the far-field (.ffe) and Q (.sph) methods.
near_field_request = "NearField1";
far_field_request  = "FarField1";

% Touchstone output folder.
snp_dir = fullfile(repo_root, 'output', 'SnP');
% =========================

%% --- Physical Constants ---
c0  = physconst('lightspeed');
mu0 = 4*pi*1e-7;          % Permeability of free space (H/m)
eps0 = 1/c0^2/mu0;        % Permittivity of free space (F/m)
eta0 = sqrt(mu0/eps0);    % Characteristic impedance of free space (Ohm)

%% --- Simulation Parameters ---
f0     = 4e9;         % Reference frequency for single-frequency plots (Hz)
fmin   = 0.1e9;       % Minimum frequency for plot x-axis (Hz)
fmax   = 6e9;         % Maximum frequency for plot x-axis (Hz)
lmax   = 11;           % Maximum spherical mode order to compute
n_dom  = 13;           % Number of dominant modes to identify and plot

% Coefficient source for S-parameters/ranking/SnP export; set automatically
% below: psi whenever near- or far-field data is available, Q only when it
% is the sole source.
snp_use_Q = false;

% Vector method: project the near-field E/H onto the vector spherical wave
% functions (compute_mode_coeffs_nearfield). The near field itself is always
% read (it is the reference trace in the field plots and gives Pr_NF).
use_vector_method = true;

% Far-field method: add a psi estimate from the .ffe far field
% (compute_mode_coeffs_farfield, which uses its own theta/phi grid).
use_FF_method = true;

% Q method: read FEKO's spherical-mode coefficients (.sph), converted to the
% psi convention as QTMconv/QTEconv. Q enters the psi average only when no
% field data (vector or FF) is available; otherwise it is used only for
% comparison plots. false skips the .sph file and everything that depends
% on it. At least one of use_vector_method, use_FF_method, use_Q_method
% must be true.
use_Q_method = true;

% SnP port ordering (ports 2..n_dom+1):
%   'peak_freq'  - by frequency of peak |T_lm|, ascending
%   'mode_index' - by FEKO mode index j = 2*(l*(l+1)+m-1)+s, s=1 TE, s=2 TM
%snp_order_method = 'peak_freq';
snp_order_method = 'mode_index';

% true: the far-field phase j^(l+1) is absorbed into psiTM/psiTE/QTMconv/
% QTEconv (via j_l_1), so the coupling formula uses jl1 = 1. false: psi/Q
% exclude it and the coupling formula applies jl1 = (1j)^(l+1).
use_jl1_convention = true;

% true: apply exp(-j*k*a) to move the S-parameter reference plane from the
% feed to r = a (antenna apothem). false: reference stays at the feed.
apply_ka_deembedding = false;

% DC (0 Hz) point of the SnP export:
%   'auto' - internal-load DC-conduction model if the .out file has an
%            internal load with nonzero resistance, else 'Rmin'
%   'Rmin' - reactive-limit approximation from the lowest-frequency Zant
%   'user' - prompt for S11 at DC
dc_resistance_method = 'auto';
%% --- Plot Settings ---
LW = 2;    % LineWidth
FS = 14;   % FontSize

set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% --- FEKO Import: Read Port Parameters from .out File ---
% V0: open-circuit voltage (complex), Z0: port reference impedance,
% R0: port location, symmetry_planes.physical_gnd: PEC ground plane flag,
% even_odd_excitation: log2(ports) if mode symmetry detected, a: antenna
% apothem (m), Zant: input impedance, internal_loads: loaded, unexcited ports
[V0, Z0, R0, symmetry_planes, even_odd_excitation, a, Zant, Freq, port_pol, internal_loads] = ReadFEKOout(fileprefix+".out");

Zant=transpose(Zant(:,1));Freq=transpose(Freq);
indicies_negR=find(real(Zant)<0);
Zindx=find(~(real(Zant)<0));Freq=Freq(Zindx);
Zant(find(real(Zant)<1e-4))=Zant(find(real(Zant)<1e-4))+1e-4;
Zant=Zant(Zindx);

if abs(Z0-eta0) > 0.01
    error("antenna port impedance should be eta0 for now")
else
    Z0 = eta0;Y0 = 1/Z0;
end


numf = length(Freq);

if fmax>Freq(numf)
    fmax=Freq(numf);
end
if fmin<Freq(1)
    fmin=Freq(1);
end

%% --- Compute Spherical-Mode Coefficients ---
% Near-field, far-field and Q extraction are separate functions; their
% results are averaged here so each source can also be inspected on its
% own. The mode-index grids depend only on lmax and the phase conventions
% and are shared by all three.
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

% The near field is always read: Eptr/Hptr and Pr_NF are used downstream
% even when the vector-method projection is disabled.
nf = compute_mode_coeffs_nearfield(fileprefix, near_field_request, Freq, Zindx, a, ...
    lval, mval, L, M, use_vector_method, l_factor);
if ~nf.available
    error('general_vector_spherical_modes: no near-field data could be read. Check that the near-field files <fileprefix>_%s.efe/.hfe exist.', near_field_request)
end

ff = compute_mode_coeffs_farfield(fileprefix, far_field_request, Freq, Zindx, lval, mval, L, M, use_FF_method);
qd = compute_mode_coeffs_Q(fileprefix, far_field_request, Freq, Zindx, lmax, lval, mval, L, M, ...
    symmetry_planes, j_l_1, CSphase, use_Q_method);

use_vector_method = nf.used_vector;
use_FF_method     = ff.available;
use_Q_method      = qd.available;
has_field_data    = use_vector_method || use_FF_method;

% Use Q for S-parameters/ranking/SnP export only when it is the sole source.
snp_use_Q = use_Q_method && ~has_field_data;

% Average psi over the available field methods; Q is included only when
% there is no field data.
psi_stack_TM = []; psi_stack_TE = [];
if use_vector_method
    psi_stack_TM = cat(4, psi_stack_TM, nf.psiTM);
    psi_stack_TE = cat(4, psi_stack_TE, nf.psiTE);
end
if use_FF_method
    psi_stack_TM = cat(4, psi_stack_TM, ff.psiTM);
    psi_stack_TE = cat(4, psi_stack_TE, ff.psiTE);
end
if use_Q_method && ~has_field_data
    psi_stack_TM = cat(4, psi_stack_TM, qd.QTMconv);
    psi_stack_TE = cat(4, psi_stack_TE, qd.QTEconv);
end
if isempty(psi_stack_TM)
    error('general_vector_spherical_modes: no spherical-mode data could be extracted. At least one of use_vector_method, use_FF_method, use_Q_method must be true and produce usable data (see the messages above).')
end
psiTM = mean(psi_stack_TM, 4, 'omitnan');
psiTE = mean(psi_stack_TE, 4, 'omitnan');
clearvars psi_stack_TM psi_stack_TE

QTMconv = qd.QTMconv; QTEconv = qd.QTEconv;
Q = qd.Q; freqs_Q = qd.freqs_Q; mmax = qd.mmax; num_modes = qd.num_modes;

theta = nf.theta; phi = nf.phi; Rs = nf.Rs;
THETA = nf.THETA; PHI = nf.PHI;
index_r = nf.index_r;
[~,index_f]   = min(abs(Freq-f0));
[~,index_phi] = min(abs(phi - 0));       % phi index closest to 0 rad
[~,index_th]  = min(abs(theta - pi/2));  % theta index closest to pi/2
k = nf.k; kr = nf.kr;
numth = length(theta); numph = length(phi);
Pr_NF = nf.Pr_NF;
jrFlm_r = nf.jrFlm_r;
Eptr = nf.Eptr; Hptr = nf.Hptr;
clearvars nf ff qd

if snp_use_Q
    fprintf('general_vector_spherical_modes: only Q data is available — using Q-derived coefficients for S-parameters/ranking/SnP export.\n');
end

%% --- Build Spherical Wave S-parameters for EVERY (l,m,type), for Ranking ---
% Every mode is ranked by its accepted-power coupling
% T_lm = S_{mode,1}/sqrt(1-|S_00|^2), which, unlike raw |psi|, accounts for
% how much power the antenna actually accepts at each frequency.
%   S_{mode,1} = sqrt(Z0/2^|even_odd_excitation|) * (1+S_00) * (value/V0) * exp(-jka)/j^(l+1)
S00_all = (Zant - Z0) ./ (Zant + Z0);                    % (1 x numf)
port_factor = sqrt(Z0 / 2^abs(even_odd_excitation));
if apply_ka_deembedding
    ka_factor = exp(-1j * k(:).' * a);                   % (1 x numf)
    ka_factor_r = reshape(ka_factor, 1, 1, numf);
else
    ka_factor = 1;                                       % no de-embedding phase
    ka_factor_r = 1;
end
if use_jl1_convention
    jl1_l = ones(L,1);                                   % phase already absorbed into psiTM/psiTE/QTMconv/QTEconv via j_l_1
else
    jl1_l = (1j).^(lval(:) + 1);                         % (L x 1)
end

one_plus_S00_r = reshape(1 + S00_all, 1, 1, numf);
jl1_r          = reshape(jl1_l, L, 1, 1);                % broadcasts over M and numf

psiTM_safe = psiTM; psiTM_safe(isnan(psiTM_safe)) = 0;
psiTE_safe = psiTE; psiTE_safe(isnan(psiTE_safe)) = 0;

S_TM_all   = port_factor .* one_plus_S00_r .* (psiTM_safe ./ V0(1)) .* ka_factor_r ./ jl1_r;   % (L x M x numf)
S_TE_all   = port_factor .* one_plus_S00_r .* (psiTE_safe ./ V0(1)) .* ka_factor_r ./ jl1_r;

T_denom_all = sqrt(1 - abs(S00_all).^2);                 % (1 x numf)
T_TM_all   = S_TM_all   ./ reshape(T_denom_all, 1, 1, numf);   % (L x M x numf)
T_TE_all   = S_TE_all   ./ reshape(T_denom_all, 1, 1, numf);

% Q-derived coupling (only needed when Q was read).
if use_Q_method
    S_TM_Q_all = port_factor .* one_plus_S00_r .* (QTMconv ./ V0(1)) .* ka_factor_r ./ jl1_r;
    S_TE_Q_all = port_factor .* one_plus_S00_r .* (QTEconv ./ V0(1)) .* ka_factor_r ./ jl1_r;
    T_TM_Q_all = S_TM_Q_all ./ reshape(T_denom_all, 1, 1, numf);
    T_TE_Q_all = S_TE_Q_all ./ reshape(T_denom_all, 1, 1, numf);
end

%% --- Identify Dominant Modal Coefficients across TE and TM jointly ---
% TE and TM modes are pooled and the n_dom largest by max |T_lm| within
% [fmin, fmax] are selected (from psi or Q according to snp_use_Q).
% fmax indices refer to fband_idx and are converted back to Freq indices.

fband_idx = find(Freq >= fmin & Freq <= fmax);
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

% Pool into one list [TM; TE]; type_tag: 0 = TM, 1 = TE. (:) keeps every
% pooled quantity a column vector.
pool_mag     = [rank_TM_mag(:);  rank_TE_mag(:)  ];
pool_fmax    = [rank_TM_fmax(:); rank_TE_fmax(:) ];
pool_type    = [zeros(L*M,1); ones(L*M,1)];
pool_lmidx   = [(1:L*M).';   (1:L*M).'  ];

% Select the top n_dom by magnitude, then order them by peak frequency.
[~, sort_pool] = sort(pool_mag, 'descend');
dom_pool_mag   = sort_pool(1:n_dom);
dom_peak_freqs = Freq(pool_fmax(dom_pool_mag));
[~, freq_order] = sort(dom_peak_freqs(:), 'ascend');
dom_pool = dom_pool_mag(freq_order);

% Decode each dominant entry into type, (lindx, mindx)
dom_type = pool_type(dom_pool);      % 0=TM, 1=TE, length n_dom
dom_lm   = pool_lmidx(dom_pool);     % linear index into (L x M), length n_dom
[dom_lidx, dom_midx] = ind2sub([L, M], dom_lm);

% Separate TM and TE index lists (used for the radiated-power masks).
dom_TM_lidx = dom_lidx(dom_type == 0);
dom_TM_midx = dom_midx(dom_type == 0);
dom_TE_lidx = dom_lidx(dom_type == 1);
dom_TE_midx = dom_midx(dom_type == 1);
n_dom_TM = length(dom_TM_lidx);
n_dom_TE = length(dom_TE_lidx);

% Report dominant modes ranked jointly — dB relative to the global dominant
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
% S_psi and S_Q are (n_dom+1) x (n_dom+1) x numf.
% Index 1: antenna port (S(1,1,:) = S_00); indices 2..n_dom+1: dominant
% modes in rank order. Values are slices of S_TM_all/S_TE_all above.

S_psi = zeros(n_dom+1, n_dom+1, numf);
S_Q   = zeros(n_dom+1, n_dom+1, numf);
S_psi(1,1,:) = S00_all;
S_Q(1,1,:)   = S00_all;

for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0   % TM
        coupling_psi = reshape(S_TM_all(li,mi,:),   1, numf);
        if use_Q_method
            coupling_Q = reshape(S_TM_Q_all(li,mi,:), 1, numf);
        end
    else                   % TE
        coupling_psi = reshape(S_TE_all(li,mi,:),   1, numf);
        if use_Q_method
            coupling_Q = reshape(S_TE_Q_all(li,mi,:), 1, numf);
        end
    end
    S_psi(kk+1, 1, :) = coupling_psi;
    S_psi(1, kk+1, :) = coupling_psi;
    if use_Q_method
        S_Q(kk+1, 1, :) = coupling_Q;
        S_Q(1, kk+1, :) = coupling_Q;
    end
end
% Sanitise: replace any remaining NaN with 0 before export
S_psi(isnan(S_psi)) = 0;
S_Q(isnan(S_Q))     = 0;

%% --- Mode labels, T_lm, and SnP port ordering ---
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

% Accepted power coupling T_lm = S_{k+1,1} / sqrt(1 - |S_{1,1}|^2)
T_psi = zeros(n_dom, numf);
T_Q   = zeros(n_dom, numf);
denom_psi = sqrt(1 - abs(reshape(S_psi(1,1,:), 1, numf)).^2);
denom_Q   = sqrt(1 - abs(reshape(S_Q(1,1,:),   1, numf)).^2);
for kk = 1:n_dom
    T_psi(kk,:) = reshape(S_psi(kk+1,1,:), 1, numf) ./ denom_psi;
    T_Q(kk,:)   = reshape(S_Q(kk+1,1,:),   1, numf) ./ denom_Q;
end

% snp_order: ordering of the n_dom dominant modes for SnP ports 2..n_dom+1.
% Shared by T plots and SnP export so legend/colors match port assignments.
if snp_use_Q
    T_for_order = T_Q;
else
    T_for_order = T_psi;
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

%% --- Plots: Impedance, Coefficients, Reconstruction, T_lm, Power ---
fignum = 0;

%% --- Plot Antenna Impedance vs Frequency ---
fignum = fignum+1; figure(fignum)
plot(Freq,real(Zant),"Color","b","LineWidth",LW)
hold on
plot(Freq,imag(Zant),"Color","r","LineWidth",LW)
hold off
title("Impedance","FontSize",FS)
legend("$R$","$X$","FontSize",FS)
xlim([fmin,fmax])
ylim([-eta0,eta0])
grid on

% Combined TE + TM figure — all n_dom dominant modes in rank order
fignum = fignum+1; figure(fignum);
ax = gca;
set(ax, 'LineStyleOrder', {'-','--',':','-.'}')
clr = get(groot,'defaultAxesColorOrder');   % capture color table once
n_clr = size(clr, 1);
if use_Q_method
    yyaxis left
end
set(ax,'YColor',get(groot,'defaultAxesXColor'))
hold on
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    c = clr(mod(kk-1, n_clr)+1, :);
    if dom_type(kk) == 0   % TM
        lbl = sprintf('TM $l=%d,m=%d$', lval(li), mval(mi));
        data = abs(reshape(psiTM(li,mi,:), 1, numf));
    else                  % TE
        lbl = sprintf('TE $l=%d,m=%d$', lval(li), mval(mi));
        data = abs(reshape(psiTE(li,mi,:), 1, numf));
    end
    plot(Freq, data, "LineWidth",LW-0.5, "LineStyle","-", "Color",c, "DisplayName",lbl)
end
hold off
ylabel("$|\psi|$","FontSize",FS)

% Right axis: Q converted to the psi convention (only when Q was read).
ymax_vals = zeros(n_dom, 2);
if use_Q_method
    yyaxis right
    set(ax,'YColor',get(groot,'defaultAxesXColor'))
    hold on
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        c = clr(mod(kk-1, n_clr)+1, :);
        if dom_type(kk) == 0   % TM
            data = abs(reshape(QTMconv(li,mi,:), 1, numf));
        else                  % TE
            data = abs(reshape(QTEconv(li,mi,:), 1, numf));
        end
        plot(Freq, data, "LineWidth",LW+0.5, "LineStyle",":", "Color",c, "HandleVisibility","off")
    end
    hold off
    ylabel("$|Q|$ converted","FontSize",FS)
end

for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0
        ymax_vals(kk,1) = max(abs(reshape(psiTM(li,mi,:),  [],1)), [], 'omitnan');
        if use_Q_method
            ymax_vals(kk,2) = max(abs(reshape(QTMconv(li,mi,:), [],1)), [], 'omitnan');
        end
    else
        ymax_vals(kk,1) = max(abs(reshape(psiTE(li,mi,:),  [],1)), [], 'omitnan');
        if use_Q_method
            ymax_vals(kk,2) = max(abs(reshape(QTEconv(li,mi,:), [],1)), [], 'omitnan');
        end
    end
end
ymax_all = max(ymax_vals(:), [], 'omitnan');
if isnan(ymax_all) || ymax_all <= 0, ymax_all = 1; end
if use_Q_method
    yyaxis left;  ylim([0 ymax_all])
    yyaxis right; ylim([0 ymax_all])
else
    ylim([0 ymax_all])
end

legend("FontSize",FS)
if use_Q_method
    title("Dominant modes: $|\psi|$ (solid) vs $|Q|$ converted (dotted)","FontSize",FS)
else
    title("Dominant modes: $|\psi|$","FontSize",FS)
end
xlim([fmin,fmax]); xlabel("Frequency (Hz)","FontSize",FS)
grid on

% Phase plot
fignum = fignum+1; figure(fignum);
ax_ph = gca;
set(ax_ph, 'LineStyleOrder', {'-','--',':','-.'}')
if use_Q_method
    yyaxis left
end
set(ax_ph,'YColor',get(groot,'defaultAxesXColor'))
hold on
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    c = clr(mod(kk-1, n_clr)+1, :);
    if dom_type(kk) == 0
        lbl = sprintf('TM $l=%d,m=%d$', lval(li), mval(mi));
        data = angle(reshape(psiTM(li,mi,:), 1, numf)) * 180/pi;
    else
        lbl = sprintf('TE $l=%d,m=%d$', lval(li), mval(mi));
        data = angle(reshape(psiTE(li,mi,:), 1, numf)) * 180/pi;
    end
    plot(Freq, data, "LineWidth",LW-0.5, "LineStyle","-", "Color",c, "DisplayName",lbl)
end
hold off
ylabel("$\angle\psi$ (deg)","FontSize",FS)

if use_Q_method
    yyaxis right
    set(ax_ph,'YColor',get(groot,'defaultAxesXColor'))
    hold on
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        c = clr(mod(kk-1, n_clr)+1, :);
        if dom_type(kk) == 0
            data = angle(reshape(QTMconv(li,mi,:), 1, numf)) * 180/pi;
        else
            data = angle(reshape(QTEconv(li,mi,:), 1, numf)) * 180/pi;
        end
        plot(Freq, data, "LineWidth",LW+0.5, "LineStyle",":", "Color",c, "HandleVisibility","off")
    end
    hold off
    ylabel("$\angle Q$ converted (deg)","FontSize",FS)
    yyaxis left;  ylim([-180 180])
    yyaxis right; ylim([-180 180])
else
    ylim([-180 180])
end

legend("FontSize",FS)
if use_Q_method
    title("Dominant modes: $\angle\psi$ (solid) vs $\angle Q$ converted (dotted)","FontSize",FS)
else
    title("Dominant modes: $\angle\psi$","FontSize",FS)
end
xlim([fmin,fmax]); xlabel("Frequency (Hz)","FontSize",FS)
grid on

%% --- Reconstruct Near-Field Fields from Coefficients, and Plot ---
% Each field component is reconstructed from exactly the n_dom dominant
% (l,m,type) contributions, plotted, and its per-mode arrays cleared to
% limit memory. TM contributes to Er and TE to Hr only; both contribute to
% the transverse components.

field_title_base = "@r="+num2str(Rs(index_r)*1000,2)+"mm, "+ ...
                   num2str(Freq(index_f)/1e9,2)+"GHz, \#dom="+num2str(n_dom);
phi_label   = "$\phi="+num2str(phi(index_phi)*180/pi,3)+"^\circ$";
theta_label = "$\theta="+num2str(theta(index_th)*180/pi,3)+"^\circ$";

%% --- Er, Hr (plus the Q-based FEKO reconstruction used as their comparison trace) ---
Er_modes = zeros(numth,numph,n_dom,numf); Hr_modes = Er_modes;
for f = 1:numf
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        if dom_type(kk) == 0   % TM -- contributes to Er only
            if isnan(psiTM(li,mi,1)), continue; end
            Er_modes(:,:,kk,f) = sqrt(eta0)*psiTM(li,mi,f)*conj(jrFlm_r(:,:,li,mi,f))/(1j*Rs(index_r));
        else                    % TE -- contributes to Hr only
            if isnan(psiTE(li,mi,1)), continue; end
            Hr_modes(:,:,kk,f) = (1/sqrt(eta0))*psiTE(li,mi,f)*conj(jrFlm_r(:,:,li,mi,f))/(1j*Rs(index_r));
        end
    end
end

% Q-based comparison trace (only when Q was read): same formula and basis
% with QTMconv/QTEconv in place of psi, at the plotted frequency index_f
% only. QTMconv/QTEconv are zero (not NaN) for invalid (l,m). Empty traces
% are omitted by plot_field_comparison_polar.
if use_Q_method
    Er_Q = zeros(numth,numph); Hr_Q = zeros(numth,numph);
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        if dom_type(kk) == 0   % TM -- contributes to Er only
            Er_Q = Er_Q + sqrt(eta0)*QTMconv(li,mi,index_f)*conj(jrFlm_r(:,:,li,mi,index_f))/(1j*Rs(index_r));
        else                    % TE -- contributes to Hr only
            Hr_Q = Hr_Q + (1/sqrt(eta0))*QTEconv(li,mi,index_f)*conj(jrFlm_r(:,:,li,mi,index_f))/(1j*Rs(index_r));
        end
    end
    Hr_FEKO_th = Hr_Q(:,index_phi); Hr_FEKO_ph = Hr_Q(index_th,:);
    Er_FEKO_th = Er_Q(:,index_phi); Er_FEKO_ph = Er_Q(index_th,:);
    clearvars Er_Q Hr_Q
else
    Hr_FEKO_th = []; Hr_FEKO_ph = [];
    Er_FEKO_th = []; Er_FEKO_ph = [];
end

clearvars jrFlm_r   % only needed for Er_modes/Hr_modes above

Hr_total = sum(Hr_modes, 3, 'omitnan');
Er_total = sum(Er_modes, 3, 'omitnan');
clearvars Er_modes Hr_modes

fignum = plot_field_comparison_polar(fignum, ...
    theta/pi, Hr_total(:,index_phi,index_f), Hr_FEKO_th, Hptr(:,index_phi,index_r,3,index_f), ...
    phi/pi,   Hr_total(index_th,:,index_f),  Hr_FEKO_ph, Hptr(index_th,:,index_r,3,index_f), ...
    eta0, '$\eta_0|H_r|$ (V/m)', {"(dominant-modes) $H_r$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Hr_total Hr_FEKO_th Hr_FEKO_ph

fignum = plot_field_comparison_polar(fignum, ...
    theta/pi, Er_total(:,index_phi,index_f), Er_FEKO_th, Eptr(:,index_phi,index_r,3,index_f), ...
    phi/pi,   Er_total(index_th,:,index_f),  Er_FEKO_ph, Eptr(index_th,:,index_r,3,index_f), ...
    1, '$|E_r|$ (V/m)', {"(dominant-modes) $E_r$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Er_total Er_FEKO_th Er_FEKO_ph

%% --- Gradient basis for the dominant modes (transverse-field plots) ---
% Angular gradient of Ylm on the near-field grid, for the n_dom dominant
% modes only.
gradYlm_th_dom = zeros(numth,numph,n_dom); gradYlm_ph_dom = gradYlm_th_dom;
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    [gradYlm_th_dom(:,:,kk), gradYlm_ph_dom(:,:,kk)] = gradYlm(lval(li), mval(mi), THETA, PHI, have_CSphase);
end

%% --- H_phi ---
Hph_TM = zeros(numth,numph,n_dom,numf); Hph_TE = Hph_TM;
for f = 1:numf
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        if dom_type(kk) == 0   % TM
            if isnan(psiTM(li,mi,1)), continue; end
            Hph_TM(:,:,kk,f) = (1/sqrt(eta0))*psiTM(li,mi,f)*sphericalRicatti_h(lval(li),2,kr(f))/j_l_1(li) ...
                                *(-gradYlm_th_dom(:,:,kk)/Rs(index_r))  ;
        else                    % TE
            if isnan(psiTE(li,mi,1)), continue; end
            Hph_TE(:,:,kk,f) = (1j/sqrt(eta0))*psiTE(li,mi,f)*sphericalRicatti_h_prime(lval(li),2,kr(f))/j_l_1(li) ...
                                *(-gradYlm_ph_dom(:,:,kk)/Rs(index_r))  ;
        end
    end
end
Hph_total = sum(Hph_TM + Hph_TE, 3, 'omitnan');
clearvars Hph_TM Hph_TE

fignum = plot_field_comparison_polar(fignum, ...
    theta/pi, Hph_total(:,index_phi,index_f), [], Hptr(:,index_phi,index_r,1,index_f), ...
    phi/pi,   Hph_total(index_th,:,index_f),  [], Hptr(index_th,:,index_r,1,index_f), ...
    eta0, '$\eta_0|H_\phi|$ (V/m)', {"(dominant-modes) $H_\phi$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Hph_total

%% --- H_theta ---
Hth_TM = zeros(numth,numph,n_dom,numf); Hth_TE = Hth_TM;
for f = 1:numf
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        if dom_type(kk) == 0   % TM
            if isnan(psiTM(li,mi,1)), continue; end
            Hth_TM(:,:,kk,f) = (1/sqrt(eta0))*psiTM(li,mi,f)*sphericalRicatti_h(lval(li),2,kr(f))/j_l_1(li) ...
                                *(+gradYlm_ph_dom(:,:,kk)/Rs(index_r))  ;
        else                    % TE
            if isnan(psiTE(li,mi,1)), continue; end
            Hth_TE(:,:,kk,f) = (1j/sqrt(eta0))*psiTE(li,mi,f)*sphericalRicatti_h_prime(lval(li),2,kr(f))/j_l_1(li) ...
                                *(-gradYlm_th_dom(:,:,kk)/Rs(index_r))  ;
        end
    end
end
Hth_total = sum(Hth_TM + Hth_TE, 3, 'omitnan');
clearvars Hth_TM Hth_TE

fignum = plot_field_comparison_polar(fignum, ...
    theta/pi, Hth_total(:,index_phi,index_f), [], Hptr(:,index_phi,index_r,2,index_f), ...
    phi/pi,   Hth_total(index_th,:,index_f),  [], Hptr(index_th,:,index_r,2,index_f), ...
    eta0, '$\eta_0|H_\theta|$ (V/m)', {"(dominant-modes) $H_\theta$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Hth_total

%% --- E_phi ---
Eph_TM = zeros(numth,numph,n_dom,numf); Eph_TE = Eph_TM;
for f = 1:numf
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        if dom_type(kk) == 1   % TE
            if isnan(psiTE(li,mi,1)), continue; end
            Eph_TE(:,:,kk,f) = sqrt(eta0)*psiTE(li,mi,f)*sphericalRicatti_h(lval(li),2,kr(f))/j_l_1(li) ...
                                *(+gradYlm_th_dom(:,:,kk)/Rs(index_r))  ;
        else                    % TM
            if isnan(psiTM(li,mi,1)), continue; end
            Eph_TM(:,:,kk,f) = 1j*sqrt(eta0)*psiTM(li,mi,f)*sphericalRicatti_h_prime(lval(li),2,kr(f))/j_l_1(li) ...
                                *(-gradYlm_ph_dom(:,:,kk)/Rs(index_r))  ;
        end
    end
end
Eph_total = sum(Eph_TM + Eph_TE, 3, 'omitnan');
clearvars Eph_TM Eph_TE

fignum = plot_field_comparison_polar(fignum, ...
    theta/pi, Eph_total(:,index_phi,index_f), [], Eptr(:,index_phi,index_r,1,index_f), ...
    phi/pi,   Eph_total(index_th,:,index_f),  [], Eptr(index_th,:,index_r,1,index_f), ...
    1, '$|E_\phi|$ (V/m)', {"(dominant-modes) $E_\phi$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Eph_total

%% --- E_theta ---
Eth_TM = zeros(numth,numph,n_dom,numf); Eth_TE = Eth_TM;
for f = 1:numf
    for kk = 1:n_dom
        li = dom_lidx(kk); mi = dom_midx(kk);
        if dom_type(kk) == 1   % TE
            if isnan(psiTE(li,mi,1)), continue; end
            Eth_TE(:,:,kk,f) = sqrt(eta0)*psiTE(li,mi,f)*sphericalRicatti_h(lval(li),2,kr(f))/j_l_1(li) ...
                                *(-gradYlm_ph_dom(:,:,kk)/Rs(index_r))  ;
        else                    % TM
            if isnan(psiTM(li,mi,1)), continue; end
            Eth_TM(:,:,kk,f) = 1j*sqrt(eta0)*psiTM(li,mi,f)*sphericalRicatti_h_prime(lval(li),2,kr(f))/j_l_1(li) ...
                                *(-gradYlm_th_dom(:,:,kk)/Rs(index_r))  ;
        end
    end
end
Eth_total = sum(Eth_TM + Eth_TE, 3, 'omitnan');
clearvars Eth_TM Eth_TE

fignum = plot_field_comparison_polar(fignum, ...
    theta/pi, Eth_total(:,index_phi,index_f), [], Eptr(:,index_phi,index_r,2,index_f), ...
    phi/pi,   Eth_total(index_th,:,index_f),  [], Eptr(index_th,:,index_r,2,index_f), ...
    1, '$|E_\theta|$ (V/m)', {"(dominant-modes) $E_\theta$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Eth_total

clearvars Hptr Eptr gradYlm_th_dom gradYlm_ph_dom

%% --- Accepted Power Coupling T_lm ---
fignum = fignum+1; figure(fignum);
ax_S = gca;
set(ax_S, 'LineStyleOrder', {'-','--',':','-.'}')
if use_Q_method
    yyaxis left
end
set(ax_S,'YColor',get(groot,'defaultAxesXColor'))
hold on
plot(Freq, db(reshape(S_psi(1,1,:),1,numf)), "LineStyle","--", "LineWidth",LW, ...
    "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{00}$")
for snp_k = 1:n_dom
    kk = snp_order(snp_k);
    c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
    plot(Freq, db(T_psi(kk,:)), "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
end
hold off
if use_Q_method
    ylabel("$|T_{lm}|$ dB ($\psi$)","FontSize",FS)
else
    ylabel("$|T_{lm}|$ dB","FontSize",FS)
end

% Right axis: Q-derived T_lm (only when Q was read; otherwise T_Q is zero
% and db() would give -Inf).
if use_Q_method
    yyaxis right
    set(ax_S,'YColor',get(groot,'defaultAxesXColor'))
    hold on
    plot(Freq, db(reshape(S_Q(1,1,:),1,numf)), "LineStyle","--", "LineWidth",LW+0.5, ...
        "Color",get(groot,'defaultAxesXColor'), "HandleVisibility","off")
    for snp_k = 1:n_dom
        kk = snp_order(snp_k);
        c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
        plot(Freq, db(T_Q(kk,:)), "LineStyle",":", "LineWidth",LW+0.5, "Color",c, "HandleVisibility","off")
    end
    hold off
    ylabel("$|T_{lm}|$ dB ($Q$)","FontSize",FS)
end

if use_Q_method
    yyaxis left;  ylim([round(min(max(db(T_psi(:,fband_idx)), [], 2)))-5 0])
    yyaxis right; ylim([round(min(max(db(T_Q(:,fband_idx)),   [], 2)))-5 0])
else
    ylim([round(min(max(db(T_psi(:,fband_idx)), [], 2)))-5 0])
end

if use_Q_method
    title("Accepted Power Coupling $T_{lm} = S_{lm}/\sqrt{1-|S_{00}|^2}$: $\psi$ (solid) vs $Q$ (dotted)","FontSize",FS)
else
    title("Accepted Power Coupling $T_{lm} = S_{lm}/\sqrt{1-|S_{00}|^2}$","FontSize",FS)
end
legend("FontSize",FS)
xlabel("Frequency (Hz)","FontSize",FS)
grid on; xlim([fmin,fmax])

% Phase plot — same snp_order
fignum = fignum+1; figure(fignum);
ax_Sph = gca;
set(ax_Sph, 'LineStyleOrder', {'-','--',':','-.'}')
if use_Q_method
    yyaxis left
end
set(ax_Sph,'YColor',get(groot,'defaultAxesXColor'))
hold on
plot(Freq, angle(reshape(S_psi(1,1,:),1,numf))*180/pi, "LineStyle","--", "LineWidth",LW, ...
    "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{00}$")
for snp_k = 1:n_dom
    kk = snp_order(snp_k);
    c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
    plot(Freq, angle(T_psi(kk,:))*180/pi, "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
end
hold off
if use_Q_method
    ylabel("$\angle T_{lm}$ deg ($\psi$)","FontSize",FS)
else
    ylabel("$\angle T_{lm}$ deg","FontSize",FS)
end

if use_Q_method
    yyaxis right
    set(ax_Sph,'YColor',get(groot,'defaultAxesXColor'))
    hold on
    plot(Freq, angle(reshape(S_Q(1,1,:),1,numf))*180/pi, "LineStyle","--", "LineWidth",LW+0.5, ...
        "Color",get(groot,'defaultAxesXColor'), "HandleVisibility","off")
    for snp_k = 1:n_dom
        kk = snp_order(snp_k);
        c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
        plot(Freq, angle(T_Q(kk,:))*180/pi, "LineStyle",":", "LineWidth",LW+0.5, "Color",c, "HandleVisibility","off")
    end
    hold off
    ylabel("$\angle T_{lm}$ deg ($Q$)","FontSize",FS)
end

if use_Q_method
    yyaxis left;  ylim([-180 180])
    yyaxis right; ylim([-180 180])
    title("Accepted Power Coupling Phase $\angle T_{lm}$: $\psi$ (solid) vs $Q$ (dotted)","FontSize",FS)
else
    ylim([-180 180])
    title("Accepted Power Coupling Phase $\angle T_{lm}$","FontSize",FS)
end
legend("FontSize",FS)
xlabel("Frequency (Hz)","FontSize",FS)
grid on; xlim([fmin,fmax])

%% --- Radiated Power Comparison ---
Pr_PSI = zeros(numf,1); Pr_PSIdom = zeros(numf,1);

mask_TM_dom = false(L, M);
mask_TE_dom = false(L, M);
for kk = 1:n_dom_TM
    mask_TM_dom(dom_TM_lidx(kk), dom_TM_midx(kk)) = true;
end
for kk = 1:n_dom_TE
    mask_TE_dom(dom_TE_lidx(kk), dom_TE_midx(kk)) = true;
end

for f = 1:numf
    Pr_PSI(f) = (1/2) * sum(abs(psiTM(:,:,f)).^2 + abs(psiTE(:,:,f)).^2, 'all', 'omitnan');
    psiTM_f = psiTM(:,:,f); psiTE_f = psiTE(:,:,f);
    Pr_PSIdom(f) = (1/2) * (sum(abs(psiTM_f(mask_TM_dom)).^2, 'omitnan') + sum(abs(psiTE_f(mask_TE_dom)).^2, 'omitnan'));
end

if use_Q_method
    Pr_Q = (1/2) * 8*pi * sum(abs(Q).^2, 1).';

    Pr_Qdom = zeros(numf,1);
    for f = 1:numf
        Pr_Qdom_f = 0;
        for kk = 1:n_dom
            li = dom_lidx(kk); mi = dom_midx(kk);
            if dom_type(kk) == 0
                j_pos = SphModeIndex(lval(li),  mval(mi), 'TM');
                Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(Q(j_pos,f))^2;
                if mval(mi) ~= 0
                    j_neg = SphModeIndex(lval(li), -mval(mi), 'TM');
                    Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(Q(j_neg,f))^2;
                end
            else
                j_pos = SphModeIndex(lval(li),  mval(mi), 'TE');
                Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(Q(j_pos,f))^2;
                if mval(mi) ~= 0
                    j_neg = SphModeIndex(lval(li), -mval(mi), 'TE');
                    Pr_Qdom_f = Pr_Qdom_f + (1/2)*8*pi*abs(Q(j_neg,f))^2;
                end
            end
        end
        Pr_Qdom(f) = Pr_Qdom_f;
    end
else
    % Q unavailable: NaN leaves a gap in the plot.
    Pr_Q    = nan(numf,1);
    Pr_Qdom = nan(numf,1);
end

Pr_Spsi = zeros(numf,1);
for f = 1:numf
    port_voltage_factor = abs(V0(1) / (1 + S_psi(1,1,f)))^2;
    Pr_Spsi(f) = (1/2) * Y0 * port_voltage_factor * 2^abs(even_odd_excitation) * ...
                 sum(abs(S_psi(2:end,1,f)).^2);
end
if use_Q_method
    Pr_SQ = zeros(numf,1);
    for f = 1:numf
        port_voltage_factor = abs(V0(1) / (1 + S_Q(1,1,f)))^2;
        Pr_SQ(f) = (1/2) * Y0 * port_voltage_factor * 2^abs(even_odd_excitation) * ...
                   sum(abs(S_Q(2:end,1,f)).^2);
    end
else
    % Q unavailable: NaN leaves a gap in the plot.
    Pr_SQ = nan(numf,1);
end

fignum = fignum+1; figure_P = figure(fignum);
hold on
plot(Freq, real(Pr_NF),   "LineWidth",LW)
plot(Freq, Pr_PSI,        "--", "LineWidth",LW+1)
plot(Freq, Pr_PSIdom,     "s",  "LineWidth",LW+1, "MarkerSize",6)
plot_labels = {"NF Poynting","$|\psi|^2$ all","$|\psi|^2$ dom"};
if use_Q_method
    plot(Freq, Pr_Q,    "--", "LineWidth",LW)
    plot(Freq, Pr_Qdom, "^",  "LineWidth",LW,   "MarkerSize",3)
    plot_labels = [plot_labels, {"$|Q|^2$ all","$|Q|^2$ dom"}];
end
plot(Freq, Pr_Spsi,       "o",  "LineWidth",LW,   "MarkerSize",4)
plot_labels = [plot_labels, {"$|S|^2$ $\psi$"}];
if use_Q_method
    plot(Freq, Pr_SQ, "d",  "LineWidth",LW,   "MarkerSize",4)
    plot_labels = [plot_labels, {"$|S|^2$ $Q$"}];
end
hold off
legend(plot_labels{:},"FontSize",FS)
title("Radiated Power: Different Methods","FontSize",FS)
ylabel("Watts","FontSize",FS); xlabel("Freq (Hz)","FontSize",FS)
grid on; xlim([fmin,fmax])

%% --- L2 Norm Error of Radiated Power vs NF Poynting ---
% Relative L2 error over [fmin, fmax] of each method against the near-field
% Poynting power Pr_NF (computed directly from the field data). Q methods
% are reported only when Q was read.
Pr_NF_band  = real(Pr_NF(fband_idx));
Pr_ref_norm = norm(Pr_NF_band);

L2_err_PSI    = norm(Pr_PSI(fband_idx)    - Pr_NF_band) / Pr_ref_norm;
L2_err_PSIdom = norm(Pr_PSIdom(fband_idx) - Pr_NF_band) / Pr_ref_norm;
L2_err_Spsi   = norm(Pr_Spsi(fband_idx)   - Pr_NF_band) / Pr_ref_norm;

fprintf('\nRadiated power L2 error vs NF Poynting (%.3g-%.3g GHz band):\n', fmin/1e9, fmax/1e9);
fprintf('  |psi|^2 all : %6.2f%%\n', 100*L2_err_PSI);
fprintf('  |psi|^2 dom : %6.2f%%\n', 100*L2_err_PSIdom);
if use_Q_method
    L2_err_Q    = norm(Pr_Q(fband_idx)    - Pr_NF_band) / Pr_ref_norm;
    L2_err_Qdom = norm(Pr_Qdom(fband_idx) - Pr_NF_band) / Pr_ref_norm;
    fprintf('  |Q|^2 all   : %6.2f%%\n', 100*L2_err_Q);
    fprintf('  |Q|^2 dom   : %6.2f%%\n', 100*L2_err_Qdom);
end
fprintf('  |S|^2 psi   : %6.2f%%\n', 100*L2_err_Spsi);
if use_Q_method
    L2_err_SQ = norm(Pr_SQ(fband_idx) - Pr_NF_band) / Pr_ref_norm;
    fprintf('  |S|^2 Q     : %6.2f%%\n', 100*L2_err_SQ);
end

%% --- Export Spherical Wave S-parameters to Touchstone File ---
% Ports 2..n_dom+1 follow snp_order (see snp_order_method).
if ~exist(snp_dir, 'dir'), mkdir(snp_dir); end

if snp_use_Q
    S_snp_base = S_Q;
    snp_label  = 'Q';
else
    S_snp_base = S_psi;
    snp_label  = 'psi';
end

snp_port_order = [1, snp_order(:).' + 1];
S_snp = S_snp_base(snp_port_order, snp_port_order, :);

% DC (0 Hz) point, per dc_resistance_method.
%
% Reactive-limit approximation: with no DC-conducting path the radiation
% resistance vanishes faster than the reactance changes, so Zin tends to a
% near-short (inductive) or near-open (capacitive). S11 is approximated to
% first order about that limit, with Rmin = real(Zant) at the lowest
% frequency keeping |S11_DC| just below 1.
Rmin = real(Zant(1));
if imag(Zant(1)) > 0
    dc_reactance_type = 'inductive';
    S11_DC_reactive = -(1 - 2*Rmin/Z0);
else
    dc_reactance_type = 'capacitive';
    S11_DC_reactive = +(1 - 2*Rmin/Z0);
end

% Internal-load DC-conduction model: internal loads (loaded, never excited
% ports; see ReadFEKOout) provide a DC path, so the DC impedance is the
% parallel combination of their resistances, Z1, and S11_DC =
% (Z1-Z0)/(Z1+Z0). Loads with zero real part are excluded (not treated as
% shorts); if none remain, 'auto' falls back to the reactive limit.
if ~isempty(internal_loads)
    Z_real_all = real([internal_loads.Z]);
    valid_mask = Z_real_all ~= 0;
    Z_loads = Z_real_all(valid_mask);
    used_load_names = {internal_loads(valid_mask).name};
    if any(~valid_mask)
        fprintf('Spherical modes: DC point — ignoring purely reactive internal load(s) with zero real part: %s.\n', ...
            strjoin({internal_loads(~valid_mask).name}, ', '));
    end
else
    Z_loads = [];
    used_load_names = {};
end
have_load_model = ~isempty(Z_loads);
if have_load_model
    Z1 = 1 / sum(1 ./ Z_loads);
    S11_DC_load = (Z1 - Z0) / (Z1 + Z0);
end

switch dc_resistance_method
    case 'user'
        S11_DC = input('Enter S11_DC (DC reflection coefficient) directly: ');
        fprintf('Spherical modes: DC point — user-specified value, S11_DC = %.4f\n', S11_DC);
    case 'Rmin'
        S11_DC = S11_DC_reactive;
        fprintf('Spherical modes: DC point — forced reactive-limit ("low-frequency") model: %s (Im[Zant]=%.3g Ohm), Rmin = %.4g Ohm, S11_DC = %.4f\n', ...
            dc_reactance_type, imag(Zant(1)), Rmin, S11_DC);
    otherwise   % 'auto'
        if have_load_model
            S11_DC = S11_DC_load;
            fprintf('Spherical modes: DC point — DC-conduction model via internal load(s) (%s), Z1 = %.4g Ohm, S11_DC = %.4f\n', ...
                strjoin(used_load_names, ', '), Z1, S11_DC);
        else
            S11_DC = S11_DC_reactive;
            fprintf('Spherical modes: DC point — no usable internal load found; %s (Im[Zant]=%.3g Ohm), Rmin = %.4g Ohm, S11_DC = %.4f\n', ...
                dc_reactance_type, imag(Zant(1)), Rmin, S11_DC);
        end
end

S_DC      = zeros(n_dom+1, n_dom+1);
S_DC(1,1) = S11_DC;

% Export only frequencies up to fmax.
export_freq_idx = find(Freq <= fmax);
Freq_export = Freq(export_freq_idx);
S_snp_export = S_snp(:,:,export_freq_idx);

snp_path = fullfile(snp_dir, "S_dom_"+snp_label+"_"+antenna_name+".s"+num2str(n_dom+1)+"p");
rfwrite(cat(3, S_DC, S_snp_export), [0, Freq_export], snp_path, "ReferenceResistance", Z0)

% Build ADS-compatible port names (no spaces/special chars): port 1 is the
% antenna feed; ports 2..n_dom+1 are named by mode type/l/m, e.g. TM_L1_M0, TE_L2_Mneg1.
port_names = strings(1, n_dom+1);
port_names(1) = "Feed";
for snp_k = 1:n_dom
    kk = snp_order(snp_k);
    li = dom_lidx(kk); mi = dom_midx(kk);
    type_str = 'TM'; if dom_type(kk) == 1, type_str = 'TE'; end
    if mval(mi) < 0
        m_str = sprintf('Mneg%d', abs(mval(mi)));
    else
        m_str = sprintf('M%d', mval(mi));
    end
    port_names(snp_k+1) = sprintf('%s_L%d_%s', type_str, lval(li), m_str);
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

% Port names as "! Port[n] name" comments, the form ADS writes in its own
% Touchstone V1.1 files. ADS shows them only when they are the only comment
% lines before the '#' option line besides the header lines above; a
% [Port Names] keyword is not part of the spec and ADS rejects it.
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

