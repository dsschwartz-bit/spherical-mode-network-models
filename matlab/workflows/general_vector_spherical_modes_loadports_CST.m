%GENERAL_VECTOR_SPHERICAL_MODES_LOADPORTS_CST  Multi-port CST antenna -> spherical-mode network model.
%
% CST counterpart of general_vector_spherical_modes_loadports_FEKO: port 1
% is the feed, ports 2..N are internal load ports. Spherical-mode
% coefficients psiTM/psiTE are computed per port with the far-field method
% only (no near field, no Q), modes are ranked by accepted-power coupling
% |T_lm| at the ranking port, and an (N+n_dom)-port Touchstone file is
% exported (N antenna ports followed by the n_dom dominant modes).
%
% Inputs (CST exports):
%   <base_directory>\<antenna_name>\SP\<sparam_file>
%       N-port Touchstone file. The antenna apothem a is read from its
%       comments ("! Parameters = {...; a=<value>}"), or prompted for.
%   <base_directory>\<antenna_name>\Export\Farfield\farfield (f=<freq>) [<port>].txt
%       one far-field text file per (port, frequency); found automatically
%       by discover_cst_farfield_files.
%
% Coupling formula: the port excitation is the incident wave
% V0plus = sqrt(2*Z0*Pin_CST) (not a Thevenin source as in FEKO), so Z0
% cancels and S_{mode,p} = psi_p / sqrt(2^(|even_odd_excitation|+1)*Pin_CST)
% * exp(-jka) / j^(l+1). This assumes CST's stimulated-power port
% convention matches the HFSS one. CST symmetry planes need no Pin scaling;
% symmetry is handled via even_odd_excitation.
%
% Outputs:
%   <snp_dir>/S_dom_psi_<antenna_name>_CST.s<N+n_dom>p with ADS-compatible
%   port names; per-port figures of impedance, far-field reconstruction,
%   T_lm and radiated power.

clear
close all
clc

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
% CST project layout: <base_directory>\<antenna_name>\SP\ holds the
% Touchstone file, <base_directory>\<antenna_name>\Export\Farfield\ the
% far-field text files.
base_directory = "C:\path\to\CST\";
antenna_name   = "circularloop_FDTD";
sparam_directory   = fullfile(base_directory, antenna_name, "SP");
farfield_directory = fullfile(base_directory, antenna_name, "Export", "Farfield");
sparam_file   = antenna_name+".s1p";

% Touchstone output folder.
snp_dir = fullfile(repo_root, 'output', 'SnP');
% =========================

%% --- Discover Far-Field Files ---
% cst_farfield_files{p}: far-field file names for port p, one per entry of
% cst_farfield_freqs (in cst_freq_unit). All ports share one frequency
% list; S-parameters are interpolated onto it (see "Read S-Parameters").
[cst_farfield_files, cst_farfield_freqs] = discover_cst_farfield_files(farfield_directory);

%% --- Unit Configuration ---
% Units of plain numbers in the CST exports (they carry no unit labels).
%   cst_length_unit  - apothem in the Touchstone "! Parameters = {...}"
%                      comment (an HFSS-style "! Variables:" block states
%                      its own unit, metres if none)
%   cst_freq_unit    - frequencies in the far-field file names (the
%                      Touchstone sweep uses the unit on its '#' line)
%   cst_voltage_unit - far-field mag(rEPhi)/mag(rETheta) columns (rE is a
%                      voltage, despite the "[V/m]" headers; see
%                      read_cst_farfield)
%   cst_field_unit   - far-field Abs(E) column (not used)
% Recognized (case-insensitive): length 'm','cm','mm','um','nm';
% frequency 'Hz','kHz','MHz','GHz'; voltage 'V','mV','kV'.
cst_length_unit  = 'mm';
cst_freq_unit    = 'GHz';
cst_voltage_unit = 'V';
cst_field_unit   = 'V/m';   % not used

%% --- Physical Constants ---
c0  = physconst('lightspeed');
mu0 = 4*pi*1e-7;          % Permeability of free space (H/m)
eps0 = 1/c0^2/mu0;        % Permittivity of free space (F/m)
eta0 = sqrt(mu0/eps0);    % Characteristic impedance of free space (Ohm)

%% --- Simulation Parameters ---
f0     = 1e9;         % Reference frequency for the reconstruction plot (Hz)
fmin   = 0.05e9;       % Minimum frequency for plot x-axis (Hz)
fmax   = 4e9;         % Maximum frequency for plot x-axis (Hz)
lmax   = 9;           % Maximum spherical mode order to compute
n_dom  = 7;           % Number of dominant modes to identify and plot

% antenna_apothem is read from the Touchstone comments if present (see
% local_read_touchstone); otherwise you are prompted for it.
antenna_apothem   = NaN;
even_odd_excitation = 0;    % symmetry-plane mode-counting exponent (2^|even_odd_excitation| excited ports per group); 0 = no symmetry reduction -- EDIT if your model uses symmetry planes

Pin_CST = 0.5;   % Watt -- CST's discrete/waveguide port stimulated power
                 % (0.5 W by definition); edit if yours differs

% true: the far-field phase j^(l+1) is absorbed into psiTM/psiTE (via
% j_l_1), so the coupling formula uses jl1 = 1. false: the coupling formula
% applies jl1 = (1j)^(l+1).
use_jl1_convention = true;

% true: apply exp(-j*k*a) to move the S-parameter reference plane from each
% port's feed to r = antenna_apothem. false: reference stays at the feeds.
apply_ka_deembedding = true;

% SnP port ordering: 'peak_freq' (ascending freq of peak |T_lm|) or
% 'mode_index' (FEKO mode index j_ = 2*(l*(l+1)+m-1)+s, s=1 TE/2 TM)
snp_order_method = 'mode_index';

% DC (0 Hz) point of the SnP export:
%   'parallel_loads' / 'series_loads' - DC-conduction model: each port's
%       load resistance combined in parallel / series. A port's resistance
%       defaults to Z0 (a port terminated in its reference impedance);
%       override it per port with known_load_resistance.
%   'Rmin' - each port's lowest-frequency input resistance, using the
%       reactive-limit approximation (radiation resistance vanishes as f->0)
%   'user' - prompt for the DC values
dc_resistance_method = 'parallel_loads';
known_load_resistance = nan(1, numel(cst_farfield_files));   % (Ohm) only used for 'parallel_loads'/'series_loads'; NaN entries default to that port's own Z0 -- e.g. known_load_resistance(2) = 1e6 to override just port 2;

%% --- Plot Settings ---
LW = 2;    % LineWidth
FS = 14;   % FontSize

set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% --- Read S-Parameters from Touchstone (.sNp) File ---
% local_read_touchstone (local function) is a self-contained reader (no RF
% Toolbox needed). num_internal is the file's port count.
sparam_path = fullfile(sparam_directory, sparam_file);
[S_dense, Freq_dense, Z0, apothem_from_file] = local_read_touchstone(sparam_path, cst_length_unit);
num_internal = size(S_dense,1);
fprintf('general_vector_spherical_modes_loadports_CST: read %d-port Touchstone file "%s": %d frequency point(s), %.4g to %.4g GHz, Z0 = %.4g Ohm.\n', ...
    num_internal, sparam_path, numel(Freq_dense), min(Freq_dense)/1e9, max(Freq_dense)/1e9, Z0);

if numel(cst_farfield_files) ~= num_internal
    error('general_vector_spherical_modes_loadports_CST: cst_farfield_files has %d entries but the Touchstone file has %d port(s) -- these must match one-to-one.', ...
        numel(cst_farfield_files), num_internal);
end
for p_idx = 1:num_internal
    if numel(cst_farfield_files{p_idx}) ~= numel(cst_farfield_freqs)
        error('general_vector_spherical_modes_loadports_CST: cst_farfield_files{%d} has %d filename(s) but cst_farfield_freqs has %d entries -- these must match one-to-one.', ...
            p_idx, numel(cst_farfield_files{p_idx}), numel(cst_farfield_freqs));
    end
end

% Working frequency axis: only the frequencies where far fields were solved.
% Far-field patterns are not interpolated in frequency; the (smooth)
% S-parameters are linearly interpolated onto these points instead.
Freq = cst_farfield_freqs(:).' * cst_unit_scale(cst_freq_unit);
numf = numel(Freq);
if min(Freq) < min(Freq_dense) || max(Freq) > max(Freq_dense)
    error('general_vector_spherical_modes_loadports_CST: cst_farfield_freqs (%.4g to %.4g GHz) falls outside the Touchstone file''s own frequency range (%.4g to %.4g GHz) -- cannot extrapolate S-parameters.', ...
        min(Freq)/1e9, max(Freq)/1e9, min(Freq_dense)/1e9, max(Freq_dense)/1e9);
end
S = zeros(num_internal, num_internal, numf);
for ii = 1:num_internal
    for jj = 1:num_internal
        S(ii,jj,:) = interp1(Freq_dense, reshape(S_dense(ii,jj,:),1,[]), Freq, 'linear');
    end
end
fprintf('general_vector_spherical_modes_loadports_CST: S-parameters interpolated from %d Touchstone point(s) down to the %d far-field frequency point(s): %s GHz.\n', ...
    numel(Freq_dense), numf, strjoin(string(Freq/1e9), ', '));

% Antenna apothem: from the Touchstone comments if found, else prompt.
if ~isnan(apothem_from_file)
    antenna_apothem = apothem_from_file;
    fprintf('general_vector_spherical_modes_loadports_CST: antenna apothem (a) found in Touchstone file comments: %.6g m.\n', antenna_apothem);
else
    fprintf('general_vector_spherical_modes_loadports_CST: WARNING - antenna apothem (a) not found in "%s"''s comments (expected an HFSS "! Variables:" block or a CST "! Parameters = {...}" line with a=<value>).\n', sparam_path);
    antenna_apothem = [];
    while isempty(antenna_apothem) || ~isnumeric(antenna_apothem) || ~isscalar(antenna_apothem) || ~isreal(antenna_apothem) || antenna_apothem <= 0
        antenna_apothem = input('Enter antenna apothem (meters): ');
    end
end

if fmax > Freq(numf), fmax = Freq(numf); end
if fmin < Freq(1),    fmin = Freq(1);    end

k = 2*pi*Freq(:).' / c0;   % (1 x numf) wavenumber vs frequency


int_name = strings(1, num_internal);
int_name(1) = "Feed";
for k_sec = 2:num_internal
    int_name(k_sec) = "load" + num2str(k_sec-1);
end
fprintf('general_vector_spherical_modes_loadports_CST: %d internal port(s): %s\n', num_internal, strjoin(int_name, ', '));

%% --- Compute Spherical-Mode Coefficients for Every Internal Port ---
% compute_port_FF_coeffs (local function) reads each port's far-field files
% and computes psiTM/psiTE with the far-field method.
port_data = cell(1, num_internal);
for p_idx = 1:num_internal
    port_files = cst_farfield_files{p_idx};
    port_paths = strings(1, numel(port_files));
    for ff_idx = 1:numel(port_files)
        port_paths(ff_idx) = fullfile(farfield_directory, port_files(ff_idx));
    end
    port_data{p_idx} = compute_port_FF_coeffs(port_paths, ...
        Freq, lmax, use_jl1_convention, cst_voltage_unit, int_name(p_idx));
end

% Dominant-mode ranking uses the Feed's own coefficients when available,
% else the first other port that has usable data.
rank_port = 1;
if ~port_data{1}.available
    rank_port = find(cellfun(@(pd) pd.available, port_data), 1);
    if isempty(rank_port)
        error('general_vector_spherical_modes_loadports_CST: no internal port has usable far-field data -- check cst_farfield_files.')
    end
    fprintf('general_vector_spherical_modes_loadports_CST: WARNING - primary port "%s" has no usable far-field data; ranking dominant modes from port "%s" instead.\n', ...
        int_name(1), int_name(rank_port));
end
pd_rank = port_data{rank_port};

% Unpack the ranking port's data.
psiTM = pd_rank.psiTM; psiTE = pd_rank.psiTE;
L = pd_rank.L; M = pd_rank.M; lval = pd_rank.lval; mval = pd_rank.mval;
l_factor = pd_rank.l_factor; j_l_1 = pd_rank.j_l_1; CSphase = pd_rank.CSphase;
Norm = pd_rank.Norm; Normp = pd_rank.Normp; midx_0 = find(mval == 0);
theta = pd_rank.theta; phi = pd_rank.phi;
[~,index_f] = min(abs(Freq-f0));

%% --- Build S-parameters for EVERY (l,m,type), for Ranking (ranking port only) ---
% Modes are ranked by their accepted-power coupling at the ranking port,
% T_lm = S_{mode,p}/sqrt(1-|S_pp|^2), which, unlike raw |psi|, accounts for
% how much power the port actually accepts at each frequency. The full
% multi-port coupling is built later for the dominant modes.
S_pp_rank = reshape(S(rank_port, rank_port, :), 1, numf);   % ranking port's own S_pp
if apply_ka_deembedding
    ka_factor = exp(-1j * k(:).' * antenna_apothem);        % (1 x numf)
else
    ka_factor = ones(1, numf);                              % no de-embedding phase
end
if use_jl1_convention
    jl1_l = ones(L,1);                                      % phase already absorbed into psiTM/psiTE via j_l_1
else
    jl1_l = (1j).^(lval(:) + 1);                            % (L x 1)
end
ka_factor_r = reshape(ka_factor, 1, 1, numf);
jl1_r    = reshape(jl1_l, L, 1, 1);                         % broadcasts over M and numf

psiTM_safe = psiTM; psiTM_safe(isnan(psiTM_safe)) = 0;
psiTE_safe = psiTE; psiTE_safe(isnan(psiTE_safe)) = 0;

S_TM_all = psiTM_safe ./ sqrt(2^(abs(even_odd_excitation)+1) * Pin_CST) .* ka_factor_r ./ jl1_r;   % (L x M x numf)
S_TE_all = psiTE_safe ./ sqrt(2^(abs(even_odd_excitation)+1) * Pin_CST) .* ka_factor_r ./ jl1_r;

T_denom_rank = local_T_denom(S_pp_rank);                    % (1 x numf); NaN where |S_pp|->1 (no accepted power)
T_TM_all = S_TM_all ./ reshape(T_denom_rank, 1, 1, numf);   % (L x M x numf)
T_TE_all = S_TE_all ./ reshape(T_denom_rank, 1, 1, numf);

%% --- Identify Dominant Modal Coefficients across TE and TM jointly ---
% TE and TM modes are pooled and the n_dom largest by max |T_lm| within
% [fmin, fmax] are selected, then ordered by peak frequency.
fband_idx = find(Freq >= fmin & Freq <= fmax);

[TM_mag, TM_fmax_band] = max(abs(reshape(T_TM_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
[TE_mag, TE_fmax_band] = max(abs(reshape(T_TE_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
TM_fmax  = fband_idx(TM_fmax_band(:)).';
TE_fmax  = fband_idx(TE_fmax_band(:)).';

TM_mag(isnan(TM_mag)) = 0;
TE_mag(isnan(TE_mag)) = 0;

pool_mag     = [TM_mag(:);  TE_mag(:)  ];
pool_fmax    = [TM_fmax(:); TE_fmax(:) ];
pool_type    = [zeros(L*M,1); ones(L*M,1)];
pool_lmidx   = [(1:L*M).';   (1:L*M).'  ];

[~, sort_pool] = sort(pool_mag, 'descend');
dom_pool_mag   = sort_pool(1:n_dom);
dom_peak_freqs = Freq(pool_fmax(dom_pool_mag));
[~, freq_order] = sort(dom_peak_freqs(:), 'ascend');
dom_pool = dom_pool_mag(freq_order);

dom_type = pool_type(dom_pool);      % 0=TM, 1=TE, length n_dom
dom_lm   = pool_lmidx(dom_pool);     % linear index into (L x M), length n_dom
[dom_lidx, dom_midx] = ind2sub([L, M], dom_lm);

dom_TM_lidx = dom_lidx(dom_type == 0);
dom_TM_midx = dom_midx(dom_type == 0);
dom_TE_lidx = dom_lidx(dom_type == 1);
dom_TE_midx = dom_midx(dom_type == 1);
n_dom_TM = length(dom_TM_lidx);
n_dom_TE = length(dom_TE_lidx);

fprintf('\n--- Top %d dominant modes (TE and TM ranked jointly by max |T_lm|, accepted-power coupling, port "%s") ---\n', n_dom, int_name(rank_port))
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0
        mode_vs_freq = abs(reshape(T_TM_all(li,mi,:), 1, numf));
        type_str = 'TM';
    else
        mode_vs_freq = abs(reshape(T_TE_all(li,mi,:), 1, numf));
        type_str = 'TE';
    end
    [~, f_max_idx_band] = max(mode_vs_freq(fband_idx));
    f_max_idx = fband_idx(f_max_idx_band);
    f_max = Freq(f_max_idx);
    all_at_f = [abs(reshape(T_TM_all(:,:,f_max_idx), [], 1)); abs(reshape(T_TE_all(:,:,f_max_idx), [], 1))];
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
% S_mode is (num_internal+n_dom) x (num_internal+n_dom) x numf.
% Indices 1..num_internal: internal ports (the port-to-port block is S).
% Indices num_internal+1..end: dominant spherical modes; each port's
% coupling to each mode uses that port's own psiTM/psiTE (zero if
% unavailable):
%   S_{mode,p} = psi_p / sqrt(2^(|even_odd_excitation|+1) * Pin_CST) * exp(-j*k*a) / j^(l+1)
% i.e. port_factor*psi/V0plus with port_factor = sqrt(Z0/2^q) and
% V0plus = sqrt(2*Z0*Pin_CST); Z0 cancels. Mode ports are not coupled to
% each other.
n_tot = num_internal + n_dom;
S_mode = zeros(n_tot, n_tot, numf);
S_mode(1:num_internal, 1:num_internal, :) = S;

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
            else                   % TE
                psi_val = NaN;
                if pd.available, psi_val = pd.psiTE(li,mi,f); end
            end
            if isnan(psi_val), psi_val = 0; end

            coupling = psi_val / sqrt(2^(abs(even_odd_excitation)+1) * Pin_CST) * exp_ka / jl1;

            S_mode(num_internal+kk, p, f) = coupling;
            S_mode(p, num_internal+kk, f) = coupling;
        end
    end
end
% Sanitise: replace any remaining NaN with 0 before export
S_mode(isnan(S_mode)) = 0;

%% --- Per-Port Plots: Impedance, Reconstruction, T_lm, Power ---
% Plots for every port with usable far-field data.
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

% snp_order: mode ordering for the SnP ports and T_lm legends, based on the
% Feed (port 1) so it is the same in every port's plots.
T_feed = zeros(n_dom, numf);
denom_feed = local_T_denom(reshape(S_mode(1,1,:), 1, numf));
for kk = 1:n_dom
    T_feed(kk,:) = reshape(S_mode(num_internal+kk,1,:), 1, numf) ./ denom_feed;
end
if strcmpi(snp_order_method, 'mode_index')
    dom_l = lval(dom_lidx);
    dom_m = mval(dom_midx);
    s_val = ones(n_dom,1); s_val(dom_type == 0) = 2;   % TM=2, TE=1
    dom_j = 2*(dom_l.*(dom_l+1) + dom_m(:) - 1) + s_val;
    [~, snp_order] = sort(dom_j, 'ascend');
else   % 'peak_freq' (default)
    [~, T_fmax_band] = max(abs(T_feed(:,fband_idx)), [], 2);
    T_fmax = fband_idx(T_fmax_band);
    [~, snp_order] = sort(T_fmax, 'ascend');
end

mask_TM_dom = false(L, M);
mask_TE_dom = false(L, M);
for kk = 1:n_dom_TM
    mask_TM_dom(dom_TM_lidx(kk), dom_TM_midx(kk)) = true;
end
for kk = 1:n_dom_TE
    mask_TE_dom(dom_TE_lidx(kk), dom_TE_midx(kk)) = true;
end

fignum = 0;
for p = 1:num_internal
    pd = port_data{p};

    if ~pd.available
        fprintf('general_vector_spherical_modes_loadports_CST: skipping plots for port "%s" -- no usable far-field data.\n', int_name(p));
        continue
    end

    % --- Input Impedance ---
    Zin_p = Z0 * (1 + reshape(S_mode(p,p,:), 1, numf)) ./ (1 - reshape(S_mode(p,p,:), 1, numf));
    fignum = fignum+1; figure(fignum)
    plot(Freq, real(Zin_p), "Color","b","LineWidth",LW)
    hold on
    plot(Freq, imag(Zin_p), "Color","r","LineWidth",LW)
    hold off
    title("Impedance ("+int_name(p)+")","FontSize",FS)
    legend("$R$","$X$","FontSize",FS)
    xlim([fmin,fmax])
    ylim([-Z0*4,Z0*4])
    grid on

    % --- Far-field reconstruction vs imported far field (polar cuts) ---
    % All-modes expansion (line) vs the imported far field (markers).
    [~,index_phi_p] = min(abs(pd.phi - 0));
    [~,index_th_p]  = min(abs(pd.theta - pi/2));
    [PHI_p, THETA_p] = meshgrid(pd.phi, pd.theta);

    [Eth_modes, Eph_modes] = reconstruct_farfield_from_modes( ...
        pd.psiTM(:,:,index_f), pd.psiTE(:,:,index_f), pd.theta, pd.phi, lval, mval, eta0);
    Eth_imp = pd.rEpt(:,:,2,index_f);   % imported: dim3 1 = rE_phi, 2 = rE_theta
    Eph_imp = pd.rEpt(:,:,1,index_f);

    ds_th = 1:numel(pd.theta);   % imported-field markers at every sample (no downsampling)
    ds_ph = 1:numel(pd.phi);

    field_title_base = num2str(Freq(index_f)/1e9,2)+"GHz, \#dom="+num2str(n_dom)+", "+int_name(p);
    phi_label_p   = "$\phi="+num2str(pd.phi(index_phi_p)*180/pi,3)+"^\circ$";
    theta_label_p = "$\theta="+num2str(pd.theta(index_th_p)*180/pi,3)+"^\circ$";

    fignum = fignum+1; figure(fignum)
    local_ff_polar(pd.theta, abs(Eth_modes(:,index_phi_p)), abs(Eth_imp(:,index_phi_p)), ds_th, ...
        {"(all-modes) $E_\theta$ far field vs $\theta$", field_title_base, phi_label_p}, "$|rE_\theta|$ (V)", LW, FS)

    fignum = fignum+1; figure(fignum)
    local_ff_polar(pd.phi, abs(Eth_modes(index_th_p,:)).', abs(Eth_imp(index_th_p,:)).', ds_ph, ...
        {"(all-modes) $E_\theta$ far field vs $\phi$", field_title_base, theta_label_p}, "$|rE_\theta|$ (V)", LW, FS)

    fignum = fignum+1; figure(fignum)
    local_ff_polar(pd.theta, abs(Eph_modes(:,index_phi_p)), abs(Eph_imp(:,index_phi_p)), ds_th, ...
        {"(all-modes) $E_\phi$ far field vs $\theta$", field_title_base, phi_label_p}, "$|rE_\phi|$ (V)", LW, FS)

    fignum = fignum+1; figure(fignum)
    local_ff_polar(pd.phi, abs(Eph_modes(index_th_p,:)).', abs(Eph_imp(index_th_p,:)).', ds_ph, ...
        {"(all-modes) $E_\phi$ far field vs $\phi$", field_title_base, theta_label_p}, "$|rE_\phi|$ (V)", LW, FS)

    clearvars Eth_modes Eph_modes Eth_imp Eph_imp

    % --- T_lm Magnitude / Phase (this port's own data) ---
    T_p = zeros(n_dom, numf);
    denom_p = local_T_denom(reshape(S_mode(p,p,:), 1, numf));
    for kk = 1:n_dom
        T_p(kk,:) = reshape(S_mode(num_internal+kk,p,:), 1, numf) ./ denom_p;
    end

    fignum = fignum+1; figure(fignum);
    hold on
    plot(Freq, db(reshape(S_mode(p,p,:),1,numf)), "LineStyle","--", "LineWidth",LW, ...
        "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{"+int_name(p)+int_name(p)+"}$")
    for snp_k = 1:n_dom
        kk = snp_order(snp_k);
        c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
        plot(Freq, db(T_p(kk,:)), "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
    end
    hold off
    ylabel("$|T_{lm}|$ dB","FontSize",FS)
    T_peak_db = max(db(T_p(:,fband_idx)), [], 2, 'omitnan');   % each mode's peak (NaN points ignored)
    T_ylo = round(min(T_peak_db(isfinite(T_peak_db)))) - 5;
    if isempty(T_ylo) || ~isfinite(T_ylo) || T_ylo >= 0
        T_ylo = -60;   % fallback if no finite peak below 0 dB
    end
    ylim([T_ylo 0])
    title("Accepted Power Coupling $T_{lm} = S_{lm}/\sqrt{1-|S_{"+int_name(p)+int_name(p)+"}|^2}$ ("+int_name(p)+"), Far-Field method","FontSize",FS)
    legend("FontSize",FS)
    xlabel("Frequency (Hz)","FontSize",FS)
    grid on; xlim([fmin,fmax])

    fignum = fignum+1; figure(fignum);
    hold on
    plot(Freq, angle(reshape(S_mode(p,p,:),1,numf))*180/pi, "LineStyle","--", "LineWidth",LW, ...
        "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{"+int_name(p)+int_name(p)+"}$")
    for snp_k = 1:n_dom
        kk = snp_order(snp_k);
        c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
        plot(Freq, angle(T_p(kk,:))*180/pi, "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
    end
    hold off
    ylabel("$\angle T_{lm}$ deg","FontSize",FS)
    ylim([-180 180])
    title("Accepted Power Coupling Phase $\angle T_{lm}$ ("+int_name(p)+"), Far-Field method","FontSize",FS)
    legend("FontSize",FS)
    xlabel("Frequency (Hz)","FontSize",FS)
    grid on; xlim([fmin,fmax])

    % --- Radiated Power: direct far-field integral, |psi|^2 all/dom, |S|^2 ---
    Pr_FF_direct_p = zeros(numf,1);
    for f = 1:numf
        U = (abs(pd.rEpt(:,:,1,f)).^2 + abs(pd.rEpt(:,:,2,f)).^2) / (2*eta0);
        Pr_FF_direct_p(f) = trapz(pd.phi, trapz(pd.theta, U.*sin(THETA_p), 1), 2);
    end

    Pr_psi_p = zeros(numf,1);
    for f = 1:numf
        Pr_psi_p(f) = (1/2) * sum(abs(pd.psiTM(:,:,f)).^2 + abs(pd.psiTE(:,:,f)).^2, 'all', 'omitnan');
    end

    Pr_psidom_p = zeros(numf,1);
    for f = 1:numf
        psiTM_f = pd.psiTM(:,:,f);
        psiTE_f = pd.psiTE(:,:,f);
        Pr_psidom_p(f) = (1/2) * ( sum(abs(psiTM_f(mask_TM_dom)).^2) + sum(abs(psiTE_f(mask_TE_dom)).^2) );
    end

    % Pr_S_p: power from S_mode; equals Pr_psidom_p, since 2^q*Pin_CST
    % (q=|even_odd_excitation|) undoes S_mode's scaling.
    Pr_S_p = zeros(numf,1);
    for f = 1:numf
        Pr_S_p(f) = 2^abs(even_odd_excitation) * Pin_CST * sum(abs(S_mode(num_internal+1:end,p,f)).^2);
    end

    fignum = fignum+1; figure(fignum);
    hold on
    plot(Freq, Pr_FF_direct_p, "-",  "LineWidth",LW+1)
    plot(Freq, Pr_psi_p,       "--", "LineWidth",LW+1)
    plot(Freq, Pr_psidom_p,    "s",  "LineWidth",LW+1, "MarkerSize",6)
    plot(Freq, Pr_S_p,         "o",  "LineWidth",LW,   "MarkerSize",4)
    hold off
    legend("Far-field (direct)","$|\psi|^2$ all","$|\psi|^2$ dom","$|S|^2$","FontSize",FS)
    title("Radiated Power ("+int_name(p)+"): Different Methods","FontSize",FS)
    ylabel("Watts","FontSize",FS); xlabel("Freq (Hz)","FontSize",FS)
    grid on; xlim([fmin,fmax])

    Pr_ref_band_p = Pr_FF_direct_p(fband_idx);
    Pr_ref_norm_p = norm(Pr_ref_band_p);
    if Pr_ref_norm_p > 0
        L2_err_psi_p    = norm(Pr_psi_p(fband_idx)    - Pr_ref_band_p) / Pr_ref_norm_p;
        L2_err_psidom_p = norm(Pr_psidom_p(fband_idx) - Pr_ref_band_p) / Pr_ref_norm_p;
        L2_err_S_p      = norm(Pr_S_p(fband_idx)      - Pr_ref_band_p) / Pr_ref_norm_p;
        fprintf('Radiated power L2 error vs direct far-field integral (%s, over [%.3g,%.3g] GHz):\n', int_name(p), fmin/1e9, fmax/1e9);
        fprintf('  |psi|^2 all: %6.2f%%\n  |psi|^2 dom: %6.2f%%\n  |S|^2:       %6.2f%%\n', 100*L2_err_psi_p, 100*L2_err_psidom_p, 100*L2_err_S_p);
    end
end

%% --- Export Spherical Wave S-parameters to Touchstone File ---
% Ports 1..num_internal are internal ports; num_internal+1..n_tot are the
% dominant modes, ordered by snp_order.
if ~exist(snp_dir, 'dir'), mkdir(snp_dir); end

valid_dc_methods = {'Rmin', 'user', 'parallel_loads', 'series_loads'};
if ~ismember(dc_resistance_method, valid_dc_methods)
    error('general_vector_spherical_modes_loadports_CST: dc_resistance_method = ''%s'' is not recognized. Valid values are: %s.', ...
        dc_resistance_method, strjoin(valid_dc_methods, ', '));
end

snp_port_order = [1:num_internal, snp_order(:).' + num_internal];
S_snp = S_mode(snp_port_order, snp_port_order, :);

% DC (0 Hz): diagonal = each port's DC reflection, off-diagonal = DC
% coupling between ports. 'Rmin' treats each port as isolated at DC
% (first-order approximation of S=(Z-Z0)/(Z+Z0) about the short/open
% limit); 'parallel_loads'/'series_loads' combine the other ports' load
% resistances (Z0 unless set in known_load_resistance).
S_DC = zeros(n_tot, n_tot);

Z_load_all = known_load_resistance;
Z_load_all(isnan(Z_load_all)) = Z0;

for p = 1:num_internal
    Sp1 = S_mode(p, p, 1);
    Zin_p1 = Z0 * (1 + Sp1) / (1 - Sp1);
    switch dc_resistance_method
        case 'user'
            Sp_DC = input(sprintf('Enter S_DC for port "%s" directly: ', int_name(p)));
            fprintf('Spherical modes: DC point (%s) — user-specified value, S_DC = %.4f\n', int_name(p), Sp_DC);
        case {'parallel_loads', 'series_loads'}
            other_idx = setdiff(1:num_internal, p);
            if isempty(other_idx)
                % No other port to combine with: the port is isolated at
                % DC, so use the open/short limit (S=-1 inductive, S=+1
                % capacitive) instead of Z1 = Inf.
                if imag(Zin_p1) > 0
                    dc_reactance_type = 'inductive';
                    Sp_DC = -1;
                else
                    dc_reactance_type = 'capacitive';
                    Sp_DC = +1;
                end
                fprintf('Spherical modes: DC point (%s) — %s, no other internal port to combine with (isolated at DC), S_DC = %.4f\n', ...
                    int_name(p), dc_reactance_type, Sp_DC);
            else
                Z_others = Z_load_all(other_idx);
                if strcmp(dc_resistance_method, 'parallel_loads')
                    Z1 = 1 / sum(1 ./ Z_others);
                else
                    Z1 = sum(Z_others);
                end
                Sp_DC = (Z1 - Z0) / (Z1 + Z0);
                fprintf('Spherical modes: DC point (%s) — DC-conduction model (%s), Z1 = %.4g Ohm, S_DC = %.4f\n', ...
                    int_name(p), dc_resistance_method, Z1, Sp_DC);
            end
        otherwise   % 'Rmin'
            Rmin_p = real(Zin_p1);
            if imag(Zin_p1) > 0
                dc_reactance_type = 'inductive';
                Sp_DC = -(1 - 2*Rmin_p/Z0);
            else
                dc_reactance_type = 'capacitive';
                Sp_DC = +(1 - 2*Rmin_p/Z0);
            end
            fprintf('Spherical modes: DC point (%s) — %s (Im[Zin]=%.3g Ohm), Rmin = %.4g Ohm, S_DC = %.4f\n', ...
                int_name(p), dc_reactance_type, imag(Zin_p1), Rmin_p, Sp_DC);
    end
    S_DC(p,p) = Sp_DC;
end

% DC coupling between port pairs. Series/parallel models (all ports'
% resistances Ri, including the Feed's, or conductances Gi = 1/Ri):
%   series:   S_nm = 2*sqrt(Rn*Rm) / sum(all Ri)
%   parallel: S_nm = 2*sqrt(Gn*Gm) / sum(all Gi)
% 'user': entered manually. 'Rmin': |S_nm| at the lowest frequency, with the
% sign taken from its imaginary part so the DC value is real.
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
        for m = (n+1):num_internal
            Snm_lowfreq = S(n, m, 1);
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

snp_path = fullfile(snp_dir, "S_dom_psi_"+antenna_name+"_CST.s"+num2str(n_tot)+"p");
rfwrite(cat(3, S_DC, S_snp_export), [0, Freq_export], snp_path, "ReferenceResistance", Z0)

% Build ADS-compatible port names (no spaces/special chars).
port_names = strings(1, n_tot);
port_names(1:num_internal) = int_name;
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
comment_lines = comment_lines + "! Written by " + string(mfilename) + ".m (Far-Field method, CST)," +...
                " a="+num2str(antenna_apothem,3)+"m, Z0="+num2str(Z0,5)+"ohms, ka_factor="+ka_factor_str + newline;
comment_lines = comment_lines + "! Date : " + string(datetime('now')) + newline;

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

function pdata = compute_port_FF_coeffs(farfield_paths, Freq, lmax, use_jl1_convention, cst_voltage_unit, port_label)
%COMPUTE_PORT_FF_COEFFS  Far-field-method psiTM/psiTE for one CST port.
%   Reads the port's far-field text files (farfield_paths(f) at Freq(f)) via
%   read_cst_farfield, stacks them as (numth x numph x 2 x numf) and
%   projects onto the vector spherical harmonics.
%
%   Returns struct pdata:
%     .available     - true if every file was read, all share one
%                      (theta,phi) grid, and it covers the full sphere
%     .psiTM, .psiTE - (L x M x numf) mode coefficients
%     .theta, .phi   - angular grid (rad)
%     .rEpt          - (numth x numph x 2 x numf) complex far field
%                      (dim3: 1 = rE_phi, 2 = rE_theta)
%     .L, .M, .lval, .mval, .l_factor, .j_l_1, .CSphase, .Norm, .Normp
%                    - mode-index grids and normalisation
c0  = physconst('lightspeed');
mu0 = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0);

lval = (1:1:lmax).';
mval = (-lmax:lmax);
L = length(lval); M = length(mval);
numf = numel(Freq);

pdata = struct('available', false, ...
    'psiTM', nan(L,M,numf), 'psiTE', nan(L,M,numf), ...
    'theta', [], 'phi', [], 'rEpt', [], ...
    'L', L, 'M', M, 'lval', lval, 'mval', mval, ...
    'l_factor', [], 'j_l_1', [], 'CSphase', [], 'Norm', [], 'Normp', []);

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
pdata.l_factor = l_factor; pdata.j_l_1 = j_l_1; pdata.CSphase = CSphase;
pdata.Norm = Norm; pdata.Normp = Normp;

if numel(farfield_paths) ~= numf
    fprintf('compute_port_FF_coeffs (%s): %d far-field file(s) given but %d frequencies expected -- this port has no spherical-mode coupling.\n', port_label, numel(farfield_paths), numf);
    return
end

% CST reports the far field at r = 1 m including the propagation phase
% e^(-jk*1). The projection expects the pattern with e^(-jkr)/r removed, so
% each frequency is multiplied by e^(+jk). This is a pure phase and does
% not change radiated power.
k = 2*pi*Freq/c0;

theta = []; phi = []; rEpt_all = [];
for f = 1:numf
    try
        [rEpt_f, theta_f, phi_f] = read_cst_farfield(farfield_paths(f), cst_voltage_unit);
    catch ME
        fprintf('compute_port_FF_coeffs (%s): could not read far-field file "%s" (%s) -- this port has no spherical-mode coupling.\n', port_label, farfield_paths(f), ME.message);
        return
    end
    if f == 1
        theta = theta_f; phi = phi_f;
        numth = numel(theta); numph = numel(phi);
        rEpt_all = zeros(numth, numph, 2, numf);
    else
        if numel(theta_f) ~= numth || numel(phi_f) ~= numph || max(abs(theta_f(:)-theta(:))) > 1e-6 || max(abs(phi_f(:)-phi(:))) > 1e-6
            fprintf('compute_port_FF_coeffs (%s): far-field file "%s" uses a different (theta,phi) grid than "%s" -- this port has no spherical-mode coupling.\n', port_label, farfield_paths(f), farfield_paths(1));
            return
        end
    end
    rEpt_all(:,:,:,f) = rEpt_f * exp(1j*k(f)*1);
end

if ~local_full_angular_range(theta, phi)
    fprintf('compute_port_FF_coeffs (%s): far-field files for this port do not cover the full angular range (theta: 0 to pi, phi: 0 to 2*pi) -- this port has no spherical-mode coupling.\n', port_label);
    return
end

numth = numel(theta); numph = numel(phi);
[PHI, THETA] = meshgrid(phi, theta);

gradYlm_th = zeros(numth,numph,L,M); gradYlm_ph = gradYlm_th;
for lindx = 1:L
    for mindx = 1:M
        if abs(mval(mindx)) > lval(lindx)
            continue
        end
        [gradYlm_th(:,:,lindx,mindx),gradYlm_ph(:,:,lindx,mindx)] = gradYlm(lval(lindx),mval(mindx),THETA,PHI,have_CSphase) ;
    end
end

rHpt_all = zeros(size(rEpt_all));
rHpt_all(:,:,1,:) =  rEpt_all(:,:,2,:) / eta0;
rHpt_all(:,:,2,:) = -rEpt_all(:,:,1,:) / eta0;

psiTM = nan(L,M,numf); psiTE = psiTM;
for lindx = 1:L
    for mindx = 1:M
        if abs(mval(mindx)) > lval(lindx)
            continue
        end
        for f = 1:numf
            psiTM(lindx,mindx,f) = 1/sqrt(eta0) * ...
                trapz(phi,...
                        trapz( theta, ( rEpt_all(:,:,1,f).*(-gradYlm_ph(:,:,lindx,mindx)) + ...
                                        rEpt_all(:,:,2,f).*(-gradYlm_th(:,:,lindx,mindx)) ...
                                       ).*sin(THETA)...
                        ,1) ...
                ,2) ;

            psiTE(lindx,mindx,f) = sqrt(eta0) * ...
                trapz(phi,...
                        trapz( theta, ( rHpt_all(:,:,1,f).*(-gradYlm_ph(:,:,lindx,mindx)) + ...
                                        rHpt_all(:,:,2,f).*(-gradYlm_th(:,:,lindx,mindx)) ...
                                       ).*sin(THETA)...
                        ,1) ...
                ,2) ;
        end
    end
end

pdata.psiTM = psiTM; pdata.psiTE = psiTE;
pdata.theta = theta; pdata.phi = phi; pdata.rEpt = rEpt_all;
pdata.available = true;
end

function local_ff_polar(ang, E_modes, E_imp, ds_idx, ttl, rlabel_str, LW, FS)
% LOCAL_FF_POLAR  One polar cut: spherical-mode expansion (line) and
% imported far field (solid black triangle markers, every sample).
ang = ang(:); E_modes = E_modes(:); E_imp = E_imp(:);
polarplot(ang, E_modes, "LineWidth", LW, "DisplayName", "Mode expansion")
hold on
polarplot(ang(ds_idx), E_imp(ds_idx), "LineStyle", "none", "Marker", "^", "MarkerSize", 4, ...
    "MarkerFaceColor", "k", "MarkerEdgeColor", "k", "DisplayName", "Imported")
hold off
title(ttl, "FontSize", FS)
rl = get(gca, 'RAxis'); rl.Label.String = rlabel_str; rl.Label.Interpreter = 'latex';
legend("Location", "southoutside", "Orientation", "horizontal", "FontSize", FS)
end

function d = local_T_denom(S_pp)
% LOCAL_T_DENOM  sqrt(1-|S_pp|^2), the accepted-power normalization of the
% coupling T_lm = S_lm/sqrt(1-|S_pp|^2). Where the port is (numerically)
% fully reflective, |S_pp| >= 1 - tiny (e.g. a lossless antenna at its very
% lowest frequency, |S11| = 1 exactly), the accepted power is zero and T_lm
% is 0/0 or Inf: those points are returned as NaN so ranking, limits and
% plots simply skip them instead of being swamped by Inf.
d = sqrt(max(1 - abs(S_pp).^2, 0));
d(d < 1e-6) = NaN;
end

function tf = local_full_angular_range(theta_chk, phi_chk)
%LOCAL_FULL_ANGULAR_RANGE  True if theta spans [0,pi] and phi spans 2*pi.
%   phi is checked by its total span, not its start value, since the first
%   sample need not lie exactly at 0. Tolerance 1 degree.
tol = 1*pi/180;   % 1 degree tolerance
theta_chk = theta_chk(:);
phi_chk   = phi_chk(:);
tf = numel(theta_chk) > 1 && numel(phi_chk) > 1 && ...
     min(theta_chk) < tol && max(theta_chk) > (pi - tol) && ...
     (max(phi_chk) - min(phi_chk)) > (2*pi - tol);
end

function [S, Freq, Z0, apothem_from_file] = local_read_touchstone(filename, cst_length_unit)
%LOCAL_READ_TOUCHSTONE  Minimal N-port Touchstone (.sNp) reader.
%   Handles "# <Hz|kHz|MHz|GHz> S <MA|DB|RI> R <refZ>" and any N (taken from
%   the file extension). 2-port files use the S11,S21,S12,S22 order; other
%   N use row-major order.
%
%   Returns:
%     S    - (N x N x numf) complex, S(i,j,f) = response at port i, port j
%            driven
%     Freq - (numf x 1) frequencies, Hz
%     Z0   - reference impedance (Ohm) from the option line, or from a
%            "! Parameters = {Z0=...}" comment when the option line gives 0
%     apothem_from_file - antenna apothem (m), from an HFSS "! Variables:"
%            line "!  a = <value> [mm|cm|um|nm|m]" (metres if no unit) or a
%            CST "! Parameters = {...; a=<value>; ...}" line (in
%            cst_length_unit). The variable must be named "a". NaN if absent.
fid = fopen(filename, 'r');
if fid < 0
    error('local_read_touchstone: cannot open file: %s', filename);
end

freq_unit = 'GHz'; fmt = 'MA'; Z0 = 50;
Z0_from_comment = NaN;   % reference impedance from CST's "! Parameters = {Z0=...}" comment
apothem_from_file = NaN;
data_lines = {};
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end
    if line(1) == '!'
        atok = regexp(line, '^!\s*a\s*=\s*([\d.eE+-]+)\s*(mm|cm|um|nm|m)?\s*$', 'tokens', 'ignorecase');
        if ~isempty(atok)
            aval = str2double(atok{1}{1});
            if numel(atok{1}) >= 2
                aunit = lower(atok{1}{2});
            else
                aunit = '';   % optional unit group didn't participate in the match at all
            end
            switch aunit
                case 'mm', aval = aval/1e3;
                case 'cm', aval = aval/1e2;
                case 'um', aval = aval/1e6;
                case 'nm', aval = aval/1e9;
                case 'm'   % already meters
                otherwise  % no unit given -- assumed already meters
            end
            apothem_from_file = aval;
        else
            % CST "! Parameters = {Z0=...; a=...}" form: Z0 is used if the
            % option line gives R 0 (see after the loop).
            z0tok = regexp(line, '[{;]\s*Z0\s*=\s*([\d.eE+-]+)\s*[;}]', 'tokens', 'ignorecase');
            if ~isempty(z0tok)
                Z0_from_comment = str2double(z0tok{1}{1});
            end
            % Apothem from the same line, in cst_length_unit.
            atok_cst = regexp(line, '[{;]\s*a\s*=\s*([\d.eE+-]+)\s*[;}]', 'tokens', 'ignorecase');
            if ~isempty(atok_cst)
                apothem_from_file = str2double(atok_cst{1}{1}) * cst_unit_scale(cst_length_unit);
            end
        end
        continue
    end
    if line(1) == '#'
        tok = regexp(line, '#\s*(\S+)\s+S\s+(\S+)\s+R\s+([\d.eE+-]+)', 'tokens', 'ignorecase');
        if ~isempty(tok)
            freq_unit = tok{1}{1};
            fmt = upper(tok{1}{2});
            Z0 = str2double(tok{1}{3});
        end
        continue
    end
    data_lines{end+1} = line; %#ok<AGROW>
end
fclose(fid);

% CST writes "R 0" when the S-parameters are not renormalized; the
% reference impedance then comes from the project's "Z0" parameter.
if ~(isfinite(Z0) && Z0 > 0)
    if isfinite(Z0_from_comment) && Z0_from_comment > 0
        fprintf('local_read_touchstone: option line gives R=%g (CST: S-parameters not renormalized) -- using Z0=%g Ohm from the "! Parameters = {Z0=...}" comment instead.\n', Z0, Z0_from_comment);
        Z0 = Z0_from_comment;
    else
        error(['local_read_touchstone: the Touchstone option line gives a reference impedance of R=%g (CST writes 0 when the ' ...
               'S-parameters are not renormalized), and no "Z0=..." entry was found in the "! Parameters = {...}" comment either. ' ...
               'Either define a "Z0" parameter in the CST project, or re-export the Touchstone file with renormalization to a ' ...
               'known impedance.'], Z0);
    end
end

switch lower(freq_unit)
    case 'hz',  freq_scale = 1;
    case 'khz', freq_scale = 1e3;
    case 'mhz', freq_scale = 1e6;
    case 'ghz', freq_scale = 1e9;
    otherwise
        error('local_read_touchstone: unrecognized frequency unit "%s"', freq_unit);
end

all_vals = [];
for i = 1:numel(data_lines)
    toks = regexp(data_lines{i}, '[-+]?[\d.]+(?:[eE][-+]?\d+)?', 'match');
    all_vals = [all_vals, str2double(toks)]; %#ok<AGROW>
end

% Port count N comes from the file extension (e.g. .s2p -> N=2); the data
% length alone can be ambiguous.
[~, ~, ext] = fileparts(filename);
ext_tok = regexp(ext, '^\.[sS](\d+)[pP]$', 'tokens');
if isempty(ext_tok)
    error('local_read_touchstone: cannot determine port count from filename extension "%s" (expected .sNp, e.g. .s2p)', ext);
end
N = str2double(ext_tok{1}{1});
per_point = 1 + 2*N^2;
if mod(numel(all_vals), per_point) ~= 0
    error('local_read_touchstone: filename indicates a %d-port file (per_point=%d values), but the data (%d total values) does not divide evenly into that -- the file may be malformed or truncated.', ...
        N, per_point, numel(all_vals));
end

numf = numel(all_vals) / per_point;
all_vals = reshape(all_vals, per_point, numf);

Freq = all_vals(1,:).' * freq_scale;
Sraw = all_vals(2:end,:);   % (2*N^2 x numf)

S = zeros(N, N, numf);
if N == 2
    % 2-port files: S11, S21, S12, S22
    order_i = [1 2 1 2];
    order_j = [1 1 2 2];
else
    % N ~= 2: row-major, S11,S12,...,S1N,S21,...
    [order_j, order_i] = meshgrid(1:N, 1:N);
    order_i = reshape(order_i.', 1, []);
    order_j = reshape(order_j.', 1, []);
end
for p = 1:N^2
    v1 = Sraw(2*p-1, :); v2 = Sraw(2*p, :);
    switch fmt
        case 'MA'
            Sval = v1 .* exp(1j*v2*pi/180);
        case 'DB'
            Sval = 10.^(v1/20) .* exp(1j*v2*pi/180);
        case 'RI'
            Sval = v1 + 1j*v2;
        otherwise
            error('local_read_touchstone: unrecognized S-parameter format "%s"', fmt);
    end
    S(order_i(p), order_j(p), :) = Sval;
end
end