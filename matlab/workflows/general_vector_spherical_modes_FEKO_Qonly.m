%GENERAL_VECTOR_SPHERICAL_MODES_FEKO_QONLY  Single-port FEKO antenna -> spherical-mode network model from Q only.
%
% Variant of general_vector_spherical_modes_FEKO that uses only FEKO's
% spherical-mode Q coefficients (.sph export); no near- or far-field files
% are read, so there is no psi estimate and every result is Q-derived.
% Modes are ranked by accepted-power coupling |T_lm| and an (n_dom+1)-port
% Touchstone file is exported (port 1 = feed, ports 2..n_dom+1 = dominant
% spherical modes).
%
% The six field components (Er, Hr, Hphi, Htheta, Ephi, Etheta) are
% reconstructed from the dominant-mode Q coefficients on a user-chosen
% angular grid and radius (recon_theta/recon_phi/recon_r), at frequency f0
% only, with no reference trace.
%
% Inputs (FEKO exports, named <fileprefix>...):
%   .out                port parameters, input impedance, internal loads
%   _FarField1.sph      spherical-mode Q coefficients
%
% Outputs:
%   <snp_dir>/S_dom_Q_<antenna_name>.s<n_dom+1>p with ADS-compatible port
%   names; figures of impedance, field reconstruction, T_lm and radiated
%   power.

clear
close all
clc

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
% FEKO project file prefix (folder + antenna name).
%antenna_name="loop_noloads";
antenna_name="bowtie_edgeport_45degrees_PW_excitation";
fileprefix="C:\path\to\FEKO\project\"+antenna_name;

% Touchstone output folder.
snp_dir = fullfile(repo_root, 'output', 'SnP');
% =========================

%% --- Physical Constants ---
c0  = physconst('lightspeed');
mu0 = 4*pi*1e-7;          % Permeability of free space (H/m)
eps0 = 1/c0^2/mu0;        % Permittivity of free space (F/m)
eta0 = sqrt(mu0/eps0);    % Characteristic impedance of free space (Ohm)

%% --- Simulation Parameters ---
f0     = 2e9;         % Reference frequency for single-frequency plots (Hz)
fmin   = 0.1e9;       % Minimum frequency for plot x-axis (Hz)
fmax   = 2.5e9;         % Maximum frequency for plot x-axis (Hz)
lmax   = 9;           % Maximum spherical mode order to compute
n_dom  = 5;           % Number of dominant modes to identify and plot

% true: the far-field phase j^(l+1) is absorbed into QTMconv/QTEconv (via
% j_l_1), so the coupling formula uses jl1 = 1. false: the coupling formula
% applies jl1 = (1j)^(l+1).
use_jl1_convention = true;

% true: apply exp(-j*k*a) to move the S-parameter reference plane from the
% feed to r = a (antenna apothem). false: reference stays at the feed.
apply_ka_deembedding = true;

% SnP port ordering (ports 2..n_dom+1):
%   'peak_freq'  - by frequency of peak |T_lm|, ascending
%   'mode_index' - by FEKO mode index j = 2*(l*(l+1)+m-1)+s, s=1 TE, s=2 TM
%snp_order_method = 'peak_freq';
snp_order_method = 'mode_index';

% DC (0 Hz) point of the SnP export:
%   'auto' - internal-load DC-conduction model if the .out file has an
%            internal load with nonzero resistance, else 'Rmin'
%   'Rmin' - reactive-limit approximation from the lowest-frequency Zant
%   'user' - prompt for S11 at DC
dc_resistance_method = 'auto';

% Angular grid for the field reconstruction plots (1-degree steps over the
% full sphere). The radius recon_r is set below to a + lambda/(2*pi) at the
% lowest simulated frequency, just outside the reactive near field.
recon_theta = (0:1:180)   * pi/180;   % (1 x 181)
recon_phi   = (-180:1:180) * pi/180;  % (1 x 361)

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

k = 2*pi*Freq/c0;

% Reconstruction radius and plot frequency/angle indices.
recon_r = a + c0/2/pi/min(Freq);
[~,index_f]   = min(abs(Freq-f0));
[~,index_phi] = min(abs(recon_phi - 0));
[~,index_th]  = min(abs(recon_theta - pi/2));
numth = length(recon_theta); numph = length(recon_phi);

%% --- Plot Antenna Impedance vs Frequency ---
fignum = 1; figure(fignum)
plot(Freq,real(Zant),"Color","b","LineWidth",LW)
hold on
plot(Freq,imag(Zant),"Color","r","LineWidth",LW)
hold off
title("Impedance","FontSize",FS)
legend("$R$","$X$","FontSize",FS)
xlim([fmin,fmax])
ylim([-eta0,eta0])
grid on

%% --- Read FEKO Q Coefficients (.sph) ---
% Mode-index grids and normalisation (independent of the Q data).
lval = (1:1:lmax).';   % mode orders l = 1..lmax
mval = (-lmax:lmax);   % m ranges from -lmax to +lmax
L    = length(lval);
M    = length(mval);
midx_0  = find(mval == 0);

l_factor = lval .* (lval + 1);
if use_jl1_convention
    j_l_1 = (1j) .^ (lval + 1);
else
    j_l_1 = ones(L, 1);
end

Norm=zeros(L,M);
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

% Read the Q coefficients (.sph) and convert them to the psi convention.
Q = []; QTEconv = zeros(L,M,numf); QTMconv = zeros(L,M,numf);
freqs_Q = Freq; mmax = lmax; num_modes = 2*(lmax*(lmax+1)+mmax-1)+2;   % safe fallbacks
use_Q_method = true;
try
    FF_request = "_FarField1.sph";
    [Q, ~, ~, mmax, lmax_Q_file, freqs_Q] = ReadSphModes_freqs(fileprefix+FF_request);
    freqs_Q=freqs_Q(Zindx);

    % The Q file must contain at least lmax mode orders; otherwise the
    % missing higher-order QTMconv/QTEconv would silently be zero.
    if lmax_Q_file < max(lval)
        error("Q file max mode order l=%d is less than the script's lmax=%d", lmax_Q_file, max(lval))
    end

    if any(freqs_Q ~= Freq)
        error("frequencies from Q file not the same as field data")
    end
    Q=Q(:,Zindx);
catch ME
    use_Q_method = false;
    fprintf('general_vector_spherical_modes_Qonly: could not read Q coefficients (%s) — no spherical-mode data available.\n', ME.message);
end

if use_Q_method
    % Conversion from FEKO Q (TICRA convention) to psi. conj(Q) converts
    % TICRA's exp(-iwt) time convention to the exp(+jwt) convention used
    % here. QTMconv/QTEconv then take the place of psiTM/psiTE throughout.
    conv_factor_TE = 1j*sqrt(8*pi) .* (1j.*j_l_1) * 2^(symmetry_planes.physical_gnd);
    conv_factor_TM =    sqrt(8*pi) .* (1j.*j_l_1) * 2^(symmetry_planes.physical_gnd);

    for lindx=1:L
        for mindx=1:M
            if abs(mval(mindx))>lval(lindx)
                continue
            else
                if mval(mindx)==0
                    QTMconv(lindx,mindx,:)=conj(Q(SphModeIndex(lval(lindx),mval(mindx),'TM'),:)) .* conv_factor_TM(lindx) ;
                    QTEconv(lindx,mindx,:)=conj(Q(SphModeIndex(lval(lindx),mval(mindx),'TE'),:)) .* conv_factor_TE(lindx) ;
                elseif mval(mindx)<0
                    QTMconv(lindx,mindx,:)=(CSphase(mindx))*(1j/sqrt(2))*...
                        conj( Q(SphModeIndex(lval(lindx),mval(mindx),'TM'),:)-(-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),-mval(mindx),'TM'),:) )...
                        .* conv_factor_TM(lindx) ;
                    QTEconv(lindx,mindx,:)=(CSphase(mindx))*(1j/sqrt(2))*...
                        conj( Q(SphModeIndex(lval(lindx),mval(mindx),'TE'),:)-(-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),-mval(mindx),'TE'),:) )...
                        .* conv_factor_TE(lindx) ;
                elseif mval(mindx)>0
                    QTMconv(lindx,mindx,:)=(CSphase(mindx))*(1/sqrt(2))*...
                        conj( (-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),mval(mindx),'TM'),:)+Q(SphModeIndex(lval(lindx),-mval(mindx),'TM'),:) )...
                        .* conv_factor_TM(lindx) ;
                    QTEconv(lindx,mindx,:)=(CSphase(mindx))*(1/sqrt(2))*...
                        conj( (-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),mval(mindx),'TE'),:)+Q(SphModeIndex(lval(lindx),-mval(mindx),'TE'),:) )...
                        .* conv_factor_TE(lindx) ;
                end
            end
        end
    end
end

if ~use_Q_method
    error('general_vector_spherical_modes_Qonly: no Q data available -- this script requires the .sph export to run at all (unlike the full script, there is no near/far-field fallback here).')
end

%% --- Build Spherical Wave S-parameters for EVERY (l,m,type), for Ranking ---
% Every mode is ranked by its accepted-power coupling
% T_lm = S_{mode,1}/sqrt(1-|S_00|^2), which, unlike raw |Q|, accounts for
% how much power the antenna actually accepts at each frequency.
%   S_{mode,1} = sqrt(Z0/2^|even_odd_excitation|) * (1+S_00) * (Q/V0) * exp(-jka)/j^(l+1)
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
    jl1_l = ones(L,1);                                   % phase already absorbed into QTMconv/QTEconv via j_l_1
else
    jl1_l = (1j).^(lval(:) + 1);                         % (L x 1)
end

one_plus_S00_r = reshape(1 + S00_all, 1, 1, numf);
jl1_r          = reshape(jl1_l, L, 1, 1);                % broadcasts over M and numf

S_TM_all = port_factor .* one_plus_S00_r .* (QTMconv ./ V0(1)) .* ka_factor_r ./ jl1_r;   % (L x M x numf)
S_TE_all = port_factor .* one_plus_S00_r .* (QTEconv ./ V0(1)) .* ka_factor_r ./ jl1_r;

T_denom_all = sqrt(1 - abs(S00_all).^2);                 % (1 x numf)
T_TM_all = S_TM_all ./ reshape(T_denom_all, 1, 1, numf); % (L x M x numf)
T_TE_all = S_TE_all ./ reshape(T_denom_all, 1, 1, numf);

%% --- Identify Dominant Modal Coefficients across TE and TM jointly ---
% TE and TM modes are pooled and the n_dom largest by max |T_lm| within
% [fmin, fmax] are selected.
fband_idx = find(Freq >= fmin & Freq <= fmax);

% Max |T_lm| within [fmin,fmax] for each (l,m) pair, for each type.
[QTM_mag, QTM_fmax_band] = max(abs(reshape(T_TM_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
[QTE_mag, QTE_fmax_band] = max(abs(reshape(T_TE_all(:,:,fband_idx), L*M, length(fband_idx))), [], 2);
QTM_fmax  = fband_idx(QTM_fmax_band(:)).';
QTE_fmax  = fband_idx(QTE_fmax_band(:)).';

rank_TM_mag = QTM_mag;   rank_TM_fmax = QTM_fmax;
rank_TE_mag = QTE_mag;   rank_TE_fmax = QTE_fmax;

% Replace NaN with 0 so invalid (|m|>l) modes sort to the bottom
rank_TM_mag(isnan(rank_TM_mag)) = 0;
rank_TE_mag(isnan(rank_TE_mag)) = 0;

% Pool into a single list: [TM entries; TE entries], tagging each with type
% type_tag: 0 = TM, 1 = TE
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
fprintf('\n--- Top %d dominant modes (TE and TM ranked jointly by max |T_lm|, accepted-power coupling) ---\n', n_dom)
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

%% --- Reconstruct Fields from Q Coefficients, and Plot ---
% Same vector spherical wave formulas as general_vector_spherical_modes_FEKO,
% with QTMconv/QTEconv in place of psiTM/psiTE, for the n_dom dominant
% modes at frequency index_f only. Each plot has a single trace (reference
% and comparison traces passed as []).
[PHI, THETA] = meshgrid(recon_phi, recon_theta);
kr_f = k(index_f) * recon_r;

gradYlm_th_dom = zeros(numth,numph,n_dom); gradYlm_ph_dom = gradYlm_th_dom;
jrFlm_r_dom = zeros(numth,numph,n_dom);
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    [gradYlm_th_dom(:,:,kk), gradYlm_ph_dom(:,:,kk)] = gradYlm(lval(li), mval(mi), THETA, PHI, have_CSphase);
    [~, ~, jrFlm_r_dom(:,:,kk)] = Spherical_Mode_Field(lval(li), mval(mi), 1, THETA, PHI, kr_f);
end

field_title_base = "@r="+num2str(recon_r*1000,2)+"mm, "+ ...
                   num2str(Freq(index_f)/1e9,2)+"GHz, \#dom="+num2str(n_dom)+", $Q$-derived";
phi_label   = "$\phi="+num2str(recon_phi(index_phi)*180/pi,3)+"^\circ$";
theta_label = "$\theta="+num2str(recon_theta(index_th)*180/pi,3)+"^\circ$";

%% --- Er, Hr ---
Er_modes = zeros(numth,numph,n_dom); Hr_modes = Er_modes;
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0   % TM -- contributes to Er only
        Er_modes(:,:,kk) = sqrt(eta0)*QTMconv(li,mi,index_f)*conj(jrFlm_r_dom(:,:,kk))/(1j*recon_r);
    else                    % TE -- contributes to Hr only
        Hr_modes(:,:,kk) = (1/sqrt(eta0))*QTEconv(li,mi,index_f)*conj(jrFlm_r_dom(:,:,kk))/(1j*recon_r);
    end
end
Hr_total = sum(Hr_modes, 3);
Er_total = sum(Er_modes, 3);
clearvars Er_modes Hr_modes

fignum = plot_field_comparison_polar(fignum, ...
    recon_theta/pi, Hr_total(:,index_phi), [], [], ...
    recon_phi/pi,   Hr_total(index_th,:),  [], [], ...
    eta0, '$\eta_0|H_r|$ (V/m)', {"(dominant-modes) $H_r$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Hr_total

fignum = plot_field_comparison_polar(fignum, ...
    recon_theta/pi, Er_total(:,index_phi), [], [], ...
    recon_phi/pi,   Er_total(index_th,:),  [], [], ...
    1, '$|E_r|$ (V/m)', {"(dominant-modes) $E_r$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Er_total

%% --- H_phi ---
Hph_TM = zeros(numth,numph,n_dom); Hph_TE = Hph_TM;
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0
        Hph_TM(:,:,kk) = (1/sqrt(eta0))*QTMconv(li,mi,index_f)*sphericalRicatti_h(lval(li),2,kr_f)/j_l_1(li) ...
                          *(-gradYlm_th_dom(:,:,kk)/recon_r);
    else
        Hph_TE(:,:,kk) = (1j/sqrt(eta0))*QTEconv(li,mi,index_f)*sphericalRicatti_h_prime(lval(li),2,kr_f)/j_l_1(li) ...
                          *(-gradYlm_ph_dom(:,:,kk)/recon_r);
    end
end
Hph_total = sum(Hph_TM + Hph_TE, 3);
clearvars Hph_TM Hph_TE

fignum = plot_field_comparison_polar(fignum, ...
    recon_theta/pi, Hph_total(:,index_phi), [], [], ...
    recon_phi/pi,   Hph_total(index_th,:),  [], [], ...
    eta0, '$\eta_0|H_\phi|$ (V/m)', {"(dominant-modes) $H_\phi$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Hph_total

%% --- H_theta ---
Hth_TM = zeros(numth,numph,n_dom); Hth_TE = Hth_TM;
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0
        Hth_TM(:,:,kk) = (1/sqrt(eta0))*QTMconv(li,mi,index_f)*sphericalRicatti_h(lval(li),2,kr_f)/j_l_1(li) ...
                          *(+gradYlm_ph_dom(:,:,kk)/recon_r);
    else
        Hth_TE(:,:,kk) = (1j/sqrt(eta0))*QTEconv(li,mi,index_f)*sphericalRicatti_h_prime(lval(li),2,kr_f)/j_l_1(li) ...
                          *(-gradYlm_th_dom(:,:,kk)/recon_r);
    end
end
Hth_total = sum(Hth_TM + Hth_TE, 3);
clearvars Hth_TM Hth_TE

fignum = plot_field_comparison_polar(fignum, ...
    recon_theta/pi, Hth_total(:,index_phi), [], [], ...
    recon_phi/pi,   Hth_total(index_th,:),  [], [], ...
    eta0, '$\eta_0|H_\theta|$ (V/m)', {"(dominant-modes) $H_\theta$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Hth_total

%% --- E_phi ---
Eph_TM = zeros(numth,numph,n_dom); Eph_TE = Eph_TM;
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 1   % TE
        Eph_TE(:,:,kk) = sqrt(eta0)*QTEconv(li,mi,index_f)*sphericalRicatti_h(lval(li),2,kr_f)/j_l_1(li) ...
                          *(+gradYlm_th_dom(:,:,kk)/recon_r);
    else                    % TM
        Eph_TM(:,:,kk) = 1j*sqrt(eta0)*QTMconv(li,mi,index_f)*sphericalRicatti_h_prime(lval(li),2,kr_f)/j_l_1(li) ...
                          *(-gradYlm_ph_dom(:,:,kk)/recon_r);
    end
end
Eph_total = sum(Eph_TM + Eph_TE, 3);
clearvars Eph_TM Eph_TE

fignum = plot_field_comparison_polar(fignum, ...
    recon_theta/pi, Eph_total(:,index_phi), [], [], ...
    recon_phi/pi,   Eph_total(index_th,:),  [], [], ...
    1, '$|E_\phi|$ (V/m)', {"(dominant-modes) $E_\phi$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Eph_total

%% --- E_theta ---
Eth_TM = zeros(numth,numph,n_dom); Eth_TE = Eth_TM;
for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 1   % TE
        Eth_TE(:,:,kk) = sqrt(eta0)*QTEconv(li,mi,index_f)*sphericalRicatti_h(lval(li),2,kr_f)/j_l_1(li) ...
                          *(-gradYlm_ph_dom(:,:,kk)/recon_r);
    else                    % TM
        Eth_TM(:,:,kk) = 1j*sqrt(eta0)*QTMconv(li,mi,index_f)*sphericalRicatti_h_prime(lval(li),2,kr_f)/j_l_1(li) ...
                          *(-gradYlm_th_dom(:,:,kk)/recon_r);
    end
end
Eth_total = sum(Eth_TM + Eth_TE, 3);
clearvars Eth_TM Eth_TE gradYlm_th_dom gradYlm_ph_dom jrFlm_r_dom

fignum = plot_field_comparison_polar(fignum, ...
    recon_theta/pi, Eth_total(:,index_phi), [], [], ...
    recon_phi/pi,   Eth_total(index_th,:),  [], [], ...
    1, '$|E_\theta|$ (V/m)', {"(dominant-modes) $E_\theta$ reconstructed", field_title_base}, FS, LW, phi_label, theta_label);
clearvars Eth_total

%% --- Build Spherical Wave S-parameters for Dominant Modes ---
% S_Q is (n_dom+1) x (n_dom+1) x numf.
% Index 1: antenna port (S(1,1,:) = S_00); indices 2..n_dom+1: dominant
% modes in rank order. Values are slices of S_TM_all/S_TE_all above.

S_Q = zeros(n_dom+1, n_dom+1, numf);
S_Q(1,1,:) = S00_all;

for kk = 1:n_dom
    li = dom_lidx(kk); mi = dom_midx(kk);
    if dom_type(kk) == 0   % TM
        coupling_Q = reshape(S_TM_all(li,mi,:), 1, numf);
    else                   % TE
        coupling_Q = reshape(S_TE_all(li,mi,:), 1, numf);
    end
    S_Q(kk+1, 1, :) = coupling_Q;
    S_Q(1, kk+1, :) = coupling_Q;
end
% Sanitise: replace any remaining NaN with 0 before export
S_Q(isnan(S_Q)) = 0;

%% --- Plot Dominant Mode S-parameters (Accepted Power Coupling T_lm) ---
% T_lm = S_{mode,1} / sqrt(1 - |S_00|^2).
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

T_Q = zeros(n_dom, numf);
denom_Q = sqrt(1 - abs(reshape(S_Q(1,1,:), 1, numf)).^2);
for kk = 1:n_dom
    T_Q(kk,:) = reshape(S_Q(kk+1,1,:), 1, numf) ./ denom_Q;
end

% snp_order: ordering of the n_dom dominant modes for SnP ports 2..n_dom+1.
if strcmpi(snp_order_method, 'mode_index')
    % FEKO mode index: j_ = 2*(l*(l+1) + m - 1) + s, s=1 for TE, s=2 for TM
    dom_l = lval(dom_lidx);
    dom_m = mval(dom_midx);
    s_val = ones(n_dom,1); s_val(dom_type == 0) = 2;   % TM=2, TE=1
    dom_j = 2*(dom_l.*(dom_l+1) + dom_m(:) - 1) + s_val;
    [~, snp_order] = sort(dom_j, 'ascend');
else   % 'peak_freq' (default)
    [~, T_fmax_band] = max(abs(T_Q(:,fband_idx)), [], 2);
    T_fmax = fband_idx(T_fmax_band);
    [~, snp_order] = sort(T_fmax, 'ascend');
end

fignum = fignum+1; figure(fignum);
hold on
plot(Freq, db(reshape(S_Q(1,1,:),1,numf)), "LineStyle","--", "LineWidth",LW, ...
    "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{00}$")
for snp_k = 1:n_dom
    kk = snp_order(snp_k);
    c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
    plot(Freq, db(T_Q(kk,:)), "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
end
hold off
ylabel("$|T_{lm}|$ dB ($Q$)","FontSize",FS)
ylim([round(min(max(db(T_Q(:,fband_idx)), [], 2)))-5 0])
title("Accepted Power Coupling $T_{lm} = S_{lm}/\sqrt{1-|S_{00}|^2}$, $Q$-derived","FontSize",FS)
legend("FontSize",FS)
xlabel("Frequency (Hz)","FontSize",FS)
grid on; xlim([fmin,fmax])

fignum = fignum+1; figure(fignum);
hold on
plot(Freq, angle(reshape(S_Q(1,1,:),1,numf))*180/pi, "LineStyle","--", "LineWidth",LW, ...
    "Color",get(groot,'defaultAxesXColor'), "DisplayName","$S_{00}$")
for snp_k = 1:n_dom
    kk = snp_order(snp_k);
    c = clr_S(mod(snp_k-1, n_clr_S)+1, :);
    plot(Freq, angle(T_Q(kk,:))*180/pi, "LineStyle","-", "LineWidth",LW, "Color",c, "DisplayName",mode_labels{kk})
end
hold off
ylabel("$\angle T_{lm}$ deg ($Q$)","FontSize",FS)
ylim([-180 180])
title("Accepted Power Coupling Phase $\angle T_{lm}$, $Q$-derived","FontSize",FS)
legend("FontSize",FS)
xlabel("Frequency (Hz)","FontSize",FS)
grid on; xlim([fmin,fmax])

%% --- Export Spherical Wave S-parameters to Touchstone File ---
% Ports 2..n_dom+1 follow snp_order (see snp_order_method).
if ~exist(snp_dir, 'dir'), mkdir(snp_dir); end

snp_port_order = [1, snp_order(:).' + 1];
S_snp = S_Q(snp_port_order, snp_port_order, :);

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

snp_path = fullfile(snp_dir, "S_dom_Q_"+antenna_name+".s"+num2str(n_dom+1)+"p");
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
comment_lines = comment_lines + "! Written by " + string(mfilename) + ".m (Q-derived)," +...
                " a="+num2str(a,3)+"m, Z0="+num2str(Z0,5)+"ohms, ka_factor="+ka_factor_str + newline;
comment_lines = comment_lines + "! Date : " + string(datetime('now')) + newline;

% Port names as "! Port[n] name" comments, the form ADS writes in its own
% Touchstone V1.1 files.
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

%% --- Radiated Power (Q-derived only) ---
% Total |Q|^2 power, dominant-mode |Q|^2 power and S-parameter power from
% S_Q. Without near-field data the reference is the all-modes |Q|^2 power.
mask_TM_dom = false(L, M);
mask_TE_dom = false(L, M);
for kk = 1:n_dom_TM
    mask_TM_dom(dom_TM_lidx(kk), dom_TM_midx(kk)) = true;
end
for kk = 1:n_dom_TE
    mask_TE_dom(dom_TE_lidx(kk), dom_TE_midx(kk)) = true;
end

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

Pr_SQ = zeros(numf,1);
for f = 1:numf
    % The (1+S_00) factor cancels the one in S_Q(mode,1,f), so Pr_SQ depends
    % only on Q and V0.
    port_voltage_factor = abs(V0(1) / (1 + S_Q(1,1,f)))^2;
    Pr_SQ(f) = (1/2) * Y0 * port_voltage_factor * 2^abs(even_odd_excitation) * ...
               sum(abs(S_Q(2:end,1,f)).^2);
end

fignum = fignum+1; figure(fignum);
hold on
plot(Freq, Pr_Q,    "--", "LineWidth",LW)
plot(Freq, Pr_Qdom, "^",  "LineWidth",LW, "MarkerSize",3)
plot(Freq, Pr_SQ,   "d",  "LineWidth",LW, "MarkerSize",4)
hold off
legend("$|Q|^2$ all","$|Q|^2$ dom","$|S|^2$ $Q$","FontSize",FS)
title("Radiated Power, $Q$-derived","FontSize",FS)
ylabel("Watts","FontSize",FS); xlabel("Freq (Hz)","FontSize",FS)
grid on; xlim([fmin,fmax])

%% --- L2 Norm Error of Dominant-Mode Power vs All-Modes Power ---
% Relative L2 error over [fmin, fmax] of the dominant-mode and
% S-parameter radiated power against the all-modes |Q|^2 power Pr_Q. This
% measures how much of the full Q expansion the n_dom truncation captures,
% not absolute accuracy (no near-field Poynting reference is available).
Pr_Q_band  = real(Pr_Q(fband_idx));
Pr_ref_norm = norm(Pr_Q_band);

L2_err_Qdom = norm(Pr_Qdom(fband_idx) - Pr_Q_band) / Pr_ref_norm;
L2_err_SQ   = norm(Pr_SQ(fband_idx)   - Pr_Q_band) / Pr_ref_norm;

fprintf('\nRadiated power L2 error vs |Q|^2 all-modes power (%.3g-%.3g GHz band):\n', fmin/1e9, fmax/1e9);
fprintf('  |Q|^2 dom : %6.2f%%\n', 100*L2_err_Qdom);
fprintf('  |S|^2 Q   : %6.2f%%\n', 100*L2_err_SQ);
