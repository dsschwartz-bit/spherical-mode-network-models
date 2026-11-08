%PLOT_CST_FARFIELDS  3D vector plot and pattern cuts of CST far-field exports.
%
% Reads the CST far-field exports of one port with the same readers used by
% general_vector_spherical_modes_loadports_CST.m (discover_cst_farfield_files,
% read_cst_farfield, cst_unit_scale) and, for each requested frequency, plots
%   1. a 3D cone plot of the far-field E vector over the sphere
%      (plot_cst_farfield_vector), and
%   2. theta and phi cuts of |E_total| = sqrt(|E_phi|^2 + |E_theta|^2)
%      (plot_cst_farfield_cuts).
%
% CST far-field convention: the export is r*E evaluated at r = 1 m with the
% e^(-jkr) phase still included. It is removed by multiplying by exp(+jk),
% as in compute_port_FF_coeffs. This only changes the phase (cone-plot
% snapshot), not the cut magnitudes.
%
% Input : CST far-field ASCII exports "<...>f=<freq>[<port>]<...>" in
%         farfield_directory (CST project folder\Export\Farfield).
% Output: figures only.

clear
close all
clc
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ (incl. matlab/io/cst) and third_party/

%% ===== USER SETTINGS =====
farfield_directory = "C:\path\to\CST\dipole_project\Export\Farfield";   % CST far-field export folder
% =========================

cst_freq_unit    = "GHz";   % units of the "f=<freq>" token in each filename
cst_voltage_unit = "V";     % units of the exported field columns (see read_cst_farfield)
port_number      = 1;       % port whose far fields are plotted (the "[<port>]" token in each filename)

requested_freqs_GHz = [1, 3, 5];   % frequencies to plot, in cst_freq_unit's units
freq_tol_GHz = 0.05;                % warn if the nearest available file is farther
                                     % from a requested frequency than this

snapshot_mode = 'real';    % 'real' (t=0, Re{phasor}) or 'negimag' (omega*t=pi/2, -Im{phasor})
n_theta_ds = 21;           % cone-plot angular down-sampling (a finer grid
n_phi_ds   = 21;           % makes the cone plot unreadable)

%% --- Plot Settings ---
LW = 2; FS = 14;
set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% --- Physical Constants ---
c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0);

%% --- Discover and Select Far-Field Files ---
[cst_farfield_files, cst_farfield_freqs_raw] = discover_cst_farfield_files(farfield_directory);
if port_number > numel(cst_farfield_files)
    error('plot_cst_farfields: port_number=%d requested, but only %d port(s) found in "%s".', ...
        port_number, numel(cst_farfield_files), farfield_directory);
end
freq_scale = cst_unit_scale(cst_freq_unit);
available_freqs_Hz = cst_farfield_freqs_raw * freq_scale;   % ascending, matches cst_farfield_files{p} order

requested_freqs_Hz = requested_freqs_GHz * 1e9;
fignum = 0;

for rf = 1:numel(requested_freqs_Hz)
    f_req = requested_freqs_Hz(rf);
    [f_err, f_idx] = min(abs(available_freqs_Hz - f_req));
    if f_err > freq_tol_GHz*1e9
        fprintf('plot_cst_farfields: WARNING - requested %.6g GHz, nearest available file is %.6g GHz (%.3g GHz away) -- plotting that instead.\n', ...
            f_req/1e9, available_freqs_Hz(f_idx)/1e9, f_err/1e9);
    end
    f_actual = available_freqs_Hz(f_idx);
    k_actual = 2*pi*f_actual/c0;

    filename = fullfile(farfield_directory, cst_farfield_files{port_number}(f_idx));
    [rEpt, theta, phi] = read_cst_farfield(filename, cst_voltage_unit);

    % remove the e^(-jk) phase of CST's r = 1 m export (see header)
    rEpt = rEpt * exp(1j*k_actual*1);

    Eph = rEpt(:,:,1);   % rEpt dim3: 1=rE_phi, 2=rE_theta (read_cst_farfield's convention)
    Eth = rEpt(:,:,2);

    title_str = sprintf('Port %d, %.4g GHz', port_number, f_actual/1e9);

    %% --- 3D Vector Plot (E-field) ---
    % field(:,:,1/2/3) = Vr/Vth/Vph (plot_cst_farfield_vector convention);
    % Vr = 0 in the far field.
    field_E = zeros(numel(theta), numel(phi), 3);
    field_E(:,:,2) = Eth;
    field_E(:,:,3) = Eph;

    fignum = plot_cst_farfield_vector(fignum, 1, theta, phi, field_E, ...
        n_theta_ds, n_phi_ds, 1, snapshot_mode, '$\vec{E}$', title_str, FS, LW);

    %% --- Theta/Phi Cuts of Total Field Magnitude ---
    [~, index_phi_cut] = min(abs(phi - 0));         % phi index closest to 0 rad
    [~, index_th_cut]  = min(abs(theta - pi/2));    % theta index closest to pi/2 (broadside)

    fignum = plot_cst_farfield_cuts(fignum, theta, phi, Eph, Eth, ...
        index_phi_cut, index_th_cut, 1, '$|E_{\mathrm{total}}|$ (V/m)', title_str, FS, LW);
end
