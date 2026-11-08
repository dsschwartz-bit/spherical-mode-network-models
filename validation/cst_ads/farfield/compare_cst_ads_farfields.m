%COMPARE_CST_ADS_FARFIELDS  Compare CST and ADS far-field patterns at selected frequencies.
%
% For each requested frequency, draws a polar figure of |E_theta| and |E_phi|
% (rows) over the requested theta and phi cuts (columns), CST against ADS.
% (A 2x2 3D cone-plot comparison, CST/ADS x t=0/quarter period, is available
% as a commented-out block in the frequency loop.)
%
% Inputs:
%   CST : far-field ASCII exports of one port in
%         <cst_root>\<cst_project>\Export\Farfield, read with
%         discover_cst_farfield_files / read_cst_farfield / cst_unit_scale.
%   ADS : DataSet CSV export <ads_dir>\<antenna_name>_FF_export.csv with
%         columns phi, theta, freq, mag/phase(Vtheta), mag/phase(Vphi),
%         read with read_ads_farfield_csv (raw columns) and reshaped here
%         into (theta x phi x freq) complex Vtheta / Vphi.
% Output: figures only.
%
% Both CST and ADS report r*E at r = 1 m with the e^(-jkr) phase still
% included; both are multiplied by exp(+jk) on reading, as in
% compute_port_FF_coeffs and plot_cst_farfields.m.
% Each source is plotted on its own angular grid (no interpolation): the
% requested cut angle is matched to each source's nearest sample.

clear
close all
clc
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ (incl. matlab/io/cst) and third_party/

%% ===== USER SETTINGS =====
cst_root = "C:\path\to\CST\";            % folder containing the CST project folders
ads_dir  = "C:\path\to\ADS\workspace\data\";   % ADS workspace data folder (CSV exports)
% =========================

antenna_name = "rectangularpatch";
cst_variant  = "symmetry_FDTD";     % CST model / solver variant
cst_project  = antenna_name + "_" + cst_variant;   % CST project folder

%% --- Configuration: CST ---
cst_farfield_directory = string(fullfile(cst_root, cst_project, "Export", "Farfield"));
cst_freq_unit    = "GHz";   % units of the "f=<freq>" token in each CST filename
cst_voltage_unit = "V";     % units of CST's exported field columns
cst_port_number  = 1;       % which CST port's far fields to compare

%% --- Configuration: ADS ---
ads_csv_file = string(fullfile(ads_dir, antenna_name + "_FF_export.csv"));

%% --- Comparison Frequencies ---
requested_freqs_GHz = [1,2,2.5,3,4];   % frequencies to compare, in GHz
freq_tol_GHz = 0.05;                % warn if either source's nearest available
                                     % frequency is farther from a request than this

n_theta_ds = 21;           % cone-plot angular down-sampling, applied independently
n_phi_ds   = 21;           % to each source's own grid (21x21 is a good start)
cst_cut_ds_step_deg = 5;   % angular spacing (degrees) between CST markers in each
                           % polar cut plot, drawn as solid black triangle markers;
                           % ADS is NOT downsampled (drawn as a continuous red line
                           % at its own full resolution)

plot_voltage_unit = "V";   % unit for the polar-cut radial axes: "V", "mV" or "uV"
                            % (scales the plotted magnitudes by 1, 1e3 or 1e6)

% Polar cuts. Each entry adds one COLUMN to the polar figure (left to right:
% all theta cuts first, then all phi cuts). Row 1 is |E_theta|, row 2 is |E_phi|.
% Each cut angle is matched to each source's own nearest sample independently.
%   theta cut = constant phi, theta varies. A theta cut at phi0 is the full
%               plane through z: phi0 is drawn at +theta and phi0+180 (wrapped
%               to 0..360) at -theta, so DON'T also list phi0+180 here.
%   phi cut   = constant theta, phi varies over 0..360.
phi_cuts_deg   = [0];       % theta cuts at these fixed phi (deg), e.g. [0 90]
theta_cuts_deg = [90 45];   % phi cuts at these fixed theta (deg)
cut_tol_deg    = 1;         % a panel is NOT plotted (text notice instead) if either
                            % source's nearest grid sample is farther than this
                            % from the requested cut angle (deg)

%% --- Plot Settings ---
LW = 2; FS = 18; legend_FS = 16;   % legend font size, separate from FS
set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

switch plot_voltage_unit
    case "V",  plot_mag_scale = 1;   plot_unit_label = "V";
    case "mV", plot_mag_scale = 1e3; plot_unit_label = "mV";
    case "uV", plot_mag_scale = 1e6; plot_unit_label = "$\mu$V";
    otherwise
        error('compare_cst_ads_farfields: plot_voltage_unit must be "V", "mV" or "uV" (got "%s").', plot_voltage_unit);
end

%% --- Physical Constants ---
c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0); %#ok<NASGU>

%% --- Discover CST Far-Field Files ---
[cst_farfield_files, cst_farfield_freqs_raw] = discover_cst_farfield_files(cst_farfield_directory);
if cst_port_number > numel(cst_farfield_files)
    error('compare_cst_ads_farfields: cst_port_number=%d requested, but only %d port(s) found in "%s".', ...
        cst_port_number, numel(cst_farfield_files), cst_farfield_directory);
end
cst_freq_scale = cst_unit_scale(cst_freq_unit);
cst_available_freqs_Hz = cst_farfield_freqs_raw * cst_freq_scale;   % ascending

%% --- Read the Full ADS Dataset (all frequencies at once) ---
% read_ads_farfield_csv returns raw (n x 1) columns; they are reshaped into
% a (theta x phi x freq) grid and mag/phase pairs combined below.
ads_data = read_ads_farfield_csv(ads_csv_file);

% Angle columns: phi_deg / theta_deg (degrees, sanitized to phideg /
% thetadeg) or phi / theta (radians); either way converted to radians in
% ads_data.phi / ads_data.theta.
if isfield(ads_data, 'phideg') && isfield(ads_data, 'thetadeg')
    ads_data.phi   = ads_data.phideg   * pi/180;
    ads_data.theta = ads_data.thetadeg * pi/180;
elseif ~(isfield(ads_data, 'phi') && isfield(ads_data, 'theta'))
    error('compare_cst_ads_farfields: "%s" has no phi/theta (radians) or phi_deg/theta_deg columns (found: %s).', ...
        ads_csv_file, strjoin(ads_data.column_names, ', '));
end

ads_phi    = unique(ads_data.phi);      % ascending
ads_theta  = unique(ads_data.theta);    % ascending
ads_freq_Hz = unique(ads_data.freq);    % ascending, already in Hz
numph_ads = numel(ads_phi); numth_ads = numel(ads_theta); numf_ads = numel(ads_freq_Hz);

if numph_ads*numth_ads*numf_ads ~= numel(ads_data.phi)
    error('compare_cst_ads_farfields: ADS data in "%s" does not form a complete (phi x theta x freq) grid: %d x %d x %d = %d combinations expected, but %d rows found.', ...
        ads_csv_file, numph_ads, numth_ads, numf_ads, numph_ads*numth_ads*numf_ads, numel(ads_data.phi));
end

[~, ads_pidx] = ismember(ads_data.phi,   ads_phi);
[~, ads_tidx] = ismember(ads_data.theta, ads_theta);
[~, ads_fidx] = ismember(ads_data.freq,  ads_freq_Hz);
ads_lin_idx = sub2ind([numth_ads, numph_ads, numf_ads], ads_tidx, ads_pidx, ads_fidx);
if numel(unique(ads_lin_idx)) ~= numel(ads_data.phi)
    error('compare_cst_ads_farfields: ADS data in "%s" has duplicate (phi,theta,freq) rows -- the grid is not well-defined.', ads_csv_file);
end

% ADS phase(...) columns are in degrees, whatever the units of the angle
% sweep.
ads_Vtheta_all = nan(numth_ads, numph_ads, numf_ads);
ads_Vphi_all   = nan(numth_ads, numph_ads, numf_ads);
ads_Vtheta_all(ads_lin_idx) = ads_data.magVtheta .* exp(1j*ads_data.phaseVtheta*pi/180);
ads_Vphi_all(ads_lin_idx)   = ads_data.magVphi   .* exp(1j*ads_data.phaseVphi*pi/180);
if any(isnan(ads_Vtheta_all(:))) || any(isnan(ads_Vphi_all(:)))
    error('compare_cst_ads_farfields: some (phi,theta,freq) combinations are missing from "%s".', ads_csv_file);
end

%% --- Compare at Each Requested Frequency ---
requested_freqs_Hz = requested_freqs_GHz * 1e9;
fignum = 0;
phi_cut_rad   = phi_cuts_deg(:).'   * pi/180;   % row vectors of cut angles (rad)
theta_cut_rad = theta_cuts_deg(:).' * pi/180;

for rf = 1:numel(requested_freqs_Hz)
    f_req = requested_freqs_Hz(rf);

    % --- CST: nearest available file ---
    [cst_err, cst_idx] = min(abs(cst_available_freqs_Hz - f_req));
    if cst_err > freq_tol_GHz*1e9
        fprintf('compare_cst_ads_farfields: WARNING - requested %.6g GHz, nearest available CST file is %.6g GHz (%.3g GHz away).\n', ...
            f_req/1e9, cst_available_freqs_Hz(cst_idx)/1e9, cst_err/1e9);
    end
    cst_f_actual = cst_available_freqs_Hz(cst_idx);
    cst_k = 2*pi*cst_f_actual/c0;

    cst_filename = fullfile(cst_farfield_directory, cst_farfield_files{cst_port_number}(cst_idx));
    [cst_rEpt, cst_theta, cst_phi] = read_cst_farfield(cst_filename, cst_voltage_unit);
    cst_rEpt = cst_rEpt * exp(1j*cst_k*1);   % remove e^(-jk) of the r = 1 m export
    cst_Eph = cst_rEpt(:,:,1);   % dim3: 1=rE_phi, 2=rE_theta (read_cst_farfield's convention)
    cst_Eth = cst_rEpt(:,:,2);

    % --- ADS: nearest available frequency slice ---
    [ads_err, ads_idx] = min(abs(ads_freq_Hz - f_req));
    if ads_err > freq_tol_GHz*1e9
        fprintf('compare_cst_ads_farfields: WARNING - requested %.6g GHz, nearest available ADS frequency is %.6g GHz (%.3g GHz away).\n', ...
            f_req/1e9, ads_freq_Hz(ads_idx)/1e9, ads_err/1e9);
    end
    ads_f_actual = ads_freq_Hz(ads_idx);
    ads_k = 2*pi*ads_f_actual/c0;

    ads_Vth = ads_Vtheta_all(:,:,ads_idx) * exp(1j*ads_k*1);   % remove e^(-jk) of the r = 1 m export
    ads_Vph = ads_Vphi_all(:,:,ads_idx)   * exp(1j*ads_k*1);

    title_str = sprintf('%.4g GHz', cst_f_actual/1e9);

    %% --- 3D Vector Plot Comparison (E-field) ---
    % Optional: uncomment to draw the 2x2 cone-plot comparison.
    % field(:,:,1/2/3) = Vr/Vth/Vph; Vr = 0 in the far field.
    % field_cst = zeros(numel(cst_theta), numel(cst_phi), 3);
    % field_cst(:,:,2) = cst_Eth;
    % field_cst(:,:,3) = cst_Eph;
    % 
    % field_ads = zeros(numel(ads_theta), numel(ads_phi), 3);
    % field_ads(:,:,2) = ads_Vth;
    % field_ads(:,:,3) = ads_Vph;
    % 
    % fignum = plot_cst_farfield_vector_compare(fignum, 1, ...
    %     cst_theta, cst_phi, field_cst, "CST", ads_theta, ads_phi, field_ads, "ADS", ...
    %     n_theta_ds, n_phi_ds, 1, '$\vec{E}$', title_str, FS, LW);

    %% --- Polar Cuts: |E_theta| and |E_phi| (2x2) ---
    % CST: black triangle markers about every cst_cut_ds_step_deg degrees;
    % ADS: continuous red line at its own resolution
    % (see plot_cst_farfield_cuts_compare).
    fignum = plot_cst_farfield_cuts_compare(fignum, ...
        cst_theta, cst_phi, cst_Eph, cst_Eth, ads_theta, ads_phi, ads_Vph, ads_Vth, ...
        phi_cut_rad, theta_cut_rad, cst_cut_ds_step_deg, plot_mag_scale, plot_unit_label, title_str, FS, legend_FS, LW, cut_tol_deg);
end