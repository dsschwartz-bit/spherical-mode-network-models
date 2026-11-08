%% RUN_TX_RX_LINK_SIMULINK
% Build the Simulink Tx-to-Rx spherical-mode link (build_tx_rx_link_simulink.m), simulate it
% (optionally for several distances r) and plot Vtheta, Vphi and the Rx feed voltage.
% The link is built on the RF Blockset platform (builder default); pass 'Platform','simscape' for the
% Simscape Electrical version.
%
% Inputs : Touchstone file(s) of the Tx / Rx spherical-mode networks (output/SnP); for
%          waveform_mode "cst" a CST ASCII time-signal export (cst_file).
% Outputs: the model and <model>_vars.mat in model_dir, figures in plot_dir, and the
%          struct array "results" in the workspace.
%
% Tx waveform: waveform_mode = "sine" | "cst" | "workspace"
%   "cst"       : CST ASCII time-signal export (importCSTSignal, as in cst_to_ads_mdf.m)
%   "workspace" : the signal defined in the USER SIGNAL section below (any [t v] data)
% After the first build set rebuild = false to reuse the model: only the workspace
% variables are updated (r, theta, phi, rotation, signal), link_update_gains refreshes the
% derived ADS VARs (prop_t, gradYlm_*, ...), and the stop time follows prop_t.
%
% EXAMPLES (copy the lines between %{ and %} into the MATLAB console):
%{
% 1) run this script with the settings below (builds, simulates, plots)
run_tx_rx_link_simulink

% 2) build only: sine source, r = 5 m, Rx 30 deg off boresight, then simulate
snp = fullfile(repo_root, 'output', 'SnP', 'S_dom_Q_dipole_cylinder_edgeport_full.s6p');
mdl = build_tx_rx_link_simulink(snp, snp, 'r', 5, 'theta_deg', 60, 'Waveform', 'sine', 'SourceFreq', 1e9);
out = sim(mdl);
v = out.yout.getElement('Vtheta').Values;  plot(v.Time*1e9, v.Data), xlabel('ns')

% 3) Tx waveform from a CST export (same file as cst_to_ads_mdf.m reads)
mdl = build_tx_rx_link_simulink(snp, snp, 'Waveform', 'cst', ...
    'CSTFile', 'C:\path\to\CST\<antenna>\Export\time_signals\port2_OC\iowaves_gaussian_2GHz_3GHz.txt', ...
    'CSTSignal', 'i1');

% 4) user-defined Tx waveform: a function handle, [t v] matrix, timeseries or struct(t,v)
mywave = @(t) exp(-((t-1.5e-9)/0.4e-9).^2) .* sin(2*pi*1e9*t);
mdl = build_tx_rx_link_simulink(snp, snp, 'Waveform', 'workspace', 'WaveformVar', 'mywave', ...
    'WaveformDuration', 3e-9, 'SourceType', 'generator');      % 'voltage' = ideal source
t = (0:20e-12:3e-9).';  wave = [t, sin(2*pi*1e9*t).*hann(numel(t))];   % [t(s) v(V)] matrix
mdl = build_tx_rx_link_simulink(snp, snp, 'Waveform', 'workspace', 'WaveformVar', 'wave');

% 5) reuse the built model: change workspace variables, refresh the derived ones, simulate
r = 20;  theta = 60*pi/180;  rot_angle = 15*pi/180;           % ADS VAR values (SI, radians)
link_update_gains                                              % prop_t, Lind, gradYlm_*, Gtx_prime, ...
t_exc = 5e-9;  set_param(mdl, 'StopTime', 'prop_t + t_exc')   % excitation duration + propagation delay
out = sim(mdl);

% 6) read results: Dataset elements are named like the nets
out.yout.getElement('FeedRx').Values                           % differential port voltage, timeseries
out.yout.getElement('Vphi').Values.Data
out.yout.getElementNames                                       % everything that was logged
%}

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
cst_file  = "C:\path\to\CST\<antenna>\Export\time_signals\port2_OC\iowaves_gaussian_2GHz_3GHz.txt";   % CST export (waveform_mode "cst")
snp_dir   = fullfile(repo_root, 'output', 'SnP');              % folder with the Touchstone files
model_dir = fullfile(repo_root, 'output', 'simulink');         % model + <model>_vars.mat
plot_dir  = fullfile(repo_root, 'output', 'simulink', 'plots');
% =========================

%% ---- Touchstone files (Tx and Rx antenna) -------------------------------
tx_snp  = fullfile(snp_dir, 'S_dom_Q_dipole_cylinder_edgeport_full.s6p');
rx_snp  = tx_snp;

%% ---- link geometry (the ADS VAR values) -----------------------------------
r_values      = 10;         % m; a vector runs one simulation per distance
theta_deg     = 90;         % Rx direction seen from Tx
phi_deg       = 0;
rot_angle_deg = 0;          % Rx rotation about axis_u
axis_u        = [1 1 1];    % rotation axis (Ux Uy Uz)

%% ---- Tx input signal -------------------------------------------------------
waveform_mode = "workspace";                 % "sine" | "cst" | "workspace"
source_type   = "generator";                 % "generator": V = incident-wave voltage V+ (open circuit 2*V+), series feed Z
                                             % "voltage"  : ideal source ("cst" derives the type from the file)
% sine
sine_freq_Hz  = 1e9;   sine_amp_V = 1;
% CST export (cst_file: USER SETTINGS)
cst_signal    = "i1";                        % signal name or index
% USER SIGNAL (waveform_mode = "workspace"): anything that ends up in Vtx_signal = [t(s) v(V)]
dt  = 20e-12;
t   = (0:dt:3e-9).';
Vtx_signal = [t, exp(-((t - 1.5e-9)/0.4e-9).^2) .* sin(2*pi*1e9*t)];

Z0_ohm = [];     % Feed source / termination and internal-port load resistance (the variable Z0);
                 % [] = the CST export's Z0 for waveform_mode "cst", else the workspace variable Z0
                 % (if it exists), else the Touchstone header

%% ---- run options ---------------------------------------------------------------
rebuild      = true;      % false: reuse the saved model, only update workspace variables
model_name   = '';        % '' = default name from the antenna names
fit_tol_dB   = -40;       % rational-fit tolerance of the S-parameter blocks
% Simulink solver (RF platform), user choice: "ode45" = non-stiff (linear network) | "ode15s" = stiff (switches, time-varying R/L/C).
% The RF circuit itself is stepped at MaxStep (default Nyquist/6) whichever solver is chosen.
rf_solver    = "ode45";
font_size    = 14;
plot_results = true;
save_plots   = true;      % write the figure as PNG (and .fig) into plot_dir

%% ---- build or update ------------------------------------------------------------
args = {'r', r_values(1), 'theta_deg', theta_deg, 'phi_deg', phi_deg, 'rot_angle_deg', rot_angle_deg, ...
    'axis', axis_u, 'FitTolerance', fit_tol_dB, 'RFSolver', char(rf_solver), 'FontSize', font_size, 'OpenModel', false, 'OutputDir', model_dir};
if ~isempty(Z0_ohm), args = [args, {'Z0', Z0_ohm}]; end
if ~isempty(model_name), args = [args, {'ModelName', model_name}]; end
switch lower(waveform_mode)
    case "sine"
        args = [args, {'Waveform', 'sine', 'SourceFreq', sine_freq_Hz, 'SourceAmplitude', sine_amp_V}];
    case "cst"
        args = [args, {'Waveform', 'cst', 'CSTFile', char(cst_file), 'CSTSignal', cst_signal}];
    case "workspace"
        Vtx_wave = Vtx_signal;                                   % N x 2 [t v] read by the builder
        args = [args, {'Waveform', 'workspace', 'WaveformVar', 'Vtx_wave', 'SourceType', char(source_type)}];
    otherwise
        error('waveform_mode must be "sine", "cst" or "workspace".');
end

if rebuild
    mdl = build_tx_rx_link_simulink(tx_snp, rx_snp, args{:});
else
    if isempty(model_name), error('rebuild = false needs model_name (the saved model to reuse).'); end
    mdl = model_name;
    if ~bdIsLoaded(mdl), load_system(fullfile(model_dir, mdl)); end
end

%% ---- simulate (one run per distance) ---------------------------------------------
results = struct('r', {}, 'prop_t', {}, 'out', {});
for i = 1:numel(r_values)
    % workspace variables the blocks read (ADS VAR block)
    assignin('base', 'r', r_values(i));
    assignin('base', 'theta', theta_deg*pi/180);  assignin('base', 'phi', phi_deg*pi/180);
    assignin('base', 'rot_angle', rot_angle_deg*pi/180);
    assignin('base', 'Ux', axis_u(1));  assignin('base', 'Uy', axis_u(2));  assignin('base', 'Uz', axis_u(3));
    link_update_gains;                                                % prop_t, Lind, gradYlm_*, Gtx_prime, ...
    pt = evalin('base', 'prop_t');
    if i > 1 || ~rebuild
        if string(waveform_mode) == "workspace"
            ts = timeseries(Vtx_signal(:, 2), Vtx_signal(:, 1), 'Name', 'Vtx_wave');
            assignin('base', 'Vtx_wave', ts);
        end
    end
    wave = evalin('base', 'Vtx_wave');
    if isa(wave, 'timeseries'), tend = wave.Time(end); else, tend = wave(end, 1); end
    assignin('base', 't_exc', tend);  set_param(mdl, 'StopTime', 'prop_t + t_exc');
    fprintf('r = %g m: prop_t = %.4g ns, simulating to %.4g ns...\n', r_values(i), pt*1e9, (pt + tend)*1e9);
    out = sim(mdl);
    results(i).r = r_values(i);  results(i).prop_t = pt;  results(i).out = out;
end

%% ---- plots ---------------------------------------------------------------------------
if plot_results
    sig = @(out, name) out.yout.getElement(name).Values;
    set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
    figure('WindowStyle', 'docked', 'Name', 'Tx-Rx link: far-field voltages');
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact');
    names = {'Vtheta', 'Vphi', 'FeedRx'};
    for k = 1:3
        nexttile; hold on
        for i = 1:numel(results)
            ts = sig(results(i).out, names{k});
            plot(ts.Time*1e9, ts.Data, 'DisplayName', sprintf('r = %g m', results(i).r));
        end
        ylabel([names{k} ' (V)']); grid on
        if k == 1, legend('Location', 'best'); end
    end
    xlabel(tl, 'Time (ns)', 'Interpreter', 'latex');
    title(tl, sprintf('%s  ($\\theta=%g^\\circ$, $\\phi=%g^\\circ$, rot $=%g^\\circ$)', strrep(mdl, '_', '\_'), theta_deg, phi_deg, rot_angle_deg), 'Interpreter', 'latex');
    fh = gcf;
    if save_plots
        if ~isfolder(plot_dir), mkdir(plot_dir); end
        exportgraphics(fh, fullfile(plot_dir, mdl + "_Vtheta_Vphi_FeedRx.png"), 'Resolution', 120);
        savefig(fh, fullfile(plot_dir, mdl + "_Vtheta_Vphi_FeedRx.fig"));
        fprintf('Plots saved in %s\n', plot_dir);
    end
end
fprintf('Results in the struct array "results"; signals: results(i).out.yout (Dataset), out.tout.\n');
