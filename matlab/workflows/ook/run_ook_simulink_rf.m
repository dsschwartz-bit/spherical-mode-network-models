%% RUN_OOK_SIMULINK_RF
% OOK direct antenna modulation with the RF Blockset version of the Tx-only model: RF S-parameters block
% (time domain), the internal (load) port terminated by an RF Variable Resistor that acts as the SPICE switch
% R = Ron*(Roff/Ron)^(1-s), s = clamp(Vs, 0, 1), driven by the same ideal square wave Vs (scale 1) as the
% CST / ADS OOK simulations; the feed is the generator Vinc (scale 2, series Rsrc_Tx = Z).
% Inputs are built as in ook_cst_to_ads_mdf.m (same CST exports, same time grid). The model and its
% _vars.mat are written to model_dir; every observation point is then simulated. This is the model used
% for the paper results (validation/simulink/ook/compare_ook_CST_ADS_simulink.m builds, simulates and
% compares in one go).

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
cst_root  = "C:\path\to\CST\";                                % folder containing the CST projects
snp_dir   = fullfile(repo_root, 'output', 'SnP');               % Touchstone files (rational fits in snp_dir/rationalfit)
model_dir = fullfile(repo_root, 'output', 'simulink', 'ook');   % model + <model>_vars.mat

% antenna: spherical-mode model with the feed (port 1) and an internal port for the switch
snp_file   = "S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p";   % in snp_dir
model_name = "Txrectangularpatch_internal_ports_OOK_RF";           % Simulink model to build

% excitation: feed incident voltage Vinc and switch control Vs
excitation_source = "user";          % "cst": CST external-port V/I exports; "user": a .mat from generate_user_excitations.m
cst_sig_dir  = fullfile(cst_root, "rectangularpatch_symmetry_DAM", "Export", "time_signals");   % CST time-signal exports
voltage_file = "excitation_filteredsine_voltage_2.45GHz.txt";   % "cst": port voltages
current_file = "excitation_filteredsine_current_2.45GHz.txt";   % "cst": port currents
rf_port = 1;  ctrl_port = 2;         % "cst": CST port numbers of the feed and the switch
user_exc_file = fullfile(repo_root, 'output', 'cst_ads', 'user_excitations', ...
    'excitation_sine2.45GHz_30ns_step66.7MHz_30ns_Vinc1_50Ohms_Vs2.mat');   % "user"

% observation points
obs_dirs = "cst";                   % "cst" = every far-field probe in ff_files, or [theta phi; ...] in deg
r_obs    = [];                      % observation distance (m); [] = the CST probe distance
ff_files = ["ff_Etheta_filteredsine_port1_2.45GHz.txt" "ff_Ephi_filteredsine_port1_2.45GHz.txt"];   % CST far-field probes in
                                    % cst_sig_dir; read for obs_dirs = "cst", r_obs = [] or compare_cst
compare_cst = true;                 % overlay and compare the CST far fields of ff_files (same excitation in CST); false: Simulink only
% =========================

%% ---- settings ---------------------------------------------------------------------
fast_restart = true;                % Simulink Fast Restart over the observation points (compile once; r, hence the delay, is fixed)
rf_fit      = "block";              % S-parameter model: "block" = the RF S-parameters block fits the Touchstone file itself (FitTol, MaxPoles,
                                    % FitOptions below); "file" = a stored rational fit, e.g. the one behind the SPICE netlist in snp_dir/rationalfit
rf_fit_file = fullfile(snp_dir, 'rationalfit', snp_file + "_Rfit.mat");   % .mat with Rfit (rationalfit,
                                    % N x N rfmodel.rational, any pole sharing) or a rational object; from generate_rationalfit_spice.m
rf_fit_z0   = [];                   % reference impedance (ohm) of that fit; [] = Z0 saved with it, else the Touchstone file's (eta0 = 376.73 here)
step_s  = 2e-12;                     % RF Blockset step (max step of the variable-step solver)
rf_fitopt = 'Share all poles';       % S-parameter fit: 'Share all poles' | 'Share poles by columns' | 'Fit individually'. If the result is inaccurate, try by columns or individually (slower), or tighten fit_tol_dB
% Simulink solver (RF platform), user choice: "ode45" = non-stiff (linear network) | "ode15s" = stiff (switches, time-varying R/L/C).
% The RF circuit itself is stepped at MaxStep (default Nyquist/6) whichever solver is chosen.
rf_solver = 'ode15s';               % ode15s here: the OOK model has a switch at the internal port
max_poles = 99;                     % clipped to 99 by the RF block
fit_tol_dB = -40;


mdl = char(model_name);
snp = char(fullfile(snp_dir, snp_file));
if ~isfile(snp), error('Touchstone file not found: %s', snp); end

%% ---- OOK inputs: CST exports (as ook_cst_to_ads_mdf.m) or the user-defined excitation --------
switch lower(string(excitation_source))       % only the selected source is read
    case "cst",  E = ook_load_excitation("cst", struct('sig_dir', cst_sig_dir, 'voltage_file', voltage_file, ...
                     'current_file', current_file, 'rf_port', rf_port, 'ctrl_port', ctrl_port));
    case "user", E = ook_load_excitation("user", struct('user_file', user_exc_file));
    otherwise,   error('excitation_source must be "cst" or "user".');
end
t = E.t;  Vinc = E.Vinc;  Vs = E.Vs;  Z = E.Z;  t0 = E.t0;
assignin('base', 'Vtx_wave', timeseries(Vinc, t - t0, 'Name', 'Vtx_wave'));
assignin('base', 'Vs_ctrl',  timeseries(Vs,   t - t0, 'Name', 'Vs_ctrl'));
assignin('base', 'Z0', Z);

%% ---- CST far field and observation points ------------------------------------------------
if ~isnumeric(obs_dirs) || isempty(r_obs) || compare_cst   % probe directions / distance / comparison
    cst = readCstFarfieldProbes(fullfile(cst_sig_dir, ff_files), ["Theta" "Phi"], 1, 1e-3);
end

if isempty(r_obs), r_probe = cst(1).r; else, r_probe = r_obs; end
if isnumeric(obs_dirs), dirs = obs_dirs; else, dirs = [[cst.theta].' [cst.phi].']; end
fprintf('%d observation point(s) at r = %g m\n', size(dirs, 1), r_probe);
th = dirs(1, 1);  ph = dirs(1, 2);

%% ---- S-parameter model: the block's own fit, or a stored rational fit (with its reference impedance) ------------
rfx = {};
if lower(string(rf_fit)) == "file"
    if ~isfile(rf_fit_file), error('Rational fit file not found: %s (save Rfit from generate_rationalfit_spice.m).', rf_fit_file); end
    F = load(rf_fit_file);  fn = fieldnames(F);
    kf = find(cellfun(@(q) isa(F.(q), 'rfmodel.rational') || isa(F.(q), 'rational'), fn), 1);
    if isempty(kf), error('%s holds no rational fit (rfmodel.rational array or rational object).', rf_fit_file); end
    assignin('base', 'ratTx_rf', F.(fn{kf}));
    if ~isempty(rf_fit_z0), z0fit = rf_fit_z0; elseif isfield(F, 'Z0'), z0fit = F.Z0; else, z0fit = sparameters(snp).Impedance; end
    rfx = {'RFRationalVar', 'ratTx_rf', 'RFRationalZ0', z0fit};
    [~, ffn, ffe] = fileparts(rf_fit_file);
    fprintf('S-parameter model: stored rational fit %s (variable %s), reference impedance %g ohm\n', [ffn ffe], fn{kf}, z0fit);
else
    fprintf('S-parameter model: fitted by the RF block (tol %g dB, max %d poles, %s)\n', fit_tol_dB, min(max_poles, 99), rf_fitopt);
end

%% ---- build (RF Blockset) and simulate -------------------------------------------------------
tBuild = tic;
build_tx_rx_link_simulink(snp, '', 'TxOnly', true, 'Platform', 'rf', 'LoadSwitch', true, 'ControlVar', 'Vs_ctrl', ...
    'Waveform', 'workspace', 'WaveformVar', 'Vtx_wave', 'SourceType', 'generator', 'Z0', Z, 'FeedZ', Z, ...
    'r', r_probe, 'theta_deg', th, 'phi_deg', ph, 'ModelName', mdl, 'MaxStep', step_s, ...
    'MaxPoles', max_poles, 'FitTolerance', fit_tol_dB, 'RFSolver', rf_solver, 'FitOptions', rf_fitopt,'OutputDir', model_dir, 'OpenModel', false, rfx{:});   % model + vars saved in model_dir
fprintf('Model build (incl. rational fit) %.1f s\n', toc(tBuild));
%% ---- simulate every observation point (one build; the direction only changes gains) --------
if compare_cst, cstCmp = cst; else, cstCmp = []; end
res = ook_observation_loop(mdl, dirs, r_probe, cstCmp, t0, 'RF Blockset', fast_restart);
