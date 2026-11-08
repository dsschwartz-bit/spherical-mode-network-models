%% RUN_TX_RADIATED_SIMULINK
% Radiated far field of the Simulink Tx-only spherical-mode model driven by a pulse, e.g. the
% Gaussian-modulated sine of a CST time-domain run. The model (build_tx_rx_link_simulink.m,
% 'TxOnly') is driven by the incident wave of a CST port-signal export, simulated at each
% observation direction, and r*E_theta(t), r*E_phi(t) are plotted.
% Results are kept in cache_file (simulinkRadiatedField.m): a rerun only simulates the
% directions not done yet, unless sim_mode = "rerun". The comparison with CST and ADS is
% validation/cst_ads/time_signals/compare_radiated_fields.m.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
cst_dir    = "C:\path\to\CST\circularloop_FDTD\Export\time_signals\";            % CST time-signal exports
tx_snp     = fullfile(repo_root, 'output', 'SnP', 'S_dom_psi_circularloop_FDTD_CST.s8p');   % spherical-mode Touchstone model
cst_io     = cst_dir + "iowaves_gaussiansine_0.2GHz_2GHz.txt";   % port-signal export: the Tx waveform
cst_signal = "i1";                    % its incident wave at the feed port
r          = [];                      % observation distance (m); [] = the probe distance of cst_ff
cst_ff     = cst_dir + "ff_Etheta_gaussiansine_0.2GHz_2GHz.txt";  % CST far-field export (only read when r = [])
directions = [0 45; 0 90; 90 90];     % [theta phi] in deg
Rload      = 1e9;                     % internal (load) port: 1e9 = open, 1e-6 = short; unused without one
cache_file = fullfile(repo_root, 'output', 'simulink', 'tx_radiated', 'circularloop_rE_gaussiansine_RF.mat');
% =========================

%% ---- Simulink model ------------------------------------------------------------------------
platform   = "rf";                    % "rf": RF Blockset (paper) | "simscape": Simscape Electrical
rf_fitopt  = "Share all poles";       % "Share all poles" | "Share poles by columns" | "Fit individually"
% Simulink solver (RF platform), user choice: "ode45" = non-stiff (linear network) | "ode15s" = stiff (switches, time-varying R/L/C).
% The RF circuit itself is stepped at MaxStep (default Nyquist/6) whichever solver is chosen.
rf_solver  = "ode45";
fit_tol_dB = -40;                     % rational-fit tolerance
sim_mode   = "auto";                  % "auto": cache + missing directions | "rerun": all again | "cache": cache only

%% ---- simulate ------------------------------------------------------------------------------
for f = [string(tx_snp), cst_io]
    if ~isfile(f), error('File not found: %s', f); end
end
if isempty(r)
    p = readCstFarfieldProbes(cst_ff, "Theta", 1, 1e-3);
    r = p(1).r;
end
if ~isfolder(fileparts(cache_file)), mkdir(fileparts(cache_file)); end
opt = struct('tx_snp', tx_snp, 'cst_io', cst_io, 'cst_signal', cst_signal, 'r', r, 'Rload', Rload, ...
    'fit_tol_dB', fit_tol_dB, 'platform', platform, 'rf_solver', rf_solver, 'rf_fitopt', rf_fitopt, ...
    'cache_file', cache_file, 'mode', sim_mode);
res = simulinkRadiatedField(directions, opt);

%% ---- plot ----------------------------------------------------------------------------------
comp = ["\theta" "\phi"];
for d = 1:numel(res)
    figure('Name', sprintf('Simulink r*E, theta %g, phi %g', res(d).theta, res(d).phi));
    for c = 1:2
        subplot(2, 1, c);
        plot(res(d).t * 1e9, res(d).e(:, c));  grid on
        ylabel(sprintf('r E_{%s} (V)', comp(c)));
        title(sprintf('\\theta = %g°, \\phi = %g°, r = %g m', res(d).theta, res(d).phi, r));
    end
    xlabel('t (ns)');
    fprintf('theta %5.1f  phi %5.1f :  peak |rE_theta| %.4g V,  peak |rE_phi| %.4g V\n', ...
        res(d).theta, res(d).phi, max(abs(res(d).e(:, 1))), max(abs(res(d).e(:, 2))));
end
