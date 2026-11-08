%% plot_RFpulse_comparison.m
% RF-pulse response of the Tx-loop / Rx-bowtie link, FEKO (.s2p block) vs ADS (rational-fit
% modes), with the generator voltage on the left axis.
%
% Input:  ADS transient export Tx_loop_Rx_bowtie_RFpulse.csv (csvFile below), read by
%         import_RFpulse_response.m: time [s], Vgen [V], FEKO_S2P [V], ADS_RationalFitModes [V].
% Output: relative MSE of ADS vs FEKO, printed; figures:
%   RFpulse_response   left axis V_gen (mV, thick semi-transparent line); right axis V_rx (uV),
%                      FEKO thinned markers, ADS line
%   TimeDomain_legend  one-row legend shared with plot_impulse_comparison.m (timeDomainLegend.m)

clear; close all;
scriptFolder = fileparts(mfilename("fullpath"));
run(fullfile(scriptFolder, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(scriptFolder, fullfile(scriptFolder, '..', 'common'));

%% ===== USER SETTINGS =====
csvFile   = "C:\path\to\ADS\data\Tx_loop_Rx_bowtie_RFpulse.csv";   % ADS transient export
% =========================

import_RFpulse_response;

%% ---- SETTINGS ---------------------------------------------------
r  = 5;  a = 0.075;  c0 = 299792458;      % distance, antenna radius [m]
tr = (r - a)/c0;  tr_ns = tr*1e9;          % first arrival from the Tx antenna surface (sets the plotted window)
tMark_ns = r/c0*1e9;  tMarkLabel = '$t = r/c_0$';   % vertical causality marker (alternative: tr_ns, '$t = (r-a)/c_0$')
obsTheta = 60;  obsPhi = 90;               % Rx direction [deg] (ADS schematic Txloop_Rxbowtie_internal_ports_transient)
tAfter_ns = 12;                            % plotted window: 0 ... tr + tAfter
markerStep_s = 0.25e-9;                    % FEKO marker every this time step ([] = every sample)
namePrefix = "RFpulse_response";
S = plotStyle();
colGen = S.colExc;                         % generator voltage (left axis): thick, 40 % transparent blue

%% ---- PLOT -------------------------------------------------------
win = time >= 0 & time <= tr + tAfter_ns*1e-9;
t   = time(win);  tns = t*1e9;
fe  = FEKO_S2P(win)*1e6;  ad = ADS_RationalFitModes(win)*1e6;  vg = Vgen(win)*1e3;
mF  = fekoMarkerIdx(t, markerStep_s);

[fig, ax] = newFigure(namePrefix);
yyaxis(ax, 'left');
h(1) = plot(ax, tns, vg, "LineStyle", "-", "Color", [colGen S.alphaExc], "LineWidth", S.lwExc, "Tag", "keepWidth");
ylabel(ax, '$V_{\mathrm{gen}}$ (mV)', 'Interpreter', 'latex', 'FontSize', S.labelSize);
ax.YColor = colGen;
yyaxis(ax, 'right');  hold(ax, 'on');
h(2) = plot(ax, tns(mF), fe(mF), "LineStyle", "none", "Marker", "d", "Color", S.colFEKO);
h(3) = plot(ax, tns, ad, "LineStyle", "-", "Color", S.colADS);
ylabel(ax, '$V_{\mathrm{rx}}$ ($\mu$V)', 'Interpreter', 'latex', 'FontSize', S.labelSize);
ax.YColor = [0 0 0];
xline(ax, tMark_ns, '-', tMarkLabel, 'Color', [0.3 0.5 0.5], 'LineWidth', 1.5, ...
    'Interpreter', 'latex', 'FontSize', S.fontSize, 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
styleAxes(ax, S, "Time (ns)");
styleTraces(h, S);                           % both y sides (findobj only sees the active one)
xlim(ax, [0, tr_ns + tAfter_ns]);
obsLabel(ax, S, obsTheta, obsPhi);      % observation point, top-right corner

%% ---- LEGEND (shared with the impulse-response figure) -------------
timeDomainLegend(S);

%% ---- RELATIVE MSE compared to FEKO: sum|x - FEKO|^2 / sum|FEKO|^2 ----------------------------
relMSE = @(x, ref) sum((x(:) - ref(:)).^2) / sum(ref(:).^2);
fprintf('rel. MSE compared to FEKO (full record | plotted window):\n');
fprintf('  ADS rational-fit modes : %.3f %% | %.3f %%\n', 100*relMSE(ADS_RationalFitModes, FEKO_S2P), 100*relMSE(ad, fe));
