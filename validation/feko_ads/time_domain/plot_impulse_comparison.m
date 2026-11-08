%% plot_impulse_comparison.m
% Impulse response s21(t) of the Tx-loop / Rx-bowtie link, FEKO (.s2p block) vs ADS (modal SnP
% blocks and rational-fit modal blocks).
%
% Input:  ADS export Tx_loop_Rx_bowtie_ts.csv (csvFile below), read by import_impulse_response.m:
%         time [s], FEKO_S2P, ADS_SnPModes, ADS_RationalFitModes.
% Output: relative MSE of each ADS model vs FEKO, printed; figures:
%   impulse_response   FEKO and ADS SnP as thinned markers, ADS rational fit as a line
%   TimeDomain_legend  one-row legend shared with plot_RFpulse_comparison.m (timeDomainLegend.m)

clear; close all;
scriptFolder = fileparts(mfilename("fullpath"));
run(fullfile(scriptFolder, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(scriptFolder, fullfile(scriptFolder, '..', 'common'));

%% ===== USER SETTINGS =====
csvFile   = "C:\path\to\ADS\data\Tx_loop_Rx_bowtie_ts.csv";        % ADS impulse-response export
% =========================

import_impulse_response;

%% ---- SETTINGS ---------------------------------------------------
r  = 5;  a = 0.075;  c0 = 299792458;      % distance, antenna radius [m]
tr = (r - a)/c0;  tr_ns = tr*1e9;          % first-arrival time (anchors the plotted window)
tMark_ns = r/c0*1e9;  tMarkLabel = '$t = r/c_0$';   % vertical causality marker (alternative: tr_ns, '$t = (r-a)/c_0$')
obsTheta = 60;  obsPhi = 90;               % Rx direction [deg] (ADS schematic Txloop_Rxbowtie_internal_ports_transient)
tBefore_ns = 5;  tAfter_ns = 12;           % plotted window around tr
yLim = [];                                 % [] = symmetric, 10 % above the largest |value| in the window (rounded up to 0.01)
markerStep_s = 0.15e-9;                    % one marker per source every this time step (staggered)
namePrefix = "impulse_response";
S = plotStyle();

%% ---- PLOT -------------------------------------------------------
win = time >= tr - tBefore_ns*1e-9 & time <= tr + tAfter_ns*1e-9;
t   = time(win);  tns = t*1e9;
mF  = fekoMarkerIdx(t, markerStep_s, 1, 2);
mS  = fekoMarkerIdx(t, markerStep_s, 2, 2);

[fig, ax] = newFigure(namePrefix);
fe = FEKO_S2P(win);  sn = ADS_SnPModes(win);  rf = ADS_RationalFitModes(win);
h(1) = plot(ax, tns(mF), fe(mF), "LineStyle", "none", "Marker", "d", "Color", S.colFEKO);
h(2) = plot(ax, tns(mS), sn(mS), "LineStyle", "none", "Marker", "v", "Color", S.colors(4,:));   % Wong bluish green
h(3) = plot(ax, tns, rf, "LineStyle", "-", "Color", S.colADS);
xline(ax, tMark_ns, '-', tMarkLabel, 'Color', [0.3 0.5 0.5], 'LineWidth', 1.5, ...
    'Interpreter', 'latex', 'FontSize', S.fontSize, 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
styleAxes(ax, S, "Time (ns)", '$s_{21}(t)$');
if isempty(yLim), yMax = ceil(1.1*max(abs([fe; sn; rf]))*100)/100;  yLim = [-yMax yMax]; end
xlim(ax, tr_ns + [-tBefore_ns tAfter_ns]);  ylim(ax, yLim);
obsLabel(ax, S, obsTheta, obsPhi);      % observation point, top-right corner

%% ---- LEGEND (shared with the RF-pulse figure) ---------------------
timeDomainLegend(S);

%% ---- RELATIVE MSE compared to FEKO: sum|x - FEKO|^2 / sum|FEKO|^2 ----------------------------
relMSE = @(x, ref) sum((x(:) - ref(:)).^2) / sum(ref(:).^2);
fprintf('rel. MSE compared to FEKO (full record | plotted window):\n');
fprintf('  ADS SnP modes          : %.3f %% | %.3f %%\n', 100*relMSE(ADS_SnPModes, FEKO_S2P), 100*relMSE(sn, fe));
fprintf('  ADS rational-fit modes : %.3f %% | %.3f %%\n', 100*relMSE(ADS_RationalFitModes, FEKO_S2P), 100*relMSE(rf, fe));
