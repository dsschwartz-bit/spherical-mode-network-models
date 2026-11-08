%% plot_current_comparison_2.m
% Dipole load current under plane-wave excitation, Z-aligned vs rotated by 45 deg (both ADS).
%
% Input:  one ADS CSV export read by import_ADS_csv.m:
%         freq [Hz], mag_rotated [A], phase_rotated [deg], mag [A], phase [deg].
% Output: figures:
%   rotated_load_current_magnitude, rotated_load_current_phase,
%   rotated_load_current_phasor (polar plot of the phasor I over frequency),
%   rotated_load_current_legend

clear; close all;
scriptFolder = fileparts(mfilename("fullpath"));
run(fullfile(scriptFolder, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(scriptFolder, fullfile(scriptFolder, '..', 'common'));

%% ===== USER SETTINGS =====
adsCsvFile = "C:\path\to\ADS\data\PW_exciation_dipole_comparison.csv";   % ADS export (see import_ADS_csv.m)
% =========================

import_ADS_csv

freqScale  = 1e-9;  freqUnit = "GHz";
namePrefix = "rotated_load_current";
S = plotStyle();
clr = [S.colAlt; S.colADS];                       % Z-aligned blue, rotated red
labels = ["Z-aligned", "Rotated $45^{\circ}$"];
ls  = ["-", "--"];                              % Z-aligned solid, rotated dashed

%% ---- MAGNITUDE --------------------------------------------------
[fig, ax] = newFigure(namePrefix + "_magnitude");
h(1) = plot(ax, freq*freqScale, mag*1e6,         "Color", clr(1,:), "LineStyle", ls(1));
h(2) = plot(ax, freq*freqScale, mag_rotated*1e6, "Color", clr(2,:), "LineStyle", ls(2));
styleAxes(ax, S, "Frequency (" + freqUnit + ")", '$|I|$ ($\mu$A)');

%% ---- PHASE ------------------------------------------------------
[fig, ax] = newFigure(namePrefix + "_phase");
plot(ax, freq*freqScale, phase,         "Color", clr(1,:), "LineStyle", ls(1));
plot(ax, freq*freqScale, phase_rotated, "Color", clr(2,:), "LineStyle", ls(2));
styleAxes(ax, S, "Frequency (" + freqUnit + ")", 'Phase ($^\circ$)');
if all(abs([phase(:); phase_rotated(:)]) <= 180), ylim(ax, [-180 180]); yticks(ax, -180:90:180); end

%% ---- PHASOR (polar) --------------------------------------------
[fig, pax] = newFigure(namePrefix + "_phasor", "polar");
polarplot(pax, deg2rad(phase),         mag*1e6,         "Color", clr(1,:), "LineStyle", ls(1));
polarplot(pax, deg2rad(phase_rotated), mag_rotated*1e6, "Color", clr(2,:), "LineStyle", ls(2));
styleAxes(pax, S);
pax.RLim = [0, 1.1*max([mag(:); mag_rotated(:)])*1e6];
pax.ThetaTickLabel = compose('$%d^\\circ$', pax.ThetaTick);
polarLabel(fig, S, '$I$ ($\mu$A)');     % phasor plot: the complex current I itself

%% ---- LEGEND -----------------------------------------------------
legendFigure(namePrefix + "_legend", h, labels, 2);
