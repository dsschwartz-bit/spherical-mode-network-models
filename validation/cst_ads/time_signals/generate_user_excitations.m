%GENERATE_USER_EXCITATIONS  User-defined two-port excitations for CST, ADS and Simulink.
%
% Builds one CST ASCII signal per port (for import as a user-defined CST
% excitation) and ONE ADS MDF file with both ports in the layout of
% cst_to_ads_mdf.m: one block, shared time column, Vinc_1 (port 1,
% generator, source impedance Z0) and Vs_2 (port 2, ideal voltage source).
% The same two signals are saved for Simulink in a .mat file next to the
% MDF: timeseries Vtx_wave (= Vinc_1) and Vs_ctrl (= Vs_2), the variable
% names of the models' From Workspace blocks:
%   evalin('base', 'load(''<file>.mat'')')   before sim(...)
%
%   port 1 : sine at f_sine, starting at t = 0 (sin(0) = 0), duration T_sine;
%            generator signal in W^(1/2) (ADS Vinc_1 = value * sqrt(Z0), source
%            impedance Z0, the convention of cst_to_ads_mdf.m)
%   port 2 : controlling port, ideal voltage source: square wave of fundamental
%            f_square, 0 V for the first half of each period and step_level V
%            for the second, over T_sine; CST unit V, ADS Vs_2
%
% Each raw signal is written as a CST time-signal file and passed through
% excitationCausalityWorkflow.m with the options of
% filter_excitation_gaussian.m: the S-parameter block of the ADS model is
% driven by it and, if the responses are not causal enough, the signal is
% filtered by a delayed Gaussian; filtering stops when the next filter
% improves it by less than MinImprovePct. Both ports are always filtered.
% Port 1: if it already meets the target, with a Gaussian band-pass, -3 dB at
% the S-parameter band edges (lowest S frequency and fmax). Port 2 (the
% step): with the Gaussian low-pass, so its edges are smooth and its DC
% level is kept. Non-causal energy is counted above fmin_check only (the low
% end contributes little to radiation and is poorly known). The causality
% target is EncTolPct for port 1 and ctrl_tol_pct for port 2.
%
% Timing between the ports: each Gaussian delays its signal by its own shift.
% After the first pass the port with the smaller delay is processed again with
% the difference as an exact extra delay (ExtraDelay, linear phase, same
% filter), so the step stays aligned with the sine in the output files too.
%
% Inputs : a CST time-signal export of the target project (ref_export, for
%          the #Parameters line / Z0 and optionally the time step) and the
%          spherical-mode Touchstone model (snp_file, output/SnP).
% Outputs: output/cst_ads/user_excitations: raw\ CST signals, filtered CST
%          signals, the combined .mdf and the Simulink .mat.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(fullfile(here, 'utils')); % excitationCausalityWorkflow

%% ===== USER SETTINGS =====
cst_project = "C:\path\to\CST\rectangularpatch_symmetry_DAM\";    % target CST project folder (DAM variant, port Z0 = 50 ohm)
% =========================

%% ---- target model and files ------------------------------------------------
ref_export  = cst_project + "Export\time_signals\iowaves_filteredsine_port1_2.45GHz.txt";  % for #Parameters (Z0) (and the time step if dt = [])
snp_file    = fullfile(repo_root, "output", "SnP", "S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p");
out_dir     = fullfile(repo_root, "output", "cst_ads", "user_excitations");
raw_dir     = fullfile(out_dir, "raw");

%% ---- signals -------------------------------------------------------------------
f_sine     = 2.45e9;      % Hz, port 1 sine
T_sine     = 30e-9;       % s, duration of the sine and of the square wave
A_sine     = 1;           % CST amplitude of the sine, W^(1/2)
step_level = 1;           % V, port 2 (controlling port, ideal voltage source) high level
f_square   = 2/T_sine;    % Hz, square-wave fundamental (2/T_sine: two low/high periods over the sine;
                          %   1/T_sine: one step up at T_sine/2 and down at T_sine)
ctrl_tol_pct = 0.01;      % %, causality target (non-causal energy) for the controlling port 2 (same as port 1)
samples_per_cycle = 100;  % sample rate = samples_per_cycle * f_sine (CST/ADS interpolate in between)
dt         = 1/(samples_per_cycle*f_sine);   % s, sample step; [] = step of ref_export

fmin_check = 100e6;       % Hz, non-causal energy counted above this only (lowest frequency of the s14p data)
opts = { ...                  % as in filter_excitation_gaussian.m
    'EncTolPct',     0.01, ...
    'SigdB',         20, ...
    'FirstTargetdB', -40, ...
    'StepdB',        10, ...
    'ShiftSigmas',   5, ...
    'MaxIter',       6, ...
    'MinImprovePct', 5, ...       % stop when a filter improves the result by less than this
    'FminCheck',     fmin_check, ...
    'WriteIfCausal', true, ...    % write files even if no filter is needed
    'PlotAlways',    true, ...
    'OutDir',        out_dir, ...
    'WriteMDF',      false};      % no single-port MDFs: one MDF with both ports is written below

%% ---- reference export: parameters line, Z0, time step ------------------------
Lr      = readlines(ref_export);
parLine = Lr(find(startsWith(strtrim(Lr), "#Parameters"), 1));
[tr, ~, infoR] = importCSTSignal(ref_export);
Z0 = infoR.Z0(find(strcmp(infoR.names, "i1"), 1));
if isempty(dt), dt = median(diff(tr{1})); end
fprintf('Target %s: Z0 = %g ohm, dt = %.4g ps\n', cst_project, Z0, dt*1e12);

%% ---- raw signals (CST unit W^(1/2)) --------------------------------------------
t  = (0:dt:T_sine).';
x1 = A_sine * sin(2*pi*f_sine*t);
sig3 = @(x) string(num2str(x, 3));                 % numbers in file names: 3 significant figures
ghz = sig3(f_sine/1e9) + "GHz";
name1 = "excitation_sine" + ghz + "_" + round(T_sine*1e9) + "ns_port1";   % <f_sine>_<duration>
name2 = "excitation_step" + sig3(f_square/1e6) + "MHz_" + round(T_sine*1e9) + "ns_port2";   % <fundamental>_<duration>
hs = @(u) double(u > dt*1e-6) + 0.5*(abs(u) <= dt*1e-6);   % unit step, 1/2 at the edge sample (edge centred on it)
Tsq = 1/f_square;  kEdge = 0:floor(T_sine*f_square + 1e-9) - 1;      % periods in the record
stepSignal = @(tt, delay) step_level * sum(hs(tt - (kEdge + 0.5)*Tsq - delay) - hs(tt - (kEdge + 1)*Tsq - delay), 2);   % volts: up at (k+1/2)Tsq, down at (k+1)Tsq
x2 = stepSignal(t, 0);

if ~isfolder(raw_dir), mkdir(raw_dir); end
raw1 = writeRaw(fullfile(raw_dir, name1 + ".txt"), parLine, t, x1, 1, "W^(1/2)");
raw2 = writeRaw(fullfile(raw_dir, name2 + ".txt"), parLine, t, x2, 2, "V");
opts2 = opts;  opts2{find(strcmp(opts2(1:2:end), 'EncTolPct'))*2} = ctrl_tol_pct;   % relaxed target for the controlling port

%% ---- causality check / filter, then align the two filter delays --------------
fprintf('\n===== port 1: %s =====\n', name1);
o1 = excitationCausalityWorkflow(raw1, 1, snp_file, opts{:}, 'ForceFilter', true);   % port 1 is always Gaussian filtered
fprintf('\n===== port 2: %s =====\n', name2);
o2 = excitationCausalityWorkflow(raw2, 2, snp_file, opts2{:}, 'ForceFilter', true, 'ForceType', "lowpass");   % step always smoothed (Gaussian low-pass)

dD = o2.delay - o1.delay;                 % > 0: the step is delayed more than the sine
if abs(dD) > dt/2                         % re-run the less delayed port with the difference as an exact extra delay
    lag = 1 + (dD < 0);                   % port to delay
    fprintf('\nfilter delays: port 1 %.1f ps, port 2 %.1f ps -> port %d gets %.1f ps extra delay\n', o1.delay*1e12, o2.delay*1e12, lag, abs(dD)*1e12);
    close(findobj(0, 'Type', 'figure', '-regexp', 'Name', ternary(lag == 1, name1, name2)));
    if lag == 1
        delete(o1.cstFile);
        o1 = excitationCausalityWorkflow(raw1, 1, snp_file, opts{:}, 'ForceFilter', true, 'ExtraDelay', dD);
    else
        delete(o2.cstFile);
        o2 = excitationCausalityWorkflow(raw2, 2, snp_file, opts2{:}, 'ForceFilter', true, 'ForceType', "lowpass", 'ExtraDelay', -dD);
    end
end

%% ---- one ADS MDF with both ports (as cst_to_ads_mdf.m: one block, shared time column) -------
% common grid: the sample step of the signals, from 0 to the later end; zero where a signal has no data
dtc = min(median(diff(o1.t)), median(diff(o2.t)));
tc  = (0:dtc:max(o1.t(end), o2.t(end))).';
V1  = interp1(o1.t, o1.y, tc, 'linear', 0) * sqrt(Z0);   % Vinc_1 = CST [W^(1/2)] * sqrt(Z0)
V2  = interp1(o2.t, o2.y, tc, 'linear', 0);              % Vs_2 in V
[~, n1] = fileparts(o1.cstFile);  [~, n2] = fileparts(o2.cstFile);
mdfAll = fullfile(out_dir, strrep(regexprep(n1, '_port1.*$', '') + "_" + regexprep(extractAfter(n2, "excitation_"), '_port2.*$', '') + ...
    "_Vinc1_" + sig3(Z0) + "Ohms_Vs2", "-", "_") + ".mdf");   % ADS cannot import MDF names with hyphens
cmt = [ ...
    "User excitations for " + string(cst_project)
    "Vinc_1: port 1, generator open-circuit voltage source with source impedance Z0 = " + string(num2str(Z0)) + " ohm; Vinc_1 = CST [W^(1/2)] * sqrt(Z0); CST file " + string(n1) + ".txt"
    "Vs_2: port 2, ideal voltage source (V); CST file " + string(n2) + ".txt"
    sprintf('filter delays (incl. alignment): port 1 %.1f ps, port 2 %.1f ps', o1.delay*1e12, o2.delay*1e12)
    sprintf('common grid dt = %.6g s, %d points; zero outside each signal''s record', dtc, numel(tc))];
generate_mdf(tc, [V1 V2], mdfAll, 'BlockName', 'BLK', 'Names', {'time', 'Vinc_1', 'Vs_2'}, 'Comments', cmt);

% the same signals for Simulink: timeseries in a .mat file next to the MDF, named as the Simulink
% models' From Workspace blocks expect them (Vtx_wave: feed generator incident voltage = Vinc_1,
% Vs_ctrl: ideal-source control = Vs_2; interpolation is on in those blocks)
Vtx_wave = timeseries(V1, tc, 'Name', 'Vtx_wave');  Vtx_wave.DataInfo.Units = 'V';
Vs_ctrl  = timeseries(V2, tc, 'Name', 'Vs_ctrl');   Vs_ctrl.DataInfo.Units  = 'V';
exc_info = struct('cst_project', string(cst_project), 'Z0', Z0, 'cst_files', [string(o1.cstFile) string(o2.cstFile)], ...
    'mdf_file', string(mdfAll), 'delay_s', [o1.delay o2.delay], 'note', ...
    "Vtx_wave = Vinc_1 (generator incident voltage, source impedance Z0); Vs_ctrl = Vs_2 (ideal voltage source)");
matAll = replace(mdfAll, ".mdf", ".mat");
save(matAll, 'Vtx_wave', 'Vs_ctrl', 'exc_info');

%% ---- summary and overview plot ----------------------------------------------------
fprintf('\n%-6s %-12s %-12s %-14s %s\n', 'port', 'non-causal', 'filter (GHz)', 'delay (ps)', 'CST file');
o = {o1, o2};
for p = 1:2
        if ~isfinite(o{p}.f3) && isempty(o{p}.band), fd = "none"; elseif ~isempty(o{p}.band), fd = sprintf('BP %.3g-%.3g', o{p}.band/1e9); else, fd = sprintf('LP %.4g', o{p}.f3/1e9); end
    [~, fc, ec] = fileparts(o{p}.cstFile);
    fprintf('%-6d %8.4f %%   %-12s %-14.1f %s\n', p, o{p}.encOut, fd, o{p}.delay*1e12, fc + ec);
end
[~, fm, em] = fileparts(mdfAll);  [~, fa, ea] = fileparts(matAll);
fprintf('ADS MDF (Vinc_1, Vs_2): %s\nSimulink (Vtx_wave, Vs_ctrl): %s\nin %s\n', fm + em, fa + ea, out_dir);

set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
try, set(groot, 'defaultFigureWindowStyle', 'docked'); catch, end
fig = figure('Name', 'User excitations');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
[~, projName] = fileparts(strip(char(cst_project), 'right', '\'));       % project folder name only (a full path has backslashes, invalid in LaTeX)
title(tl, sprintf('User excitations for %s', strrep(projName, '_', '\_')), 'Interpreter', 'latex');
ax = nexttile(tl);  plot(ax, o1.t*1e9, o1.y*sqrt(Z0));  grid(ax, 'on');  box(ax, 'on');
ylabel(ax, '$V_{\mathrm{inc},1}$ (V)');  title(ax, sprintf('port 1: %g GHz sine, %g ns', f_sine/1e9, T_sine*1e9));
ax2 = nexttile(tl);  plot(ax2, o2.t*1e9, o2.y);  grid(ax2, 'on');  box(ax2, 'on');
ylabel(ax2, '$V_{\mathrm{s},2}$ (V)');  xlabel(ax2, 'Time (ns)');  title(ax2, sprintf('port 2: $0 \\rightarrow %g$ V step at the middle of the sine', step_level));
linkaxes([ax ax2], 'x');

%% ---- local functions --------------------------------------------------------------
function f = writeRaw(f, parLine, t, x, port, unit)
% CST time-signal layout (as excitationCausalityWorkflow writes it), time in ns
fid = fopen(f, 'w');
if ~isempty(parLine) && ~ismissing(parLine), fprintf(fid, '%s\n', parLine); end
fprintf(fid, '#"Time / ns"\t"i%d [Real / %s]"\n', port, unit);   % W^(1/2): generator, V: ideal voltage source
fprintf(fid, '#---------------------------------\n');
fprintf(fid, '%.13g\t%.13e\n', [t*1e9, x].');
fclose(fid);
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
