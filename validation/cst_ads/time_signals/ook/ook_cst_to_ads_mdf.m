%OOK_CST_TO_ADS_MDF  ADS MDF excitation for on-off keying (OOK) direct antenna modulation from CST port V/I.
% On-off-keying direct antenna modulation (DAM): build the ADS excitation MDF
% from the total voltages and currents at the CST schematic's EXTERNAL ports.
%
%   port 1 (RF)      : external port connected straight to the 3D model's port 1.
%                      Kurokawa power waves from the total V and I,
%                          a1 = (V1 + Z I1) / (2 sqrt(Z)),  b1 = (V1 - Z I1) / (2 sqrt(Z))
%                      (Z = port reference impedance, the export's Z0), and the
%                      incident wave becomes a generator excitation with source
%                      impedance Z: ADS Vinc_1 = a1 * sqrt(Z)  (= (V1 + Z I1)/2, the
%                      convention of cst_to_ads_mdf.m). Checked against i1 / o1 of the
%                      CST i/o-wave export (iowaves_file), which must agree.
%   port 2 (control) : external port drives internal schematic circuitry, not the
%                      3D model's port 2 directly, so there is no i2 / S22 wave at the
%                      3D port. Its total voltage is imposed by an ideal voltage source:
%                      ADS Vs_2 = V2 (CST values in V, used as-is). The port current I2
%                      is only reported (an ideal control input draws none).
%
% Both columns share one time grid (the finer of the two CST grids, linear
% interpolation, zero outside a signal's range) in one MDF, so the ports stay
% time-aligned in ADS. MDF names have no hyphens (ADS rejects them).
%
% Helpers: importCSTSignal.m, generate_mdf.m (matlab/io).
% Inputs : CST external-port voltage and current exports (and optionally the
%          i/o-wave export) in signal_directory.
% Output : one MDF in output/cst_ads/mdf_out (mdf_directory).

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
base_directory   = "C:\path\to\CST\";   % folder containing the CST project folders
mdf_directory    = fullfile(repo_root, "output", "cst_ads", "mdf_out");   % MDF output folder
% =========================

%% ---- project paths -------------------------------------------------------
antenna_name     = "rectangularpatch_symmetry_DAM";
signal_directory = fullfile(base_directory, antenna_name, "Export", "time_signals");

%% ---- user settings ----------------------------------------------------
voltage_file   = "excitation_filteredsine_voltage_2.45GHz.txt";  % external-port voltages "PortN [Real / V]"
current_file   = "";      % "" = voltage_file with _voltage_ -> _current_ ("PortN [Real / A]")
iowaves_file   = "iowaves_filteredsine_port1_2.45GHz.txt";      % i1/o1 of port 1 for the check; "" = skip
rf_port        = 1;       % generator port (power waves -> Vinc_<port>)
ctrl_port      = 2;       % controlling port (total voltage -> Vs_<port>)
port_Z         = [];      % RF port reference impedance (ohm); [] = Z0 from the export's #Parameters
blockName      = "BLK";
prefix_port    = "Vinc_";
prefix_voltage = "Vs_";
check_tol_pct  = 0.1;     % warn if a1/b1 differ from i1/o1 by more than this (relative RMS, %)
makePlot       = true;

%% ---- import ------------------------------------------------------------
vFile = fullfile(signal_directory, voltage_file);
if strlength(current_file) > 0, iFile = fullfile(signal_directory, current_file);
else, iFile = fullfile(signal_directory, strrep(voltage_file, "_voltage_", "_current_")); end
if ~isfile(vFile), error('Voltage export not found: %s', vFile); end
if ~isfile(iFile), error('Current export not found: %s', iFile); end

[tV, vV, infoV] = importCSTSignal(vFile, []);
[tI, vI, infoI] = importCSTSignal(iFile, []);
kV = @(p) find(infoV.port == p, 1);   kI = @(p) find(infoI.port == p, 1);
if isempty(kV(rf_port)) || isempty(kV(ctrl_port)) || isempty(kI(rf_port)) || isempty(kI(ctrl_port))
    error('Ports %d and %d must both be in %s and %s', rf_port, ctrl_port, vFile, iFile);
end
if isempty(port_Z), Z = infoV.Z0(kV(rf_port)); else, Z = port_Z; end
if isempty(Z) || isnan(Z), error('No Z0 in %s: set port_Z', vFile); end

t1 = tV{kV(rf_port)};    V1 = vV{kV(rf_port)};    I1 = interp1(tI{kI(rf_port)},   vI{kI(rf_port)},   t1, 'linear', 0);
t2 = tV{kV(ctrl_port)};  V2 = vV{kV(ctrl_port)};  I2 = interp1(tI{kI(ctrl_port)}, vI{kI(ctrl_port)}, t2, 'linear', 0);

fprintf('Antenna : %s\nVoltages: %s\nCurrents: %s\n', antenna_name, vFile, iFile);
fprintf('Port %d (RF)     : %d points, dt %.4g ps, %.4g-%.4g ns, Z = %g ohm\n', rf_port, numel(t1), median(diff(t1))*1e12, t1(1)*1e9, t1(end)*1e9, Z);
fprintf('Port %d (control): %d points, dt %.4g ps, %.4g-%.4g ns\n', ctrl_port, numel(t2), median(diff(t2))*1e12, t2(1)*1e9, t2(end)*1e9);

%% ---- port 1: power waves -------------------------------------------------
a1 = (V1 + Z*I1) / (2*sqrt(Z));      % W^(1/2)
b1 = (V1 - Z*I1) / (2*sqrt(Z));
relRMS = @(x, ref) 100 * norm(x - ref) / norm(ref);
chk = "";
if strlength(iowaves_file) > 0
    wFile = fullfile(signal_directory, iowaves_file);
    [tW, vW, infoW] = importCSTSignal(wFile, Z);
    ki = find(strcmp(infoW.names, "i" + rf_port), 1);  ko = find(strcmp(infoW.names, "o" + rf_port), 1);
    if isempty(ki), error('No i%d in %s', rf_port, wFile); end
    i1 = vW{ki} / sqrt(infoW.Z0(ki));                    % importCSTSignal scales W^(1/2) by sqrt(Z0)
    eA = relRMS(interp1(t1, a1, tW{ki}, 'linear', 0), i1);
    fprintf('Check vs %s:\n  a%d vs i%d: relative RMS difference %.3g %%\n', iowaves_file, rf_port, rf_port, eA);
    chk = string(sprintf('a%d vs CST i%d (%s): relative RMS difference %.3g %%', rf_port, rf_port, iowaves_file, eA));
    if ~isempty(ko)
        o1 = vW{ko} / sqrt(infoW.Z0(ko));
        eB = relRMS(interp1(t1, b1, tW{ko}, 'linear', 0), o1);
        fprintf('  b%d vs o%d: relative RMS difference %.3g %%\n', rf_port, rf_port, eB);
        chk = chk + sprintf('; b%d vs o%d: %.3g %%', rf_port, rf_port, eB);
    end
    if eA > check_tol_pct
        warning('Port %d incident wave from V/I differs from the CST i/o-wave export by %.3g %% (> %g %%): check port_Z and the exports.', ...
            rf_port, eA, check_tol_pct);
    end
end

%% ---- port 2: total voltage ---------------------------------------------
fprintf('Port %d: V %.4g to %.4g V, max |I| %.3g A (ideal voltage source, current not used)\n', ...
    ctrl_port, min(V2), max(V2), max(abs(I2)));

%% ---- common grid and MDF ------------------------------------------------
dtc = min(median(diff(t1)), median(diff(t2)));
t   = (min(t1(1), t2(1)) : dtc : max(t1(end), t2(end))).';
Vinc = interp1(t1, a1*sqrt(Z), t, 'linear', 0);       % ADS generator: Vinc = a * sqrt(Z)
Vs   = interp1(t2, V2,         t, 'linear', 0);
names = [prefix_port + rf_port, prefix_voltage + ctrl_port];

[~, vn] = fileparts(vFile);
baseName = regexprep(strrep(vn, "_voltage", ""), '(_[\d.]+[kMGT]?Hz)+$', '');
Tdur = t(end) - t(1);
outFile = fullfile(mdf_directory, strrep(baseName + "_OOK_" + names(1) + "_" + names(2) + "_" + ...
    round(Tdur*1e9, 2, 'significant') + "ns", "-", "_") + ".mdf");

comments = [
    "Antenna: " + antenna_name
    "OOK direct antenna modulation: external-port total voltages and currents"
    "Voltages: " + vFile
    "Currents: " + iFile
    "Time in seconds"
    names(1) + ": generator incident wave at port " + rf_port + ", a = (V + Z I)/(2 sqrt(Z)), Vinc = a sqrt(Z), source impedance Z = " + Z + " ohm"
    names(2) + ": total voltage at external port " + ctrl_port + ", ideal (zero-impedance) voltage source, CST values in V used as-is"
    sprintf('Common grid: dt = %.6g s, %d points; zero outside each signal''s range', dtc, numel(t))
    ];
if chk ~= "", comments(end+1) = chk; end

if ~isfolder(mdf_directory), mkdir(mdf_directory); end
generate_mdf(t, [Vinc Vs], outFile, 'BlockName', blockName, 'Names', cellstr(["time" names]), 'Comments', comments);
fprintf('Wrote %s\n  %d points, %.4g to %.4g ns; %s peak %.4g V, %s %.4g to %.4g V\n', outFile, numel(t), ...
    t(1)*1e9, t(end)*1e9, names(1), max(abs(Vinc)), names(2), min(Vs), max(Vs));
fprintf(['  ADS: port %d generator (source impedance %g ohm) with Vinc = %s; port %d ideal voltage source with V = %s,\n' ...
    '       driving the same internal control circuitry as in the CST schematic\n'], rf_port, Z, names(1), ctrl_port, names(2));

%% ---- plots (docked tabs) -------------------------------------------------
if makePlot
    set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
    try, set(groot, 'defaultFigureWindowStyle', 'docked'); catch, end
    co = orderedcolors("gem");
    figure('Name', "Port " + rf_port + " waves");
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile; hold on; box on; grid on;
    plot(t1*1e9, a1, 'Color', co(1,:), 'LineWidth', 1.4, 'DisplayName', sprintf('$a_{%d}$ from $V$, $I$', rf_port));
    if exist('i1', 'var'), plot(tW{ki}*1e9, i1, 'k--', 'LineWidth', 1, 'DisplayName', sprintf('CST $i_{%d}$', rf_port)); end
    ylabel('$\mathrm{W}^{1/2}$'); title(sprintf('Incident wave, port %d', rf_port)); legend('Location', 'best');
    nexttile; hold on; box on; grid on;
    plot(t1*1e9, b1, 'Color', co(2,:), 'LineWidth', 1.4, 'DisplayName', sprintf('$b_{%d}$ from $V$, $I$', rf_port));
    if exist('o1', 'var'), plot(tW{ko}*1e9, o1, 'k--', 'LineWidth', 1, 'DisplayName', sprintf('CST $o_{%d}$', rf_port)); end
    xlabel('Time (ns)'); ylabel('$\mathrm{W}^{1/2}$'); title(sprintf('Reflected wave, port %d', rf_port)); legend('Location', 'best');

    figure('Name', "Port " + ctrl_port + " control");
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile; plot(t2*1e9, V2, 'Color', co(1,:), 'LineWidth', 1.4); grid on; box on;
    ylabel(sprintf('$V_{%d}$ (V)', ctrl_port)); title(sprintf('Total voltage, external port %d', ctrl_port));
    nexttile; plot(t2*1e9, I2*1e3, 'Color', co(2,:), 'LineWidth', 1.4); grid on; box on;
    xlabel('Time (ns)'); ylabel(sprintf('$I_{%d}$ (mA)', ctrl_port)); title(sprintf('Current, external port %d', ctrl_port));

    figure('Name', 'MDF columns');
    yyaxis left;  plot(t*1e9, Vinc, 'LineWidth', 1.2); ylabel(sprintf('$V_{\\mathrm{inc},%d}$ (V)', rf_port));
    yyaxis right; plot(t*1e9, Vs,   'LineWidth', 1.4); ylabel(sprintf('$V_{\\mathrm{s},%d}$ (V)', ctrl_port));
    grid on; box on; xlabel('Time (ns)'); title('ADS excitation (common time grid)');
end
