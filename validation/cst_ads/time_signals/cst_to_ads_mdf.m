%CST_TO_ADS_MDF  Convert CST port time signals into a Keysight ADS MDF excitation file.
% Imports the port time signals of a CST export in the antenna's Export
% folder and writes them into one ADS MDF file (output/cst_ads/mdf_out by
% default). Helpers: importCSTSignal.m, generate_mdf.m, cstVItoWaves.m
% (matlab/io); excitationCausalityWorkflow.m, plotCSTSignal.m (utils/).
%
% The MDF has one shared time column (seconds) and one voltage column per
% port, named by excitation type and port number:
%   PORT (generator)       : Vinc_<port number>  (i1 -> Vinc_1, ...)
%   VOLTAGE (ideal source) : Vs_<port number>    (i1 -> Vs_1, ...)
% In ADS, each source picks its port by using the matching name as its
% expression.
%
% The CST signal unit decides the conversion to volts:
%   W^(1/2) : generator open-circuit voltage with source impedance = port
%             Z0; V = value * sqrt(Z0), Z0 from the export's #Parameters
%             line (or Z0default). In ADS use a source with that impedance.
%   V       : ideal voltage source (zero source impedance); values used
%             as-is. In ADS use an ideal (zero-impedance) source.
%
% CST can write each port on its own time grid; the signals are resampled
% (linear interpolation) onto one common grid with the finest time step,
% spanning all signals, and set to zero where a signal has no data.
%
% Output name (the excitation is not tied to one simulation):
%   <CST export name>_<PORT|VOLTAGE>_<duration>_<f_lower>_<f_upper>.mdf
% e.g. excitation_gaussiansine_port_PORT_20ns_2.4GHz_2.6GHz.mdf
%   PORT    : CST signal in W^(1/2) (generator with source impedance Z0)
%   VOLTAGE : CST signal in V (ideal voltage source)
%   duration: signal length
%   f_lower, f_upper: edges of the relative bandwidth, where the spectrum
%   magnitude is bw_level_dB (20 dB) below its maximum, always in GHz
%   (up to 2 decimals); duration to 2 significant figures. A trailing frequency tag on the CST
%   export name (e.g. _2.4GHz_2.6GHz) is dropped.
%
% Port-voltage exports (CST 'Port voltages': "Port N [exc] [Real / V]") hold
% the total port voltage sqrt(Z0)*(a+b) for each excitation. By the
% substitution theorem an ideal voltage source imposing that voltage gives
% the same fields, so these are written as VOLTAGE excitations (Vs_<port>),
% one MDF per excitation, tagged _exc<excitation> (e.g. _exc1,
% _exccosim_Tran1_1). pv_ports = "driven" keeps only the driven port's
% voltage (the other ports stay terminated in Z0 in ADS, as in CST); "all"
% writes every port's voltage (all ports as ideal sources).
% If the matching port-current export exists ("Port N [exc] [Real / A]",
% file name with _current_ in place of _voltage_) and vi_mode = "incident",
% the Kurokawa power waves a = (V + Z I)/(2 sqrt(Z)), b = (V - Z I)/(2 sqrt(Z))
% (Z = each port's reference impedance: port_Z, else impedance_file, else Z0)
% are computed (cstVItoWaves.m), written to <mdf_directory>\from_VI\<name _iowaves_>.txt
% in CST layout, and the incident wave of the driven port becomes a PORT
% excitation (Vinc_<port>), one MDF per excitation.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(fullfile(here, 'utils')); % excitationCausalityWorkflow, plotCSTSignal

%% ===== USER SETTINGS =====
base_directory   = "C:\path\to\CST\";   % folder containing the CST project folders
mdf_directory    = fullfile(repo_root, "output", "cst_ads", "mdf_out");   % MDF output folder
% =========================

%% ---- project paths -------------------------------------------------------
antenna_name     = "rectangularpatch_symmetry_DAM";
sparam_directory = fullfile(base_directory, antenna_name, "SP");
signal_directory = fullfile(base_directory, antenna_name, "Export", "time_signals");

%% ---- user settings ----------------------------------------------------
signal_file    = "excitation_filteredsine_voltage_2.45GHz.txt";    % time-signal export in signal_directory
bw_level_dB    = 20;       % band edges: spectrum this many dB below its max
Z0default      = [];       % ohms; only used if the export has no Z0 parameter
blockName      = "BLK";    % name after BEGIN in the MDF file
prefix_port    = "Vinc_";  % column name for PORT (generator) signals = prefix + port number
prefix_voltage = "Vs_";    % column name for VOLTAGE (ideal source) signals
makePlot       = true;     % quick look at the imported signals
plot_band_dB   = 40;       % zoomed spectrum plot: band edges this many dB below the peak (shaded on
                           % the 0-to-Nyquist plot); [] = use bw_level_dB
plot_frange_GHz = [];      % zoomed spectrum x-limits [fmin fmax] in GHz, e.g. [1.5 3.5]; overrides plot_band_dB
check_causality = true;   % drive the ADS model's S-parameter block with each excitation; if it misses
                          % the causality target, show the check and write a filtered excitation
snp_check_file  = fullfile(repo_root, "output", "SnP", "S_dom_psi_circularloop_FDTD_CST.s8p");   % S-parameter model of the ADS circuit
causality_opts  = {'EncTolPct', 0.01, 'SigdB', 20, 'FirstTargetdB', -40, 'StepdB', 10, 'ShiftSigmas', 5};
                          % an excitation that meets the target only gets a one-line report; the full
                          % workflow (figures, filter, filtered CST/MDF files) runs only when it fails.
                          % For the S-parameter analysis and a filter on demand: filter_excitation_gaussian.m
% port-voltage exports only ("Port N [exc]")
vi_mode        = "incident";  % "incident": with the matching current export, power waves
                              %   a = (V + Z0 I)/(2 sqrt(Z0)) -> PORT MDF (Vinc_<port>)
                              % "voltage": the port voltage itself -> VOLTAGE MDF (Vs_<port>)
current_file   = "";          % current export; "" = signal_file with _voltage_ -> _current_
port_Z         = [];          % constant reference impedance per port (ohm), indexed by port number,
                              %   e.g. [50 75]; a scalar applies to all ports; [] = Z0 from the
                              %   export's #Parameters line for all ports
impedance_file = "";          % optional, used when port_Z is []: CST port-impedance export in
                              %   signal_directory, e.g. "cosim_port_impedances.txt"
pv_excitations = strings(0);  % excitation tags to write, e.g. ["1" "cosim Tran1 1"]; empty = all
pv_ports       = "driven";    % "driven": only the driven port's signal; "all": every port

if ~isfolder(mdf_directory), mkdir(mdf_directory); end

%% ---- import ------------------------------------------------------------
cstFile = fullfile(signal_directory, signal_file);
if ~isfile(cstFile)
    error('Time-signal file not found: %s', cstFile);
end

% each signal on its own grid
[tAll, vAll, infoAll] = importCSTSignal(cstFile, Z0default);

% port voltages + port currents -> power waves (written as a CST-layout file)
if any(strcmp(infoAll.sourceType, 'port voltage')) && vi_mode == "incident"
    if strlength(current_file) > 0
        curFile = fullfile(signal_directory, current_file);
    else
        curFile = fullfile(signal_directory, strrep(signal_file, "_voltage_", "_current_"));
    end
    if isfile(curFile) && ~strcmp(curFile, cstFile)
        [~, vn] = fileparts(cstFile);
        waveFile = fullfile(mdf_directory, "from_VI", strrep(vn, "_voltage_", "_iowaves_") + ".txt");
        if strcmp(waveFile, fullfile(mdf_directory, "from_VI", vn + ".txt"))
            waveFile = fullfile(mdf_directory, "from_VI", vn + "_iowaves.txt");
        end
        fprintf('Port voltages %s\n  + currents %s\n', cstFile, curFile);
        if ~isempty(port_Z)
            Zref = port_Z;                                       % user values
        elseif strlength(impedance_file) > 0
            Zref = fullfile(signal_directory, impedance_file);   % CST port-impedance export
            if ~isfile(Zref), error('Port-impedance file not found: %s', Zref); end
        else
            Zref = [];                                           % Z0 from #Parameters, all ports
        end
        cstVItoWaves(cstFile, curFile, waveFile, Zref);
        cstFile = waveFile;
        [tAll, vAll, infoAll] = importCSTSignal(cstFile, Z0default);
    else
        warning(['No current export %s: the port voltages are written as ideal voltage-source ' ...
            'excitations (vi_mode = "voltage").'], curFile);
    end
end
nAll = numel(infoAll.names);
[srcDir, srcName, srcExt] = fileparts(cstFile);

% excitations: incident waves i1, i2, ... (CST can also export output
% signals o1,1, ... in the same file), or, for a port-voltage export, the
% port voltages of each excitation
isInc = ~cellfun(@isempty, regexp(infoAll.names, '^i\d+( \[.*\])?$', 'once'));
isPV  = strcmp(infoAll.sourceType, 'port voltage');
groups = {};  gtags = strings(0);
if any(isInc) && all(cellfun(@isempty, infoAll.excitation(isInc)))
    groups = {find(isInc)};  gtags = "";                  % plain i1, i2, ...: one MDF
elseif any(isInc) || any(isPV)
    if any(isInc), isSel = isInc; else, isSel = isPV; end  % one MDF per excitation
    tags = unique(string(infoAll.excitation(isSel)), 'stable');
    if ~isempty(pv_excitations)
        miss = setdiff(string(pv_excitations), tags);
        if ~isempty(miss)
            error('Excitation(s) %s not in %s (available: %s)', strjoin(miss, ', '), cstFile, strjoin(tags, ', '));
        end
        tags = string(pv_excitations);
    end
    for tg = tags
        idx = find(isSel & strcmp(infoAll.excitation, tg));
        if pv_ports == "driven"
            tok = regexp(char(tg), '(\d+)\s*$', 'tokens', 'once');      % "1" -> 1, "cosim Tran1 2" -> 2
            if isempty(tok), drv = NaN; else, drv = str2double(tok{1}); end
            if isnan(drv) || ~any(infoAll.port(idx) == drv)
                error('Cannot tell the driven port of excitation "%s"; set pv_ports = "all".', tg);
            end
            idx = idx(infoAll.port(idx) == drv);
        end
        groups{end+1} = idx;                                                      %#ok<SAGROW>
        gtags(end+1)  = "_exc" + regexprep(regexprep(tg, '[^\w]+', '_'), '^_|_$', ''); %#ok<SAGROW>
    end
else
    error('No incident signals (i1, i2, ...) or port voltages (Port N [exc]) in %s', cstFile);
end
used = unique([groups{:}]);
if numel(used) < nAll
    fprintf('Skipping signals: %s\n', strjoin(infoAll.names(setdiff(1:nAll, used)), ', '));
end

fprintf('Antenna : %s\n', antenna_name);
fprintf('Source  : %s\n', cstFile);

S = struct('antenna_name', antenna_name, 'cstFile', cstFile, 'srcDir', srcDir, 'srcName', srcName, ...
    'srcExt', srcExt, 'bw_level_dB', bw_level_dB, 'mdf_directory', mdf_directory, 'Z0default', Z0default, ...
    'blockName', blockName, 'prefix_port', prefix_port, 'prefix_voltage', prefix_voltage, ...
    'makePlot', makePlot, 'plot_frange', plot_frange_GHz*1e9, 'plot_band_dB', plot_band_dB, 'check_causality', check_causality, 'snp_check_file', snp_check_file, ...
    'causality_opts', {causality_opts});
for g = 1:numel(groups)
    writeExcitationMDF(S, tAll, vAll, infoAll, groups{g}, gtags(g));
end

%% ========================================================================
function writeExcitationMDF(S, tAll, vAll, infoAll, sigIdx, gtag)
% one MDF from the signals sigIdx (one excitation)
nAll = numel(infoAll.names);
tRaw = tAll(sigIdx);  vRaw = vAll(sigIdx);
info = keepSignals(infoAll, sigIdx, nAll);
nSig = numel(sigIdx);
isPV = strcmp(info.sourceType, 'port voltage');
if gtag ~= "", fprintf('\nExcitation %s\n', extractAfter(gtag, "_exc")); end
fprintf('Signals : %d (%s)\n', nSig, strjoin(info.names, ', '));

% all signals on one shared grid (finest step, full span, zero outside a signal's range)
sameGrid = all(cellfun(@(x) isequal(size(x), size(tRaw{1})) && ...
    max(abs(x - tRaw{1})) <= 1e-9 * max(abs(tRaw{1})), tRaw));
if sameGrid
    t = tRaw{1};  V = [vRaw{:}];  resampled = false;
else
    dtc = min(cellfun(@(x) median(diff(x)), tRaw));
    t = (min(cellfun(@min, tRaw)) : dtc : max(cellfun(@max, tRaw))).';
    V = zeros(numel(t), nSig);
    for k = 1:nSig
        [tu_, iu] = unique(tRaw{k});
        V(:, k) = interp1(tu_, vRaw{k}(iu), t, 'linear', 0);
    end
    resampled = true;
end

%% ---- column names from port numbers ----------------------------------
varNames = strings(1, nSig);
for k = 1:nSig
    portNum = regexp(info.names{k}, '\d+', 'match', 'once');   % 'i1' -> '1', 'v2 [1]' -> '2'
    if isempty(portNum)
        portNum = matlab.lang.makeValidName(info.names{k});
        warning('Signal "%s" has no port number; its column uses the full name.', info.names{k});
    end
    if strcmp(info.sourceType{k}, 'generator')
        varNames(k) = S.prefix_port + portNum;                  % 'Vinc_1'
    else
        varNames(k) = S.prefix_voltage + portNum;               % 'Vs_1' (voltage source or port voltage)
    end
end
if numel(unique(varNames)) ~= nSig
    error('Two signals map to the same column name: %s', strjoin(varNames, ', '));
end

%% ---- comments ------------------------------------------------------------
comments = [
    "Antenna: "          + S.antenna_name
    "Source directory: " + S.srcDir
    "Source file: "      + S.srcName + S.srcExt
    "Time in seconds (CST export in " + string(info.timeUnit{1}) + ")"
    ];
for k = 1:nSig
    comments(end+1, 1) = varNames(k) + ": CST signal " + string(info.rawNames{k}) + ...
        " - " + signalDescription(info, k); %#ok<AGROW>
end
if resampled
    comments(end+1, 1) = sprintf(['Signals resampled onto a common grid: ' ...
        'dt = %.6g s, %d points; zero outside each signal''s range'], t(2) - t(1), numel(t));
end

%% ---- excitation type, duration and band -----------------------------------
if all(strcmp(info.sourceType, 'generator'))
    excType = "PORT";
elseif ~any(strcmp(info.sourceType, 'generator'))
    excType = "VOLTAGE";
else
    excType = "MIXED";
end
% band edges = frequencies where the spectrum magnitude is bw_level_dB
% below its maximum; with several signals, those of the signal with the
% largest spectral peak
Tdur = 0;  Xpeak = 0;  fpk = 0;  flo = 0;  fhi = 0;
for k = 1:nSig
    tk = tRaw{k}(:);  vk = vRaw{k}(:);
    Tdur = max(Tdur, tk(end) - tk(1));
    [fl, fu, fp, Xm] = spectrumBand(tk, vk, S.bw_level_dB);
    fprintf('  %s: peak %.4g GHz, -%g dB band %.4g - %.4g GHz\n', ...
        info.names{k}, fp/1e9, S.bw_level_dB, fl/1e9, fu/1e9);
    comments(end+1, 1) = sprintf('%s: spectral peak %.4g GHz; -%g dB band %.4g to %.4g GHz', ...
        varNames(k), fp/1e9, S.bw_level_dB, fl/1e9, fu/1e9); %#ok<AGROW>
    if Xm > Xpeak, Xpeak = Xm;  fpk = fp;  flo = fl;  fhi = fu; end
end
comments(end+1, 1) = sprintf(['Excitation type %s; duration %.4g ns; -%g dB bandwidth ' ...
    '%.4g to %.4g GHz (peak %.4g GHz)'], excType, Tdur*1e9, S.bw_level_dB, flo/1e9, fhi/1e9, fpk/1e9);

%% ---- causality check against the S-parameter block ---------------------------
if S.check_causality && isfile(S.snp_check_file)
    for kk = 1:nSig
        pn = info.port(kk);
        wfArgs = [S.causality_opts, {'Signal', info.names{kk}, 'OutDir', fullfile(S.mdf_directory, "filtered"), ...
            'PrefixPort', S.prefix_port, 'PrefixVoltage', S.prefix_voltage}];
        res = excitationCausalityWorkflow(S.cstFile, pn, S.snp_check_file, wfArgs{:}, ...
            'CheckOnly', true, 'Verbose', false);                 % quick pass/fail check only
        if res.passed
            it = find(strcmp(wfArgs(1:2:end), 'EncTolPct'), 1);
            if isempty(it), tol = 0.01; else, tol = wfArgs{2*it}; end   % workflow default 0.01 %
            fprintf('  %s meets the causality target: worst non-causal energy %.2g %% (target %g %%), no filter\n', ...
                info.names{kk}, res.encIn, tol);
            comments(end+1, 1) = sprintf('%s: causal with %s (worst non-causal energy %.2g %%)', ...
                varNames(kk), S.snp_check_file, res.encIn); %#ok<AGROW>
        else
            res = excitationCausalityWorkflow(S.cstFile, pn, S.snp_check_file, wfArgs{:});   % full workflow
            fprintf(2, ['  %s does NOT meet the causality target (%.3g %% non-causal energy): use the filtered excitation\n' ...
                '    CST: %s\n    ADS: %s\n'], info.names{kk}, res.encIn, res.cstFile, res.mdfFile);
        end
    end
elseif S.check_causality
    warning('S-parameter file for the causality check not found: %s', S.snp_check_file);
end

%% ---- write the MDF --------------------------------------------------------
if ~isfolder(S.mdf_directory), mkdir(S.mdf_directory); end
baseName = regexprep(S.srcName, '(_[\d.]+[kMGT]?Hz)+$', '');   % drop trailing frequency tag
outFile = fullfile(S.mdf_directory, baseName + gtag + "_" + excType + "_" + ...
    engStr(Tdur, "s") + "_" + ghzStr(flo) + "_" + ghzStr(fhi) + ".mdf");
generate_mdf(t, V, outFile, ...
    'BlockName', S.blockName, ...
    'Names',     cellstr(["time", varNames]), ...
    'Comments',  comments);

fprintf('Wrote %s\n', outFile);
fprintf('  %d points, %.4g to %.4g ns\n', numel(t), t(1)*1e9, t(end)*1e9);
for k = 1:nSig
    fprintf('  %-8s <- %-20s  %s, peak %.4g V\n', ...
        varNames(k), info.names{k}, signalDescription(info, k), max(abs(V(:, k))));
end
if any(isPV)
    fprintf(['  ADS: ideal (zero-impedance) voltage source(s) on port(s) %s; ' ...
        'terminate the other ports in Z0 = %g ohm\n'], strjoin(string(info.port), ', '), info.Z0(1));
end

%% ---- optional plot: time signal, spectrum magnitude and phase -----------
if S.makePlot
    plotBand = S.plot_band_dB;  if isempty(plotBand), plotBand = S.bw_level_dB; end
    for k = 1:nSig
        [~, ~, hp] = plotCSTSignal(S.cstFile, 'Signal', sigIdx(k), 'Z0', S.Z0default, ...
            'FRange', S.plot_frange, 'BandLevel', plotBand);
        hp.figure.Name = char(varNames(k) + " (" + info.names{k} + ")");
    end
end
end

function s = engStr(x, unit)
% value with SI prefix, 2 significant figures: 20.1e-9,'s' -> "20ns"
pre = ["f" "p" "n" "u" "m" "" "k" "M" "G" "T"];
if x == 0, s = "0" + unit; return, end
x = round(x, 2, 'significant');
e = floor(log10(abs(x)) / 3) * 3;
e = min(max(e, -15), 12);
s = string(num2str(x / 10^e)) + pre(e/3 + 6) + unit;
end

function s = ghzStr(f)
% frequency always in GHz, up to 2 decimals: 2.4e9 -> "2.4GHz", 8e8 -> "0.8GHz"
s = string(num2str(round(f/1e9, 2))) + "GHz";
end

function s = keepSignals(s, idx, nAll)
% keep entries idx of every per-signal field (length nAll) of info struct s
fn = fieldnames(s);
for j = 1:numel(fn)
    x = s.(fn{j});
    if ~ischar(x) && ~isstruct(x) && numel(x) == nAll && nAll > 1
        s.(fn{j}) = x(idx);
    end
end
end

function [flo, fhi, fpk, Xpk] = spectrumBand(t, v, levelDB)
% spectrum peak and the band edges where |X(f)| is levelDB below the peak
% (nonuniform DFT with trapezoid weights; coarse scan, then refined)
dt = diff(t);
w  = [dt(1); dt(1:end-1) + dt(2:end); dt(end)] / 2;
X  = @(f) abs(nufft(v .* w, t, f(:)));
fNyq = 1 / (2 * median(dt));
df   = 1 / (2 * (t(end) - t(1)));
fk = (0:df:fNyq).';  Xk = X(fk);
[~, ip] = max(Xk);
ff = linspace(max(fk(ip) - 2*df, 0), min(fk(ip) + 2*df, fNyq), 2001).';
[Xpk, jp] = max(X(ff));  fpk = ff(jp);
thr = Xpk * 10^(-levelDB/20);
il = find(Xk(1:ip) < thr, 1, 'last');
iu = ip - 1 + find(Xk(ip:end) < thr, 1, 'first');
if isempty(il), flo = 0;    else, flo = refine(fk(il), fk(il+1)); end
if isempty(iu), fhi = fNyq; else, fhi = refine(fk(iu-1), fk(iu)); end
    function fc = refine(fa, fb)
        g  = linspace(fa, fb, 2001).';
        Xg = X(g) - thr;
        j  = find(sign(Xg(1:end-1)) ~= sign(Xg(2:end)), 1, 'first');
        if isempty(j), fc = (fa + fb)/2; return, end
        fc = g(j) - Xg(j) * (g(j+1) - g(j)) / (Xg(j+1) - Xg(j));
    end
end

function s = signalDescription(info, k)
% one-line description of how signal k was converted to volts
switch info.sourceType{k}
    case 'generator'
        if strcmp(info.Z0source{k}, 'file'), zsrc = "from CST export"; else, zsrc = "Z0default setting"; end
        s = "generator open-circuit voltage, source impedance = port Z0 = " + ...
            sprintf('%g', info.Z0(k)) + " ohm (" + zsrc + "); V = CST [" + ...
            string(info.signalUnit{k}) + "] * sqrt(Z0)";
    case 'port voltage'
        s = "total port voltage sqrt(Z0)*(a+b) of a port excitation (Z0 = " + sprintf('%g', info.Z0(k)) + ...
            " ohm), imposed by an ideal (zero-impedance) voltage source; CST values in " + ...
            string(info.signalUnit{k}) + ", used as-is";
    otherwise
        s = "ideal voltage source (zero source impedance); CST values in " + ...
            string(info.signalUnit{k}) + ", used as-is";
end
end
