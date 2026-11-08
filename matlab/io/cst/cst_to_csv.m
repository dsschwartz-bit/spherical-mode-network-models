function outFile = cst_to_csv(cstFile, outFile, varargin)
%CST_TO_CSV  Convert one CST ASCII time signal to a two-column CSV for FEKO.
%
%   cst_to_csv(cstFile)
%   cst_to_csv(cstFile, outFile)
%   cst_to_csv(cstFile, outFile, 'Signal', 2, 'Z0', 50, 'Header', true)
%   outFile = cst_to_csv(...)
%
%   Writes a CSV with column 1 = time in ns and column 2 = voltage in
%   volts, for use as a FEKO time signal. Uses importCSTSignal.m to read
%   the file; time is converted to ns whatever unit CST used. The CST
%   signal unit decides the voltage conversion:
%     W^(1/2) : generator open-circuit voltage (source impedance = port
%               Z0); V = value * sqrt(Z0), Z0 from the CST header
%     V       : ideal voltage source; values used as-is, Z0 not needed
%
%   Inputs
%     cstFile : CST ASCII export (e.g. excitation_port1.txt)
%     outFile : (optional) output CSV name or path. A bare file name
%               (no folder) is written to <repo>/output/cst_ads/csv_out,
%               which is created if needed; a name with a folder is
%               written exactly there. Default: <cst name>.csv, or
%               <cst name>_<signal>.csv when the CST file holds more than
%               one signal. Pass [] or "" for the default when also
%               giving options.
%
%   Options
%     'Signal'    : which signal to write when the file holds several.
%                   An index (1, 2, ...) or a CST signal name ('i2').
%                   Default 1 (the first signal).
%     'Z0'        : reference impedance in ohms, used only if the CST file
%                   has no Z0 in its #Parameters line
%     'Header'    : true writes a first row "time_ns,voltage_V"
%                   (default false: numbers only)
%     'TimeUnit'  : time unit of the CSV, 'ns' (default), 'ps', 'us'
%                   or 's'
%     'Precision' : significant digits (default 12)
%
%   Example
%     cst_to_csv("C:\path\to\CST\project\Export\excitation_port1.txt");
%     cst_to_csv(cstFile, "port2.csv", 'Signal', 'i2');

if nargin < 2, outFile = []; end
p = inputParser;
addParameter(p, 'Signal',    1);
addParameter(p, 'Z0',        []);
addParameter(p, 'Header',    false);
addParameter(p, 'Precision', 12);
addParameter(p, 'TimeUnit',  'ns');
parse(p, varargin{:});
o = p.Results;

switch lower(char(o.TimeUnit))
    case 's',  tFactor = 1;
    case 'ms', tFactor = 1e3;
    case 'us', tFactor = 1e6;
    case 'ns', tFactor = 1e9;
    case 'ps', tFactor = 1e12;
    otherwise
        error('cst_to_csv:timeUnit', 'TimeUnit must be s, ms, us, ns or ps.');
end
tUnit = lower(char(o.TimeUnit));

% ---- read ----------------------------------------------------------------
cstFile = char(cstFile);
if ~isfile(cstFile)
    error('cst_to_csv:notFound', 'CST file not found: %s', cstFile);
end
try
    [t, v, info] = importCSTSignal(cstFile, o.Z0);
catch err
    if strcmp(err.identifier, 'importCSTSignal:noZ0')
        error('cst_to_csv:noZ0', ...
            'No Z0 in %s. Pass it with ''Z0'', e.g. cst_to_csv(file, [], ''Z0'', 50).', cstFile);
    end
    rethrow(err);
end
nSig = numel(t);

% ---- pick the signal -----------------------------------------------------
sel = o.Signal;
if ischar(sel) || isstring(sel)
    k = find(strcmpi(info.names, char(sel)), 1);
    if isempty(k)
        error('cst_to_csv:signal', 'Signal "%s" not in file. Available: %s', ...
            char(sel), strjoin(info.names, ', '));
    end
else
    k = sel;
    if ~(isscalar(k) && k == round(k) && k >= 1 && k <= nSig)
        error('cst_to_csv:signal', 'Signal index must be 1 to %d. Available: %s', ...
            nSig, strjoin(info.names, ', '));
    end
end
if nSig > 1
    fprintf('%s holds %d signals (%s); writing %s.\n', cstFile, nSig, ...
        strjoin(info.names, ', '), info.names{k});
end

% ---- output name ---------------------------------------------------------
if isempty(outFile) || strlength(string(outFile)) == 0
    [~, srcName] = fileparts(cstFile);
    if nSig > 1
        srcName = [srcName '_' matlab.lang.makeValidName(info.names{k})];
    end
    outFile = [srcName '.csv'];
end
outFile = char(outFile);
[outDir, ~, ext] = fileparts(outFile);
if isempty(outDir)                       % bare file name -> output/cst_ads/csv_out
    outDir  = fullfile(repo_root, 'output', 'cst_ads', 'csv_out');
    outFile = fullfile(outDir, outFile);
end
if isempty(ext)
    outFile = [outFile '.csv'];
end
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end

% ---- write ---------------------------------------------------------------
fid = fopen(outFile, 'w');
if fid < 0
    error('cst_to_csv:open', 'Could not open %s for writing.', outFile);
end
closer = onCleanup(@() fclose(fid));

if o.Header
    fprintf(fid, 'time_%s,voltage_V\n', tUnit);
end
fmt = sprintf('%%.%dg,%%.%dg\n', o.Precision, o.Precision);
fprintf(fid, fmt, [t{k} * tFactor, v{k}].');   % importCSTSignal gives seconds

if strcmp(info.sourceType{k}, 'generator')
    how = sprintf('generator signal [%s] * sqrt(Z0), Z0 = %g ohm', info.signalUnit{k}, info.Z0(k));
else
    how = sprintf('voltage source, values in %s used as-is', info.signalUnit{k});
end
fprintf('Wrote %s  (%s, %d pts, %.4g to %.4g %s; %s)\n', outFile, ...
    info.names{k}, numel(t{k}), t{k}(1)*tFactor, t{k}(end)*tFactor, tUnit, how);

if nargout == 0, clear outFile; end
end
