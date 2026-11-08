function outFile = writeCSTSignal(t, v, outFile, varargin)
%WRITECSTSIGNAL  Write a time signal as a CST ASCII signal file (for import
%as a user-defined excitation signal).
%
%   writeCSTSignal(t, v, outFile)
%   writeCSTSignal(dt, v, outFile)            % uniform sampling, starts at 0
%   writeCSTSignal(t, v, outFile, 'TimeUnit', 'ns', 'Name', 'i1', ...
%                  'Unit', 'V', 'Parameters', "Z0=50; a=1")
%
%   The file uses the same layout CST writes for its own signal exports:
%       #Parameters = {Z0=50; a=1}            (only if 'Parameters' given)
%       #"Time / ns"	"i1 [Real / V]"
%       #---------------------------------
%       0.000000000000	0.000000000000
%       ...
%   Lines starting with # are comments for CST; the data are two columns,
%   time and amplitude.
%
%   Inputs
%     t       : N x 1 time vector in SECONDS, or a scalar time step dt in
%               seconds (then t = (0:N-1)*dt)
%     v       : N x 1 signal amplitude
%     outFile : output file name. A bare name (no folder) goes to
%               <repo>/output/cst_ads/cst_out (created if needed);
%               '.txt' is added if there is no extension.
%
%   Options
%     'TimeUnit'   : unit of the time column, which must match the time
%                    unit of the CST project: 'ns' (default), 'ps', 'us',
%                    'ms', 's', 'fs'
%     'Name'       : signal name in the header (default 'i1')
%     'Unit'       : amplitude unit in the header (default 'V').
%                    CST port signals are 'V' (voltage source) or
%                    'W^(1/2)' (generator signal); this function only
%                    writes the label, the values are written as given.
%     'Parameters' : text for a #Parameters = {...} line, e.g. "Z0=50"
%                    (default: no parameter line)
%     'Precision'  : significant digits (default 12)
%
%   Notes for CST
%     - The signal should start at t = 0.
%     - CST resamples the signal to its own solver time step, so the
%       sampling here only needs to resolve the signal's bandwidth.
%     - A signal that does not end at zero is cut off abruptly when the
%       simulation passes its last sample; pad with zeros if needed
%       (see zeroPadCSTSignal).
%
%   Example
%     dt = 1e-12;  t = (0:4999)'*dt;
%     v  = exp(-((t-0.5e-9)/0.1e-9).^2) .* sin(2*pi*3e9*t);
%     writeCSTSignal(t, v, "gauss_sine_3GHz");

p = inputParser;
addParameter(p, 'TimeUnit',   'ns');
addParameter(p, 'Name',       'i1');
addParameter(p, 'Unit',       'V');
addParameter(p, 'Parameters', '');
addParameter(p, 'Precision',  12);
parse(p, varargin{:});
o = p.Results;

% ---- inputs ------------------------------------------------------------
v = v(:);
if isscalar(t)
    t = (0:numel(v)-1).' * t;
else
    t = t(:);
end
if numel(t) ~= numel(v)
    error('writeCSTSignal:size', 't and v must have the same length (%d vs %d).', ...
        numel(t), numel(v));
end
% values below double precision relative to the peak (e.g. 1e-300 in a
% Gaussian tail) are written as exact zeros so CST parses them cleanly
v(abs(v) < eps * max(abs(v))) = 0;
if any(diff(t) <= 0)
    error('writeCSTSignal:time', 'Time must be strictly increasing.');
end
if abs(t(1)) > 1e-3 * median(diff(t))
    warning('writeCSTSignal:start', ...
        'Signal starts at t = %g s, not 0. CST excitation signals normally start at 0.', t(1));
end
if abs(v(end)) > 1e-3 * max(abs(v))
    warning('writeCSTSignal:end', ...
        ['Signal does not end at zero (last value %g, peak %g); CST will cut it ' ...
         'off abruptly after the last sample.'], v(end), max(abs(v)));
end

tu = lower(char(o.TimeUnit));
switch tu
    case 's',  tf = 1;
    case 'ms', tf = 1e3;
    case 'us', tf = 1e6;
    case 'ns', tf = 1e9;
    case 'ps', tf = 1e12;
    case 'fs', tf = 1e15;
    otherwise
        error('writeCSTSignal:timeUnit', 'TimeUnit must be s, ms, us, ns, ps or fs.');
end

% ---- output name -------------------------------------------------------
outFile = char(outFile);
[outDir, ~, ext] = fileparts(outFile);
if isempty(outDir)
    outDir  = fullfile(repo_root, 'output', 'cst_ads', 'cst_out');
    outFile = fullfile(outDir, outFile);
end
if isempty(ext)
    outFile = [outFile '.txt'];
end
if ~isfolder(outDir)
    mkdir(outDir);
end

% ---- write ---------------------------------------------------------------
fid = fopen(outFile, 'w');
if fid < 0
    error('writeCSTSignal:open', 'Could not open %s for writing.', outFile);
end
closer = onCleanup(@() fclose(fid));

params = strtrim(char(o.Parameters));
if ~isempty(params)
    params = regexprep(params, '^\{|\}$', '');      % accept with or without braces
    fprintf(fid, '#Parameters = {%s}\n', params);
end
fprintf(fid, '#"Time / %s"\t"%s [Real / %s]"\n', tu, char(o.Name), char(o.Unit));
fprintf(fid, '#---------------------------------\n');
fmt = sprintf('%%.%dg\\t%%.%dg\\n', o.Precision, o.Precision);
fprintf(fid, fmt, [t * tf, v].');

fprintf('Wrote %s  (%d pts, %.4g to %.4g %s)\n', outFile, numel(t), ...
    t(1)*tf, t(end)*tf, tu);

if nargout == 0, clear outFile; end
end
