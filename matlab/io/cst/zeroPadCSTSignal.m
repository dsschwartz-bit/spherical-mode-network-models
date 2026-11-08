function outFile = zeroPadCSTSignal(cstFile, tEnd, outFile, varargin)
%ZEROPADCSTSIGNAL  Extend a CST time signal with zeros to a later end time
%and write it back in CST ASCII format (for re-import as an excitation).
%
%   zeroPadCSTSignal(cstFile, tEnd)
%   zeroPadCSTSignal(cstFile, tEnd, outFile)
%   zeroPadCSTSignal(cstFile, tEnd, outFile, 'Signal', 2, 'TimeUnit', 'ps')
%
%   The original samples are kept exactly. New samples are appended after
%   the last one at the signal's median time step, with amplitude 0, up
%   to tEnd. The amplitude is written in its original CST unit (V or
%   W^(1/2), no sqrt(Z0) scaling) with the original signal name and
%   #Parameters line, so the padded file describes the same excitation.
%
%   Inputs
%     cstFile : CST ASCII signal export
%     tEnd    : new end time, in 'TimeUnit' (by default the time unit of
%               the CST file, normally ns), e.g. 20 for 20 ns
%     outFile : (optional) output name. Default:
%                 <cst name>_zeropad_<tEnd><unit>.txt
%               in the same folder as the CST file, with tEnd rounded to
%               3 significant figures (e.g. excitation_zeropad_20ns.txt).
%               A bare name (no folder) is also written to the CST
%               file's folder; a name with a folder is written there.
%               Pass [] or "" for the default when giving options.
%
%   Options
%     'Signal'   : which signal if the file holds several (index or CST
%                  name, default 1); one signal is written per file
%     'TimeUnit' : time unit of tEnd AND of the output file's time
%                  column ('s','ms','us','ns','ps','fs'). Default: the
%                  time unit of the CST file.
%     'Z0'       : only needed if a W^(1/2) file has no Z0 parameter
%                  (used while reading; the output stays in W^(1/2))
%
%   Example
%     zeroPadCSTSignal("C:\path\to\CST\project\Export\time_signals\excitation_sinestep.txt", 20);   % to 20 ns

p = inputParser;
addParameter(p, 'Signal',   1);
addParameter(p, 'TimeUnit', '');
addParameter(p, 'Z0',       []);
parse(p, varargin{:});
o = p.Results;
if nargin < 3, outFile = []; end

% ---- read -----------------------------------------------------------------
[tc, vc, info] = importCSTSignal(cstFile, o.Z0);
k = o.Signal;
if ischar(k) || isstring(k)
    k = find(strcmpi(info.names, char(k)), 1);
    if isempty(k)
        error('zeroPadCSTSignal:signal', 'Signal not found. Available: %s', ...
            strjoin(info.names, ', '));
    end
end
t   = tc{k}(:);
raw = vc{k}(:) / info.scale(k);       % back to the CST unit (V or W^(1/2))

% ---- tEnd in seconds -------------------------------------------------------
tu = char(o.TimeUnit);
if isempty(tu), tu = info.timeUnit{k}; end
tuScale = struct('s',1,'ms',1e-3,'us',1e-6,'ns',1e-9,'ps',1e-12,'fs',1e-15);
if ~isfield(tuScale, lower(tu))
    error('zeroPadCSTSignal:timeUnit', 'TimeUnit must be s, ms, us, ns, ps or fs.');
end
tEndUser = tEnd;
tEnd = tEnd * tuScale.(lower(tu));    % seconds from here on

% ---- pad ------------------------------------------------------------------
dt = median(diff(t));
nPad = (tEnd - t(end)) / dt;
if nPad > 5e7
    error('zeroPadCSTSignal:tooLong', ...
        ['Padding to %g %s would add %.3g samples (dt = %.4g ps). tEnd is in %s; ' ...
         'check the value or the ''TimeUnit'' option.'], ...
        tEndUser, tu, nPad, dt*1e12, tu);
end
if tEnd <= t(end)
    warning('zeroPadCSTSignal:noPad', ...
        'tEnd (%g s) is not after the last sample (%g s); writing the signal unchanged.', ...
        tEnd, t(end));
    tPad = zeros(0, 1);
else
    tPad = (t(end) + dt : dt : tEnd).';
    if isempty(tPad) || tPad(end) < tEnd - 1e-6*dt
        tPad(end+1, 1) = tEnd;
    end
end
tNew = [t; tPad];
vNew = [raw; zeros(numel(tPad), 1)];

if abs(raw(end)) > 1e-3 * max(abs(raw))
    warning('zeroPadCSTSignal:step', ...
        ['The signal is not near zero at its last sample (%g, peak %g), so ' ...
         'padding with zeros adds a step at t = %g s.'], raw(end), max(abs(raw)), t(end));
end

% ---- header details from the original file --------------------------------
params = '';
blk = info.blocks(min(k, numel(info.blocks)));
for j = 1:numel(blk.header)
    tok = regexp(blk.header{j}, 'Parameters\s*=\s*\{(.*)\}', 'tokens', 'once');
    if ~isempty(tok), params = tok{1}; break; end
end

% ---- output name ----------------------------------------------------------
[srcDir, srcName] = fileparts(char(cstFile));
if isempty(srcDir), srcDir = pwd; end
if isempty(outFile) || strlength(string(outFile)) == 0
    tName = num2str(round(tEndUser, 3, 'significant'));   % 3 significant figures
    outFile = sprintf('%s_zeropad_%s%s.txt', srcName, tName, lower(tu));
end
outFile = char(outFile);
if isempty(fileparts(outFile))           % bare name -> next to the CST file
    outFile = fullfile(srcDir, outFile);
end

outFile = writeCSTSignal(tNew, vNew, outFile, 'TimeUnit', tu, ...
    'Name', info.names{k}, 'Unit', info.signalUnit{k}, 'Parameters', params);
fprintf('  padded %s from %.4g to %.4g ns (%d zero samples added, dt = %.4g ps)\n', ...
    info.names{k}, t(end)*1e9, tNew(end)*1e9, numel(tPad), dt*1e12);

if nargout == 0, clear outFile; end
end
