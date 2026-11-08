function [time, voltage, info] = importCSTSignal(filename, Z0, varargin)
%IMPORTCSTSIGNAL  Import CST time-domain port signals (one or many) from an ASCII export.
%
%   [time, voltage, info] = importCSTSignal(filename)
%   [time, voltage, info] = importCSTSignal(filename, Z0)
%   [time, voltage, info] = importCSTSignal(filename, Z0, 'TimeUnit', 'ns')
%
%   The export may contain any number of signals. CST writes each one as
%   its own block with its own header, e.g.
%       #Parameters = {Z0=376.73; a=39.71; ...}
%       #"Time / ns"	"i1 [Real / W^(1/2)]"
%       #---------------------------------
%       0.000   -3.27e-08
%       ...
%   Blocks with several signal columns side by side are also handled.
%
%   The signal unit tells which kind of excitation it is:
%     W^(1/2) : generator open-circuit voltage with source impedance equal
%               to the port impedance. Converted to volts as
%               value * sqrt(Z0); Z0 is required.
%     V       : ideal voltage source (zero source impedance). The values
%               are already volts and are used as-is; Z0 is not needed.
%     V on a "Port N [exc]" column : total port voltage sqrt(Z0)*(a+b) of a
%               port excitation (CST 'Port voltages' export), sourceType
%               'port voltage', name 'vN [exc]'. Usable as an ideal-voltage-
%               source excitation (substitution theorem).
%     A on a "Port N [exc]" column : port current into the port (CST 'Port
%               currents' export), sourceType 'port current', name 'cN [exc]'.
%               With the port voltage it gives the power waves (cstVItoWaves.m).
%
%   Inputs
%     filename : path to the exported text file
%     Z0       : (optional) port impedance in ohms. Used only for W^(1/2)
%                signals, and only if the file has no Z0 in its
%                #Parameters line. If the file does contain Z0, the file
%                value is used (with a warning if the value you passed
%                differs).
%
%   Name/value options
%     'TimeUnit'   : force the time unit ('s','ms','us','ns','ps','fs')
%                    if the header is missing or wrong
%     'CommonTime' : false (default) keeps each signal on its own time
%                    grid. true resamples all signals onto one shared grid
%                    (finest step, full span, zero outside a signal's
%                    range) and returns time as a vector and voltage as an
%                    N x K matrix.
%
%   Outputs (default, CommonTime = false)
%     time    : 1 x K cell array, time{k} is an Nk x 1 vector in SECONDS
%     voltage : 1 x K cell array, voltage{k} in volts
%     info    : struct with
%                 .names       signal names, e.g. {'i1','i2','o1,1'}; port-wave names are
%                              returned in lower case (co-simulation exports write I1, O1,1)
%                 .rawNames    names exactly as written in the file
%                 .port        port number of each signal (NaN if none)
%                 .excitation  excitation tag ('' if none), e.g. '1', '2',
%                              'cosim Tran1 1'; tagged names read 'v1 [cosim Tran1 1]'
%                 .signalUnit  signal unit for each signal ('W^(1/2)', 'V')
%                 .sourceType  'generator' (W^(1/2)), 'voltage source' (V) or
%                              'port voltage' (V, total voltage at a port) or
%                              'port current' (A, current into a port)
%                 .scale       factor applied to the CST values
%                              (sqrt(Z0) for generator, 1 for voltage source)
%                 .Z0          Z0 used for each signal (NaN if not needed)
%                 .Z0source    'file', 'input' or 'not used'
%                 .timeUnit    time unit found for each signal
%                 .parameters  struct of #Parameters values (first block)
%                 .resampled   true if CommonTime had to resample
%                 .blocks      per-block details (parameters, header)
%
%   Example
%     [t, v, info] = importCSTSignal('excitation_ports.txt');
%     figure; hold on
%     for k = 1:numel(t)
%         plot(t{k}*1e9, v{k}, 'DisplayName', info.names{k});
%     end
%     xlabel('Time (ns)'); ylabel('Voltage (V)'); legend

% ---- inputs -------------------------------------------------------------
if nargin < 2, Z0 = []; end
if ~isempty(Z0)
    validateattributes(Z0, {'numeric'}, {'scalar', 'real', 'positive', 'finite'}, ...
        mfilename, 'Z0');
end
p = inputParser;
addParameter(p, 'TimeUnit',   '',    @(s) ischar(s) || isstring(s));
addParameter(p, 'CommonTime', false, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
forcedUnit = char(p.Results.TimeUnit);

% ---- read file and split into header/data blocks -----------------------
txt = fileread(filename);
lines = regexp(txt, '\r?\n', 'split');
lines = strtrim(lines);
lines = lines(~cellfun(@isempty, lines));
if isempty(lines)
    error('importCSTSignal:empty', 'File %s is empty.', filename);
end

isHdr = startsWith(lines, '#');
% a block starts at a header line that follows a data line (or the file start)
blockStart = find(isHdr & [true, ~isHdr(1:end-1)]);
if isempty(blockStart) || blockStart(1) ~= 1
    blockStart = [1, blockStart];   % data before any header
end
blockStart = unique(blockStart);
blockEdges = [blockStart, numel(lines) + 1];

time = {}; voltage = {};
names = {}; rawNames = {}; Z0used = []; Z0src = {}; tUnits = {}; sUnits = {};
srcTypes = {}; scales = []; portNums = []; excTags = {};
blocks = struct('parameters', {}, 'header', {}, 'columnNames', {});

for b = 1:numel(blockStart)
    L = lines(blockEdges(b):blockEdges(b+1)-1);
    hdr = L(startsWith(L, '#'));
    dat = L(~startsWith(L, '#'));
    if isempty(dat), continue; end

    % --- numeric data
    nCol = numel(sscanf(dat{1}, '%f'));
    if nCol < 2
        error('importCSTSignal:format', ...
            'Block %d: expected at least 2 numeric columns.', b);
    end
    vals = sscanf(strjoin(dat, ' '), '%f');
    if mod(numel(vals), nCol) ~= 0
        error('importCSTSignal:format', ...
            'Block %d: rows have inconsistent numbers of columns.', b);
    end
    M = reshape(vals, nCol, []).';

    % --- parameters and column names
    params   = parseParameters(hdr);
    colNames = {};
    for k = 1:numel(hdr)
        tok = regexp(hdr{k}, '"([^"]*)"', 'tokens');
        if numel(tok) >= 2
            colNames = cellfun(@(c) c{1}, tok, 'UniformOutput', false);
            break
        end
    end

    % --- time unit
    if ~isempty(forcedUnit)
        tu = forcedUnit;
    elseif ~isempty(colNames)
        u = regexp(colNames{1}, '(?:/|\[)\s*([^\]\s]+)', 'tokens', 'once');
        if isempty(u)
            error('importCSTSignal:noTimeUnit', ...
                'Block %d: no time unit in "%s". Use ''TimeUnit''.', b, colNames{1});
        end
        tu = u{1};
    else
        error('importCSTSignal:noTimeUnit', ...
            'Block %d has no header. Pass the unit with ''TimeUnit''.', b);
    end
    t = M(:, 1) * unitToSeconds(tu);

    % --- one output per signal column
    for c = 2:nCol
        tag = '';
        if numel(colNames) >= c
            % label = header without the trailing "[Real / unit]"; an excitation
            % tag before it ("Port 1 [cosim Tran1 1] [Real / V]") is kept apart
            lbl = strtrim(regexprep(colNames{c}, '\[[^\]]*/[^\]]*\]\s*$', ''));
            tg  = regexp(lbl, '\[([^\]]*)\]\s*$', 'tokens', 'once');
            if ~isempty(tg)
                tag = strtrim(tg{1});
                lbl = strtrim(regexprep(lbl, '\[[^\]]*\]\s*$', ''));
            end
            nm = strtrim(regexprep(lbl, '\[.*$', ''));
            u  = regexp(colNames{c}, '/\s*([^\]]+)\]', 'tokens', 'once');
            if isempty(u), u = regexp(colNames{c}, '\[([^\]]+)\]', 'tokens', 'once'); end
            if isempty(u), su = ''; else, su = strtrim(u{1}); end
        else
            nm = sprintf('signal%d', numel(names) + 1);  su = '';
        end
        pnum = regexp(nm, '^Port\s*(\d+)$', 'tokens', 'once');   % "Port 1" = port voltage export
        rawNm = nm;  if ~isempty(tag), rawNm = [nm ' [' tag ']']; end

        % --- excitation type from the unit
        suc = lower(strrep(su, ' ', ''));
        if any(strcmp(suc, {'w^(1/2)', 'w^0.5', 'sqrt(w)', 'w^1/2'}))
            stype = 'generator';
            [z, zsrc] = pickZ0(params, Z0, b);
            scale = sqrt(z);
        elseif any(strcmp(suc, {'v', 'volt', 'volts'})) && ~isempty(pnum)
            % "Port N [exc] [Real / V]": total port voltage sqrt(Z0)*(a+b) of a
            % port excitation, NOT the incident wave
            stype = 'port voltage';
            z = NaN;  zsrc = 'not used';  scale = 1;
            if isfield(params, 'Z0') && isnumeric(params.Z0), z = params.Z0;  zsrc = 'file'; end
            nm = ['v' pnum{1}];
        elseif any(strcmp(suc, {'a', 'amp', 'ampere', 'amperes'})) && ~isempty(pnum)
            % "Port N [exc] [Real / A]": port current (into the port)
            stype = 'port current';
            z = NaN;  zsrc = 'not used';  scale = 1;
            if isfield(params, 'Z0') && isnumeric(params.Z0), z = params.Z0;  zsrc = 'file'; end
            nm = ['c' pnum{1}];
        elseif any(strcmp(suc, {'v', 'volt', 'volts'}))
            stype = 'voltage source';
            z = NaN;  zsrc = 'not used';  scale = 1;
        elseif isempty(su)
            error('importCSTSignal:noSignalUnit', ...
                ['Signal "%s" has no unit in the header, so it is unclear whether ' ...
                 'it is a W^(1/2) generator signal or a voltage.'], nm);
        else
            error('importCSTSignal:signalUnit', ...
                'Signal "%s" has unit "%s"; expected W^(1/2) (generator) or V (voltage source).', ...
                nm, su);
        end

        time{end+1}     = t;                        %#ok<AGROW>
        voltage{end+1}  = M(:, c) * scale;          %#ok<AGROW>
        rawNames{end+1} = rawNm;                    %#ok<AGROW>
        if ~isempty(regexp(nm, '^[IO]\d', 'once')), nm = [lower(nm(1)) nm(2:end)]; end   % I1 -> i1, O1,1 -> o1,1
        if ~isempty(pnum), pn = str2double(pnum{1});
        else, pn = str2double(regexp(nm, '^[io](\d+)', 'tokens', 'once')); end
        if isempty(pn), pn = NaN; end
        portNums(end+1) = pn;                       %#ok<AGROW>
        excTags{end+1}  = tag;                      %#ok<AGROW>
        if ~isempty(tag), nm = [nm ' [' tag ']']; end   % v1 [1], v1 [cosim Tran1 1]
        names{end+1}    = nm;                       %#ok<AGROW>
        sUnits{end+1}   = su;                       %#ok<AGROW>
        srcTypes{end+1} = stype;                    %#ok<AGROW>
        scales(end+1)   = scale;                    %#ok<AGROW>
        Z0used(end+1)   = z;                        %#ok<AGROW>
        Z0src{end+1}    = zsrc;                     %#ok<AGROW>
        tUnits{end+1}   = tu;                       %#ok<AGROW>
    end

    blocks(end+1) = struct('parameters', params, 'header', {hdr}, ...
                           'columnNames', {colNames}); %#ok<AGROW>
end

if isempty(time)
    error('importCSTSignal:noData', 'No numeric data found in %s.', filename);
end

% ---- optional resampling onto one shared grid --------------------------
resampled = false;
if p.Results.CommonTime
    sameGrid = all(cellfun(@(x) isequal(size(x), size(time{1})) && ...
        max(abs(x - time{1})) <= 1e-9 * max(abs(time{1})), time));
    if sameGrid
        tq = time{1};
    else
        dt = min(cellfun(@(x) median(diff(x)), time));
        tq = (min(cellfun(@min, time)) : dt : max(cellfun(@max, time))).';
        resampled = true;
    end
    V = zeros(numel(tq), numel(voltage));
    for k = 1:numel(voltage)
        if sameGrid
            V(:, k) = voltage{k};
        else
            [tu_, iu] = unique(time{k});
            V(:, k) = interp1(tu_, voltage{k}(iu), tq, 'linear', 0);
        end
    end
    time = tq;  voltage = V;
end

info = struct( ...
    'names',      {names}, ...
    'rawNames',   {rawNames}, ...
    'port',       portNums, ...
    'excitation', {excTags}, ...
    'signalUnit', {sUnits}, ...
    'sourceType', {srcTypes}, ...
    'scale',      scales, ...
    'Z0',         Z0used, ...
    'Z0source',   {Z0src}, ...
    'timeUnit',   {tUnits}, ...
    'parameters', blocks(1).parameters, ...
    'resampled',  resampled, ...
    'blocks',     blocks);
end

% ========================================================================
function [z, zsrc] = pickZ0(params, Z0, b)
% Z0 for a W^(1/2) signal: file value first, then the Z0 argument
if isfield(params, 'Z0') && isnumeric(params.Z0) && params.Z0 > 0
    z = params.Z0;  zsrc = 'file';
    if ~isempty(Z0) && abs(Z0 - z) > 1e-9 * z
        warning('importCSTSignal:Z0mismatch', ...
            'Block %d: file Z0 = %g ohm overrides the %g ohm you passed.', b, z, Z0);
    end
elseif ~isempty(Z0)
    z = Z0;  zsrc = 'input';
else
    error('importCSTSignal:noZ0', ...
        ['Block %d: signal is in W^(1/2) but no Z0 was found in the file. ' ...
         'Pass it explicitly, e.g. importCSTSignal(file, 50).'], b);
end
end

% ========================================================================
function params = parseParameters(hdr)
% parse "#Parameters = {a=39.71; feedW=1; ...}" into a struct
params = struct();
for k = 1:numel(hdr)
    body = regexp(hdr{k}, 'Parameters\s*=\s*\{(.*)\}', 'tokens', 'once');
    if isempty(body), continue; end
    pairs = strsplit(body{1}, ';');
    for j = 1:numel(pairs)
        kv = strsplit(strtrim(pairs{j}), '=');
        if numel(kv) ~= 2, continue; end
        name = matlab.lang.makeValidName(strtrim(kv{1}));
        val  = str2double(strtrim(kv{2}));
        if isnan(val), val = strtrim(kv{2}); end
        params.(name) = val;
    end
    return
end
end

% ========================================================================
function s = unitToSeconds(u)
switch lower(strtrim(u))
    case {'s', 'sec', 'second', 'seconds'}
        s = 1;
    case {'ms'}
        s = 1e-3;
    case {'us', [char(181) 's'], [char(956) 's']}   % us, µs, μs
        s = 1e-6;
    case {'ns'}
        s = 1e-9;
    case {'ps'}
        s = 1e-12;
    case {'fs'}
        s = 1e-15;
    otherwise
        error('importCSTSignal:unknownUnit', ...
            'Unrecognized time unit "%s". Use ''TimeUnit'' to specify it.', u);
end
end
