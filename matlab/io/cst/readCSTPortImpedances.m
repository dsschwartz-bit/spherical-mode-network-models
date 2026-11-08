function P = readCSTPortImpedances(file)
%READCSTPORTIMPEDANCES  Port reference impedances from a CST 1D export.
%
%   P = readCSTPortImpedances(file)
%
%   Reads blocks such as
%       #"Frequency / GHz"   "Z1 [Real / Ohm]"
%       #"Frequency / GHz"   "Z1 [Imaginary / Ohm]"     (optional)
%   one per port and part. A port without an imaginary block is taken as
%   purely real.
%
%   P : struct array, one element per port, with fields
%       port  port number
%       f     frequency (Hz), column
%       Z     complex impedance (ohm) at f
%       Zc    the single value to use in the time domain (mean of Z)
%       isConstant  true if Z varies by less than 0.1 % over f
%       isReal      true if |Im Z| < 1e-6 |Z|

txt = fileread(file);
blk = regexp(txt, '#"([^"]*)"\s*"([^"]*)"\s*\r?\n#-+\s*\r?\n([^#]*)', 'tokens');
if isempty(blk), error('No impedance blocks in %s.', file); end

S = struct('port', {}, 'part', {}, 'f', {}, 'x', {});
for k = 1:numel(blk)
    xname = blk{k}{1};  yname = blk{k}{2};
    pn = regexp(yname, '^\s*Z\s*(\d+)', 'tokens', 'once');
    if isempty(pn), continue; end
    u = regexp(xname, '/\s*(\w+)', 'tokens', 'once');
    fs = 1;
    if ~isempty(u)
        switch lower(u{1})
            case 'hz', fs = 1;  case 'khz', fs = 1e3;  case 'mhz', fs = 1e6;
            case 'ghz', fs = 1e9;  case 'thz', fs = 1e12;
        end
    end
    if contains(yname, 'Imag', 'IgnoreCase', true), part = "imag"; else, part = "real"; end
    M = sscanf(blk{k}{3}, '%f');  M = reshape(M, 2, []).';
    S(end+1) = struct('port', str2double(pn{1}), 'part', part, 'f', M(:,1)*fs, 'x', M(:,2)); %#ok<AGROW>
end

ports = unique([S.port]);
P = struct('port', {}, 'f', {}, 'Z', {}, 'Zc', {}, 'isConstant', {}, 'isReal', {});
for p = ports
    r = S([S.port] == p & [S.part] == "real");
    if isempty(r), error('Port %d: no real part in %s.', p, file); end
    f = r(1).f;  Z = r(1).x;
    m = S([S.port] == p & [S.part] == "imag");
    if ~isempty(m), Z = Z + 1j * interp1(m(1).f, m(1).x, f, 'linear', 'extrap'); end
    Zc = mean(Z);
    P(end+1) = struct('port', p, 'f', f, 'Z', Z, 'Zc', Zc, ...
        'isConstant', max(abs(Z - Zc)) <= 1e-3 * abs(Zc), ...
        'isReal', all(abs(imag(Z)) <= 1e-6 * abs(Z))); %#ok<AGROW>
end
end
