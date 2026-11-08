function W = cstVItoWaves(vFile, iFile, outFile, Z)
%CSTVITOWAVES  Power waves from CST port voltage and current exports.
%
%   W = cstVItoWaves(vFile, iFile)
%   W = cstVItoWaves(vFile, iFile, outFile)       also writes a CST-layout file
%   W = cstVItoWaves(vFile, iFile, outFile, Z)    port impedance(s) Z
%
%   vFile : CST 'Port voltages' export, columns "Port N [exc] [Real / V]"
%   iFile : CST 'Port currents' export, columns "Port N [exc] [Real / A]"
%           (current into the port), same ports and excitations
%   Z     : port reference impedance(s): a scalar (all ports), a vector
%           indexed by port number, or a CST port-impedance export
%           ("Z1 [Real / Ohm]", ..., see readCSTPortImpedances). Default:
%           Z0 from the #Parameters line for all ports. Each port's
%           impedance must be real and frequency-independent for
%           time-domain waves.
%
%   Kurokawa power waves, for each port N and excitation:
%       a = (V + Z  I) / (2 sqrt(|Re Z|))     incident (forward)
%       b = (V - Z* I) / (2 sqrt(|Re Z|))     outgoing (reflected)
%   in W^(1/2). In the time domain Z must be real (a frequency-independent
%   resistance), so Z* = Z.
%
%   The incident wave on the driven port is the CST excitation; on the
%   other ports it is ~0 (matched). Excitation tags ending in a port number
%   ("1", "cosim Tran1 2") give the driven port, so the outgoing wave of
%   port N is named oN,<driven port>, as in CST.
%
%   outFile (optional): written in the layout of a CST time-signal export
%   (one block per signal, time in ns), signals "iN [exc]" and
%   "oN,M [exc]" in W^(1/2), readable by importCSTSignal. Each block's
%   #Parameters line carries Z0 = that port's impedance, so importCSTSignal
%   converts each wave to volts with its own port's sqrt(Z0).
%
%   W : struct array with fields port, excitation, drivenPort, t (s),
%       a, b (W^(1/2)), iName, oName, Z (port impedance used), Zsource

[tv, vv, iv] = importCSTSignal(vFile);
[ti, ci, ii] = importCSTSignal(iFile);
if ~all(strcmp(iv.sourceType, 'port voltage'))
    error('%s is not a port-voltage export ("Port N [exc] [Real / V]").', vFile);
end
if ~all(strcmp(ii.sourceType, 'port current'))
    error('%s is not a port-current export ("Port N [exc] [Real / A]").', iFile);
end
if nargin < 3, outFile = ""; end
if nargin < 4, Z = []; end
% port impedance lookup: Zp(port) -> ohm
ports = unique(iv.port);
if ischar(Z) || isstring(Z)
    PZ = readCSTPortImpedances(Z);
    [~, zn, ze] = fileparts(Z);  zsrc = string(zn) + ze;
    Zp = nan(1, max(ports));
    for p = ports
        q = PZ([PZ.port] == p);
        if isempty(q), error('Port %d not in %s.', p, Z); end
        if ~q.isReal || ~q.isConstant
            error(['Port %d impedance in %s varies with frequency or is complex (%.4g%+.4gj to %.4g%+.4gj ohm);\n' ...
                'time-domain power waves need a constant real impedance.'], p, Z, ...
                real(min(q.Z)), imag(min(q.Z)), real(max(q.Z)), imag(max(q.Z)));
        end
        Zp(p) = real(q.Zc);
    end
elseif isempty(Z)
    if ~isfield(iv.parameters, 'Z0'), error('No Z0 in %s; pass Z.', vFile); end
    Zp = repmat(iv.parameters.Z0, 1, max(ports));  zsrc = "#Parameters Z0";
elseif isscalar(Z)
    Zp = repmat(Z, 1, max(ports));  zsrc = "input";
else
    if numel(Z) < max(ports)
        error('Z has %d values but the export has port %d; give one impedance per port.', numel(Z), max(ports));
    end
    Zp = Z(:).';  zsrc = "input";
end
if any(~isreal(Zp(ports))) || any(Zp(ports) <= 0) || any(isnan(Zp(ports)))
    error('Port impedances must be positive real resistances for time-domain power waves.');
end

W = struct('port', {}, 'excitation', {}, 'drivenPort', {}, 't', {}, 'a', {}, 'b', {}, ...
    'iName', {}, 'oName', {}, 'Z', {}, 'Zsource', {});
for k = 1:numel(tv)
    j = find(ii.port == iv.port(k) & strcmp(ii.excitation, iv.excitation{k}), 1);
    if isempty(j)
        error('No current for port %d, excitation "%s" in %s.', iv.port(k), iv.excitation{k}, iFile);
    end
    t = tv{k};  V = vv{k};  I = ci{j};
    if numel(ti{j}) ~= numel(t) || max(abs(ti{j} - t)) > 1e-6 * max(abs(t))
        [tu, iu] = unique(ti{j});
        I = interp1(tu, I(iu), t, 'pchip', 0);               % different grid: resample the current
    end
    Zk = Zp(iv.port(k));
    a = (V + Zk*I) / (2*sqrt(abs(real(Zk))));
    b = (V - conj(Zk)*I) / (2*sqrt(abs(real(Zk))));
    tag = iv.excitation{k};
    tok = regexp(tag, '(\d+)\s*$', 'tokens', 'once');
    if isempty(tok), drv = NaN; else, drv = str2double(tok{1}); end
    p = iv.port(k);
    if isempty(tag), sfx = ''; else, sfx = [' [' tag ']']; end
    iName = sprintf('i%d%s', p, sfx);
    if isnan(drv), oName = sprintf('o%d%s', p, sfx); else, oName = sprintf('o%d,%d%s', p, drv, sfx); end
    W(end+1) = struct('port', p, 'excitation', tag, 'drivenPort', drv, 't', t, 'a', a, 'b', b, ...
        'iName', iName, 'oName', oName, 'Z', Zk, 'Zsource', zsrc); %#ok<AGROW>
end

% ---- optional CST-layout file -----------------------------------------------
if strlength(outFile) > 0
    od = fileparts(outFile);
    if ~isempty(od) && ~isfolder(od), mkdir(od); end
    L = readlines(vFile);
    parLine = L(find(startsWith(strtrim(L), "#Parameters"), 1));
    fid = fopen(outFile, 'w');
    cleaner = onCleanup(@() fclose(fid));
    [~, vn, ve] = fileparts(vFile);  [~, in, ie] = fileparts(iFile);
    for pass = 1:2                                          % all incident waves, then all outgoing
        for k = 1:numel(W)
            if pass == 1, nm = W(k).iName; x = W(k).a; else, nm = W(k).oName; x = W(k).b; end
            if ~isempty(parLine)                            % Z0 of this block = this port's impedance
                pl = regexprep(parLine, 'Z0\s*=\s*[^;}]*', sprintf('Z0=%.12g', W(k).Z));
                if ~contains(pl, 'Z0='), pl = regexprep(pl, '\}\s*$', sprintf('; Z0=%.12g}', W(k).Z)); end
                fprintf(fid, '%s\n', pl);
            else
                fprintf(fid, '#Parameters = {Z0=%.12g}\n', W(k).Z);
            end
            fprintf(fid, '#"Time / ns"\t"%s [Real / W^(1/2)]"\n', nm);
            fprintf(fid, '#---------------------------------\n');
            fprintf(fid, '%.13g\t%.13e\n', [W(k).t*1e9, x].');
            fprintf(fid, '\n');
        end
    end
    clear cleaner
    fprintf('Power waves (Z = %s ohm, from %s) from %s and %s written to\n  %s\n', ...
        strjoin("port " + ports + ": " + string(Zp(ports)), ', '), zsrc, ...
        [char(vn) char(ve)], [char(in) char(ie)], char(outFile));
end
end
