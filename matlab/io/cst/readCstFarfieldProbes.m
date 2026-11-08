function cst = readCstFarfieldProbes(file, components, excitation, length_unit)
%READCSTFARFIELDPROBES  Far-field probe time signals from a CST ASCII export.
%
%   cst = readCstFarfieldProbes(file, components, excitation, length_unit)
%
%   file        : CST farfield time-signal export, or a string array of
%                 exports (CST writes at most 25 curves per ASCII export,
%                 so many probes may be split over several files); the
%                 probes of all files are merged. Each block has a header
%                 like  #"Time / ns"  "E_Field (Farfield) (Spherical)
%                 (theta phi r)(Theta) [1] [Real / V/m]"
%   components  : string array of block names, e.g. ["Theta" "Phi"]
%   excitation  : excitation (port) index in the [n] tag to keep, e.g. 1.
%                 Co-simulation tags [cosim Tran1 n] use the trailing n.
%                 Blocks with another tag are skipped; blocks with no
%                 number in the tag (e.g. [cosim Tran1]) are kept.
%                 [] keeps every block.
%   length_unit : CST length unit in metres (1e-3 for mm), for r
%
%   cst : struct array, one entry per probe position, with fields
%         theta, phi (deg), r (m), t{c} (s) and e{c} = r*E (V) for each
%         component c (empty when that block is missing).
%   A block cut short at the end of the file is kept up to its last full
%   (time, value) pair; blocks much shorter than the rest give a warning.
%
%   CST accepts probes with theta > 180 (e.g. (225, 0)). These are mapped to
%   the standard range: (theta, phi) -> (360 - theta, phi + 180), and since
%   the theta and phi unit vectors of CST's (theta > 180, phi) point opposite
%   to those of the mapped direction, E_theta and E_phi change sign
%   (CST (225, 0) E_theta = -(135, 180) E_theta). A mapped probe that
%   duplicates one already read is skipped.

file = string(file);
nC = numel(components);
cst = struct('theta', {}, 'phi', {}, 'r', {}, 't', {}, 'e', {});
for fi = 1:numel(file)
C = strtrim(readlines(file(fi)));
C = C(strlength(C) > 0);
nHdr = nnz(startsWith(C, '#') & contains(C, '(Spherical)'));
if nHdr >= 25
    warning(['%s holds %d curves; CST ASCII exports stop at 25, so some probes ' ...
        'may be missing. Export fewer curves per file (e.g. only [%d] Theta/Phi) ' ...
        'and list all files in cst_file.'], file(fi), nHdr, excitation);
end
for c = 1:nC
    hdr = find(startsWith(C, '#') & contains(C, "(" + components(c) + ")"));
    for h = hdr.'
        % excitation tag: [1], or a trailing number after a space, [cosim Tran1 2] -> 2;
        % [cosim Tran1] (one excitation) counts as untagged and is kept
        tag = regexp(C(h), '\)\s*\[(?:[^\]]*\s)?(\d+)\s*\]', 'tokens', 'once');
        if ~isempty(excitation) && ~isempty(tag) && str2double(tag{1}) ~= excitation
            continue
        end
        pr = regexp(C(h), '\(Spherical\)\s*\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)', ...
            'tokens', 'once');
        if isempty(pr)
            error('No (Spherical) (theta phi r) probe position in the CST header:\n%s', C(h));
        end
        probe = str2double(pr);                        % [theta phi r]
        sgn = 1;
        th = mod(probe(1), 360);
        if th > 180                                    % map to theta in 0..180
            th = 360 - th;  probe(2) = probe(2) + 180;
            if any(strcmpi(components(c), ["Theta" "Phi"])), sgn = -1; end
        end
        probe(1) = th;  probe(2) = mod(probe(2), 360);
        tu = regexp(C(h), '"Time / (\w+)"', 'tokens', 'once');
        tScale = struct('s',1,'ms',1e-3,'us',1e-6,'ns',1e-9,'ps',1e-12).(tu{1});
        k = h + 1;  while k <= numel(C) && startsWith(C(k), '#'), k = k + 1; end
        if k > numel(C), continue, end
        e = k;      while e < numel(C) && ~startsWith(C(e+1), '#'), e = e + 1; end
        v = sscanf(strjoin(C(k:e), ' '), '%f');
        if mod(numel(v), 2)
            v = v(1:end-1);
        end
        M = reshape(v, 2, []).';
        rk = probe(3) * length_unit;
        n = find(abs([cst.theta] - probe(1)) < 1e-9 & abs([cst.phi] - probe(2)) < 1e-9, 1);
        if ~isempty(n) && ~isempty(cst(n).e{c})
            continue                                   % duplicate direction
        end
        if isempty(n)
            n = numel(cst) + 1;
            cst(n).theta = probe(1);  cst(n).phi = probe(2);  cst(n).r = rk;
            cst(n).t = cell(1, nC);   cst(n).e = cell(1, nC);
        end
        cst(n).t{c} = M(:, 1) * tScale;
        cst(n).e{c} = sgn * M(:, 2) * rk;              % r*E in V
    end
end
end
if isempty(cst)
    error('No probe blocks for excitation [%d] in %s', excitation, strjoin(file, ', '));
end
% flag blocks much shorter than the rest (an export cut short)
len = cellfun(@numel, [cst.t]);
nTyp = median(len(len > 0));
for n = 1:numel(cst)
    for c = 1:nC
        if ~isempty(cst(n).t{c}) && numel(cst(n).t{c}) < 0.9 * nTyp
            warning('CST (%s) block at theta = %g, phi = %g has %d of %d points (export cut short?).', ...
                components(c), cst(n).theta, cst(n).phi, numel(cst(n).t{c}), nTyp);
        end
    end
end
end
