function q = compute_port_mode_coeffs_Q(fileprefix, far_name, far_occ, far_type, Freq, Zindx, lmax, ...
        lval, mval, L, M, symmetry_planes, j_l_1, CSphase, use_Q_method, port_label)
%COMPUTE_PORT_MODE_COEFFS_Q  Per-port TM/TE mode coefficients from a FEKO .sph export.
%   q = compute_port_mode_coeffs_Q(fileprefix, far_name, far_occ, far_type,
%       Freq, Zindx, lmax, lval, mval, L, M, symmetry_planes, j_l_1,
%       CSphase, use_Q_method, port_label)
%   Multi-port form of compute_mode_coeffs_Q: the .sph file of one internal
%   port is located with resolve_field_file.
%
%   FEKO's Q (TICRA, exp(-iwt)) are converted to the psi convention of
%   this toolbox (exp(+jwt)): conj(Q), scaled by
%       TM: sqrt(8 pi) * j * j_l_1 * 2^physical_gnd
%       TE: j * sqrt(8 pi) * j * j_l_1 * 2^physical_gnd
%   and the complex e^(jm phi) pairs (m, -m) combined into the real
%   cos/sin (tesseral) basis used by gradYlm, times CSphase(m).
%   No angular grid or projection is needed.
%
%   fileprefix   : string, path prefix of the FEKO output files
%   far_name, far_occ, far_type : the port's request name, occurrence
%                  count and reuse type (see resolve_field_file); if
%                  far_name is empty, q.available is false
%   Freq, Zindx  : (1 x numf) frequencies and the indices into the file's
%                  frequency list selected by the caller; the selected .sph
%                  frequencies must equal Freq
%   lmax         : maximum degree; an error (caught, method skipped) if the
%                  .sph file holds fewer degrees
%   lval, mval, L, M : mode-index grids built by the caller from lmax
%   symmetry_planes : struct with field .physical_gnd (from ReadFEKOout*);
%                  the coefficients are multiplied by 2^physical_gnd
%   j_l_1, CSphase : (L x 1) and (1 x M) convention factors from the caller
%   use_Q_method : if false, return immediately with q.available = false
%   port_label   : port name used in printed messages
%
% Returns struct q:
%   .available         - true if the .sph file was read
%   .QTMconv, .QTEconv - (L x M x numf) converted coefficients; zero for
%                        |m| > l; [] if not available
%   .Q, .freqs_Q, .mmax, .num_modes - raw Q (num_modes x numf) and its
%                        metadata, e.g. for radiated power from Q
%
% See also ReadSphModes_freqs, SphModeIndex.

numf = numel(Freq);
num_modes_fallback = 2*(lmax*(lmax+1)+lmax-1)+2;

q = struct('available', false, 'QTMconv', [], 'QTEconv', [], ...
    'Q', [], 'freqs_Q', Freq, 'mmax', lmax, 'num_modes', num_modes_fallback);

if ~use_Q_method || strlength(far_name) == 0
    return
end

try
    [sph_file, ~] = resolve_field_file(fileprefix, far_name, far_occ, far_type, "sph");
    if strlength(sph_file) == 0
        error('compute_port_mode_coeffs_Q:sph_not_found', 'no .sph file found for request "%s" (occurrence %d)', far_name, far_occ);
    end
    [Q, ~, ~, mmax, lmax_Q_file, freqs_Q] = ReadSphModes_freqs(sph_file);
    freqs_Q = freqs_Q(Zindx);

    if lmax_Q_file < max(lval)
        error("Q file max mode order l=%d is less than the script's lmax=%d", lmax_Q_file, max(lval))
    end

    if any(freqs_Q ~= Freq)
        error("frequencies from Q file not the same as field data")
    end
    Q = Q(:,Zindx);
    num_modes = 2*(lmax*(lmax+1)+mmax-1)+2;
catch ME
    fprintf('compute_port_mode_coeffs_Q (%s): could not read Q coefficients (%s) — skipping Q method for this port.\n', port_label, ME.message);
    return
end

% conj(Q): TICRA exp(-iwt) -> exp(+jwt).
QTEconv = zeros(L,M,numf); QTMconv = zeros(L,M,numf);

conv_factor_TE = 1j*sqrt(8*pi) .* (1j.*j_l_1) * 2^(symmetry_planes.physical_gnd);
conv_factor_TM =    sqrt(8*pi) .* (1j.*j_l_1) * 2^(symmetry_planes.physical_gnd);

for lindx=1:L
    for mindx=1:M
        if abs(mval(mindx))>lval(lindx)
            % |m| > l: no such mode, left at zero
            continue
        else
            if mval(mindx)==0
                QTMconv(lindx,mindx,:)=conj(Q(SphModeIndex(lval(lindx),mval(mindx),'TM'),:)) .* conv_factor_TM(lindx) ;
                QTEconv(lindx,mindx,:)=conj(Q(SphModeIndex(lval(lindx),mval(mindx),'TE'),:)) .* conv_factor_TE(lindx) ;
            elseif mval(mindx)<0
                QTMconv(lindx,mindx,:)=(CSphase(mindx))*(1j/sqrt(2))*...
                    conj( Q(SphModeIndex(lval(lindx),mval(mindx),'TM'),:)-(-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),-mval(mindx),'TM'),:) )...
                    .* conv_factor_TM(lindx) ;
                QTEconv(lindx,mindx,:)=(CSphase(mindx))*(1j/sqrt(2))*...
                    conj( Q(SphModeIndex(lval(lindx),mval(mindx),'TE'),:)-(-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),-mval(mindx),'TE'),:) )...
                    .* conv_factor_TE(lindx) ;
            elseif mval(mindx)>0
                QTMconv(lindx,mindx,:)=(CSphase(mindx))*(1/sqrt(2))*...
                    conj( (-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),mval(mindx),'TM'),:)+Q(SphModeIndex(lval(lindx),-mval(mindx),'TM'),:) )...
                    .* conv_factor_TM(lindx) ;
                QTEconv(lindx,mindx,:)=(CSphase(mindx))*(1/sqrt(2))*...
                    conj( (-1)^mval(mindx)*Q(SphModeIndex(lval(lindx),mval(mindx),'TE'),:)+Q(SphModeIndex(lval(lindx),-mval(mindx),'TE'),:) )...
                    .* conv_factor_TE(lindx) ;
            end
        end
    end
end

q.available = true;
q.QTMconv = QTMconv; q.QTEconv = QTEconv;
q.Q = Q; q.freqs_Q = freqs_Q; q.mmax = mmax; q.num_modes = num_modes;

end
