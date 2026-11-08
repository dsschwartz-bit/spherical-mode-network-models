function result = compute_spherical_mode_coefficients_Q(fileprefix, lmax, fmin, fmax)
%COMPUTE_SPHERICAL_MODE_COEFFICIENTS_Q  TM/TE mode coefficients from a
%FEKO spherical-wave (.sph) export.
%   Reads the FEKO Q coefficients and converts them to the psi convention
%   (same Q-to-psi conversion as
%   matlab/workflows/general_vector_spherical_modes_FEKO.m).
%
%   result = COMPUTE_SPHERICAL_MODE_COEFFICIENTS_Q(fileprefix, lmax)
%   result = COMPUTE_SPHERICAL_MODE_COEFFICIENTS_Q(fileprefix, lmax, fmin, fmax)
%
%   fileprefix : string, path prefix of the FEKO output files, e.g.
%                fullfile(feko_dir, "cylinder_rotated_yaxis_45deg");
%                <prefix>_FarField1.sph and <prefix>.out must exist
%   lmax       : maximum spherical mode order (l = 1..lmax; there is no
%                l=0 TE/TM mode). Must not exceed the order stored in the
%                .sph file (error otherwise).
%   fmin, fmax : optional frequency range in Hz, applied after removing
%                points with negative real port impedance (Zindx). The
%                .sph frequencies must match the .out file exactly.
%                Omitted or empty: full range.
%
%   Returns a struct with fields:
%     QTMconv, QTEconv : (L x M x numf) complex coefficient arrays in the
%                    same psi convention and normalisation as psiTM/psiTE
%                    from the near-field and far-field functions in this
%                    folder. FEKO's Q (TICRA, exp(-i*omega*t)) is
%                    conjugated and scaled by conv_factor_TE/TM.
%                    Row lindx <-> l = lval(lindx) = lindx (l = 1..lmax).
%                    Column mindx <-> m = mval(mindx), mval = -lmax:lmax,
%                    fixed across all rows (m=0 always at column lmax+1).
%                    Entries with |m| > l are zero (the field-based
%                    functions use NaN instead).
%     lval, mval   : as above
%     L, M         : size(QTMconv,1), size(QTMconv,2)
%     Freq         : (1 x numf) frequency vector in Hz, after Zindx
%                    filtering and fmin/fmax
%     mmax, num_modes : from the .sph file (see ReadSphModes_freqs), for
%                    re-indexing into FEKO's raw Q vector via SphModeIndex
%     fileprefix   : echoed back
%
%   No ka de-embedding is applied; that belongs to the S-parameter step.

if nargin < 3 || isempty(fmin)
    fmin = -Inf;
end
if nargin < 4 || isempty(fmax)
    fmax = Inf;
end

fileprefix = string(fileprefix);

%% --- Physical Constants ---
c0   = physconst('lightspeed'); %#ok<NASGU>
mu0  = 4*pi*1e-7;
eps0 = 1/c0^2/mu0;
eta0 = sqrt(mu0/eps0);

%% --- FEKO Import: Read Port Parameters from .out File ---
[~, Z0, ~, symmetry_planes, ~, ~, Zant, Freq, ~] = ReadFEKOout(fileprefix+".out");

Zant = transpose(Zant(:,1)); Freq = transpose(Freq);
Zindx = find(~(real(Zant)<0));
Freq  = Freq(Zindx);

if abs(Z0-eta0) > 0.01
    error("antenna port impedance should be eta0 for now")
end

%% --- Mode index setup ---
lval = (1:1:lmax).';   % mode orders l = 1..lmax
mval = (-lmax:lmax);   % m ranges from -lmax to +lmax
L    = length(lval);
M    = length(mval);

use_jl1_convention = true;
if use_jl1_convention
    j_l_1 = (1j) .^ (lval + 1);   % (L x 1) j^(l+1) far-field phase normalisation
else
    j_l_1 = ones(L, 1);           % (L x 1) no phase normalisation
end

have_CSphase = false;
if have_CSphase
    CSphase = (-1).^mval;
else
    CSphase = ones(size(mval));   % (1 x M) no Condon-Shortley phase
end

%% --- Read Q Coefficients from .sph File ---
FF_request = "_FarField1.sph";
[Q, ~, ~, mmax, lmax_Q_file, freqs_Q] = ReadSphModes_freqs(fileprefix+FF_request);
freqs_Q = freqs_Q(Zindx);

% The Q file must contain all requested mode orders; otherwise the
% missing orders would silently come out as zero.
if lmax_Q_file < lmax
    error('compute_spherical_mode_coefficients_Q:lmaxExceedsFile', ...
          'Q file max mode order l=%d is less than the requested lmax=%d', lmax_Q_file, lmax)
end

if any(freqs_Q ~= Freq)
    error('compute_spherical_mode_coefficients_Q:freqMismatch', ...
          'frequencies from Q file are not the same as the .out file''s own impedance sweep')
end
Q = Q(:,Zindx);

%% --- Restrict to requested frequency range ---
frange_idx = find(Freq >= fmin & Freq <= fmax);
if isempty(frange_idx)
    error('compute_spherical_mode_coefficients_Q:noFreqInRange', ...
          'No frequency points found in [%.4g, %.4g] Hz (available range: [%.4g, %.4g] Hz)', ...
          fmin, fmax, min(Freq), max(Freq));
end
Freq = Freq(frange_idx);
Q    = Q(:,frange_idx);

numf = length(Freq);
num_modes = 2*(lmax*(lmax+1)+mmax-1)+2;

%% --- Convert Q (TICRA convention) to psi convention ---
% conj(Q) converts from the TICRA exp(-i*omega*t) convention to the
% exp(+j*omega*t) convention used throughout this toolbox.
QTEconv = zeros(L,M,numf); QTMconv = zeros(L,M,numf);

conv_factor_TE = 1j*sqrt(8*pi) .* (1j.*j_l_1) * 2^(symmetry_planes.physical_gnd);
conv_factor_TM =    sqrt(8*pi) .* (1j.*j_l_1) * 2^(symmetry_planes.physical_gnd);

for lindx=1:L
    for mindx=1:M
        if abs(mval(mindx))>lval(lindx)
            % |m| > l: mode does not exist, leave QTEconv/QTMconv as zero
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

%% --- Package result ---
result.QTMconv = QTMconv;
result.QTEconv = QTEconv;
result.lval  = lval;
result.mval  = mval;
result.L     = L;
result.M     = M;
result.Freq  = Freq;
result.mmax  = mmax;
result.num_modes = num_modes;
result.fileprefix = fileprefix;

end
