function [Q_TM, Q_TE] = ExtractLMmodes(Q, lmax)
% ExtractLMmodes  Rearrange spherical mode coefficients from the compressed
% single-index format into 2D (l, m) matrices, separately for TM and TE modes.
%
% The output matrices are indexed as:
%   Row    : l = 1..lmax  (polar order)
%   Column : m + lmax + 1 (m = -lmax -> col 1, m = 0 -> col lmax+1,
%                          m = +lmax -> col 2*lmax+1; see getmindx)
%
% Entries where |m| > l are zero (not physical).
%
% INPUTS:
%   Q     - Full complex mode coefficient array, size [num_j_indices x Nfreq],
%           as returned by ReadSphModes_freqs.
%   lmax  - Maximum polar mode order (nmax from ReadSphModes_freqs).
%
% OUTPUTS:
%   Q_TM  - TM (s=2) coefficients, size [lmax x (2*lmax+1) x Nfreq].
%   Q_TE  - TE (s=1) coefficients, size [lmax x (2*lmax+1) x Nfreq].
%
% INDEXING EXAMPLE (lmax=3):
%   Q_TM(l, m+4, f) is the TM coefficient of mode (l, m) at frequency f,
%   e.g. Q_TM(2, 2, f) -> l=2, m=-2;  Q_TM(2, 4, f) -> l=2, m=0.
%
% USAGE:
%   [Q, j_, j3D, mmax, nmax, freqs] = ReadSphModes_freqs('myantenna.sph');
%   [Q_TM, Q_TE] = ExtractLMmodes(Q, nmax);

Nfreq   = size(Q, 2);
ncols   = 2*lmax + 1;

Q_TM = zeros(lmax, ncols, Nfreq);
Q_TE = zeros(lmax, ncols, Nfreq);

for l = 1:lmax
    for m = -l:l
        col = m + lmax + 1;   % column index: m=-lmax -> 1, m=0 -> lmax+1

        jidx_TE = 2*(l*(l+1) + m - 1) + 1;   % s=1: TE
        jidx_TM = 2*(l*(l+1) + m - 1) + 2;   % s=2: TM

        Q_TE(l, col, :) = Q(jidx_TE, :);
        Q_TM(l, col, :) = Q(jidx_TM, :);
    end
end

end