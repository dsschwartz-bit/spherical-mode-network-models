function [Q_m0, j_m0, j3D_m0] = ExtractM0Modes(Q, nmax, j3D)
% ExtractM0Modes  Extract the m=0 spherical mode coefficients from the
% full Q array returned by ReadSphModes_freqs.
%
% For m=0 the compressed mode index reduces to:
%   j_ = 2*(n*(n+1) - 1) + s,   n = 1..nmax,  s = 1 (TE), 2 (TM)
%
% INPUTS:
%   Q     - Full complex mode coefficient array, size [num_j_indices x Nfreq],
%             as returned by ReadSphModes_freqs.
%   nmax  - Maximum elevation mode order, as returned by ReadSphModes_freqs.
%   j3D   - String label array indexed by compressed mode index, as returned
%             by ReadSphModes_freqs. (Optional - pass [] to skip.)
%
% OUTPUTS:
%   Q_m0  - Complex m=0 mode coefficients, size [2*nmax x Nfreq].
%             Rows are ordered as: (s=1,n=1), (s=2,n=1),
%                                  (s=1,n=2), (s=2,n=2), ...
%                                  (s=1,n=nmax), (s=2,n=nmax).
%             i.e. TE then TM for each n, stepping n from 1 to nmax.
%   j_m0  - Vector of compressed mode indices for the extracted modes,
%             length 2*nmax. (Named j_m0 to match the j_ convention.)
%   j3D_m0 - String labels for the extracted modes, length 2*nmax.
%             Empty if j3D input was [].
%
% USAGE:
%   [Q, j_, j3D, mmax, nmax, freqs] = ReadSphModes_freqs('myantenna.sph');
%   [Q_m0, j_m0, j3D_m0] = ExtractM0Modes(Q, nmax, j3D);

m = 0;
num_m0_modes = 2 * nmax;   % 2 modes (TE+TM) per n value
j_m0   = zeros(1, num_m0_modes);
Nfreq  = size(Q, 2);
Q_m0   = zeros(num_m0_modes, Nfreq);

row = 1;
for n = 1:nmax
    for s = 1:2   % s=1: TE,  s=2: TM
        jidx = 2*(n*(n+1) + m - 1) + s;   % compressed index for this (s,m=0,n)
        j_m0(row)    = jidx;
        Q_m0(row, :) = Q(jidx, :);
        row = row + 1;
    end
end

% Build label output
if ~isempty(j3D)
    j3D_m0 = j3D(j_m0);
else
    j3D_m0 = [];
end

end
