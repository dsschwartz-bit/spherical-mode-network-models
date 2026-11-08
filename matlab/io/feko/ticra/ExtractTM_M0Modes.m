function [Q_TM_m0, j_TM_m0, j3D_TM_m0] = ExtractTM_M0Modes(Q, nmax, j3D)
% ExtractTM_M0Modes  Extract the TM, m=0 spherical mode coefficients from
% the full Q array returned by ReadSphModes_freqs.
%
% TM modes have s=2. For m=0 the compressed mode index is:
%   j_ = 2*(n*(n+1) - 1) + 2,   n = 1..nmax
%
% INPUTS:
%   Q     - Full complex mode coefficient array, size [num_j_indices x Nfreq],
%             as returned by ReadSphModes_freqs.
%   nmax  - Maximum elevation mode order, as returned by ReadSphModes_freqs.
%   j3D   - String label array indexed by compressed mode index, as returned
%             by ReadSphModes_freqs. (Optional - pass [] to skip.)
%
% OUTPUTS:
%   Q_TM_m0  - Complex TM, m=0 mode coefficients, size [nmax x Nfreq].
%                Row n contains the coefficient for (s=2, m=0, n=n).
%   j_TM_m0  - Vector of compressed mode indices, length nmax.
%   j3D_TM_m0 - String labels '2 0 n' for each mode. Empty if j3D input was [].
%
% USAGE:
%   [Q, j_, j3D, mmax, nmax, freqs] = ReadSphModes_freqs('myantenna.sph');
%   [Q_TM_m0, j_TM_m0, j3D_TM_m0] = ExtractTM_M0Modes(Q, nmax, j3D);

m = 0;
s = 2;  % TM modes
Nfreq    = size(Q, 2);
Q_TM_m0  = zeros(nmax, Nfreq);
j_TM_m0  = zeros(1, nmax);

for n = 1:nmax
    jidx = 2*(n*(n+1) + m - 1) + s;   % compressed index for (s=2, m=0, n)
    j_TM_m0(n)    = jidx;
    Q_TM_m0(n, :) = Q(jidx, :);
end

% Build label output
if ~isempty(j3D)
    j3D_TM_m0 = j3D(j_TM_m0);
else
    j3D_TM_m0 = [];
end

end
