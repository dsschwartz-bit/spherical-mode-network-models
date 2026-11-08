function RotM = build_shell_rotation_matrices(R, max_l)
%BUILD_SHELL_ROTATION_MATRICES  Rotation matrices for real spherical
%harmonics, via direct recursion.
%
%   RotM = BUILD_SHELL_ROTATION_MATRICES(R, max_l) returns a
%   containers.Map, keyed by angular momentum l = 0 .. max_l, whose
%   values are the (2l+1)x(2l+1) rotation matrices for real spherical
%   harmonics of that angular momentum, given the 3x3 real-space
%   rotation matrix R.
%
%   Ported from OpenMolcas's Tools/mort/mort.py
%   (Copyright (C) 2021,2024-2026, Ignacio Fdez. Galvan; LGPL-2.1),
%   part of https://gitlab.com/Molcas/OpenMolcas.
%   Translated from Python to MATLAB by David Schwartz, 2026. This file is
%   distributed under the GNU LGPL 2.1; see COPYING.LGPL-2.1.txt in this folder.
%
%   Implements the recursive algorithm of:
%
%     J. Ivanic and K. Ruedenberg,
%     "Rotation Matrices for Real Spherical Harmonics.
%      Direct Determination by Recursion",
%     J. Phys. Chem. 100, 6342-6347 (1996). doi:10.1021/jp953350u
%
%     J. Ivanic and K. Ruedenberg,
%     "Rotation Matrices for Real Spherical Harmonics. Direct
%      Determination by Recursion (Additions and Corrections)",
%     J. Phys. Chem. A 102, 9099-9100 (1998). doi:10.1021/jp9833350
%
%   The recursion (Eq. 8.1 of the 1996 paper) is:
%
%       R^l_{m,m'} = u * U^l_{m,m'} + v * V^l_{m,m'} + w * W^l_{m,m'}
%
%   where u, v, w are numerical coefficients (Table 1) and U, V, W are
%   built from R^{l-1} and R^1 via the functions iP_{mu,m'} (Table 2).
%
%   Two additional corrections, beyond the official 1998 erratum, are
%   applied here as identified by the OpenMolcas developers (see NB
%   comments inside SPHERICAL_SHELL_ROTATION below):
%     1. In the m' = +-l cases, the original formulas had m' -> l-1
%        where they should not.
%     2. In the correction for m < 0, a factor sqrt(1-delta) should
%        be sqrt(1+delta).
%
%   Example:
%       R = eye(3);                       % identity rotation
%       RotM = build_shell_rotation_matrices(R, 4);
%       Rd = RotM(2);                      % 5x5 rotation matrix for d shell

R = normalize_rot(R, true);

RotM = containers.Map('KeyType', 'double', 'ValueType', 'any');

% Seed with the trivial Cartesian s and p shells
RotM(0) = cartesian_shell_rotation(0, R);
RotM(1) = cartesian_shell_rotation(1, R);
% Reorder the p rotation to match the m ordering (+1, 0, -1) = (x,y,z)
R1 = RotM(1);
RotM(1) = R1([2, 3, 1], [2, 3, 1]);

% Recursively build l = 2 .. max_l
for l = 2:max_l
    RotM(l) = spherical_shell_rotation(l, RotM); %#ok<NASGU>
end

end


% ------------------------------------------------------------------
% l = 0, 1 base cases (Cartesian shells; trivial but built generically)
% ------------------------------------------------------------------

function Rot = normalize_rot(A, test)
%NORMALIZE_ROT  Find the closest unitary matrix to A (cleans up
%numerical noise).
if nargin < 2
    test = true;
end
[U, S, V] = svd(A);
s = diag(S);
if test
    if ~all(abs(s - 1.0) < 1e-3)
        error('build_shell_rotation_matrices:notUnitary', ...
              'Not a unitary rotation');
    end
end
Rot = U * V';
end


function perms = permute_unique(nums)
%PERMUTE_UNIQUE  Unique permutations (permutations with repeated
%elements), returned as a cell array of row vectors.
perms = {zeros(1, 0)};
for k = 1:length(nums)
    n = nums(k);
    new_perm = {};
    for p = 1:length(perms)
        perm = perms{p};
        for i = 0:length(perm)
            new_perm{end+1} = [perm(1:i), n, perm(i+1:end)]; %#ok<AGROW>
            if i < length(perm) && perm(i+1) == n
                break;
            end
        end
    end
    perms = new_perm;
end
end


function f = factors(l, m)
%FACTORS  Cartesian factors (0=x, 1=y, 2=z) for angular momentum l,
%shell index m (m need not be a valid magnetic quantum number; it is
%used purely as a linear index over the shell's functions).
ly = floor((sqrt(8 * (m + l) + 1) - 1) / 2);
lz = m + l - ly * (ly + 1) / 2;
lx = l - ly;
ly = ly - lz;
f = [zeros(1, lx), ones(1, ly), 2 * ones(1, lz)];
end


function Rot = cartesian_shell_rotation(l, R)
%CARTESIAN_SHELL_ROTATION  Rotation matrix for the functions in a
%Cartesian shell of angular momentum l, given the real-space 3x3
%rotation matrix R. Used here only for l = 0, 1 to seed the recursion.
if l < 0
    error('build_shell_rotation_matrices:badL', 'l must be >= 0');
end
num = round((l + 1) * (l + 2) / 2);
Rot = zeros(num, num);

f = cell(1, num);
idx = 1;
for m = -l:(num - l - 1)
    f{idx} = factors(l, m);
    idx = idx + 1;
end

for m1 = 1:num
    perms = permute_unique(f{m1});
    for m2 = 1:num
        for pidx = 1:length(perms)
            p = perms{pidx};
            K = 1.0;
            for i = 1:length(p)
                % p(i), f{m2}(i) are 0-based Cartesian indices (0=x,1=y,2=z)
                K = K * R(p(i) + 1, f{m2}(i) + 1);
            end
            Rot(m2, m1) = Rot(m2, m1) + K;
        end
    end
end
end


% ------------------------------------------------------------------
% l > 1: the Ivanic-Ruedenberg recursion, with corrections
% ------------------------------------------------------------------

function Rot = spherical_shell_rotation(l, RotM)
%SPHERICAL_SHELL_ROTATION  Rotation matrix for the functions in a real
%spherical harmonics shell of angular momentum l (l > 1), built
%recursively from the l=1 rotation matrix RotM(1) and the (l-1)
%rotation matrix RotM(l-1).
%
%   RotM: containers.Map, l -> rotation matrix; must already contain
%         RotM(1) (reordered so column/row 1 <-> m=+1, 2 <-> m=0,
%         3 <-> m=-1) and RotM(l-1).
%
%   Note: this builds the inverse/transpose, so indexing may look
%   reversed relative to the paper's R_{m,m'} convention.
if l <= 1
    error('build_shell_rotation_matrices:badL', 'l must be > 1');
end

num = 2 * l + 1;
Rot = zeros(num, num);

for m2 = -l:l
    if abs(m2) == l
        den = 2 * l * (2 * l - 1);
    else
        den = (l + m2) * (l - m2);
    end
    for m1 = -l:l
        u = sqrt((l + m1) * (l - m1) / den);
        rU = rP(0, m1, m2);
        if m1 == 0
            v = -0.5 * sqrt(2 * (l - 1) * l / den);
            rV = rP(1, 1, m2) + rP(-1, -1, m2);
            w = 0;
            rW = 0;
        else
            if m1 > 0
                s = 1;
            else
                s = -1;
            end
            absm = abs(m1);
            v = 0.5 * sqrt((l + absm - 1) * (l + absm) / den);
            % NB: The case m < 0 is wrong in the 1998 correction
            % (sqrt(1-delta) -> sqrt(1+delta)); corrected here.
            if absm == 1
                rV = sqrt(2) * rP(s, absm - 1, m2);
            else
                rV = rP(1, m1 - s, m2) - s * rP(-1, -m1 + s, m2);
            end
            w = -0.5 * sqrt((l - absm - 1) * (l - absm) / den);
            rW = rP(1, m1 + s, m2) + s * rP(-1, -m1 - s, m2);
        end
        Rot(l + m2 + 1, l + m1 + 1) = u * rU + v * rV + w * rW;
    end
end

    function val = rP(i, mu, m2)
        % NB: The m2 == +-l cases are wrong in the original 1996 paper
        % (m' -> l-1); corrected here.
        if abs(mu) > l - 1
            val = 0;
            return;
        end
        R1 = RotM(1);
        Rlm1 = RotM(l - 1);
        if m2 == l
            val = R1(3, i + 2) * Rlm1(2*l - 1, l + mu) ...
                - R1(1, i + 2) * Rlm1(1, l + mu);
        elseif m2 == -l
            val = R1(3, i + 2) * Rlm1(1, l + mu) ...
                + R1(1, i + 2) * Rlm1(2*l - 1, l + mu);
        else
            val = R1(2, i + 2) * Rlm1(l + m2, l + mu);
        end
    end

end
