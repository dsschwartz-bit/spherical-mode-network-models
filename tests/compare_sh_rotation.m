%COMPARE_SH_ROTATION  Cross-check two real-SH rotation implementations.
%   Compares two independent implementations of the Ivanic-Ruedenberg
%   real spherical harmonic rotation recursion:
%
%   (A) build_shell_rotation_matrices.m  - ported from OpenMolcas's
%       Tools/mort/mort.py (quantum chemistry)
%   (B) SHRotate.m + rotate_coeff.m      - independent toolbox
%       (computer graphics / SH lighting lineage)
%
% Both are applied to the same set of real spherical harmonic
% coefficients under the same 3x3 rotation matrix (built from an
% axis+angle via rotationmat3D.m). The rotated mode coefficients from
% each method are printed side by side, and the resulting function is
% plotted on the unit sphere (r=1) over a THETA/PHI meshgrid, before
% and after rotation.
%
% No external data needed (all functions are in matlab/ and third_party/).

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

%% ---------------- Input: SH coefficient matrix ----------------
% coeff_matrix is (L+1) x (2L+1), matching the psiTM/psiTE layout in the
% workflow scripts (matlab/workflows):
%   row l+1  -> angular momentum l  (l = 0 .. L)
%   col j    -> magnetic number m = j - L - 1, FIXED for every row
%               (columns span the full -L..L range regardless of that
%               row's own l), so m=0 always sits at the same column
%               (L+1) for every row.
% For a given row l, only the middle (2l+1) columns -- m = -l .. l --
% are valid; entries outside that range (|m| > l) are structurally
% unused and must be 0 or NaN (NaN is treated as 0 below, matching the
% nan(...) preallocation of psiTM/psiTE).
%
% Three ways to provide input -- use ONE of these:
%
% OPTION A: a single pure mode Y_l0^m0. Just change l0/m0 -- the
% matrix is sized directly from them, nothing else to keep in sync.
l0 = 3; m0 = 2;
if abs(m0) > l0
    error('m0 (%d) is not a valid mode for l0 (%d): need |m0| <= l0', m0, l0);
end
coeff_matrix = zeros(l0+1, 2*l0+1);
coeff_matrix(l0+1, l0+1+m0) = 1;

% OPTION B: supply your own full matrix instead (any (L+1)x(2L+1)
% matrix in the layout described above -- L is inferred automatically
% from its size, further down). Comment out Option A above and use:
% coeff_matrix = my_data;

% OPTION C: an arbitrary superposition across bands, e.g.:
% Lmax = 4; rand('seed', 1);
% coeff_matrix = zeros(Lmax+1, 2*Lmax+1);
% for l = 0:Lmax
%     col_range = (Lmax+1-l):(Lmax+1+l);   % m = -l .. l
%     coeff_matrix(l+1, col_range) = randn(1, 2*l+1);
% end

rot_angle = deg2rad(70);     % rotation angle
rot_axis  = [1, 1, 0.3];     % rotation axis (need not be pre-normalized)

grid_n_theta = 80;           % meshgrid resolution (colatitude)
grid_n_phi   = 160;          % meshgrid resolution (longitude)

%% ---------------- Sanitize and validate ----------------
coeff_matrix(isnan(coeff_matrix)) = 0;   % NaN (invalid slot) treated as absent

% lmax_input is derived from the size of coeff_matrix.
lmax_input = size(coeff_matrix, 1) - 1;
if size(coeff_matrix, 2) ~= 2*lmax_input + 1
    error('coeff_matrix must be (L+1) x (2L+1): got %d rows but %d columns (expected %d)', ...
          size(coeff_matrix,1), size(coeff_matrix,2), 2*lmax_input+1);
end

% Reject nonzero values in invalid (|m| > l) slots.
for l = 0:lmax_input
    invalid_cols = [1:(lmax_input-l), (lmax_input+l+2):(2*lmax_input+1)];
    if any(coeff_matrix(l+1, invalid_cols) ~= 0)
        error('coeff_matrix row l=%d has a nonzero entry outside |m|<=%d', l, l);
    end
end

%% ---------------- Determine lmax from the data itself ----------------
% Rotation is block-diagonal in l (rotating Y_l^m only ever produces a
% superposition of Y_l^m' within the SAME l), so there is never a
% reason to build rotation matrices or evaluate bands above the
% highest l that actually has a nonzero coefficient. Trim to that.
active_rows = find(any(coeff_matrix ~= 0, 2));
if isempty(active_rows)
    lmax = 0;
else
    lmax = max(active_rows) - 1;
end
% Re-center: keep the SAME global column alignment (m=0 at column
% lmax_input+1) while dropping unused higher-l rows.
coeff_matrix = coeff_matrix(1:lmax+1, :);

fprintf('Input populated up to l = %d (of %d allocated) -- using lmax = %d\n\n', ...
        lmax, lmax_input, lmax);

%% ---------------- Convert matrix -> flat vector ----------------
% Flat layout matches real_spherical_harmonics.m / rotate_coeff.m:
% blocks of size (2l+1) for l = 0,1,2,..., each ordered m = -l .. l,
% i.e. global 1-based index = l^2 + (m + l) + 1. Note this flat layout
% is left-aligned per block (unlike coeff_matrix, which stays centered
% on a fixed global m=0 column) -- the conversion below re-indexes
% between the two.
ncoef = (lmax + 1)^2;
coeff = zeros(1, ncoef);
for l = 0:lmax
    col_range = (lmax_input+1-l) : (lmax_input+1+l);   % m = -l .. l in coeff_matrix
    coeff(l^2+1 : l^2+2*l+1) = coeff_matrix(l+1, col_range);
end

%% ---------------- Build rotation matrix ----------------
R = rotationmat3D(rot_angle, rot_axis);

fprintf('Rotation axis (normalized) = [%.3f %.3f %.3f], angle = %.2f deg\n', ...
        rot_axis / norm(rot_axis), rad2deg(rot_angle));
fprintf('det(R) = %.10f   max|R*R^T - I| = %.3e\n\n', ...
        det(R), max(max(abs(R * R' - eye(3)))));

%% ---------------- Rotate coefficients: two independent methods ----------------
% Method A: build_shell_rotation_matrices.m
RotM = build_shell_rotation_matrices(R, lmax);
dest_ours = zeros(size(coeff));
for l = 0:lmax
    sp = l^2;
    dest_ours(sp+1:sp+2*l+1) = coeff(sp+1:sp+2*l+1) * RotM(l);
end

% Method B: SHRotate.m / rotate_coeff.m (independent toolbox)
dest_toolbox = rotate_coeff(coeff, R, lmax);

% Rotated results in the centred matrix layout (row = l+1, col = m+lmax+1)
dest_ours_matrix    = zeros(lmax+1, 2*lmax+1);
dest_toolbox_matrix = zeros(lmax+1, 2*lmax+1);
for l = 0:lmax
    sp = l^2;
    col_range = (lmax+1-l) : (lmax+1+l);
    dest_ours_matrix(l+1, col_range)    = dest_ours(sp+1:sp+2*l+1);
    dest_toolbox_matrix(l+1, col_range) = dest_toolbox(sp+1:sp+2*l+1);
end

%% ---------------- Print mode coefficients ----------------
fprintf('%-4s %-4s | %-12s | %-14s %-14s %-10s\n', ...
        'l', 'm', 'orig', 'rot (ours)', 'rot (toolbox)', 'diff');
fprintf('%s\n', repmat('-', 1, 64));
k = 1;
for l = 0:lmax
    for m = -l:l
        fprintf('%-4d %-4d | %-12.6f | %-14.6f %-14.6f %-10.2e\n', ...
                l, m, coeff(k), dest_ours(k), dest_toolbox(k), ...
                abs(dest_ours(k) - dest_toolbox(k)));
        k = k + 1;
    end
end
fprintf('\nMax |ours - toolbox| over all coefficients: %.3e\n\n', ...
        max(abs(dest_ours - dest_toolbox)));

%% ---------------- Evaluate on a THETA/PHI grid, r = 1 ----------------
theta = linspace(0, pi, grid_n_theta);
phi   = linspace(0, 2*pi, grid_n_phi);
[PHI, THETA] = meshgrid(phi, theta);

F_orig    = zeros(size(THETA));
F_ours    = zeros(size(THETA));
F_toolbox = zeros(size(THETA));

for i = 1:numel(THETA)
    F_orig(i)    = real_spherical_harmonics([THETA(i), PHI(i)], coeff,        lmax);
    F_ours(i)    = real_spherical_harmonics([THETA(i), PHI(i)], dest_ours,    lmax);
    F_toolbox(i) = real_spherical_harmonics([THETA(i), PHI(i)], dest_toolbox, lmax);
end

r = 1;
X = r * sin(THETA) .* cos(PHI);
Y = r * sin(THETA) .* sin(PHI);
Z = r * cos(THETA);

%% ---------------- Plot ----------------
clim = max(abs([F_orig(:); F_ours(:); F_toolbox(:)])) * [-1 1];
if clim(2) == 0
    clim = [-1 1];
end
n = 128;
r_ = [linspace(0.10, 1, n), linspace(1, 0.70, n)];
g_ = [linspace(0.30, 1, n), linspace(1, 0.10, n)];
b_ = [linspace(0.75, 1, n), linspace(1, 0.10, n)];
cmap = [r_', g_', b_'];

fig = figure('Color', 'w', 'Position', [100 100 1300 900]);

ax1 = subplot(2,2,1);
surf(X, Y, Z, F_orig, 'EdgeColor', 'none'); axis equal off; view(135,20);
title(ax1, 'Original function');
colormap(ax1, cmap); caxis(ax1, clim); colorbar(ax1);

ax2 = subplot(2,2,2);
surf(X, Y, Z, F_ours, 'EdgeColor', 'none'); axis equal off; view(135,20);
title(ax2, 'Rotated: build\_shell\_rotation\_matrices');
colormap(ax2, cmap); caxis(ax2, clim); colorbar(ax2);

ax3 = subplot(2,2,3);
surf(X, Y, Z, F_toolbox, 'EdgeColor', 'none'); axis equal off; view(135,20);
title(ax3, 'Rotated: SHRotate / rotate\_coeff');
colormap(ax3, cmap); caxis(ax3, clim); colorbar(ax3);

ax4 = subplot(2,2,4);
Diff = F_ours - F_toolbox;
surf(X, Y, Z, Diff, 'EdgeColor', 'none'); axis equal off; view(135,20);
title(ax4, sprintf('Difference (max = %.1e)', max(abs(Diff(:)))));
colormap(ax4, cmap); colorbar(ax4);

annotation('textbox', [0 0.95 1 0.05], 'String', ...
    sprintf('SH coefficients up to l = %d, rotated %.0f deg about axis [%.2f %.2f %.2f]', ...
            lmax, rad2deg(rot_angle), rot_axis / norm(rot_axis)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

%% ---------------- Per-band decomposition ----------------
% Rotation is block-diagonal in l, so it only redistributes amplitude
% among Y_l^m of the same l. For each non-empty band, evaluate each
% weighted rotated term, check that their sum equals a direct evaluation
% of that band, then check that the sum over all bands reproduces F_ours.

fprintf('Per-band decomposition of the rotated function:\n');
F_band_sum = zeros(size(THETA));   % accumulates the full reconstruction across all bands

for l = 0:lmax
    sp = l^2;
    nterms = 2*l + 1;
    orig_block = coeff(sp+1:sp+nterms);
    rot_block  = dest_ours(sp+1:sp+nterms);
    if all(orig_block == 0) && all(rot_block == 0)
        continue;   % nothing in this band, skip it
    end

    fprintf('\n l = %d:\n', l);
    for jj = 1:nterms
        mm = jj - l - 1;
        fprintf('   m=%3d : orig=%+9.6f   rotated=%+9.6f\n', mm, orig_block(jj), rot_block(jj));
    end

    % Evaluate each individual weighted rotated term on the grid
    F_terms = zeros(nterms, numel(THETA));
    for jj = 1:nterms
        coeff_term = zeros(1, ncoef);
        coeff_term(sp+jj) = rot_block(jj);
        for i = 1:numel(THETA)
            F_terms(jj,i) = real_spherical_harmonics([THETA(i), PHI(i)], coeff_term, lmax);
        end
    end
    F_band = reshape(sum(F_terms, 1), size(THETA));
    F_band_sum = F_band_sum + F_band;

    % sanity check: evaluating this block alone (via real_spherical_harmonics'
    % own internal summation) must agree with our manual term-by-term sum
    coeff_block_only = zeros(1, ncoef);
    coeff_block_only(sp+1:sp+nterms) = rot_block;
    F_band_direct = zeros(size(THETA));
    for i = 1:numel(THETA)
        F_band_direct(i) = real_spherical_harmonics([THETA(i), PHI(i)], coeff_block_only, lmax);
    end
    fprintf('   max |band superposition - direct eval of this block| = %.3e\n', ...
            max(abs(F_band(:) - F_band_direct(:))));

    % Plot this band: each weighted term, tiled in a
    % ceil(sqrt(2l+1)) x ceil(sqrt(2l+1)) grid (no band-sum panel)
    ntile = ceil(sqrt(nterms));
    figure('Color', 'w', 'Position', [50 50 220*ntile 220*ntile]);
    termclim = max(abs(F_terms(:))) * [-1 1];
    if termclim(2) == 0
        termclim = [-1 1];
    end
    for jj = 1:nterms
        mm = jj - l - 1;
        axj = subplot(ntile, ntile, jj);
        Fj = reshape(F_terms(jj,:), size(THETA));
        surf(X, Y, Z, Fj, 'EdgeColor', 'none'); axis equal off; view(135,20);
        title(axj, sprintf('l=%d m=%d\ncoeff=%.3f', l, mm, rot_block(jj)));
        colormap(axj, cmap); caxis(axj, termclim);
    end

    annotation('textbox', [0 0.95 1 0.05], 'String', ...
        sprintf('Band l=%d: rotated superposition of Y_{%d,m} terms', l, l), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);
end

%% ---------------- Full reconstruction across all bands ----------------
recon_err = max(abs(F_band_sum(:) - F_ours(:)));
fprintf('\nMax |sum of ALL band contributions - full rotated evaluation| = %.3e\n', recon_err);

figure('Color', 'w', 'Position', [100 100 950 320]);

axr1 = subplot(1,3,1);
surf(X, Y, Z, F_band_sum, 'EdgeColor', 'none'); axis equal off; view(135,20);
title(axr1, 'Reconstructed: sum of all bands');
colormap(axr1, cmap); caxis(axr1, clim); colorbar(axr1);

axr2 = subplot(1,3,2);
surf(X, Y, Z, F_ours, 'EdgeColor', 'none'); axis equal off; view(135,20);
title(axr2, 'Direct evaluation (F\_ours)');
colormap(axr2, cmap); caxis(axr2, clim); colorbar(axr2);

axr3 = subplot(1,3,3);
Drec = F_band_sum - F_ours;
surf(X, Y, Z, Drec, 'EdgeColor', 'none'); axis equal off; view(135,20);
title(axr3, sprintf('Difference (max = %.1e)', max(abs(Drec(:)))));
colormap(axr3, cmap); colorbar(axr3);

annotation('textbox', [0 0.92 1 0.06], 'String', ...
    'Full reconstruction: summing every band''s rotated superposition reproduces the complete rotated function', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);