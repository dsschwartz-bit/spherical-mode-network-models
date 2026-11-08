function plot_stem_type(type_name, dom_lm, lval, mval, psi_rot, psi_via_matrix, f0_idx, f0_actual, rot_angle)
%PLOT_STEM_TYPE  Stem plot of the dominant TM or TE modes at one frequency.
%   Magnitude (dB) and phase (deg) of the direct and rotation-matrix
%   coefficients, at frequency index f0_idx, on a compact single-index
%   x-axis (l = 1..lmax, m = -l..l within each l). If dom_lm is empty a
%   note is printed and no figure is created.
%   (Project helper for tests/compare_rotated_antenna_coefficients.m; not
%   part of the third-party real_spherical_harmonics package.)
%
%   plot_stem_type(type_name, dom_lm, lval, mval, psi_rot, ...
%                   psi_via_matrix, f0_idx, f0_actual, rot_angle)
%
%   type_name       : 'TM' or 'TE' (used only for labeling)
%   dom_lm          : Nx2 array of [lindx, mindx] pairs identified as
%                     dominant for THIS type (e.g. dom_lindx_mindx_TM)
%   lval, mval      : mode index vectors, as returned by
%                     compute_spherical_mode_coefficients
%   psi_rot         : (L x M x numf) direct (rotated-model) coefficients
%   psi_via_matrix  : (L x M x numf) rotation-matrix-applied coefficients
%   f0_idx          : frequency index (into psi_rot's 3rd dim) to
%                     evaluate every mode at
%   f0_actual       : the actual frequency (Hz) at f0_idx, for the title
%   rot_angle       : rotation angle in radians, for the title

if isempty(dom_lm)
    fprintf('  %s: no dominant modes -- skipping %s stem plot\n', type_name, type_name);
    return
end

lmax = max(lval);

% Sort the dominant (lindx,mindx) pairs into l-then-m order for the x-axis
[~, sort_idx] = sortrows(dom_lm, [1 2]);
dom_lm = dom_lm(sort_idx, :);

n_modes = size(dom_lm, 1);
mode_l = zeros(n_modes,1); mode_m = zeros(n_modes,1);
direct_db = zeros(n_modes,1); matrix_db = zeros(n_modes,1);
direct_ph = zeros(n_modes,1); matrix_ph = zeros(n_modes,1);

for kk = 1:n_modes
    lindx = dom_lm(kk,1); mindx = dom_lm(kk,2);
    mode_l(kk) = lval(lindx); mode_m(kk) = mval(mindx);

    vd = psi_rot(lindx, mindx, f0_idx);
    vm = psi_via_matrix(lindx, mindx, f0_idx);
    direct_db(kk) = 20*log10(max(abs(vd), realmin));
    matrix_db(kk) = 20*log10(max(abs(vm), realmin));
    direct_ph(kk) = rad2deg(angle(vd));
    matrix_ph(kk) = rad2deg(angle(vm));
end

idx = (1:n_modes).';
l_starts = [1; find(diff(mode_l) ~= 0) + 1];
l_labels = mode_l(l_starts);

figure('Color','w','Position',[100 100 max(500, 80*n_modes+200) 700]);

ax1 = subplot(2,1,1);
stem(idx, direct_db, 'b', 'filled', 'MarkerSize', 5); hold on
stem(idx, matrix_db, 'r', 'MarkerSize', 5);
hold off
ylabel(sprintf('$|\\psi_{%s}|$ (dB)', type_name));
title(sprintf('%s dominant modes, magnitude', type_name));
legend('direct (rotated model)','via rotation matrix','Location','best','Interpreter','latex');
grid on

ax2 = subplot(2,1,2);
stem(idx, direct_ph, 'b', 'filled', 'MarkerSize', 5); hold on
stem(idx, matrix_ph, 'r', 'MarkerSize', 5);
hold off
ylabel(sprintf('$\\angle\\psi_{%s}$ (deg)', type_name));
title(sprintf('%s dominant modes, phase', type_name));
ylim([-180 180]); grid on

for ax = [ax1, ax2]
    axes(ax); %#ok<LAXES>
    xlabel(sprintf('mode index (single-indexing: $l=1..%d$, $m=-l..l$ within each $l$)', lmax));
    xlim([0.5, n_modes+0.5]);
    yl = ylim;
    for bb = 1:length(l_starts)
        line([l_starts(bb)-0.5, l_starts(bb)-0.5], yl, 'Color', [0.8 0.8 0.8], 'HandleVisibility','off');
        text(l_starts(bb), yl(2), sprintf('$l=%d$',l_labels(bb)), ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', ...
            'FontSize', 8, 'Color', [0.4 0.4 0.4], 'Interpreter', 'latex');
    end
end

annotation('textbox',[0 0.96 1 0.04], 'String', ...
    sprintf('%s dominant modes at $f_0 = %.4g$ Hz: direct vs via rotation matrix ($%.0f^\\circ$ about $y$-axis)', ...
            type_name, f0_actual, rad2deg(rot_angle)), ...
    'EdgeColor','none', 'HorizontalAlignment','center', 'FontSize', 12, 'Interpreter', 'latex');

end
