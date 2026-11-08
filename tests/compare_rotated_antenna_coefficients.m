%COMPARE_ROTATED_ANTENNA_COEFFICIENTS  Validate SHRotate on FEKO data.
%   Compares two routes to the mode coefficients of an antenna rotated by
%   45 deg about the y-axis:
%     direct     : psiTM/psiTE from the FEKO model that is already rotated
%                  (cylinder_rotated_yaxis_45deg);
%     via matrix : psiTM/psiTE from the unrotated FEKO model
%                  (dipole_cylinder_edgeport_symmetry_full), rotated with
%                  the real spherical-harmonic rotation matrices (SHRotate).
%   For each l, TM and TE coefficients each transform with the same
%   (2l+1)x(2l+1) block as the real harmonics Y_lm, and never mix across l.
%   The two routes should agree if both models differ only by the rotation
%   and the harmonicY real basis convention matches SHRotate's
%   (Ivanic-Ruedenberg eq. 2.3-2.5); a clean per-m sign flip or per-l scale
%   factor points to a convention mismatch (check l = 1 first).
%
%   Inputs : FEKO near-field exports (<prefix>_NearField1.efe/.hfe and
%            <prefix>.out) for both models, in base_dir (USER SETTINGS).
%   Outputs: table of dominant modes, relative L2 errors, and magnitude,
%            phase and stem plots.

clear; clc; close all;

%% ===== USER SETTINGS =====
base_dir = "C:\path\to\FEKO\dipole";   % unrotated model here; rotated model in base_dir\rotated
% =========================

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(here);                           % compute_spherical_mode_coefficients_* helpers

% Dock all figures in one tabbed window and default all text to the
% LaTeX interpreter.
set(0, 'DefaultFigureWindowStyle', 'docked');
try, set(0, 'DefaultTextInterpreter', 'latex'); catch, end
try, set(0, 'DefaultLegendInterpreter', 'latex'); catch, end
try, set(0, 'DefaultAxesTickLabelInterpreter', 'latex'); catch, end

%% ---------------- Configuration ----------------
lmax = 9;   % maximum spherical mode order

fmin = [];   % Hz, lower frequency bound; leave empty [] to use the full range in the files
fmax = [];   % Hz, upper frequency bound; leave empty [] to use the full range in the files

fileprefix_unrot = fullfile(base_dir, "dipole_cylinder_edgeport_symmetry_full");
fileprefix_rot   = fullfile(base_dir, "rotated", "cylinder_rotated_yaxis_45deg");

rot_axis  = [0, 1, 0];     % y-axis
rot_angle = deg2rad(45);   % 45 degrees, matching the antenna's actual rotation

n_dom = 5;   % number of dominant (l,m,type) modes to report/plot in detail

f0 = 4e9;   % Hz, frequency of interest for the stem plots (nearest available
            % point in Freq_common is used -- see "Stem plot" section below)

%% ---------------- Compute coefficients for both models ----------------
fprintf('Computing spherical mode coefficients: unrotated model...\n');
res_unrot = compute_spherical_mode_coefficients_nearfield(fileprefix_unrot, lmax, fmin, fmax);

fprintf("Computing spherical mode coefficients: rotated model "+ num2str(rad2deg(rot_angle)) + "deg about ["+ num2str(rot_axis)+ "] ...\n");
res_rot = compute_spherical_mode_coefficients_nearfield(fileprefix_rot, lmax, fmin, fmax);

lval = res_unrot.lval; mval = res_unrot.mval; L = res_unrot.L; M = res_unrot.M;
if ~isequal(lval, res_rot.lval) || ~isequal(mval, res_rot.mval)
    error('lval/mval mismatch between the two models -- were they computed with the same lmax?');
end

%% ---------------- Interpolate rotated model onto the unrotated model's frequency grid ----------------
% The two FEKO sweeps may not have been solved at identical frequency
% points, so rather than requiring an exact match, interpolate the
% rotated model's (complex) psiTM/psiTE onto the unrotated model's
% Freq grid. interp1 interpolates real and imaginary parts directly,
% which is the right behavior for smoothly-varying frequency-domain data.
Freq_common = res_unrot.Freq;
numf = length(Freq_common);

if min(Freq_common) < min(res_rot.Freq) || max(Freq_common) > max(res_rot.Freq)
    warning('compare_rotated_antenna_coefficients:extrapolation', ...
        ['Unrotated model frequency grid extends beyond the rotated model''s ' ...
         'range [%.4g, %.4g] Hz -- some points will be extrapolated (NaN).'], ...
        min(res_rot.Freq), max(res_rot.Freq));
end

psiTM_unrot = res_unrot.psiTM; psiTE_unrot = res_unrot.psiTE;

% reshape (L x M x numf_rot) -> (numf_rot x L*M) for a single vectorized
% interp1 call across all (l,m) at once, then reshape back
psiTM_rot_2d = reshape(res_rot.psiTM, L*M, []).';
psiTE_rot_2d = reshape(res_rot.psiTE, L*M, []).';
psiTM_rot_interp_2d = interp1(res_rot.Freq, psiTM_rot_2d, Freq_common, 'linear');
psiTE_rot_interp_2d = interp1(res_rot.Freq, psiTE_rot_2d, Freq_common, 'linear');
psiTM_rot = reshape(psiTM_rot_interp_2d.', L, M, numf);
psiTE_rot = reshape(psiTE_rot_interp_2d.', L, M, numf);

fprintf('Interpolated rotated-model coefficients onto %d unrotated-model frequency points (rotated sweep had %d points)\n', ...
        numf, length(res_rot.Freq));

%% ---------------- Build the SH rotation matrix ----------------
R = rotationmat3D(rot_angle, rot_axis);
fprintf('Rotation: axis=[%.0f %.0f %.0f], angle=%.0f deg -- det(R)=%.10f, max|R*R^T-I|=%.2e\n', ...
        rot_axis, rad2deg(rot_angle), det(R), max(max(abs(R*R'-eye(3)))));

Rtb = SHRotate(R, lmax);   % Rtb{l+1} = (2l+1)x(2l+1) real orthogonal rotation matrix for degree l

%% ---------------- Rotate the unrotated model's coefficients, band by band ----------------
% For each l, take the (2l+1)-wide m = -l..l slice (centred layout, as in
% compare_sh_rotation.m), rotate it at every frequency, and write it back.
psiTM_via_matrix = nan(L, M, numf);
psiTE_via_matrix = nan(L, M, numf);

for lindx = 1:L
    l = lval(lindx);
    col_range = (lmax+1-l) : (lmax+1+l);   % m = -l .. l within the fixed M=2*lmax+1 width
    Rl = Rtb{l+1};                          % (2l+1)x(2l+1) real orthogonal

    TM_slice = squeeze(psiTM_unrot(lindx, col_range, :));   % (2l+1) x numf, complex
    TE_slice = squeeze(psiTE_unrot(lindx, col_range, :));

    % dest(:,f) = Rl.' * src(:,f) so that, equivalently, for each f the
    % row vector src(:,f).' is right-multiplied by Rl: src(:,f).' * Rl
    psiTM_via_matrix(lindx, col_range, :) = Rl.' * TM_slice;
    psiTE_via_matrix(lindx, col_range, :) = Rl.' * TE_slice;
end

%% ---------------- Find dominant modes (from the directly-computed rotated model) ----------------
% "Dominant" = largest peak |coefficient| over the common frequency range,
% pooling TM and TE (no frequency-band restriction).
pool_mag  = [max(abs(reshape(psiTM_rot, L*M, numf)), [], 2); ...
             max(abs(reshape(psiTE_rot, L*M, numf)), [], 2)];
pool_type = [zeros(L*M,1); ones(L*M,1)];             % 0=TM, 1=TE
pool_lm   = [(1:L*M).'; (1:L*M).'];

valid = ~isnan(pool_mag) & pool_mag > 0;
pool_mag(~valid) = -Inf;

[~, order] = sort(pool_mag, 'descend');
order = order(1:min(n_dom, sum(valid)));

fprintf('\n%-4s %-4s %-4s | %-12s %-12s | %-10s %-10s\n', ...
        'l','m','type','|direct|','|via matrix|','rel.err','phase diff(deg)');
fprintf('%s\n', repmat('-',1,70));

dom_results = struct('l',{},'m',{},'type',{},'mag_direct',{},'mag_matrix',{}, ...
                      'rel_err',{},'phase_diff_deg',{},'lindx',{},'mindx',{});

for kk = 1:length(order)
    idx = order(kk);
    [lindx, mindx] = ind2sub([L,M], pool_lm(idx));
    l = lval(lindx); m = mval(mindx);
    is_TE = pool_type(idx) == 1;

    if is_TE
        vec_direct = squeeze(psiTE_rot(lindx,mindx,:));
        vec_matrix = squeeze(psiTE_via_matrix(lindx,mindx,:));
        type_str = 'TE';
    else
        vec_direct = squeeze(psiTM_rot(lindx,mindx,:));
        vec_matrix = squeeze(psiTM_via_matrix(lindx,mindx,:));
        type_str = 'TM';
    end

    [mag_direct, fpk] = max(abs(vec_direct));
    mag_matrix = abs(vec_matrix(fpk));
    rel_err = abs(vec_matrix(fpk) - vec_direct(fpk)) / max(mag_direct, eps);
    phase_diff = rad2deg(angle(vec_matrix(fpk)) - angle(vec_direct(fpk)));
    phase_diff = mod(phase_diff + 180, 360) - 180;   % wrap to [-180,180]

    fprintf('%-4d %-4d %-4s | %-12.4e %-12.4e | %-10.2e %-10.2f\n', ...
            l, m, type_str, mag_direct, mag_matrix, rel_err, phase_diff);

    dom_results(kk) = struct('l',l,'m',m,'type',type_str, ...
        'mag_direct',mag_direct,'mag_matrix',mag_matrix, ...
        'rel_err',rel_err,'phase_diff_deg',phase_diff, ...
        'lindx',lindx,'mindx',mindx);
end

% (l,m) pairs in dom_results, pooled across TM and TE. Not used for the
% error metric below, which needs per-type masks.
dom_lindx_mindx = unique([[dom_results.lindx]', [dom_results.mindx]'], 'rows');

%% ---------------- Overall agreement metric (dominant modes only) ----------------
% Relative L2-norm error over the dominant modes of each type only.
% Non-dominant coefficients sit at the numerical noise floor of the
% near-field projection, so their relative error is meaningless. TM and TE
% masks are built separately from the dom_results entries of that type,
% so a weak TE coefficient is not included just because the same (l,m) is
% dominant in TM.

valid_mask   = ~isnan(psiTM_rot) & ~isnan(psiTM_via_matrix);
valid_mask_E = ~isnan(psiTE_rot) & ~isnan(psiTE_via_matrix);

is_TM = strcmp({dom_results.type}, 'TM');
is_TE = strcmp({dom_results.type}, 'TE');
dom_lindx_mindx_TM = unique([[dom_results(is_TM).lindx]', [dom_results(is_TM).mindx]'], 'rows');
dom_lindx_mindx_TE = unique([[dom_results(is_TE).lindx]', [dom_results(is_TE).mindx]'], 'rows');

dom_mask_TM = false(L,M);
for kk = 1:size(dom_lindx_mindx_TM,1)
    dom_mask_TM(dom_lindx_mindx_TM(kk,1), dom_lindx_mindx_TM(kk,2)) = true;
end
dom_mask_TE = false(L,M);
for kk = 1:size(dom_lindx_mindx_TE,1)
    dom_mask_TE(dom_lindx_mindx_TE(kk,1), dom_lindx_mindx_TE(kk,2)) = true;
end
dom_mask_TM = repmat(dom_mask_TM, 1, 1, numf) & valid_mask;
dom_mask_TE = repmat(dom_mask_TE, 1, 1, numf) & valid_mask_E;

num_TM_dom = norm(psiTM_via_matrix(dom_mask_TM) - psiTM_rot(dom_mask_TM));
den_TM_dom = norm(psiTM_rot(dom_mask_TM));
num_TE_dom = norm(psiTE_via_matrix(dom_mask_TE) - psiTE_rot(dom_mask_TE));
den_TE_dom = norm(psiTE_rot(dom_mask_TE));

n_TM_total = nnz(~isnan(psiTM_rot(:,:,1)));
n_TE_total = nnz(~isnan(psiTE_rot(:,:,1)));

fprintf('\nOverall relative L2 error, dominant modes only:\n');
if den_TM_dom > 0
    fprintf('  TM (%d of %d (l,m) pairs): ||via_matrix - direct|| / ||direct|| = %.3e\n', ...
            size(dom_lindx_mindx_TM,1), n_TM_total, num_TM_dom/den_TM_dom);
else
    fprintf('  TM: N/A -- no TM modes among the %d dominant entries\n', length(dom_results));
end
if den_TE_dom > 0
    fprintf('  TE (%d of %d (l,m) pairs): ||via_matrix - direct|| / ||direct|| = %.3e\n', ...
            size(dom_lindx_mindx_TE,1), n_TE_total, num_TE_dom/den_TE_dom);
else
    fprintf('  TE: N/A -- no TE modes among the %d dominant entries\n', length(dom_results));
end

%% ---------------- Plot: dominant modes, magnitude vs frequency (separate figure) ----------------
n_plot = min(6, length(dom_results));
ntile = ceil(sqrt(n_plot));
figure('Color','w','Position',[100 100 350*ntile 300*ntile]);
for kk = 1:n_plot
    d = dom_results(kk);
    if strcmp(d.type,'TM')
        vd = squeeze(psiTM_rot(d.lindx,d.mindx,:));
        vm = squeeze(psiTM_via_matrix(d.lindx,d.mindx,:));
    else
        vd = squeeze(psiTE_rot(d.lindx,d.mindx,:));
        vm = squeeze(psiTE_via_matrix(d.lindx,d.mindx,:));
    end

    subplot(ntile, ntile, kk);
    plot(Freq_common/1e9, abs(vd), 'b-', 'LineWidth', 2); hold on
    plot(Freq_common/1e9, abs(vm), 'r--', 'LineWidth', 2);
    hold off
    title(sprintf('%s, $l=%d$, $m=%d$', d.type, d.l, d.m));
    xlabel('Freq (GHz)'); ylabel('$|\psi|$');
    legend('direct (rotated model)', 'via rotation matrix', 'Location','best', 'Interpreter','latex');
    grid on
end
sgtitle_str = sprintf('Dominant mode magnitude: direct vs rotation-matrix-applied ($%.0f^\\circ$ about $y$-axis)', rad2deg(rot_angle));
annotation('textbox',[0 0.96 1 0.04],'String',sgtitle_str,'EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',12,'Interpreter','latex');

%% ---------------- Plot: dominant modes, phase vs frequency (separate figure) ----------------
figure('Color','w','Position',[100 100 350*ntile 300*ntile]);
for kk = 1:n_plot
    d = dom_results(kk);
    if strcmp(d.type,'TM')
        vd = squeeze(psiTM_rot(d.lindx,d.mindx,:));
        vm = squeeze(psiTM_via_matrix(d.lindx,d.mindx,:));
    else
        vd = squeeze(psiTE_rot(d.lindx,d.mindx,:));
        vm = squeeze(psiTE_via_matrix(d.lindx,d.mindx,:));
    end

    subplot(ntile, ntile, kk);
    plot(Freq_common/1e9, rad2deg(angle(vd)), 'b-', 'LineWidth', 2); hold on
    plot(Freq_common/1e9, rad2deg(angle(vm)), 'r--', 'LineWidth', 2);
    hold off
    title(sprintf('%s, $l=%d$, $m=%d$', d.type, d.l, d.m));
    xlabel('Freq (GHz)'); ylabel('$\angle\psi$ (deg)');
    legend('direct (rotated model)', 'via rotation matrix', 'Location','best', 'Interpreter','latex');
    ylim([-180 180]); grid on
end
sgtitle_str = sprintf('Dominant mode phase: direct vs rotation-matrix-applied ($%.0f^\\circ$ about $y$-axis)', rad2deg(rot_angle));
annotation('textbox',[0 0.96 1 0.04],'String',sgtitle_str,'EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',12,'Interpreter','latex');

%% ---------------- Stem plot: DOMINANT modes only, dB magnitude and phase ----------------
% Compact x-axis: l = 1..lmax, m = -l..l within each l (not the FEKO
% j = 2*(l*(l+1)+m-1)+s numbering). TM and TE get separate figures using
% their own dominant-mode sets; a type with no dominant modes is skipped.
% All modes are evaluated at the frequency in Freq_common nearest to f0.

[~, f0_idx] = min(abs(Freq_common - f0));
fprintf('\nStem plots evaluated at f0 = %.4g Hz (nearest available: %.4g Hz)\n', ...
        f0, Freq_common(f0_idx));

plot_stem_type('TM', dom_lindx_mindx_TM, lval, mval, psiTM_rot, psiTM_via_matrix, ...
                f0_idx, Freq_common(f0_idx), rot_angle);
plot_stem_type('TE', dom_lindx_mindx_TE, lval, mval, psiTE_rot, psiTE_via_matrix, ...
                f0_idx, Freq_common(f0_idx), rot_angle);