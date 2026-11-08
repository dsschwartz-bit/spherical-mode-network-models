function fignum = plot_field_comparison_polar_OMNI(fignum, ...
                                         th_sweep, field_psi_th, field_FEKO_th, field_computed_th, ...
                                         mag_scale, field_label, title_str, FS, LW, phi_label)
%PLOT_FIELD_COMPARISON_POLAR_OMNI  Polar plot of reconstructed vs simulated field vs theta only.
%   fignum = plot_field_comparison_polar_OMNI(fignum, th_sweep, field_psi_th,
%       field_FEKO_th, field_computed_th, mag_scale, field_label, title_str,
%       FS, LW, phi_label)
%   For omnidirectional (TM_l0-only) antennas, which have no phi variation
%   and only E_r, E_theta and H_phi. Arguments as in
%   plot_field_comparison_polar (theta cut only); pass [] for
%   field_FEKO_th to omit the FEKO Q trace. title_str may be a string
%   (lines split at \n) or a cell array of lines.
%
%   See also plot_field_comparison_polar.

narginchk(11, 11);
fignum = fignum + 1;
figure(fignum);

ax = polaraxes();

angle_rad   = th_sweep(:) * pi;
angle_clean = angle_rad;

% Pole mask (sin(theta) = 0) for the reconstructed traces.
ispole = abs(sin(angle_clean)) < eps;

mag_psi_clean = mag_scale * abs(field_psi_th(:));
mag_psi_clean(isnan(mag_psi_clean) | ispole) = NaN;

h_psi = polarplot(ax, angle_clean, mag_psi_clean, "LineWidth", LW+2, "LineStyle", "-");
hold(ax, 'on')

h_FEKO = [];
if ~isempty(field_FEKO_th)
    mag_FEKO_clean = mag_scale * abs(field_FEKO_th(:));
    mag_FEKO_clean(isnan(mag_FEKO_clean) | ispole) = NaN;
    h_FEKO = polarplot(ax, angle_clean, mag_FEKO_clean, "LineWidth", LW-0.5, "LineStyle", "--");
end

h_comp = [];
if ~isempty(field_computed_th)
    valid_c = ~isnan(angle_clean) & ~isnan(mag_scale*abs(field_computed_th(:)));
    h_comp = polarplot(ax, angle_clean(valid_c), mag_scale*abs(field_computed_th(valid_c)), ...
                       "LineStyle", "none", "Marker", "*", "LineWidth", LW-1, "Color", [0 0 0]);
end
hold(ax, 'off')

% Restrict angular axis to data range
ax.ThetaLim = [min(angle_clean)*180/pi, max(angle_clean)*180/pi];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir          = 'clockwise';
ax.FontSize = FS - 2;

% Title lines from title_str, with field_label/phi_label as the last line.
sub_label = [field_label ' vs $\theta$, ' char(phi_label)];
if iscell(title_str)
    title_lines = [title_str(:)', {sub_label}];
else
    parts = strsplit(char(title_str), '\n');
    title_lines = [parts(:)', {sub_label}];
end
title(ax, title_lines, "FontSize", FS, "Interpreter", "latex")
grid(ax, 'on')

if ~isempty(h_FEKO) && ~isempty(h_comp)
    legend([h_psi, h_FEKO, h_comp], {"Reconstructed $\psi$", "Reconstructed FEKO $Q$", "Computed"}, ...
        "FontSize", FS, "Interpreter", "latex")
elseif ~isempty(h_FEKO)
    legend([h_psi, h_FEKO], {"Reconstructed $\psi$", "Reconstructed FEKO $Q$"}, ...
        "FontSize", FS, "Interpreter", "latex")
elseif ~isempty(h_comp)
    legend([h_psi, h_comp], {"Reconstructed $\psi$", "Computed"}, ...
        "FontSize", FS, "Interpreter", "latex")
else
    legend(h_psi, {"Reconstructed $\psi$"}, "FontSize", FS, "Interpreter", "latex")
end

end