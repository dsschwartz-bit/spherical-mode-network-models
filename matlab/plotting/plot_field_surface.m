function fignum = plot_field_surface(fignum, ...
                                      theta, phi, field, ...
                                      th_sweep, field_psi_th, field_FEKO_th, field_computed_th, ...
                                      ph_sweep, field_psi_ph, field_FEKO_ph, field_computed_ph, ...
                                      mag_scale, field_label, title_str, FS, LW, phi_label, theta_label)
%PLOT_FIELD_SURFACE  3-D pattern surface of a field plus magnitude/phase cuts.
%   fignum = plot_field_surface(fignum, theta, phi, field, th_sweep, ...
%       field_psi_th, field_FEKO_th, field_computed_th, ph_sweep, ...
%       field_psi_ph, field_FEKO_ph, field_computed_ph, mag_scale, ...
%       field_label, title_str, FS, LW, phi_label, theta_label)
%   Top: surface with radius mag_scale*|field| and colour = phase (cyclic
%   hsv map). Bottom left/right: magnitude (left axis) and phase (right
%   axis) vs theta at fixed phi and vs phi at fixed theta.
%
%   INPUTS:
%   fignum          - current figure number (incremented and returned)
%   theta           - theta vector (numth x 1), radians
%   phi             - phi vector (1 x numph), radians
%   field           - complex field at one frequency (numth x numph)
%   th_sweep        - theta/pi vector for phase plot x-axis
%   field_psi_th    - reconstructed field vs theta (for phase plot)
%   field_FEKO_th   - FEKO Q reconstruction vs theta ([] to omit)
%   field_computed_th - computed reference vs theta ([] to omit)
%   ph_sweep        - phi/pi vector for phase plot x-axis
%   field_psi_ph    - reconstructed field vs phi
%   field_FEKO_ph   - FEKO Q reconstruction vs phi ([] to omit)
%   field_computed_ph - computed reference vs phi ([] to omit)
%   mag_scale       - scalar magnitude multiplier (e.g. eta0 or 1)
%   field_label     - magnitude label string (LaTeX)
%   title_str       - overall figure title
%   FS              - FontSize
%   LW              - LineWidth base
%   phi_label       - fixed phi label string for theta phase subplot
%   theta_label     - fixed theta label string for phi phase subplot

fignum = fignum + 1;
figure(fignum);

%% --- Spherical surface plot (top, spanning full width) ---
ax_surf = subplot(2, 2, [1 2]);

[PHI, THETA] = meshgrid(phi, theta);

mag_field   = mag_scale * abs(field);
phase_field = angle(field) * 180/pi;   % degrees, in (-180, 180]

X = mag_field .* sin(THETA) .* cos(PHI);
Y = mag_field .* sin(THETA) .* sin(PHI);
Z = mag_field .* cos(THETA);

surf(ax_surf, X, Y, Z, phase_field, ...
     'EdgeColor', 'none', 'FaceColor', 'interp')
% Cyclic colormap, so +-180 deg get the same colour.
colormap(ax_surf, hsv(256))
clim(ax_surf, [-180 180])
cb = colorbar(ax_surf);
cb.Label.Interpreter = 'latex';
cb.Label.String = 'phase (deg)';
cb.Label.FontSize = FS;
xlabel(ax_surf, '$x$', 'FontSize', FS, 'Interpreter', 'latex')
ylabel(ax_surf, '$y$', 'FontSize', FS, 'Interpreter', 'latex')
zlabel(ax_surf, '$z$', 'FontSize', FS, 'Interpreter', 'latex')
title(ax_surf, field_label+" (radius = magnitude, color = phase)", 'FontSize', FS, 'Interpreter', 'latex')
axis(ax_surf, 'equal')
view(ax_surf, -45, 30)
grid(ax_surf, 'on')

baseColor = theme().BaseColorStyle;
if strcmp(baseColor, 'dark'),       computed_color = [1 1 1];
elseif strcmp(baseColor, 'light'),  computed_color = [0 0 0];
else,                               computed_color = [0.5 0.5 0.5];
end

%% --- Magnitude + Phase vs theta (bottom left) ---
ax_th = subplot(2, 2, 3);

yyaxis(ax_th, 'left')
set(ax_th, 'YColor', get(groot,'defaultAxesXColor'))
colororder(ax_th, get(groot,'defaultAxesColorOrder'))
set(ax_th, 'ColorOrderIndex', 1)
hold(ax_th, 'on')
plot(ax_th, th_sweep, mag_scale*abs(field_psi_th(:)),  "LineWidth", LW+1, "LineStyle", "-")
if ~isempty(field_FEKO_th)
    plot(ax_th, th_sweep, mag_scale*abs(field_FEKO_th(:)), "LineWidth", LW-0.5, "LineStyle", "-")
end
if ~isempty(field_computed_th)
    plot(ax_th, th_sweep, mag_scale*abs(field_computed_th(:)), ...
         "LineStyle", "none", "Marker", "*", "LineWidth", LW-1, "Color", computed_color)
end
hold(ax_th, 'off')
ylabel(ax_th, field_label, 'FontSize', FS, 'Interpreter', 'latex')
ylim(ax_th, [0 inf])

yyaxis(ax_th, 'right')
set(ax_th, 'YColor', get(groot,'defaultAxesXColor'))
colororder(ax_th, get(groot,'defaultAxesColorOrder'))
set(ax_th, 'ColorOrderIndex', 1)
hold(ax_th, 'on')
plot(ax_th, th_sweep, angle(field_psi_th(:))*180/pi,  "LineWidth", LW+1, "LineStyle", ":")
if ~isempty(field_FEKO_th)
    plot(ax_th, th_sweep, angle(field_FEKO_th(:))*180/pi, "LineWidth", LW-0.5, "LineStyle", ":")
end
if ~isempty(field_computed_th)
    plot(ax_th, th_sweep, angle(field_computed_th(:))*180/pi, ...
         "LineStyle", "none", "Marker", "^", "LineWidth", LW-1, "Color", computed_color)
end
hold(ax_th, 'off')
ylim(ax_th, [-180 180])
ylabel(ax_th, 'phase (deg)', 'FontSize', FS, 'Interpreter', 'latex')
xlabel(ax_th, ['$\theta/\pi$, ' char(phi_label)], 'FontSize', FS, 'Interpreter', 'latex')
grid(ax_th, 'on')

%% --- Magnitude + Phase vs phi (bottom right) ---
ax_ph = subplot(2, 2, 4);

yyaxis(ax_ph, 'left')
set(ax_ph, 'YColor', get(groot,'defaultAxesXColor'))
colororder(ax_ph, get(groot,'defaultAxesColorOrder'))
set(ax_ph, 'ColorOrderIndex', 1)
hold(ax_ph, 'on')
plot(ax_ph, ph_sweep, mag_scale*abs(field_psi_ph(:)),  "LineWidth", LW+1, "LineStyle", "-")
if ~isempty(field_FEKO_ph)
    plot(ax_ph, ph_sweep, mag_scale*abs(field_FEKO_ph(:)), "LineWidth", LW-0.5, "LineStyle", "-")
end
if ~isempty(field_computed_ph)
    plot(ax_ph, ph_sweep, mag_scale*abs(field_computed_ph(:)), ...
         "LineStyle", "none", "Marker", "*", "LineWidth", LW-1, "Color", computed_color)
end
hold(ax_ph, 'off')
ylabel(ax_ph, field_label, 'FontSize', FS, 'Interpreter', 'latex')
ylim(ax_ph, [0 inf])

yyaxis(ax_ph, 'right')
set(ax_ph, 'YColor', get(groot,'defaultAxesXColor'))
colororder(ax_ph, get(groot,'defaultAxesColorOrder'))
set(ax_ph, 'ColorOrderIndex', 1)
hold(ax_ph, 'on')
plot(ax_ph, ph_sweep, angle(field_psi_ph(:))*180/pi,  "LineWidth", LW+1, "LineStyle", ":")
if ~isempty(field_FEKO_ph)
    plot(ax_ph, ph_sweep, angle(field_FEKO_ph(:))*180/pi, "LineWidth", LW-0.5, "LineStyle", ":")
end
if ~isempty(field_computed_ph)
    plot(ax_ph, ph_sweep, angle(field_computed_ph(:))*180/pi, ...
         "LineStyle", "none", "Marker", "^", "LineWidth", LW-1, "Color", computed_color)
end
hold(ax_ph, 'off')
ylim(ax_ph, [-180 180])
ylabel(ax_ph, 'phase (deg)', 'FontSize', FS, 'Interpreter', 'latex')
xlabel(ax_ph, ['$\phi/\pi$, ' char(theta_label)], 'FontSize', FS, 'Interpreter', 'latex')
grid(ax_ph, 'on')

%% --- Shared legend ---
yyaxis(ax_th, 'left')
leg_handles = flipud(ax_th.Children);

if ~isempty(field_FEKO_th) && ~isempty(field_computed_th)
    leg_labels = {"Reconstructed $\psi$", "Reconstructed FEKO $Q$", "Computed"};
elseif ~isempty(field_FEKO_th)
    leg_labels = {"Reconstructed $\psi$", "Reconstructed FEKO $Q$"};
elseif ~isempty(field_computed_th)
    leg_labels = {"Reconstructed $\psi$", "Computed"};
else
    leg_labels = {"Reconstructed $\psi$"};
end

legend(ax_th, leg_handles(1:numel(leg_labels)), leg_labels, "FontSize", FS-2, "Interpreter", "latex", "Location", "best");

%% --- Overall title ---
sgtitle(title_str, "FontSize", FS+4, "FontWeight", "bold", "Interpreter", "latex")

end