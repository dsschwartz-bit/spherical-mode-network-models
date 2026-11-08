function fignum = plot_field_comparison(fignum, ...
                                         th_sweep, field_psi_th, field_FEKO_th, field_computed_th, ...
                                         ph_sweep, field_psi_ph, field_FEKO_ph, field_computed_ph, ...
                                         mag_scale, field_label, title_str, FS, LW, phi_label, theta_label)
%PLOT_FIELD_COMPARISON  Magnitude/phase of reconstructed vs simulated field along two cuts.
%   fignum = plot_field_comparison(fignum, th_sweep, field_psi_th, field_FEKO_th,
%       field_computed_th, ph_sweep, field_psi_ph, field_FEKO_ph, field_computed_ph,
%       mag_scale, field_label, title_str, FS, LW, phi_label, theta_label)
%   Top subplot: vs theta; bottom: vs phi. Magnitude on the left axis
%   (solid lines), phase in degrees on the right axis (dotted lines);
%   the simulator field is drawn with markers.
%
%   fignum        figure counter; the plot uses fignum+1, which is returned
%   th_sweep      theta values in units of pi (theta/pi)
%   field_psi_th  field reconstructed from the psi coefficients vs theta
%   field_FEKO_th field reconstructed from FEKO Q vs theta, or [] to omit
%   field_computed_th  simulator field vs theta (markers)
%   ph_sweep, field_*_ph  the same along the phi cut (phi/pi)
%   mag_scale     factor applied to the magnitudes
%   field_label   LaTeX label of the field (e.g. '$|E_r|$')
%   title_str     figure title (LaTeX)
%   FS, LW        font size and base line width
%   phi_label, theta_label  LaTeX text naming the fixed angle of each cut
%
% See also plot_field_comparison_polar.

fignum = fignum + 1;
figure(fignum);

for sp = 1:2
    subplot(2, 1, sp);

    if sp == 1
        sweep  = th_sweep;
        f_psi  = field_psi_th;
        f_FEKO = field_FEKO_th;
        f_comp = field_computed_th;
        xlab = ['$\theta/\pi$, ' char(phi_label)];
    else
        sweep  = ph_sweep;
        f_psi  = field_psi_ph;
        f_FEKO = field_FEKO_ph;
        f_comp = field_computed_ph;
        xlab = ['$\phi/\pi$, ' char(theta_label)];
    end

    yyaxis left
    set(gca, 'YColor', get(groot,'defaultAxesXColor'))
    colororder(gca, get(groot,'defaultAxesColorOrder'))
    set(gca, 'ColorOrderIndex', 1)
    hold on
    plot(sweep, mag_scale*abs(f_psi(:)),  "LineWidth", LW+2,   "LineStyle", "-")
    if ~isempty(f_FEKO)
        plot(sweep, mag_scale*abs(f_FEKO(:)), "LineWidth", LW-0.5, "LineStyle", "-")
    end
    baseColor = theme().BaseColorStyle;
    if strcmp(baseColor, 'dark')
        computed_color = [1 1 1];
    elseif strcmp(baseColor, 'light')
        computed_color = [0 0 0];
    else
        computed_color = [0.5 0.5 0.5];
    end
    plot(sweep, mag_scale*abs(f_comp(:)), "*", "LineStyle", "none", "LineWidth", LW-1, "Color", computed_color)
    hold off
    ylabel(field_label, "FontSize", FS, "Interpreter", "latex")

    yyaxis right
    set(gca, 'YColor', get(groot,'defaultAxesXColor'))
    colororder(gca, get(groot,'defaultAxesColorOrder'))
    set(gca, 'ColorOrderIndex', 1)
    hold on
    plot(sweep, angle(f_psi(:))*180/pi,  "LineWidth", LW+1,   "LineStyle", ":")
    if ~isempty(f_FEKO)
        plot(sweep, angle(f_FEKO(:))*180/pi, "LineWidth", LW-0.5, "LineStyle", ":")
    end
    plot(sweep, angle(f_comp(:))*180/pi, "^", "LineStyle", "none", "LineWidth", LW-1, "Color", computed_color)
    hold off
    ylim([-180 180])
    ylabel("phase (deg)", "FontSize", FS, "Interpreter", "latex")

    xlabel(xlab, "FontSize", FS, "Interpreter", "latex")
    grid on

    if sp == 1
        legend_handles = flipud(get(gca, 'Children'));
    end
end

% Shared legend to the right of both subplots.
if ~isempty(field_FEKO_th)
    leg_labels = {"Reconstructed $\psi$", "Reconstructed FEKO $Q$", "Computed"};
else
    leg_labels = {"Reconstructed $\psi$", "Computed"};
end
lg = legend(legend_handles(1:numel(leg_labels)), leg_labels, "FontSize", FS, "Interpreter", "latex");
lg.Orientation = 'vertical';
lg.Units       = 'normalized';
lg.Position(1) = 0.88;
lg.Position(2) = 0.5 - lg.Position(4)/2;

% Shrink subplots to make room for legend
for sp = 1:2
    ax = subplot(2, 1, sp);
    ax.Position(3) = 0.72;
end

sgtitle(title_str, "FontSize", FS, "Interpreter", "latex")

end