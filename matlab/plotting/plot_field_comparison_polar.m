function fignum = plot_field_comparison_polar(fignum, ...
                                         th_sweep, field_psi_th, field_FEKO_th, field_computed_th, ...
                                         ph_sweep, field_psi_ph, field_FEKO_ph, field_computed_ph, ...
                                         mag_scale, field_label, title_str, FS, LW, phi_label, theta_label)
%PLOT_FIELD_COMPARISON_POLAR  Polar plots of reconstructed vs simulated field magnitude.
%   fignum = plot_field_comparison_polar(...)  takes the same 16 arguments
%   as plot_field_comparison. Left: vs theta (0 at top, clockwise); right:
%   vs phi. Reconstructed traces are blanked where sin(theta) = 0 on the
%   theta cut; the simulator field (markers) is not.
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
%   See also plot_field_comparison, plot_field_comparison_polar_OMNI.

narginchk(16, 16);
fignum = fignum + 1;
figure(fignum);

% Axes positions leave a bottom margin for the legend.
positions = [0.05 0.12 0.38 0.62;
             0.46 0.12 0.38 0.62];

for sp = 1:2

    ax = polaraxes('Position', positions(sp,:));

    if sp == 1
        angle_rad = th_sweep * pi;
        f_psi     = field_psi_th;
        f_FEKO    = field_FEKO_th;
        f_comp    = field_computed_th;
        sub_label = [field_label ' vs $\theta$, ' char(phi_label)];
    else
        angle_rad = ph_sweep * pi;
        f_psi     = field_psi_ph;
        f_FEKO    = field_FEKO_ph;
        f_comp    = field_computed_ph;
        sub_label = [field_label ' vs $\phi$, ' char(theta_label)];
    end

    angle_clean = angle_rad(:);

    % Pole mask (theta cut only) for the reconstructed traces.
    if sp == 1
        ispole = abs(sin(angle_clean)) < eps;
    else
        ispole = false(size(angle_clean));
    end

    mag_psi_clean = mag_scale * abs(f_psi(:));
    mag_psi_clean(isnan(mag_psi_clean) | ispole) = NaN;

    h_psi = polarplot(ax, angle_clean, mag_psi_clean, "LineWidth", LW+2, "LineStyle", "-");
    hold(ax, 'on')
    if ~isempty(f_FEKO)
        mag_FEKO_clean = mag_scale * abs(f_FEKO(:));
        mag_FEKO_clean(isnan(mag_FEKO_clean) | ispole) = NaN;
        h_FEKO = polarplot(ax, angle_clean, mag_FEKO_clean, "LineWidth", LW-0.5, "LineStyle", "--");
    else
        h_FEKO = [];
    end
    baseColor = theme().BaseColorStyle;
    if strcmp(baseColor, 'dark')
        computed_color = [1 1 1];
    elseif strcmp(baseColor, 'light')
        computed_color = [0 0 0];
    else
        computed_color = [0.5 0.5 0.5];
    end
    % Simulator field: no pole mask.
    valid_c = ~isnan(angle_clean);
    if ~isempty(f_comp)
        valid_c = valid_c & ~isnan(mag_scale*abs(f_comp(:)));
    end
    if ~isempty(f_comp)
        h_comp = polarplot(ax, angle_clean(valid_c), mag_scale*abs(f_comp(valid_c)), ...
                           "LineStyle", "none", "Marker", "*", "LineWidth", LW-1, "Color", computed_color);
    else
        h_comp = [];
    end
    hold(ax, 'off')

    % Limit the angular axis to the data range so the arc is not closed.
    valid_angles = angle_clean(~isnan(angle_clean));
    if length(valid_angles) >= 2 && (max(valid_angles) - min(valid_angles)) > 0
        ax.ThetaLim = [min(valid_angles)*180/pi, max(valid_angles)*180/pi];
    end
    hold(ax, 'off')

    % Set after plotting, since polarplot resets these.
    if sp == 1
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir          = 'clockwise';
    else
        ax.ThetaZeroLocation = 'right';
        ax.ThetaDir          = 'counterclockwise';
    end

    ax.FontSize = FS - 2;
    title(ax, sub_label, "FontSize", FS, "Interpreter", "latex")
    grid(ax, 'on')

    if sp == 1
        if ~isempty(h_FEKO) && ~isempty(h_comp)
            leg_handles = [h_psi, h_FEKO, h_comp];
            leg_labels  = {"Reconstructed $\psi$", "Reconstructed FEKO $Q$", "Computed"};
        elseif ~isempty(h_FEKO)
            leg_handles = [h_psi, h_FEKO];
            leg_labels  = {"Reconstructed $\psi$", "Reconstructed FEKO $Q$"};
        elseif ~isempty(h_comp)
            leg_handles = [h_psi, h_comp];
            leg_labels  = {"Reconstructed $\psi$", "Computed"};
        else
            leg_handles = h_psi;
            leg_labels  = {"Reconstructed $\psi$"};
        end
    end
end

% Horizontal legend centred at the bottom of the figure.
lg = legend(leg_handles, leg_labels, "FontSize", FS, "Interpreter", "latex");
lg.Orientation = 'horizontal';
lg.Units         = 'normalized';
lg.Position(1)   = 0.435 - lg.Position(3)/2;
lg.Position(2)   = 0.01;                       % bottom

sg = sgtitle(title_str, "FontSize", FS+4, "FontWeight", "bold", "Interpreter", "latex");

end