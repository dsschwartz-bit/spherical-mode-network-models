function fignum = plot_cst_farfield_cuts(fignum, theta, phi, Eph, Eth, ...
                                          index_phi_cut, index_th_cut, ...
                                          mag_scale, field_label, title_str, FS, LW)
%PLOT_CST_FARFIELD_CUTS  Theta and phi cuts of the total far-field magnitude
% |E_total| = sqrt(|Eph|^2 + |Eth|^2) of one source, as two subplots.
%
%   fignum        : last figure number used (incremented and returned)
%   theta, phi    : angular grid (rad)
%   Eph, Eth      : (numth x numph) complex components at one frequency
%   index_phi_cut : phi index held fixed for the theta cut
%   index_th_cut  : theta index held fixed for the phi cut
%   mag_scale     : multiplier on the magnitude (e.g. 1/eta0 to show H)
%   field_label   : LaTeX y-label, e.g. '$|E_{\mathrm{total}}|$'
%   title_str     : figure title (string or cell array)
%   FS, LW        : font size, line width
%
% See also plot_cst_farfield_vector, plot_field_comparison.

E_total = sqrt(abs(Eph).^2 + abs(Eth).^2);

phi_label   = "$\phi="   + string(round(phi(index_phi_cut)*180/pi,1))   + "^\circ$";
theta_label = "$\theta=" + string(round(theta(index_th_cut)*180/pi,1)) + "^\circ$";

fignum = fignum + 1;
figure(fignum);

% --- Theta cut (fixed phi) ---
subplot(2,1,1)
plot(theta*180/pi, mag_scale*E_total(:,index_phi_cut), "LineWidth", LW, "Color", "b")
xlabel("$\theta$ (deg), " + phi_label, "FontSize", FS, "Interpreter", "latex")
ylabel(field_label, "FontSize", FS, "Interpreter", "latex")
xlim([0 180])
ylim([0 inf])
grid on
title("vs $\theta$", "FontSize", FS, "Interpreter", "latex")

% --- Phi cut (fixed theta) ---
subplot(2,1,2)
plot(phi*180/pi, mag_scale*E_total(index_th_cut,:), "LineWidth", LW, "Color", "r")
xlabel("$\phi$ (deg), " + theta_label, "FontSize", FS, "Interpreter", "latex")
ylabel(field_label, "FontSize", FS, "Interpreter", "latex")
xlim([0 360])
ylim([0 inf])
grid on
title("vs $\phi$", "FontSize", FS, "Interpreter", "latex")

sgtitle(title_str, "FontSize", FS+2, "FontWeight", "bold", "Interpreter", "latex")

end
