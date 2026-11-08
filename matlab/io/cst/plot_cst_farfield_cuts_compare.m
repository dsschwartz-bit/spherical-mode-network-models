function fignum = plot_cst_farfield_cuts_compare(fignum, ...
        theta_cst, phi_cst, Eph_cst, Eth_cst, theta_ads, phi_ads, Eph_ads, Eth_ads, ...
        phi_cut_rad, theta_cut_rad, cst_ds_step_deg, mag_scale, unit_label, title_str, FS, legend_FS, LW, cut_tol_deg)
% PLOT_CST_FARFIELD_CUTS_COMPARE   Polar plots comparing CST and ADS (2 rows x N cut columns):
%
%                 theta cut(s) at fixed phi      phi cut(s) at fixed theta
%     |E_theta|    [ phi_1 ... phi_n ]       [ theta_1 ... theta_m ]
%     |E_phi|      [ phi_1 ... phi_n ]       [ theta_1 ... theta_m ]
%
% Each requested cut adds one column (left to right: all theta cuts, then
% all phi cuts), so the figure has 2 rows x (n+m) columns.
%
% Each polar axis plots the COMPLEX magnitude abs() of ONE field component
% (not the total), scaled by mag_scale (1, 1e3, 1e6 for V, mV, uV). CST and ADS are
% each taken on their own native grid; the cut angle is matched to each
% source's nearest sample independently (no interpolation).
%
% CST: solid black triangles, downsampled to ~cst_ds_step_deg spacing.
% ADS: continuous red line at its own resolution.
% Theta cut (constant phi0): zenith at top, increasing clockwise, -180..180
%   deg -- phi0 at +theta, phi0+180 at -theta (the full plane through z).
%   Don't request both phi0 and phi0+180: they are the same plane.
% Phi cut: 0 at right, increasing counterclockwise, 0-360 deg.
%
%   theta_*,phi_*,Eph_*,Eth_* : angular grid (theta numth x 1, phi 1 x numph)
%                          and (numth x numph) complex components, per source
%   phi_cut_rad   : vector of fixed phi angles (rad), one THETA cut each
%   theta_cut_rad : vector of fixed theta angles (rad), one PHI cut each
%   cst_ds_step_deg        : nominal CST marker spacing (degrees)
%   mag_scale              : multiplier on the magnitudes (1 -> V)
%   unit_label            : LaTeX unit string shown in panel titles, e.g. "mV"
%   title_str              : overall figure title (e.g. the frequency)
%   FS, legend_FS, LW      : axes/title font size, legend font size, ADS line width
%   cut_tol_deg            : (optional, default 1) max allowed distance between a
%                            requested cut angle and a source's nearest grid
%                            sample. If either source is farther off, that panel
%                            is NOT plotted; a text notice is shown in its place.

if nargin < 19 || isempty(cut_tol_deg), cut_tol_deg = 1; end

theta_cst = theta_cst(:); phi_cst = phi_cst(:).';
theta_ads = theta_ads(:); phi_ads = phi_ads(:).';
phi_cut_rad = phi_cut_rad(:).'; theta_cut_rad = theta_cut_rad(:).';

% Column list: all theta cuts first, then all phi cuts.
col_type  = [repmat('t', 1, numel(phi_cut_rad)), repmat('p', 1, numel(theta_cut_rad))];
col_angle = [phi_cut_rad, theta_cut_rad];
ncol = numel(col_type);
if ncol == 0
    error('plot_cst_farfield_cuts_compare: no cuts requested (both cut-angle lists are empty).');
end

fignum = fignum + 1;
figure(fignum);
% The figure is left docked (MATLAB default: one window, one tab per figure),
% so its size is not set here. Layout is done in PIXELS from the actual
% figure size: each polar axes is a square of side h*(figure height), and
% columns are spaced col_pitch_factor x that diameter (room for the angle
% labels), capped so all ncol columns fit across the figure.
h = 0.30;                             % axes height (normalized); vertical layout unchanged
row_y = [0.60 0.12];                  % vertical layout: gap must hold bottom angle labels + next title
col_pitch_factor = 1.40;              % column spacing / axes diameter (lower = more compact)

set(gcf, 'Units', 'pixels');  drawnow;
fp = get(gcf, 'Position');            % [x y width height] of the figure tab, pixels
set(gcf, 'Units', 'normalized');
diam_pix  = h * fp(4);
pitch_pix = col_pitch_factor * diam_pix;
Wpix  = fp(3);
pitch = min(pitch_pix/Wpix, 0.96/ncol);     % normalized column spacing
w     = min(diam_pix/Wpix, 0.9*pitch);      % square axes box in pixels
left  = (1 - ncol*pitch)/2;                 % centre the block of columns

comp_names = {'$|E_{\theta}|$', '$|E_{\phi}|$'};
cst_comp = {Eth_cst, Eph_cst};
ads_comp = {Eth_ads, Eph_ads};

leg_handles = []; leg_labels = {"CST", "ADS"};

% All panels (both components, every cut, both sources) share ONE radial
% scale, so a component that is essentially zero shows as a dot at the
% centre instead of being auto-scaled up into a misleading full-size pattern.
all_ax = gobjects(0);
r_max  = 0;

for row = 1:2
    for col = 1:ncol
        x0 = left + (col-0.5)*pitch - w/2;
        ax = polaraxes('Position', [x0 row_y(row) w h]);

        if col_type(col) == 't'
            % theta cut: the full plane through z at phi = phi0. The phi0
            % half-plane is drawn at +theta, the phi0+180 half-plane at
            % -theta (theta = 0..180 on each side, so -180..180 overall).
            ph0 = col_angle(col);  ph1 = ph0 + pi;
            [E_cst, angle_cst, got_c] = theta_plane(cst_comp{row}, theta_cst, phi_cst, ph0, ph1, mag_scale);
            [E_ads, angle_ads, got_a] = theta_plane(ads_comp{row}, theta_ads, phi_ads, ph0, ph1, mag_scale);
            want = [ph0 ph1];  sym = '\phi';
            full_range_deg = [-180 180];
            sub_label = sprintf('%s (%s), $\\phi=%.0f^\\circ/%.0f^\\circ$', comp_names{row}, char(unit_label), ...
                mod(ph0*180/pi, 360), mod(ph1*180/pi, 360));
        else
            % phi cut at fixed theta
            [~, ic] = min(abs(theta_cst - col_angle(col)));
            [~, ia] = min(abs(theta_ads - col_angle(col)));
            E_cst = mag_scale * abs(cst_comp{row}(ic, :)).';
            E_ads = mag_scale * abs(ads_comp{row}(ia, :)).';
            angle_cst = phi_cst(:); angle_ads = phi_ads(:);
            got_c = theta_cst(ic);  got_a = theta_ads(ia);
            want = col_angle(col);  sym = '\theta';
            full_range_deg = [0 360];
            sub_label = sprintf('%s (%s), $\\theta=%.0f^\\circ$', comp_names{row}, char(unit_label), col_angle(col)*180/pi);
        end

        % Grid-mismatch check: if either source's nearest sample is more than
        % cut_tol_deg from the requested cut angle, don't plot this panel;
        % show a text notice in its place instead.
        msg = [mismatch_lines('CST', sym, want, got_c, cut_tol_deg), ...
               mismatch_lines('ADS', sym, want, got_a, cut_tol_deg)];
        if ~isempty(msg)
            delete(ax);
            ax = axes('Position', [x0 row_y(row) w h]);
            axis(ax, 'off');
            text(ax, 0.5, 0.5, [{sprintf('\\textbf{Not plotted:}~grid mismatch $>%g^\\circ$', cut_tol_deg)}, msg], ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Interpreter', 'latex', 'FontSize', FS - 4);
            title(ax, sub_label, "FontSize", FS, "Interpreter", "latex")
            if row == 1
                fprintf('plot_cst_farfield_cuts_compare (%s): column %d not plotted -- %s\n', ...
                    strrep(char(title_str), '$', ''), col, strjoin(strrep(strrep(msg, '$', ''), '^\circ', ' deg'), '; '));
            end
            continue
        end

        % Downsample CST only (after extracting the cut at full resolution)
        target_angles_rad = (full_range_deg(1):cst_ds_step_deg:full_range_deg(2)) * pi/180;
        ds_idx = zeros(size(target_angles_rad));
        for ii = 1:numel(target_angles_rad)
            [~, ds_idx(ii)] = min(abs(angle_cst - target_angles_rad(ii)));
        end
        ds_idx = unique(ds_idx, 'stable');

        r_max = max([r_max; E_ads(:); E_cst(ds_idx(:))]);
        all_ax(end+1) = ax; %#ok<AGROW>

        hold(ax, 'on')
        h_ads = polarplot(ax, angle_ads, E_ads, "LineStyle", "-", "Color", "r", "LineWidth", LW, "DisplayName", "ADS");
        h_cst = polarplot(ax, angle_cst(ds_idx), E_cst(ds_idx), "LineStyle", "none", "Marker", "^", "MarkerSize", 4, ...
                          "MarkerFaceColor", "k", "MarkerEdgeColor", "k", "DisplayName", "CST");
        hold(ax, 'off')

        ax.ThetaLim = full_range_deg;
        if col_type(col) == 't'
            ax.ThetaZeroLocation = 'top';
            ax.ThetaDir          = 'clockwise';
        else
            ax.ThetaZeroLocation = 'right';
            ax.ThetaDir          = 'counterclockwise';
        end
        ax.FontSize = FS - 4;   % angle/radial tick labels a bit smaller than titles (FS) to avoid overlap
        title(ax, sub_label, "FontSize", FS, "Interpreter", "latex")
        grid(ax, 'on')

        if isempty(leg_handles)                 % first panel actually plotted
            leg_handles = [h_cst, h_ads];
        end
    end
end

if r_max > 0
    set(all_ax, 'RLim', [0 r_max]);
end

if ~isempty(leg_handles)
    lg = legend(leg_handles, leg_labels, "FontSize", legend_FS, "Interpreter", "latex");
    lg.Orientation = 'horizontal';
    lg.Units       = 'normalized';
    lg.Position(1) = 0.5 - lg.Position(3)/2;
    lg.Position(2) = 0.01;
end

% Overall title (e.g. the frequency) in the top-left corner rather than
% centred, so it does not collide with the subplot titles of the top row.
annotation(gcf, 'textbox', [0.005 0.90 0.25 0.095], 'String', title_str, ...
    'EdgeColor', 'none', 'Interpreter', 'latex', 'FontSize', FS, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

end


function [E, ang, got] = theta_plane(C, theta, phi, ph0, ph1, mag_scale)
% Full theta cut through z: |C| along theta (0..pi) at the phi sample
% nearest ph0 (returned at +theta) and at the phi sample nearest ph1 =
% ph0+pi (returned at -theta, in increasing angle order). phi is matched
% with wrap-around, so e.g. a cut at phi=180 finds a grid that stops at
% 355 deg or one that runs -180..180. got = [phi(i0) phi(i1)], the phi
% samples actually used (rad), for the caller's grid-mismatch check.
wrapd = @(a) abs(angle(exp(1j*a)));        % angular distance in [0, pi]
[~, i0] = min(wrapd(phi - ph0));
[~, i1] = min(wrapd(phi - ph1));
got = [phi(i0) phi(i1)];
theta = theta(:);
Ep = mag_scale * abs(C(:, i0));
En = mag_scale * abs(C(:, i1));
ang = [-flipud(theta); theta];
E   = [flipud(En); Ep];
end


function lines = mismatch_lines(src, sym, want, got, tol_deg)
% One LaTeX text line per requested angle whose nearest grid sample (got)
% is more than tol_deg away (with wrap-around); empty cell if all are OK.
% want/got in rad; sym is the LaTeX angle symbol, e.g. '\phi'.
lines = {};
for k = 1:numel(want)
    d = abs(angle(exp(1j*(got(k) - want(k))))) * 180/pi;
    if d > tol_deg
        lines{end+1} = sprintf('%s: nearest $%s=%.1f^\\circ$ (want $%.0f^\\circ$)', ...
            src, sym, mod(got(k)*180/pi, 360), mod(want(k)*180/pi, 360)); %#ok<AGROW>
    end
end
end
