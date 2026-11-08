function fignum = plot_cst_farfield_vector_compare(fignum, r, ...
        theta1, phi1, field1, label1, theta2, phi2, field2, label2, ...
        n_theta_ds, n_phi_ds, vec_scale, field_label, title_str, FS, LW)
% PLOT_CST_FARFIELD_VECTOR_COMPARE   2x2 3D vector (cone) plots of a far
% field from two sources (e.g. CST and ADS):
%
%                  real part (t=0)     -Imag part (omega*t = pi/2)
%     source 1          [1]                     [2]
%     source 2          [3]                     [4]
%
% The second column is the NEGATIVE of the imaginary part, i.e. the field
% a quarter period LATER (v(t) = Re{V exp(j w t)} at w t = pi/2 is -Im{V}).
% Each source stays on its own native angular grid (no interpolation), and
% all four tiles share ONE color scale (cone color = Cartesian vector
% magnitude of the plotted real-valued snapshot, in the field's units).
%
%   fignum               : last figure number used (incremented and returned)
%   r                    : sphere radius used for cone placement
%   theta*,phi*,field*   : angular grid and (numth x numph x 3) complex
%                          Vr/Vth/Vph field (Vr = 0 for a far field)
%   label1, label2       : LaTeX source names, e.g. "CST", "ADS"
%   n_theta_ds, n_phi_ds : cone down-sampling, applied to each source's grid
%   vec_scale            : multiplier on the fields before plotting
%   field_label          : LaTeX field name, e.g. '$\vec{E}$'
%   title_str            : overall figure title (e.g. the frequency)
%   FS, LW               : font size, line width
%
% See also plot_cst_farfield_vector.

thetas = {theta1(:), theta2(:)};
phis   = {phi1(:).', phi2(:).'};
fields = {field1, field2};
labels = {label1, label2};
modes  = {'real', 'negimag'};
mode_titles = {"real part", "$-$imag part"};

Xc = cell(2,2); Yc = cell(2,2); Zc = cell(2,2);
Vxc = cell(2,2); Vyc = cell(2,2); Vzc = cell(2,2); magc = cell(2,2);
maxmag = 0;
for s = 1:2
    numth = numel(thetas{s}); numph = numel(phis{s});
    th_idx = unique(round(linspace(1, numth, min(n_theta_ds, numth))));
    ph_idx = unique(round(linspace(1, numph, min(n_phi_ds,   numph))));
    [PHI_grid, TH_grid] = meshgrid(phis{s}(ph_idx), thetas{s}(th_idx));
    [X, Y, Z] = local_sph_pos(r, TH_grid, PHI_grid);
    for m = 1:2
        Vr  = local_time_snapshot(vec_scale*fields{s}(th_idx,ph_idx,1), modes{m});
        Vth = local_time_snapshot(vec_scale*fields{s}(th_idx,ph_idx,2), modes{m});
        Vph = local_time_snapshot(vec_scale*fields{s}(th_idx,ph_idx,3), modes{m});
        [Vx,Vy,Vz] = local_sph_vec2cart(TH_grid, PHI_grid, Vr, Vth, Vph);
        mag = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
        Xc{s,m}=X; Yc{s,m}=Y; Zc{s,m}=Z; Vxc{s,m}=Vx; Vyc{s,m}=Vy; Vzc{s,m}=Vz; magc{s,m}=mag;
        smax = max(mag(:), [], 'omitnan');
        if smax > maxmag && isfinite(smax)
            maxmag = smax;
        end
    end
end
if ~(maxmag > 0)
    maxmag = 1;
end

coneLength      = 0.2*r;
baseRadiusRatio = 0.25;
lim_val         = 1.25*r;   % tighter axes limits -> bigger sphere in each tile
cmap            = local_whiteYellowOrangeRed(256);

[Xs, Ys, Zs] = sphere(20);
Xs = r*Xs; Ys = r*Ys; Zs = r*Zs;

fignum = fignum + 1;
figure(fignum);
% The figure is left docked (MATLAB default: one window, one tab per figure);
% setting its Position would undock it into a separate window.
tl = tiledlayout(2, 2, "TileSpacing", "tight", "Padding", "tight");

for s = 1:2
    for m = 1:2
        ax = nexttile; hold(ax, "on")
        surf(ax, Xs, Ys, Zs, "FaceColor", [0.85 0.85 0.85], "FaceAlpha", 0.9, ...
            "EdgeColor", [0.55 0.55 0.55], "LineWidth", 0.5, "HandleVisibility", "off")

        X=Xc{s,m}; Y=Yc{s,m}; Z=Zc{s,m}; Vx=Vxc{s,m}; Vy=Vyc{s,m}; Vz=Vzc{s,m}; mag=magc{s,m};
        for i = 1:numel(X)
            if ~isfinite(Vx(i)) || ~isfinite(Vy(i)) || ~isfinite(Vz(i)) || mag(i) == 0
                continue
            end
            colorIdx = max(1, min(size(cmap,1), round((mag(i)/maxmag) * (size(cmap,1)-1)) + 1));
            local_plotCone(ax, [X(i) Y(i) Z(i)], [Vx(i) Vy(i) Vz(i)], ...
                coneLength, coneLength*baseRadiusRatio, cmap(colorIdx,:));
        end
        hold(ax, "off")

        colormap(ax, cmap);
        clim(ax, [0, maxmag]);
        axis(ax, "equal")
        xlim(ax, [-lim_val lim_val]); ylim(ax, [-lim_val lim_val]); zlim(ax, [-lim_val lim_val]);
        grid(ax, "on"); view(ax, 45, 30)
        xlabel(ax, "$x$", "FontSize", FS, "Interpreter", "latex")
        ylabel(ax, "$y$", "FontSize", FS, "Interpreter", "latex")
        zlabel(ax, "$z$", "FontSize", FS, "Interpreter", "latex")
        ax.FontSize = FS;
        title(ax, labels{s} + ": " + mode_titles{m}, "FontSize", FS, "Interpreter", "latex")

        if s == 2 && m == 2
            cb = colorbar(ax);
            cb.Layout.Tile = "east";
            cb.Label.Interpreter = "latex";
            cb.Label.String = "Cone color: vector magnitude of plotted " + field_label;
            cb.Label.FontSize = FS;
            cb.FontSize = FS;
        end
    end
end

title(tl, title_str, "FontSize", FS, "Interpreter", "latex")

end

% ------------------------------------------------------------------------
function [x,y,z] = local_sph_pos(r, th, ph)
x = r.*sin(th).*cos(ph);
y = r.*sin(th).*sin(ph);
z = r.*cos(th);
end

function [Vx,Vy,Vz] = local_sph_vec2cart(th, ph, Vr, Vth, Vph)
Vx = Vr.*sin(th).*cos(ph) + Vth.*cos(th).*cos(ph) - Vph.*sin(ph);
Vy = Vr.*sin(th).*sin(ph) + Vth.*cos(th).*sin(ph) + Vph.*cos(ph);
Vz = Vr.*cos(th)          - Vth.*sin(th);
end

function Vout = local_time_snapshot(V, mode)
switch mode
    case 'real'
        Vout = real(V);
    case 'negimag'
        Vout = -imag(V);
    otherwise
        error('plot_cst_farfield_vector_compare: unknown snapshot mode "%s"', mode)
end
end

function local_plotCone(ax, origin, direction, coneLength, baseRadius, faceColor)
nsides = 12;
[cx, cy, cz] = cylinder([baseRadius 0], nsides);
cz = cz * coneLength;

dirUnit = direction(:) / norm(direction);
zaxis = [0; 0; 1];
v = cross(zaxis, dirUnit);
s = norm(v);
c = dot(zaxis, dirUnit);

if s < 1e-10
    if c > 0
        R = eye(3);
    else
        R = diag([1, -1, -1]);
    end
else
    vx = [   0    -v(3)   v(2);
           v(3)     0    -v(1);
          -v(2)    v(1)     0  ];
    R = eye(3) + vx + vx^2 * ((1 - c) / s^2);
end

pts = R * [cx(:)'; cy(:)'; cz(:)'];
cxr = reshape(pts(1,:), size(cx)) + origin(1);
cyr = reshape(pts(2,:), size(cy)) + origin(2);
czr = reshape(pts(3,:), size(cz)) + origin(3);

surf(ax, cxr, cyr, czr, "FaceColor", faceColor, "EdgeColor", "none", "HandleVisibility", "off")
end

function cmap = local_whiteYellowOrangeRed(n)
if nargin < 1
    n = 256;
end
keypoints = [1   1   1;
             1   1   0;
             1   0.5 0;
             1   0   0];
nk = size(keypoints, 1);
segLen = floor(n / (nk-1));
cmap = [];
for k = 1:nk-1
    seg = [linspace(keypoints(k,1), keypoints(k+1,1), segLen)', ...
           linspace(keypoints(k,2), keypoints(k+1,2), segLen)', ...
           linspace(keypoints(k,3), keypoints(k+1,3), segLen)'];
    cmap = [cmap; seg]; %#ok<AGROW>
end
if size(cmap,1) < n
    cmap = [cmap; repmat(cmap(end,:), n-size(cmap,1), 1)];
elseif size(cmap,1) > n
    cmap = cmap(1:n,:);
end
end