function fignum = plot_field_comparison_vector(fignum, r, theta, phi, ...
                                                field_psi, field_Q, field_computed, ...
                                                n_theta_ds, n_phi_ds, ...
                                                vec_scale, snapshot_mode, field_label, title_str, FS, LW)
%PLOT_FIELD_COMPARISON_VECTOR  Cone-glyph plot of a vector field on a sphere, one tile per source.
%   fignum = plot_field_comparison_vector(fignum, r, theta, phi, field_psi,
%       field_Q, field_computed, n_theta_ds, n_phi_ds, vec_scale,
%       snapshot_mode, field_label, title_str, FS, LW)
%   Draws the instantaneous field at radius r on a down-sampled
%   (theta,phi) grid as fixed-length cones (no stems) pointing along the
%   field, coloured by magnitude (white -> yellow -> orange -> red) on a
%   colour scale shared by all tiles, over a grey reference sphere. Each
%   source gets its own tile of a 1 x n tiledlayout.
%
%   theta (numth x 1), phi (1 x numph) : full angular grids, radians
%   field_psi, field_Q, field_computed : (numth x numph x 3) complex
%       phasors, pages = [V_r, V_theta, V_phi]: reconstruction from psi,
%       reconstruction from FEKO Q ([] to omit), simulator field
%   n_theta_ds, n_phi_ds : number of theta/phi samples kept for the cones
%       (about 21 x 21 is readable)
%   vec_scale     : factor applied to all fields (e.g. eta0 for H)
%   snapshot_mode : 'real' (t = 0, Re{V}) or 'negimag' (omega t = pi/2,
%       -Im{V}), for the exp(+j omega t) convention
%   field_label   : LaTeX name of the field, e.g. '$\vec{E}$'
%   title_str     : subtitle line(s) (string or cell array)
%   LW is not used.
%
%   See also plot_field_comparison, plot_field_comparison_polar.

theta = theta(:);      % numth x 1
phi   = phi(:).';      % 1 x numph
numth = numel(theta);
numph = numel(phi);

src_names = {"Reconstructed $\psi$", "Reconstructed $Q$", "Computed"};
src_data  = {field_psi, field_Q, field_computed};

keep      = ~cellfun(@isempty, src_data);
src_names = src_names(keep);
src_data  = src_data(keep);
n_src     = numel(src_data);

% --- Down-sampled grid for cone glyphs ---
th_idx = unique(round(linspace(1, numth, min(n_theta_ds, numth))));
ph_idx = unique(round(linspace(1, numph, min(n_phi_ds,   numph))));
[PHI_grid, TH_grid] = meshgrid(phi(ph_idx), theta(th_idx));   % (numth_ds x numph_ds)

[X, Y, Z] = local_sph_pos(r, TH_grid, PHI_grid);

% --- Cartesian components on the down-sampled grid and a common magnitude scale ---
Vxyz    = cell(n_src,1);
magGrid = cell(n_src,1);
maxmag  = 0;
for si = 1:n_src
    Vr  = local_time_snapshot(vec_scale*src_data{si}(th_idx,ph_idx,1), snapshot_mode);
    Vth = local_time_snapshot(vec_scale*src_data{si}(th_idx,ph_idx,2), snapshot_mode);
    Vph = local_time_snapshot(vec_scale*src_data{si}(th_idx,ph_idx,3), snapshot_mode);
    [Vx,Vy,Vz] = local_sph_vec2cart(TH_grid, PHI_grid, Vr, Vth, Vph);
    Vxyz{si}    = cat(3, Vx, Vy, Vz);
    magGrid{si} = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
    maxmag = max([maxmag; magGrid{si}(:)], [], 'omitnan');
end
if ~(maxmag > 0) || isnan(maxmag)
    maxmag = 1;
end

coneLength      = 0.15*r;   % fixed cone length -- does not scale with magnitude
baseRadiusRatio = 0.15;     % cone base radius, as a fraction of its length
lim_val         = 1.5*r;
cmap            = local_whiteYellowOrangeRed(256);

[Xs, Ys, Zs] = sphere(20);
Xs = r*Xs; Ys = r*Ys; Zs = r*Zs;

fignum = fignum + 1;
figure(fignum);
tl = tiledlayout(1, n_src, "TileSpacing", "compact", "Padding", "compact");

for si = 1:n_src
    ax = nexttile; hold(ax, "on")

    % Grey reference sphere
    surf(ax, Xs, Ys, Zs, "FaceColor", [0.85 0.85 0.85], "FaceAlpha", 0.9, ...
        "EdgeColor", [0.55 0.55 0.55], "LineWidth", 0.5, "HandleVisibility", "off")

    Vx = Vxyz{si}(:,:,1); Vy = Vxyz{si}(:,:,2); Vz = Vxyz{si}(:,:,3);
    magv = magGrid{si};
    for i = 1:numel(X)
        if ~isfinite(Vx(i)) || ~isfinite(Vy(i)) || ~isfinite(Vz(i)) || magv(i) == 0
            continue
        end
        colorIdx = max(1, min(size(cmap,1), round((magv(i)/maxmag) * (size(cmap,1)-1)) + 1));
        coneColor = cmap(colorIdx,:);
        local_plotCone(ax, [X(i) Y(i) Z(i)], [Vx(i) Vy(i) Vz(i)], ...
            coneLength, coneLength*baseRadiusRatio, coneColor);
    end
    hold(ax, "off")

    colormap(ax, cmap);
    clim(ax, [0, maxmag]);

    axis(ax, "equal")
    xlim(ax, [-lim_val lim_val]); ylim(ax, [-lim_val lim_val]); zlim(ax, [-lim_val lim_val]);
    grid(ax, "on"); view(ax, 45, 30)
    xlabel(ax, "$x$ (m)", "FontSize", FS-4, "Interpreter", "latex")
    ylabel(ax, "$y$ (m)", "FontSize", FS-4, "Interpreter", "latex")
    zlabel(ax, "$z$ (m)", "FontSize", FS-4, "Interpreter", "latex")
    title(ax, src_names{si}, "FontSize", FS, "Interpreter", "latex")

    if si == n_src
        cb = colorbar(ax);
        cb.Layout.Tile = "east";
        cb.Label.Interpreter = "latex";
        cb.Label.String = field_label + " magnitude";
        cb.Label.FontSize = FS;
    end
end

title(tl, field_label + " --- " + string(snapshot_label(snapshot_mode)), "FontSize", FS, "Interpreter", "latex")
subtitle(tl, title_str, "FontSize", FS, "Interpreter", "latex")

end

% ------------------------------------------------------------------------
function [x,y,z] = local_sph_pos(r, th, ph)
% Position of a point on a sphere of radius r at colatitude th, azimuth ph.
x = r.*sin(th).*cos(ph);
y = r.*sin(th).*sin(ph);
z = r.*cos(th);
end

function [Vx,Vy,Vz] = local_sph_vec2cart(th, ph, Vr, Vth, Vph)
% Convert spherical vector components (Vr,Vth,Vph), evaluated at (th,ph),
% into Cartesian components using the standard spherical unit-vector basis:
%   r_hat  = [ sin(th)cos(ph),  sin(th)sin(ph),  cos(th)]
%   th_hat = [ cos(th)cos(ph),  cos(th)sin(ph), -sin(th)]
%   ph_hat = [-sin(ph),         cos(ph),         0      ]
Vx = Vr.*sin(th).*cos(ph) + Vth.*cos(th).*cos(ph) - Vph.*sin(ph);
Vy = Vr.*sin(th).*sin(ph) + Vth.*cos(th).*sin(ph) + Vph.*cos(ph);
Vz = Vr.*cos(th)          - Vth.*sin(th);
end

function Vout = local_time_snapshot(V, mode)
% Instantaneous field v(t) = Re{V exp(j omega t)} at the chosen instant.
%   mode = 'real'    -> t = 0            -> Re{V}
%   mode = 'negimag' -> omega*t = pi/2   -> Re{V*exp(j*pi/2)} = -Im{V}
switch mode
    case 'real'
        Vout = real(V);
    case 'negimag'
        Vout = -imag(V);
    otherwise
        error('plot_field_comparison_vector: unknown snapshot mode "%s"', mode)
end
end

function lbl = snapshot_label(mode)
switch mode
    case 'real'
        lbl = "$t=0$ (real part)";
    case 'negimag'
        lbl = "$\omega t=\pi/2$ ($-\mathrm{Im}$ of phasor)";
    otherwise
        error('plot_field_comparison_vector: unknown snapshot mode "%s"', mode)
end
end

function local_plotCone(ax, origin, direction, coneLength, baseRadius, faceColor)
% Solid cone with its base at origin, pointing along direction (only the
% direction is used), of length coneLength and base radius baseRadius.
% A +z cone from cylinder() is rotated (Rodrigues) and translated.
nsides = 12;
[cx, cy, cz] = cylinder([baseRadius 0], nsides);   % unit-height cone along +z, base at z=0
cz = cz * coneLength;

dirUnit = direction(:) / norm(direction);
zaxis = [0; 0; 1];
v = cross(zaxis, dirUnit);
s = norm(v);
c = dot(zaxis, dirUnit);

if s < 1e-10
    if c > 0
        R = eye(3);                % direction already along +z
    else
        R = diag([1, -1, -1]);     % direction along -z: 180-degree flip
    end
else
    vx = [   0    -v(3)   v(2);
           v(3)     0    -v(1);
          -v(2)    v(1)     0  ];
    R = eye(3) + vx + vx^2 * ((1 - c) / s^2);   % Rodrigues' rotation formula
end

pts = R * [cx(:)'; cy(:)'; cz(:)'];
cxr = reshape(pts(1,:), size(cx)) + origin(1);
cyr = reshape(pts(2,:), size(cy)) + origin(2);
czr = reshape(pts(3,:), size(cz)) + origin(3);

surf(ax, cxr, cyr, czr, "FaceColor", faceColor, "EdgeColor", "none", "HandleVisibility", "off")
end

function cmap = local_whiteYellowOrangeRed(n)
% Sequential n x 3 colormap white -> yellow -> orange -> red.
if nargin < 1
    n = 256;
end
keypoints = [1   1   1;    % white
             1   1   0;    % yellow
             1   0.5 0;    % orange
             1   0   0];   % red
nk = size(keypoints, 1);
segLen = floor(n / (nk-1));
cmap = [];
for k = 1:nk-1
    seg = [linspace(keypoints(k,1), keypoints(k+1,1), segLen)', ...
           linspace(keypoints(k,2), keypoints(k+1,2), segLen)', ...
           linspace(keypoints(k,3), keypoints(k+1,3), segLen)'];
    cmap = [cmap; seg]; %#ok<AGROW>
end
% Pad/trim to exactly n rows.
if size(cmap,1) < n
    cmap = [cmap; repmat(cmap(end,:), n-size(cmap,1), 1)];
elseif size(cmap,1) > n
    cmap = cmap(1:n,:);
end
end