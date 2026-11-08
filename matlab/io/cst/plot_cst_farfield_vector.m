function fignum = plot_cst_farfield_vector(fignum, r, theta, phi, field, ...
                                            n_theta_ds, n_phi_ds, ...
                                            vec_scale, snapshot_mode, field_label, title_str, FS, LW)
%PLOT_CST_FARFIELD_VECTOR  3D cone plot of one far-field quantity (E or H)
% on a sphere of radius r, at a single time instant.
%
% Fixed-length cones (no stems) are drawn at down-sampled (theta, phi)
% points on a grey reference sphere; cone colour (white-yellow-orange-red)
% encodes the vector magnitude of the plotted real field.
% Single-source counterpart of plot_field_comparison_vector.
%
%   r             : plotting radius (only scales the picture; 1 is natural
%                   for r*E far-field data)
%   theta, phi    : angular grid (rad)
%   field         : (numth x numph x 3) complex, pages Vr/Vth/Vph (Vr = 0
%                   for a far field); note this order differs from
%                   read_cst_farfield's (phi, theta) pages
%   n_theta_ds, n_phi_ds : number of theta/phi samples kept for cones
%                   (about 21 x 21 is readable)
%   vec_scale     : multiplier on the field (e.g. eta0 to show H in V/m)
%   snapshot_mode : 'real'    -> t = 0          -> Re{phasor}
%                   'negimag' -> omega*t = pi/2 -> -Im{phasor}
%                   (exp(+j*omega*t) time convention)
%   field_label   : LaTeX field name, e.g. '$\vec{E}$'
%   title_str     : subtitle (string or cell array)
%   FS, LW        : font size, line width
%
% See also plot_field_comparison_vector, plot_cst_farfield_cuts.

theta = theta(:);      % numth x 1
phi   = phi(:).';      % 1 x numph
numth = numel(theta);
numph = numel(phi);

% --- Down-sampled grid for cone glyphs ---
th_idx = unique(round(linspace(1, numth, min(n_theta_ds, numth))));
ph_idx = unique(round(linspace(1, numph, min(n_phi_ds,   numph))));
[PHI_grid, TH_grid] = meshgrid(phi(ph_idx), theta(th_idx));   % (numth_ds x numph_ds)

[X, Y, Z] = local_sph_pos(r, TH_grid, PHI_grid);

Vr  = local_time_snapshot(vec_scale*field(th_idx,ph_idx,1), snapshot_mode);
Vth = local_time_snapshot(vec_scale*field(th_idx,ph_idx,2), snapshot_mode);
Vph = local_time_snapshot(vec_scale*field(th_idx,ph_idx,3), snapshot_mode);
[Vx,Vy,Vz] = local_sph_vec2cart(TH_grid, PHI_grid, Vr, Vth, Vph);
magGrid = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
maxmag = max(magGrid(:), [], 'omitnan');
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
ax = gca; hold(ax, "on")

% Opaque gray reference sphere with visible grid lines
surf(ax, Xs, Ys, Zs, "FaceColor", [0.85 0.85 0.85], "FaceAlpha", 0.9, ...
    "EdgeColor", [0.55 0.55 0.55], "LineWidth", 0.5, "HandleVisibility", "off")

for i = 1:numel(X)
    if ~isfinite(Vx(i)) || ~isfinite(Vy(i)) || ~isfinite(Vz(i)) || magGrid(i) == 0
        continue
    end
    colorIdx = max(1, min(size(cmap,1), round((magGrid(i)/maxmag) * (size(cmap,1)-1)) + 1));
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
xlabel(ax, "$x$", "FontSize", FS-2, "Interpreter", "latex")
ylabel(ax, "$y$", "FontSize", FS-2, "Interpreter", "latex")
zlabel(ax, "$z$", "FontSize", FS-2, "Interpreter", "latex")

cb = colorbar(ax);
cb.Label.Interpreter = "latex";
cb.Label.String = magnitude_label(field_label, snapshot_mode);
cb.Label.FontSize = FS;

title(ax, field_label + " --- " + string(snapshot_label(snapshot_mode)), "FontSize", FS, "Interpreter", "latex")
subtitle(ax, title_str, "FontSize", FS-2, "Interpreter", "latex")

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
% Extract the physical (real, time-domain) field at a given instant from
% a complex phasor V, assuming the +j*omega*t time convention used
% throughout this toolbox: v(t) = Re{ V * exp(j*omega*t) }.
%   mode = 'real'    -> t = 0            -> Re{V}
%   mode = 'negimag' -> omega*t = pi/2   -> Re{V*exp(j*pi/2)} = -Im{V}
switch mode
    case 'real'
        Vout = real(V);
    case 'negimag'
        Vout = -imag(V);
    otherwise
        error('plot_cst_farfield_vector: unknown snapshot mode "%s"', mode)
end
end

function lbl = snapshot_label(mode)
switch mode
    case 'real'
        lbl = "$t=0$ (real part)";
    case 'negimag'
        lbl = "$\omega t=\pi/2$ ($-\mathrm{Im}$ of phasor)";
    otherwise
        error('plot_cst_farfield_vector: unknown snapshot mode "%s"', mode)
end
end

function lbl = magnitude_label(field_label, mode)
% MAGNITUDE_LABEL  Colour-bar label: the norm of the real snapshot,
% ||Re[F]|| or ||Im[F]|| (= ||-Im[F]||), not the phasor magnitude. The "$"
% delimiters of field_label are stripped so the label is one math block.
field_symbol = strrep(field_label, "$", "");
switch mode
    case 'real'
        lbl = "$||\mathrm{Re}[" + field_symbol + "]||$";
    case 'negimag'
        lbl = "$||\mathrm{Im}[" + field_symbol + "]||$";
    otherwise
        error('plot_cst_farfield_vector: unknown snapshot mode "%s"', mode)
end
end

function local_plotCone(ax, origin, direction, coneLength, baseRadius, faceColor)
% LOCAL_PLOTCONE  Draws a single solid 3D cone glyph on axes ax: base at
% 'origin', pointing along 'direction' (any nonzero length; only its
% direction is used), extending for 'coneLength', tapering from
% 'baseRadius' at the base to a point at the tip.
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
% LOCAL_WHITEYELLOWORANGERED  Sequential colormap: white (small) ->
% yellow -> orange -> red (large).
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
if size(cmap,1) < n
    cmap = [cmap; repmat(cmap(end,:), n-size(cmap,1), 1)];
elseif size(cmap,1) > n
    cmap = cmap(1:n,:);
end
end