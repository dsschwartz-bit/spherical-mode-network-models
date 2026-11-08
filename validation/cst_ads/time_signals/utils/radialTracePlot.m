function [ax, h] = radialTracePlot(theta_deg, t, V, varargin)
%RADIALTRACEPLOT  Polar "spoke" plot of time signals.
%
%   Each signal is drawn along a radial spoke at angle theta. Distance from
%   the origin is time, and the voltage is drawn perpendicular to the spoke
%   (positive voltage = counter-clockwise side of the spoke).
%
%   radialTracePlot(theta_deg, t, V)
%     theta_deg : 1xK vector of spoke angles in degrees
%     t         : Nx1 time vector (shared by all signals)
%     V         : NxK matrix, one column per angle (KxN is also accepted)
%
%   Name/value options:
%     'VoltScale'    : radial units per volt. Default sizes the largest
%                      trace to ~45% of the gap between spokes at the rim.
%     'RadialOffset' : radius where t(1) sits (default 0). Set > 0 to open
%                      a hole in the middle so traces don't pile up there.
%     'Fill'         : shade between trace and spoke (default true)
%     'FillAlpha'    : shading transparency (default 0.25)
%     'Colors'       : Kx3 color matrix (default lines(K))
%     'NumRings'     : approximate number of time rings (default 5). Rings
%                      are placed at round times (multiples of 1, 2 or 5
%                      x 10^k in TimeUnit), so the labels are round numbers
%     'TimeUnit'     : label appended to ring values, e.g. 'ms' (default 's')
%     'TimeFactor'   : multiply t by this for ring labels, e.g. 1e3 (default 1)
%     'Parent'       : axes to draw into (default gca)
%     'AngleGrid'    : spacing in degrees of radial grid lines drawn all the
%                      way round, e.g. 30 (default [] = none). Grid lines
%                      that do not coincide with a data spoke are labelled.
%     'AngleLabel'   : function handle mapping a plot angle (deg) to its
%                      label text, used for spokes and grid lines
%                      (default: the angle with a degree sign)
%     'Interpreter'  : text interpreter for all labels, 'latex' (default),
%                      'tex' or 'none'. With 'latex', AngleLabel strings
%                      must be valid LaTeX (e.g. '$45^\circ$')
%     'Grid'         : draw rings, angle grid, spokes, labels and scale bar
%                      (default true). Set false to overlay further traces
%                      on an existing plot (pass the same 'VoltScale', and
%                      a t that starts at the same t(1)).
%     'LineStyle'    : trace line style (default '-')
%     'Marker'       : trace marker (default 'none'), e.g. 'o' with
%                      'LineStyle','none' for markers only
%     'MarkerSize'   : marker size (default 4)
%     'RingLabelAngle': plot angle (deg) along which the time-ring labels
%                      are written (default [] = halfway between the first
%                      two spokes)
%     'ScaleRef'     : voltage shown by the scale bar (default [] = a round
%                      number near half the largest |V|). Pass the same
%                      value to panels that share a 'VoltScale'.
%     'FontSize'     : font size of the ring and grid labels (default 8);
%                      spoke labels and the scale bar label are 2 larger
%     'Sector'       : [a1 a2] plot-angle range in degrees drawn (rings, angle
%                      grid, axes limits), e.g. [0 90] for a quarter or
%                      [0 180] for a half polar plot (default [0 360]). With a
%                      partial sector the ring labels sit just below the a1
%                      edge and the scale bar at the lower right.
%     'SweepLabel'   : label (e.g. '$\theta$') of a curved arrow outside the
%                      rim showing the swept angle (default [] = none)
%     'SweepDir'     : +1 counter-clockwise (default), -1 clockwise
%     'SweepAngle'   : plot angle (deg) at which the arrow is centred (15)
%
%   Outputs: ax, and h with fields traces, spokes, fills, rings, angleGrid,
%   scaleBar (graphics handles).
%
%   Run with no arguments for a demo:  radialTracePlot

if nargin == 0
    demo();
    return
end

p = inputParser;
addParameter(p, 'VoltScale',    []);
addParameter(p, 'RadialOffset', 0);
addParameter(p, 'Fill',         true);
addParameter(p, 'FillAlpha',    0.25);
addParameter(p, 'Colors',       []);
addParameter(p, 'NumRings',     5);
addParameter(p, 'TimeUnit',     's');
addParameter(p, 'TimeFactor',   1);
addParameter(p, 'Parent',       []);
addParameter(p, 'AngleGrid',    []);
addParameter(p, 'AngleLabel',   []);
addParameter(p, 'Interpreter',  'latex');
addParameter(p, 'Grid',         true);
addParameter(p, 'LineStyle',    '-');
addParameter(p, 'Marker',       'none');
addParameter(p, 'MarkerSize',   4);
addParameter(p, 'FontSize',     8);
addParameter(p, 'RingLabelAngle', []);
addParameter(p, 'ScaleRef',     []);
addParameter(p, 'SweepLabel',   []);     % e.g. '$\theta$': arrow outside the rim showing the swept angle
addParameter(p, 'SweepDir',     1);      % +1 counter-clockwise, -1 clockwise (direction of increasing angle)
addParameter(p, 'SweepAngle',   15);     % plot angle (deg) at which the sweep arrow is centred
addParameter(p, 'Sector',       [0 360]);
parse(p, varargin{:});
o = p.Results;
isLatex = strcmpi(o.Interpreter, 'latex');
if isempty(o.AngleLabel)
    if isLatex
        o.AngleLabel = @(a) sprintf('$%g^\\circ$', a);
    else
        o.AngleLabel = @(a) sprintf('%g%c', a, char(176));
    end
end
if isLatex, pad = '$\;$'; pad2 = '$\;\;$'; else, pad = ' '; pad2 = '  '; end   % leading space survives LaTeX
txtArgs = {'Interpreter', o.Interpreter};

% ---- shape checks -------------------------------------------------------
theta_deg = theta_deg(:).';
K = numel(theta_deg);
t = t(:);
N = numel(t);
if size(V,1) ~= N && size(V,2) == N
    V = V.';
end
assert(isequal(size(V), [N K]), 'V must be N x K (N = numel(t), K = numel(theta)).');
if isempty(o.Colors), o.Colors = lines(K); end

% ---- geometry -----------------------------------------------------------
r    = o.RadialOffset + (t - t(1));   % radius for each sample
r0   = r(1);
rmax = r(end);

if isempty(o.VoltScale)
    if K > 1
        th  = sort(mod(theta_deg(:), 360));
        dth = min(diff([th; th(1) + 360]));
    else
        dth = 90;
    end
    gapAtRim    = rmax * sind(min(dth, 90));
    vmax        = max(abs(V(:)));
    if vmax == 0, vmax = 1; end
    o.VoltScale = 0.45 * gapAtRim / vmax;
end

% ---- axes ---------------------------------------------------------------
ax = o.Parent;
if isempty(ax), ax = gca; end
holdState = ishold(ax);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');

gridCol = [0.85 0.85 0.85];
txtCol  = [0.35 0.35 0.35];
sec     = o.Sector;
fullC   = diff(sec) >= 360 - 1e-9;
phi     = linspace(sec(1), sec(2), 361) * pi/180;

% ---- time rings ---------------------------------------------------------
% rings at round times in display units (TimeFactor applied)
tsc  = [t(1) t(end)] * o.TimeFactor;
step = niceNumber((tsc(2) - tsc(1)) / o.NumRings);
tk   = (ceil(tsc(1)/step)*step : step : tsc(2)) / o.TimeFactor;
if o.RadialOffset == 0, tk = tk(tk > t(1)); end   % skip the zero-radius ring
rk = o.RadialOffset + (tk - t(1));

% put ring labels halfway between the first two spokes
if K > 1
    th   = sort(mod(theta_deg, 360));
    lblA = th(1) + min(diff([th th(1)+360])) / 2;
else
    lblA = theta_deg(1) + 45;
end
if ~isempty(o.RingLabelAngle), lblA = o.RingLabelAngle; end

h.rings     = gobjects(0, 1);
h.angleGrid = gobjects(0, 1);
h.scaleBar  = gobjects(0, 1);
if o.Grid
    h.rings = gobjects(numel(rk), 1);
    for i = 1:numel(rk)
        h.rings(i) = plot(ax, rk(i)*cos(phi), rk(i)*sin(phi), '-', 'Color', gridCol);
        tv = round(tk(i)*o.TimeFactor/step)*step;
        if isLatex, rlbl = sprintf('%s\\textbf{%g %s}', pad, tv, o.TimeUnit);
        else,       rlbl = sprintf('%s%g %s', pad, tv, o.TimeUnit); end
        if fullC
            text(ax, rk(i)*cosd(lblA), rk(i)*sind(lblA), rlbl, ...
                txtArgs{:}, 'Color', txtCol, 'FontSize', o.FontSize, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        else                                           % just outside the a1 edge of the sector
            if isLatex, rlbl = sprintf('\\textbf{%g %s}', tv, o.TimeUnit); else, rlbl = sprintf('%g %s', tv, o.TimeUnit); end
            text(ax, rk(i)*cosd(sec(1)) + 0.02*rmax*sind(sec(1)), rk(i)*sind(sec(1)) - 0.02*rmax*cosd(sec(1)), rlbl, ...
                txtArgs{:}, 'Color', txtCol, 'FontSize', o.FontSize, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Rotation', sec(1));
        end
    end

    % ---- angle grid -----------------------------------------------------
    if ~isempty(o.AngleGrid)
        if fullC, ga = 0:o.AngleGrid:360 - o.AngleGrid/2;
        else,     ga = ceil(sec(1)/o.AngleGrid)*o.AngleGrid : o.AngleGrid : sec(2); end
        h.angleGrid = gobjects(numel(ga), 1);
        for i = 1:numel(ga)
            a = ga(i);
            h.angleGrid(i) = plot(ax, [r0 rmax]*cosd(a), [r0 rmax]*sind(a), ...
                '-', 'Color', gridCol);
            onSpoke = any(abs(mod(a - theta_deg + 180, 360) - 180) < 1e-6);
            if ~onSpoke                 % data spokes get their own bold label
                rl = rmax * 1.08;
                text(ax, rl*cosd(a), rl*sind(a), o.AngleLabel(a), ...
                    txtArgs{:}, 'Color', txtCol, 'FontSize', o.FontSize, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end
    end
end

% ---- spokes, traces, fills ---------------------------------------------
h.spokes = gobjects(K, 1);
h.traces = gobjects(K, 1);
h.fills  = gobjects(K, 1);

for k = 1:K
    c = cosd(theta_deg(k));
    s = sind(theta_deg(k));

    % baseline point on the spoke, and offset along the normal (-s, c)
    bx = r * c;
    by = r * s;
    d  = o.VoltScale * V(:, k);
    x  = bx - d * s;
    y  = by + d * c;

    if o.Fill
        h.fills(k) = patch(ax, [x; flipud(bx)], [y; flipud(by)], o.Colors(k,:), ...
            'FaceAlpha', o.FillAlpha, 'EdgeColor', 'none');
    end
    if o.Grid
        h.spokes(k) = plot(ax, [r0 rmax]*c, [r0 rmax]*s, '-', 'Color', [0.5 0.5 0.5]);
    end
    h.traces(k) = plot(ax, x, y, 'LineStyle', o.LineStyle, 'Marker', o.Marker, ...
        'MarkerSize', o.MarkerSize, 'Color', o.Colors(k,:), 'LineWidth', 1.4);

    % angle label just past the rim
    if o.Grid
        rl = rmax * 1.08;
        lbl = o.AngleLabel(theta_deg(k));
        if isLatex, lbl = ['\boldmath ' lbl]; end
        text(ax, rl*c, rl*s, lbl, txtArgs{:}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', 'FontSize', o.FontSize + 2);
    end
end

if o.Grid
    % ---- voltage scale bar ----------------------------------------------
    % Axes are equal-aspect, so a vertical bar of length VoltScale*vref is
    % the correct voltage length for every spoke orientation.
    if isempty(o.ScaleRef), vref = niceNumber(max(abs(V(:))) / 2); else, vref = o.ScaleRef; end
    L    = o.VoltScale * vref;
    if fullC
        xb = rmax * 1.15;  yb = -rmax * 1.05;
    else                                               % sector extents (with the origin)
        as  = [linspace(sec(1), sec(2), 181) sec(1) sec(2)];
        sx  = [0 rmax*cosd(as)];  sy = [0 rmax*sind(as)];
        xb  = max(sx) + 0.12*rmax;  yb = min(sy) + 0.02*rmax;
    end
    h.scaleBar = plot(ax, [xb xb], [yb yb + L], 'k-', 'LineWidth', 2);
    text(ax, xb, yb + L/2, sprintf('%s%g V', pad2, vref), txtArgs{:}, ...
        'VerticalAlignment', 'middle', ...
        'FontSize', o.FontSize + 2);

    % ---- swept-angle indicator (curved arrow + label outside the rim) ----
    if ~isempty(o.SweepLabel)
        a0 = o.SweepAngle;  span = 16;  ra = rmax * 1.17;  sd = sign(o.SweepDir);
        aa = linspace(a0 - sd*span/2, a0 + sd*span/2, 40);
        h.sweep = plot(ax, ra*cosd(aa), ra*sind(aa), 'k-', 'LineWidth', 1.3, 'Clipping', 'off');
        ae = aa(end);  tx = -sind(ae)*sd;  ty = cosd(ae)*sd;  nx = cosd(ae);  ny = sind(ae);
        hl = 0.05*rmax;  hw = 0.022*rmax;  px = ra*nx;  py = ra*ny;
        patch(ax, [px, px - hl*tx + hw*nx, px - hl*tx - hw*nx], ...
                  [py, py - hl*ty + hw*ny, py - hl*ty - hw*ny], 'k', 'EdgeColor', 'none', 'Clipping', 'off');
        slbl = o.SweepLabel;  if isLatex, slbl = ['\boldmath ' slbl]; end
        rl = rmax * 1.27;
        text(ax, rl*cosd(a0), rl*sind(a0), slbl, txtArgs{:}, 'FontSize', o.FontSize + 4, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Clipping', 'off');
    end

    if fullC
        lim = rmax * 1.25;
        xlim(ax, [-lim lim]);
        ylim(ax, [-lim lim]);
    else
        xlim(ax, [min(sx) - 0.14*rmax, max(sx) + 0.3*rmax]);
        ylim(ax, [min(sy) - 0.14*rmax, max(sy) + 0.2*rmax]);
    end
end

if ~holdState, hold(ax, 'off'); end
if nargout == 0, clear ax h; end
end

% ========================================================================
function v = niceNumber(x)
% round x to 1, 2 or 5 times a power of ten
if x <= 0, v = 1; return, end
e = 10^floor(log10(x));
m = x / e;
if m < 1.5,     v = 1*e;
elseif m < 3.5, v = 2*e;
elseif m < 7.5, v = 5*e;
else,           v = 10*e;
end
end

% ========================================================================
function demo()
% 8 coarse angles, each with a damped pulse that arrives later as the
% angle increases.
theta = 0:45:315;
t     = linspace(0, 10e-3, 600).';          % 10 ms
V     = zeros(numel(t), numel(theta));
for k = 1:numel(theta)
    t0     = 2e-3 + 5e-3 * k / numel(theta);  % arrival time
    tau    = 0.8e-3;
    f      = 1.5e3;
    env    = exp(-((t - t0) / tau).^2);
    V(:,k) = (0.5 + 0.5*cosd(theta(k))) * env .* sin(2*pi*f*(t - t0)) ...
             + 0.02 * randn(size(t));
end

figure('Color', 'w');
radialTracePlot(theta, t, V, 'RadialOffset', 1e-3, ...
    'TimeUnit', 'ms', 'TimeFactor', 1e3);
title('Voltage vs. time along each spoke');
end
