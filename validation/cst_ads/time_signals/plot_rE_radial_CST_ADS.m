%PLOT_RE_RADIAL_CST_ADS  Radial (spoke) plots of CST and ADS radiated far-field time signals.
% Radiated far-field time signals r*E_theta(t) and r*E_phi(t) from CST and
% ADS, drawn with radialTracePlot in a 1 x 2 tiled layout (one panel per
% component): one spoke per observation direction in a cut, time along the
% spoke, voltage across it. CST solid with shading, ADS dotted; radial grid
% lines every angle_grid degrees. Both panels share one volt scale, so the
% two components can be compared directly.
%
% The files may hold many (theta, phi) directions; the cut set below picks
% which of them are drawn:
%  cut = "theta" : plane through the z axis at phi = phi_cut. This plane
%            also holds the half-plane phi = phi_cut + 180 (wrapped to
%            0..360), which is drawn on the left as negative theta: a probe
%            at (theta, phi_cut+180) sits at -theta. z is up (plot angle
%            90 - theta). The field is plotted as given (no sign change).
%  cut = "phi"   : cone theta = theta_cut (the xy plane for 90). Spokes at
%            phi over the full 0..360, phi = 0 (x axis) to the right,
%            increasing counter-clockwise.
% Directions not in the cut are skipped. With spokes = "any" a spoke is drawn
% for every direction in the cut found in either file (a source that lacks it
% has no trace on that spoke); with spokes = "both" only directions present in
% both files are drawn. Errors are computed where both files have data.
%
% CST : farfield probe export, blocks (Abs)/(Phi)/(Theta) per probe and per
%       excitation [n], header "(Spherical) (theta phi r)" with r in the CST
%       length unit (mm). Only the blocks of excitation cst_excitation are
%       used (readCstFarfieldProbes.m).
% ADS : Data Display CSV with columns time [s], Vtheta_prime, Vphi_prime
%       (r*E in V) and theta_deg [deg] (or theta [rad]), optionally phi_deg
%       / phi. Without a phi column the sweep is taken at phi = ads_phi
%       (readAdsFarfieldSweeps.m).
%
% Errors (where both have the direction), relative to CST over the overlap of
% the records: relative MSE = integral (ADS - CST)^2 dt / integral CST^2 dt,
% per component and re the CST total field. Components that neither file
% holds are dropped from a cut's panels. Output: figures and console tables.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(fullfile(here, 'utils')); % radialTracePlot

%% ===== USER SETTINGS =====
cst_root = "C:\path\to\CST\";                    % folder containing the CST project folders
ads_dir  = "C:\path\to\ADS\workspace\data\";   % ADS workspace data folder (CSV exports)
% =========================

%% ---- antenna, excitation and folders ------------------------------------
antenna_name = "rectangularpatch";
cst_variant  = "symmetry_FDTD_2.5GHz";   % CST model / solver variant
excitation   = "gaussian";
band         = "2GHz_3GHz";         % excitation band, as in the exported file names
port         = 1;                    % excited port
cst_project  = antenna_name + "_" + cst_variant;   % CST project folder name (some projects
                                     % add "_" + excitation + "_<fc>")
cst_subdir = "port2_SC";            % subfolder of Export\time_signals for this run ("" = none)
cst_dir  = string(fullfile(cst_root, cst_project, "Export", "time_signals", cst_subdir)) + filesep;

%% ---- files and settings ---------------------------------------------------
ads_file = string(fullfile(ads_dir, antenna_name + "_SC_rE_" + excitation + "_" + band + ".csv"));
cst_file = cst_dir + ["ff_Etheta_" + excitation + "_port" + port + "_" + band + ".txt" ...
                      "ff_Ephi_"   + excitation + "_port" + port + "_" + band + ".txt"];   % [1] Theta / Phi blocks
cst_file = cst_file(isfile(cst_file));        % whichever components were exported
if isempty(cst_file)                          % else any far-field export ff_E*_<band>.txt
    L = dir(cst_dir + "ff_E*_" + band + ".txt");
    cst_file = string(fullfile({L.folder}, {L.name}));
    if isempty(cst_file), error('No CST far-field export ff_E*_%s.txt in %s', band, cst_dir); end
    fprintf('CST far-field file(s): %s\n', strjoin(string({L.name}), ', '));
end
% CST ASCII exports hold at most 25 curves: split many probes over several
% exports and list them all, e.g. cst_file = ["...part1.txt" "...part2.txt"]
model_name = "";          % plot title; "" = built from excitation, cst_subdir (port termination) and band
if model_name == ""
    if contains(excitation, "cosim"), exTxt = "co-simulation";
    else, exTxt = upper(extractBefore(excitation, 2)) + extractAfter(excitation, 1) + " excitation"; end
    parts = ["Rectangular patch" exTxt];
    tk = regexp(cst_subdir, 'port(\d+)_(OC|SC)', 'tokens', 'once');
    if ~isempty(tk)
        term = struct('OC', "open", 'SC', "shorted");
        parts(end+1) = "port " + tk{1} + " " + term.(tk{2});
    end
    parts(end+1) = regexprep(band, '([\d.]+)GHz_([\d.]+)GHz', '$1--$2 GHz');
    model_name = strjoin(parts, ", ");
end

% ---- cuts to plot: one figure each ----
%  cut "theta": plane phi = angle (+ angle + 180 drawn as negative theta)
%  cut "phi"  : cone theta = angle, phi 0..360
%  components : panels of that figure (CST block, ADS V<comp>_prime)
cuts = struct( ...
    'cut',        {"theta",                    "phi"}, ...
    'angle',      {90,                          45}, ... % theta cut: phi = angle; phi cut: theta = angle
    'components', {["Theta" "Phi"],            ["Theta" "Phi"]}, ...
    'only',       {[], []});   % spokes to keep (cut coordinate, deg); [] = all
spokes    = "both";       % "both": only directions present in both files (for comparison);
                          %         a cut with none is skipped with a message
                          % "any" : every direction in the cut found in either file
                          %         (a source without it has no trace there)

allComp        = ["Theta" "Phi"];   % components read from the files
cst_excitation = 1;       % CST excitation [n] to use (ADS excites port 1 only)
ads_theta      = 90;      % deg, theta of the ADS sweep if the file has no theta column
ads_phi        = 180;     % deg, phi of the ADS sweep if the file has no phi column
ring_label_angle = 255;   % plot angle (deg) of the time-ring labels; pick an empty sector
angle_grid = 30;          % radial grid line spacing [deg]
font_size  = 12;          % ring/grid labels; spoke labels, scale bar, legend
                          % and titles are 2 larger
fig_size   = [1300 760];  % figure width and height [pixels] (not used: figures open docked as tabs)
t_start    = [];          % ns, time at the inner end of each spoke; [] = (r - a)/c0 - 0.5 ns
                          % (before (r - a)/c0 the field is zero, so this drops the
                          % empty early part and spreads the signals out)
a_enclosing = [];         % m, enclosing-sphere radius; [] = parameter a of the CST export
radial_hole = 0.12;       % fraction of the radius left empty at the centre so
                          % the traces don't crowd together there
cst_length_unit = 1e-3;   % CST length unit [m]
angle_tol = 0.1;          % deg, direction matching

c0 = 299792458;

%% ---- read both files --------------------------------------------------------
if isempty(a_enclosing)                                % enclosing-sphere radius from CST
    hdr = readlines(cst_file(1));
    tk  = regexp(hdr(find(startsWith(strtrim(hdr), "#Parameters"), 1)), '[{;]\s*a\s*=\s*([-+0-9.eE]+)', 'tokens', 'once');
    if isempty(tk), a_enclosing = 0; else, a_enclosing = str2double(tk{1}) * cst_length_unit; end
end
cst = readCstFarfieldProbes(cst_file, allComp, cst_excitation, cst_length_unit);
ads = readAdsFarfieldSweeps(ads_file, allComp, ads_phi);
for k = 1:numel(ads)                                   % no theta column: use ads_theta
    if isnan(ads(k).theta), ads(k).theta = ads_theta; end
end

for ic = 1:numel(cuts)                                 % ===== one figure per cut =====
cut = cuts(ic).cut;  components = cuts(ic).components;
% keep a component only if both CST and ADS hold it
hasC = arrayfun(@(s) any(arrayfun(@(p) ~isempty(p.e{allComp == s}), cst)), components);
hasA = arrayfun(@(s) any(arrayfun(@(p) ~all(isnan(p.e(:, allComp == s))), ads)), components);
srcs = ["CST" "ADS"];
for c = find(~(hasC & hasA)), fprintf('E_%s not in the %s data: panel dropped\n', lower(components(c)), strjoin(srcs([~hasC(c) ~hasA(c)]), ' or ')); end
hasC = hasC & hasA;
components = components(hasC);  nC = numel(components);
if nC == 0, fprintf('no component of this cut is in the files; figure skipped\n'); continue, end
cIdx = arrayfun(@(s) find(allComp == s), components);  % index into the read components
if cut == "theta", phi_cut = cuts(ic).angle; theta_cut = NaN;
else,               theta_cut = cuts(ic).angle; phi_cut = NaN; end
inCut = @(th, ph) cutCoord(th, ph, cut, phi_cut, theta_cut, angle_tol);
if cut == "theta", cutTxt = sprintf('theta cut at phi = %g/%g deg', phi_cut, mod(phi_cut + 180, 360));
else,               cutTxt = sprintf('phi cut at theta = %g deg', theta_cut); end
fprintf('\n===== %s (%s) =====\n', cutTxt, strjoin("E_" + lower(components), ', '));

% cut coordinate (signed theta, or phi) of every direction; NaN = not in the cut
uC = arrayfun(@(p) inCut(p.theta, p.phi), cst);
uA = arrayfun(@(p) inCut(p.theta, p.phi), ads);

u = sort(uniquetol([uC(~isnan(uC)) uA(~isnan(uA))], angle_tol, 'DataScale', 1));
if cut == "theta", u(abs(u + 180) < angle_tol) = 180; u = unique(u); end
if isempty(u), fprintf('no CST probe or ADS sweep lies in this cut; figure skipped\n'); continue, end
K = numel(u);
jC = zeros(1, K);  jA = zeros(1, K);      % index into cst / ads, 0 = none
for k = 1:K
    m = find(abs(uC - u(k)) < angle_tol, 1);  if ~isempty(m), jC(k) = m; end
    m = find(abs(uA - u(k)) < angle_tol, 1);  if ~isempty(m), jA(k) = m; end
end

if cut == "theta"
    cutName = sprintf('held constant: $\\phi = %g^\\circ$ (right half) and $\\phi = %g^\\circ$ (left half, drawn as negative $\\theta$)', ...
        phi_cut, mod(phi_cut + 180, 360));
    sweepName = 'spoke angle: $\theta$ from the $z$ axis (up); time along each spoke';  axName = 'theta';
    panelCut = sprintf('$\\phi = %g^\\circ$ / $%g^\\circ$', phi_cut, mod(phi_cut + 180, 360));
    ang = 90 - u;                                     % z (theta = 0) up
    angLabel = @(a) sprintf('$%g^\\circ$', mod(90 - a + 180, 360) - 180);
else
    cutName = sprintf('held constant: $\\theta = %g^\\circ$', theta_cut);
    sweepName = 'spoke angle: $\phi$ from the $x$ axis; time along each spoke';  axName = 'phi';
    panelCut = sprintf('$\\theta = %g^\\circ$', theta_cut);
    ang = u;                                          % phi = 0 to the right
    angLabel = @(a) sprintf('$%g^\\circ$', mod(a, 360));
end
yn = '-x';
fprintf('\ndirections in the cut (x = data present)\n%7s  CST  ADS\n', axName);
for k = 1:K
    fprintf('%7g   %c    %c\n', u(k), yn((jC(k) > 0) + 1), yn((jA(k) > 0) + 1));
end
% optionally keep only the directions both files have
if spokes == "both"
    both = jC > 0 & jA > 0;
    if ~any(both)
        lC = strjoin(compose('%g', u(jC > 0)), ', ');  if strlength(lC) == 0, lC = "none"; else, lC = lC + " deg"; end
        lA = strjoin(compose('%g', u(jA > 0)), ', ');  if strlength(lA) == 0, lA = "none"; else, lA = lA + " deg"; end
        fprintf(['no direction in this cut is in both files; figure skipped\n' ...
            '  CST has %s = %s; ADS has %s = %s\n'], axName, lC, axName, lA);
        continue
    end
    u = u(both);  jC = jC(both);  jA = jA(both);  ang = ang(both);  K = numel(u);
end
only = cuts(ic).only;
if ~isempty(only)
    keep = arrayfun(@(x) any(abs(x - only) < angle_tol), u);
    u = u(keep);  jC = jC(keep);  jA = jA(keep);  ang = ang(keep);  K = numel(u);
    if K == 0, fprintf('none of the requested spokes is available; figure skipped\n'); continue, end
end
fprintf('plotted: %s = %s deg\n', axName, strjoin(compose('%g', u), ', '));

%% ---- common time grid -------------------------------------------------------
if any(jC), r_m = cst(jC(find(jC, 1))).r; else, r_m = cst(1).r; end
tEnds = [];  dts = [];
for k = 1:K
    if jC(k)
        for c = cIdx
            if ~isempty(cst(jC(k)).t{c})
                tEnds(end+1) = cst(jC(k)).t{c}(end); dts(end+1) = median(diff(cst(jC(k)).t{c})); %#ok<SAGROW>
            end
        end
    end
    if jA(k), tEnds(end+1) = ads(jA(k)).t(end); dts(end+1) = median(diff(ads(jA(k)).t)); end %#ok<SAGROW>
end
dt = min(dts);  tend = max(tEnds);                    % every record shown to its own end
tg = (0:dt:tend).';
Vc = nan(numel(tg), K, nC);  Va = Vc;
for k = 1:K
    for c = 1:nC
        cc = cIdx(c);                                  % component in the read data
        if jC(k) && ~isempty(cst(jC(k)).t{cc})
            tk = cst(jC(k)).t{cc};               % zero before the record, NaN (no trace) after it
            Vc(:, k, c) = interp1(tk, cst(jC(k)).e{cc}, tg, 'linear', NaN);  Vc(tg < tk(1), k, c) = 0;
        end
        if jA(k) && ~all(isnan(ads(jA(k)).e(:, cc)))
            tk = ads(jA(k)).t;
            Va(:, k, c) = interp1(tk, ads(jA(k)).e(:, cc), tg, 'linear', NaN);  Va(tg < tk(1), k, c) = 0;
        end
    end
end

% draw only what both sources hold: a trace pair is cut to the overlap of
% the two records; a component one source lacks at a spoke is not drawn
ends = zeros(0, 3);
for k = 1:K
    for c = 1:nC
        hc = ~isnan(Vc(:, k, c));  ha = ~isnan(Va(:, k, c));
        if ~any(hc) || ~any(ha)
            if any(hc) || any(ha)
                miss = "ADS";  if ~any(hc), miss = "CST"; end
                fprintf('E_%s at %s = %g deg: %s has no data; not drawn\n', lower(components(c)), axName, u(k), miss);
            end
            Vc(:, k, c) = NaN;  Va(:, k, c) = NaN;  continue
        end
        both = hc & ha;
        if any(hc & ~ha) || any(ha & ~hc)
            ends(end+1, :) = 1e9*[tg(find(hc, 1, 'last')), tg(find(ha, 1, 'last')), tg(find(both, 1, 'last'))]; %#ok<AGROW>
        end
        Vc(~both, k, c) = NaN;  Va(~both, k, c) = NaN;
    end
end
if ~isempty(ends)   % one line per cut rather than per trace
    ends = uniquetol(ends, 1e-3, 'ByRows', true, 'DataScale', 1);
    for e = ends', fprintf('records: CST to %.2f ns, ADS to %.2f ns; drawn to %.2f ns\n', e); end
end
% a component that is numerically zero in both files (e.g. E_theta on the
% z axis at phi = 90) is still drawn, so the viewer sees the null, but it
% gets no error: relative errors against ~1e-16 V are meaningless
% (same rule as compare_rE_time_at_probe_CST_ADS)
pkAll = max(abs([Vc(:); Va(:)]), [], 'omitnan');
isZero = false(K, nC);
for k = 1:K
    for c = 1:nC
        pk = max(abs([Vc(:, k, c); Va(:, k, c)]), [], 'omitnan');
        if ~isempty(pk) && ~isnan(pk) && pk < 1e-9 * pkAll
            fprintf('E_%s at %s = %g deg is zero in both files (peak %.2g V); drawn, no error computed\n', lower(components(c)), axName, u(k), pk);
            isZero(k, c) = true;
        end
    end
end

%% ---- errors re CST over the overlap of the records, where both have the direction
% tg runs from 0 to the latest record end; each error uses the samples where
% both traces of that direction exist (the overlap of the two records)
E2 = @(x, m) trapz(tg(m), x(m).^2);
fprintf('\n  %5s  comp   peak CST     peak ADS     rel. MSE (re CST comp.)   MSE re CST total field\n', axName);
for k = find(jC & jA)
    have = squeeze(any(~isnan(Vc(:, k, :)), 1) & any(~isnan(Va(:, k, :)), 1)).';
    cOK  = find(have & ~isZero(k, :));                  % components both sources hold (non-zero)
    if isempty(cOK), continue, end
    m    = all(~isnan(Vc(:, k, cOK)) & ~isnan(Va(:, k, cOK)), 3);   % common overlap window
    nTot = sum(arrayfun(@(c) E2(Vc(:, k, c), m), cOK));
    for c = cOK
        nd = E2(Va(:, k, c) - Vc(:, k, c), m);
        nc = E2(Vc(:, k, c), m);
        if nc > 0, eC = sprintf('%9.3g %%', 100*nd/nc); else, eC = '      n/a'; end
        fprintf('  %5g  %-5s  %9.4g V  %9.4g V   %s                 %9.3g %%\n', u(k), components(c), ...
            max(abs(Vc(:, k, c))), max(abs(Va(:, k, c))), eC, 100*nd/nTot);
    end
end
if ~any(jC & jA), fprintf('  (no direction in the cut is in both files)\n'); end

%% ---- radial plots, one panel per component ----------------------------------
% time axis: start just before the field can arrive
if isempty(t_start), ts = (r_m - a_enclosing)/c0 - 0.5e-9; else, ts = t_start*1e-9; end
ip   = tg >= max(ts, 0);
tp   = tg(ip);  Vc = Vc(ip, :, :);  Va = Va(ip, :, :);
span = tp(end) - tp(1);
hole = radial_hole / (1 - radial_hole) * span;        % radius of t(1)
rmax = hole + span;
vmax = max(abs([Vc(:); Va(:)]), [], 'omitnan');
sa   = sort(mod(ang, 360));
if K > 1, dth = min(diff([sa, sa(1) + 360])); else, dth = 90; end
vs   = 0.45 * rmax * sind(min(dth, 90)) / vmax;       % shared volt scale
e10  = 10^floor(log10(vmax/2));  m = vmax/2/e10;      % shared scale-bar value
vref = e10 * ((m < 1.5)*1 + (m >= 1.5 & m < 3.5)*2 + (m >= 3.5 & m < 7.5)*5 + (m >= 7.5)*10);
cols = [0 0 0; 0.85 0 0];                              % CST black, ADS red
if cut == "theta", tabName = sprintf('rE polar %s=%g', char(966), phi_cut);
else,               tabName = sprintf('rE polar %s=%g', char(952), theta_cut); end
fig = figure('Color', 'w', 'NumberTitle', 'off', 'Name', tabName);   % docked: a tab in the figure window
tl  = tiledlayout(fig, 1, nC, 'TileSpacing', 'none', 'Padding', 'tight');
common = {'VoltScale', vs, 'ScaleRef', vref, 'TimeUnit', 'ns', 'TimeFactor', 1e9, ...
    'AngleLabel', angLabel, 'FontSize', font_size, 'RadialOffset', hole};
for c = 1:nC
    ax = nexttile(tl);
    [~, hc] = radialTracePlot(ang, tp, Vc(:, :, c), 'Parent', ax, common{:}, ...
        'Colors', repmat(cols(1,:), K, 1), 'AngleGrid', angle_grid, ...
        'RingLabelAngle', ring_label_angle);
    [~, ha] = radialTracePlot(ang, tp, Va(:, :, c), 'Parent', ax, common{:}, ...
        'Colors', repmat(cols(2,:), K, 1), 'Grid', false, 'Fill', false, 'LineStyle', ':');
    % swept-angle indicator: arc + arrowhead between the first two angle
    % labels, pointing in the direction of increasing angle
    if cut == "theta", a1 = 83; a2 = 67; swName = '$\theta$';
    else,               a1 = 7;  a2 = 23; swName = '$\phi$'; end
    ra = 1.08 * rmax;  aa = linspace(a1, a2, 40);  hold(ax, 'on');
    plot(ax, ra*cosd(aa), ra*sind(aa), 'k-', 'LineWidth', 1.2);
    tng = sign(a2 - a1) * [-sind(a2) cosd(a2)];  nrm = [cosd(a2) sind(a2)];
    p = ra * nrm;  hl = 0.045*rmax;  hw = 0.02*rmax;
    patch(ax, [p(1)+tng(1)*hl, p(1)+nrm(1)*hw, p(1)-nrm(1)*hw], ...
              [p(2)+tng(2)*hl, p(2)+nrm(2)*hw, p(2)-nrm(2)*hw], 'k', 'EdgeColor', 'none');
    am = (a1 + a2)/2;
    text(ax, 1.17*rmax*cosd(am), 1.17*rmax*sind(am), swName, 'Interpreter', 'latex', ...
        'FontSize', font_size + 4, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    % trim the empty margin around the circle (labels reach ~1.15 rmax,
    % the scale bar and its label ~1.3 rmax on the right)
    xlim(ax, [-1.18 1.32] * rmax);  ylim(ax, [-1.16 1.23] * rmax);
    title(ax, {sprintf('$rE_\\%s(t)$', lower(components(c))), panelCut}, ...
        'Interpreter', 'latex', 'FontSize', font_size + 2);
    if c == 1
        lgd = legend(ax, [hc.traces(find(jC, 1)) ha.traces(find(jA, 1))], {'CST', 'ADS'}, ...
            'Location', 'southwest', 'FontSize', font_size + 2, 'Interpreter', 'latex');
        lgd = nudgeLegend(lgd, ax);                    % keep it bottom-left, clear of the angle labels
    end
end
title(tl, {sprintf('%s, $r = %g$ m', model_name, r_m), cutName, sweepName}, ...
    'FontSize', font_size + 2, 'Interpreter', 'latex');
end                                                    % ===== end of cut loop =====

%% ---- local functions --------------------------------------------------------
function lgd = nudgeLegend(lgd, ax)
% Keep the legend in the lower-left corner but move it down and left (not
% past the figure edge) until it no longer covers any text label of ax.
drawnow;
fig = ancestor(ax, 'figure');  set([lgd ax], 'Units', 'pixels');
pa = ax.Position;  xl = ax.XLim;  yl = ax.YLim;
toPix = @(x, y) [pa(1) + (x - xl(1))/diff(xl)*pa(3), pa(2) + (y - yl(1))/diff(yl)*pa(4)];
boxes = zeros(0, 4);
for tx = findobj(ax, 'type', 'text')'
    tx.Units = 'data';  ex = tx.Extent;
    p1 = toPix(ex(1), ex(2));  p2 = toPix(ex(1) + ex(3), ex(2) + ex(4));
    boxes(end+1, :) = [p1, p2 - p1]; %#ok<AGROW>
end
pl = lgd.Position;  gap = 10;                          % px clearance from the labels
for it = 1:20
    hit = boxes(rectint(pl, boxes) > 0, :);
    if isempty(hit), break, end
    dx = max(pl(1) + pl(3) - hit(:, 1)) + gap;           % shift left to clear the labels
    dy = max(pl(2) + pl(4) - hit(:, 2)) + gap;           % or shift down
    if pl(1) - dx >= 2, pl(1) = pl(1) - dx;
    elseif pl(2) - dy >= 2, pl(2) = pl(2) - dy;
    else, pl(1:2) = [2 2]; break, end
end
lgd.Position = pl;  set([lgd ax], 'Units', 'normalized');
if ~isempty(boxes(rectint(pl, boxes) > 0, :)), lgd.Location = 'southwestoutside'; end   % last resort
end

function u = cutCoord(th, ph, cut, phi_cut, theta_cut, tol)
% Position of direction (th, ph) [deg] in the cut: signed theta for a theta
% cut (negative in the half-plane phi_cut + 180), phi in 0..360 for a phi
% cut; NaN when the direction is not in the cut.
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
u = NaN;
if cut == "theta"
    if th < tol || th > 180 - tol
        % on the z axis the theta/phi unit vectors still rotate with phi, so
        % only the sample at phi = phi_cut is used (frame of the right half)
        if angd(ph, phi_cut) < tol, u = th; end
    elseif angd(ph, phi_cut) < tol
        u = th;
    elseif angd(ph, phi_cut + 180) < tol
        u = -th;
    end
elseif cut == "phi"
    if abs(th - theta_cut) < tol
        u = mod(ph, 360);
    end
else
    error('cut must be "theta" or "phi".');
end
end
