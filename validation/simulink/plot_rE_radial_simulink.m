function [fig, res, leg] = plot_rE_radial_simulink(sim_sig, cst, ads, varargin)
%PLOT_RE_RADIAL_SIMULINK  Radial (spoke) plots of the radiated far-field time signals
%   r*E_theta(t), r*E_phi(t): Simulink link results against CST and ADS, drawn with
%   radialTracePlot (one spoke per observation direction of a cut, time along the spoke,
%   voltage across it; one panel per component, shared volt scale).  Counterpart of
%   validation/cst_ads/time_signals/plot_rE_radial_CST_ADS.m with the Simulink traces added.
%   Simulink: light blue, wide (drawn first, underneath).  CST: black, shaded.  ADS: red, dotted.
%
%   [fig, res, leg] = plot_rE_radial_simulink(sim_sig, cst, ads, Name, Value, ...)
%     sim_sig : struct array, one element per simulated direction, fields
%               theta, phi (deg), t (s), e (N x 2: [rE_theta rE_phi] in V), or []
%     cst     : struct array from readCstFarfieldProbes (theta, phi, r, t{c}, e{c}), or []
%     ads     : struct array from readAdsFarfieldSweeps (theta, phi, t, e), or []
%
%   'Sources' ["Simulink" "CST" "ADS"]: which of the three are drawn (data of the others is ignored).
%   Errors are relative to the first of CST, ADS, Simulink that is drawn.
%
%   Cut (which directions form the spokes: every direction of the cut in the sources drawn):
%     'Cut'    "theta" (default): plane through the z axis at phi = Angle; the half-plane
%              phi = Angle + 180 is drawn on the left as negative theta (z up, plot angle
%              90 - theta).   "phi": cone theta = Angle (phi 0..360, x axis to the right,
%              counter-clockwise).
%     'Angle'  phi (theta cut) or theta (phi cut) in degrees (default 90)
%   Other options:
%     'Components' ["Theta" "Phi"]   panels      'AllComponents' ["Theta" "Phi"]  column order of e
%     'Spokes' []                    only these spokes: signed theta for a theta cut (negative = half-plane phi + 180),
%                                    phi for a phi cut
%     'r' 1 (m), 'a' 0 (m)           time axis starts at (r-a)/c0 - 0.5 ns unless 'TStart' (ns)
%     'RingLabelAngle' 255, 'AngleGrid' 30, 'FontSize' 12, 'Hole' 0.12 (fraction of the radius
%     left empty at the centre), 'AngleTol' 0.1, 'Title' (figure/sheet title),
%     'SimLabel' 'Simulink', 'SaveFile' (PNG path; a .fig is written next to it)
%     'FigureSize' [] ([w h] cm: undocked figure of that size; [] = docked tab),
%     'Subtitle' "top" (cut label under the component title) | "grid" (on the polar grid, at
%     'SubtitleAngle' 135 deg and 'SubtitleRadius' 0.95 of the outer radius) | "none",
%     'LineWidths' [] ([Simulink CST ADS] in pt; default [3.5 1 1.4]),
%     'TimeShift' 0 (s, subtracted from every time axis: de-embedded propagation delay, e.g. r/c0),
%     'TimeArg' 't' (label argument: rE_phi(TimeArg), e.g. 't + r/c_0'),
%     'LegendLocation' 'southwest', 'ShowLegend' true, 'CompLabel' "top" | "bottom" | "left" | "none" (component label),
%     'Parent' [] (an existing axes, e.g. a tile of a multi-panel figure, to draw into; one component only)
%     'VMax' [] (V: the volt scale is set for this peak, e.g. the same for several plots; [] = largest |V| drawn)
%     'SimAlpha' 1 (opacity of the Simulink traces), 'SimOnTop' false (Simulink in front of CST / ADS),
%     'LineStyles' [] ([Simulink CST ADS], default ["-" "-" ":"]),
%     'Mirror' false (each trace also drawn negated: a two-sided +/- envelope when the records are envelopes)
%     'OnTop' [] (sources brought to the front in this order, last frontmost: "Simulink" | "CST" | "ADS")
%     'Envelope' [] (struct t (s, the time base after 'TimeShift'), env (peak 1), col, alpha, name): a shaded
%                    +/- envelope on every spoke, scaled to the reference trace (least squares on |hilbert|)
%
%   leg: struct with fields h (line handles) and names, for a legend shared by several panels.
%
%   res: struct array (spoke, component) with the relative MSE of each drawn source against
%   the reference over the overlap of the records.

p = inputParser;
addParameter(p, 'Cut', "theta");
addParameter(p, 'Angle', 90);
addParameter(p, 'Components', ["Theta" "Phi"]);
addParameter(p, 'AllComponents', ["Theta" "Phi"]);
addParameter(p, 'AngleTol', 0.1);
addParameter(p, 'RingLabelAngle', 255);
addParameter(p, 'AngleGrid', 30);
addParameter(p, 'FontSize', 12);
addParameter(p, 'Hole', 0.12);
addParameter(p, 'TStart', []);
addParameter(p, 'r', 1);
addParameter(p, 'a', 0);
addParameter(p, 'Title', '');
addParameter(p, 'SimLabel', 'Simulink');
addParameter(p, 'SaveFile', '');
addParameter(p, 'Sources', ["Simulink" "CST" "ADS"]);
addParameter(p, 'FigureSize', []);         % [width height] in cm: undocked figure of that size (e.g. for publication); [] = docked tab
addParameter(p, 'Subtitle', "top");       % cut label (e.g. phi = 0): "top" under the component title, "grid" on the polar grid, "none"
addParameter(p, 'SubtitleAngle', 135);     % deg, position of the label on the grid ("grid"), counter-clockwise from the x axis
addParameter(p, 'SubtitleRadius', 0.95);   % fraction of the outer radius ("grid")
addParameter(p, 'LegendLocation', 'southwest');
addParameter(p, 'ShowLegend', true);
addParameter(p, 'Verbose', true);          % false: no console table (the errors are still returned in res)
addParameter(p, 'TimeShift', 0);           % s, taken off every record (de-embedded propagation delay, e.g. r/c0)
addParameter(p, 'TimeArg', 't');           % argument in the component label, e.g. 't + r/c_0'
addParameter(p, 'NumRings', 5);            % approximate number of time rings (radialTracePlot)
addParameter(p, 'LineWidths', []);         % [Simulink CST ADS] line widths (pt); [] = [3.5 1 1.4]
addParameter(p, 'CompLabel', "top");      % component label rE_theta(t) / rE_phi(t): "top" (title), "bottom" (under the grid), "left" (vertical, left of the grid), "none"
addParameter(p, 'Parent', []);             % existing axes to draw into (one component), e.g. a tile of a multi-panel figure
addParameter(p, 'VMax', []);              % V, peak used for the volt scale ([] = largest |V| drawn)
addParameter(p, 'Envelope', []);          % struct t, env, col, alpha, name: shaded envelope on every spoke
addParameter(p, 'LineStyles', strings(0)); % [Simulink CST ADS] line styles; empty = ["-" "-" ":"]
addParameter(p, 'Mirror', false);         % true: each trace also drawn negated (two-sided envelope)
addParameter(p, 'SimAlpha', 1);           % opacity of the Simulink traces (< 1: drawn as transparent patch edges)
addParameter(p, 'SimOnTop', false);       % true: Simulink traces brought in front of CST and ADS
addParameter(p, 'OnTop', strings(0));     % sources brought to the front in this order (last = frontmost), e.g. ["Simulink" "CST"]
addParameter(p, 'Spokes', []);    % only these spokes: signed theta (deg; negative = half-plane phi + 180) for a theta cut, phi for a phi cut; [] = all
parse(p, varargin{:});
o = p.Results;
cut = string(o.Cut);  components = string(o.Components);  allComp = string(o.AllComponents);
c0 = 299792458;  tol = o.AngleTol;
if exist('radialTracePlot', 'file') ~= 2
    tsdir = fullfile(fileparts(mfilename('fullpath')), '..', 'cst_ads', 'time_signals', 'utils');
    if isfolder(tsdir), addpath(tsdir); else, error('radialTracePlot.m not found on the path.'); end
end
srcSel = string(o.Sources);
if isempty(sim_sig) || ~any(srcSel == "Simulink"), sim_sig = struct('theta', {}, 'phi', {}, 't', {}, 'e', {}); end
if isempty(cst) || ~any(srcSel == "CST"), cst = struct('theta', {}, 'phi', {}, 'r', {}, 't', {}, 'e', {}); end
if isempty(ads) || ~any(srcSel == "ADS"), ads = struct('theta', {}, 'phi', {}, 't', {}, 'e', {}); end

% de-embedded propagation delay: every record shifted by -TimeShift
if o.TimeShift ~= 0
    for j = 1:numel(sim_sig), sim_sig(j).t = sim_sig(j).t - o.TimeShift; end
    for j = 1:numel(cst), for q_ = 1:numel(cst(j).t), cst(j).t{q_} = cst(j).t{q_} - o.TimeShift; end, end
    for j = 1:numel(ads), ads(j).t = ads(j).t - o.TimeShift; end
end

%% ---- spokes of the cut --------------------------------------------------------
if cut == "theta", phi_cut = o.Angle; theta_cut = NaN; else, theta_cut = o.Angle; phi_cut = NaN; end
if cut == "theta"     % z-axis points given on the far half-plane (e.g. theta = -180 / phi + 180): one spoke, not two
    sim_sig = axisToCutHalf(sim_sig, phi_cut, tol, false);
    cst     = axisToCutHalf(cst, phi_cut, tol, true);
    ads     = axisToCutHalf(ads, phi_cut, tol, false);
end
cc = @(s) cutCoord(s.theta, s.phi, cut, phi_cut, theta_cut, tol);
uS = arrayfun(cc, sim_sig);  uC = arrayfun(cc, cst);  uA = arrayfun(cc, ads);
uAll = [uS(:); uC(:); uA(:)];
u = sort(uniquetol(uAll(~isnan(uAll)), tol, 'DataScale', 1)).';
if ~isempty(o.Spokes)
    if cut == "theta", keep = arrayfun(@(x) any(abs(x - o.Spokes) < tol), u);
    else,              keep = arrayfun(@(x) any(abs(mod(x - o.Spokes + 180, 360) - 180) < tol), u); end
    u = u(keep);
end
if isempty(u), error('No direction of the sources drawn lies in this cut.'); end
K = numel(u);
jS = zeros(1, K); jC = zeros(1, K); jA = zeros(1, K);
for k = 1:K
    m = find(abs(uS - u(k)) < tol, 1);  if ~isempty(m), jS(k) = m; end
    m = find(abs(uC - u(k)) < tol, 1);  if ~isempty(m), jC(k) = m; end
    m = find(abs(uA - u(k)) < tol, 1);  if ~isempty(m), jA(k) = m; end
end
if cut == "theta"
    axName = 'theta';  ang = 90 - u;
    panelCut = sprintf('$\\phi = %g^\\circ$ / $%g^\\circ$', phi_cut, mod(phi_cut + 180, 360));
    if all(u > -tol), panelCut = sprintf('$\\phi = %g^\\circ$', phi_cut); end   % only the phi_cut half-plane drawn
    angLabel = @(a) sprintf('$%g^\\circ$', mod(90 - a + 180, 360) - 180);
    swName = '$\theta$';  a1 = 83;  a2 = 67;
else
    axName = 'phi';  ang = u;
    panelCut = sprintf('$\\theta = %g^\\circ$', theta_cut);
    angLabel = @(a) sprintf('$%g^\\circ$', mod(a, 360));
    swName = '$\phi$';  a1 = 7;  a2 = 23;
end
if o.Verbose
    fprintf('radial plot, %s: %s = %s deg (%s: %d, CST: %d, ADS: %d of %d)\n', panelCut, axName, ...
        strjoin(compose('%g', u), ', '), o.SimLabel, nnz(jS), nnz(jC), nnz(jA), K);
end

%% ---- common time grid, all sources interpolated onto it ------------------------------
dts = [];  tEnds = [];
for k = 1:K
    if jS(k), ts_ = sim_sig(jS(k)).t(:);  dts(end+1) = median(diff(ts_));  tEnds(end+1) = ts_(end); end %#ok<AGROW>
    if jC(k)
        for c = 1:numel(cst(jC(k)).t)
            if ~isempty(cst(jC(k)).t{c}), dts(end+1) = median(diff(cst(jC(k)).t{c})); tEnds(end+1) = cst(jC(k)).t{c}(end); end %#ok<AGROW>
        end
    end
    if jA(k), dts(end+1) = median(diff(ads(jA(k)).t)); tEnds(end+1) = ads(jA(k)).t(end); end %#ok<AGROW>
end
dt = min(dts);  tg = (0:dt:min(tEnds)).';     % plots end at the minimum of the last times of all records
cIdx = arrayfun(@(s) find(allComp == s), components);
nC = numel(components);
Vs = nan(numel(tg), K, nC);  Vc = Vs;  Va = Vs;
for k = 1:K
    for c = 1:nC
        q = cIdx(c);
        if jS(k), Vs(:, k, c) = onGrid(sim_sig(jS(k)).t, sim_sig(jS(k)).e(:, q), tg); end
        if jC(k) && numel(cst(jC(k)).e) >= q && ~isempty(cst(jC(k)).e{q})
            Vc(:, k, c) = onGrid(cst(jC(k)).t{q}, cst(jC(k)).e{q}, tg);
        end
        if jA(k) && ~all(isnan(ads(jA(k)).e(:, q)))
            Va(:, k, c) = onGrid(ads(jA(k)).t, ads(jA(k)).e(:, q), tg);
        end
    end
end

%% ---- errors against the reference (first of CST, ADS, Simulink drawn) over the overlap of the records ----
if ~all(isnan(Vc(:))),     refName = "CST";               Vr = Vc;
elseif ~all(isnan(Va(:))), refName = "ADS";               Vr = Va;
else,                      refName = string(o.SimLabel);  Vr = Vs; end
res = struct('axis', {}, 'u', {}, 'comp', {}, 'ref', {}, 'peak_cst', {}, 'peak_ads', {}, 'peak_sim', {}, ...
    'relmse_cst', {}, 'relmse_ads', {}, 'relmse_sim', {}, 'magerr_cst', {}, 'magerr_ads', {}, 'magerr_sim', {}, ...
    'eref', {}, 'ediff', {});      % eref: integral ref^2 dt; ediff [CST ADS Sim]: integral (x - ref)^2 dt (for a cut total)
if o.Verbose
    fprintf('  %5s  comp   peak CST     peak ADS     peak %s   relMSE CST   relMSE ADS   relMSE %s (compared to %s)\n', ...
        axName, o.SimLabel, o.SimLabel, refName);
end
pkAll = max(abs([Vs(:); Vc(:); Va(:)]), [], 'omitnan');
for k = 1:K
    for c = 1:nC
        pk = max(abs([Vs(:, k, c); Vc(:, k, c); Va(:, k, c)]), [], 'omitnan');
        if pk < 1e-4*pkAll || ~any(~isnan(Vr(:, k, c))), continue, end    % null component (80 dB down: numerical noise) / no reference
        E2 = @(x, m) trapz(tg(m), x(m).^2);
        vr = Vr(:, k, c);  eR = nan(1, 3);  eD = nan(1, 3);  V3 = {Vc(:, k, c), Va(:, k, c), Vs(:, k, c)};
        eDiff = nan(1, 3);  eRef = E2(vr, ~isnan(vr));
        for j = 1:3
            m_ = ~isnan(V3{j}) & ~isnan(vr);
            if any(m_) && ~isequaln(V3{j}, vr)
                eDiff(j) = E2(V3{j} - vr, m_);
                eR(j) = eDiff(j) / E2(vr, m_);
                eD(j) = bandMagErr(vr(m_), V3{j}(m_), 20);    % mean |dB| error, -20 dB band of the reference
            end
        end
        res(end+1) = struct('axis', axName, 'u', u(k), 'comp', components(c), 'ref', refName, ...
            'peak_cst', max(abs(Vc(:, k, c)), [], 'omitnan'), 'peak_ads', max(abs(Va(:, k, c)), [], 'omitnan'), ...
            'peak_sim', max(abs(Vs(:, k, c)), [], 'omitnan'), 'relmse_cst', eR(1), 'relmse_ads', eR(2), 'relmse_sim', eR(3), ...
            'magerr_cst', eD(1), 'magerr_ads', eD(2), 'magerr_sim', eD(3), 'eref', eRef, 'ediff', eDiff); %#ok<AGROW>
        if o.Verbose
            fprintf('  %5g  %-5s  %9.4g V  %9.4g V  %9.4g V   %7.2f %%     %7.2f %%     %7.2f %%\n', u(k), components(c), ...
                res(end).peak_cst, res(end).peak_ads, res(end).peak_sim, 100*eR);
        end
    end
end

%% ---- radial plots, one panel per component ---------------------------------------------------------
if isempty(o.TStart), tst = (o.r - o.a)/c0 - o.TimeShift - 0.5e-9; else, tst = o.TStart*1e-9; end
ip = tg >= max(tst, 0);
tp = tg(ip);  Vs = Vs(ip, :, :);  Vc = Vc(ip, :, :);  Va = Va(ip, :, :);  Vr = Vr(ip, :, :);
span = tp(end) - tp(1);
hole = o.Hole/(1 - o.Hole)*span;
rmax = hole + span;
vmax = max(abs([Vs(:); Vc(:); Va(:)]), [], 'omitnan');
if ~isempty(o.VMax), vmax = o.VMax; end
sa = uniquetol(mod(ang, 360), tol/360);     % distinct spoke angles (coincident spokes would give a zero volt scale)
if numel(sa) > 1, dth = min(diff([sa(:).', sa(1) + 360])); else, dth = 90; end
vs = 0.45*rmax*sind(min(dth, 90))/vmax;
e10 = 10^floor(log10(vmax/2));  m = vmax/2/e10;
vref = e10*((m < 1.5)*1 + (m >= 1.5 & m < 3.5)*2 + (m >= 3.5 & m < 7.5)*5 + (m >= 7.5)*10);
% draw order: Simulink (wide, underneath), then CST (thin, shaded), then ADS (dotted) on top, so all stay visible
srcs = struct('name', {o.SimLabel, 'CST', 'ADS'}, 'V', {Vs, Vc, Va}, 'col', {[0.55 0.75 1], [0 0 0], [0.85 0 0]}, ...
    'ls', {'-', '-', ':'}, 'fill', {false, true, false}, 'lw', {3.5, 1, 1.4});
if ~isempty(o.LineWidths), for s_ = 1:3, srcs(s_).lw = o.LineWidths(s_); end, end
if ~isempty(o.LineStyles), for s_ = 1:3, srcs(s_).ls = char(string(o.LineStyles(s_))); end, end
srcTag = strjoin(srcSel, '/');
if cut == "theta", tabName = sprintf('rE polar %s %s=%g', srcTag, char(966), phi_cut);
else,              tabName = sprintf('rE polar %s %s=%g', srcTag, char(952), theta_cut); end
tl = [];
if ~isempty(o.Parent)
    if nC > 1, error('''Parent'': one component only (got %d).', nC); end
    fig = ancestor(o.Parent, 'figure');
else
    if isempty(o.FigureSize)
        fig = figure('Color', 'w', 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', tabName);
    else
        fig = figure('Color', 'w', 'WindowStyle', 'normal', 'NumberTitle', 'off', 'Name', tabName, ...
            'Units', 'centimeters', 'Position', [2 2 o.FigureSize(:).'], 'PaperPositionMode', 'auto');
    end
    tl = tiledlayout(fig, 1, nC, 'TileSpacing', 'none', 'Padding', 'tight');
end
leg = struct('h', gobjects(0), 'names', strings(0));
common = {'VoltScale', vs, 'ScaleRef', vref, 'TimeUnit', 'ns', 'TimeFactor', 1e9, 'AngleLabel', angLabel, ...
    'FontSize', o.FontSize, 'RadialOffset', hole, 'NumRings', o.NumRings};
for c = 1:nC
    if isempty(o.Parent), ax = nexttile(tl); else, ax = o.Parent; end
    first = true;  lh = gobjects(0);  ln = strings(0);  hSrc = cell(1, numel(srcs));  hSim = gobjects(0);
    for s = 1:numel(srcs)
        V = srcs(s).V(:, :, c);
        k1 = find(any(~isnan(V), 1), 1);
        if isempty(k1), continue, end
        args = {'Colors', repmat(srcs(s).col, K, 1), 'LineStyle', srcs(s).ls, 'Fill', srcs(s).fill};
        if first
            args = [args, {'AngleGrid', o.AngleGrid, 'RingLabelAngle', o.RingLabelAngle}]; %#ok<AGROW>
        else
            args = [args, {'Grid', false}]; %#ok<AGROW>
        end
        hT = gobjects(0);  hT1 = [];
        for sg = [1, -ones(1, o.Mirror)]     % 'Mirror': the trace and its negative (two-sided envelope)
        [~, h] = radialTracePlot(ang, tp, sg*V, 'Parent', ax, common{:}, args{:});
        set(h.traces, 'LineWidth', srcs(s).lw);
        if isempty(hT1), hT1 = h.traces; end
        hT = [hT(:); h.traces(:)];
        args = [args, {'Grid', false}]; %#ok<AGROW>
        end
        if s == 1 && o.SimAlpha < 1          % Simulink: transparent (patch edges, so the alpha survives vector export)
            for i_ = 1:numel(hT)
                if ~isgraphics(hT(i_)), continue, end
                x_ = hT(i_).XData(:);  y_ = hT(i_).YData(:);
                pT = patch(ax, [x_; NaN], [y_; NaN], 'k', 'FaceColor', 'none', 'EdgeColor', srcs(s).col, ...
                    'EdgeAlpha', o.SimAlpha, 'LineWidth', srcs(s).lw, 'LineStyle', srcs(s).ls);
                delete(hT(i_));  hT(i_) = pT;
            end
        end
        if s == 1, hSim = hT(isgraphics(hT)); end
        hSrc{s} = hT(isgraphics(hT)); %#ok<AGROW>
        lh(end+1) = hT1(k1);  ln(end+1) = string(srcs(s).name); %#ok<AGROW>
        first = false;
    end
    if o.SimOnTop && ~isempty(hSim), uistack(hSim, 'top'); end
    for nm_ = string(o.OnTop)                % e.g. ["Simulink" "CST"]: CST frontmost
        s_ = find([string(srcs(1).name) "CST" "ADS"] == nm_ | (["Simulink" "CST" "ADS"] == nm_), 1);
        if ~isempty(s_) && ~isempty(hSrc{s_}), uistack(hSrc{s_}, 'top'); end
    end
    if ~isempty(o.Envelope)              % shaded envelope on every spoke, scaled to the reference trace
        En = o.Envelope;  en = interp1(En.t(:), En.env(:), tp, 'linear', 0);  rr = hole + (tp - tp(1));  hE = gobjects(0);
        for k = 1:K
            vr = Vr(:, k, c);  m = ~isnan(vr);
            if ~any(m) || max(abs(vr(m))) < 1e-4*pkAll, continue, end
            ar = zeros(size(vr));  ar(m) = abs(hilbert(vr(m)));
            g = (ar(m)'*en(m)) / max(en(m)'*en(m), realmin);
            d = vs*g*en;  d(~m) = 0;  cA = cosd(ang(k));  sA = sind(ang(k));
            hE = patch(ax, [rr*cA - d*sA; flipud(rr*cA + d*sA)], [rr*sA + d*cA; flipud(rr*sA - d*cA)], En.col, ...
                'FaceAlpha', En.alpha, 'EdgeColor', 'none');
        end
        if ~isempty(hE), lh(end+1) = hE;  ln(end+1) = string(En.name); end %#ok<AGROW>
    end
    hold(ax, 'on');
    ra = 1.08*rmax;  aa = linspace(a1, a2, 40);
    plot(ax, ra*cosd(aa), ra*sind(aa), 'k-', 'LineWidth', 1.2);
    tng = sign(a2 - a1)*[-sind(a2) cosd(a2)];  nrm = [cosd(a2) sind(a2)];
    pt = ra*nrm;  hl = 0.045*rmax;  hw = 0.02*rmax;
    patch(ax, [pt(1) + tng(1)*hl, pt(1) + nrm(1)*hw, pt(1) - nrm(1)*hw], ...
              [pt(2) + tng(2)*hl, pt(2) + nrm(2)*hw, pt(2) - nrm(2)*hw], 'k', 'EdgeColor', 'none');
    am = (a1 + a2)/2;
    text(ax, 1.17*rmax*cosd(am), 1.17*rmax*sind(am), swName, 'Interpreter', 'latex', ...
        'FontSize', o.FontSize + 4, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim(ax, [-1.18 1.32]*rmax);  ylim(ax, [-1.16 1.23]*rmax);
    compTitle = sprintf('$rE_\\%s(%s)$', lower(components(c)), o.TimeArg);
    ttl = {};
    if string(o.CompLabel) == "top", ttl{end+1} = compTitle; end
    if string(o.Subtitle) == "top",  ttl{end+1} = panelCut; end
    if ~isempty(ttl), title(ax, ttl, 'Interpreter', 'latex', 'FontSize', o.FontSize + 2); end
    if string(o.CompLabel) == "left"         % vertical, left of the grid (like a y label)
        % left of the leftmost angle label (e.g. -90 deg / 180 deg) so the two do not overlap
        drawnow;  ta = findall(ax, 'Type', 'text');  xlft = -1.18*rmax;
        if ~isempty(ta), ex = cell2mat(get(ta, 'Extent'));  xlft = min([xlft; ex(:, 1)]); end
        text(ax, xlft - 0.1*rmax, 0, compTitle, 'Interpreter', 'latex', 'FontSize', o.FontSize + 4, 'Rotation', 90, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    elseif string(o.CompLabel) == "bottom"   % under the grid, below the bottom angle label
        text(ax, 0, -1.28*rmax, compTitle, 'Interpreter', 'latex', 'FontSize', o.FontSize + 4, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
    if string(o.Subtitle) == "grid"           % cut label on the polar grid
        rs = o.SubtitleRadius*rmax;
        text(ax, rs*cosd(o.SubtitleAngle), rs*sind(o.SubtitleAngle), panelCut, 'Interpreter', 'latex', ...
            'FontSize', o.FontSize + 2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'BackgroundColor', 'w', 'Margin', 2);
    end
    % widen the limits so that no label (e.g. the volt scale with large fonts) is cut off
    drawnow;  tx = findall(ax, 'Type', 'text');
    if ~isempty(tx)
        ex = cell2mat(get(tx, 'Extent'));  xl = xlim(ax);  yl = ylim(ax);
        xlim(ax, [min(xl(1), min(ex(:, 1))), max(xl(2), max(ex(:, 1) + ex(:, 3)))]);
        ylim(ax, [min(yl(1), min(ex(:, 2))), max(yl(2), max(ex(:, 2) + ex(:, 4)))]);
    end
    if c == 1 && ~isempty(lh)
        leg = struct('h', lh, 'names', ln);
        if o.ShowLegend
            legend(ax, lh, cellstr(ln), 'Location', o.LegendLocation, 'FontSize', o.FontSize + 2, 'Interpreter', 'latex');
        end
    end
end
if ~isempty(o.Title) && ~isempty(tl), title(tl, o.Title, 'Interpreter', 'latex', 'FontSize', o.FontSize + 3); end
if ~isempty(o.SaveFile)
    sd = fileparts(o.SaveFile);  if ~isempty(sd) && ~isfolder(sd), mkdir(sd); end
    ax_ = findall(fig, 'Type', 'axes');  for a_ = ax_.', if ~isempty(a_.Toolbar), a_.Toolbar.Visible = 'off'; end, end   % no axes toolbar in the PNG
    exportgraphics(fig, o.SaveFile, 'Resolution', 120);
    [d_, n_] = fileparts(o.SaveFile);  savefig(fig, fullfile(d_, [n_ '.fig']));
end
end

function e = bandMagErr(xr, x, lev)
% mean |20 log10(|X|/|Xr|)| over the band where the reference spectrum is within lev dB of its peak
% (uniformly sampled records of equal length; zero-padded FFT)
Nf = 2^nextpow2(8*numel(xr));
Xr = abs(fft(xr, Nf));  Xr = Xr(1:Nf/2 + 1);  X = abs(fft(x, Nf));  X = X(1:Nf/2 + 1);
b = Xr >= max(Xr)*10^(-lev/20);
e = mean(abs(20*log10(max(X(b), realmin) ./ Xr(b))));
end

function S = axisToCutHalf(S, phi_cut, tol, isCell)
% On the z axis (theta = 0 or 180) the point given at phi_cut + 180 is the same as at phi_cut, with the
% theta and phi unit vectors reversed: move it to phi_cut and negate both components.
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
for j = 1:numel(S)
    if (S(j).theta < tol || abs(S(j).theta - 180) < tol) && angd(S(j).phi, phi_cut + 180) < tol
        S(j).phi = phi_cut;
        if isCell, S(j).e = cellfun(@(x) -x, S(j).e, 'UniformOutput', false); else, S(j).e = -S(j).e; end
    end
end
end

function u = cutCoord(th, ph, cut, phi_cut, theta_cut, tol)
% coordinate of a direction along the cut (signed theta, or phi); NaN if it is not in the cut
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
if cut == "theta"
    if angd(ph, phi_cut) < tol, u = th;
    elseif angd(ph, phi_cut + 180) < tol, u = -th;
    else, u = NaN; end
else
    if abs(th - theta_cut) < tol, u = mod(ph, 360); else, u = NaN; end
end
end

function v = onGrid(t, x, tg)
% x(t) on the grid tg: zero before the record starts, NaN (no trace) after it ends
t = t(:);  x = x(:);
[t, iu] = unique(t, 'last');  x = x(iu);
v = interp1(t, x, tg, 'linear', NaN);
v(tg < t(1)) = 0;
end
