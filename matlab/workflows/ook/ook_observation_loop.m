function res = ook_observation_loop(mdl, dirs, r, cst, t0, label, fastRestart, makePlots, stopTime)
%OOK_OBSERVATION_LOOP  Simulate an already built / loaded Tx-only OOK model at a list of
%   observation points and compare each with the CST far-field probe in that direction.
%
%   res = ook_observation_loop(mdl, dirs, r, cst, t0, label)
%     mdl   : model name (loaded, base-workspace variables set: Vtx_wave, Vs_ctrl, link_cfg, ...)
%     dirs  : [theta phi] in deg, one row per observation point
%     r     : observation distance (m)
%     cst   : struct array from readCstFarfieldProbes (theta, phi, t{c}, e{c}), or [] (no comparison)
%     t0    : start time of the CST excitation grid (s): CST time = Simulink time + t0
%     label : name of the Simulink curve in plots / console (e.g. 'RF Blockset', 'Simscape')
%     fastRestart : true (default) = Simulink Fast Restart: the model is compiled once and the
%             direction gains are updated between runs (r, hence the propagation delay, is the same
%             for every point); false = full compile for every point. Restored to off afterwards.
%     makePlots : true (default) = one figure per point (and the CST comparison); false = simulate only
%     stopTime  : Simulink stop time (s, Simulink time base); [] (default) = the end of the CST records
%             when cst is given (latest probe sample - t0), else prop_t + t_exc
%
%   The direction only changes gains and the propagation delay (link_update_gains), so the
%   model is not rebuilt between points.  One figure (docked tab) per point: rE_theta and
%   rE_phi, Simulink light blue (wide), CST black.  A summary table is printed.
%
%   res : struct array, one element per point: theta, phi, r, t (s, CST time base),
%         e (N x 2: [rE_theta rE_phi] in V), relmse (1 x 2, NaN without a CST probe),
%         peak_sim, peak_cst (1 x 2, V), out (Simulink output, first point only)

if nargin < 6 || isempty(label), label = 'Simulink'; end
if nargin < 7 || isempty(fastRestart), fastRestart = true; end
if nargin < 8 || isempty(makePlots), makePlots = true; end
if nargin < 9, stopTime = []; end
if isempty(stopTime) && ~isempty(cst)            % run to the same final time as the CST records
    stopTime = max(arrayfun(@(c) max(cellfun(@(x) max([x(:); -Inf]), c.t)), cst)) - t0;
end
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);
tol  = 0.1;
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
try, set(groot, 'defaultFigureWindowStyle', 'docked'); catch, end
cn = {'\theta', '\phi'};  nm = {'Vtheta', 'Vphi'};
colSim = [0.55 0.75 1];

nD = size(dirs, 1);
res = struct('theta', {}, 'phi', {}, 'r', {}, 't', {}, 'e', {}, 'relmse', {}, 'err_energy', {}, 'ref_energy', {}, ...
    'peak_sim', {}, 'peak_cst', {}, 'sim_time', {}, 'out', {});
tAll = tic;
% stop time = prop_t + t_exc: prop_t depends on r only, so it is fixed for the whole loop
assignin('base', 'r', r);  assignin('base', 'theta', dirs(1, 1)*pi/180);  assignin('base', 'phi', dirs(1, 2)*pi/180);
evalin('base', 'link_update_gains');
set_param(mdl, 'FastRestart', 'off');
if isempty(stopTime)
    set_param(mdl, 'StopTime', 'prop_t + t_exc');
    fprintf('%s: stop time prop_t + t_exc = %.4g ns\n', label, evalin('base', 'prop_t + t_exc')*1e9);
else
    set_param(mdl, 'StopTime', num2str(stopTime, '%.10g'));
    fprintf('%s: stop time %.4g ns (end of the CST records)\n', label, stopTime*1e9);
end
if fastRestart
    try
        set_param(mdl, 'FastRestart', 'on');
    catch ME
        warning('Fast Restart not available for %s (%s): full compile per point.', mdl, ME.message);  fastRestart = false;
    end
end
cleanupFR = onCleanup(@() set_param(mdl, 'FastRestart', 'off'));   % leave the model editable
for d = 1:nD
    th = dirs(d, 1);  ph = dirs(d, 2);
    assignin('base', 'theta', th*pi/180);  assignin('base', 'phi', ph*pi/180);
    evalin('base', 'link_update_gains');            % direction gains (tunable); prop_t unchanged (same r)
    try
        out = sim(mdl);
    catch ME
        if ~fastRestart, rethrow(ME); end
        warning('Fast Restart run failed (%s): switching it off.', ME.message);
        set_param(mdl, 'FastRestart', 'off');  fastRestart = false;  out = sim(mdl);
    end
    ti  = out.SimulationMetadata.TimingInfo;        % Simulink's own timing: solve only, no plotting
    tSim = [ti.InitializationElapsedWallTime, ti.ExecutionElapsedWallTime, ti.TotalElapsedWallTime];
    vt = out.yout.getElement(nm{1}).Values;  vp = out.yout.getElement(nm{2}).Values;
    ts = vt.Time(:);  es = [vt.Data(:) vp.Data(:)];

    % CST probe in this direction (the z axis: any phi is the same point, but the components
    % differ, so only an exact phi match is compared there)
    j = [];
    if ~isempty(cst)
        j = find(arrayfun(@(c) abs(c.theta - th) < tol && angd(c.phi, ph) < tol, cst), 1);
    end
    rm = nan(1, 2);  eErr = nan(1, 2);  eRef = nan(1, 2);  pkC = nan(1, 2);  pkS = max(abs(es), [], 1);
    if makePlots
    fig = figure('Name', sprintf('%s th %g ph %g', label, th, ph), 'NumberTitle', 'off', 'Color', 'w');
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
    axs = gobjects(1, 2);
    for c = 1:2
        ax = nexttile(tl);  hold(ax, 'on');  axs(c) = ax;
        h = plot(ax, ts*1e9, es(:, c), '-', 'Color', colSim, 'LineWidth', 3, 'DisplayName', label);
        names = {label};
        if ~isempty(j) && ~isempty(cst(j).e{c})
            tc = cst(j).t{c}(:) - t0;  ec = cst(j).e{c}(:);  [tc, iu] = unique(tc);  ec = ec(iu);
            h(2) = plot(ax, tc*1e9, ec, 'k-', 'LineWidth', 1.1, 'DisplayName', 'CST');  names{2} = 'CST';
            pkC(c) = max(abs(ec));
            if pkC(c) > 1e-3
                tB = max(tc(1), ts(1));  tE = min(tc(end), ts(end));  tg = (tB:2e-12:tE).';
                xs = interp1(ts, es(:, c), tg, 'pchip', 0);  xc = interp1(tc, ec, tg, 'pchip', 0);
                eErr(c) = trapz(tg, (xs - xc).^2);  eRef(c) = trapz(tg, xc.^2);   % error / reference energy (V^2 s)
                rm(c) = eErr(c) / eRef(c);
            end
        end
        grid(ax, 'on');  box(ax, 'on');  ylabel(ax, ['$rE_{' cn{c} '}$ (V)']);
        legend(ax, h, names, 'Location', 'best');
        if ~isnan(rm(c)), title(ax, sprintf('rel. MSE compared to CST %.3g\\%%', 100*rm(c))); end
    end
    % same scale for rE_theta and rE_phi: symmetric y limits from the larger of the two
    ymax = max(arrayfun(@(a) max(abs([a.YLim])), axs));
    linkaxes(axs, 'xy');  ylim(axs(1), [-1 1]*ymax);
    xlabel(tl, 'Time (ns)', 'Interpreter', 'latex');      % (tiledlayout labels ignore the groot default)
    title(tl, sprintf('OOK direct modulation, %s: $\\theta = %g^\\circ$, $\\phi = %g^\\circ$, $r = %g$ m', label, th, ph, r), 'Interpreter', 'latex');
    end
    o = [];  if d == 1, o = out; end
    res(d) = struct('theta', th, 'phi', ph, 'r', r, 't', ts + t0, 'e', es, 'relmse', rm, 'err_energy', eErr, 'ref_energy', eRef, ...
        'peak_sim', pkS, 'peak_cst', pkC, 'sim_time', tSim, 'out', o); %#ok<AGROW>
    fprintf('  %d/%d  theta %g, phi %g: Simulink init %.2f s + solve %.2f s = %.2f s\n', d, nD, th, ph, tSim);
end
st = vertcat(res.sim_time);
if fastRestart, frTxt = 'Fast Restart on'; else, frTxt = 'Fast Restart off'; end
fprintf('%s (%s): %d observation point(s); Simulink init %.1f s + solve %.1f s = %.1f s (mean solve %.2f s per point); loop incl. plotting %.1f s\n', ...
    label, frTxt, nD, sum(st(:, 1)), sum(st(:, 2)), sum(st(:, 3)), mean(st(:, 2)), toc(tAll));

if ~makePlots, return, end          % simulate-only use (e.g. compare_ook_CST_ADS_simulink.m): no table
% ---- summary ----
if isempty(cst), fprintf('\n%s (no CST comparison)\n', label); else, fprintf('\n%s compared to CST\n', label); end
fprintf(' theta   phi   peak rE_theta (V)      peak rE_phi (V)        rel. MSE\n');
fprintf('             Simulink  CST        Simulink  CST        rE_theta   rE_phi\n');
pct = @(x) ternary(isnan(x), '     -  ', sprintf('%7.3f%%', 100*x));
for d = 1:nD
    q = res(d);
    fprintf(' %5g %5g   %8.3g  %8.3g    %8.3g  %8.3g    %s  %s\n', q.theta, q.phi, q.peak_sim(1), q.peak_cst(1), ...
        q.peak_sim(2), q.peak_cst(2), pct(q.relmse(1)), pct(q.relmse(2)));
end
% total rel. MSE over every point and component compared (null CST components, peak <= 1 mV, are excluded):
% sum of error energies / sum of reference energies
Ee = vertcat(res.err_energy);  Er = vertcat(res.ref_energy);
if any(~isnan(Er(:)))
    tot = @(m) sum(Ee(m)) / sum(Er(m));
    okT = ~isnan(Er(:, 1));  okP = ~isnan(Er(:, 2));  okA = ~isnan(Er);
    fprintf('Total rel. MSE compared to CST (%d field(s), zero fields excluded): %.3f %%', nnz(okA), 100*tot(okA));
    if any(okT), fprintf('; rE_theta %.3f %% (%d)', 100*sum(Ee(okT, 1))/sum(Er(okT, 1)), nnz(okT)); end
    if any(okP), fprintf('; rE_phi %.3f %% (%d)', 100*sum(Ee(okP, 2))/sum(Er(okP, 2)), nnz(okP)); end
    fprintf('\n');
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
