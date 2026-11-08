function sim = simulinkRadiatedField(dirs, opt)
%SIMULINKRADIATEDFIELD  Simulink radiated far-field time signals r*E_theta(t), r*E_phi(t)
%   of the Tx-only spherical-mode model (matlab/simulink/build_tx_rx_link_simulink.m) at the
%   observation directions dirs = [theta phi] (deg, one row each), taken from a cache file
%   and, for directions not in it, simulated automatically (the model is built first) and
%   added to the cache.
%
%   sim = simulinkRadiatedField(dirs, opt)
%     sim : struct array (one element per row of dirs): theta, phi (deg), t (s), e (N x 2:
%           [rE_theta rE_phi] in V)
%     opt : struct with fields
%       tx_snp      Touchstone file of the Tx antenna
%       cst_io      CST port-signal export used as the Tx waveform
%       cst_signal  its signal name (e.g. "i1")
%       r           observation distance (m)
%       Rload       internal-port load resistance (Ohm; unused without internal ports)
%       fit_tol_dB  rational-fit tolerance (dB), platform ("simscape" | "rf"), rf_solver, rf_fitopt
%       cache_file  .mat file of the simulated signals
%       mode        "auto"  : use the cache, simulate only the missing directions (default)
%                   "rerun" : simulate every direction again (and rewrite the cache)
%                   "cache" : cache only; error if a direction is missing
%       angle_tol   deg (default 0.1)
%
%   The cache is only used while its settings (Touchstone and waveform files and their
%   dates, signal, r, load, fit tolerance, platform) are those of opt; otherwise it is rebuilt.

if ~isfield(opt, 'mode') || strlength(opt.mode) == 0, opt.mode = "auto"; end
if ~isfield(opt, 'angle_tol'), opt.angle_tol = 0.1; end
opt.mode = string(opt.mode);
tol = opt.angle_tol;
angd = @(a, b) abs(mod(a - b + 180, 360) - 180);

% settings that identify the simulation (a change invalidates the cache)
fs = dir(opt.tx_snp);  fw = dir(opt.cst_io);
key = struct('tx_snp', string(opt.tx_snp), 'snp_date', fs.datenum, 'cst_io', string(opt.cst_io), ...
    'io_date', fw.datenum, 'cst_signal', string(opt.cst_signal), 'r', opt.r, 'Rload', opt.Rload, ...
    'fit_tol_dB', opt.fit_tol_dB, 'platform', string(opt.platform), 'rf_fitopt', string(opt.rf_fitopt));

cache = struct('theta', {}, 'phi', {}, 't', {}, 'e', {});
if isfile(opt.cache_file) && opt.mode ~= "rerun"
    S = load(opt.cache_file, 'key', 'cache');
    if isequaln(S.key, key)
        cache = S.cache;
    else
        fprintf('Simulink cache %s was made with other settings: simulating again.\n', opt.cache_file);
    end
end

inCache = @(C, th, ph) find(arrayfun(@(s) abs(s.theta - th) < tol && angd(s.phi, ph) < tol, C), 1);
miss = false(size(dirs, 1), 1);
for d = 1:size(dirs, 1), miss(d) = isempty(inCache(cache, dirs(d, 1), dirs(d, 2))); end

if any(miss)
    if opt.mode == "cache"
        error('Simulink data missing in %s for (theta, phi) = %s; use mode "auto" to simulate them.', ...
            opt.cache_file, mat2str(dirs(miss, :)));
    end
    fprintf('Simulink: %d of %d direction(s) not simulated yet: building the model and simulating.\n', nnz(miss), numel(miss));
    todo = dirs(miss, :);
    mdl = build_tx_rx_link_simulink(char(opt.tx_snp), '', 'TxOnly', true, 'Waveform', 'cst', ...
        'CSTFile', char(opt.cst_io), 'CSTSignal', char(opt.cst_signal), 'r', opt.r, ...
        'theta_deg', todo(1, 1), 'phi_deg', todo(1, 2), 'TxLoadR', opt.Rload, 'FitTolerance', opt.fit_tol_dB, ...
        'OpenModel', false, 'Platform', char(opt.platform), 'RFSolver', char(opt.rf_solver), ...
        'FitOptions', char(opt.rf_fitopt));
    assignin('base', 'Rload_Tx', opt.Rload);
    wave = evalin('base', 'Vtx_wave');  tend = wave.Time(end);
    for d = 1:size(todo, 1)
        cache(end+1) = simulateDirection(mdl, todo(d, 1), todo(d, 2), opt.r, tend); %#ok<AGROW>
    end
    cdir = fileparts(opt.cache_file);  if ~isempty(cdir) && ~isfolder(cdir), mkdir(cdir); end
    save(opt.cache_file, 'key', 'cache');
    fprintf('Simulink results saved in %s (%d directions)\n', opt.cache_file, numel(cache));
else
    fprintf('Simulink: all %d direction(s) taken from %s\n', size(dirs, 1), opt.cache_file);
end

sim = struct('theta', {}, 'phi', {}, 't', {}, 'e', {});
for d = 1:size(dirs, 1), sim(end+1) = cache(inCache(cache, dirs(d, 1), dirs(d, 2))); end %#ok<AGROW>
end

function s = simulateDirection(mdl, th, ph, r, tend)
% simulate the Tx-only model for the direction (theta, phi) at distance r: r*E_theta, r*E_phi in V
assignin('base', 'r', r);
assignin('base', 'theta', th*pi/180);
assignin('base', 'phi', ph*pi/180);
evalin('base', 'link_update_gains');
pt = evalin('base', 'prop_t');
assignin('base', 't_exc', tend);  set_param(mdl, 'StopTime', 'prop_t + t_exc');   % excitation duration + propagation delay
fprintf('  Simulink theta = %g, phi = %g: simulating to %.4g ns...\n', th, ph, (pt + tend)*1e9);
out = sim(mdl);
vth = out.yout.getElement('Vtheta').Values;
vph = out.yout.getElement('Vphi').Values;
s = struct('theta', th, 'phi', ph, 't', vth.Time, 'e', [vth.Data(:), vph.Data(:)]);
end
