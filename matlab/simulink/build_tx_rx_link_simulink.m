function mdl = build_tx_rx_link_simulink(tx_snp_path, rx_snp_path, varargin)
%BUILD_TX_RX_LINK_SIMULINK  Build a Simulink model of the Tx-to-Rx spherical-mode link.
%   Simulink counterpart (RF Blockset by default, or Simscape) of the ADS link built by
%   python/ads/create_tx_rx_link_internal_ports_TRANSIENT.py.
%
%   mdl = build_tx_rx_link_simulink(tx_snp_path, rx_snp_path, Name, Value, ...)
%
%   tx_snp_path, rx_snp_path : Touchstone files of the Tx / Rx spherical-mode networks (port names
%   in the "! Port[N] Name" header: Feed, TM_L<l>_M<m> / TE_..., internal ports).  The model and
%   <model>_vars.mat (all workspace variables it needs; loaded by the model's PreLoadFcn) are saved
%   in 'OutputDir' (default <repo>/output/simulink).
%
%   Topology (same as the ADS schematic):
%     Tx antenna : Simscape "S-Parameter Rational Fit" block, one +/- terminal
%                  pair per port (OneRefPinPerPort).
%       Feed       : Rsrc_Tx in series with the Tx waveform source, across the port.
%       mode port  : eta0 to ground, reference grounded, voltage sensed.
%       other port : (internal ports, e.g. load1) Z0 across the port.
%     Tx modes -> Transport Delay(prop_t) -> Gtx (VCVS gains, from the gradYlm_*
%       variables) -> Rrot (VectorFieldRotate) -> Vtheta, Vphi
%     Vtheta,Vphi -> 1/Lind + integrator (the mu0*r inductor) -> Itheta, Iphi
%     Itheta,Iphi -> Grx (CCCS gains) -> current sources into the Rx mode ports.
%     Rx antenna : second S-Parameter block; Feed/internal ports = Z0 across the port.
%   Feed and internal ports are NOT referenced to ground: the load/source sits
%   between the port's + and - terminals and the - terminal floats.  An optional
%   gmin (Rgmin, default 1e12 Ohm, as in the ADS rational-fit netlists) to
%   ground keeps the solver well posed; RefGmin = Inf leaves the node fully floating.
%
%   The Tx mode voltages (v_Tx_<mode>) are wired to the Mux and the Rx mode currents (i_Rx_<mode>)
%   from the Demux; wires from the far edge of an S-parameter block run over the top, and the routes are
%   chosen so that no wires cross.
%
%   TX WAVEFORM (comes from the base workspace, variable Vtx_wave):
%     'Waveform','sine'       sine source: Vamp, f0 (default)
%     'Waveform','cst'        CST ASCII export read with importCSTSignal.m (as in
%                             cst_to_ads_mdf.m): 'CSTFile', 'CSTSignal' (index or
%                             name like 'i1'; default 1).  W^(1/2) exports become a
%                             generator: incident-wave voltage V+ = a*sqrt(Z0 of the
%                             export) (the ADS Vinc / Vi_1), source impedance = that Z0
%                             (50 Ohm), open-circuit voltage 2*V+.  V / port-voltage
%                             exports become an ideal voltage source (series R = 1e-6 Ohm).
%     'Waveform','workspace'  user-defined signal in the variable named by
%                             'WaveformVar' (default 'Vtx_wave'), given as a
%                             timeseries, an N x 2 matrix [t(s) v(V)], a struct with
%                             fields t, v (optional type), or a function handle
%                             @(t) in volts (sampled over 'WaveformDuration' s at
%                             'WaveformDt').  'SourceType' = 'generator' (the signal is
%                             the incident-wave voltage V+; open-circuit voltage 2*V+,
%                             series resistance FeedZ; default) or 'voltage' (ideal source).
%   The waveform is stored as the timeseries Vtx_wave (From Workspace block), the series
%   resistance as Rsrc_Tx, the Rx Feed termination as Rfeed_Rx and the factor between the
%   waveform and the source open-circuit voltage as Vsrc_scale (2 generator, 1 voltage);
%   edit them in the workspace and re-run.  'FeedZ' (Ohm) sets Rsrc_Tx and Rfeed_Rx: default
%   the Z0 of the CST export (50 Ohm) for 'cst', else the Touchstone Z0.  NOTE the Touchstone
%   reference impedance (Z0 = eta0 = 376.73) is a property of the S-parameter block; the Feed
%   source/termination impedance is separate.
%
%   WORKSPACE VARIABLES (the ADS VAR block, SI units; ADS mu0=1256.6 nH/m is
%   1256.6e-9 H/m here): eps0, mu0, eta0, Z0, r, theta, phi, rot_angle, Ux, Uy, Uz,
%   a_Tx, a_Rx, plus Rgmin, c0, f0, Vamp, Vtx_wave, Rsrc_Tx and the derived ADS VARs
%   prop_t, Lind (=mu0*r), Ftheta_l<l>_m<m> / Fphi_l<l>_m<m> (+ _prime; the ADS Ftheta_Fphi_vars entries, read from the
%   Wolfram-exported FTheta_ADS_Discrete.txt / FPhi_ADS_Discrete.txt in the folder given by the workspace variable Ylm_data_dir,
%   default <repo>/mathematica/tables), theta_prime, phi_prime, R11..R33, magU,
%   gradYlm_theta[_prime]_l<l>_m<m>, gradYlm_phi[_prime]_l<l>_m<m>, Rrot.  After
%   changing r/theta/phi/rot_angle/Ux/Uy/Uz/a_Tx/a_Rx/mu0/c0 run link_update_gains.
%   The Gain blocks use expressions that reference the gradYlm_* variables by name
%   (like the ADS VCVS/CCCS G fields).
%
%   Nets are named like the ADS nets: FeedTx, FeedRx, Tx_<mode>, Rx_<mode>,
%   Tx_<port>, Rx_<port>.
%   Logged outputs (Outport blocks, out.yout): Vtheta, Vphi; the differential port
%   voltages FeedTx, FeedRx, Tx_<port>, Rx_<port> (V(+)-V(-)) for the feed/internal
%   ports (the absolute voltage of the - terminal is not logged: only the differential voltage is used).
%
%   INTERNAL-PORT TERMINATIONS: the resistance across every Tx / Rx internal port (e.g. load1) is
%   the workspace variable Rload_Tx / Rload_Rx (default Z0, as the ADS Term).  'TxLoadR', 'RxLoadR'
%   (Ohm) set them: 1e9 ~ open circuit, 1e-6 ~ short circuit.
%
%   Z0 (priority): 'Z0' option > the CST export (only when replicating CST/ADS results, i.e. with
%   'Waveform','cst' or 'CSTParamFile') > the base-workspace variable Z0 (the normal case with a
%   workspace waveform: set Z0 yourself) > the Touchstone header.  Z0 is a free variable and
%   need not agree with the Touchstone.
%
%   Z0 AND a FROM CST: the '#Parameters = {... a=39.71; Z0=50 ...}' line of a CST export
%   ('CSTFile', or 'CSTParamFile' for any CST time-signal / far-field export, e.g. with a
%   workspace waveform) gives the port impedance Z0 (workspace variable Z0, the default
%   Feed source/termination and internal-port load resistance) and the enclosing-sphere
%   radius a (a_Tx; 'CSTLengthUnit' converts it to metres, default 1e-3).  The Touchstone
%   header 'a=...m' must agree with it (error if they differ by > 0.5 %, i.e. beyond the header's
%   rounding); Z0 is a free variable and need not match the Touchstone.  Without a CST file the
%   header values are used.  The S-parameter block's reference impedance
%   (Z0 = eta0 = 376.73 in the .snp) is a property of the Touchstone file and unchanged.
%
%   Other options:
%     'ModelName' (default Tx<ant>_Rx<ant>, plus _internal_ports if the antennas have internal ports: the ADS cell name), 'r' 10, 'theta_deg' 90,
%     'phi_deg' 0, 'rot_angle_deg' 0, 'axis' [1 1 1], 'SourceFreq' (Hz, sine),
%     'SourceAmplitude' (V, sine), 'StopTime', 'MaxStep' (default Nyquist/6 = 1/(12*fmax)), 'FitTolerance' (dB, -40),
%     'RefGmin' (Ohm, 1e12; Inf = none), 'OutputDir' (default <repo>/output/simulink, created if missing), 'OpenModel' (true),
%     'FontSize' (14), 'ShowCaptions' (true), 'RFSolver' ('auto': RF platform uses variable-step ode45, or ode15s
%     when 'LoadSwitch' is true because switches make the network stiff; also any variable-step solver name or 'fixed').
%     Fit accuracy: 'FitTolerance' (dB, e.g. -40: the fit stops at the fewest poles reaching it) with 'MaxPoles' as the cap.
%     'Platform' ('rf' default = RF Blockset circuit with an S-parameters block, suffix _RF, the setup used for
%       the paper results; 'simscape' = Simscape Electrical with the S-Parameter Rational Fit block).
%     'ImportFit' (Simscape platform): a rational object, or a .mat file holding one, used as is (nothing is fitted), as ADS imports
%     a fit. Note: the rational object shares one pole set across all entries; non-shared fits ('FitOptions' 'Fit individually' or
%     'Share poles by columns') exist only in the RF platform's S-parameters block.
%     'SParamBlock' (Simscape platform; the RF platform always uses 'rf'): 'rational' (default) = Simscape S-Parameter Rational Fit fed by a rational fit made here with
%       'FitTolerance' (-40 dB) and 'MaxPoles' (100); 'rf' = RF Blockset S-parameters block reading the Touchstone file
%       itself, with 'SParamModeling' ('Time domain' | 'Frequency domain'), 'FitTolerance', 'MaxPoles', 'FitOptions'
%       ('Share all poles' | 'Share poles by columns' | 'Fit individually')).
%     'TxOnly' (false): build only the Tx side (no Rx S-parameter block, 1/Lind, integrator, Grx):
%       Tx modes -> delay (r - a_Tx)/c0 -> Gtx -> Rrot -> Vtheta, Vphi = r*E_theta, r*E_phi.  Used to
%       compare radiated fields with CST/ADS probes (no Rx de-embedding delay a_Rx/c0).  The Rx file
%       argument may be '' ; default ModelName Tx<ant>[_internal_ports] (the ADS name without the Rx antenna), e.g.
%       mdl = build_tx_rx_link_simulink(snp, '', 'TxOnly', true, 'Waveform', 'cst', ...);
%
%   EXAMPLES (copy the lines between %{ and %} into the MATLAB console):
%{
snp = fullfile(repo_root, 'output', 'SnP', 'S_dom_Q_dipole_cylinder_edgeport_full.s6p');
mdl = build_tx_rx_link_simulink(snp, snp);                                  % sine source, r = 10 m
mdl = build_tx_rx_link_simulink(snp, snp, 'r', 5, 'theta_deg', 60, 'rot_angle_deg', 15);
mdl = build_tx_rx_link_simulink(snp, snp, 'Waveform', 'cst', ...
    'CSTFile', 'C:\path\to\CST\<antenna>\Export\time_signals\port2_OC\iowaves_gaussian_2GHz_3GHz.txt', 'CSTSignal', 'i1');
mywave = @(t) exp(-((t-1.5e-9)/0.4e-9).^2) .* sin(2*pi*1e9*t);
mdl = build_tx_rx_link_simulink(snp, snp, 'Waveform', 'workspace', 'WaveformVar', 'mywave', 'WaveformDuration', 3e-9);
out = sim(mdl);                                                              % logged signals in out.yout
v = out.yout.getElement('FeedRx').Values;  plot(v.Time*1e9, v.Data)
%}

p = inputParser;
addParameter(p, 'ModelName', '');
addParameter(p, 'r', 10);
addParameter(p, 'theta_deg', 90);
addParameter(p, 'phi_deg', 0);
addParameter(p, 'rot_angle_deg', 0);
addParameter(p, 'axis', [1 1 1]);
addParameter(p, 'SourceFreq', []);
addParameter(p, 'SourceAmplitude', 1);
addParameter(p, 'StopTime', []);
addParameter(p, 'MaxStep', []);
addParameter(p, 'FitTolerance', -40);
addParameter(p, 'RefGmin', 1e12);
addParameter(p, 'OutputDir', '');
addParameter(p, 'OpenModel', true);
addParameter(p, 'Waveform', 'sine');
addParameter(p, 'CSTFile', '');
addParameter(p, 'CSTSignal', 1);
addParameter(p, 'WaveformVar', 'Vtx_wave');
addParameter(p, 'SourceType', 'generator');
addParameter(p, 'WaveformDuration', []);
addParameter(p, 'WaveformDt', []);
addParameter(p, 'Z0', []);
addParameter(p, 'FeedZ', []);
addParameter(p, 'CSTParamFile', '');
addParameter(p, 'CSTLengthUnit', 1e-3);
addParameter(p, 'TxLoadR', []);
addParameter(p, 'RxLoadR', []);
addParameter(p, 'FontSize', 14);
addParameter(p, 'ShowCaptions', true);
addParameter(p, 'TxOnly', false);
addParameter(p, 'SParamBlock', 'rational');
addParameter(p, 'SParamModeling', 'Time domain');
addParameter(p, 'MaxPoles', 100);
addParameter(p, 'FitOptions', 'Share all poles');   % RF platform: 'Share all poles' (default, fewest poles) | 'Share poles by columns' | 'Fit individually'; if the simulation is inaccurate, try by columns or individually (or tighten FitTolerance)
addParameter(p, 'RFSolver', 'auto');          % RF platform solver: 'auto' = variable-step ode45 (ode15s when LoadSwitch is true: switches make it stiff), or any variable-step solver name, or 'fixed' (FixedStepDiscrete). MaxStep (RF step) is the max step.
addParameter(p, 'ImportFit', []);              % Simscape platform: rational object (or .mat file holding one) used as is instead of fitting, as ADS imports a fit; else 'FitTolerance' (dB) / 'MaxPoles' set the fit accuracy
addParameter(p, 'RFRationalVar', '');         % RF platform: name of a base-workspace rational fit fed to the RF block ('Rational model') instead of its own fit:
                                              % a rational object (rational) or an N x N rfmodel.rational array (rationalfit, any pole sharing)
addParameter(p, 'RFRationalZ0', 50);          % RF platform with RFRationalVar: reference impedance (ohm) of that rational fit. The block's rational model is
                                              % always referenced to 50 ohm, so for any other value an ideal transformer (turns ratio n = sqrt(Z0/50)) is put
                                              % on every port: the network then has the fitted S-parameters referenced to Z0 (e.g. eta0 = 376.73 ohm)
addParameter(p, 'LoadSwitch', false);       % RF platform, Tx internal ports: voltage-controlled switch (RF Variable Resistor) instead of the fixed load
addParameter(p, 'ControlVar', 'Vs_ctrl');    % workspace variable (timeseries / [t v]) of the switch control voltage (ideal source, scale 1)
addParameter(p, 'Platform', 'rf');           % 'rf' (RF Blockset circuit, S-parameters block; used for the paper) | 'simscape' (Simscape Electrical)
addParameter(p, 'SParamSubsystem', []);      % true: each S-parameter block wrapped in a subsystem (SnP_Tx / SnP_Rx, terminals named
                                              % <port>+ / <port>-); [] = true on the RF platform, false on Simscape
parse(p, varargin{:});
o = p.Results;
o.Waveform = lower(char(o.Waveform));
o.Platform = lower(char(o.Platform));
if strcmp(o.Platform, 'rf'), o.SParamBlock = 'rf'; end
if strcmp(o.Platform, 'rf') && ~strcmpi(o.SParamModeling, 'Time domain')
    % the block's frequency-domain model does not give usable results for a baseband (0 Hz carrier) transient
    warning('build_tx_rx_link_simulink:RFFreqDomain', 'RF Blockset frequency-domain S-parameter modelling is not used (unstable baseband transient): using Time domain.');
    o.SParamModeling = 'Time domain';
end
o.SParamBlock = lower(char(o.SParamBlock));
useRF = strcmp(o.SParamBlock, 'rf');
isRFplat = strcmp(o.Platform, 'rf');
if o.TxOnly || isempty(rx_snp_path), rx_snp_path = tx_snp_path; end     % Tx-only: no Rx file needed (pass '')

here = fileparts(mfilename('fullpath'));
if exist('repo_root', 'file') ~= 2, run(fullfile(here, '..', '..', 'setup.m')); end   % toolbox path (matlab/sphmodes: gradYlm)
if isempty(o.OutputDir), o.OutputDir = fullfile(repo_root, 'output', 'simulink'); end
o.OutputDir = abs_path(char(o.OutputDir));
if ~isfolder(o.OutputDir), mkdir(o.OutputDir); end
% the RF Blockset S-parameters block stores the Touchstone path: make it absolute so the model
% simulates from any current folder
tx_snp_path = abs_path(char(tx_snp_path));
if ~isempty(rx_snp_path), rx_snp_path = abs_path(char(rx_snp_path)); end

%% ---------------------------------------------------------- read files
tx = read_snp_info(tx_snp_path);
rx = read_snp_info(rx_snp_path);
fprintf('Tx: %d ports (%d modes);  Rx: %d ports (%d modes).\n', ...
    numel(tx.names), sum(strcmp(tx.kind, 'mode')), numel(rx.names), sum(strcmp(rx.kind, 'mode')));
if ~any(strcmp(tx.kind, 'feed')) || ~any(strcmp(rx.kind, 'feed'))
    error('Both Touchstone files need a port named "Feed" in their "! Port[N] Name" header.');
end
if isempty(o.ModelName)
    % ADS cell naming: Tx<ant>[_Rx<ant>][_internal_ports]; the suffix only when internal ports are present
    sfx = '';
    if any(strcmp(tx.kind, 'internal')) || (~o.TxOnly && any(strcmp(rx.kind, 'internal'))), sfx = '_internal_ports'; end
    if o.TxOnly
        o.ModelName = matlab.lang.makeValidName(sprintf('Tx%s%s', tx.ant, sfx));
    else
        o.ModelName = matlab.lang.makeValidName(sprintf('Tx%s_Rx%s%s', tx.ant, rx.ant, sfx));
    end
    if strcmp(o.Platform, 'rf'), o.ModelName = [o.ModelName '_RF']; end     % RF Blockset version
end
mdl = o.ModelName;
fmax = min(tx.fmax, rx.fmax);  fmin = max(tx.fmin, rx.fmin);
if isempty(o.SourceFreq), o.SourceFreq = min(max(1e9, fmin), fmax); end

%% ------------------------------------------ base-workspace variables (ADS VAR block)
Z0val = 376.73;
if ~isnan(tx.Z0), Z0val = tx.Z0; elseif ~isnan(rx.Z0), Z0val = rx.Z0; end
pfile = char(o.CSTParamFile);  if isempty(pfile), pfile = char(o.CSTFile); end
cstp = struct('Z0', NaN, 'a', NaN);
if ~isempty(pfile)
    assert(isfile(pfile), 'CST file not found: %s', pfile);
    cstp = read_cst_params(pfile, o.CSTLengthUnit);
    fprintf('CST parameters (%s): Z0 = %g Ohm, a = %g m\n', pfile, cstp.Z0, cstp.a);
end
z0src = 'Touchstone header';
if ~isempty(o.Z0)
    Z0val = o.Z0;  z0src = 'option Z0';
elseif ~isnan(cstp.Z0)
    Z0val = cstp.Z0;  z0src = 'CST export';
elseif evalin('base', 'exist(''Z0'', ''var'') && isnumeric(Z0) && isscalar(Z0)')
    Z0val = evalin('base', 'Z0');  z0src = 'workspace variable Z0';
end
fprintf('Z0 = %g Ohm (%s)\n', Z0val, z0src);
a_Tx = tx.a;  if isnan(a_Tx), a_Tx = 0.075; end
if ~isnan(cstp.a)
    if ~isnan(tx.a) && abs(tx.a - cstp.a) > 0.005*cstp.a
        error(['Touchstone header a = %g m (%s) and CST export a = %g m (%s) disagree: the spherical-mode ' ...
            'data and the CST signals must use the same enclosing-sphere radius.'], tx.a, tx_file_name(tx_snp_path), ...
            cstp.a, tx_file_name(pfile));
    end
    a_Tx = cstp.a;                       % the CST value is the more precise (the header is rounded)
end
a_Rx = rx.a;  if isnan(a_Rx), a_Rx = a_Tx;  end
if strcmp(rx_snp_path, tx_snp_path), a_Rx = a_Tx; end
vars = struct('eps0', 8.8542e-12, 'eta0', 376.73, 'Z0', Z0val, 'Rgmin', o.RefGmin, 'c0', 299792458, ...
    'mu0', 1256.6e-9, 'r', o.r, 'theta', o.theta_deg*pi/180, 'phi', o.phi_deg*pi/180, ...
    'rot_angle', o.rot_angle_deg*pi/180, 'Ux', o.axis(1), 'Uy', o.axis(2), 'Uz', o.axis(3), ...
    'a_Tx', a_Tx, 'a_Rx', a_Rx, 'f0', o.SourceFreq, 'Vamp', o.SourceAmplitude);
fn = fieldnames(vars);
for i = 1:numel(fn), assignin('base', fn{i}, vars.(fn{i})); end
cfg.tx_modes = tx.modes;  cfg.rx_modes = rx.modes;
cfg.tx_needs_corr = ~tx.ka_is_one;  cfg.rx_needs_corr = ~rx.ka_is_one && ~o.TxOnly;   % Tx-only: delay (r - a_Tx)/c0
cfg.tx_only = o.TxOnly;
assignin('base', 'link_cfg', cfg);
link_update_gains;
prop_t = evalin('base', 'prop_t');
fprintf('prop_t = %.4g ns (Tx ka_factor=1: %d, Rx ka_factor=1: %d)\n', ...
    prop_t*1e9, tx.ka_is_one, rx.ka_is_one);

%% ------------------------------------------------------- Tx waveform
[wave, srcType, waveDesc, feedZ_wave] = get_waveform(o, Z0val, fmax);     % wave = timeseries or []
if isempty(o.FeedZ), Zfeed = feedZ_wave; else, Zfeed = o.FeedZ; end
if isnan(Zfeed), Zfeed = Z0val; end
if strcmp(srcType, 'generator'), Rsrc = Zfeed; Vscale = 2; else, Rsrc = 1e-6; Vscale = 1; end
assignin('base', 'Rsrc_Tx', Rsrc);  assignin('base', 'Rfeed_Rx', Zfeed);  assignin('base', 'Vsrc_scale', Vscale);
if isempty(o.TxLoadR), assignin('base', 'Rload_Tx', Z0val); else, assignin('base', 'Rload_Tx', o.TxLoadR); end
if isempty(o.RxLoadR), assignin('base', 'Rload_Rx', Z0val); else, assignin('base', 'Rload_Rx', o.RxLoadR); end
if ~isempty(wave)
    assignin('base', 'Vtx_wave', wave);
    tend = wave.Time(end);
    dtw = median(diff(wave.Time));
else
    tend = 10e-9;  dtw = Inf;
end
assignin('base', 'SW_Ron', 0.01);  assignin('base', 'SW_Roff', 1e6);      % switch on / off resistance (SPICE switch RON, ROFF)
assignin('base', 't_exc', tend);                         % excitation duration at the feed
if isempty(o.StopTime), stopStr = 'prop_t + t_exc'; else, stopStr = num2str(o.StopTime, '%.6g'); end   % stop time = excitation duration + propagation delay (workspace variables)
if isempty(o.MaxStep)
    % Nyquist step of the highest Touchstone frequency, divided by 6: caps the variable-step solver and, on the
    % RF platform, is the Configuration block's time step (so it must resolve the waveform)
    o.MaxStep = min(1/(2*fmax)/6, dtw);
end
fprintf('Waveform: %s\n', waveDesc);

%% --------------------------------------------------------- rational fit
if useRF
    fprintf('S-parameter block: RF Blockset S-parameters (%s, fit tol %g dB, max %d poles, %s).\n', ...
        o.SParamModeling, o.FitTolerance, min(o.MaxPoles, 99), o.FitOptions);
    ZrefTx = NaN;  ZrefRx = NaN;
else
    if ~isempty(o.ImportFit)
        % imported fit (as in ADS): used as is, nothing is fitted here
        ratTx = import_fit(o.ImportFit);
        ZrefTx = sparameters(tx_snp_path).Impedance;
        fprintf('Rational fit (Tx): imported, %d poles, %d ports, ErrDB %.1f.\n', ratTx.NumPoles, ratTx.NumPorts, ratTx.ErrDB);
    else
        fprintf('Rational fit (Tx), tolerance %g dB, max %d poles...\n', o.FitTolerance, o.MaxPoles);
        [ratTx, ZrefTx] = fit_rational(tx_snp_path, o.FitTolerance, o.MaxPoles);
        fprintf('  -> %d poles, ErrDB %.1f.\n', ratTx.NumPoles, ratTx.ErrDB);
    end
    if strcmp(rx_snp_path, tx_snp_path)
        ratRx = ratTx;  ZrefRx = ZrefTx;
    else
        fprintf('Rational fit (Rx)...\n');
        [ratRx, ZrefRx] = fit_rational(rx_snp_path, o.FitTolerance, o.MaxPoles);
    end
    assignin('base', 'ratTx', ratTx);  assignin('base', 'ratRx', ratRx);
end
rfo = [];
if useRF
    mp = o.MaxPoles;
    if mp > 99                                   % the RF Blockset block accepts at most 99 poles
        warning('build_tx_rx_link_simulink:MaxPoles', 'RF Blockset S-parameters block allows at most 99 poles: using 99 instead of %d.', mp);
        mp = 99;
    end
    rfo = struct('file', '', 'mode', o.SParamModeling, 'tol', o.FitTolerance, 'maxpoles', mp, 'fitopt', o.FitOptions, 'ratvar', o.RFRationalVar, ...
        'ratz0', o.RFRationalZ0);
end

%% ----------------------------------------------------------- model
if bdIsLoaded(mdl), close_system(mdl, 0); end
% the .slx of the previous build is on the path under the same name: that is expected, no warning
wShadow = warning('off', 'Simulink:Engine:MdlFileShadowing');
new_system(mdl);
warning(wShadow);
if o.OpenModel, open_system(mdl); end

fs = o.FontSize;
ROW = max(160, 11*fs);  Y0 = 400;  W_SP = 140;
TAGW = 18*fs;                                  % width of Goto/From tag blocks
SUPW = 190 + TAGW;                             % horizontal room used by one port's support circuit
X_TXSP = 1000 + SUPW;
nTx = numel(tx.names);  nRx = numel(rx.names);
opt = struct('fs', fs, 'TAGW', TAGW, 'gmin', isfinite(o.RefGmin), 'wave', ~isempty(wave), 'ROW', ROW, ...
    'rf', isRFplat, 'step', o.MaxStep, 'sw', o.LoadSwitch && isRFplat, 'ctrlvar', o.ControlVar);

rfoT = rfo;  if useRF, rfoT.file = tx_snp_path; end
sTx = add_sparam(mdl, 'SnP1', X_TXSP, Y0, W_SP, ROW, nTx, 'ratTx', ZrefTx, rfoT);
TxS = build_side(mdl, 'Tx', tx, sTx, opt);
yTopCh = Y0 - 150;                             % first channel above the S-parameter blocks for wires that go over the top
ymid = Y0 + sTx.height/2;
PM = max(34, 2.4*fs);                          % pitch of the Mux / Demux pins
nT = numel(TxS.mode_names);
Ptx = port_positions(TxS.ports);
Ltx = hub_layout(Ptx, TxS.sgn, +1, yTopCh, model_edge(mdl, +1));    % wires Tx mode ports -> Mux, routed without crossings
xmux = Ltx.xpin_min + 60;
ytop_t = ymid - nT*PM/2;
mux = addb(mdl, 'Mux_TxModes', 'simulink/Signal Routing/Mux', [xmux ytop_t xmux+8 ytop_t+nT*PM], 'Inputs', num2str(nT));
pinsT = port_positions(arrayfun(@(k) pH(mux).Inport(k), (1:nT).'));
WT = hub_waypoints(Ltx, pinsT, +1);
for k = 1:nT
    pn = Ltx.order(k);
    wire(mdl, TxS.ports{pn}, pH(mux).Inport(k), ['v_Tx_' TxS.mode_names{pn}], WT{pn}(2:end-1, :));
end
nFarRx = 0;
if ~o.TxOnly
    CH_W = 1900;                               % width reserved for the signal chain
    nRm = sum(strcmp(rx.kind, 'mode'));
    X_RXSP = xmux + CH_W + SUPW + 40*nRm + 250;
    rfoR = rfo;  if useRF, rfoR.file = rx_snp_path; end
    sRx = add_sparam(mdl, 'SnP2', X_RXSP, Y0, W_SP, ROW, nRx, 'ratRx', ZrefRx, rfoR);
    RxS = build_side(mdl, 'Rx', rx, sRx, opt);
    nR = numel(RxS.mode_names);
    Lrx = hub_layout(port_positions(RxS.ports), RxS.sgn, -1, yTopCh, model_edge(mdl, -1));   % wires Demux -> Rx mode ports
    nFarRx = Lrx.nF;
end
xx = xmux + 8 + 150;
td = addb(mdl, 'PropDelay', 'simulink/Continuous/Transport Delay', [xx ymid-25 xx+100 ymid+25], ...
    'DelayTime', 'prop_t', 'InitialOutput', '0', 'BufferSize', '65536');
xx = xx + 190;
rx_modes_ord = rx.modes;  if ~o.TxOnly, rx_modes_ord = rx.modes(Lrx.order); end     % Demux output order
[gtx_expr, grx_expr] = gain_expressions(tx.modes(Ltx.order), rx_modes_ord);                % Mux input order
gtx = addb(mdl, 'VCVS_TxGain_inst', 'simulink/Math Operations/Gain', [xx ymid-25 xx+120 ymid+25], ...
    'Gain', gtx_expr, 'Multiplication', 'Matrix(K*u)');
xx = xx + 190;
rot = addb(mdl, 'VectorFieldRotate_inst', 'simulink/Math Operations/Gain', [xx ymid-25 xx+120 ymid+25], ...
    'Gain', 'Rrot', 'Multiplication', 'Matrix(K*u)');
xx = xx + 190;
dmx = addb(mdl, 'Demux_VthVph', 'simulink/Signal Routing/Demux', [xx ymid-40 xx+8 ymid+40], 'Outputs', '2');
oVth = addb(mdl, 'Vtheta', 'simulink/Sinks/Out1', [xx+90 ymid-110 xx+120 ymid-96]);
oVph = addb(mdl, 'Vphi',   'simulink/Sinks/Out1', [xx+90 ymid+96 xx+120 ymid+110]);
xx = xx + 190;
if ~o.TxOnly
    gL = addb(mdl, 'L1_L2_1overLind', 'simulink/Math Operations/Gain', [xx ymid-25 xx+90 ymid+25], 'Gain', '1/Lind');
    xx = xx + 150;
    integ = addb(mdl, 'L1_L2_Integrator', 'simulink/Continuous/Integrator', [xx ymid-25 xx+50 ymid+25]);
    xx = xx + 150;
    grx = addb(mdl, 'CCCS_RxAperture_inst', 'simulink/Math Operations/Gain', [xx ymid-25 xx+120 ymid+25], ...
        'Gain', grx_expr, 'Multiplication', 'Matrix(K*u)');
    xx = xx + 190;
end

l = add_line(mdl, pH(mux).Outport(1), pH(td).Inport(1));      set_param(l, 'Name', 'Tx_modes');
add_line(mdl, pH(td).Outport(1),  pH(gtx).Inport(1));
l = add_line(mdl, pH(gtx).Outport(1), pH(rot).Inport(1));     set_param(l, 'Name', 'Vtheta_Vphi_prime');
l = add_line(mdl, pH(rot).Outport(1), pH(dmx).Inport(1));     set_param(l, 'Name', 'Vtheta_Vphi');
l = add_line(mdl, pH(dmx).Outport(1), pH(oVth).Inport(1)); set_param(l, 'Name', 'Vtheta');
l = add_line(mdl, pH(dmx).Outport(2), pH(oVph).Inport(1)); set_param(l, 'Name', 'Vphi');
if ~o.TxOnly
add_line(mdl, pH(rot).Outport(1), pH(gL).Inport(1));
add_line(mdl, pH(gL).Outport(1),  pH(integ).Inport(1));
l = add_line(mdl, pH(integ).Outport(1), pH(grx).Inport(1));   set_param(l, 'Name', 'Itheta_Iphi');

ytop_r = ymid - nR*PM/2;
if nR == 1
    pinsR = port_positions(pH(grx).Outport(1));
else
    dmr = addb(mdl, 'Demux_RxModes', 'simulink/Signal Routing/Demux', [xx ytop_r xx+8 ytop_r+nR*PM], 'Outputs', num2str(nR));
    add_line(mdl, pH(grx).Outport(1), pH(dmr).Inport(1));
    pinsR = port_positions(arrayfun(@(k) pH(dmr).Outport(k), (1:nR).'));
end
if -pinsR(1, 1) < Lrx.xpin_min
    warning('build_tx_rx_link_simulink:RxCorridor', 'Rx wire corridor is narrow: some Rx wires may touch.');
end
WR = hub_waypoints(Lrx, pinsR, -1);
for k = 1:nR
    pn = Lrx.order(k);
    if nR == 1, src = pH(grx).Outport(1); else, src = pH(dmr).Outport(k); end
    wire(mdl, src, RxS.ports{pn}, ['i_Rx_' RxS.mode_names{pn}], flipud(WR{pn}(2:end-1, :)));
end
end   % ~TxOnly

%% ----------------------------------------------------- solver settings
if isRFplat
    % RF Blockset: one Configuration block (baseband, carrier 0 Hz) on the Tx network (and one per separate Rx network)
    add_rf_cfg(mdl, TxS.gnd_for_solver, 'Config_Tx', o.MaxStep);
    if ~o.TxOnly, add_rf_cfg(mdl, RxS.gnd_for_solver, 'Config_Rx', o.MaxStep); end
    rfs = char(o.RFSolver);
    if strcmpi(rfs, 'auto'), rfs = 'ode45'; if o.LoadSwitch, rfs = 'ode15s'; end, end   % stiff solver when the network has a switch
    if strcmpi(rfs, 'fixed')
        set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepDiscrete', 'FixedStep', num2str(o.MaxStep, '%.6g'));
    else
        set_param(mdl, 'SolverType', 'Variable-step', 'Solver', rfs, 'MaxStep', num2str(o.MaxStep, '%.6g'));
    end
    set_param(mdl, 'StopTime', stopStr, ...
        'SaveOutput', 'on', 'SaveFormat', 'Dataset', 'ReturnWorkspaceOutputs', 'on', 'SaveTime', 'on');
else
    add_solver_cfg(mdl, TxS.gnd_for_solver, 'Solver_Tx');
    if ~o.TxOnly, add_solver_cfg(mdl, RxS.gnd_for_solver, 'Solver_Rx'); end
    set_param(mdl, 'Solver', 'ode23t', 'StopTime', stopStr, ...
        'MaxStep', num2str(o.MaxStep, '%.6g'), 'RelTol', '1e-4', ...
        'SaveOutput', 'on', 'SaveFormat', 'Dataset', 'ReturnWorkspaceOutputs', 'on', ...
        'SaveTime', 'on');
end

%% ------------------- antenna subsystems: S-parameter block + its port circuits
% everything on the antenna's side of the signal chain goes into one subsystem per antenna: the S-parameter block,
% the mode-port resistors, grounds, gmin resistors, the feed source (and its scaling), the load / switch (and its
% control mapping), the sensors and the solver / Configuration block.  The workspace sources (From Workspace), the
% Mux / Demux of the signal chain and the root Outports stay outside.
wrapSP = o.SParamSubsystem;  if isempty(wrapSP), wrapSP = isRFplat; end
txParent = mdl;  rxParent = mdl;
if wrapSP
    txParent = wrap_antenna(mdl, 'SnP1', 'Antenna_Tx', fs, tx);
    if ~o.TxOnly, rxParent = wrap_antenna(mdl, 'SnP2', 'Antenna_Rx', fs, rx); end
    % createSubsystem drops the names of lines that now leave a subsystem: restore each root Outport's signal
    % name (FeedTx, FeedRx, Tx_<port>, ...) so out.yout.getElement(<name>) works as on the unwrapped model
    for ob = reshape(find_system(mdl, 'SearchDepth', 1, 'BlockType', 'Outport'), 1, [])
        pob = get_param(ob{1}, 'PortHandles');  lob = get_param(pob.Inport(1), 'Line');
        if lob > 0 && isempty(get_param(lob, 'Name')), set_param(lob, 'Name', get_param(ob{1}, 'Name')); end
    end
end

%% -------------------------------------------- fonts, captions, labels
xTxL = X_TXSP - SUPW;  xTxR = X_TXSP + W_SP + SUPW;
cy = Y0 - 36*max(Ltx.nF, nFarRx) - 90;        % captions sit above the wire channels
if wrapSP                                      % compact top level: captions just above the blocks
    bp_ = get_param(find_system(mdl, 'SearchDepth', 1, 'Type', 'Block'), 'Position');
    bp_ = vertcat(bp_{:});  cy = min(bp_(:, 2)) - 60;
    xTxL = min(bp_(:, 1));  xmux = xTxL + 900;
end
if o.ShowCaptions
    note(mdl, sprintf('Tx antenna: %s   (%s)', tx.ant, tx_file_name(tx_snp_path)), xTxL, cy-300, fs+6, true);
    if ~o.TxOnly
        note(mdl, sprintf('Rx antenna: %s   (%s)', rx.ant, tx_file_name(rx_snp_path)), X_RXSP-SUPW, cy-300, fs+6, true);
    end
    note(mdl, sprintf(['Workspace variables (ADS VAR block): eps0, mu0, eta0, Z0, r, theta, phi, rot_angle, Ux, Uy, Uz, a_Tx, a_Rx,\n' ...
        'prop_t, Lind, theta_prime, phi_prime, R11..R33, magU, gradYlm_*; run link_update_gains after changing r/theta/phi/...']), ...
        xTxL, cy-200, fs, false);
    note(mdl, ['Tx waveform: ' waveDesc ' (variable Vtx_wave; series R Rsrc_Tx = ' srcType ', Vsrc_scale = ' num2str(Vscale) ')'], xTxL, cy-130, fs, false);
    note(mdl, sprintf(['Mode ports: eta0 to ground, reference grounded (Tx senses voltage, Rx injects current).\n' ...
        'Feed / internal ports: Z0 across the port (source on Tx Feed), reference floating (gmin to ground).']), ...
        xTxL, cy-70, fs, false);
    if o.TxOnly
        note(mdl, 'TX-ONLY model (radiated-field comparison):  v_Tx modes -> delay prop_t = (r - a_Tx)/c0 -> Gtx -> Rrot -> Vtheta, Vphi  (= r*E_theta, r*E_phi)', ...
            xmux, cy-130, fs+2, true);
    else
        note(mdl, 'Signal chain:  v_Tx modes -> delay prop_t -> Gtx -> Rrot -> Vtheta, Vphi -> 1/Lind -> integral -> Grx -> i_Rx modes', ...
            xmux, cy-130, fs+2, true);
    end
end
label_terminals(txParent, sTx, tx, fs);
if ~o.TxOnly, label_terminals(rxParent, sRx, rx, fs); end
blks = find_system(mdl, 'SearchDepth', 1, 'Type', 'Block');
for i = 1:numel(blks)
    try, set_param(blks{i}, 'FontSize', fs); catch, end %#ok<NOCOM>
end
lns = find_system(mdl, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line');
for i = 1:numel(lns)
    try, set_param(lns(i), 'FontSize', fs); catch, end %#ok<NOCOM>
end

%% ---------------------------------------------- save model + variables
matfile = fullfile(o.OutputDir, [mdl '_vars.mat']);
wsvars = evalin('base', 'who');
keepNames = {'eps0','eta0','Z0','Rgmin','c0','mu0','r','theta','phi','rot_angle','Ux','Uy','Uz','magU', ...
    'a_Tx','a_Rx','f0','Vamp','link_cfg','prop_t','Lind','Gtx_prime','Rrot','Grx','theta_prime','phi_prime', ...
    'ratTx','ratRx','Rsrc_Tx','Rfeed_Rx','Vsrc_scale','Rload_Tx','Rload_Rx','Vtx_wave','t_exc','SW_Ron','SW_Roff'};
keep = wsvars(ismember(wsvars, keepNames) | startsWith(wsvars, 'gradYlm_') | startsWith(wsvars, 'Ftheta') | startsWith(wsvars, 'Fphi') | ...
    ~cellfun('isempty', regexp(wsvars, '^R[123][123]$', 'once')));
for i = 1:numel(keep)
    eval([keep{i} ' = evalin(''base'', ''' keep{i} ''');']); %#ok<EVLDIR>
end
save(matfile, keep{:});
% PreLoadFcn: load <model>_vars.mat from the folder the model file is in (no absolute path stored, so the
% model and its .mat can be moved together); which() is the fallback if FileName is not yet set at preload
matname = [mdl '_vars.mat'];
set_param(mdl, 'PreLoadFcn', sprintf(['try, load(fullfile(fileparts(get_param(bdroot, ''FileName'')), ''%s'')); ' ...
    'catch, load(fullfile(fileparts(which(''%s'')), ''%s'')); end'], matname, mdl, matname));
wShadow = warning('off', 'Simulink:Engine:MdlFileShadowing');
save_system(mdl, fullfile(o.OutputDir, mdl));
warning(wShadow);
fprintf('Built %s (saved in %s). Variables stored in %s.\n', mdl, o.OutputDir, matfile);
fprintf('Run:  out = sim(''%s'');   Edit workspace vars (+ link_update_gains) to change r, theta, ...\n', mdl);
end

%% ===================================================================
%% ========================== local functions =========================
%% ===================================================================

function y = xt(ytop, i, pitch, fs)
y = ytop + (i-1)*pitch + (pitch-1.4*fs)/2;
end

function n = tx_file_name(fp)
[~, s, e] = fileparts(fp);  n = [s e];
end

function p = abs_path(p)
% absolute form of a file / folder path (relative paths are taken from the current folder)
if isempty(p) || ~isempty(regexp(p, '^([A-Za-z]:[\\/]|[\\/]|~)', 'once')), return, end
p = fullfile(pwd, p);
end

function note(mdl, txt, x, y, fs, bold)
a = Simulink.Annotation(mdl, txt);
a.Position = [x y];
a.FontSize = fs;
if bold, a.FontWeight = 'bold'; end
end

function label_terminals(mdl, sp, info, fs)
% Touchstone-header port names next to every +/- terminal of the S-parameter block
for k = 1:numel(info.names)
    s = sp.s(k);
    for t = 1:2
        if t == 1, pos = sp.sig(k, :);  sgn = '+'; else, pos = sp.ref(k, :);  sgn = '-'; end
        txt = sprintf('%d: %s (%s)', k, info.names{k}, sgn);
        w = 0.62*fs*numel(txt);
        if s > 0, x = pos(1) + 8;  else, x = pos(1) - 8 - w; end
        note(mdl, txt, x, pos(2) - 2.4*fs, fs, false);
    end
end
end

function [wave, srcType, desc, feedZ] = get_waveform(o, Z0val, fmax)
% Tx waveform timeseries (volts) from a CST export, the workspace, or [] for the sine source.
wave = [];  srcType = 'generator';  feedZ = NaN;
switch o.Waveform
    case 'sine'
        desc = sprintf('sine, Vamp = %g V at f0 = %.4g GHz', o.SourceAmplitude, o.SourceFreq/1e9);
    case 'cst'
        assert(~isempty(o.CSTFile) && isfile(o.CSTFile), 'CSTFile not found: %s', char(o.CSTFile));
        [t, v, info] = importCSTSignal(char(o.CSTFile), Z0val);   % the file's own Z0 (from #Parameters) wins
        if isnumeric(o.CSTSignal)
            k = o.CSTSignal;
        else
            k = find(strcmpi(info.names, char(o.CSTSignal)), 1);
            assert(~isempty(k), 'Signal "%s" not in %s (available: %s)', char(o.CSTSignal), ...
                char(o.CSTFile), strjoin(info.names, ', '));
        end
        zNote = '';
        if strcmp(info.sourceType{k}, 'generator')
            srcType = 'generator';
            feedZ = info.Z0(k);                       % CST port impedance: source resistance and Rx termination
            zNote = sprintf(', V+ = a*sqrt(%g), Rs = %g Ohm, Voc = 2*V+', feedZ, feedZ);
        else
            srcType = 'voltage';
        end
        wave = make_ts(t{k}, v{k});
        [~, fnm, fext] = fileparts(char(o.CSTFile));
        desc = sprintf('CST %s%s, signal %s (%s source%s)', fnm, fext, info.names{k}, srcType, zNote);
    case 'workspace'
        srcType = lower(char(o.SourceType));
        assert(ismember(srcType, {'generator', 'voltage'}), 'SourceType must be ''generator'' or ''voltage''.');
        nm = char(o.WaveformVar);
        assert(evalin('base', sprintf('exist(''%s'', ''var'')', nm)) == 1, ...
            'Workspace variable "%s" with the Tx waveform does not exist.', nm);
        w = evalin('base', nm);
        if isa(w, 'timeseries')
            wave = make_ts(w.Time, squeeze(w.Data));
        elseif isstruct(w)
            wave = make_ts(w.t, w.v);
            if isfield(w, 'type') && ismember(lower(char(w.type)), {'generator', 'voltage'})
                srcType = lower(char(w.type));
            end
        elseif isa(w, 'function_handle')
            dur = o.WaveformDuration;  dt = o.WaveformDt;
            assert(~isempty(dur), 'A function-handle waveform needs ''WaveformDuration''.');
            if isempty(dt), dt = 1/(10*fmax); end
            t = (0:dt:dur).';
            wave = make_ts(t, w(t));
        elseif isnumeric(w) && size(w, 2) == 2
            wave = make_ts(w(:, 1), w(:, 2));
        else
            error('Waveform variable "%s" must be a timeseries, N x 2 [t v], struct(t, v) or @(t) handle.', nm);
        end
        desc = sprintf('workspace variable %s (%s source)', nm, srcType);
    otherwise
        error('Waveform must be ''sine'', ''cst'' or ''workspace''.');
end
end

function ts = make_ts(t, v)
t = t(:);  v = v(:);
[t, iu] = unique(t);  v = v(iu);
if t(1) > 0, t = [0; t];  v = [0; v]; end
ts = timeseries(v, t, 'Name', 'Vtx_wave');
end

function [gtx, grx] = gain_expressions(txm, rxm)
% Gain-block expressions that reference the ADS-style gradYlm_* variables by
% name (the same formulas as the ADS VCVS/CCCS G fields).
vn = @(prefix, md) sprintf('%s_l%d_%s', prefix, md.l, ternary(md.m < 0, sprintf('mneg%d', abs(md.m)), sprintf('m%d', md.m)));
r1 = cell(1, numel(txm));  r2 = r1;
for i = 1:numel(txm)
    md = txm(i);
    if strcmp(md.type, 'TM')
        r1{i} = ['-' vn('gradYlm_theta_prime', md)];  r2{i} = ['-' vn('gradYlm_phi_prime', md)];
    else
        r1{i} = ['-' vn('gradYlm_phi_prime', md)];    r2{i} = vn('gradYlm_theta_prime', md);
    end
end
gtx = sprintf('[%s; %s]', strjoin(r1, ', '), strjoin(r2, ', '));
parts = cell(1, numel(rxm));
for i = 1:numel(rxm)
    md = rxm(i);
    par = ternary(mod(md.l, 2) == 0, '-1', '1');
    if strcmp(md.type, 'TM')
        parts{i} = sprintf('%s*4*pi*(-%s), %s*4*pi*(-%s)', par, vn('gradYlm_theta', md), par, vn('gradYlm_phi', md));
    else   % TE: extra polarity flip (ADS reverses the CCCS output terminals)
        parts{i} = sprintf('-%s*4*pi*(-%s), -%s*4*pi*(%s)', par, vn('gradYlm_phi', md), par, vn('gradYlm_theta', md));
    end
end
grx = ['[' strjoin(parts, '; ') ']'];
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function S = build_side(mdl, side, info, sp, opt)
% Circuitry around one S-parameter block. The block puts the + and - terminal
% of each port on the SAME edge (ports 1..N/2 left, the rest right), so every
% port's support circuit sits just outside its own edge and nothing wraps
% around the block.  The mode voltage / current pins are returned in S.ports (wired by the caller).
S = struct('mode_names', {{}}, 'gnd_for_solver', [], 'ports', {{}}, 'sgn', []);
n = numel(info.names);
fs = opt.fs;  TAGW = opt.TAGW;
if opt.rf, S = build_side_rf(mdl, side, info, sp, opt); return, end
for k = 1:n
    nm = info.names{k};
    y1 = sp.sig(k, 2);  y2 = sp.ref(k, 2);  xe = sp.sig(k, 1);
    s = sp.s(k);  if s == 0, s = 1; end
    if s > 0, ori = 'right'; flip = 'left'; else, ori = 'left'; flip = 'right'; end
    X = @(d, w) bx(xe, s, d, w);
    sigH = sp.sigH(k);  refH = sp.refH(k);
    kind = info.kind{k};
    if strcmp(kind, 'feed'), net = ['Feed' side]; else, net = [side '_' nm]; end
    refnet = [net '_ref'];
    hh = 1.4*fs;                                    % height of the small converter / tag blocks

    if strcmp(kind, 'mode')
        S.mode_names{end+1} = nm;
        xr1 = X(40, 40);
        Rb = addb(mdl, ['R_' net], 'fl_lib/Electrical/Electrical Elements/Resistor', [xr1(1) y1-15 xr1(2) y1+15], ...
            'R', 'eta0', 'Orientation', ori);
        ph = pH(Rb);  wire(mdl, sigH, ph.LConn(1), net);
        xg = X(120, 24);
        G1 = gnd(mdl, ['GND_R_' net], xg(1), y1-12, ori);
        wire(mdl, ph.RConn(1), pH(G1).LConn(1));
        xp = X(110, 40);  xt1 = X(190, TAGW);
        if strcmp(side, 'Tx')
            Sb = addb(mdl, ['Vsens_' net], 'fl_lib/Electrical/Electrical Sensors/Voltage Sensor', ...
                [xr1(1) y1+35 xr1(2) y1+75], 'Orientation', ori);
            ps = pH(Sb);  wire(mdl, sigH, ps.LConn(1), net);
            G2 = gnd(mdl, ['GND_V_' net], xg(1), y1+62, ori);
            wire(mdl, ps.RConn(2), pH(G2).LConn(1));
            Pb = addb(mdl, ['PS2S_' net], 'nesl_utility/PS-Simulink Converter', [xp(1) y1+40 xp(2) y1+40+hh], 'Orientation', ori);
            wire(mdl, ps.RConn(1), pH(Pb).LConn(1));
            S.ports{end+1} = pH(Pb).Outport(1);  S.sgn(end+1) = s;      % Tx: mode voltage leaves here (wired to the Mux)
        else
            Ib = addb(mdl, ['Isrc_' net], 'fl_lib/Electrical/Electrical Sources/Controlled Current Source', ...
                [xr1(1) y1+35 xr1(2) y1+75], 'Orientation', ori);
            pi_ = pH(Ib);  wire(mdl, sigH, pi_.LConn(1), net);   % + terminal injects into the mode node
            G2 = gnd(mdl, ['GND_I_' net], xg(1), y1+62, ori);
            wire(mdl, pi_.RConn(2), pH(G2).LConn(1));
            Sc = addb(mdl, ['S2PS_' net], 'nesl_utility/Simulink-PS Converter', [xp(1) y1+40 xp(2) y1+40+hh], 'Orientation', flip);
            wire(mdl, pH(Sc).RConn(1), pi_.RConn(1));
            S.ports{end+1} = pH(Sc).Inport(1);  S.sgn(end+1) = s;      % Rx: mode current enters here (wired from the Demux)
        end
        % reference of every mode port: straight to ground, just outside the block
        xgr = X(30, 24);
        Gr = gnd(mdl, ['GND_ref_' net], xgr(1), y2-12, ori);
        wire(mdl, refH, pH(Gr).LConn(1), refnet);
        if isempty(S.gnd_for_solver), S.gnd_for_solver = Gr; end
    else
        % ---- Feed / internal port: Z0 across the port, reference floating ----
        isTxFeed = strcmp(kind, 'feed') && strcmp(side, 'Tx');
        xx = X(40, 40);
        if isTxFeed
            rname = 'Rsrc_Tx';
        elseif strcmp(kind, 'feed')
            rname = 'Rfeed_Rx';
        else
            rname = ['Rload_' side];
        end
        if isTxFeed, bname = 'Term1';                       % ADS instance names
        elseif strcmp(kind, 'feed'), bname = 'Term2';
        else, bname = ['Term_' side '_' nm]; end
        Rb = addb(mdl, bname, 'fl_lib/Electrical/Electrical Elements/Resistor', [xx(1) y1-15 xx(2) y1+15], ...
            'R', rname, 'Orientation', ori);
        ph = pH(Rb);  wire(mdl, sigH, ph.LConn(1), net);
        if isTxFeed
            xv = X(130, 40);
            Vs = addb(mdl, 'VtSine1', 'fl_lib/Electrical/Electrical Sources/Controlled Voltage Source', ...
                [xv(1) y1-20 xv(2) y1+20], 'Orientation', ori);
            pv = pH(Vs);  wire(mdl, ph.RConn(1), pv.LConn(1), net);
            xc = X(200, 40);
            S2 = addb(mdl, 'S2PS_VtSine1', 'nesl_utility/Simulink-PS Converter', [xc(1) y1-1.4*fs-20 xc(2) y1-20], 'Orientation', flip);
            wire(mdl, pH(S2).RConn(1), pv.RConn(1));
            xg = X(270, 70);
            Gs = addb(mdl, 'Vsrc_gain', 'simulink/Math Operations/Gain', [xg(1) y1-1.6*fs-26 xg(2) y1-26], ...
                'Gain', 'Vsrc_scale', 'Orientation', flip);
            set_param(Gs, 'ShowName', 'off');
            add_line(mdl, pH(Gs).Outport(1), pH(S2).Inport(1));
            xz = X(370, 90);
            if opt.wave
                Sn = addb(mdl, 'VtSine1_wave', 'simulink/Sources/From Workspace', [xz(1) y1-1.6*fs-26 xz(2) y1-26], ...
                    'VariableName', 'Vtx_wave', 'Interpolate', 'on', 'OutputAfterFinalValue', 'Setting to zero', ...
                    'SampleTime', '0', 'Orientation', flip);
            else
                Sn = addb(mdl, 'VtSine1_sine', 'simulink/Sources/Sine Wave', [xz(1) y1-1.6*fs-26 xz(2) y1-26], ...
                    'Amplitude', 'Vamp', 'Frequency', '2*pi*f0', 'SampleTime', '0', 'Orientation', flip);
            end
            add_line(mdl, pH(Sn).Outport(1), pH(Gs).Inport(1));
            far = pv.RConn(2);
        else
            far = ph.RConn(1);
        end
        farP = get_param(far, 'Position');
        % return path to the - terminal (one row below), just outside the elements
        wire(mdl, far, refH, refnet, [farP(1)+s*30 farP(2); farP(1)+s*30 y2; xe+s*15 y2]);
        % differential port voltage V(+)-V(-): sensor between the two terminals of the port
        xf = X(40, 40);  xq = X(110, 40);  xo = X(190, 40);
        Sf = addb(mdl, ['Vsens_' net], 'fl_lib/Electrical/Electrical Sensors/Voltage Sensor', ...
            [xf(1) y1+35 xf(2) y1+75], 'Orientation', ori);
        pf = pH(Sf);  wire(mdl, sigH, pf.LConn(1));
        wire(mdl, pf.RConn(2), refH);
        Pf = addb(mdl, ['PS2S_' net], 'nesl_utility/PS-Simulink Converter', [xq(1) y1+40 xq(2) y1+40+hh], 'Orientation', ori);
        wire(mdl, pf.RConn(1), pH(Pf).LConn(1));
        Of = addb(mdl, net, 'simulink/Sinks/Out1', [xo(1) y1+40 xo(2) y1+40+hh], 'Orientation', ori);
        l = add_line(mdl, pH(Pf).Outport(1), pH(Of).Inport(1));  set_param(l, 'Name', net);
        % reference node: only the differential port voltage (above) is used; optional gmin gives the floating node a DC path
        a = X(60, 40);  g = X(140, 24);
        if opt.gmin
            Rg = addb(mdl, ['Rgmin_' refnet], 'fl_lib/Electrical/Electrical Elements/Resistor', [a(1) y2+15 a(2) y2+45], ...
                'R', 'Rgmin', 'Orientation', ori);
            wire(mdl, refH, pH(Rg).LConn(1));
            Gg = gnd(mdl, ['GND_gmin_' refnet], g(1), y2+17, ori);
            wire(mdl, pH(Rg).RConn(1), pH(Gg).LConn(1));
        end
    end
end
if isempty(S.gnd_for_solver)
    error('%s antenna has no mode port to provide the ground reference.', side);
end
end

function S = build_side_rf(mdl, side, info, sp, opt)
% RF Blockset version of build_side, same topology as the Simscape one:
%   mode ports: reference (-) terminal to an RF ground, eta0 to ground, sensed by an RF Outport (Tx, ideal voltage)
%     or driven by an RF Inport (Rx, ideal current);
%   Feed / internal ports: the termination (Z0, Rload_*, or the switch) and, for the Tx Feed, the source (RF Inport,
%     ideal voltage, in series with Rsrc_Tx) sit BETWEEN the port's + and - terminals; the - terminal floats (optional
%     gmin Rgmin to ground); the port voltage V(+)-V(-) is sensed differentially.
% Baseband, carrier 0 Hz: the signals are real time-domain voltages / currents.
S = struct('mode_names', {{}}, 'gnd_for_solver', [], 'ports', {{}}, 'sgn', []);
fs = opt.fs;  hh = 1.4*fs;
RFR = 'simrfV2elements/R';  RFG = 'simrfV2elements/Gnd';
for k = 1:numel(info.names)
    nm = info.names{k};
    y1 = sp.sig(k, 2);  y2 = sp.ref(k, 2);  xe = sp.sig(k, 1);
    s = sp.s(k);  if s == 0, s = 1; end
    if s > 0, ori = 'right'; flip = 'left'; else, ori = 'left'; flip = 'right'; end
    X = @(d, w) bx(xe, s, d, w);
    sigH = sp.sigH(k);  refH = sp.refH(k);
    kind = info.kind{k};
    if strcmp(kind, 'feed'), net = ['Feed' side]; else, net = [side '_' nm]; end
    isTxFeed = strcmp(kind, 'feed') && strcmp(side, 'Tx');
    isRxMode = strcmp(kind, 'mode') && strcmp(side, 'Rx');
    isMode   = strcmp(kind, 'mode');
    refnet   = [net '_ref'];
    if isMode
        % mode port: reference terminal straight to ground
        xgr = X(30, 24);
        Gr = addb(mdl, ['GND_ref_' net], RFG, [xgr(1) y2-12 xgr(1)+24 y2+12], 'Orientation', 'right');
        wire(mdl, refH, pH(Gr).LConn(1));
        if isempty(S.gnd_for_solver), S.gnd_for_solver = Gr; end
    elseif opt.gmin
        % Feed / internal port: floating reference, gmin to ground gives it a DC path
        a = X(60, 40);  g = X(140, 24);
        Rg = addb(mdl, ['Rgmin_' refnet], RFR, [a(1) y2+15 a(2) y2+45], 'Resistance', 'Rgmin', 'Orientation', ori);
        wire(mdl, refH, pH(Rg).LConn(1), refnet);
        Gg = addb(mdl, ['GND_gmin_' refnet], RFG, [g(1) y2+17 g(1)+24 y2+41], 'Orientation', ori);
        wire(mdl, pH(Rg).RConn(1), pH(Gg).LConn(1));
    end
    % termination resistor: to ground (mode port) or across the port, to its - terminal (Feed / internal port)
    if ~isTxFeed
        if strcmp(kind, 'mode'), rname = 'eta0';
        elseif strcmp(kind, 'feed'), rname = 'Rfeed_Rx';
        else, rname = ['Rload_' side]; end
        if strcmp(kind, 'feed'), bname = 'Term2';
        elseif strcmp(kind, 'internal'), bname = ['Term_' side '_' nm];
        else, bname = ['R_' net]; end
        xr = X(40, 40);
        if opt.sw && strcmp(kind, 'internal') && strcmp(side, 'Tx')
            % voltage-controlled switch: RF Variable Resistor, R = Ron*(Roff/Ron)^(1 - s), s = clamp(Vctrl, 0, 1) (SPICE switch,
            % VOFF = 0 V, VON = 1 V); the control is the ideal source Vs (From Workspace, scale 1)
            Rb = addb(mdl, bname, 'simrfV2elements/Variable Resistor', [xr(1) y1-15 xr(2) y1+15], 'Orientation', ori);
            xs1 = X(40, 70);  xs2 = X(140, 60);  xs3 = X(230, 90);  ys = y1 - 90;
            Sc = addb(mdl, 'VctrlSwitch', 'simulink/Sources/From Workspace', [xs3(1) ys xs3(2) ys+30], ...
                'VariableName', opt.ctrlvar, 'Interpolate', 'on', 'OutputAfterFinalValue', 'Setting to zero', ...
                'SampleTime', '0', 'Orientation', flip);
            Sa = addb(mdl, 'VctrlSat', 'simulink/Discontinuities/Saturation', [xs2(1) ys xs2(2) ys+30], ...
                'UpperLimit', '1', 'LowerLimit', '0', 'Orientation', flip);
            Sf = addb(mdl, 'Rswitch_map', 'simulink/User-Defined Functions/Fcn', [xs1(1) ys xs1(2) ys+30], ...
                'Expr', 'SW_Ron*(SW_Roff/SW_Ron)^(1-u)', 'Orientation', flip);
            add_line(mdl, pH(Sc).Outport(1), pH(Sa).Inport(1));
            add_line(mdl, pH(Sa).Outport(1), pH(Sf).Inport(1));
            add_line(mdl, pH(Sf).Outport(1), pH(Rb).Inport(1));
        else
            Rb = addb(mdl, bname, RFR, [xr(1) y1-15 xr(2) y1+15], 'Resistance', rname, 'Orientation', ori);
        end
        wire(mdl, sigH, pH(Rb).LConn(1), net);
        if isMode
            xg = X(120, 24);
            G1 = addb(mdl, ['GND_R_' net], RFG, [xg(1) y1-12 xg(1)+24 y1+12], 'Orientation', ori);
            wire(mdl, pH(Rb).RConn(1), pH(G1).LConn(1));
        else
            wire(mdl, pH(Rb).RConn(1), refH, refnet);           % across the port: + ... R ... -
        end
    end
    if isTxFeed
        xv = X(130, 60);
        % the Inport's own source impedance is ignored for an ideal-voltage source: series resistor Rsrc_Tx
        xr = X(40, 40);
        Rs = addb(mdl, 'Term1', RFR, [xr(1) y1-15 xr(2) y1+15], 'Resistance', 'Rsrc_Tx', 'Orientation', ori);
        wire(mdl, sigH, pH(Rs).LConn(1), net);
        Ib = addb(mdl, 'VtSine1', 'simrfV2util1/Inport', [xv(1) y1-20 xv(2) y1+20], 'Orientation', flip, ...
            'SimulinkInputSignalType', 'Ideal voltage', 'CarrierFreq', '0', 'ZS', '0', 'InternalGrounding', 'off');
        wire(mdl, pH(Ib).RConn(1), pH(Rs).RConn(1));            % source + to Rsrc_Tx
        wire(mdl, pH(Ib).RConn(2), refH, refnet);              % source - to the port's - terminal
        xg = X(270, 70);
        Gs = addb(mdl, 'Vsrc_gain', 'simulink/Math Operations/Gain', [xg(1) y1-1.6*fs-26 xg(2) y1-26], ...
            'Gain', 'Vsrc_scale', 'Orientation', flip);
        set_param(Gs, 'ShowName', 'off');
        xz = X(370, 90);
        if opt.wave
            Sn = addb(mdl, 'VtSine1_wave', 'simulink/Sources/From Workspace', [xz(1) y1-1.6*fs-26 xz(2) y1-26], ...
                'VariableName', 'Vtx_wave', 'Interpolate', 'on', 'OutputAfterFinalValue', 'Setting to zero', ...
                'SampleTime', '0', 'Orientation', flip);
        else
            Sn = addb(mdl, 'VtSine1_sine', 'simulink/Sources/Sine Wave', [xz(1) y1-1.6*fs-26 xz(2) y1-26], ...
                'Amplitude', 'Vamp', 'Frequency', '2*pi*f0', 'SampleTime', '0', 'Orientation', flip);
        end
        add_line(mdl, pH(Sn).Outport(1), pH(Gs).Inport(1));
        add_line(mdl, pH(Gs).Outport(1), pH(Ib).Inport(1));
    end
    if isRxMode
        xi = X(110, 40);
        Ic = addb(mdl, ['Isrc_' net], 'simrfV2util1/Inport', [xi(1) y1+35 xi(2) y1+75], 'Orientation', flip, ...
            'SimulinkInputSignalType', 'Ideal current', 'CarrierFreq', '0', 'ZS', 'Inf', 'InternalGrounding', 'on');
        wire(mdl, pH(Ic).RConn(1), sigH, net);
        % the RF ideal-current Inport draws a positive input current OUT of its node (the Simscape Controlled
        % Current Source pushes it in): negate so the mode current is injected into the mode node, as on Simscape
        xn = X(190, 40);
        Ng = addb(mdl, ['Ineg_' net], 'simulink/Math Operations/Gain', [xn(1) y1+40 xn(2) y1+40+hh], ...
            'Gain', '-1', 'Orientation', flip);
        add_line(mdl, pH(Ng).Outport(1), pH(Ic).Inport(1));
        S.mode_names{end+1} = nm;
        S.ports{end+1} = pH(Ng).Inport(1);  S.sgn(end+1) = s;
    else
        % voltage sensor: mode port to ground; Feed / internal port differential, V(+) - V(-)
        xq = X(110, 40);
        Ob = addb(mdl, ['Vsens_' net], 'simrfV2util1/Outport', [xq(1) y1+35 xq(2) y1+75], 'Orientation', ori, ...
            'SensorType', 'Ideal voltage', 'OutputFormat', 'Real Passband', 'CarrierFreq', '0', 'ZL', 'Inf', ...
            'InternalGrounding', ternary(isMode, 'on', 'off'));
        wire(mdl, sigH, pH(Ob).LConn(1), net);
        if ~isMode, wire(mdl, pH(Ob).LConn(2), refH, refnet); end
        if strcmp(kind, 'mode')
            S.mode_names{end+1} = nm;
            S.ports{end+1} = pH(Ob).Outport(1);  S.sgn(end+1) = s;
        else
            xo = X(190, 40);
            Of = addb(mdl, net, 'simulink/Sinks/Out1', [xo(1) y1+40 xo(2) y1+40+hh], 'Orientation', ori);
            l = add_line(mdl, pH(Ob).Outport(1), pH(Of).Inport(1));  set_param(l, 'Name', net);
        end
    end
end
if isempty(S.mode_names) && strcmp(side, 'Tx'), error('Tx antenna has no mode port.'); end
end

function add_rf_cfg(mdl, gblk, name, step)
b = addb(mdl, name, 'simrfV2util1/Configuration', [0 0 60 40], 'StepSize', num2str(step, '%.6g'), 'StepSize_unit', 's', ...
    'AutoFreq', 'on', 'AddNoise', 'off');
g = get_param(gblk, 'Position');
set_param(b, 'Position', [g(1) g(2)-70 g(1)+60 g(2)-30]);
wire(mdl, pH(b).LConn(1), pH(gblk).LConn(1));
end

function sp = add_sparam(mdl, name, x, y, w, row, nports, ratname, zref, rfo)
% S-parameter block: Simscape "S-Parameter Rational Fit" (rational object from the workspace), or, when rfo is a
% struct, the RF Blockset "S-parameters" block (data file, own fit / frequency-domain model), ungrounded terminals.
if nargin < 10, rfo = []; end
if isempty(rfo)
    blk = addb(mdl, name, 'ee_lib/Passive/S-Parameter Rational Fit', [x y x+w y+100], ...
        'numPorts', num2str(nports), 'rfmodelRationalObj', ratname, 'referenceImpedance', num2str(zref));
    ph = pH(blk);
    T = [ph.LConn(:); ph.RConn(:)];                 % terminal t: p_k = 2k-1, n_k = 2k
    assert(numel(T) == 2*nports, 'Unexpected terminal count on the S-parameter block.');
    sigT = T(1:2:end);  refT = T(2:2:end);
    nL = numel(ph.LConn);
else
    blk = addb(mdl, name, 'simrfV2elements/S-parameters', [x y x+w y+100]);
    if isempty(rfo.ratvar)
        set_param(blk, 'DataSource', 'Data file', 'File', rfo.file, 'InternalGrounding', 'off', ...
            'SparamRepresentation', rfo.mode, 'FitOpt', rfo.fitopt, 'FitTol', num2str(rfo.tol), ...
            'MaxPoles', num2str(rfo.maxpoles));
    elseif evalin('base', sprintf('isa(%s, ''rfmodel.rational'')', rfo.ratvar))
        % rationalfit output (N x N rfmodel.rational): the block takes N x N cells of poles, residues and direct terms
        R = evalin('base', rfo.ratvar);
        assert(isequal(size(R), [nports nports]), 'RFRationalVar %s is %s, expected %d x %d.', rfo.ratvar, mat2str(size(R)), nports, nports);
        assert(all(arrayfun(@(q) q.Delay == 0, R(:))), 'RFRationalVar %s: rfmodel.rational with Delay ~= 0 is not supported.', rfo.ratvar);
        Pc = arrayfun(@(q) q.A(:).', R, 'uni', 0);  Rc = arrayfun(@(q) q.C(:).', R, 'uni', 0);  Dc = arrayfun(@(q) q.D, R, 'uni', 0);
        e0 = cellfun(@isempty, Pc);  Pc(e0) = {-2*pi*1e12};  Rc(e0) = {0};   % entries without poles: one pole, zero residue
        vP = [rfo.ratvar '_poles'];  vR = [rfo.ratvar '_residues'];  vD = [rfo.ratvar '_direct'];
        assignin('base', vP, Pc);  assignin('base', vR, Rc);  assignin('base', vD, Dc);
        set_param(blk, 'DataSource', 'Rational model', 'isRationalObj', 'off', 'Poles', vP, 'Residues', vR, 'DF', vD, ...
            'InternalGrounding', 'off', 'SparamRepresentation', rfo.mode);
        fprintf('RF S-parameters block: rationalfit model %s (%d poles per entry, max).\n', rfo.ratvar, max(cellfun(@numel, Pc(:))));
    elseif evalin('base', sprintf('isa(%s, ''rational'') && iscell(%s.Poles)', rfo.ratvar, rfo.ratvar))
        % rational object with PoleSharing 'Column' / 'Individual' (N x N cells of poles and residues): the block's
        % RationalObject field takes shared-pole objects only, so it gets the N x N cells instead
        R = evalin('base', rfo.ratvar);
        Pc = cellfun(@(v) v(:).', R.Poles, 'uni', 0);  Rc = cellfun(@(v) v(:).', R.Residues, 'uni', 0);  Dc = num2cell(R.DirectTerm);
        % entries without poles (all-zero data, e.g. with PoleSharing 'Individual'): the block's network does not compile with
        % empty entries, so they get one fast stable pole with zero residue (contributes nothing)
        e0 = cellfun(@isempty, Pc);  Pc(e0) = {-2*pi*1e12};  Rc(e0) = {0};
        vP = [rfo.ratvar '_poles'];  vR = [rfo.ratvar '_residues'];  vD = [rfo.ratvar '_direct'];
        assignin('base', vP, Pc);  assignin('base', vR, Rc);  assignin('base', vD, Dc);
        set_param(blk, 'DataSource', 'Rational model', 'isRationalObj', 'off', 'Poles', vP, 'Residues', vR, 'DF', vD, ...
            'InternalGrounding', 'off', 'SparamRepresentation', rfo.mode);
        fprintf('RF S-parameters block: rational model %s, non-shared poles (%d per entry, max).\n', rfo.ratvar, max(cellfun(@numel, Pc(:))));
    else
        set_param(blk, 'DataSource', 'Rational model', 'isRationalObj', 'on', 'RationalObject', rfo.ratvar, 'InternalGrounding', 'off', ...
            'SparamRepresentation', rfo.mode);
    end
    ph = pH(blk);
    % terminal order from the mask's port labels ('k+' / 'k-' on LConn / RConn)
    md = get_param(blk, 'MaskDisplay');
    tk = regexp(md, 'port_label\(''(LConn|RConn)'',\s*(\d+),\s*''(\d+)([+-])''\)', 'tokens');
    sigT = zeros(nports, 1);  refT = zeros(nports, 1);
    for q = 1:numel(tk)
        t = tk{q};
        if strcmp(t{1}, 'LConn'), hh = ph.LConn(str2double(t{2})); else, hh = ph.RConn(str2double(t{2})); end
        k = str2double(t{3});
        if strcmp(t{4}, '+'), sigT(k) = hh; else, refT(k) = hh; end
    end
    assert(all(sigT ~= 0) && all(refT ~= 0), 'Could not map the RF Blockset S-parameters block terminals.');
    nL = numel(ph.LConn);
end
h = nL * row;
set_param(blk, 'Position', [x y x+w y+h]);
sp.blk = blk;  sp.top = y;  sp.height = h;
sp.sigH = sigT;  sp.refH = refT;
sp.sig = zeros(nports, 2);  sp.ref = zeros(nports, 2);  sp.s = zeros(nports, 1);
for k = 1:nports
    sp.sig(k, :) = get_param(sp.sigH(k), 'Position');
    sp.ref(k, :) = get_param(sp.refH(k), 'Position');
    sp.s(k) = sign(sp.sig(k, 1) - (x + w/2));
end
% RF rational model referenced to Z0 ~= 50 ohm: the block's rational model is referenced to 50 ohm, so each port gets an
% ideal transformer (RF Blockset 'Ideal Transformer', turns ratio n = sqrt(Z0/50), block on the secondary side): a network
% whose 50-ohm S-parameters are S, seen through it, has the S-parameters S referenced to Z0 (agrees with the block's own
% fit of the Touchstone file). Each block-side loop is isolated by the transformer, so its reference is tied to ground
% through a large resistor (no current flows in it).
if ~isempty(rfo) && ~isempty(rfo.ratvar) && isfield(rfo, 'ratz0') && abs(rfo.ratz0 - 50) > 1e-9
    nStr = sprintf('sqrt(%.10g/50)', rfo.ratz0);
    for k = 1:nports
        p = sp.sig(k, :);  s = sp.s(k);
        tname = sprintf('%s_Z0xfmr%d', name, k);
        if s < 0, pos = [p(1)-95 p(2)-13 p(1)-50 p(2)+38]; ori = 'right'; else, pos = [p(1)+50 p(2)-13 p(1)+95 p(2)+38]; ori = 'left'; end
        tb = addb(mdl, tname, 'simrfV2elements/Ideal Transformer', pos, 'n', nStr, 'Orientation', ori);
        th = pH(tb);
        wire(mdl, th.RConn(1), sp.sigH(k));  wire(mdl, th.RConn(2), sp.refH(k));       % secondary: block side
        rg = addb(mdl, sprintf('Rgmin_%s_xf%d', name, k), 'simrfV2elements/R', [pos(1) pos(4)+15 pos(3) pos(4)+35], ...
            'Resistance', '1e12', 'Orientation', 'down');
        gg = addb(mdl, sprintf('GND_%s_xf%d', name, k), 'simrfV2elements/Gnd', [pos(1)+10 pos(4)+50 pos(1)+34 pos(4)+74]);
        wire(mdl, pH(rg).LConn(1), sp.refH(k));  wire(mdl, pH(rg).RConn(1), pH(gg).LConn(1));
        sp.sigH(k) = th.LConn(1);  sp.refH(k) = th.LConn(2);                             % primary: the network's ports
        sp.sig(k, :) = get_param(sp.sigH(k), 'Position');  sp.ref(k, :) = get_param(sp.refH(k), 'Position');
    end
    fprintf('RF rational model referenced to %g ohm: ideal transformers (n = %s) on all %d ports.\n', rfo.ratz0, nStr, nports);
end
end

function add_solver_cfg(mdl, gblk, name)
b = addb(mdl, name, 'nesl_utility/Solver Configuration', [0 0 60 40]);
g = get_param(gblk, 'Position');
set_param(b, 'Position', [g(1) g(2)-70 g(1)+60 g(2)-30]);
wire(mdl, pH(b).RConn(1), pH(gblk).LConn(1));
end

function ss = wrap_antenna(mdl, spname, newname, fs, info)
% put the S-parameter block <mdl>/<spname> and its port circuits (mode-port eta0 resistors and grounds, mode-port
% sensors, gmin resistors, solver / Configuration block) into one subsystem <newname>.  Outside stay: the feed
% generator (source, its scaling and series resistance Term1), the Feed / internal-port terminations (Term*, the
% load switch and its control chain), the differential Feed / internal-port sensors (placed south of the block, their
% outputs pointing down), the workspace sources, the Mux / Demux of the signal chain and the root Outports.
% The Feed / internal ports leave the subsystem as RF terminals <port>+ / <port>- on its west side.
stopType = {'FromWorkspace', 'Outport', 'Inport', 'Mux', 'Demux', 'TransportDelay', 'Integrator'};
stopName = {'VCVS_TxGain_inst', 'VectorFieldRotate_inst', 'CCCS_RxAperture_inst', 'L1_L2_1overLind'};
isOutside = @(h) any(strcmp(get_param(h, 'BlockType'), stopType)) || any(strcmp(get_param(h, 'Name'), stopName)) || ...
    ~isempty(regexp(get_param(h, 'Name'), '^(Term|VtSine1|Vsrc_gain|Vctrl|Rswitch_map)', 'once'));
sp0 = get_param([mdl '/' spname], 'Handle');
keep = sp0;  todo = sp0;
while ~isempty(todo)
    b = todo(1);  todo(1) = [];
    pc = get_param(b, 'PortConnectivity');
    nb = [];
    for j = 1:numel(pc), nb = [nb, pc(j).SrcBlock(:).', pc(j).DstBlock(:).']; end %#ok<AGROW>
    nb = unique(nb(nb > 0));
    for h = nb
        if any(keep == h) || ~strcmp(get_param(get_param(h, 'Parent'), 'Type'), 'block_diagram'), continue, end
        if isOutside(h), continue, end
        keep(end+1) = h;  todo(end+1) = h; %#ok<AGROW>
    end
end
pos = get_param(sp0, 'Position');
h0 = find_system(mdl, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'SubSystem');
Simulink.BlockDiagram.createSubsystem(keep);
h1 = find_system(mdl, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'SubSystem');
hs = setdiff(h1, h0);
assert(numel(hs) == 1, 'Could not create the subsystem for %s.', spname);
set_param(hs, 'Name', newname, 'ShowName', 'on', 'FontSize', fs);
ss = [mdl '/' newname];
inner = [ss '/' spname];

% ---- names: signal outputs after the sensor inside (Vsens_<net> -> <net>); RF terminals after the S-parameter port
inO = find_system(ss, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'Outport');
inI = find_system(ss, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'Inport');
pio = find_system(ss, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'PMIOPort');
allp = [inO(:); inI(:); pio(:)];
for i = 1:numel(allp), set_param(allp(i), 'Name', sprintf('tmp_port_%d', i)); end
for i = 1:numel(inO)
    pc = get_param(inO(i), 'PortConnectivity');  src = pc(1).SrcBlock;
    if ~isempty(src) && src > 0, set_unique_name(ss, inO(i), regexprep(get_param(src, 'Name'), '^Vsens_', '')); end
end
php = get_param(ss, 'PortHandles');
for i = 1:numel(inI)
    k = str2double(get_param(inI(i), 'Port'));  ln = get_param(php.Inport(k), 'Line');  nm = '';
    if ln > 0
        sb = get_param(ln, 'SrcBlockHandle');
        if strcmp(get_param(sb, 'BlockType'), 'FromWorkspace'), nm = get_param(sb, 'VariableName'); else, nm = get_param(sb, 'Name'); end
    end
    if ~isempty(nm), set_unique_name(ss, inI(i), nm); end
end
sph = get_param(inner, 'PortHandles');                   % S-parameter terminal labels 'k+' / 'k-' (mask), by handle
tk = regexp(get_param(inner, 'MaskDisplay'), 'port_label\(''(LConn|RConn)'',\s*(\d+),\s*''(\d+)([+-])''\)', 'tokens');
lab = containers.Map('KeyType', 'double', 'ValueType', 'any');
for j = 1:numel(tk)
    t = tk{j};
    if strcmp(t{1}, 'LConn'), hh = sph.LConn(str2double(t{2})); else, hh = sph.RConn(str2double(t{2})); end
    lab(hh) = {str2double(t{3}), t{4}};
end
% electrical nets inside the subsystem, found geometrically (wire segments sharing a point are one net; a
% terminal belongs to the net of the wire point at its position) - physical branch segments do not report ends
ils = find_system(ss, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line');
P = cell(numel(ils), 1);  for i = 1:numel(ils), P{i} = get_param(ils(i), 'Points'); end
comp = 1:numel(ils);
for i = 1:numel(ils)
    for j = i+1:numel(ils)
        if any(ismember(P{i}, P{j}, 'rows'))
            a = comp(j);  b = comp(i);  comp(comp == a) = b;
        end
    end
end
allT = [];
for bh = reshape(find_system(ss, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'block'), 1, [])
    ph_ = get_param(bh, 'PortHandles');  allT = [allT; ph_.LConn(:); ph_.RConn(:)]; %#ok<AGROW>
end
nets = {};
for c = unique(comp)
    pts = vertcat(P{comp == c});  e = [];
    for t_ = allT.'
        tp_ = get_param(t_, 'Position');  on_ = false;
        for k_ = find(comp == c)
            Q = P{k_};
            for s_ = 1:size(Q, 1) - 1          % terminal on (or within 1 px of) a segment of this wire
                a_ = Q(s_, :);  b_ = Q(s_ + 1, :);
                if tp_(1) >= min(a_(1), b_(1)) - 1 && tp_(1) <= max(a_(1), b_(1)) + 1 && ...
                        tp_(2) >= min(a_(2), b_(2)) - 1 && tp_(2) <= max(a_(2), b_(2)) + 1 && ...
                        abs((b_(1) - a_(1))*(tp_(2) - a_(2)) - (b_(2) - a_(2))*(tp_(1) - a_(1))) <= norm(b_ - a_) + 1
                    on_ = true;  break
                end
            end
            if on_, break, end
        end
        if ~on_, ph_ = get_param(t_, 'Line');  on_ = ph_ > 0 && any(ils(comp == c) == ph_); end
        if on_, e(end+1, 1) = t_; end %#ok<AGROW>
    end
    nets{end+1} = e; %#ok<AGROW>
end
pnames = cell(numel(pio), 1);  pkey = inf(numel(pio), 1);
for i = 1:numel(pio)
    pp = get_param(pio(i), 'PortHandles');  pp = [pp.LConn(:); pp.RConn(:)];
    ends = [];
    for c = 1:numel(nets), if any(ismember(pp, nets{c})), ends = [ends; nets{c}]; end, end %#ok<AGROW>
    ends = ends(ends > 0 & isKey(lab, num2cell(ends)));
    if isempty(ends), continue, end
    v = lab(ends(1));
    pkey(i) = 2*v{1} - strcmp(v{2}, '+');
    pnames{i} = v;
end
[~, ord] = sort(pkey);
for i = 1:numel(pio), set_param(pio(i), 'Side', 'Left'); end
for n = 1:numel(ord)
    i = ord(n);
    set_param(pio(i), 'Port', num2str(n));
    v = pnames{i};
    if ~isempty(v)
        if v{1} <= numel(info.names), nm = [info.names{v{1}} v{2}]; else, nm = sprintf('%d%s', v{1}, v{2}); end
        set_unique_name(ss, pio(i), nm);
    end
end

% ---- outputs to the hub (Mux of the Tx modes / Demux of the Rx modes) in the hub's pin order
php = get_param(ss, 'PortHandles');  hub = [];
for i = 1:numel(php.Outport)
    ln = get_param(php.Outport(i), 'Line');  if ln < 0, continue, end
    d = get_param(ln, 'DstBlockHandle');  d = d(d > 0);
    if ~isempty(d) && any(strcmp(get_param(d(1), 'BlockType'), {'Mux', 'Demux'})), hub = d(1); break, end
end
side = 'out';
if isempty(hub)
    for i = 1:numel(php.Inport)
        ln = get_param(php.Inport(i), 'Line');  if ln < 0, continue, end
        s_ = get_param(ln, 'SrcBlockHandle');
        if s_ > 0 && any(strcmp(get_param(s_, 'BlockType'), {'Mux', 'Demux'})), hub = s_; side = 'in'; break, end
    end
end
if ~isempty(hub)
    hph = get_param(hub, 'PortHandles');
    if strcmp(side, 'out'), blkT = 'Outport'; hp = php.Outport; pins = hph.Inport; else, blkT = 'Inport'; hp = php.Inport; pins = hph.Outport; end
    key = 1e6*ones(numel(hp), 1);  innerP = zeros(numel(hp), 1);
    pb = find_system(ss, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', blkT);
    for i = 1:numel(pb)
        k = str2double(get_param(pb(i), 'Port'));  innerP(k) = pb(i);
        ln = get_param(hp(k), 'Line');  if ln < 0, continue, end
        if strcmp(side, 'out'), pe = get_param(ln, 'DstPortHandle'); else, pe = get_param(ln, 'SrcPortHandle'); end
        j = find(ismember(pins, pe), 1);  if ~isempty(j), key(k) = j; end
    end
    [~, o2] = sort(key + (1:numel(key)).'/1e3);
    for n = 1:numel(o2), set_param(innerP(o2(n)), 'Port', num2str(n)); end
end

% ---- compact block next to the hub, at its pin pitch
php = get_param(ss, 'PortHandles');
nO = numel(php.Outport);  nI = numel(php.Inport);  W = 12*fs;
if ~isempty(hub)
    hpos = get_param(hub, 'Position');
    if strcmp(side, 'out'), nh = numel(hph.Inport); n = nO; else, nh = numel(hph.Outport); n = nI; end
    pitch = (hpos(4) - hpos(2))/max(nh, 1);  H = pitch*max(n, 2);  top = hpos(2);
    if strcmp(side, 'out'), x1 = hpos(1) - 260 - W; else, x1 = hpos(3) + 260; end
else
    H = max(120, 2.4*fs*max(nO, nI));  top = (pos(2) + pos(4))/2 - H/2;  x1 = pos(1);
end
set_param(ss, 'Position', round([x1 top x1 + W top + H]));
php = get_param(ss, 'PortHandles');
ppos = @(h) get_param(h, 'Position');

% ---- outside blocks: generator / terminations to the west of their terminals, differential sensors to the south
term = php.LConn;                                        % west terminals, in port order (Feed+, Feed-, load1+, ...)
moved = []; %#ok<NASGU>
nS = 0;  ybot = top + H;
for t = 1:2:numel(term)
    tp = ppos(term(t));  yP = tp(2);
    % blocks connected to this terminal pair (outside)
    nbr = [];
    pcs_ = get_param(ss, 'PortConnectivity');
    for q = t:min(t + 1, numel(term))
        k_ = find(strcmp({pcs_.Type}, sprintf('LConn%d', q)), 1);  if isempty(k_), continue, end
        nbr = [nbr; pcs_(k_).SrcBlock(:); pcs_(k_).DstBlock(:)]; %#ok<AGROW>
    end
    nbr = unique(nbr(nbr > 0 & nbr ~= get_param(ss, 'Handle')));
    for h = nbr(:).'
        nm = get_param(h, 'Name');
        if startsWith(nm, 'Vsens_')                      % differential sensor: south of the block, output down
            nS = nS + 1;  xs = x1 + 30 + (nS - 1)*130;
            set_param(h, 'Orientation', 'down', 'Position', round([xs, ybot + 70, xs + 40, ybot + 110]));
            moved(end+1) = h; %#ok<AGROW>
            php_ = get_param(h, 'PortHandles');  ln = get_param(php_.Outport(1), 'Line');
            if ln > 0
                d = get_param(ln, 'DstBlockHandle');
                for db = d(d > 0).'
                    if strcmp(get_param(db, 'BlockType'), 'Outport')
                        set_param(db, 'Orientation', 'down', 'Position', round([xs + 13, ybot + 190, xs + 27, ybot + 190 + 1.4*fs]));
                    end
                end
            end
        elseif startsWith(nm, 'Term')                    % series source resistance / termination, west of the terminal
            set_param(h, 'Orientation', 'right', 'Position', round([x1 - 150, yP - 15, x1 - 110, yP + 15]));
            moved(end+1) = h; %#ok<AGROW>
            % generator (RF Inport source) and the signal chains feeding it / the variable resistor
            chain_west(h, x1 - 150, yP);
        elseif strcmp(get_param(h, 'BlockType'), 'SubSystem') && contains(get_param(h, 'ReferenceBlock'), 'Inport')
            set_param(h, 'Orientation', 'right', 'Position', round([x1 - 270, yP - 20, x1 - 210, yP + 20]));
            chain_west(h, x1 - 270, yP);
        end
    end
end

% ---- measurement outputs (Feed / internal-port voltages, sensed inside the block) not going to the hub: root
%      Outports south of the block, pointing down
php = get_param(ss, 'PortHandles');
for i = 1:numel(php.Outport)
    ln = get_param(php.Outport(i), 'Line');  if ln < 0, continue, end
    d = get_param(ln, 'DstBlockHandle');  d = d(d > 0);
    for db = d(:).'
        if ~strcmp(get_param(db, 'BlockType'), 'Outport'), continue, end
        nS = nS + 1;  xs = x1 + 20 + (nS - 1)*110;
        set_param(db, 'Orientation', 'down', 'Position', round([xs, ybot + 90, xs + 14, ybot + 90 + 1.4*fs]));
    end
end
% ---- redraw every top-level wire touching the subsystem or a block placed here (autorouted)
lns = find_system(mdl, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line');
sh = get_param(ss, 'Handle');
for i = 1:numel(lns)
    if ~ishandle(lns(i)) || get_param(lns(i), 'LineParent') ~= -1, continue, end
    try
        src = get_param(lns(i), 'SrcPortHandle');  dst = get_param(lns(i), 'DstPortHandle');
    catch
        continue
    end
    ends = unique([src(:); dst(:); line_ends(lns(i))], 'stable');  ends = ends(ends > 0);
    if numel(ends) < 2, continue, end
    % redraw every top-level wire (blocks were moved; stale bend points would leave diagonal wires)
    nm_ = get_param(lns(i), 'Name');
    delete_line(lns(i));
    for q = 2:numel(ends)
        try
            l2 = add_line(mdl, ends(1), ends(q), 'autorouting', 'smart');
            if ~isempty(nm_), set_param(l2, 'Name', nm_); end
        catch
        end
    end
end

% ---- physical branch segments keep stale bend points: replace any diagonal route by an orthogonal one
lns = find_system(mdl, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line');
for i = 1:numel(lns)
    if ~ishandle(lns(i)), continue, end
    P = get_param(lns(i), 'Points');
    if size(P, 1) < 2, continue, end
    dP = diff(P, 1, 1);
    if all(dP(:, 1) == 0 | dP(:, 2) == 0), continue, end
    a = P(1, :);  b = P(end, :);
    try, set_param(lns(i), 'Points', [a; a(1) b(2); b]); catch, end
end
    function chain_west(h, xr, yc)
        % move the signal blocks feeding block h (Gain, Saturation, Fcn, From Workspace, RF Inport) to its west, in a row
        x = xr;  cur = h;
        for step = 1:5
            pc_ = get_param(cur, 'PortConnectivity');  src_ = [];
            for jj = 1:numel(pc_)
                if ~isempty(pc_(jj).SrcBlock) && pc_(jj).SrcBlock > 0 && isempty(pc_(jj).DstBlock), src_ = pc_(jj).SrcBlock; break, end
            end
            if isempty(src_), break, end          % no further Simulink signal source (e.g. an RF block's physical input)
            bp_ = get_param(src_, 'Position');  w_ = bp_(3) - bp_(1);  hh_ = bp_(4) - bp_(2);
            x = x - 50 - w_;
            set_param(src_, 'Orientation', 'right', 'Position', round([x, yc - 60 - hh_/2, x + w_, yc - 60 + hh_/2]));
            cur = src_;
        end
    end
end


function e = line_ends(ln)
% every port handle on a wire and its branches
e = [];  todo = ln;  seen = [];
while ~isempty(todo)
    l = todo(1);  todo(1) = [];
    if any(seen == l) || l <= 0, continue, end
    seen(end+1) = l; %#ok<AGROW>
    e = [e; reshape(get_param(l, 'SrcPortHandle'), [], 1); reshape(get_param(l, 'DstPortHandle'), [], 1)]; %#ok<AGROW>
    ch = get_param(l, 'LineChildren');  todo = [todo, reshape(ch, 1, [])]; %#ok<AGROW>
end
end

function set_unique_name(ss, h, nm)
nm = strtrim(regexprep(char(nm), '[\n\r/]', ' '));
if isempty(nm), return, end
base = nm;  n = 1;
while ~isempty(find_system(ss, 'SearchDepth', 1, 'Name', nm)), n = n + 1; nm = sprintf('%s_%d', base, n); end
set_param(h, 'Name', nm);
end

function wrap_sparam(mdl, name, newname, info)
% put the S-parameter block <mdl>/<name> into a subsystem <newname> at the same place and size (all terminal
% connections kept); the subsystem's terminals are named <port name>+ / <port name>- from the Touchstone header
blk = [mdl '/' name];
pos = get_param(blk, 'Position');
h0 = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'SubSystem');
Simulink.BlockDiagram.createSubsystem(get_param(blk, 'Handle'));
h1 = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'SubSystem');
ss = setdiff(h1, h0);
assert(numel(ss) == 1, 'Could not create the subsystem for %s.', blk);
set_param(ss{1}, 'Name', newname);
ss = [mdl '/' newname];
set_param(ss, 'Position', pos, 'ShowName', 'on');
% terminal labels of the S-parameter block ('k+' / 'k-' on LConn / RConn, from its mask), by port handle
inner = [ss '/' name];  ph = get_param(inner, 'PortHandles');
tk = regexp(get_param(inner, 'MaskDisplay'), 'port_label\(''(LConn|RConn)'',\s*(\d+),\s*''(\d+)([+-])''\)', 'tokens');
lab = containers.Map('KeyType', 'double', 'ValueType', 'char');
for j = 1:numel(tk)
    t = tk{j};  k = str2double(t{3});
    if strcmp(t{1}, 'LConn'), hh = ph.LConn(str2double(t{2})); else, hh = ph.RConn(str2double(t{2})); end
    if k <= numel(info.names), lab(hh) = [info.names{k} t{4}]; else, lab(hh) = [t{3} t{4}]; end
end
pio = find_system(ss, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'PMIOPort');   % handles (names change below)
for i = 1:numel(pio)                 % temporary names first (the final names may already be in use)
    set_param(pio(i), 'Name', sprintf('tmp_term_%d', i));
end
for i = 1:numel(pio)
    pp = get_param(pio(i), 'PortHandles');  pp = [pp.LConn(:); pp.RConn(:)];
    ln = get_param(pp(1), 'Line');
    if ln < 0, continue, end
    ends = [get_param(ln, 'SrcPortHandle'); get_param(ln, 'DstPortHandle')];
    ends = ends(ends > 0 & isKey(lab, num2cell(ends)));
    if ~isempty(ends), set_param(pio(i), 'Name', lab(ends(1))); end
end
end

function xr = bx(anchor, s, d, w)
% x-range of a block of width w, offset d from anchor in direction s
if s > 0, xr = [anchor + d, anchor + d + w]; else, xr = [anchor - d - w, anchor - d]; end
end

function b = addb(mdl, name, src, pos, varargin)
b = [mdl '/' name];
add_block(src, b, 'Position', pos, varargin{:});
if ~isempty(regexp(name, '^(R_|GND|PS2S|Vsens|S2PS|Rgmin|Isrc)', 'once'))
    set_param(b, 'ShowName', 'off');
end
end

function g = gnd(mdl, name, x, y, ori)
if nargin < 5, ori = 'right'; end
g = addb(mdl, name, 'fl_lib/Electrical/Electrical Elements/Electrical Reference', [x y x+24 y+24], 'Orientation', ori);
end

function h = pH(blk)
h = get_param(blk, 'PortHandles');
end

function l = wire(mdl, a, b, name, via)
% physical/signal line from port a to port b, optionally through waypoints
if nargin >= 5 && ~isempty(via)
    pa = get_param(a, 'Position');  pb = get_param(b, 'Position');
    l = add_line(mdl, [pa; via; pb]);
else
    l = add_line(mdl, a, b);
end
if nargin >= 4 && ~isempty(name), set_param(l, 'Name', name); end
end

function info = read_snp_info(fp)
% Header parse: port names, a, Z0, ka_factor; frequency range from the data.
info = struct('names', {{}}, 'kind', {{}}, 'modes', struct('type', {}, 'l', {}, 'm', {}), ...
    'a', NaN, 'Z0', NaN, 'ka_is_one', false, 'fmin', NaN, 'fmax', NaN, 'ant', '');
fid = fopen(fp, 'r');
assert(fid ~= -1, 'Could not open %s', fp);
c = onCleanup(@() fclose(fid));
pm = containers.Map('KeyType', 'double', 'ValueType', 'any');
ln = fgetl(fid);
while ischar(ln)
    t = strtrim(ln);
    if ~startsWith(t, '!')
        if ~isempty(t) && ~startsWith(t, '#'), break; end
        ln = fgetl(fid);  continue
    end
    tk = regexp(t, '^!\s*Port\[(\d+)\]\s+(\S+)', 'tokens', 'once');
    if ~isempty(tk), pm(str2double(tk{1})) = tk{2}; end
    tk = regexpi(t, 'a=([\d.eE+-]+)m.*?Z0=([\d.eE+-]+)ohms', 'tokens', 'once');
    if ~isempty(tk) && isnan(info.a), info.a = str2double(tk{1}); info.Z0 = str2double(tk{2}); end
    tk = regexpi(t, 'ka_factor\s*=\s*([^\s,]+)', 'tokens', 'once');
    if ~isempty(tk), info.ka_is_one = strcmp(strtrim(tk{1}), '1'); end
    ln = fgetl(fid);
end
assert(pm.Count > 0, 'No "! Port[N] Name" comments found in %s', fp);
npt = max(cell2mat(keys(pm)));
info.names = cell(1, npt);  info.kind = cell(1, npt);
for i = 1:npt
    if isKey(pm, i), info.names{i} = pm(i); else, info.names{i} = sprintf('Port%d', i); end
    tk = regexp(info.names{i}, '^(TM|TE)_L(\d+)_M(neg)?(\d+)$', 'tokens', 'once');
    if ~isempty(tk)
        info.kind{i} = 'mode';
        m = str2double(tk{4}); if ~isempty(tk{3}), m = -m; end
        info.modes(end+1) = struct('type', tk{1}, 'l', str2double(tk{2}), 'm', m); %#ok<AGROW>
    elseif strcmp(info.names{i}, 'Feed')
        info.kind{i} = 'feed';
    else
        info.kind{i} = 'internal';
    end
end
sp = sparameters(fp);
f = sp.Frequencies(:);  f = f(f > 0);
info.fmin = f(1);  info.fmax = f(end);
[~, stem] = fileparts(fp);
tk = regexp(stem, 'S_dom_(?:psi|Q)_([A-Za-z0-9]+)_', 'tokens', 'once');
if isempty(tk), info.ant = stem; else, info.ant = tk{1}; end
end

function c = read_cst_params(fp, lenUnit)
% Z0 (Ohm) and a (m) from the '#Parameters = {...}' line of a CST ASCII export (NaN if absent)
c = struct('Z0', NaN, 'a', NaN);
fid = fopen(fp, 'r');  assert(fid ~= -1, 'Could not open %s', fp);
cl = onCleanup(@() fclose(fid));
for i = 1:40
    ln = fgetl(fid);
    if ~ischar(ln), break, end
    if startsWith(strtrim(ln), '#Parameters')
        tk = regexp(ln, '[{;]\s*Z0\s*=\s*([-+0-9.eE]+)', 'tokens', 'once');
        if ~isempty(tk), c.Z0 = str2double(tk{1}); end
        tk = regexp(ln, '[{;]\s*a\s*=\s*([-+0-9.eE]+)', 'tokens', 'once');
        if ~isempty(tk), c.a = str2double(tk{1}) * lenUnit; end
        break
    end
end
end

function rat = import_fit(src)
% rational object, or path to a .mat file holding one (first rational variable found)
if isa(src, 'rational'), rat = src; return, end
src = char(src);
assert(isfile(src), 'ImportFit: file not found: %s', src);
S = load(src);  fn = fieldnames(S);  rat = [];
for k = 1:numel(fn), if isa(S.(fn{k}), 'rational'), rat = S.(fn{k}); break, end, end
assert(~isempty(rat), 'ImportFit: no rational object in %s', src);
end

function [rat, zref] = fit_rational(fp, tol, maxPoles)
sp = sparameters(fp);
zref = sp.Impedance;
try
    rat = rational(sp, 'Tolerance', tol, 'MaxPoles', maxPoles);
catch
    rat = rational(sp);
end
end

function P = port_positions(h)
% n x 2 matrix of the [x y] positions of port handles (cell array or vector)
if iscell(h), h = [h{:}]; end
P = zeros(numel(h), 2);
for k = 1:numel(h), P(k, :) = get_param(h(k), 'Position'); end
end

function L = hub_layout(P, sg, hubdir, yTop, xEdge)
% Routing plan for wires between the mode ports (positions P, side sg = +1 right edge / -1 left edge of
% the S-parameter block) and a hub (Mux / Demux) standing in direction hubdir (+1 right, -1 left).
% Computed in a frame mirrored so the hub is on the right.  Ports on the hub-facing edge ('near') go
% straight across; ports on the other edge ('far') leave outward, rise to a channel above the block
% and cross over the top (rainbow nesting: the lowest port takes the outermost route), so nothing crosses.
sg = sg(:);
Pm = [hubdir*P(:, 1), P(:, 2)];
far = (sg*hubdir) < 0;
F = find(far);   [~, i] = sort(Pm(F, 2), 'descend');  F = F(i);     % lowest row first = outermost
N = find(~far);  [~, i] = sort(Pm(N, 2), 'ascend');   N = N(i);
nF = numel(F);  nN = numel(N);
if nN > 0, xN0 = max(Pm(N, 1)) + 50; else, xN0 = max(Pm(:, 1)) + 50; end
xF0 = xN0 + 40*nN + 40;
q = (1:nF).';
L = struct('Pm', Pm, 'F', F, 'N', N, 'nF', nF, 'nN', nN, 'order', [F; N], 'xN0', xN0);
if nF > 0
    L.xup = min(xEdge, min(Pm(F, 1))) - 50*(nF - q + 1);   % left of every block on that side
    L.yc  = yTop - 36*(nF - q);
    L.xdF = xF0 + 40*(nF - q);
else
    L.xup = [];  L.yc = [];  L.xdF = [];
end
L.xpin_min = xF0 + 40*nF + 60;
end

function W = hub_waypoints(L, pins, hubdir)
% polyline (port -> hub pin, actual coordinates) for every port; pins: n x 2 in the order L.order
pm = [hubdir*pins(:, 1), pins(:, 2)];
W = cell(size(L.Pm, 1), 1);
far = cell(L.nF, 1);
for q = 1:L.nF
    f = L.F(q);
    far{q} = [L.Pm(f, :); L.xup(q) L.Pm(f, 2); L.xup(q) L.yc(q); L.xdF(q) L.yc(q); L.xdF(q) pm(q, 2); pm(q, :)];
end
nN = L.nN;  best = cell(nN, 1);  bestc = inf;
if nN > 0
    if nN <= 7, pr = perms(1:nN); else, pr = 1:nN; end
    for ip = 1:size(pr, 1)
        xd = L.xN0 + 40*(pr(ip, :) - 1);
        polys = cell(nN, 1);
        for j = 1:nN
            nn = L.N(j);  k = L.nF + j;
            polys{j} = [L.Pm(nn, :); xd(j) L.Pm(nn, 2); xd(j) pm(k, 2); pm(k, :)];
        end
        c = count_crossings([far; polys]);
        if c < bestc, bestc = c;  best = polys;  end
        if c == 0, break, end
    end
else
    bestc = count_crossings(far);
end
if bestc > 0
    warning('build_tx_rx_link_simulink:WireCrossings', '%d wire crossing(s) in a mode-wire group.', bestc);
end
for q = 1:L.nF, W{L.F(q)} = far{q}; end
for j = 1:nN, W{L.N(j)} = best{j}; end
for i = 1:numel(W)
    if isempty(W{i}), continue, end
    pts = W{i};
    keep = [true; any(diff(pts) ~= 0, 2)];
    pts = pts(keep, :);
    pts(:, 1) = hubdir*pts(:, 1);
    W{i} = pts;
end
end

function c = count_crossings(polys)
c = 0;
for a = 1:numel(polys)
    for b = a+1:numel(polys)
        A = polys{a};  B = polys{b};
        for i = 1:size(A, 1)-1
            for j = 1:size(B, 1)-1
                if seg_inter(A(i, :), A(i+1, :), B(j, :), B(j+1, :)), c = c + 1; end
            end
        end
    end
end
end

function tf = seg_inter(a, b, c, d)
cr = @(u, v) u(1)*v(2) - u(2)*v(1);
o1 = cr(b - a, c - a);  o2 = cr(b - a, d - a);  o3 = cr(d - c, a - c);  o4 = cr(d - c, b - c);
tf = (o1*o2 < 0) && (o3*o4 < 0);
if ~tf
    on = @(p, q, r) min(p(1), q(1)) <= r(1) && r(1) <= max(p(1), q(1)) && min(p(2), q(2)) <= r(2) && r(2) <= max(p(2), q(2));
    tf = (o1 == 0 && on(a, b, c)) || (o2 == 0 && on(a, b, d)) || (o3 == 0 && on(c, d, a)) || (o4 == 0 && on(c, d, b));
end
end

function x = model_edge(mdl, hubdir)
% extreme block edge in the frame mirrored so the hub is on the right: the leftmost left edge of any block
% (hubdir = +1) or minus the rightmost right edge (hubdir = -1); the wires that go over the top stay outside it
b = find_system(mdl, 'SearchDepth', 1, 'Type', 'Block');
pos = zeros(numel(b), 4);
for k = 1:numel(b), pos(k, :) = get_param(b{k}, 'Position'); end
if hubdir > 0, x = min(pos(:, 1)); else, x = -max(pos(:, 3)); end
end
