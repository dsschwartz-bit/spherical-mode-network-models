function E = ook_load_excitation(source, opt)
%OOK_LOAD_EXCITATION  OOK inputs for the Simulink Tx-only models: generator incident voltage
%   Vinc (feed, port rf_port) and ideal control voltage Vs (switch, port ctrl_port), on one
%   common time grid.
%
%   E = ook_load_excitation("cst", opt)    from the CST external-port voltage / current exports
%       opt.sig_dir, opt.voltage_file, opt.current_file, opt.rf_port, opt.ctrl_port
%       a1 = (V1 + Z I1)/(2 sqrt(Z)), Vinc = a1 sqrt(Z) = (V1 + Z I1)/2, Vs = V2 (as
%       ook_cst_to_ads_mdf.m), Z = the export's Z0.
%   E = ook_load_excitation("user", opt)   from a user-excitation .mat (generate_user_excitations.m)
%       opt.user_file: .mat with timeseries Vtx_wave (= Vinc_1) and Vs_ctrl (= Vs_2) and the
%       struct exc_info (Z0, ...).
%
%   E : struct with t (s, column, starting at t0), Vinc, Vs (V), Z (ohm), t0 (s, start of the
%       grid: Simulink time = t - t0), desc (text for the console / plots)

source = lower(string(source));
switch source
    case "cst"
        [tV, vV, infoV] = importCSTSignal(char(fullfile(opt.sig_dir, opt.voltage_file)), []);
        [tI, vI, infoI] = importCSTSignal(char(fullfile(opt.sig_dir, opt.current_file)), []);
        kV = @(p) find(infoV.port == p, 1);   kI = @(p) find(infoI.port == p, 1);
        Z  = infoV.Z0(kV(opt.rf_port));
        t1 = tV{kV(opt.rf_port)};    V1 = vV{kV(opt.rf_port)};
        I1 = interp1(tI{kI(opt.rf_port)}, vI{kI(opt.rf_port)}, t1, 'linear', 0);
        t2 = tV{kV(opt.ctrl_port)};  V2 = vV{kV(opt.ctrl_port)};
        a1 = (V1 + Z*I1) / (2*sqrt(Z));
        dtc = min(median(diff(t1)), median(diff(t2)));
        t   = (min(t1(1), t2(1)) : dtc : max(t1(end), t2(end))).';
        Vinc = interp1(t1, a1*sqrt(Z), t, 'linear', 0);
        Vs   = interp1(t2, V2,         t, 'linear', 0);
        desc = "CST exports " + opt.voltage_file + " / " + opt.current_file;
    case "user"
        if ~isfile(opt.user_file), error('User excitation file not found: %s', opt.user_file); end
        S = load(opt.user_file, 'Vtx_wave', 'Vs_ctrl', 'exc_info');
        if ~isfield(S, 'Vtx_wave') || ~isfield(S, 'Vs_ctrl')
            error('%s must hold the timeseries Vtx_wave and Vs_ctrl (generate_user_excitations.m).', opt.user_file);
        end
        t1 = S.Vtx_wave.Time(:);  x1 = double(S.Vtx_wave.Data(:));
        t2 = S.Vs_ctrl.Time(:);   x2 = double(S.Vs_ctrl.Data(:));
        if isfield(S, 'exc_info') && isfield(S.exc_info, 'Z0'), Z = S.exc_info.Z0; else, Z = 50; end
        if isequal(t1, t2)
            t = t1;  Vinc = x1;  Vs = x2;
        else                                                 % common grid: the finer of the two
            dtc = min(median(diff(t1)), median(diff(t2)));
            t   = (min(t1(1), t2(1)) : dtc : max(t1(end), t2(end))).';
            Vinc = interp1(t1, x1, t, 'linear', 0);  Vs = interp1(t2, x2, t, 'linear', 0);
        end
        [~, fn, fe] = fileparts(opt.user_file);
        desc = "user excitation " + fn + fe;
    otherwise
        error('excitation_source must be "cst" or "user" (got "%s").', source);
end
E = struct('t', t, 'Vinc', Vinc, 'Vs', Vs, 'Z', Z, 't0', t(1), 'desc', desc);
fprintf('Excitation: %s\n  Z = %g ohm; %.4g-%.4g ns, %d points; Vinc peak %.4g V; Vs %.3g to %.3g V\n', ...
    desc, Z, t(1)*1e9, t(end)*1e9, numel(t), max(abs(Vinc)), min(Vs), max(Vs));
end
