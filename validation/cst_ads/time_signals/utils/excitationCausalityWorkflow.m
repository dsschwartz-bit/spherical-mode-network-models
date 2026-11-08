function out = excitationCausalityWorkflow(in_file, port, snp_file, varargin)
%EXCITATIONCAUSALITYWORKFLOW  Check an excitation against the ADS model's
%S-parameter block and, if needed, low-pass it and write CST / ADS files.
%
%   out = excitationCausalityWorkflow(in_file, port, snp_file, Name, Value, ...)
%
%   in_file  : CST time-signal export holding the incident wave i<port>
%   port     : port number (signal i<port>)
%   snp_file : Touchstone file of the S-parameter block used in ADS
%
%   1. The S-parameter block is driven by the excitation: every response
%      y_ij = IFFT(S_ij * X), with S known only up to its highest frequency
%      fmax, as the ADS transient sees it. The non-causal energy of a
%      response is its energy before the excitation starts (t < 0), as a
%      percentage of its own energy; only responses within SigdB of the
%      strongest count.
%   2. If the worst one exceeds EncTolPct, the excitation is low-passed
%      with a Gaussian, delayed by ShiftSigmas standard deviations of its
%      impulse response (H = exp(-(f/f3)^2 3ln10/20) exp(-j2pi f shift)):
%      to FirstTargetdB at fmax, or StepdB below the current level if it is
%      already lower. Step 1 is repeated, up to MaxIter filters.
%   3. Causality figures are shown when the input fails (PlotIfFail) and
%      the filtered excitation is written as a CST ASCII file (layout of a
%      CST time-signal export) and an ADS MDF file (Vinc_<port> / Vs_<port>,
%      conventions of cst_to_ads_mdf.m).
%
%   out: struct with fields passed (input already causal), encIn, encOut
%   (worst significant non-causal energy, %), levelIn / levelOut (spectrum
%   at fmax, dB re peak), f3, shift, sigma_t, cstFile, mdfFile, t, y (final
%   signal, CST unit), unit.
%
%   Options (defaults): EncTolPct 0.01, SigdB 20, FirstTargetdB -40,
%   StepdB 10, ShiftSigmas 5, MaxIter 6, DfCheck 5e6, WriteIfCausal false,
%   NormalizePeak false, TailSigmas 10, BwLeveldB 20, OutDir
%   <repo>/output/cst_ads/filtered_excitations, PrefixPort "Vinc_", PrefixVoltage
%   "Vs_", BlockName "BLK", Z0default [], PlotIfFail true, PlotAlways false,
%   Verbose true, ForceFilter false (apply at least one filter even when the
%   input meets the target; then the filter is a Gaussian band-pass, -3 dB at
%   the ForceBand edges, default [lowest S-parameter frequency, fmax]; with
%   ForceType "lowpass" the usual Gaussian low-pass, e.g. to smooth a step that must
%   keep its DC level), ForceType "bandpass", ForceBand [], CheckOnly false (only
%   the pass/fail check of the input: no filter, plots or files),
%   MinImprovePct 5 (stop filtering when the next filter improves the worst
%   non-causal energy by less than this, and keep the previous one), ExtraDelay 0 (s, exact extra delay added to the
%   output, e.g. to line up several ports' filter delays), FminCheck 0 (Hz:
%   only the spectral content above FminCheck of a response and of its part
%   before t = 0 is counted, e.g. to ignore a poorly known low end), Signal "" (signal name; default "i<port>", e.g. "v1 [1]" for a
%   port-voltage export, which is used as an ideal voltage source), LowBandWarnPct []
%   (warn when the non-causal energy from frequencies below the lowest S-parameter
%   point, DC excluded, exceeds this %; [] = EncTolPct/10), WriteMDF true (false:
%   write the CST ASCII file only).
%
%   Requires RF Toolbox (sparameters), importCSTSignal.m, generate_mdf.m.

ip = inputParser;
ip.addParameter('EncTolPct', 0.01);      ip.addParameter('SigdB', 20);
ip.addParameter('FirstTargetdB', -40);   ip.addParameter('StepdB', 10);
ip.addParameter('ShiftSigmas', 5);       ip.addParameter('MaxIter', 6);
ip.addParameter('DfCheck', 5e6);         ip.addParameter('WriteIfCausal', false);
ip.addParameter('NormalizePeak', false); ip.addParameter('TailSigmas', 10);
ip.addParameter('BwLeveldB', 20);
ip.addParameter('OutDir', fullfile(repo_root, 'output', 'cst_ads', 'filtered_excitations'));
ip.addParameter('PrefixPort', "Vinc_");  ip.addParameter('PrefixVoltage', "Vs_");
ip.addParameter('BlockName', "BLK");     ip.addParameter('Z0default', []);
ip.addParameter('PlotIfFail', true);     ip.addParameter('PlotAlways', false);
ip.addParameter('Verbose', true);        ip.addParameter('Signal', "");
ip.addParameter('ForceFilter', false);   % apply at least one Gaussian even if the input is already causal
ip.addParameter('ExtraDelay', 0);        % s, exact extra delay of the output (linear phase, e.g. to align ports)
ip.addParameter('WriteMDF', true);       % false: CST ASCII only (e.g. user-defined CST excitations; MDF made later from the CST export)
ip.addParameter('ForceType', "bandpass");  % forced filter on an input that already meets the target: "bandpass" or "lowpass"
ip.addParameter('ForceBand', []);        % Hz [flo fhi], band-pass edges of a forced filter; [] = [lowest S frequency, fmax]
ip.addParameter('MinImprovePct', 5);     % %, stop when a filter improves the worst non-causal energy by less than this
ip.addParameter('LowBandWarnPct', []);   % %, warn if the non-causal energy below the lowest S frequency exceeds this; [] = EncTolPct/10
ip.addParameter('FminCheck', 0);         % Hz, count non-causal energy only above this frequency (0 = all)
ip.addParameter('CheckOnly', false);     % true: only the causality check of the input (no filter, plots or files)
ip.parse(varargin{:});  o = ip.Results;
say = @(varargin) o.Verbose && fprintf(varargin{:});

%% ---- S-parameters on the check grid ------------------------------------------
[Su, fu, tc, toTime, nameS, fData] = snpOnGrid(snp_file, o.DfCheck);
fmax = fu(end);
say('Causality check against %s (0-%.4g GHz, %d non-zero entries)\n', snp_file, fmax/1e9, numel(Su));

%% ---- input excitation --------------------------------------------------------
[tcs, vcs, info] = importCSTSignal(in_file, o.Z0default);   % seconds, volts
if strlength(o.Signal) > 0, sigName = string(o.Signal); else, sigName = "i" + port; end
k = find(strcmp(info.names, sigName), 1);
if isempty(k), error('No signal %s in %s (signals: %s).', sigName, in_file, strjoin(info.names, ', ')); end
t0 = tcs{k}(:);  vV = vcs{k}(:);
isGen = strcmp(info.sourceType{k}, 'generator');   % port voltages and voltage sources -> VOLTAGE
if ~strcmp(sigName, "i" + port)                     % e.g. "v1 [cosim Tran1 1]": keep the tag in file names
    xtag = regexprep(strtrim(regexprep(char(sigName), '^\S+', '')), '[^\w]+', '_');
    xtag = regexprep(xtag, '^_|_$', '');
else
    xtag = '';
end
if isGen, scale = sqrt(info.Z0(k)); else, scale = 1; end
x0 = vV / scale;                                        % CST unit (W^(1/2) or V)
unitCST = string(info.signalUnit{k});
Lin = readlines(in_file);  parLine = Lin(find(startsWith(strtrim(Lin), "#Parameters"), 1));
[~, inName] = fileparts(in_file);
if ~isempty(xtag), inName = string(inName) + "_exc" + xtag; end
spec = @(tt, vv, ff) nufft(vv .* [diff(tt(1:2)); diff(tt(1:end-1)) + diff(tt(2:end)); diff(tt(end-1:end))]/2, tt, ff);
levelAt = @(X) 20*log10(abs(X(end)) / max(abs(X)));

%% ---- check, filter, re-check ------------------------------------------------------
dt = median(diff(t0));
tu = t0;  y = x0;  f3 = Inf;  fc = [];  shift = 0;  sigma_t = 0;  filtered = false;
X  = spec(t0 - t0(1), x0, fu);  X0 = X;           % input spectrum: every filter is applied to the input
dfc = fu(2) - fu(1);
Cin = checkCausal(Su, X, tc, toTime, o.SigdB, o.FminCheck, dfc, fData);  C = Cin;
lev = levelAt(X);  levIn = lev;
if o.CheckOnly                                   % check only: pass/fail of the unfiltered input
    out = struct('passed', Cin.worst <= o.EncTolPct, 'encIn', Cin.worst, 'encOut', Cin.worst, 'levelIn', lev, ...
        'levelOut', lev, 'f3', Inf, 'band', [], 'encLow', Cin.worstLow, 'shift', 0, 'delay', 0, 'sigma_t', 0, ...
        'cstFile', "", 'mdfFile', "", 't', t0, 'y', x0, 'unit', unitCST);
    return
end
say('  excitation %s (%s): spectrum at %.3g GHz %.1f dB re peak; worst non-causal energy %.4f %% (target %g %%%s)\n', ...
    inName, sigName, fmax/1e9, lev, C.worst, o.EncTolPct, ternary(o.FminCheck > 0, sprintf(', above %.4g MHz', o.FminCheck/1e6), ''));
if lev > o.FirstTargetdB, target = o.FirstTargetdB; else, target = lev - o.StepdB; end
nFilt = 0;
for it = 1:o.MaxIter
    meets = C.worst <= o.EncTolPct;
    if meets && ~(o.ForceFilter && it == 1), break, end
    if meets && strcmpi(o.ForceType, "bandpass")   % forced filter on an input that already meets the target:
        fb = o.ForceBand;                % Gaussian band-pass, -3 dB at both band edges
        if isempty(fb), fb = [fData fmax]; end
        fcN = fb;  f3N = NaN;            % band edges; sigma_t from the impulse-response width
    else                                 % Gaussian low-pass to the target level at fmax
        fcN = [];
        lvl = @(f) levelAt(X0 .* gaussLP(fu, f, 0));   % level of the filtered input (filters do not stack)
        if levIn > -150                  % spectrum at fmax resolvable: solve for the filtered level
            f3N = 10^fzero(@(lf) lvl(10.^lf) - target, log10([0.05 50]*fmax));
        else                             % at the numerical floor (e.g. a whole-cycle sine): the Gaussian
            f3N = fmax * sqrt(3 / (levIn - target));   % itself gives the extra attenuation at fmax (H dB = -3 (f/f3)^2)
        end
    end
    if isempty(fcN), sigN = sqrt(2*log(10)*3/20) / (2*pi*f3N);   % Gaussian low-pass: std of its impulse response
    else,            sigN = bpSigma(fcN); end                    % band-pass: equivalent std (rms width)
    shN = o.ShiftSigmas * sigN;
    [tuN, yN] = applyGauss(t0, x0, dt, f3N, shN + o.ExtraDelay, o.TailSigmas*sigN, fcN);
    XN = spec(tuN - t0(1), yN, fu);
    CN = checkCausal(Su, XN, tc, toTime, o.SigdB, o.FminCheck, dfc, fData);
    levN = levelAt(XN);
    if ~isempty(fcN), fdesc = sprintf('Gaussian band-pass -3 dB at %.4g and %.4g GHz', fcN/1e9);
    else,       fdesc = sprintf('Gaussian -3 dB at %.4g GHz', f3N/1e9); end
    say('    filter %d: %s, shift %.1f ps -> %.1f dB at %.3g GHz, worst non-causal energy %.4f %%\n', ...
        it, fdesc, shN*1e12, levN, fmax/1e9, CN.worst);
    if filtered && CN.worst > C.worst * (1 - o.MinImprovePct/100)
        say('      less than %g %% better than filter %d: keeping filter %d\n', o.MinImprovePct, nFilt, nFilt);
        break
    end
    tu = tuN;  y = yN;  X = XN;  C = CN;  lev = levN;  f3 = f3N;  fc = fcN;  sigma_t = sigN;  shift = shN;
    filtered = true;  nFilt = it;
    target = lev - o.StepdB;
end
passed = ~filtered;
if ~filtered && o.ExtraDelay > 0       % delay of the unfiltered input: shift its time axis (samples unchanged;
    tu = [t0(1); t0 + o.ExtraDelay];  y = [0; x0];   % a band-limited delay would ring at sharp edges), zero before it
    say('  delayed by %.1f ps (ExtraDelay)\n', o.ExtraDelay*1e12);
end
if filtered && C.worst > o.EncTolPct
    warning('Still %.3f %% non-causal energy after %d filter(s).', C.worst, nFilt);
end
if isempty(o.LowBandWarnPct), lowThr = o.EncTolPct/10; else, lowThr = o.LowBandWarnPct; end
say('  non-causal energy from 0 < f < %.4g MHz (below the lowest S-parameter point): %.4f %% (warning above %g %%)\n', fData/1e6, C.worstLow, lowThr);
if C.worstLow > lowThr
    warning(['Frequencies below the lowest S-parameter point (%.4g MHz, DC excluded), where %s has no data, ' ...
        'contribute %.4f %% non-causal energy (worst significant response; warning level %g %%).'], fData/1e6, snp_file, C.worstLow, lowThr);
end
if passed, say('  meets the %g %% causality target: no filter needed\n', o.EncTolPct); end
if o.NormalizePeak && filtered, y = y * max(abs(x0)) / max(abs(y)); end

out = struct('passed', passed, 'encIn', Cin.worst, 'encOut', C.worst, 'levelIn', levIn, 'levelOut', lev, ...
    'f3', f3, 'band', fc, 'encLow', C.worstLow, 'shift', shift, 'delay', shift*filtered + o.ExtraDelay, 'sigma_t', sigma_t, 'cstFile', "", 'mdfFile', "", 't', tu, 'y', y, 'unit', unitCST);

%% ---- plots -------------------------------------------------------------------------
if (o.PlotIfFail && ~passed) || o.PlotAlways
    plotCausality(inName, port, unitCST, t0, x0, tu, y, fu, Cin, C, nameS, o, filtered, f3, shift, spec, fc);
end

%% ---- write CST ASCII and ADS MDF ----------------------------------------------
if passed && ~o.WriteIfCausal, return, end
fs = linspace(0, min(1/(2*dt), 12e9), 4801).';
[flo, fhi] = bandEdges(fs, abs(spec(tu, y, fs)), o.BwLeveldB);
if ~isfolder(o.OutDir), mkdir(o.OutDir); end
if filtered && ~isempty(fc), tag = "_gBP" + ghzStr(fc(1)) + "_" + ghzStr(fc(2));
elseif filtered,       tag = "_gLP" + ghzStr(f3);
else,                  tag = ""; end
cstFile  = fullfile(o.OutDir, inName + tag + ".txt");
baseName = regexprep(inName, '(_[\d.]+[kMGT]?Hz)+$', '') + tag;
if isGen, excType = "PORT"; mdfName = o.PrefixPort + port;
else,     excType = "VOLTAGE"; mdfName = o.PrefixVoltage + port; end
Tdur    = tu(end) - tu(1);
if isGen, srcTag = "_PORT_" + string(num2str(info.Z0(k))) + "Ohms";   % generator with source impedance Z0
else,     srcTag = "_VOLTAGE"; end                                % ideal voltage source
mdfFile = fullfile(o.OutDir, strrep(baseName + srcTag, "-", "_") + ".mdf");   % ADS cannot import MDF names with hyphens

fid = fopen(cstFile, 'w');                              % CST time-signal layout, time in ns
if ~isempty(parLine), fprintf(fid, '%s\n', parLine); end
fprintf(fid, '#"Time / ns"\t"i%d [Real / %s]"\n', port, unitCST);
fprintf(fid, '#---------------------------------\n');
fprintf(fid, '%.13g\t%.13e\n', [tu*1e9, y].');
fclose(fid);

comments = "Excitation: " + string(inName) + " (CST signal " + sigName + ")";
if filtered
    if ~isempty(fc), fd = sprintf('Gaussian band-pass in f^2, -3 dB at %.4g and %.4g GHz (input already causal, ForceFilter)', fc/1e9);
    else,      fd = sprintf('Gaussian low-pass -3 dB at %.4g GHz', f3/1e9); end
    comments(end+1) = sprintf(['%s, sigma_t %.4g ps, shift %.4g ps ' ...
        '(%g sigma_t, linear phase), for causality with %s'], fd, sigma_t*1e12, shift*1e12, o.ShiftSigmas, snp_file);
end
if o.ExtraDelay > 0, comments(end+1) = sprintf('Extra delay %.4g ps (exact, linear phase) to line up with other ports', o.ExtraDelay*1e12); end
comments(end+1) = sprintf('Worst non-causal energy of the S-parameter responses: %.4g %% (input %.4g %%)', C.worst, Cin.worst);
comments(end+1) = "CST ASCII file with the same signal: " + string(cstFile);
if isGen
    comments(end+1) = sprintf('%s = CST [%s] * sqrt(Z0), Z0 = %g ohm: generator open-circuit voltage, source impedance Z0', ...
        mdfName, unitCST, info.Z0(k));
else
    comments(end+1) = sprintf('%s: ideal voltage source, CST values in %s used as-is', mdfName, unitCST);
end
comments(end+1) = sprintf('Excitation type %s; duration %.4g ns; -%g dB bandwidth %.4g to %.4g GHz', ...
    excType, Tdur*1e9, o.BwLeveldB, flo/1e9, fhi/1e9);
out.cstFile = string(cstFile);
if o.WriteMDF
    generate_mdf(tu, y*scale, mdfFile, 'BlockName', o.BlockName, 'Names', {'time', char(mdfName)}, 'Comments', comments.');
    out.mdfFile = string(mdfFile);
    say('  wrote %s\n        %s\n', cstFile, mdfFile);
else
    say('  wrote %s\n', cstFile);
end
end

%% ========================================================================
function [Su, fu, t, toTime, names, fData] = snpOnGrid(snp_file, df)
% non-zero S-parameter entries on a uniform grid 0..fmax (cached per file)
persistent cache
key = string(snp_file) + "|" + df;
if ~isempty(cache) && cache.key == key
    Su = cache.Su;  fu = cache.fu;  t = cache.t;  toTime = cache.toTime;  names = cache.names;  fData = cache.fData;  return
end
S  = sparameters(snp_file);  fS = S.Frequencies(:);  P = S.Parameters;  nP = size(P, 1);
fData = min(fS(fS > 0));                 % lowest S-parameter frequency other than DC
fu = (0:df:fS(end)).';  N = 2*(numel(fu) - 1);  dt = 1/(N*df);
t  = ((0:N-1).' - N/2) * dt;
toTime = @(X) fftshift(real(ifft([X; conj(X(end-1:-1:2))])));
pn = "P" + (1:nP);
tok = regexp(readlines(snp_file), '^!\s*Port\[(\d+)\]\s*(.+)$', 'tokens', 'once');
for q = find(~cellfun(@isempty, tok)).', pn(str2double(tok{q}{1})) = strtrim(string(tok{q}{2})); end
Su = {};  names = strings(0);
for i = 1:nP
    for j = 1:nP
        s = squeeze(P(i, j, :));
        if all(s == 0), continue, end
        x = interp1(fS, real(s), fu, 'pchip') + 1j*interp1(fS, imag(s), fu, 'pchip');
        x(1) = real(x(1));  Su{end+1} = x; %#ok<AGROW>
        names(end+1) = pn(i) + "," + pn(j); %#ok<AGROW>
    end
end
cache = struct('key', key, 'Su', {Su}, 'fu', fu, 't', t, 'toTime', toTime, 'names', names, 'fData', fData);
end

function C = checkCausal(Su, X, t, toTime, sig_dB, fmin, df, fLow)
% responses y = IFFT(S * X); non-causal energy = energy before t = 0 as % of
% each response's own energy; worst over responses within sig_dB of the strongest
% encLow: non-causal energy in 0 < f < fLow (below the lowest S-parameter point,
% DC excluded), as % of the same response energy
n = numel(Su);  E = zeros(n, 1);  Eneg = E;  Elow = E;  Y = cell(n, 1);
N = numel(t);  fk = min((0:N-1).', N - (0:N-1).') * df;  low = fk > 0 & fk < fLow;
for q = 1:n
    Y{q} = toTime(Su{q} .* X);
    if fmin > 0                          % energies above fmin only (Parseval on the FFT bins)
        N = numel(t);  fk = min((0:N-1).', N - (0:N-1).') * df;  hi = fk >= fmin;
        yn = Y{q} .* (t < 0);
        E(q) = sum(abs(fft(Y{q})).^2 .* hi);  Eneg(q) = sum(abs(fft(yn)).^2 .* hi);
    else
        E(q) = sum(Y{q}.^2);  Eneg(q) = sum(Y{q}(t < 0).^2);
    end
    Ynf = fft(Y{q} .* (t < 0));  Ef = sum(abs(fft(Y{q})).^2 .* (fk >= fmin));   % FFT-domain scale
    Elow(q) = E(q) * sum(abs(Ynf).^2 .* low) / Ef;
end
C.enc = 100 * Eneg ./ E;
C.sig = E >= max(E) * 10^(-sig_dB/10);
C.worst = max(C.enc(C.sig));
C.encLow = 100 * Elow ./ E;  C.worstLow = max(C.encLow(C.sig));
C.Y = Y;  C.t = t;
end

function plotCausality(inName, port, unitCST, t0, x0, tu, y, fu, Cin, C, names, o, filtered, f3, shift, spec, fc)
% excitation, its spectrum against the S-parameter band, and the non-causal
% energy of the S-parameter responses for the input and the filtered excitation
esc = @(s) strrep(strrep(string(s), '_', '\_'), '%', '\%');        % plain text -> LaTeX-safe text
uL  = regexprep(unitCST, '\^\((.*?)\)', '\$^{$1}\$');               % W^(1/2) -> W$^{1/2}$
fig = figure('Color', 'w', 'NumberTitle', 'off', 'Name', "Excitation causality " + inName);
setLatex(fig);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
if filtered && ~isempty(fc), ft = sprintf('filtered: Gaussian band-pass $-3$ dB at %.3g and %.3g GHz, shift %.0f ps', fc/1e9, shift*1e12);
elseif filtered, ft = sprintf('filtered: Gaussian low-pass $-3$ dB at %.3g GHz, shift %.0f ps', f3/1e9, shift*1e12);
else,        ft = 'meets the target, no filter'; end
title(tl, {sprintf('Causality of the S-parameter responses to %s ($i_{%d}$)', esc(inName), port), ...
    sprintf('worst non-causal energy %.4f\\%% (input) $\\rightarrow$ %.4f\\%%; target %g\\%%; %s', Cin.worst, C.worst, o.EncTolPct, ft)}, ...
    'Interpreter', 'latex');
dBn = @(Z) 20*log10(abs(Z) / max(abs(Z)) + realmin);
fs  = linspace(0, 1.5*fu(end), 3001).';
ax = nexttile(tl);  plot(ax, t0*1e9, x0, 'LineWidth', 1.1);  hold(ax, 'on');  grid(ax, 'on');
if filtered, plot(ax, tu*1e9, y, 'LineWidth', 1.1); legend(ax, 'input', 'filtered', 'Location', 'best'); end
xlabel(ax, 'Time (ns)');  ylabel(ax, "$i_{" + port + "}$ (" + uL + ")");  title(ax, 'Excitation');
ax = nexttile(tl);  plot(ax, fs/1e9, dBn(spec(t0, x0, fs)), 'LineWidth', 1.1);  hold(ax, 'on');  grid(ax, 'on');
if filtered, plot(ax, fs/1e9, dBn(spec(tu, y, fs)), 'LineWidth', 1.1); legend(ax, 'input', 'filtered', 'Location', 'best'); end
xline(ax, fu(end)/1e9, 'k--', 'S-parameter $f_{\max}$', 'Interpreter', 'latex', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
ylim(ax, [-100 3]);  xlabel(ax, 'Frequency (GHz)');  ylabel(ax, 'dB re peak');  title(ax, 'Excitation spectrum');
ax = nexttile(tl);
qs = find(Cin.sig);  [~, ord] = sort(Cin.enc(qs), 'descend');  qs = qs(ord);
semilogy(ax, 1:numel(qs), max(Cin.enc(qs), 1e-7), 'o-', 'LineWidth', 1.2);  hold(ax, 'on');  grid(ax, 'on');
if filtered, semilogy(ax, 1:numel(qs), max(C.enc(qs), 1e-7), 's-', 'LineWidth', 1.2); end
yline(ax, o.EncTolPct, '--', sprintf('%g\\%%', o.EncTolPct), 'Color', [0 0.6 0], 'LineWidth', 1.3, ...
    'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
set(ax, 'XTick', 1:numel(qs), 'XTickLabel', esc(names(qs)), 'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 60);
ylabel(ax, 'non-causal energy (\%)');  title(ax, sprintf('Responses within %g dB of the strongest', o.SigdB));
if filtered, legend(ax, 'input excitation', 'filtered excitation', 'Location', 'northeast'); end
ax = nexttile(tl);  [~, m] = max(Cin.enc .* Cin.sig);
yi = Cin.Y{m};  plot(ax, Cin.t*1e9, 20*log10(abs(yi)/max(abs(yi)) + realmin), 'LineWidth', 1.1);  hold(ax, 'on');  grid(ax, 'on');
lg = sprintf('input: %.4f\\%% before $t = 0$', Cin.enc(m));
if filtered
    yf = C.Y{m};  plot(ax, C.t*1e9, 20*log10(abs(yf)/max(abs(yf)) + realmin), 'LineWidth', 1.1);
    lg = {lg, sprintf('filtered: %.4f\\%% before $t = 0$', C.enc(m))};
end
xline(ax, 0, 'k-', 'HandleVisibility', 'off');  ylim(ax, [-100 3]);  xlim(ax, [-5 20]);
xlabel(ax, 'Time (ns)');  ylabel(ax, 'dB re peak');  legend(ax, lg, 'Location', 'northeast');
title(ax, "Worst response S(" + esc(names(m)) + ")", 'Interpreter', 'latex');

%% second figure: grid of the time-domain response of every significant entry
nq  = numel(qs);  nc = ceil(sqrt(nq));  nr = ceil(nq/nc);
fig2 = figure('Color', 'w', 'NumberTitle', 'off', 'Name', "Causality responses " + inName);
setLatex(fig2);
tl2  = tiledlayout(fig2, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, {sprintf('Time-domain responses $y = \\mathrm{IFFT}(S\\,X)$ to %s ($i_{%d}$)', esc(inName), port), ...
    sprintf('responses within %g dB of the strongest, sorted by non-causal energy; target %g\\%%', o.SigdB, o.EncTolPct)}, ...
    'Interpreter', 'latex');
for r = 1:nq
    q  = qs(r);  ax = nexttile(tl2);
    yi = Cin.Y{q};  yd = 20*log10(abs(yi)/max(abs(yi)) + realmin);
    plot(ax, Cin.t*1e9, yd, 'LineWidth', 1);  hold(ax, 'on');  grid(ax, 'on');
    if filtered
        yf = C.Y{q};  plot(ax, C.t*1e9, 20*log10(abs(yf)/max(abs(yf)) + realmin), 'LineWidth', 1);
    end
    tEnd = Cin.t(find(yd > -100, 1, 'last'))*1e9;
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    ylim(ax, [-120 3]);  xlim(ax, [-5 max(20, tEnd)]);
    if filtered
        tt = sprintf('S(%s): %.2g\\%% $\\rightarrow$ %.2g\\%%', esc(names(q)), Cin.enc(q), C.enc(q));
    else
        tt = sprintf('S(%s): %.2g\\%%', esc(names(q)), Cin.enc(q));
    end
    title(ax, tt, 'Interpreter', 'latex', 'FontSize', 9, 'FontWeight', 'normal');
    if r == 1 && filtered, legend(ax, 'input', 'filtered', 'Location', 'northeast'); end
end
xlabel(tl2, 'Time (ns)', 'Interpreter', 'latex');  ylabel(tl2, 'dB re peak', 'Interpreter', 'latex');
end

function setLatex(fig)
% LaTeX interpreter (Computer Modern fonts) for every text, tick label, legend and colorbar of fig
set(fig, 'DefaultTextInterpreter', 'latex', 'DefaultAxesTickLabelInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', 'DefaultColorbarTickLabelInterpreter', 'latex');
end

function H = gaussLP(f, f3, shift, band)
% linear-phase delay shift times
%   low-pass  (band = []):        Gaussian, -3 dB at f3
%   band-pass (band = [flo fhi]): Gaussian in f^2, -3 dB at flo and fhi,
%             exp(-k ((f^2 - f0^2)/B)^2), f0^2 = (fhi^2 + flo^2)/2, B = (fhi^2 - flo^2)/2
%             (even and smooth in f, so its impulse response is compact)
k = log(10)*3/20;
if nargin < 4 || isempty(band)
    H = exp(-k*(f/f3).^2);
else
    f0s = (band(2)^2 + band(1)^2)/2;  B = (band(2)^2 - band(1)^2)/2;
    H = exp(-k*((f.^2 - f0s)/B).^2);
end
H = H .* exp(-1j*2*pi*f*shift);
end

function sg = bpSigma(band)
% equivalent impulse-response std of the band-pass: sqrt(2) x the rms duration of
% |h|^2, sqrt(int H'^2 df / int H^2 df)/(2 pi) (equal to sigma_t for the low-pass)
f  = linspace(0, 3*band(2), 200001).';  H = gaussLP(f, Inf, 0, band);  dH = gradient(H, f);
sg = sqrt(2) * sqrt(trapz(f, dH.^2) / trapz(f, H.^2)) / (2*pi);
end

function [tu, y] = applyGauss(t, x, dt, f3, shift, tailT, fc)
if nargin < 7, fc = []; end
tu = (t(1) : dt : t(end) + shift + tailT).';
xu = interp1(t, x, tu, 'linear', 0);
Nf = 2^nextpow2(2*numel(tu));
fk = (0:Nf-1).' / (Nf*dt);  fk(fk >= 1/(2*dt)) = fk(fk >= 1/(2*dt)) - 1/dt;
y  = real(ifft(fft(xu, Nf) .* gaussLP(fk, f3, shift, fc)));
y  = y(1:numel(tu));
end

function [flo, fhi] = bandEdges(f, X, levelDB)
thr = max(X) * 10^(-levelDB/20);  [~, ip] = max(X);
il = find(X(1:ip) < thr, 1, 'last');  iu = ip - 1 + find(X(ip:end) < thr, 1, 'first');
if isempty(il), flo = 0; else, flo = interp1(X(il:il+1), f(il:il+1), thr); end
if isempty(iu), fhi = f(end); else, fhi = interp1(X(iu-1:iu), f(iu-1:iu), thr); end
end

function s = ghzStr(f)
s = string(num2str(round(f/1e9, 2))) + "GHz";
end

function s = engStr(x, unit)
pre = ["f" "p" "n" "u" "m" "" "k" "M" "G" "T"];
if x == 0, s = "0" + unit; return, end
x = round(x, 2, 'significant');
e = floor(log10(abs(x)) / 3) * 3;  e = min(max(e, -15), 12);
s = string(num2str(x / 10^e)) + pre(e/3 + 6) + unit;
end


function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
