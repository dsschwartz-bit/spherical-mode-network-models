function [f, X, h] = plotCSTSignal(cstFile, varargin)
%PLOTCSTSIGNAL  Plot a CST time signal and its spectrum (magnitude and phase).
%
%   plotCSTSignal(cstFile)
%   plotCSTSignal(cstFile, 'FMax', 3e9, 'Signal', 'i2', 'Unwrap', true)
%   [f, X] = plotCSTSignal(...)
%
%   Reads the signal with importCSTSignal.m (so W^(1/2) generator signals
%   are scaled by sqrt(Z0) and V signals are used as-is), and plots v(t),
%   the spectrum magnitude in dB and the spectrum phase. The spectrum is
%   single-sided. The zoomed plots (magnitude and phase) cover the band
%   where |X| is within 'BandLevel' dB (20) of its peak, or 'FRange' if
%   given; a second magnitude plot runs from 0 to Fs/2, the Nyquist
%   frequency (Fs = 1/median time step).
%
%       X(f) = integral v(t) exp(-j 2 pi f t) dt      [V/Hz]
%
%   evaluated directly at the requested frequencies. The time samples do
%   not need to be evenly spaced:
%     'nufft'  (default) non-uniform DFT, sum_n v(t_n) w_n exp(-j2pi f t_n),
%              with trapezoidal weights w_n from the actual sample spacing.
%              Exact for the samples given; no interpolation.
%     'interp' resample onto a uniform grid (pchip) and use fft (for
%              comparison).
%
%   Options
%     'Signal'    : signal index or CST name when the file holds several (1)
%     'Z0'        : port impedance if a W^(1/2) file has no Z0 parameter
%     'FMax'      : highest frequency plotted in Hz, or 'nyquist' for
%                   Fs/2. Default: the effective maximum frequency (see
%                   'Threshold'), rounded up to a round number.
%     'FRange'    : [fmin fmax] in Hz, the plotted frequency range; overrides
%                   FMax. Useful when broadband noise (e.g. co-simulation
%                   exports) defeats the automatic limit. Default [] (auto).
%     'BandLevel' : dB below the spectrum peak (default 20, as bw_level_dB in
%                   cst_to_ads_mdf). Unless FRange is given, the zoomed range
%                   is the band at that level: the first frequencies on each
%                   side of the peak where |X| drops BandLevel dB below it.
%                   [] falls back to FMax / Threshold.
%
%   Layout (2 x 2): time signal | phase (zoomed range)
%                   magnitude (zoomed range) | magnitude 0 to Nyquist, with
%                   the zoomed range shaded.
%     'Threshold' : dB below the spectrum peak that counts as negligible
%                   when finding the effective maximum frequency (60)
%     'NumFreq'   : number of frequency points over the plotted range. Default:
%                   at least 1001, and at least a spacing of 1/(2T),
%                   T = signal duration.
%     'Method'    : 'nufft' (default) or 'interp'
%     'Normalize' : false (default) plots 20*log10|X| in dB re 1 V/Hz;
%                   true plots dB relative to the spectrum peak
%     'Unwrap'    : false (default) plots phase wrapped to +-180 deg;
%                   true unwraps it
%     'PhaseDelay': delay tau removed from the phase plot, i.e. the phase of
%                   X(f) exp(+j 2 pi f tau) is plotted, taking out the linear
%                   phase of a time shift. 'fit' (default): tau from a
%                   least-squares fit of the unwrapped phase over the zoomed
%                   range, weighted by |X|^2 (the |X|^2-weighted mean group
%                   delay). A number: tau in seconds. 'none' or 0: raw phase
%                   (relative to t = 0). The delay is marked on the time plot.
%                   With a delay removed the phase is plotted unwrapped.
%
%   Outputs
%     f : frequency vector (Hz), X : complex spectrum (V/Hz),
%     h : handles (time, phase, mag, magFull, figure), plus h.fEff
%         (effective maximum frequency, Hz), h.fNyquist (Hz) and h.tau (s)
%
%   Example
%     plotCSTSignal("C:\path\to\CST\project\Export\time_signals\excitation1.txt");

p = inputParser;
addParameter(p, 'Signal',    1);
addParameter(p, 'Z0',        []);
addParameter(p, 'FMax',      []);
addParameter(p, 'FRange',    []);
addParameter(p, 'Threshold', 60);
addParameter(p, 'BandLevel', 20);
addParameter(p, 'NumFreq',   []);
addParameter(p, 'Method',    'nufft');
addParameter(p, 'Normalize', false);
addParameter(p, 'Unwrap',    false);
addParameter(p, 'PhaseDelay', 'fit');
parse(p, varargin{:});
o = p.Results;

% ---- read ----------------------------------------------------------------
[tc, vc, info] = importCSTSignal(cstFile, o.Z0);
k = o.Signal;
if ischar(k) || isstring(k)
    k = find(strcmpi(info.names, char(k)), 1);
    if isempty(k)
        error('plotCSTSignal:signal', 'Signal not found. Available: %s', strjoin(info.names, ', '));
    end
end
t = tc{k}(:);  v = vc{k}(:);
[t, iu] = unique(t);  v = v(iu);           % sorted, no repeated times
[~, srcName, srcExt] = fileparts(char(cstFile));

% ---- effective maximum frequency (searched up to Fs/2) ---------------------
dtMed = median(diff(t));
Fs    = 1 / dtMed;
fNyq  = Fs / 2;
T     = t(end) - t(1);
fAll  = linspace(0, fNyq, ceil(fNyq * 2*T) + 1).';    % 1/(2T) spacing
XAll  = spectrum(t, v, fAll, 'nufft');
magAll = abs(XAll);
% highest frequency (anywhere up to Nyquist) still within Threshold dB of the peak
above  = find(magAll >= max(magAll) * 10^(-o.Threshold/20), 1, 'last');
fEff   = fAll(min(above + 1, numel(fAll)));

% ---- frequency grid (single-sided, 0 .. FMax) ------------------------------
fMin = 0;
if ~isempty(o.FRange)
    if numel(o.FRange) ~= 2 || o.FRange(2) <= o.FRange(1) || o.FRange(1) < 0
        error('plotCSTSignal:frange', 'FRange must be [fmin fmax] in Hz with 0 <= fmin < fmax.');
    end
    fMin = o.FRange(1);  fMax = o.FRange(2);
elseif ~isempty(o.BandLevel) && (isempty(o.FMax) || ~ismember('BandLevel', p.UsingDefaults))
    [fMin, fMax] = bandEdges(t, v, fAll, magAll, o.BandLevel);     % an explicit FMax wins over the default band
elseif isempty(o.FMax)
    fMax = min(niceCeil(fEff), fNyq);
elseif (ischar(o.FMax) || isstring(o.FMax)) && strcmpi(o.FMax, 'nyquist')
    fMax = fNyq;
else
    fMax = o.FMax;
end
if fMax > fNyq * (1 + 1e-9)
    warning('plotCSTSignal:aboveNyquist', ...
        'FMax = %g Hz is above Nyquist (%g Hz); the spectrum there is aliased.', fMax, fNyq);
end
if isempty(o.NumFreq)
    nF = max(1001, ceil((fMax - fMin) * 2*T) + 1);
else
    nF = o.NumFreq;
end
bandMode = isempty(o.FRange) && ~isempty(o.BandLevel) && (isempty(o.FMax) || ~ismember('BandLevel', p.UsingDefaults));
f = linspace(fMin, fMax, nF).';                      % zoomed range
X = spectrum(t, v, f, lower(o.Method));
fF = unique([fAll; f]);                               % extended range 0 .. Nyquist
XF = spectrum(t, v, fF, 'nufft');
if ~isempty(o.FRange)
    how = 'set by FRange';
elseif bandMode
    how = sprintf('-%g dB band edges', o.BandLevel);
elseif ~isempty(o.FMax)
    how = 'set by FMax';
else
    how = sprintf('automatic: highest frequency within %g dB of the peak %.4g GHz', o.Threshold, fEff/1e9);
end
fprintf('%s: Nyquist %.4g GHz; zoomed spectrum %.4g to %.4g GHz (%s)\n', ...
    info.names{k}, fNyq/1e9, fMin/1e9, fMax/1e9, how);

% ---- plot: [time | phase (zoomed)] / [magnitude (zoomed) | magnitude to Nyquist]
fig = figure('Color', 'w', 'Name', [srcName ' ' info.names{k}]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
tOpt = {'Interpreter', 'latex', 'FontSize', 11};   % axes titles
title(tl, {texEsc([srcName srcExt]), texEsc(sprintf('%s (%s)', info.names{k}, info.sourceType{k}))}, ...
    'Interpreter', 'latex', 'FontSize', 12);
if o.Normalize
    dB = @(Y) 20*log10(abs(Y) / max(abs(XF)));  ylab = '|X(f)| (dB re peak)';
else
    dB = @(Y) 20*log10(abs(Y));                  ylab = '|X(f)| (dB re 1 V/Hz)';
end
[fs, fu] = freqScale(fMax);

h.time = nexttile(tl, 1);
plot(t * 1e9, v, 'LineWidth', 1.1);
grid on; box on;
xlabel('Time (ns)'); ylabel('Voltage (V)');  title('Time signal', tOpt{:});

% delay to de-embed from the phase (linear phase of a time shift)
pd = o.PhaseDelay;
if (ischar(pd) || isstring(pd)) && strcmpi(pd, 'none'), pd = 0; end
if (ischar(pd) || isstring(pd)) && strcmpi(pd, 'fit')
    w   = abs(X).^2;  ph = unwrap(angle(X));
    A   = [ones(size(f)), -2*pi*f] .* sqrt(w);
    cfit = A \ (ph .* sqrt(w));               % ph ~ c - 2 pi f tau, |X|^2-weighted
    tau = cfit(2);  tauHow = 'fit';
else
    tau = double(pd);  tauHow = 'given';
end
h.tau = tau;
if tau ~= 0, xline(h.time, tau * 1e9, '--', sprintf('$\\tau = %.4g$ ns', tau*1e9), 'Color', [0.85 0.33 0.10], 'Interpreter', 'latex', 'FontSize', 11, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'right'); end

h.phase = nexttile(tl, 2);
Xd = X .* exp(1j*2*pi*f*tau);                      % delay removed
if o.Unwrap || tau ~= 0
    phd = unwrap(angle(Xd));
    [~, ip] = max(abs(Xd));  phd = phd - 2*pi*round(phd(ip) / (2*pi));   % peak phase within +-180
    plot(f / fs, phd * 180/pi, 'LineWidth', 1.1);
else
    plot(f / fs, angle(Xd) * 180/pi, '.', 'MarkerSize', 4);
    ylim([-180 180]); yticks(-180:90:180);
end
grid on; box on;
xlabel(['Frequency (' fu ')']); ylabel('Phase (deg)');
if tau ~= 0
    title({'Phase of $X(f)\,\exp(j 2\pi f \tau)$', sprintf('$\\tau = %.4g$ ns (%s)', tau*1e9, tauHow)}, tOpt{:});
else
    title('Phase (re $t = 0$)', tOpt{:});
end
xlim([fMin, fMax] / fs);

h.mag = nexttile(tl, 3);
yz = dB(X);
plot(f / fs, yz, 'LineWidth', 1.1);
grid on; box on;
xlabel(['Frequency (' fu ')']); ylabel(ylab);
title({sprintf('Magnitude, %.4g to %.4g %s', fMin/fs, fMax/fs, fu), ['(' texEsc(how) ')']}, tOpt{:});
xlim([fMin, fMax] / fs);
ok = isfinite(yz);
if any(ok), ylim([min(yz(ok)) - 3, max(yz(ok)) + 3]); end
linkaxes([h.mag, h.phase], 'x');

h.magFull = nexttile(tl, 4);
[fsF, fuF] = freqScale(fNyq);
plot(fF / fsF, dB(XF), 'LineWidth', 1.0);  hold on
yl = ylim;
patch([fMin fMax fMax fMin] / fsF, yl([1 1 2 2]), [0.85 0.33 0.10], ...
    'FaceAlpha', 0.2, 'EdgeColor', [0.85 0.33 0.10], 'EdgeAlpha', 0.6, 'HandleVisibility', 'off');   % zoomed range
ylim(yl);  hold off
grid on; box on;
xlabel(['Frequency (' fuF ')']); ylabel(ylab);
title({sprintf('Magnitude, 0 to Nyquist (%.4g %s)', fNyq/fsF, fuF), ['shaded: zoomed range (' texEsc(how) ')']}, tOpt{:});
xlim([0, fNyq] / fsF);
h.figure   = fig;
h.fEff     = fEff;
h.fNyquist = fNyq;

if nargout == 0, clear f X h; end
end

% ========================================================================
function X = spectrum(t, v, f, method)
switch method
    case 'nufft'
        dt = diff(t);
        w  = [dt(1); dt(1:end-1) + dt(2:end); dt(end)] / 2;   % trapezoid weights
        if exist('nufft', 'file')
            X = nufft(v .* w, t, f);                  % sum v w exp(-j2pi t f)
        else
            X = exp(-2j*pi * f * t.') * (v .* w);     % direct sum (slower)
        end
    case 'interp'
        dt = median(diff(t));
        tu = (t(1):dt:t(end)).';
        vu = interp1(t, v, tu, 'pchip');
        N  = numel(tu);
        df = max(f(2) - f(1), eps);
        nfft = 2^nextpow2(max(N, round(1 / (dt * df))));   % bins at least as fine as f
        Xu = fft(vu, nfft) * dt .* exp(-2j*pi * (0:nfft-1).' / (nfft*dt) * tu(1));
        fu = (0:nfft-1).' / (nfft * dt);
        X  = interp1(fu, Xu, f);                    % onto the requested grid
    otherwise
        error('plotCSTSignal:method', 'Method must be ''nufft'' or ''interp''.');
end
X = X(:);
end

function [flo, fhi] = bandEdges(t, v, fc, mc, levelDB)
% band edges where |X| first drops levelDB below its peak on each side of
% the peak (coarse grid fc/mc, then refined on a fine grid)
[pk, ip] = max(mc);  thr = pk * 10^(-levelDB/20);
il = find(mc(1:ip) < thr, 1, 'last');
iu = ip - 1 + find(mc(ip:end) < thr, 1, 'first');
if isempty(il), flo = 0;      else, flo = refine(fc(il), fc(il+1)); end
if isempty(iu), fhi = fc(end); else, fhi = refine(fc(iu-1), fc(iu)); end
    function fx = refine(fa, fb)
        g = linspace(fa, fb, 401).';
        d = abs(spectrum(t, v, g, 'nufft')) - thr;
        j = find(sign(d(1:end-1)) ~= sign(d(2:end)), 1, 'first');
        if isempty(j), fx = (fa + fb)/2; return, end
        fx = g(j) - d(j) * (g(j+1) - g(j)) / (d(j+1) - d(j));
    end
end

function y = niceCeil(x)
% round x up to 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8 or 10 x a power of ten
if x <= 0, y = x; return, end
e = 10^floor(log10(x));
steps = [1 1.2 1.5 2 2.5 3 4 5 6 8 10];
y = steps(find(steps * e >= x * (1 - 1e-12), 1)) * e;
end

function [s, u] = freqScale(fMax)
if fMax >= 1e9,     s = 1e9; u = 'GHz';
elseif fMax >= 1e6, s = 1e6; u = 'MHz';
elseif fMax >= 1e3, s = 1e3; u = 'kHz';
else,               s = 1;   u = 'Hz';
end
end


function s = texEsc(s)
%TEXESC  Escape LaTeX special characters for a text-mode title.
s = regexprep(char(s), '([_%&#$])', '\\$1');
end
