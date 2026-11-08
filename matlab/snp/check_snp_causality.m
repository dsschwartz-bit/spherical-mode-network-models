%% CHECK_SNP_CAUSALITY
% Time-domain causality check of every entry S_ij(f) of a Touchstone file.
%
% Each S_ij is put on a uniform frequency grid from DC to fmax, band
% limited (see taper: Gaussian-sine, Gaussian, Kaiser or Bessel) and
% inverse-transformed to an impulse response h_ij(t) on t = -T/2 .. T/2
% (T = 1/df). A causal entry has no response before t = 0. Two things limit
% what "before t = 0" means here:
%   - the band limit: even a perfect delta is spread over a few 1/(2 fmax).
%     The band-limited delta is computed with the same settings and its
%     onset time t_ref is the resolution limit (leads are relative to it).
%   - the frequency step: the record is periodic with T = 1/df, so a
%     response still ringing at T/2 wraps round and appears at negative
%     time. The tail level near +T/2 is reported to flag this.
%
% Reported per entry (only entries that are not identically zero), for the
% raw data and for the band-limited data:
%   lead   : onset of |h| above onset_dB re its peak, minus t_ref.
%            lead < 0 means the response starts before t = 0 by that much
%            beyond the band-limit spread, i.e. it needs an advance.
%   Enc    : non-causal energy (before t = 0 for a causal taper, before
%            -|t_ref| for a zero-phase one) as a percentage of the total.
%   tail   : peak |h| over the last 10% before +T/2, re peak (dB);
%            high values mean wrap-around (time aliasing) is likely.
% plus the bulk delay that makes the significant entries (energy within
% sig_dB of the strongest) causal:
%   tau0 = max(-lead) + margin
% (delay S by exp(-j*w*tau0) and take tau0 off a delay elsewhere, e.g.
% the far-field VCVS delay r/c0 - tau0).
%
% Cross-check: RF Toolbox iscausal (R2023a+) on the same data, giving a
% per-entry causality metric (tolerance iscausal_tol dB). iscausal says
% whether an entry is causal; the impulse-response check says by how much
% it leads and what delay fixes it.
%
% The file's own frequency points are used where they lie on the uniform
% grid; only grid points missing from the file (e.g. between an added DC
% point and the first measured frequency) are filled by interpolation.
%
% With taper "gaussian" the filter parameters are saved as snp_taper.mat
% next to the SnP file, for shaping the excitation of later simulations.
%
% Requires RF Toolbox (sparameters, iscausal) and Signal Processing
% Toolbox (kaiser, besself).

% When called from another file (e.g. validation/cst_ads/time_signals/
% filter_excitation_gaussian.m) with snp_check_external set to a Touchstone
% path, that file is checked instead of snp_file below.
if ~exist('snp_check_external', 'var'), clear; clc; end

%% ===== USER SETTINGS =====
% Touchstone file to check (default: <repo>/output/SnP; this file is in <repo>/matlab/snp)
snp_file    = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'output', 'SnP', ...
    "S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p");
% =========================
if exist('snp_check_external', 'var') && strlength(snp_check_external) > 0
    snp_file = string(snp_check_external);
end

%% ---- settings ---------------------------------------------------------------
df_uniform  = 5e6;     % Hz, uniform grid step (record 1/df); [] = smallest step in the file
kaiser_beta = 6;      % Kaiser taper (0 = rectangular; larger = less ringing, wider spread).
                       % Entries still large at fmax (higher-order modes) need a strong
                       % taper, otherwise their band-edge truncation shows up as
                       % false non-causal leads.
taper       = "gaussiansine"; % band limiting before the inverse transform:
                       %  "gaussiansine": excitation-shaped Gaussian band-pass, set by gs_upper_dB at fmax,
                       %            gs_lower_dB at the lowest measured frequency (image term included)
                       %            and a delay of gs_shift_sigmas std. devs. of its envelope
                       %  "gaussian": zero-phase Gaussian (no dispersion; the limit of high-order Bessel)
                       %  "kaiser": half Kaiser window (zero phase, spreads to t < 0)
                       %  "bessel": analog Bessel low-pass (causal, nearly flat group
                       %            delay in band; the reference becomes the filter's own
                       %            impulse response, so leads are relative to it)
gs_upper_dB     = -40;  % gaussiansine: level at the highest frequency in the file (re peak)
gs_lower_dB     = -20;  % gaussiansine: level at gs_fmin (re peak), image term included
gs_fmin         = [];   % Hz, gaussiansine lower reference; [] = lowest non-zero frequency in the file
gs_shift_sigmas = 5;    % gaussiansine: delay in std. devs. of its time envelope (0 = zero phase)
gauss_f3        = "auto"; % Hz, Gaussian -3 dB frequency, or "auto": highest one meeting enc_tol_pct and sig_dB
band_GHz        = [2 3];  % excitation band (GHz) used to report the Gaussian's in-band amplitude tilt
gauss_shift_sigmas = 5;  % Gaussian delayed by this many std. devs. of its impulse response (causal filter)
gauss_require_lead = false; % also require no lead in entries within sig_dB (else energy criterion only)
gauss_energy_ref   = "entry"; % "entry": enc_tol_pct of each significant entry's own energy; "strongest": of the strongest entry
bessel_order    = "auto"; % Bessel filter order, or "auto": smallest order meeting enc_tol_pct and sig_dB
bessel_max_order = 20;  % upper limit of the automatic order search
enc_tol_pct     = 0.01;  % % non-causal energy allowed (before t = 0) in each entry within sig_dB of the strongest
sig_dB          = 20;    % dB: entries with energy within this of the strongest are the significant ones
bessel_fc_frac  = 2/3; % Bessel cutoff as a fraction of the highest frequency in the file
                       % (besself cutoff: group delay approximately flat up to it)
onset_dB    = -30;     % onset threshold re each entry's peak
test_delay  = 0;       % s, extra delay exp(-j*w*test_delay) applied to every
                       % entry before the check (try the recommended tau0)
margin      = 50e-12;  % s, margin added to the recommended tau0
plot_cols   = [1 2];   % source ports whose impulse responses are plotted
                       % (all destination ports of each)
iscausal_tol = -60;    % dB, relative tolerance for RF Toolbox iscausal
n_worst      = 6;      % worst causality violators shown in their own figure

%% ---- read -------------------------------------------------------------------
S  = sparameters(snp_file);
f  = S.Frequencies(:);  P = S.Parameters;  nP = size(P, 1);
names = "P" + (1:nP);                                   % port names from "! Port[N] Name" comments
L = readlines(snp_file);
tok = regexp(L, '^!\s*Port\[(\d+)\]\s*(.+)$', 'tokens', 'once');
for k = find(~cellfun(@isempty, tok)).'
    names(str2double(tok{k}{1})) = strtrim(string(tok{k}{2}));
end
fprintf('%s\n  %d ports, %d frequencies, %.4g-%.4g GHz, step %.3g-%.3g MHz\n', ...
    snp_file, nP, numel(f), f(1)/1e9, f(end)/1e9, min(diff(f))/1e6, max(diff(f))/1e6);
if f(1) > 0
    fprintf('  no DC point: DC taken as the real part of the lowest frequency\n');
end

%% ---- uniform grid, window, time axis ---------------------------------------
if isempty(df_uniform), df_uniform = min(diff(f)); end
fu = (0:df_uniform:f(end)).';  Nu = numel(fu);
N  = 2*(Nu - 1);  T = 1/df_uniform;  dt = T/N;
t  = ((0:N-1).' - N/2) * dt;                            % -T/2 .. T/2 (after fftshift)
toTime = @(X) fftshift(real(ifft([X; conj(X(end-1:-1:2))])));
fprintf('  uniform grid %.4g MHz (T = %.4g ns), fmax %.4g GHz (dt = %.4g ps)\n', ...
    df_uniform/1e6, T*1e9, fu(end)/1e9, dt*1e12);
filled = fu(min(abs(fu - f.'), [], 2) > 1e-6*df_uniform);
if isempty(filled)
    fprintf('  all grid points taken from the file\n');
else
    fprintf('  %d grid point(s) not in the file, interpolated: %s MHz\n', numel(filled), mat2str(filled.'/1e6, 5));
end

% every non-zero entry on the uniform grid (with the optional test delay)
Su = cell(nP);
for i = 1:nP
    for j = 1:nP
        s = squeeze(P(i, j, :));
        if all(abs(s) == 0), continue, end
        su = interp1(f, real(s), fu, 'pchip', 'extrap') + 1j*interp1(f, imag(s), fu, 'pchip', 'extrap');
        su = su .* exp(-1j*2*pi*fu*test_delay);
        su(1) = real(su(1));                                     % DC is real
        Su{i, j} = su;
    end
end

%% ---- band limiting and check of every entry ---------------------------------
% raw data: no band limiting (abrupt cut at fmax)
Rraw = evalTaper(Su, ones(Nu, 1), t, dt, onset_dB, false, sig_dB);
fprintf('\n  RAW S-parameters (no filter): worst non-causal energy of entries within %g dB of the strongest %.4f %% (own energy), %.4f %% of the strongest; worst significant lead %.0f ps\n', ...
    sig_dB, Rraw.EncSig, Rraw.EncMaxRel, Rraw.leadSig*1e12);

if taper == "bessel"
    fcB = bessel_fc_frac * fu(end);
    if isstring(bessel_order) || ischar(bessel_order)          % "auto": smallest sufficient order
        fprintf(['  automatic Bessel order (cutoff %.4g GHz): smallest order with non-causal energy <= %g %% ' ...
            'of the strongest entry and no lead in entries within %g dB of it\n'], fcB/1e9, enc_tol_pct, sig_dB);
        for nB = 1:bessel_max_order
            w = besselResp(nB, fcB, fu);
            R = evalTaper(Su, w, t, dt, onset_dB, true, sig_dB);
            fprintf('    order %2d: |H| at fmax %6.1f dB, at 2.5 GHz %5.1f dB; worst non-causal energy %.4f %%, worst significant lead %5.0f ps\n', ...
                nB, 20*log10(abs(w(end))), 20*log10(abs(interp1(fu, w, min(2.5e9, fu(end))))), R.EncMaxRel, R.leadSig*1e12);
            if R.EncMaxRel <= enc_tol_pct && R.leadSig >= -dt/2, break, end
        end
        if ~(R.EncMaxRel <= enc_tol_pct && R.leadSig >= -dt/2)
            fprintf('    not reached by order %d; using it anyway\n', bessel_max_order);
        end
        nOrd = nB;
    else
        nOrd = bessel_order;  w = besselResp(nOrd, fcB, fu);
        R = evalTaper(Su, w, t, dt, onset_dB, true, sig_dB);
    end
    tauB = -gradient(unwrap(angle(w)), 2*pi*fu);
    fprintf('  Bessel low-pass: order %d, cutoff %.4g GHz; |H| at fmax %.1f dB; group delay %.1f ps (DC) .. %.1f ps (cutoff)\n', ...
        nOrd, fcB/1e9, 20*log10(abs(w(end))), tauB(1)*1e12, interp1(fu, tauB, fcB)*1e12);
    tapTxt = sprintf('Bessel LP order %d, $f_c$ %.3g GHz, %.0f dB at $f_{\\max}$', nOrd, fcB/1e9, 20*log10(abs(w(end))));
elseif taper == "gaussiansine"
    if isempty(gs_fmin), gs_fmin = min(f(f > 0)); end
    fgrid = linspace(0, fu(end), 4001).';
    Gs  = @(ff, fc, sf) exp(-(ff - fc).^2/(2*sf^2)) + exp(-(ff + fc).^2/(2*sf^2));   % with image term
    lvl = @(ff, fc, sf) 20*log10(Gs(ff, fc, sf) / max(Gs(fgrid, fc, sf)));
    cost = @(p) (lvl(gs_fmin, p(1)*1e9, p(2)*1e9) - gs_lower_dB)^2 + (lvl(fu(end), p(1)*1e9, p(2)*1e9) - gs_upper_dB)^2;
    p0 = [(gs_fmin + fu(end))/2, (fu(end) - gs_fmin)/5] / 1e9;
    p  = fminsearch(cost, p0, optimset('TolX', 1e-9, 'TolFun', 1e-12, 'MaxFunEvals', 4000));
    fcG = p(1)*1e9;  sfG = p(2)*1e9;  stG = 1/(2*pi*sfG);  tauG = gs_shift_sigmas*stG;
    w = Gs(fu, fcG, sfG) .* exp(-1j*2*pi*fu*tauG);  w = w / max(abs(w));  w(1) = real(w(1));
    Lg = lvl(fgrid, fcG, sfG);  [~, ipk] = max(Lg);
    e20 = [interp1(Lg(1:ipk), fgrid(1:ipk), -20), interp1(Lg(ipk:end), fgrid(ipk:end), -20)];
    fprintf(['  Gaussian-sine filter: %.1f dB at %.3g GHz, %.1f dB at %.3g GHz -> f_c %.4g GHz, sigma_f %.4g GHz; ' ...
        '-20 dB band %.3g - %.3g GHz; envelope sigma_t %.1f ps, delay %.1f ps (%g sigma_t)\n'], ...
        lvl(gs_fmin, fcG, sfG), gs_fmin/1e9, lvl(fu(end), fcG, sfG), fu(end)/1e9, fcG/1e9, sfG/1e9, e20/1e9, ...
        stG*1e12, tauG*1e12, gs_shift_sigmas);
    R = evalTaper(Su, w, t, dt, onset_dB, gs_shift_sigmas >= 3, sig_dB);
    tapTxt = sprintf('Gaussian sine $f_c$ %.3g GHz, %.0f dB at %.2g GHz, %.0f dB at %.2g GHz, delay %g$\\sigma$', ...
        fcG/1e9, gs_lower_dB, gs_fmin/1e9, gs_upper_dB, fu(end)/1e9, gs_shift_sigmas);
elseif taper == "gaussian"
    % Gaussian low-pass, delayed by gauss_shift_sigmas standard deviations of its
    % impulse response so that the filter itself is causal to that accuracy
    sigt = @(f3) 1 ./ (2*pi * f3/sqrt(2*log(10)*3/20));     % std. dev. of the impulse response
    gW   = @(f3) exp(-(fu/f3).^2 * log(10)*3/20) .* exp(-1j*2*pi*fu*gauss_shift_sigmas*sigt(f3));
    tlt  = @(f3) 3*((band_GHz(2)*1e9/f3)^2 - (band_GHz(1)*1e9/f3)^2);   % in-band tilt, dB
    if isstring(gauss_f3) || ischar(gauss_f3)              % "auto": highest feasible f3
        f3s = (0.25:0.025:fu(end)/1e9) * 1e9;  okF = false(size(f3s));
        A = zeros(size(f3s));  E = A;  Ls = A;
        for q = 1:numel(f3s)
            wq = gW(f3s(q));  wq(1) = real(wq(1));  Rq = evalTaper(Su, wq, t, dt, onset_dB, true, sig_dB);
            A(q) = 20*log10(abs(wq(end)));  Ls(q) = Rq.leadSig;
            if gauss_energy_ref == "entry", E(q) = Rq.EncSig; else, E(q) = Rq.EncMaxRel; end
            okF(q) = E(q) <= enc_tol_pct && (~gauss_require_lead || Ls(q) >= -dt/2);
        end
        qb = find(okF, 1, 'last');
        if isempty(qb), error('No Gaussian cutoff meets the causality conditions.'); end
        if gauss_energy_ref == "entry"
            crit = sprintf('energy before t = 0 <= %g %% of its own energy in every entry within %g dB of the strongest', enc_tol_pct, sig_dB);
        else
            crit = sprintf('energy before t = 0 <= %g %% of the strongest entry', enc_tol_pct);
        end
        if gauss_require_lead, crit = [crit sprintf(' and no lead in entries within %g dB of it', sig_dB)]; end
        fprintf('  automatic Gaussian cutoff (shift %g sigma_t): highest -3 dB frequency with %s\n', gauss_shift_sigmas, crit);
        fprintf('     f3 GHz   sigma_t   shift    |H| at fmax   tilt %g-%g GHz   worst Enc %%   worst significant lead\n', band_GHz);
        for q = max(1, qb-4):min(numel(f3s), qb+4)
            mk = '';  if q == qb, mk = '  <- chosen'; elseif ~okF(q), mk = '  (fails)'; end
            fprintf('     %5.3f   %5.1f ps  %5.0f ps  %8.1f dB    %8.2f dB     %8.4f    %6.0f ps%s\n', f3s(q)/1e9, ...
                sigt(f3s(q))*1e12, gauss_shift_sigmas*sigt(f3s(q))*1e12, A(q), tlt(f3s(q)), E(q), Ls(q)*1e12, mk);
        end
        f3 = f3s(qb);
    else
        f3 = gauss_f3;
    end
    w = gW(f3);  w(1) = real(w(1));  R = evalTaper(Su, w, t, dt, onset_dB, true, sig_dB);
    fprintf(['  Gaussian taper: -3 dB at %.4g GHz; |H| at fmax %.1f dB; tilt over %g-%g GHz %.2f dB; ' ...
        'sigma_t %.1f ps, shift %.1f ps (%g sigma_t); no dispersion (linear phase)\n'], f3/1e9, 20*log10(abs(w(end))), ...
        band_GHz, tlt(f3), sigt(f3)*1e12, gauss_shift_sigmas*sigt(f3)*1e12, gauss_shift_sigmas);
    tapTxt = sprintf('Gaussian, $-3$ dB at %.3g GHz, shift %.0f ps, %.0f dB at $f_{\\max}$', f3/1e9, gauss_shift_sigmas*sigt(f3)*1e12, 20*log10(abs(w(end))));
    % save the filter parameters (for shaping CST / ADS excitations)
    taper_params = struct('type', "gaussian", 'f3', f3, 'sigma_t', sigt(f3), ...
        'shift', gauss_shift_sigmas*sigt(f3), 'shift_sigmas', gauss_shift_sigmas, ...
        'att_fmax_dB', 20*log10(abs(w(end))), 'fmax', fu(end), 'snp_file', snp_file, ...
        'enc_tol_pct', enc_tol_pct, 'sig_dB', sig_dB, 'created', string(datetime('now')));
    save(fullfile(fileparts(snp_file), 'snp_taper.mat'), 'taper_params');
    fprintf('  filter saved to %s\n', fullfile(fileparts(snp_file), 'snp_taper.mat'));
else
    wk = kaiser(2*Nu - 1, kaiser_beta);  w = wk(Nu:end);    % half window, 1 at DC
    R = evalTaper(Su, w, t, dt, onset_dB, false, sig_dB);
    tapTxt = sprintf('Kaiser $\\beta = %g$', kaiser_beta);
end
lead = R.lead;  Enc = R.Enc;  EncRel = R.EncRel;  tail = R.tail;  H = R.H;  href = R.href;  t_ref = R.t_ref;
fprintf('\n  %-34s %14s %18s\n', 'causality of entries within sig_dB', 'impulse', 'filtered');
fprintf('  %-34s %13.4f%% %17.4f%%\n', 'worst non-causal energy (own)', Rraw.EncSig, R.EncSig);
fprintf('  %-34s %13.4f%% %17.4f%%\n', 'worst non-causal energy (re max)', Rraw.EncMaxRel, R.EncMaxRel);
fprintf('  %-34s %11.0f ps %15.0f ps\n', 'worst lead', Rraw.leadSig*1e12, R.leadSig*1e12);
fprintf('  %-34s %14d %18d\n', 'entries leading > 1 time step', nnz(Rraw.lead(:) < -dt & Rraw.sig(:)), nnz(R.lead(:) < -dt & R.sig(:)));
fprintf('  reference (onset of band-limited delta at %g dB): t_ref = %.1f ps\n', onset_dB, t_ref*1e12);
fprintf('  worst non-causal energy %.4f %% of the strongest entry; worst lead of entries within %g dB of it: %.0f ps\n', ...
    R.EncMaxRel, sig_dB, R.leadSig*1e12);

%% ---- report ---------------------------------------------------------------------
[~, ord] = sort(lead(:));  ord = ord(isfinite(lead(ord)));
fprintf('\nentries ordered by lead (negative = starts before t = 0 beyond the window spread)\n');
fprintf('  %-22s %9s %9s %10s %12s %9s\n', 'S(to, from)', 'lead ps', 'Enc %', 'impulse Enc %', 'Enc re max %', 'tail dB');
for q = ord.'
    [i, j] = ind2sub([nP nP], q);
    flag = '';
    if lead(i, j) < -dt && R.sig(i, j), flag = '  <- non-causal';
    elseif lead(i, j) < -dt, flag = sprintf('  (leads, but insignificant: below %g dB of the strongest)', sig_dB); end
    if tail(i, j) > onset_dB, flag = [flag '  (tail high: wrap-around possible)']; end %#ok<AGROW>
    fprintf('  S(%-9s,%-9s) %9.1f %9.4f %10.4f %12.4f %9.1f%s\n', names(i), names(j), lead(i, j)*1e12, Enc(i, j), Rraw.Enc(i, j), EncRel(i, j), tail(i, j), flag);
end
tau0 = max(0, -R.leadSig) + margin;                % significant entries only
fprintf('\nlargest lead of entries within %g dB of the strongest: %.1f ps -> bulk delay to make them causal: tau0 = %.0f ps (incl. %.0f ps margin)\n', sig_dB, ...
    -R.leadSig*1e12, tau0*1e12, margin*1e12);
if test_delay > 0, fprintf('(results include test_delay = %.0f ps)\n', test_delay*1e12); end

%% ---- cross-check: RF Toolbox iscausal ---------------------------------------------
Pd = P .* reshape(exp(-1j*2*pi*f*test_delay), 1, 1, []);        % same test delay
try
    [isC, resmat] = iscausal(f, Pd, iscausal_tol);   % isC: overall result; resmat: metric per entry
    nz = any(P ~= 0, 3);                              % entries with data
    fprintf('\nRF Toolbox iscausal (tol %g dB): overall %s\n', iscausal_tol, string(all(isC(:))));
    if isequal(size(resmat), [nP nP])
        fprintf('  metric range over non-zero entries: %.1f to %.1f dB\n', min(resmat(nz)), max(resmat(nz)));
        bad = find(nz & ~(resmat <= iscausal_tol));
        [~, o] = sort(resmat(bad), 'descend');  bad = bad(o);
        for q = bad.'
            [i, j] = ind2sub([nP nP], q);
            fprintf('  S(%-11s,%-11s) metric %7.1f dB  (time-domain lead %7.1f ps)\n', ...
                names(i), names(j), resmat(i, j), lead(i, j)*1e12);
        end
    else
        fprintf('  metric array size %s\n', mat2str(size(resmat)));
    end
catch err
    fprintf('\nRF Toolbox iscausal failed: %s\n', err.message);
end

%% ---- figures ----------------------------------------------------------------------
[~, fn, fe] = fileparts(snp_file);
esc = @(s) strrep(strrep(string(s), '_', '\_'), '%', '\%');     % plain text -> LaTeX-safe text
namesL = esc(names);                                            % port names for LaTeX labels
% LaTeX interpreter (Computer Modern fonts) for all text in the figures made below
ltxDefaults = {'DefaultTextInterpreter', 'latex', 'DefaultAxesTickLabelInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', 'DefaultColorbarTickLabelInterpreter', 'latex'};
% Figures: 1 raw (+ worst entries), 3 raw vs filtered comparison, 2 filtered.
% Each lead map has its own symmetric colour scale.
cases = {Rraw, 'impulse responses of the S-parameter block (no excitation shaping)', 'SnP causality raw'; ...
         R, ['filtered responses of the S-parameter block (' tapTxt ')'], 'SnP causality filtered response'};
for cs = [1 3 2]
    if cs == 3        % comparison figure
    figC = figure('Color', 'w', 'NumberTitle', 'off', 'Name', 'SnP causality impulse vs filtered response', ltxDefaults{:});
    tlC  = tiledlayout(figC, 2, 2, 'TileSpacing', 'compact');
    title(tlC, {sprintf('Causality of %s%s: impulse responses vs filtered responses (%s)', esc(fn), fe, tapTxt), ...
        sprintf('entries within %g dB of the strongest; non-causal energy as \\%% of each entry''s own energy', sig_dB)}, 'Interpreter', 'latex');
    % filter magnitude
    ax = nexttile(tlC, 1);  plot(ax, fu/1e9, 20*log10(abs(w)/max(abs(w)) + realmin), 'LineWidth', 1.3);  grid(ax, 'on');  hold(ax, 'on');
    xline(ax, fu(end)/1e9, 'k--', '$f_{\max}$', 'Interpreter', 'latex', 'HandleVisibility', 'off');
    if exist('gs_fmin', 'var') && ~isempty(gs_fmin), xline(ax, gs_fmin/1e9, 'k--', '$f_{\min}$', 'Interpreter', 'latex', 'HandleVisibility', 'off'); end
    ylim(ax, [-80 3]);  xlabel(ax, 'Frequency (GHz)');  ylabel(ax, '$|H|$ (dB re peak)');  title(ax, 'Filter (applied to the input signal)');
    % non-causal energy per significant entry, raw vs filtered (one of each reciprocal pair)
    ax = nexttile(tlC, 2);
    qs = find(Rraw.sig & triu(true(nP)));  [~, o3] = sort(Rraw.Enc(qs), 'descend');  qs = qs(o3);
    lbl = strings(numel(qs), 1);  for m = 1:numel(qs), [ii, jj] = ind2sub([nP nP], qs(m)); lbl(m) = namesL(ii) + "," + namesL(jj); end
    semilogy(ax, 1:numel(qs), max(Rraw.Enc(qs), 1e-6), 'o-', 1:numel(qs), max(R.Enc(qs), 1e-6), 's-', 'LineWidth', 1.2);  grid(ax, 'on');  hold(ax, 'on');
    yline(ax, enc_tol_pct, '--', sprintf('%g\\%%', enc_tol_pct), 'Color', [0 0.6 0], 'LineWidth', 1.3, 'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    set(ax, 'XTick', 1:numel(qs), 'XTickLabel', lbl, 'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 60);
    ylabel(ax, 'non-causal energy (\%)');  legend(ax, 'impulse response', 'filtered response', 'Location', 'northeast');  title(ax, 'Per entry');
    % impulse response of the worst raw entry, raw vs filtered
    [~, m] = max(Rraw.Enc(qs));  [ii, jj] = ind2sub([nP nP], qs(m));
    ax = nexttile(tlC, [1 2]);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
    hr = Rraw.H{ii, jj};  hf = R.H{ii, jj};
    plot(ax, t*1e9, 20*log10(abs(hr)/max(abs(hr)) + realmin), 'LineWidth', 1.1);
    plot(ax, t*1e9, 20*log10(abs(hf)/max(abs(hf)) + realmin), 'LineWidth', 1.1);
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');  ylim(ax, [-100 3]);  xlim(ax, [-5 10]);
    xlabel(ax, 'Time (ns)');  ylabel(ax, '$|h(t)|$ (dB re peak)');
    legend(ax, sprintf('impulse response: %.3f\\%% before $t = 0$', Rraw.Enc(ii, jj)), sprintf('filtered response: %.4f\\%% before $t = 0$', R.Enc(ii, jj)), 'Location', 'northeast');
    title(ax, sprintf('Worst entry S(%s, %s): impulse response vs filtered response', namesL(ii), namesL(jj)), 'Interpreter', 'latex');
        continue
    end
    Rc = cases{cs, 1};
    fig = figure('Color', 'w', 'NumberTitle', 'off', 'Name', cases{cs, 3}, ltxDefaults{:});
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    title(tl, {sprintf('Causality of %s%s: %s', esc(fn), fe, cases{cs, 2}), ...
        sprintf('onset at %g dB; worst non-causal energy of entries within %g dB of the strongest %.4f\\%%; worst lead %.0f ps', ...
        onset_dB, sig_dB, Rc.EncSig, Rc.leadSig*1e12)}, 'Interpreter', 'latex');
    ax = nexttile(tl);
    imagesc(ax, Rc.lead*1e12, 'AlphaData', isfinite(Rc.lead));
    axis(ax, 'image');  cb = colorbar(ax);  cb.Label.String = 'lead (ps), $<0$ non-causal';  cb.Label.Interpreter = 'latex';
    mLead = max(abs(Rc.lead(:))*1e12, [], 'omitnan');  if ~(mLead > 0), mLead = 1; end
    clim(ax, [-mLead mLead]);  colormap(ax, flipud(redblue(256)));
    set(ax, 'XTick', 1:nP, 'YTick', 1:nP, 'XTickLabel', namesL, 'YTickLabel', namesL, ...
        'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 60, 'Color', [0.9 0.9 0.9]);
    xlabel(ax, 'from port $j$');  ylabel(ax, 'to port $i$');
    title(ax, {'Onset lead of $h_{ij}(t)$', 'grey: $S_{ij}$ identically zero in the file (no coupling)'}, 'Interpreter', 'latex');

    ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
    lg = strings(0);
    for j = plot_cols
        for i = 1:nP
            if isempty(Rc.H{i, j}) || i == j, continue, end
            plot(ax, t*1e9, 20*log10(abs(Rc.H{i, j})/max(abs(Rc.H{i, j})) + realmin), 'LineWidth', 1);
            lg(end+1) = sprintf('S(%s, %s)', namesL(i), namesL(j)); %#ok<SAGROW>
        end
    end
    plot(ax, t*1e9, 20*log10(Rc.href + realmin), 'k--', 'LineWidth', 1.3);
    if cs == 1, lg(end+1) = "band-limited delta (reference)"; else, lg(end+1) = "filter, time domain (reference)"; end
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');  xline(ax, Rc.t_ref*1e9, 'k:', 'HandleVisibility', 'off');
    yline(ax, onset_dB, 'k:', 'HandleVisibility', 'off');
    xlim(ax, [min(-1.5, (Rc.t_ref + min(Rc.lead(:)))*1e9 - 0.5), max(3, Rc.t_ref*1e9 + 3)]);  ylim(ax, [-80 3]);
    xlabel(ax, 'Time (ns)');  ylabel(ax, '$|h(t)|$ (dB re peak)');
    if cs == 1, ttl = 'Impulse responses (dotted: $t_{\mathrm{ref}}$ and onset level)';
    else,       ttl = 'Filtered responses (dotted: $t_{\mathrm{ref}}$ and onset level)'; end
    title(ax, ttl, 'Interpreter', 'latex');
    legend(ax, lg, 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 7);
    if cs == 1   % worst raw entries, shown right after the raw figure
    % The n_worst significant entries with the most negative lead (shown even when all are
    % causal within the resolution); only one of each reciprocal pair. The S-parameter
    % block itself is never filtered; any filter is applied to the excitation.
    ok = find(isfinite(Rraw.lead(:)) & Rraw.sig(:));       % significant entries; ties: most non-causal energy first
    [~, o2] = sortrows([Rraw.lead(ok) -Rraw.Enc(ok)]);  cand = ok(o2).';  pick = zeros(0, 2);
    nViol = nnz(Rraw.lead(:) < -dt & Rraw.sig(:));          % significant entries only
    for q = cand
        [i, j] = ind2sub([nP nP], q);
        if any(pick(:, 1) == j & pick(:, 2) == i) && max(abs(Rraw.H{i, j} - Rraw.H{j, i})) < 1e-2*max(abs(Rraw.H{i, j}))   % reciprocal within 1 %
            continue                                       % reciprocal duplicate
        end
        pick(end+1, :) = [i j]; %#ok<AGROW>
        if size(pick, 1) == n_worst, break, end
    end
    if ~isempty(pick)
        nW  = size(pick, 1);  nc = min(3, nW);
        fig = figure('Color', 'w', 'NumberTitle', 'off', 'Name', 'SnP worst causality', ltxDefaults{:});
        tl  = tiledlayout(fig, ceil(nW/nc), nc, 'TileSpacing', 'compact');
        if nViol > 0
            verdict = sprintf('%d significant entries (within %g dB of the strongest) lead t = 0 by more than one time step (%.0f ps)', nViol, sig_dB, dt*1e12);
        else
            verdict = sprintf('all entries within %g dB of the strongest causal within one time step (%.0f ps); worst non-causal energy %.4f\\%% of the strongest', sig_dB, dt*1e12, Rraw.EncMaxRel);
        end
        title(tl, {sprintf(['Worst-case significant entries of %s%s: impulse responses ' ...
            '(raw data, no filter; any filter is applied to the excitation, not to the block)'], esc(fn), fe), ...
            verdict}, 'Interpreter', 'latex');
        for p = 1:nW
            i = pick(p, 1);  j = pick(p, 2);  h = Rraw.H{i, j};  a = abs(h)/max(abs(h));
            ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
            yyaxis(ax, 'left');
            plot(ax, t*1e9, 20*log10(a + realmin), '-', 'LineWidth', 1.2);
            plot(ax, t*1e9, 20*log10(Rraw.href + realmin), 'k--', 'LineWidth', 1);
            yline(ax, onset_dB, 'k:', 'HandleVisibility', 'off');  ylim(ax, [-80 3]);  ylabel(ax, '$|h|$ (dB re peak)');
            yyaxis(ax, 'right');
            plot(ax, t*1e9, h/max(abs(h)), '-', 'LineWidth', 0.8);  ylim(ax, [-1.05 1.05]);
            ylabel(ax, '$h / \max|h|$');
            xline(ax, 0, 'k-', 'HandleVisibility', 'off');  xline(ax, Rraw.t_ref*1e9, 'k:', 'HandleVisibility', 'off');
            xline(ax, (Rraw.t_ref + Rraw.lead(i, j))*1e9, 'r-', 'onset', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
            xlim(ax, [min(-1, (Rraw.t_ref + Rraw.lead(i, j))*1e9 - 0.5), 3]);  xlabel(ax, 'Time (ns)');
            [Smax, kmax] = max(abs(squeeze(P(i, j, :))));                % largest |S_ij| over the file's band
            title(ax, {sprintf('S(%s, %s): lead %.0f ps', namesL(i), namesL(j), Rraw.lead(i, j)*1e12), ...
                sprintf('non-causal energy %.2f\\%%, peak coupling %.3g at %.3g GHz', Rraw.Enc(i, j), Smax, f(kmax)/1e9)}, 'Interpreter', 'latex');
        end
        lgd = legend(ax, {'$|h|$ (dB)', 'band-limited delta (reference)', '$h$ (linear)'}, 'NumColumns', 3);
        lgd.Layout.Tile = 'south';
    end
    
    end
end

%% ========================================================================
function w = besselResp(n, fc, fu)
% analog Bessel low-pass of order n (besself cutoff fc) on the grid fu
[b, a] = besself(n, 2*pi*fc);
w = freqs(b, a, 2*pi*fu);  w = w(:);  w(1) = real(w(1));
end

function R = evalTaper(Su, w, t, dt, onset_dB, causalTaper, sig_dB)
% Impulse response of every entry with band limiting w, and the causality
% measures: lead (onset re the band-limited delta), Enc (% of the entry's
% energy before t = 0, or before -|t_ref| for a zero-phase taper), EncRel
% (the same energy as % of the strongest entry's energy), tail level (dB),
% sig (entries within sig_dB of the strongest), and their worst leadSig,
% EncSig, EncMaxRel.
toTime = @(X) fftshift(real(ifft([X; conj(X(end-1:-1:2))])));
href = abs(toTime(w));  R.href = href / max(href);
R.t_ref = t(find(20*log10(R.href + realmin) > onset_dB, 1));
if causalTaper, iNeg = t < 0; else, iNeg = t < -abs(R.t_ref); end
iTail = t > 0.4*(numel(t)*dt);
nP = size(Su, 1);  R.lead = nan(nP);  R.Enc = nan(nP);  R.tail = nan(nP);  R.H = cell(nP);
E = nan(nP);  Eneg = nan(nP);
for q = find(~cellfun(@isempty, Su)).'
    h = toTime(Su{q} .* w);  R.H{q} = h;  a = abs(h) / max(abs(h));
    R.lead(q) = t(find(20*log10(a + realmin) > onset_dB, 1)) - R.t_ref;
    E(q) = sum(h.^2);  Eneg(q) = sum(h(iNeg).^2);
    R.Enc(q)  = 100 * Eneg(q) / E(q);
    R.tail(q) = 20*log10(max(a(iTail)) + realmin);
end
Emax = max(E(:));  R.EncRel = 100 * Eneg / Emax;  R.EncMaxRel = max(R.EncRel(:));
R.sig = E >= Emax * 10^(-sig_dB/10);
R.leadSig = min([0; R.lead(R.sig)]);
R.EncSig = max([0; R.Enc(R.sig)]);                    % worst own-energy %, significant entries
end

function c = redblue(n)
% blue (negative) - white - red (positive)
x = linspace(-1, 1, n).';
c = [min(1, 1 + x), 1 - abs(x), min(1, 1 - x)];
c = max(min(c, 1), 0);
end
