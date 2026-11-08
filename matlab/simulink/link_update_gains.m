function link_update_gains()
%LINK_UPDATE_GAINS  Recompute the link gains and propagation delay in the base workspace.
%   Updates the derived variables of the Simulink model built by
%   build_tx_rx_link_simulink.m after r, theta, phi, rot_angle, Ux/Uy/Uz,
%   a_Tx, a_Rx, mu0 or c0 have changed, so the model needs no rebuild:
%
%       r = 5;  theta = 60*pi/180;   % edit any of these in the base workspace
%       link_update_gains;           % refresh Gtx_prime, Rrot, Grx, prop_t, Lind
%       sim('<model name>');
%
%   Reads : link_cfg (mode lists / ka_factor flags, set by the builder),
%           r, theta, phi, rot_angle, Ux, Uy, Uz, a_Tx, a_Rx, c0, mu0
%   Writes (the ADS VAR-block variables, as base-workspace variables):
%           theta_prime, phi_prime, R11..R33, magU,
%           gradYlm_theta_prime_l<l>_m<m>, gradYlm_phi_prime_l<l>_m<m>   (Tx modes)
%           gradYlm_theta_l<l>_m<m>,       gradYlm_phi_l<l>_m<m>         (Rx modes)
%           (negative m is spelled mneg<|m|>, as in ADS), prop_t, Lind,
%           plus the matrices Gtx_prime (2 x nTx), Rrot (2 x 2), Grx (nRx x 2)
%           assembled from them (kept for convenience / plotting).
%
%   These mirror the ADS cells (python/ads/create_vcvs_tx_gain.py, create_cccs_rx_aperture.py,
%   create_vector_field_rotation_vars.py):
%     Tx  TM: [-gth'; -gph']      TE: [-gph'; +gth']   (gradYlm at body-frame angles)
%     Rx  TM: par(l)*4*pi*[-gth, -gph]   TE: -par(l)*4*pi*[-gph, +gth]
%         (par(l) = +1 odd l, -1 even l; TE carries the extra polarity flip
%          the ADS CCCS gets from reversing its output terminals)
%     Rrot = Cart2Sph(theta,phi) * R * Sph2Cart(theta',phi')

here = fileparts(mfilename('fullpath'));
if exist('repo_root', 'file') ~= 2, run(fullfile(here, '..', '..', 'setup.m')); end   % adds matlab/ (sphmodes: gradYlm)
assert(exist('gradYlm', 'file') == 2, ...
    'gradYlm.m not found -- run setup.m in the repository root.');

cfg = evalin('base', 'link_cfg');
% Spherical-harmonic gradients: read the same Wolfram-exported symbolic expressions ADS uses
% (FTheta_ADS_Discrete.txt / FPhi_ADS_Discrete.txt in mathematica/tables, see create_gradYlm_vars.py); the
% base-workspace variable Ylm_data_dir overrides the folder.  gradYlm.m (numeric) is used instead when the
% files are missing or an order l is not in them.
ddir = fullfile(repo_root, 'mathematica', 'tables');
if evalin('base', 'exist(''Ylm_data_dir'', ''var'')'), ddir = evalin('base', 'Ylm_data_dir'); end
tab = ads_ylm_tables(ddir);
v = @(name) evalin('base', name);
r = v('r'); theta = v('theta'); phi = v('phi'); rot_angle = v('rot_angle');
U = [v('Ux') v('Uy') v('Uz')];
c0 = v('c0'); mu0 = v('mu0'); a_Tx = v('a_Tx'); a_Rx = v('a_Rx');

% ---- rotation (Rodrigues), body-frame sample direction ----
u = U(:) / norm(U);
K = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
R = eye(3) + sin(rot_angle) * K + (1 - cos(rot_angle)) * (K * K);
% poles: gradYlm locates theta = 0 / pi by exact equality, so snap near-polar angles to the
% exact pole (a rounding residue such as sin(pi) = 1.2e-16 would otherwise give 0*Inf noise)
poletol = 1e-9;
if abs(sin(theta)) < poletol, theta = round(theta/pi)*pi; end
n  = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
nb = R.' * n;
that = @(t, p) [cos(t)*cos(p); cos(t)*sin(p); -sin(t)];
phat = @(t, p) [-sin(p); cos(p); 0];
if hypot(nb(1), nb(2)) < poletol
    % body-frame direction is a pole: phi' is undefined, so take the azimuth the observation
    % direction's own theta-hat arrives with (keeps Rrot a proper, continuous rotation)
    th_p = (nb(3) < 0) * pi;
    tb = R.' * that(theta, phi);
    if nb(3) >= 0, ph_p = atan2(tb(2), tb(1)); else, ph_p = atan2(-tb(2), -tb(1)); end
else
    th_p = acos(max(-1, min(1, nb(3))));
    ph_p = atan2(nb(2), nb(1));
end
Rrot = [that(theta, phi).'; phat(theta, phi).'] * R * [that(th_p, ph_p), phat(th_p, ph_p)];

setv = @(name, val) assignin('base', name, val);
setv('magU', norm(U));
setv('theta_prime', th_p);  setv('phi_prime', ph_p);
for ii = 1:3, for jj = 1:3, setv(sprintf('R%d%d', ii, jj), R(ii, jj)); end, end
vn = @(prefix, l, m) sprintf('%s_l%d_%s', prefix, l, ternary(m < 0, sprintf('mneg%d', abs(m)), sprintf('m%d', m)));

% ---- Tx gains (body-frame angles) ----
nT = numel(cfg.tx_modes);
Gtx_prime = zeros(2, nT);
for i = 1:nT
    md = cfg.tx_modes(i);
    [gt, gp] = ylm_grad(md.l, md.m, th_p, ph_p, tab, setv, '_prime');
    setv(vn('gradYlm_theta_prime', md.l, md.m), gt);
    setv(vn('gradYlm_phi_prime',   md.l, md.m), gp);
    if strcmp(md.type, 'TM'), Gtx_prime(:, i) = [-gt; -gp];
    else,                     Gtx_prime(:, i) = [-gp;  gt]; end
end

% ---- Rx gains (unprimed angles) ----
nR = numel(cfg.rx_modes);
Grx = zeros(nR, 2);
for i = 1:nR
    md = cfg.rx_modes(i);
    [gt, gp] = ylm_grad(md.l, md.m, theta, phi, tab, setv, '');
    setv(vn('gradYlm_theta', md.l, md.m), gt);
    setv(vn('gradYlm_phi',   md.l, md.m), gp);
    par = 1; if mod(md.l, 2) == 0, par = -1; end
    if strcmp(md.type, 'TM'), Grx(i, :) = par * 4*pi * [-gt, -gp];
    else,                     Grx(i, :) = -par * 4*pi * [-gp, gt]; end
end

% ---- delay: whole propagation delay on the Tx side (transient causality) ----
corr = 0;
if cfg.tx_needs_corr, corr = corr + a_Tx; end
if cfg.rx_needs_corr, corr = corr + a_Rx; end
prop_t = (r - corr) / c0;
if prop_t < 0
    error('link_update_gains:NegativeDelay', ...
        'Propagation delay (r - a_Tx - a_Rx)/c0 is negative (r too small).');
end
Lind = mu0 * r;     % mu0 in H/m

for nm = {'Gtx_prime', 'Rrot', 'Grx', 'prop_t', 'Lind'}
    assignin('base', nm{1}, eval(nm{1}));
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function [gt, gp] = ylm_grad(l, m, th, ph, tab, setv, sfx)
% gradYlm_theta / gradYlm_phi exactly as the ADS VAR formulas (create_gradYlm_vars[_prime].py):
%   Nl = 1/sqrt(l(l+1)),  PhiNorm = 1/sqrt(2 pi) (m = 0) or 1/sqrt(pi),  and the library entries
%   Ftheta_l<l>_m<|m|>(theta), Fphi_l<l>_m<|m|>(theta) (always the positive-m entry)
mm = abs(m);
idx = l*l + l + mm + 1;                           % table index l^2 + l + m (+1 for MATLAB)
% numeric fallback (gradYlm.m): expression files not found, or this order l is beyond the file's lMax
if isempty(tab) || idx > numel(tab.ft) || isempty(tab.ft{idx}) || (m ~= 0 && (idx > numel(tab.fp) || isempty(tab.fp{idx})))
    if ~isempty(tab)
        warning('off', 'backtrace');
        warning('link_update_gains:YlmOrder', 'l = %d is not in the ADS expression files (lMax = %d): using gradYlm.m for it.', l, floor(sqrt(numel(tab.ft)) - 1));
        warning('on', 'backtrace');
    end
    [gt, gp] = gradYlm(l, m, th, ph);
    return
end
nl = 1/sqrt(l*(l+1));
if m == 0, pn = 1/sqrt(2*pi); else, pn = 1/sqrt(pi); end
Ft = tab.ft{idx}(th);
setv(sprintf('Ftheta%s_l%d_m%d', sfx, l, mm), Ft);          % the ADS Ftheta_Fphi_vars entries
if m == 0
    gt = nl*pn*Ft;  gp = 0;
else
    Fp = tab.fp{idx}(th);
    setv(sprintf('Fphi%s_l%d_m%d', sfx, l, mm), Fp);
    if m < 0
        gt = nl*pn*Ft*sin(mm*ph);   gp = nl*pn*mm*cos(mm*ph)*Fp;
    else
        gt = nl*pn*Ft*cos(mm*ph);   gp = -nl*pn*mm*sin(mm*ph)*Fp;
    end
end
end

function tab = ads_ylm_tables(ddir)
% parse the ADS DSCRDATA files into function handles of theta (cached between calls)
persistent cache
if ~isempty(cache) && strcmp(cache.dir, ddir), tab = cache.tab; return, end
tab = [];
fT = fullfile(ddir, 'FTheta_ADS_Discrete.txt');  fP = fullfile(ddir, 'FPhi_ADS_Discrete.txt');
if ~isfile(fT) || ~isfile(fP)
    warning('link_update_gains:NoYlmFiles', 'ADS spherical-harmonic files not found in %s: using gradYlm.m.', ddir);
    cache = struct('dir', ddir, 'tab', tab);  return
end
tab = struct('ft', {parse_dscr(fT)}, 'fp', {parse_dscr(fP)}, 'dir', ddir);
cache = struct('dir', ddir, 'tab', tab);
end

function f = parse_dscr(fn)
% 'index formula' lines between BEGIN/END DSCRDATA; ADS expression -> MATLAB anonymous function of theta
L = splitlines(fileread(fn));
f = {};
inblk = false;
for k = 1:numel(L)
    ln = strtrim(L{k});
    if isempty(ln) || startsWith(ln, 'REM') || startsWith(ln, '%'), continue, end
    if startsWith(upper(ln), 'BEGIN'), inblk = true; continue, end
    if startsWith(upper(ln), 'END'), inblk = false; continue, end
    if ~inblk, continue, end
    t = regexp(ln, '^(\d+)\s+(.*)$', 'tokens', 'once');
    expr = strrep(t{2}, '**', '^');
    f{str2double(t{1}) + 1} = str2func(['@(theta) ' expr]); %#ok<AGROW>
end
end
