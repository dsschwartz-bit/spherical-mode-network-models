function [cl_conv, nu_n, dnu_n, nu_approx_idx, dnu_approx_idx] = Bicon_Coeff_Symbolic_2(a, theta0_sym, Freq, lmax, nu_max,max_percent_error, conversion)
%BICON_COEFF_SYMBOLIC_2  Biconical antenna TM mode coefficients with exact conical eigenvalues.
%   [cl_conv, nu_n, dnu_n, nu_approx_idx, dnu_approx_idx] = Bicon_Coeff_Symbolic_2(a,
%       theta0_sym, Freq, lmax, nu_max, max_percent_error, conversion)
%   Same mode-matching solution as Bicon_Coeff, but the conical eigenvalues
%   nu_n (roots of coneHarmonic(nu,theta0,"odd") = 0) are found with
%   vpasolve (Symbolic Math Toolbox), one root per unit interval
%   [k, k+1], k = 0..floor(nu_max)-1, and d(nu_n)/d(theta0) by a finite
%   difference in theta0 (step eps^(1/3)). For each root found, the
%   conical harmonic of the exact and of the approximate nu_n
%   (Janaswamy 2022, Eq. 11-12; as in Bicon_Coeff) are compared on a
%   3-degree grid; once the worst-case difference, as a percentage of the
%   exact one, falls below max_percent_error, the approximation is used for
%   all remaining modes. A failed root search uses the approximate root in
%   that interval; if the whole search fails, all modes use the
%   approximation. A figure of the per-mode error is drawn.
%
%   Differences from Bicon_Coeff_Symbolic: N = max(round(nu_max),
%   round(Nfac*k_max a)) modes are sought; the interior log-derivative is
%   computed as (nu+1)/ka - j_(nu+1)(ka)/j_nu(ka), replaced by its
%   asymptote (nu+1)/ka when nu/ka > 10, with messages reporting the
%   frequencies where the asymptote was used; G is always solved with
%   backslash (a message is printed if rcond(G) < eps).
%
%   The acceleration terms t2/t3 of the G series depend only on n,
%   alpha0, beta_m, beta_ell and ka, not on nu_n.
%
%   Inputs:
%     a                 arm length / exterior sphere radius (m)
%     theta0_sym        cone half-angle (rad), numeric or symbolic
%     Freq              (1 x numf) frequencies (Hz)
%     lmax              highest (odd) degree; L = (lmax+1)/2
%     nu_max            upper end of the eigenvalue search
%     max_percent_error switch-over threshold in % (default 10)
%     conversion        true: scale to the psi convention (see Bicon_Coeff);
%                       false (default): raw mode-matching coefficients
%   Outputs:
%     cl_conv        (numf x L) TM_l0 coefficients, l = 1,3,...,lmax
%     nu_n, dnu_n    (1 x N) eigenvalues and d(nu_n)/d(theta0) used
%     nu_approx_idx, dnu_approx_idx  indices where the approximation was used
%
%   See also Bicon_Coeff, coneHarmonic, Legendre_Pnu_m.

if mod(lmax, 2) == 0
    error("lmax should be odd")
end
if nargin < 6 || isempty(max_percent_error)
    max_percent_error = 10;
end
if nargin < 7 || isempty(conversion)
    conversion = false;
end

%% --- Physical Constants ---
c0   = 299792458;       % m/s
mu0  = 4*pi*1e-7;       % H/m
eps0 = 1/c0^2/mu0;      % F/m
eta0 = sqrt(mu0/eps0);  % Ohm
theta0=double(theta0_sym);

%% --- Sweep and Discretization Parameters ---
numf     = length(Freq);
freq_max = max(Freq);
kaMx     = 2*pi*a*freq_max/c0;

L     = (lmax+1)/2;
Mifac = floor((L-1)/2);
Nfac  = Mifac;
N     = max(round(nu_max), round(Nfac*kaMx));  % number of conical modes to find
n     = (1:N);                     % Integer summation indices (used in t2/t3 acceleration)

lval  = (1:2:lmax);
Normp = sqrt( lval.*(lval+1) .* 2*pi * 2./(2*lval+1) );

%% --- Antenna and Geometry Parameters ---
alp0 = 1/(0.5 - theta0/pi);   % alpha_0

%% --- Approximate Eigenvalues (Janaswamy 2022, Eq. 11-12) ---
nu_n_approx  = sqrt((n*alp0).^2 - 0.25) - 0.5;
dnu_n_approx = n.^2*alp0^3 ./ (pi*(nu_n_approx + 0.5));

%% --- Exact eigenvalues via vpasolve; switch to the approximation once it is accurate ---
syms nu_sym
dtheta     = sym(eps)^(1/3);          % finite-difference step in theta0

th_grid = (theta0 : 3*pi/180 : pi/2);   % Theta grid for error check (radians)

Nu_search     = floor(nu_max);
search_start = (0:1:Nu_search);

nu_n_exact  = zeros(1, N);
dnu_n_exact = zeros(1, N);
per_mode_err_pct  = NaN(1, N);
vpasolve_ok       = true;
first_below_order = [];
nu_approx_idx     = [];    % Indices where nu_n fell back to approximate
dnu_approx_idx    = [];    % Indices where dnu_n fell back to approximate

num_found = 0;
try
    for m_idx = 1:(length(search_start)-1)
        if num_found >= N; break; end

        lo = search_start(m_idx);
        hi = search_start(m_idx+1);

        % One root per interval; on failure use the approximate root in it.
        try
            sym_sol = vpasolve(coneHarmonic(nu_sym, theta0_sym,"odd") == 0, nu_sym, [lo, hi]);
            nu_sol  = double(sym_sol);
        catch ME
            warning(ME.identifier, '%s — using approx nu (interval [%.3g, %.3g])', ME.message, lo, hi)
            approx_in_interval = find(nu_n_approx > lo & nu_n_approx < hi, 1);
            if ~isempty(approx_in_interval)
                nu_candidate = nu_n_approx(approx_in_interval);
                % Skip if already present (duplicate check)
                if ~any(abs(nu_n_exact - nu_candidate) < 1e-6)
                    num_found                 = num_found + 1;
                    nu_n_exact(num_found)     = nu_candidate;
                    dnu_n_exact(num_found)    = dnu_n_approx(approx_in_interval);
                    nu_approx_idx(end+1)  = num_found; %#ok<AGROW>
                    dnu_approx_idx(end+1) = num_found; %#ok<AGROW>
                end
            end
            continue
        end

        if isempty(nu_sol) || mod(nu_sol, 1) == 0 || any(abs(nu_n_exact - nu_sol) < 1e-6)
            continue
        end

        num_found             = num_found + 1;
        nu_n_exact(num_found) = nu_sol;

        % d(nu)/d(theta0) by finite difference; approximation on failure.
        try
            sym_pert           = vpasolve(coneHarmonic(nu_sym, theta0_sym+dtheta,"odd") == 0, nu_sym, [lo, hi]);
            dnu_n_exact(num_found) = (double(sym_pert) - nu_sol) / dtheta;
        catch ME
            warning(ME.identifier, '%s — using approx dnu (interval [%.3g, %.3g])', ME.message, lo, hi)
            [~, idx]              = min(abs(nu_n_approx - nu_sol));
            dnu_n_exact(num_found)    = dnu_n_approx(idx);
            dnu_approx_idx(end+1) = num_found; %#ok<AGROW>
        end

        % --- Per-mode check: exact vs approximate conical harmonic ---
        % num_found counts accepted roots, so nu_n_approx(num_found) is the
        % same-n counterpart of nu_n_exact(num_found).
        if ~any(nu_approx_idx == num_found)
            Lnu_exact  = arrayfun(@(t) coneHarmonic(nu_n_exact(num_found),t,"odd"), th_grid);
            Lnu_approx = arrayfun(@(t) coneHarmonic(nu_n_approx(num_found),t,"odd"), th_grid);
            scale = max(abs(Lnu_exact));
            if scale > 0
                per_mode_err_pct(num_found) = max(abs(Lnu_exact - Lnu_approx)) / scale * 100;

                % Below threshold: use the approximation for the remaining modes.
                if per_mode_err_pct(num_found) < max_percent_error && isempty(first_below_order)
                    first_below_order = num_found;
                    remaining = (num_found+1):N;
                    nu_n_exact(remaining)  = nu_n_approx(remaining);
                    dnu_n_exact(remaining) = dnu_n_approx(remaining);
                    nu_approx_idx  = [nu_approx_idx,  remaining]; %#ok<AGROW>
                    dnu_approx_idx = [dnu_approx_idx, remaining]; %#ok<AGROW>
                    num_found = N;   % Force loop exit
                    break
                end
            end
        end
    end

    if num_found < N
        warning("Bicon_Coeff_Symbolic: only found %d of %d requested eigenvalues.", num_found, N)
        nu_n_exact   = nu_n_exact(1:num_found);
        dnu_n_exact  = dnu_n_exact(1:num_found);
        nu_n_approx  = nu_n_approx(1:num_found);
        dnu_n_approx = dnu_n_approx(1:num_found);
        per_mode_err_pct = per_mode_err_pct(1:num_found);
        N = num_found;  n = (1:N);
    end

catch ME
    warning(ME.identifier, '%s — falling back to approximate eigenvalues', ME.message)
    vpasolve_ok = false;
end

%% --- Select Eigenvalues and Report ---
if ~vpasolve_ok
    nu_n           = nu_n_approx;
    dnu_n          = dnu_n_approx;
    nu_approx_idx  = 1:N;
    dnu_approx_idx = 1:N;
else
    nu_n  = nu_n_exact;
    dnu_n = dnu_n_exact;

    if ~isempty(first_below_order)
        fprintf("Bicon_Coeff_Symbolic: approximation error drops below %.3g%% at mode n=%d (nu_n=%.6g) — using approximate nu_n / dnu_n for remaining %d modes\n", ...
                max_percent_error, first_below_order, nu_n_exact(first_below_order), N - first_below_order)
    end
end

%% --- Plot per-mode errors (if any were computed) ---
if any(~isnan(per_mode_err_pct))
    figure()
    plot(1:N, per_mode_err_pct, '-', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6)
    hold on
    yline(max_percent_error, '--r', sprintf('Threshold %.3g%%', max_percent_error), ...
          'LabelHorizontalAlignment', 'left')
    hold off
    xlabel('Mode index n')
    ylabel('Max shape error (%)', 'Interpreter', 'none')
    title(sprintf('Mnu shape error: exact vs approximate nu_n  (theta_0=%.4g rad)', theta0), 'Interpreter', 'none')
    xticks(1:N)
    grid on
end

%% --- Main Frequency Sweep Loop ---
cl_conv = zeros(numf, L);
c       = cl_conv;
asymp_last_f   = [];
asymp_last_ka  = [];
asymp_last_num = [];
prev_evanescent = false(size(nu_n));

for f = 1:numf
    ka = 2*pi*a*Freq(f)/c0;

    %% --- Bessel Function Ratios ---
    % JpJ = (x j_nu)'/(x j_nu) = (nu+1)/x - j_(nu+1)/j_nu at x = ka;
    % for nu/ka > 10 the ratio j_(nu+1)/j_nu is negligible.
    evanescent = (nu_n / ka) > 10;
    j_nu   = sphericalBesselj(nu_n, ka);
    j_nu1  = sphericalBesselj(nu_n + 1, ka);
    JpJ    = zeros(size(nu_n));
    JpJ( evanescent) = (nu_n(evanescent) + 1) / ka;
    JpJ(~evanescent) = (nu_n(~evanescent) + 1) / ka - j_nu1(~evanescent) ./ j_nu(~evanescent);

    % Report modes that switched from asymptotic to exact at this frequency.
    if f > 1
        newly_exact = ~evanescent & prev_evanescent;
        for nn = find(newly_exact)
            fprintf('Bicon_Coeff_Symbolic_2: mode n=%d (nu=%.4g) last used asymptotic JpJ=(nu+1)/ka at f=%.6g Hz\n', ...
                nn, nu_n(nn), Freq(f-1))
        end
    end
    prev_evanescent = evanescent;

    if any(evanescent)
        asymp_last_f   = Freq(f);
        asymp_last_ka  = ka;
        asymp_last_num = sum(evanescent);
    end

    %% --- Mode Coupling Prefactor ---
    Fac1 = dnu_n .* (2*nu_n+1) ./ (nu_n.*(nu_n+1));
    B    = Fac1 .* JpJ;

    %% --- Double Mode Index Loop: Build G Matrix ---
    gml = zeros(L,L);
    for mi = 1:L
        m      = 2*mi - 1;
        y      = legendre(m, cos(theta0));
        p(mi)  = y(1);

        for elli = 1:mi
            ell  = 2*elli - 1;
            bm   = sqrt(m*(m+1))/alp0;
            bell = sqrt(ell*(ell+1))/alp0;

            nu_sq = nu_n .* (nu_n + 1);
            t1 = B ./ ((nu_sq - m*(m+1)) .* (nu_sq - ell*(ell+1)));

            t2 = 1./(pi*alp0^2*(bm+bell)*ka) * ...
                 (1./((n-bm).*(n-bell)) - 1./((n+bm).*(n+bell)));

            if m ~= ell
                t3 = 1./(pi*alp0^2*(bm+bell)*ka) * ...
                     (pi*sin(pi*(bm-bell))/((bm-bell)*sin(pi*bm)*sin(pi*bell)) ...
                     + 1/(bm*bell) - 2*(psi(bm)-psi(bell))/(bm-bell));
            else
                t3 = 1./(pi*alp0^2*(bm+bell)*ka) * ...
                     ((pi/sin(pi*bm))^2 + 1/bm^2 - 2*psi(1,bm));
            end

            gml(mi, elli) = sum(t1 - t2) + t3;

            y2 = legendre(ell, cos(theta0));
            gml(mi, elli) = gml(mi, elli) * sin(theta0) * y(1) * y2(1) * m*(m+1) * ell*(ell+1);

            gml(elli, mi) = gml(mi, elli);

            if elli == mi
                H1 = besselh(m - 1/2, 2, ka);
                H2 = besselh(m + 3/2, 2, ka);
                H3 = besselh(m + 1/2, 2, ka);
                H  = 1/ka + (H1-H2)/H3;
                gml(mi, elli) = -0.5*H*m*(m+1)/(2*m+1) + gml(mi, elli);
            end

        end  % elli loop
    end  % mi loop

    r = rcond(gml);
    if r < eps
        fprintf("rcond(gml) below eps")
        c(f,:) = -1j * lsqminnorm(gml, p.');
    else
        c(f,:) = -1j * gml \ p.';
    end

    if conversion
        conv_factor  = sqrt(eta0)/(2*pi) * (1j.^(lval+1)) .* Normp ./ sphericalRicatti_h(lval, 2, ka);
        cl_conv(f,:) = c(f,:) .* conv_factor;
    else
        cl_conv(f,:) = c(f,:);
    end

end  % frequency sweep loop

if ~isempty(asymp_last_f)
    if asymp_last_f == Freq(end)
        for nn = find(prev_evanescent)
            fprintf('Bicon_Coeff_Symbolic_2: mode n=%d (nu=%.4g): asymptotic JpJ=(nu+1)/ka used at all frequencies\n', ...
                nn, nu_n(nn))
        end
    else
        fprintf('Bicon_Coeff_Symbolic_2: asymptotic JpJ = (nu+1)/ka last used at f=%.6g Hz\n', ...
            asymp_last_f)
    end
else
    fprintf('Bicon_Coeff_Symbolic_2: asymptotic JpJ = (nu+1)/ka never used — all modes have nu/ka <= 10 at all frequencies\n')
end
fprintf('\n')
end