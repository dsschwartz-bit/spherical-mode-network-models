%CHECK_LOGDERIV  Riccati-Bessel logarithmic derivative vs its asymptote.
%   Compares the exact logarithmic derivative J_hat'_nu(ka)/J_hat_nu(ka)
%   with the large-order asymptote (nu+1)/ka over a range of nu at fixed
%   ka, and reports the order at which J_hat_nu underflows
%   (|J_hat_nu| < eps*max|J_hat_nu|) together with the asymptote error there.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

%% --- Parameters ---
ka    = 3.48;          % electrical size at which to evaluate
nu_vec = 0.5:0.5:50; % orders; range extends far enough to reach underflow

%% --- Exact logarithmic derivative ---
J_hat       = sphericalRicatti_j(nu_vec, ka);
J_hat_nm1   = sphericalRicatti_j(nu_vec-1, ka);
j_nu        = sphericalBesselj(nu_vec, ka);

% d/dz [z j_nu(z)] = J_hat_{nu-1}(z) - nu*j_nu(z)
dJ_hat      = J_hat_nm1 - nu_vec .* j_nu;

logderiv_exact = dJ_hat ./ J_hat;

%% --- Asymptote ---
logderiv_asymp = (nu_vec + 1) / ka;

%% --- Ratio j_{nu-1}/j_nu via recurrence (avoids underflow) ---
% j_{nu-1}/j_nu = (2*nu+1)/ka - j_{nu+1}/j_nu
% For large nu/ka, j_{nu+1}/j_nu -> 0, so j_{nu-1}/j_nu -> (2*nu+1)/ka
j_ratio_recurrence = (2*nu_vec + 1)/ka;  % leading term only
logderiv_recurrence = j_ratio_recurrence - nu_vec/ka;  % = (nu+1)/ka, the asymptote

%% --- Plot ---
figure()
plot(nu_vec, real(logderiv_exact),  'b-',  'LineWidth', 2, 'DisplayName', 'Exact')
hold on
plot(nu_vec, logderiv_asymp,        'r--', 'LineWidth', 1.5, 'DisplayName', '$(\nu+1)/ka$')
hold off
xlabel('$\nu$', 'Interpreter', 'latex')
ylabel('$\mathrm{Re}(\hat{J}''_\nu / \hat{J}_\nu)$', 'Interpreter', 'latex')
title(sprintf('Logarithmic derivative vs asymptote  ($ka = %.3g$)', ka), 'Interpreter', 'latex')
legend('Interpreter', 'latex')
grid on

thresh = eps * max(abs(J_hat));
figure()
semilogy(nu_vec, abs(J_hat),  'b-',  'LineWidth', 2, 'DisplayName', '$|\hat{J}_\nu(ka)|$')
hold on
semilogy(nu_vec, abs(dJ_hat), 'r--', 'LineWidth', 1.5, 'DisplayName', '$|\hat{J}''_\nu(ka)|$')
yline(thresh, 'm-', 'LineWidth', 1.5)
text(nu_vec(end)*0.5, thresh*3, sprintf('$\\epsilon \\cdot \\max|\\hat{J}_\\nu|$ = %.2e', thresh), ...
    'Interpreter', 'latex', 'FontSize', 11, 'Color', 'm')
hold off
xlabel('$\nu$', 'Interpreter', 'latex')
ylabel('magnitude', 'Interpreter', 'latex')
title('Bessel function magnitudes (underflow check)', 'Interpreter', 'latex')
legend('Interpreter', 'latex')
grid on

%% --- Print where exact underflows and asymptote takes over ---
underflow_idx = find(abs(J_hat) < eps * max(abs(J_hat)));
if ~isempty(underflow_idx)
    fprintf('J_hat underflows (< eps*max) for nu >= %.3g\n', nu_vec(underflow_idx(1)))
    fprintf('Asymptote error at that nu: %.4g%%\n', ...
        abs(logderiv_asymp(underflow_idx(1)) - real(logderiv_exact(underflow_idx(1)-1))) ...
        / abs(real(logderiv_exact(underflow_idx(1)-1))) * 100)
else
    fprintf('No underflow detected over tested nu range\n')
end