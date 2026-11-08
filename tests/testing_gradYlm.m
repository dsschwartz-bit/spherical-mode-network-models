%TESTING_GRADYLM  Check gradYlm against the Mathematica-derived tables.
%   For l = lmax and all m, compares the theta and phi components of the
%   angular gradient from gradYlm.m with analytic expressions built from
%   mathematica/tables/FTheta_table.txt and FPhi_table.txt (exported by
%   mathematica/LegendreFThetaFPhi.wl). Plots 3-D patterns and angular cuts
%   (theta-cuts at phi = 0, pi/2; phi-cut at theta = pi/2) for each (l,m).

clear
close all
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/
%%
theta=linspace(0,pi,75).';
phi=linspace(0,2*pi,75);
lmax=2;

numth = length(theta);
numph = length(phi);

lval = (lmax).';   % single mode order l = lmax
mval = (-lmax:lmax);   % m ranges from -lmax to +lmax
L    = length(lval);
M    = length(mval);

% Preallocate angular basis functions
Ylm = zeros(numth,numph,L,M); gradYlm_th=Ylm; gradYlm_ph=Ylm;
anlytc_gradYlm_th = gradYlm_th; anlytc_gradYlm_ph = gradYlm_ph;

% Pre-compute l-dependent quantities
l_factor = lval .* (lval + 1);   % (L x 1) — l(l+1) for each mode order

have_CSphase = false;

% Build 2D theta/phi grid for spherical harmonic evaluation
[PHI, THETA] = meshgrid(phi, theta);   % both (numth x numph)

%% ---- Load the Mathematica-exported Legendre tables as function handles ----
% The tables must cover l = 1..lmax (set lValues in
% mathematica/LegendreFThetaFPhi.wl and re-export); otherwise EvalLegendreF
% errors with "No table entry found".
dataDir = fullfile(repo_root, 'mathematica', 'tables');
[lListW, mListW, FThetaFun, FPhiFun] = LoadLegendreFHandles(dataDir);

for lindx = 1:L
    for mindx = 1:M

        if abs(mval(mindx)) > lval(lindx)
            continue
        end

        % Normalised angular basis: Ylm = Y_l_m(theta,phi) / sqrt(l(l+1))
        % tesseral harmonic (real-valued spherical harmonics), power normalized
        Ylm(:,:,lindx,mindx) = harmonicY(lval(lindx),mval(mindx), THETA, PHI, ...
            'type','real','norm',true,'phase',have_CSphase) * sqrt(2^(~(mval(mindx)==0))) / sqrt(l_factor(lindx));
        [gradYlm_th(:,:,lindx,mindx),gradYlm_ph(:,:,lindx,mindx)] = gradYlm(lval(lindx),mval(mindx),THETA,PHI,have_CSphase) ;

        % ---- Analytic (Wolfram-derived) cross-check of gradYlm_th/gradYlm_ph ----
        % F_theta = -sin(theta)*d(P_bar_l^|m|)/dx |_{x=cos(theta)}
        %           (the minus sign comes from the chain rule,
        %           d/dtheta = d/d(cos theta) * (-sin theta), and is
        %           baked directly into the Wolfram export, so no
        %           extra explicit "-" is needed in the formulas below)
        % F_phi   = P_bar_l^|m|(cos theta) / sin(theta)
        % Table lookup always uses abs(m), matching the convention used
        % internally by gradYlm.m / legendre_derivative_modified.m for
        % tesseral harmonics (Legendre part depends only on |m|; the
        % sign of m only selects/differentiates the azimuthal factor).
        m_here = mval(mindx);
        mm     = abs(m_here);

        FthetaVal = EvalLegendreF(FThetaFun, lListW, mListW, lval(lindx), mm, THETA);
        FphiVal   = EvalLegendreF(FPhiFun,   lListW, mListW, lval(lindx), mm, THETA);

        % Normalization split into two factors:
        %   Nl      = power normalization only (depends on l alone),
        %             1/sqrt(l(l+1))
        %   PhiNorm = azimuthal eigenfunction normalization,
        %             1/sqrt(2*pi) for m=0, 1/sqrt(pi) for m~=0,
        %             so that Phi_m(phi) itself satisfies
        %             int_0^2pi Phi_m(phi)^2 dphi = 1 on its own.
        Nl = 1 / sqrt(l_factor(lindx));
        if m_here == 0
            PhiNorm = 1 / sqrt(2*pi);
        else
            PhiNorm = 1 / sqrt(pi);
        end

        if m_here == 0
            anlytc_gradYlm_th(:,:,lindx,mindx) = Nl * PhiNorm * FthetaVal;
            anlytc_gradYlm_ph(:,:,lindx,mindx) = zeros(size(THETA));
        elseif m_here < 0
            anlytc_gradYlm_th(:,:,lindx,mindx) = Nl * PhiNorm * FthetaVal .* sin(mm*PHI);
            anlytc_gradYlm_ph(:,:,lindx,mindx) = Nl * PhiNorm * (mm * cos(mm*PHI)) .* FphiVal;
        else
            anlytc_gradYlm_th(:,:,lindx,mindx) = Nl * PhiNorm * FthetaVal .* cos(mm*PHI);
            anlytc_gradYlm_ph(:,:,lindx,mindx) = -Nl * PhiNorm * (mm * sin(mm*PHI)) .* FphiVal;
        end

        %% ---- Comparison plots: gradYlm (numeric) vs analytic (Wolfram-derived) ----
        plotGradientComparison(theta, phi, THETA, PHI, ...
            gradYlm_th(:,:,lindx,mindx), gradYlm_ph(:,:,lindx,mindx), ...
            anlytc_gradYlm_th(:,:,lindx,mindx), anlytc_gradYlm_ph(:,:,lindx,mindx), ...
            lval(lindx), m_here);
    end
end

%% ================= Local functions ================= %%

function plotGradientComparison(theta, phi, THETA, PHI, ...
        gradTh, gradPh, anlyticTh, anlyticPh, lVal, mVal)
%PLOTGRADIENTCOMPARISON 3D surface plots + rectangular angular cuts
%   comparing the numeric gradYlm output against the analytic
%   (Wolfram-derived) anlytc_gradYlm_th/anlytc_gradYlm_ph output, for a
%   single (l,m) pair.
%
%   Cuts produced: theta-cuts at phi=0 and phi=pi/2 combined into one
%   plot per component (with legend distinguishing phi value and
%   method), and a phi-cut at theta=pi/2, for both the theta- and
%   phi-components.

    % Sanitize Inf -> NaN for plotting only (phi-component is
    % mathematically Inf at the poles for m~=0; this is expected, see
    % gradYlm.m's own comments -- just leave a gap in the plots).
    gradPhPlot    = gradPh;    gradPhPlot(~isfinite(gradPhPlot))       = NaN;
    anlyticPhPlot = anlyticPh; anlyticPhPlot(~isfinite(anlyticPhPlot)) = NaN;

    %% ---- 3D spherical (radiation-pattern-style) comparison ----
    % Radius = magnitude of the function value; angular position comes
    % from (theta,phi); color encodes the SIGNED value (diverging
    % colormap centered at zero), so sign information isn't lost by
    % taking the magnitude for the radius.
    figure('Name', sprintf('Gradient comparison (3D), l=%d, m=%d', lVal, mVal));

    subplot(2,2,1);
    plotSphericalMagnitude(THETA, PHI, gradTh, 'gradYlm, $\theta$-component');

    subplot(2,2,2);
    plotSphericalMagnitude(THETA, PHI, anlyticTh, '\nabla Y_{lm} (analytic), $\theta$-component');

    subplot(2,2,3);
    plotSphericalMagnitude(THETA, PHI, gradPhPlot, 'gradYlm, $\phi$-component');

    subplot(2,2,4);
    plotSphericalMagnitude(THETA, PHI, anlyticPhPlot, '\nabla Y_{lm} (analytic), $\phi$-component');

    sgtitle(sprintf('Gradient comparison (3D), l = %d, m = %d', lVal, mVal));

    %% ---- Rectangular (Cartesian) cuts ----
    % Cuts are interpolated rather than taken at the nearest grid node, so
    % they sit exactly at phi = 0, pi/2 and theta = pi/2 for any grid
    % resolution (values that should be exactly zero then stay zero).
    gradTh_phi0        = interp1(phi, gradTh.',        0,    'linear').';
    gradTh_phihalfpi   = interp1(phi, gradTh.',        pi/2, 'linear').';
    anlyticTh_phi0      = interp1(phi, anlyticTh.',      0,    'linear').';
    anlyticTh_phihalfpi = interp1(phi, anlyticTh.',      pi/2, 'linear').';

    gradPh_phi0        = interp1(phi, gradPhPlot.',    0,    'linear').';
    gradPh_phihalfpi   = interp1(phi, gradPhPlot.',    pi/2, 'linear').';
    anlyticPh_phi0       = interp1(phi, anlyticPhPlot.', 0,    'linear').';
    anlyticPh_phihalfpi  = interp1(phi, anlyticPhPlot.', pi/2, 'linear').';

    gradTh_thetahalfpi    = interp1(theta, gradTh,        pi/2, 'linear');
    anlyticTh_thetahalfpi = interp1(theta, anlyticTh,      pi/2, 'linear');
    gradPh_thetahalfpi    = interp1(theta, gradPhPlot,    pi/2, 'linear');
    anlyticPh_thetahalfpi = interp1(theta, anlyticPhPlot, pi/2, 'linear');

    figure('Name', sprintf('Angular cuts, l=%d, m=%d', lVal, mVal));

    % --- theta-cuts (phi=0 and phi=pi/2 combined) ---
    subplot(2,2,1);
    plot(theta, gradTh_phi0,      'b-',  'LineWidth',2.5); hold on;
    plot(theta, anlyticTh_phi0,         'r-',  'LineWidth',1.5);
    plot(theta, gradTh_phihalfpi, 'b-.', 'LineWidth',2.5);
    plot(theta, anlyticTh_phihalfpi,    'r-.', 'LineWidth',1.5);
    yline(0, 'k:'); grid on;
    xlabel('$\theta$ (rad)'); title('$\theta$-component (theta-cuts)');
    legend('$gradYlm_{th}$, $\phi$=0$','$\nabla Y_{lm}\cdot\hat{$\theta$}$, $$\phi$=0$', ...
        '$gradYlm_{th}$, $\phi=\pi/2$','$\nabla Y_{lm}\cdot\hat{\theta}$, $\phi=\pi/2$', ...
        'Location','best', 'Interpreter','latex');

    subplot(2,2,3);
    plot(theta, gradPh_phi0,      'b-',  'LineWidth',2.5); hold on;
    plot(theta, anlyticPh_phi0,         'r-',  'LineWidth',1.5);
    plot(theta, gradPh_phihalfpi, 'b-.', 'LineWidth',2.5);
    plot(theta, anlyticPh_phihalfpi,    'r-.', 'LineWidth',1.5);
    yline(0, 'k:'); grid on;
    xlabel('$\theta$ (rad)','Interpreter','latex'); title('$\phi$-component (theta-cuts)','Interpreter','latex');
    legend('$gradYlm_{ph}$, $\phi=0$','$\nabla Y_{lm}\cdot\hat{$\phi$}$, $\phi=0$', ...
        '$gradYlm_{ph}$, $\phi$=\pi/2$','$\nabla Y_{lm}\cdot\hat{\phi}$, $\phi=\pi/2$', ...
        'Location','best', 'Interpreter','latex');

    % --- phi-cut at theta = pi/2 ---
    subplot(2,2,2);
    plot(phi, gradTh_thetahalfpi, 'b-',  'LineWidth',2.5); hold on;
    plot(phi, anlyticTh_thetahalfpi,    'r-.', 'LineWidth',1.5);
    yline(0, 'k:'); grid on;
    xlabel('$\phi$ (rad)','Interpreter','latex'); title('$\theta$-component, $\theta$=\pi/2','Interpreter','latex');
    legend('$gradYlm_{th}$','$\nabla Y_{lm}\cdot\hat{$\theta$}$', ...
        'Location','best', 'Interpreter','latex');

    subplot(2,2,4);
    plot(phi, gradPh_thetahalfpi, 'b-',  'LineWidth',2.5); hold on;
    plot(phi, anlyticPh_thetahalfpi,    'r-.', 'LineWidth',1.5);
    yline(0, 'k:'); grid on;
    xlabel('$\phi$ (rad)','Interpreter','latex'); title('$\phi$-component, $\theta$=\pi/2','Interpreter','latex');
    legend('$gradYlm_{ph}$','$\nabla Y_{lm}\cdot\hat{$\phi$}$', ...
        'Location','best', 'Interpreter','latex');

    sgtitle(sprintf('Angular cuts, l = %d, m = %d', lVal, mVal));
end

function plotSphericalMagnitude(THETA, PHI, value, labelStr)
%PLOTSPHERICALMAGNITUDE Radiation-pattern-style plot: radius = |value|,
%   angular position from (THETA,PHI), color = signed value (diverging
%   colormap centered at zero so the sign is not lost).
    
    r = abs(value);

    X = r .* sin(THETA) .* cos(PHI);
    Y = r .* sin(THETA) .* sin(PHI);
    Z = r .* cos(THETA);

    surf(X, Y, Z, value, 'EdgeColor','none');
    colormap(gca, blueWhiteRed(256));

    cmax = max(abs(value(:)), [], 'omitnan');
    if cmax > 0
        clim([-cmax, cmax]);
    end
    colorbar;

    axis equal; view(45,30); grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title(labelStr);
end

function cmap = blueWhiteRed(n)
%BLUEWHITERED Simple diverging colormap: blue (negative) -> white
%   (zero) -> red (positive). Used so positive/negative regions of the
%   plotted value are visually distinguishable at a glance.
    if nargin < 1
        n = 256;
    end
    half = floor(n/2);
    blueToWhite = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    whiteToRed  = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
    cmap = [blueToWhite; whiteToRed];
end