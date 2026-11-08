%TESTING_GRADYLM_VECTORS  Vector-field view of the gradYlm check.
%   Same comparison as testing_gradYlm.m (gradYlm.m vs the analytic
%   expressions from mathematica/tables/FTheta_table.txt and
%   FPhi_table.txt), shown as tangential vector fields on the unit sphere:
%   one cone glyph per grid point, coloured by magnitude.

clear
close all
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

% Coarse grid: the vector plot becomes unreadable with fine sampling.
theta=linspace(0,pi,21).';
phi=linspace(0,2*pi,21);
lmax=1;

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

        %% ---- Comparison plot: gradYlm (numeric) vs analytic (Wolfram-derived) ----
        plotGradientVectorComparison(THETA, PHI, ...
            gradYlm_th(:,:,lindx,mindx), gradYlm_ph(:,:,lindx,mindx), ...
            anlytc_gradYlm_th(:,:,lindx,mindx), anlytc_gradYlm_ph(:,:,lindx,mindx), ...
            lval(lindx), m_here);
    end
end

%% ================= Local functions ================= %%

function plotGradientVectorComparison(THETA, PHI, gradTh, gradPh, anlyticTh, anlyticPh, lVal, mVal)
%PLOTGRADIENTVECTORCOMPARISON 3D vector field plots (quiver3) comparing
%   the numeric gradYlm output against the analytic (Wolfram-derived)
%   anlytc_gradYlm_th/anlytc_gradYlm_ph output, for a single (l,m) pair.
%
%   Each vector is drawn at its (theta,phi) position on the unit sphere,
%   with direction/length given by converting the (theta-component,
%   phi-component) tangential vector to Cartesian (see
%   sphVecToCartVec). Two side-by-side subplots: gradYlm (numeric) vs
%   the analytic result.

    % Sanitize Inf -> NaN for plotting only (phi-component is
    % mathematically Inf at the poles for m~=0; this is expected, see
    % gradYlm.m's own comments -- just drop those arrows).
    gradPhPlot    = gradPh;    gradPhPlot(~isfinite(gradPhPlot))       = NaN;
    anlyticPhPlot = anlyticPh; anlyticPhPlot(~isfinite(anlyticPhPlot)) = NaN;

    figure('Name', sprintf('Gradient vector field comparison, l=%d, m=%d', lVal, mVal));

    subplot(1,2,1);
    plotVectorFieldOnSphere(THETA, PHI, gradTh, gradPhPlot, 'gradYlm (numeric)');

    subplot(1,2,2);
    plotVectorFieldOnSphere(THETA, PHI, anlyticTh, anlyticPhPlot, '$\nabla Y_{lm}$ (analytic)');

    sgtitle(sprintf('Gradient vector field, l = %d, m = %d', lVal, mVal));
end

function plotVectorFieldOnSphere(THETA, PHI, Vtheta, Vphi, labelStr)
%PLOTVECTORFIELDONSPHERE Draws the tangential vector field (Vtheta,Vphi)
%   at each (THETA,PHI) grid point, positioned on the unit sphere, as
%   solid cone glyphs (arrowheads only, no connecting stem line -- see
%   plotCone). An opaque gray reference sphere is drawn first for
%   visual context.
%
%   Each cone has a fixed length and is coloured by the vector
%   magnitude.

    % Position of each grid point on the unit sphere (r=1)
    X = sin(THETA) .* cos(PHI);
    Y = sin(THETA) .* sin(PHI);
    Z = cos(THETA);

    % Convert the tangential (theta,phi) vector components to Cartesian.
    % Standard closed-form spherical-to-Cartesian unit vector
    % transformation, restricted to zero radial component (the angular
    % gradient has no radial part):
    %   Vx = Vtheta*cos(theta)*cos(phi) - Vphi*sin(phi)
    %   Vy = Vtheta*cos(theta)*sin(phi) + Vphi*cos(phi)
    %   Vz = -Vtheta*sin(theta)
    [Vx, Vy, Vz] = sphVecToCartVec(Vtheta, Vphi, THETA, PHI);

    mag = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
    magMax = max(mag(:), [], 'omitnan');
    if magMax == 0 || isnan(magMax)
        magMax = 1;   % avoid divide-by-zero, e.g. the all-zero m=0 phi-component
    end

    % Fixed cone length for every vector; direction and colour carry the
    % vector information.
    fixedConeLength = 0.15;
    baseRadiusRatio = 0.15;   % cone base radius, as a fraction of its length

    cmap = whiteYellowOrangeRed(256);

    % Grey reference sphere with latitude/longitude grid lines.
    [Xs, Ys, Zs] = sphere(20);
    surf(Xs, Ys, Zs, 'FaceColor',[0.85 0.85 0.85], 'FaceAlpha',0.9, ...
        'EdgeColor',[0.55 0.55 0.55], 'LineWidth',0.5); hold on;

    % coneplot needs a volumetric grid, so each vector is drawn as an
    % individual solid cone (plotCone), coloured by its magnitude.
    for i = 1:numel(X)
        if ~isfinite(Vx(i)) || ~isfinite(Vy(i)) || ~isfinite(Vz(i)) || mag(i) == 0
            continue
        end
        colorIdx = max(1, min(size(cmap,1), round((mag(i)/magMax) * (size(cmap,1)-1)) + 1));
        coneColor = cmap(colorIdx,:);
        plotCone([X(i) Y(i) Z(i)], [Vx(i) Vy(i) Vz(i)], fixedConeLength, ...
            fixedConeLength*baseRadiusRatio, coneColor);
    end

    colormap(gca, cmap);
    clim([0, magMax]);
    colorbar;

    axis equal; view(45,30); grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    % LaTeX interpreter is needed for \nabla.
    title(labelStr, 'Interpreter','latex');
end

function [Vx, Vy, Vz] = sphVecToCartVec(Vtheta, Vphi, THETA, PHI)
%SPHVECTOCARTVEC Converts a tangential spherical vector (theta- and
%   phi-components, zero radial component) to Cartesian components, at
%   each point specified by THETA, PHI.
    Vx = Vtheta .* cos(THETA) .* cos(PHI) - Vphi .* sin(PHI);
    Vy = Vtheta .* cos(THETA) .* sin(PHI) + Vphi .* cos(PHI);
    Vz = -Vtheta .* sin(THETA);
end

function plotCone(origin, direction, coneLength, baseRadius, faceColor)
%PLOTCONE Draws a single solid 3D cone glyph: base at 'origin', pointing
%   along 'direction' (any nonzero length; only its direction is used),
%   extending for 'coneLength', tapering from 'baseRadius' at the base
%   to a point at the tip. Used as a standalone "arrowhead" glyph, with
%   no connecting stem line.
%
%   Construction: build a unit cone aligned with +z (via cylinder, base
%   radius tapering to 0), then rotate it to align with 'direction'
%   using Rodrigues' rotation formula (the standard closed-form
%   solution for the rotation matrix that maps one unit vector onto
%   another), then translate to 'origin'.

    nsides = 12;
    [cx, cy, cz] = cylinder([baseRadius 0], nsides);   % unit-height cone along +z, base at z=0
    cz = cz * coneLength;

    dirUnit = direction(:) / norm(direction);
    zaxis = [0; 0; 1];
    v = cross(zaxis, dirUnit);
    s = norm(v);
    c = dot(zaxis, dirUnit);

    if s < 1e-10
        if c > 0
            R = eye(3);                % direction already along +z
        else
            R = diag([1, -1, -1]);     % direction along -z: 180-degree flip
        end
    else
        vx = [   0    -v(3)   v(2);
               v(3)     0    -v(1);
              -v(2)    v(1)     0  ];
        R = eye(3) + vx + vx^2 * ((1 - c) / s^2);   % Rodrigues' rotation formula
    end

    pts = R * [cx(:)'; cy(:)'; cz(:)'];
    cxr = reshape(pts(1,:), size(cx)) + origin(1);
    cyr = reshape(pts(2,:), size(cy)) + origin(2);
    czr = reshape(pts(3,:), size(cz)) + origin(3);

    surf(cxr, cyr, czr, 'FaceColor', faceColor, 'EdgeColor', 'none');
end

function cmap = whiteYellowOrangeRed(n)
%WHITEYELLOWORANGERED Sequential colormap: white (small) -> yellow ->
%   orange -> red (large), so small values fade into the background.
    if nargin < 1
        n = 256;
    end
    keypoints = [1   1   1;    % white
                 1   1   0;    % yellow
                 1   0.5 0;    % orange
                 1   0   0];   % red
    nk = size(keypoints, 1);
    segLen = floor(n / (nk-1));
    cmap = [];
    for k = 1:nk-1
        seg = [linspace(keypoints(k,1), keypoints(k+1,1), segLen)', ...
               linspace(keypoints(k,2), keypoints(k+1,2), segLen)', ...
               linspace(keypoints(k,3), keypoints(k+1,3), segLen)'];
        cmap = [cmap; seg]; %#ok<AGROW>
    end
    % Pad/trim to exactly n rows in case of rounding from the floor().
    if size(cmap,1) < n
        cmap = [cmap; repmat(cmap(end,:), n-size(cmap,1), 1)];
    elseif size(cmap,1) > n
        cmap = cmap(1:n,:);
    end
end