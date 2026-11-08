%TESTING_RECURRENCE_RELATIONS  Check recurrences for rhat*Y_lm.
%   For real (tesseral) harmonics of order l = lmax and all m, rebuilds
%   x*Y_lm, y*Y_lm and z*Y_lm on the unit sphere from the l+1 / l-1
%   associated Legendre recurrence relations and compares them with the
%   direct products sin(theta)cos(phi)*Y_lm, sin(theta)sin(phi)*Y_lm and
%   cos(theta)*Y_lm, one figure per (l,m).

clear
close all
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

%% setup

numth=20;numph=20;


theta=linspace(0,pi,numth);
phi=linspace(0,2*pi,numph);
[PHI,THETA]=meshgrid(phi,theta);

lmax = 4;
lval = lmax;   % single mode order l = lmax
mval = (-lmax:lmax);   % m ranges from -lmax to +lmax
L=length(lval);M=length(mval);

Ylm=zeros(numth,numph,L,M);
factor_rho1=Ylm;factor_rho2=Ylm;factor_z1=Ylm;factor_z2=Ylm;
Fx_lm=Ylm;Fy_lm=Ylm;Fz_lm=Ylm;

for lindx=1:L
    for mindx=1:M
        l = lval(lindx);
        m=mval(mindx);
        Ylm(:,:,lindx,mindx) = TesseralHarmonic(l,m,THETA,PHI);

         
        factor_rho1(:,:,lindx,mindx) = (-1)^(m==0)...
                                *sqrt( ((2*l+1)/(2*l+3)) * (l-abs(m)+1)*(l-abs(m)+2) )...
                                *Tesseral_Plm(l+1,abs(m)-1,THETA);

        factor_rho2(:,:,lindx,mindx) = (-1)^(m==0)...
                                *sqrt( ((2*l+1)/(2*l-1)) * (l+abs(m))*(l+abs(m)-1) )...
                                *Tesseral_Plm(l-1,abs(m)-1,THETA);
        

        factor_z1(:,:,lindx,mindx) = sqrt( ((2*l+1)/(2*l+3)) * (l+abs(m)+1)*(l-abs(m)+1) )...
                                        *Tesseral_Plm(l+1,m,THETA);

        factor_z2(:,:,lindx,mindx) = sqrt( ((2*l+1)/(2*l-1)) * (l+abs(m))*(l-abs(m)) )...
                                        *Tesseral_Plm(l-1,m,THETA);
       if m==0

             Fx_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   .*cos(PHI)/sqrt(2*pi) ;
             Fy_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   .*sin(PHI)/sqrt(2*pi) ;
       elseif m==1
           Fx_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   *(1/2).*( Tesseral_Phim(sign(m)*(abs(m)+1),PHI)...
                                            + 1/sqrt(pi) ) ;

            Fy_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   *(1/2).*( Tesseral_Phim(-sign(m)*(abs(m)+1),PHI)...
                                           )*sign(m) ;
       elseif m==-1

           Fx_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   *(1/2).*( Tesseral_Phim(sign(m)*(abs(m)+1),PHI)...
                                           ) ;

            Fy_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   *(1/2).*( Tesseral_Phim(-sign(m)*(abs(m)+1),PHI)...
                                            -1/sqrt(pi) )*sign(m) ;
       else
            Fx_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   *(1/2).*( Tesseral_Phim(sign(m)*(abs(m)+1),PHI)...
                                            +Tesseral_Phim(sign(m)*(abs(m)-1),PHI) ) ;

            Fy_lm(:,:,lindx,mindx) = (-1/(2*l+1))*...
                                 ( factor_rho1(:,:,lindx,mindx)-factor_rho2(:,:,lindx,mindx) )...
                                   *(1/2).*( Tesseral_Phim(-sign(m)*(abs(m)+1),PHI)...
                                            -Tesseral_Phim(-sign(m)*(abs(m)-1),PHI) )*sign(m) ;
       
       end
        Fz_lm(:,:,lindx,mindx) = (1/(2*l+1))*( factor_z1(:,:,lindx,mindx)+factor_z2(:,:,lindx,mindx) )...
                                   .*Tesseral_Phim(m,PHI) ;
    end
end

Fx_lm_direct = sin(THETA).*cos(PHI).*Ylm;
Fy_lm_direct = sin(THETA).*sin(PHI).*Ylm;
Fz_lm_direct = cos(THETA).*Ylm;

%% ---- Comparison plot: recurrence relation vs direct (r_hat*Y_lm), per (l,m) ----
% Fx_lm, Fy_lm, Fz_lm are scalar functions of (theta,phi) (x*Y_lm, y*Y_lm,
% z*Y_lm), so they are compared as coloured surfaces, not as a vector field.
for lindx = 1:L
    for mindx = 1:M
        plotRecurrenceComparisonSurfaces(PHI, THETA, ...
            Fx_lm(:,:,lindx,mindx),        Fy_lm(:,:,lindx,mindx),        Fz_lm(:,:,lindx,mindx), ...
            Fx_lm_direct(:,:,lindx,mindx), Fy_lm_direct(:,:,lindx,mindx), Fz_lm_direct(:,:,lindx,mindx), ...
            lval(lindx), mval(mindx));
    end
end

%% ================= Local functions ================= %%

function plotRecurrenceComparisonSurfaces(PHI, THETA, Fx_num, Fy_num, Fz_num, Fx_dir, Fy_dir, Fz_dir, lVal, mVal)
%PLOTRECURRENCECOMPARISONSURFACES  Tiled scalar surface comparison of the
%   recurrence-relation reconstruction against the direct computation,
%   for a single (l,m) pair. One row per scalar quantity (Fx, Fy, Fz),
%   recurrence in the left column and direct in the right column, with a
%   shared color/height scale within each row so the two panels are
%   directly comparable.

    figure('Name', sprintf('Recurrence relation scalar comparison, l=%d, m=%d', lVal, mVal));
    tl = tiledlayout(3, 2, 'TileSpacing','compact', 'Padding','compact');

    plotScalarSurfPair(PHI, THETA, Fx_num, Fx_dir, '$F_x$')
    plotScalarSurfPair(PHI, THETA, Fy_num, Fy_dir, '$F_y$')
    plotScalarSurfPair(PHI, THETA, Fz_num, Fz_dir, '$F_z$')

    title(tl, sprintf('$F_{lm}$ components, $l=%d$, $m=%d$', lVal, mVal), 'Interpreter','latex')
end

function plotScalarSurfPair(PHI, THETA, Vnum, Vdir, compLabel)
%PLOTSCALARSURFPAIR  One row of a tiledlayout: recurrence-relation
%   surface (left) vs direct surface (right) for a single scalar
%   quantity, mapped onto the unit sphere (radius = 1). Since radius is
%   fixed, the scalar value is shown via face color instead of height,
%   with a shared color scale within each row so the two panels are
%   directly comparable.

    X = sin(THETA) .* cos(PHI);
    Y = sin(THETA) .* sin(PHI);
    Z = cos(THETA);

    zmax = max(abs([Vnum(:); Vdir(:)]), [], 'omitnan');
    if zmax == 0 || isnan(zmax)
        zmax = 1;   % avoid a degenerate zero-range scale for an all-zero field
    end

    nexttile;
    surf(X, Y, Z, Vnum, 'EdgeColor','none', 'FaceColor','interp')
    colormap(gca, 'parula'); clim([-zmax zmax])
    axis equal; view(45,30); grid on
    xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$z$','Interpreter','latex')
    title(compLabel + " (recurrence relation)", 'Interpreter','latex')
    colorbar

    nexttile;
    surf(X, Y, Z, Vdir, 'EdgeColor','none', 'FaceColor','interp')
    colormap(gca, 'parula'); clim([-zmax zmax])
    axis equal; view(45,30); grid on
    xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$z$','Interpreter','latex')
    title(compLabel + " (direct, $\hat{r}\,Y_{lm}$)", 'Interpreter','latex')
    colorbar
end