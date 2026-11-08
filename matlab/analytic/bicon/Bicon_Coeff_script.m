%% Bicon_Coeff_script.m
% Script version of Bicon_Coeff for a single frequency: biconical
% antenna with cone half-angle th0 = pi/4 and arm length a = lambda0/8
% (ka = pi/4 at f0 = 1 GHz). Builds the mode-coupling matrix G (and the
% K-matrix variants Kml, Kml2, Kml3), solves for the TM_l0 coefficients
% c, converts them to the psi convention (psi_TM), computes the
% normalised input admittance KYin, and plots |c_l| against the thin
% dipole coefficients cd. The admittance and psi_TM lines index c, KY_t
% and Gam as whole arrays, so they assume numf = 1.
% No input files.

clear
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

c0=physconst('lightspeed');
mu0=4*pi*1e-7 ; %H/m
eps0=1/c0^2/mu0;
eta0=sqrt(mu0/eps0);

%% --- Sweep and Discretization Parameters ---
numf = 1;                          % number of frequency points
freq_min=1e9;
freq_max=1e9;
Freq=linspace(freq_min,freq_max,numf);
KY_t=zeros(numf,1);Gam=KY_t;KYin=KY_t;
f0=1e9;k0=2*pi*f0/c0;a=pi/4/k0;
kaMn = 2*pi*a*freq_min/c0; kaMx = 2*pi*a*freq_max/c0;        % range of ka
Nth = 91;                        % theta points (th is not used further)
Mifac = 13; Nfac = 13;          % number of modes L and series terms N scale with k_max*a
th = linspace(0,pi/2,Nth);
%% --- Antenna and Geometry Parameters ---
th0 = pi/4                      % cone half-angle (rad)
alp0 = 1/(0.5-th0/pi);           % alpha_0
j = 1i;
K1 = 1/log(cot(th0/2));          % 1/ln(cot(th0/2)): biconical line impedance is (eta0/pi)/K1


L = ceil(Mifac*kaMx)+1;       % number of odd exterior modes
lmax=2*L-1;
lval=(1:2:lmax);
Normp = sqrt( lval.*(lval+1) .* 2*pi * 2./(2*lval+1) ) ;
n = (1:round(Nfac*kaMx));      % conical-mode index of the truncated series
N=n(end);
psi_TM = zeros(numf,L);c=psi_TM;
nu_n=zeros(numf,N);dnu_n=nu_n;t1=nu_n;
%% --- Main Frequency Sweep Loop ---
for f = 1:numf
    ka = 2*pi*a*Freq(f)/c0;

    %% --- Conical eigenvalues (Janaswamy 2022 approximation, see Bicon_Coeff) ---
    nu_n(f,:) = sqrt((n*alp0).^2-0.25) - 0.5;          % approximate roots of coneHarmonic(nu,th0,"odd")
    dnu_n(f,:) = n.^2*alp0^3./(pi*(nu_n(f,:)+0.5));          % d(nu_n)/d(th0)

    %% --- Interior radial functions ---
    B1 = besselj(nu_n(f,:)-1/2,ka);
    B2 = besselj(nu_n(f,:)+3/2,ka);
    B3 = besselj(nu_n(f,:)+1/2,ka);
    JpJ = 0.5*(1/ka + (B1-B2)./B3);               % (x j_nu(x))'/(x j_nu(x)) at x = ka

    %% --- Mode Coupling Prefactors ---
    Fac1 = dnu_n(f,:).*(2*nu_n(f,:)+1)./(nu_n(f,:).*(nu_n(f,:)+1));

    % Series kernels: B for G; Bp, Bp2, Bp3 (= Bp -/+ JpJ/ka) for the K variants.
    B   = Fac1.*JpJ;
    Bp  = Fac1.*(1+JpJ.^2-nu_n(f,:).*(nu_n(f,:)+1)./ka^2);
    Bp2 = Fac1.*(1+JpJ.^2-nu_n(f,:).*(nu_n(f,:)+1)/ka^2 + JpJ/ka);
    Bp3 = Fac1.*(1+JpJ.^2-nu_n(f,:).*(nu_n(f,:)+1)/ka^2 - JpJ/ka);

    %% --- Double Mode Index Loop: Build G and K Matrices ---
    for mi = 1:L
        m = 2*mi-1;                          % odd exterior degree (a degree, not an order)
        y = legendre(m,cos(th0));
        p(mi) = y(1);                        % excitation vector: P_m(cos th0)

        for elli = 1:mi
            ell = 2*elli-1;
            bm   = sqrt(m*(m+1))/alp0;
            bell = sqrt(ell*(ell+1))/alp0;

            %% --- Convergence-Accelerated Series for G Matrix ---
            % t1: series over conical modes; t2: its large-n asymptote
            % (subtracted); t3: closed-form sum of t2 (digamma psi).
            t1(f,:) = B./((nu_n(f,:).*(nu_n(f,:)+1)-m*(m+1)).*(nu_n(f,:).*(nu_n(f,:)+1)-ell*(ell+1)));
            t2 = 1./(pi*alp0^2*(bm+bell)*ka)*(1./((n-bm).*(n-bell))-1./((n+bm).*(n+bell)));

            if m ~= ell
                t3 = 1./(pi*alp0^2*(bm+bell)*ka)*(pi*sin(pi*(bm-bell))/((bm-bell)*sin(pi*bm)*sin(pi*bell)) ...
                    +1/(bm*bell) - 2*(psi(bm)-psi(bell))/(bm-bell));
            else  % m == ell: limiting form
                t3 = 1./(pi*alp0^2*(bm+bell)*ka)*((pi/sin(pi*bm))^2 + 1/bm^2 - 2*psi(1,bm));
            end

            gml(mi,elli) = sum(t1(f,:)-t2)+t3;

            %% --- K matrix variants (same series with kernels Bp, Bp2, Bp3) ---
            t4 = Bp ./((nu_n(f,:).*(nu_n(f,:)+1)-m*(m+1)).*(nu_n(f,:).*(nu_n(f,:)+1)-ell*(ell+1)));
            t5 = Bp2./((nu_n(f,:).*(nu_n(f,:)+1)-m*(m+1)).*(nu_n(f,:).*(nu_n(f,:)+1)-ell*(ell+1)));
            t6 = Bp3./((nu_n(f,:).*(nu_n(f,:)+1)-m*(m+1)).*(nu_n(f,:).*(nu_n(f,:)+1)-ell*(ell+1)));

            y2 = legendre(ell,cos(th0));

            % Angular weights P_m P_ell m(m+1) ell(ell+1) (times sin(th0) for G only)
            gml(mi,elli) = gml(mi,elli)*sin(th0)*y(1)*y2(1)*m*(m+1)*ell*(ell+1);
            Kml (mi,elli) = sum(t4)*y(1)*y2(1)*m*(m+1)*ell*(ell+1);
            Kml2(mi,elli) = sum(t5)*y(1)*y2(1)*m*(m+1)*ell*(ell+1);
            Kml3(mi,elli) = sum(t6)*y(1)*y2(1)*m*(m+1)*ell*(ell+1);

            % Symmetric in m and ell
            gml(elli,mi)  = gml(mi,elli);
            Kml (elli,mi) = Kml (mi,elli);
            Kml2(elli,mi) = Kml2(mi,elli);
            Kml3(elli,mi) = Kml3(mi,elli);

            %% --- Diagonal: exterior (outgoing) radiation term ---
            if elli == mi
                H1 = besselh(m-1/2,2,ka);
                H2 = besselh(m+3/2,2,ka);
                H3 = besselh(m+1/2,2,ka);
                H  = 1/ka + (H1-H2)/H3;      % 2 * (x h_m^(2)(x))'/(x h_m^(2)(x)) at x = ka
                gml(mi,elli) = -0.5*H*m*(m+1)/(2*m+1) + gml(mi,elli);
            end

        end  % elli loop
    end  % mi loop

    %% --- Solve for Biconical Antenna Mode Coefficients ---
    c(f,:) = -j*gml\p.';      % G c = -j p
    
    conv_factor = sqrt(eta0)/(2*pi) * (1j.^(lval+1)) .* Normp ./ sphericalRicatti_h(lval,2,ka) ;
    psi_TM(f,:) = c.*conv_factor ;
    %% --- Input admittance via the reflection coefficient ---
    KY_t(f) = K1*p*(c.');                              % normalised admittance at r = a
    Gam(f)  = (1-KY_t)/(1+KY_t)*exp(-j*2*ka);          % reflection coefficient moved to the feed (round trip 2ka)
    KYin(f) = (1-Gam)/(1+Gam);                         % normalised input admittance

    %% --- Reference: thin-dipole modal coefficients ---
    cd = sqrt((2*lval+1)/pi).*sphericalBesselj(lval,ka).*sphericalHankel(lval,2,ka)./(lval.*(lval+1)) ;
    

end  % frequency loop

%% --- Plot Mode Coefficient Magnitudes: Biconical vs Dipole ---
figure
semilogy(lval, abs(c(numf,:)/c(numf,1)), '-or', lval, abs(cd/cd(1)), '-og', 'markersize', 5)
axis([0 25 5e-4 2])
legend('Biconical', 'Dipole')
text(12, 0.1, '$\theta_0 = \pi/4, kL = \pi/4$', 'Interpreter','latex','FontSize',12)