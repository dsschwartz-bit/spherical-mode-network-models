function [cl, nu_n, dnu_n] = Bicon_Coeff(a, theta0, Freq, lmax, conversion)
%BICON_COEFF  Analytic TM spherical-mode coefficients of a biconical antenna.
%   [cl, nu_n, dnu_n] = Bicon_Coeff(a, theta0, Freq, lmax, conversion)
%   Mode-matching solution for a symmetric biconical antenna (cones of
%   half-angle theta0, arm length a), fed by the TEM mode. The interior
%   region is expanded in conical modes of non-integer degree nu_n (zeros
%   of coneHarmonic(nu,theta0,"odd")), here from the closed-form
%   approximation nu_n = sqrt((n alpha0)^2 - 1/4) - 1/2,
%   alpha0 = 1/(1/2 - theta0/pi). The exterior uses odd-l TM_l0 spherical
%   modes. The mode-coupling matrix G is built with the convergence-
%   accelerated series t1 - t2 + t3 and G c = -j p is solved per frequency.
%
%   Inputs:
%     a          arm length / exterior sphere radius (m)
%     theta0     cone half-angle (rad)
%     Freq       (1 x numf) frequencies (Hz)
%     lmax       highest (odd) degree; L = (lmax+1)/2 modes l = 1,3,...,lmax
%     conversion true: scale to the psi (Hertz-potential) convention of
%                this toolbox, sqrt(eta0)/(2pi) j^(l+1) Normp / (ka h_l^(2)(ka));
%                false (default): return the raw mode-matching c_l
%   Outputs:
%     cl     (numf x L) coefficients
%     nu_n   (1 x N) conical eigenvalues, N = round(floor((L-1)/2) * k_max a)
%     dnu_n  (1 x N) d(nu_n)/d(theta0)
%
%   References:
%     Janaswamy (2022), Eq. 11-12 (approximate nu_n and dnu_n).
%
%   See also Bicon_Coeff_Symbolic, coneHarmonic.
if mod(lmax,2)==0
    error("lmax should be odd")
end
if nargin < 5 || isempty(conversion)
    conversion = false;
end
c0=299792458; %m/s
mu0=4*pi*1e-7 ; %H/m
eps0=1/c0^2/mu0; %F/m
eta0=sqrt(mu0/eps0); %Ohm

%% --- Sweep and Discretization Parameters ---
numf = length(Freq);                          
freq_max=max(Freq);

kaMx = 2*pi*a*freq_max/c0;        
              
L=(lmax+1)/2;
Mifac = floor((L-1)/2); 
Nfac = Mifac;          % series length N scales with k_max*a
n = (1:round(Nfac*kaMx));      % conical-mode index of the truncated series
N=length(n);

%% --- Antenna and Geometry Parameters ---
alp0 = 1/(0.5-theta0/pi);           % alpha_0

lval=(1:2:lmax);
Normp = sqrt( lval.*(lval+1) .* 2*pi * 2./(2*lval+1) ) ;

%% --- Conical Waveguide Mode Eigenvalues (geometry only, no frequency dependence) ---
nu_n  = sqrt((n*alp0).^2-0.25) - 0.5;          % approximate roots of coneHarmonic(nu,theta0,"odd")
dnu_n = n.^2*alp0^3./(pi*(nu_n+0.5));           % d(nu_n)/d(theta0) of that approximation

cl = zeros(numf,L); c = cl;

%% --- Main Frequency Sweep Loop ---
for f = 1:numf
    ka = 2*pi*a*Freq(f)/c0;

    %% --- Interior radial functions ---
    B1 = besselj(nu_n-1/2,ka);
    B2 = besselj(nu_n+3/2,ka);
    B3 = besselj(nu_n+1/2,ka);
    JpJ = 0.5*(1/ka + (B1-B2)./B3);               % (x j_nu(x))'/(x j_nu(x)) at x = ka

    %% --- Mode Coupling Prefactors ---
    Fac1 = dnu_n.*(2*nu_n+1)./(nu_n.*(nu_n+1));

    B   = Fac1.*JpJ;                               % series kernel of G

    %% --- Double Mode Index Loop: Build G Matrix ---
    for mi = 1:L
        m = 2*mi-1;                          % odd exterior degree (here m is a degree, not an order)
        y = legendre(m,cos(theta0));
        p(mi) = y(1);                        % excitation vector: P_m(cos theta0)

        for elli = 1:mi
            ell = 2*elli-1;
            bm   = sqrt(m*(m+1))/alp0;
            bell = sqrt(ell*(ell+1))/alp0;

            %% --- Convergence-Accelerated Series for G Matrix ---
            % t1: series over conical modes; t2: its large-n asymptote
            % (subtracted); t3: closed-form sum of t2 (digamma psi).
            t1 = B./((nu_n.*(nu_n+1)-m*(m+1)).*(nu_n.*(nu_n+1)-ell*(ell+1)));
            t2 = 1./(pi*alp0^2*(bm+bell)*ka)*(1./((n-bm).*(n-bell))-1./((n+bm).*(n+bell)));

            if m ~= ell
                t3 = 1./(pi*alp0^2*(bm+bell)*ka)*(pi*sin(pi*(bm-bell))/((bm-bell)*sin(pi*bm)*sin(pi*bell)) ...
                    +1/(bm*bell) - 2*(psi(bm)-psi(bell))/(bm-bell));
            else  % m == ell: limiting form
                t3 = 1./(pi*alp0^2*(bm+bell)*ka)*((pi/sin(pi*bm))^2 + 1/bm^2 - 2*psi(1,bm));
            end

            gml(mi,elli) = sum(t1-t2)+t3;

            y2 = legendre(ell,cos(theta0));

            % Angular weight sin(theta0) P_m P_ell m(m+1) ell(ell+1); G is symmetric.
            gml(mi,elli) = gml(mi,elli)*sin(theta0)*y(1)*y2(1)*m*(m+1)*ell*(ell+1);
            gml(elli,mi)  = gml(mi,elli);

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
    c(f,:) = -1j*gml\p.';      % G c = -j p

    if conversion
        conv_factor  = sqrt(eta0)/(2*pi) * (1j.^(lval+1)) .* Normp ./ sphericalRicatti_h(lval,2,ka);
        cl(f,:) = c(f,:) .* conv_factor;
    else
        cl(f,:) = c(f,:);
    end

end  % frequency sweep loop

end