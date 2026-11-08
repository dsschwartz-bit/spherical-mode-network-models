%BICON_TESTING  Cone-harmonic eigenvalues for a biconical antenna.
%   For cone half-angle t0 = pi/4, finds the first Nmax non-integer orders
%   nu_n with L_nu(t0) = 0, where L_nu is the odd cone harmonic
%   (coneHarmonic), and their sensitivity d(nu_n)/d(t0) by finite
%   difference; compares both with the asymptotic approximations
%     nu_n ~ sqrt((n*alpha0)^2 - 1/4) - 1/2,  alpha0 = pi/(pi/2 - t0).
%   Then evaluates the coupling terms Tlnu between odd spherical modes l
%   and cone modes nu_n, and the cone-mode norms Normp2_nu.
%   Requires Symbolic Math Toolbox.

clear
close all
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/
nu1=3.462;m=0;
th=linspace(0,pi,201).';

Lnu=coneHarmonic(nu1,th,"odd");
plot(th,Lnu)
xline(pi/4)
xline(3*pi/4)
yline(0)
%% Zeros of L_nu(t0) in nu
syms nu theta x Lnu temp
t0=sym(pi/4)
x(theta)=cos(theta)
x0=x(t0)
Lnu(nu,theta)=coneHarmonic(nu,theta,"odd")


max_nu=20;inc=1;
start=(0:inc:max_nu);   % search intervals [start(m), start(m+1)] for vpasolve

Nmax=3;
nu_n=zeros(1,Nmax);
nu_n_approx=zeros(1,Nmax);

dnu_dt_n=zeros(1,Nmax);
dnu_dt_n_approx=zeros(1,Nmax);

n=1;
alpha0=pi/(pi/2-t0);
dt=eps;   % finite-difference step in t0

for m=1:(length(start)-1)
    
    nu_n_approx(n)=sqrt((n*alpha0)^2-(1/4))-(1/2);
    
    
    sym1=vpasolve(Lnu(nu,t0),nu,[start(m),start(m+1)]);
   
    
    nu_sol=double(sym1);

    if ~(~mod(nu_sol,1) || any(nu_n==nu_sol))
            nu_n(n)=nu_sol ;
            dnu_dt_n_approx(n)=alpha0*(n*alpha0)^2/(pi*(nu_n(n)+1/2));
            
            sym2=vpasolve(Lnu(nu,t0+dt),nu,[start(m),start(m+1)]);
            sym3=(sym2-sym1)/dt;
            dnu_dt_n(n)=double(sym3);

            n=n+1;
            if n>Nmax
                break
            end
    end
end

nu_n=nonzeros(nu_n);nu_n=reshape(nu_n,1,length(nu_n))
nu_n_approx

dnu_dt_n=nonzeros(dnu_dt_n);dnu_dt_n=reshape(dnu_dt_n,1,length(dnu_dt_n))
dnu_dt_n_approx

N=length(nu_n);
figure
yyaxis left
plot(nu_n,"o")
hold on
plot(nu_n_approx(1:N),"--")
hold off
yyaxis right
plot(dnu_dt_n,"^")
hold on
plot(dnu_dt_n_approx(1:N),"--")
hold off


%% L_nu over the gap region, exact vs approximate nu_1
th=linspace(double(t0),double(pi-t0),101).';
th0=double(t0);
figure
plot(th/pi,Lnu(double(nu_n(1)),th),"LineWidth",2)
hold on
plot(th/pi,Lnu(double(nu_n_approx(1)),th),"--")
hold off
xlabel("$\theta/\pi$","Interpreter","latex")
xlim([double(t0),double(pi-t0)]/pi)
grid
%% Coupling between spherical modes (odd l) and cone modes nu_n
syms dLnu_dt test
dLnu_dt(nu,theta)=diff(Lnu(nu,theta),theta);
lval=1:2:2*N-1;
Tlnu=zeros(length(lval),N);inttest=Tlnu;
Normp2_nu=zeros(1,N);
for lindx=1:N
    for n=1:N
        lfactor=lval(lindx)*(lval(lindx)+1);
        nu_factor=nu_n(n)*(nu_n(n)+1);
        Normp2_nu(n) = 2*pi*double((2*nu_factor/(2*nu_n(n)+1))*sin(t0)*dLnu_dt(nu_n(n),t0)^2/dnu_dt_n(n));
        Tlnu(lindx,n) = 2*subs(sin(theta)*FastLegendre_Plm(lval(lindx),0,x)*dLnu_dt(nu_n(n),theta)/(nu_factor-lfactor),theta,t0) ;
    end
end
figure
plot(th/pi,dLnu_dt(double(nu_n(1)),th),"LineWidth",2)
hold on
plot(th/pi,dLnu_dt(double(nu_n_approx(1)),th),"--")
hold off
xlabel("$\theta/\pi$","Interpreter","latex")
xlim([double(t0),double(pi-t0)]/pi)
grid