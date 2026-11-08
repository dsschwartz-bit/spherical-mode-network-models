%% Bicon_Internal_External.m
% Interior (conical-mode) and exterior (spherical-mode) fields of a FEKO
% biconical antenna compared with the analytic mode-matching solution
% (Janaswamy 2022).
%   Exterior (r > a): psiTM_l from projecting FEKO E_r onto Y_l0 (odd l);
%     analytic c_l from Bicon_Coeff_Symbolic_2.
%   Interior (theta0 < theta < pi - theta0, r < a): chiTM_n from
%     projecting FEKO E_r onto the conical harmonics L_nu_n(theta)
%     (coneHarmonic, non-integer degrees nu_n); analytic a_nu_n from
%     c_l via Eq. (26), (28), (30); TEM term chi0 / a_0.
% E_r and H_phi are rebuilt from both coefficient sets and plotted
% against the FEKO fields at f0, and the coefficients are plotted vs
% frequency. The analytic results are scaled by I0 = chi0/a_0_conv.
%
% Requires, in the FEKO project folder: <antenna_name>.out,
% <antenna_name>_NearField_Exterior_2.efe/.hfe and
% <antenna_name>_NearField_Interior_2.efe/.hfe (port reference impedance
% eta0). Uses the Symbolic Math Toolbox.

clear
close all
clc

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
data_dir     = "C:\path\to\FEKO\bicon";    % FEKO project folder with the exports
antenna_name = "bicon_edgeport_symmetry";    % FEKO model name (file prefix)
% =========================

%% --- Plot Settings ---
LW = 2;    % LineWidth
FS = 14;   % FontSize

baseColor=theme().BaseColorStyle;
if strcmp(baseColor,"dark")
    pointColor=[1,1,1];
elseif strcmp(baseColor,"light")
    pointColor=[0,0,0];
else
    pointColor=[0.5,0.5,0.5];
end
clearvars baseColor

set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% --- Physical Constants ---
c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;       % Permeability of free space (H/m)
eps0 = 1/c0^2/mu0;      % Permittivity of free space (F/m)
eta0 = sqrt(mu0/eps0);  % Characteristic impedance of free space (Ohm)
%% --- Cone half-angle ---
theta0_deg = 45;   
theta0=theta0_deg*pi/180;
theta0_sym=sym(theta0_deg/180*pi);
%% --- Simulation Parameters ---
lmax = 41;              % Maximum spherical mode order
f0   = 1.5e9;             % Reference frequency for single-frequency plots (Hz)
fmin = 0.3e9;           % Minimum frequency for plot x-axis (Hz)
fmax = 5e9;             % Maximum frequency for plot x-axis (Hz)

%% --- File Paths ---
fileprefix   = fullfile(data_dir, antenna_name);

%% --- Import FEKO Port Parameters ---
% V0: source voltage, Z0: port reference impedance, R0: port location,
% has_gnd_plane: PEC ground-plane flag, even_odd_excitation: symmetry flag,
% a: antenna apothem (NaN if not in the .out file), Zant: input impedance.
[V0, Z0, R0, has_gnd_plane, even_odd_excitation, a,Zant,Freq] = ReadFEKOout(fileprefix+".out");

% The analytic comparison assumes a port reference impedance of eta0.
if abs(Z0-eta0) > 0.01
    error("antenna port impedance should be eta0 for now")
else
    Z0 = eta0;
end
Y0 = 1/Z0;
Yant = 1./Zant;
numf = length(Freq);

if fmax>Freq(numf)
    fmax=Freq(numf);
end
if fmin<Freq(1)
    fmin=Freq(1);
end

%% --- Plot Antenna Admittance vs Frequency ---
fignum = 1; figure(fignum) ; fignum=fignum+1;
plot(Freq,real(Yant),"Color","b","LineWidth",LW)
hold on
plot(Freq,imag(Yant),"Color","r","LineWidth",LW)
hold off
title("Admittance","FontSize",FS)
legend("$G$","$B$","FontSize",FS)
xlim([fmin,fmax])
grid on

%% --- Import FEKO Near Fields ---
% E/H arrays: (numth x numph x numr x 3 x numf); dim 4 = [phi, theta, r].
% Exterior (r > a): spherical modes. Interior (r < a, between the cones): conical modes.
[Eptr_ext, r_ext, theta_ext, phi_ext, ~] = ReadFE_spherical(fileprefix+"_NearField_Exterior_2.efe");
[Eptr_int, r_int, theta_int, phi_int, ~] = ReadFE_spherical(fileprefix+"_NearField_Interior_2.efe");

[Hptr_ext, ~, ~, ~, ~] = ReadFE_spherical(fileprefix+"_NearField_Exterior_2.hfe");
[Hptr_int, ~, ~, ~, ~] = ReadFE_spherical(fileprefix+"_NearField_Interior_2.hfe");

% Set E_phi and E_r to zero for theta <= theta0 (on/inside the PEC cone).
index_temp=find((theta_int-theta0)<1e-6);
Eptr_int(index_temp,:,1:(end),1,:)=0;
Eptr_int(index_temp,:,1:(end),3,:)=0;
clearvars index_temp

% Limit lmax to the Nyquist limit of the exterior theta sampling.
delta_theta_ext=theta_ext(2)-theta_ext(1);
if lmax>((2*pi)/delta_theta_ext/2)
    lmax=floor((2*pi)/delta_theta_ext/2);
end
% Highest conical degree searched / kept.
nu_max=lmax;

numth_ext = length(theta_ext); numph_ext = length(phi_ext); numr_ext = length(r_ext);
numth_int = length(theta_int); numph_int = length(phi_int); numr_int = length(r_int);
% Apothem estimate if not available from the .out file
if isnan(a)
    a = r_ext(end) - c0/Freq(1)/2/pi;
end

% Frequency index; the outermost exported radius of each region is used.
[~,index_f] = min(abs(Freq-f0));
index_r_ext=length(r_ext);
index_r_int=length(r_int);
index_phi=1;

% Symmetry factor: 2 if only the upper hemisphere was exported, else 1
if max(theta_ext) > (pi/2)
    sym_factor_ext = 1;
else
    sym_factor_ext = 2;
end
if max(theta_int) > (pi/2)
    sym_factor_int = 1;
else
    sym_factor_int = 2;
end

k = 2*pi*Freq/c0;
ka = k*a;
kr_ext=k*r_ext(index_r_ext);
kr_int=k*r_int(index_r_int);


%% --- Determine Spherical Mode Coefficients (psiTM) ---
% Project the exterior E_r onto Yl0_pnorm = Y_l^0(theta)/sqrt(l(l+1))
% (odd l, m = 0) to get the TM Hertz-potential coefficients psi_l.

lval = (1:2:lmax);   % l = 1, 3, 5, ...
mval = 0;            % m = 0 only
L    = length(lval);
l_factor=lval.*(lval+1);
j_l_1=((1j).^(lval+1));
Yl0_pnorm    = zeros(numth_ext,numph_ext,L);
psiTM  = nan(numf,L); psiTE = psiTM;boundary_norm_l=psiTM;prop_factor_ext=psiTM;

% Squared norms l(l+1) * 2pi * 2/(2l+1) of the exterior modes
Normp2_l = l_factor .* (2*pi * 2 ./ (2*lval+1));

for lindx = 1:L

    % Y_l^0(theta)/sqrt(l(l+1)) (harmonicY, third_party)
    Yl0_pnorm(:,:,lindx) = harmonicY(lval(lindx),0,theta_ext,phi_ext*ones(size(theta_ext)), ...
        'type','real','norm',true,'phase',true) * sqrt(2^(mval~=0)) / sqrt(l_factor(lindx));

    for f = 1:numf
        % No boundary normalisation (coefficients referred to the far field).
        boundary_norm_l(f,lindx)=1;

        % Radial factor h_l^(2)(kr)/(kr j^(l+1)) at the exterior extraction radius
        prop_factor_ext(f,lindx) = sphericalHankel(lval(lindx),2,kr_ext(f)) / boundary_norm_l(f,lindx) / kr_ext(f) / j_l_1(lindx);

        % psiTM: projection of E_r onto Yl0_pnorm, divided by -jw and the radial factor
        psiTM(f,lindx) = (sqrt(eta0)/mu0) * 2*pi * sym_factor_ext * ...
            trapz(theta_ext, Eptr_ext(:,:,index_r_ext,3,f).*(Yl0_pnorm(:,:,lindx)).*sin(theta_ext), 1) ...
            / prop_factor_ext(f,lindx) / (-1j*2*pi*Freq(f));
    end
    
end
clearvars lindx
%% --- Analytic Biconical Antenna Coefficients ---
% Mode-matching solution (Bicon_Coeff_Symbolic_2, 1 % switch-over
% threshold): raw c_l and the conical eigenvalues nu with d(nu)/d(theta0).
% Lnu_sym = coneHarmonic(nu,theta,"odd") and its theta-derivative
% dLnu_dth_sym are kept symbolic for T_lnu and Normp2_nu below.
converted=false;
syms nu_sym theta_sym Lnu_sym(nu_sym,theta_sym) dLnu_dth(nu_sym,theta_sym)

Lnu_sym(nu_sym,theta_sym)=coneHarmonic(nu_sym,theta_sym,"odd");
dLnu_dth_sym(nu_sym,theta_sym)=diff(Lnu_sym(nu_sym,theta_sym),theta_sym);

[c_l, nu, dnu_dth, nu_approx_idx, dnu_approx_idx] = Bicon_Coeff_Symbolic_2(a, theta0_sym, Freq, lmax,nu_max,1,converted);fignum=fignum+1;


nu=nu(find(nu<nu_max));
N = length(nu);
dnu_dth=dnu_dth(1:N);

nu_factor=nu.*(nu+1);

conv_factor_l=zeros(numf,L);
conv_factor_nu=zeros(numf,N);


% Normp2_nu: squared norm of conical mode n over theta0..pi-theta0 and phi,
%   2pi * 2nu(nu+1)/(2nu+1) * sin(theta0) * L_nu'(theta0)^2 / (dnu/dtheta0)
%   (Janaswamy 2022, Eq. 28), including the nu(nu+1) energy factor.
Normp2_nu = 2*pi*double((2*nu_factor./(2*nu+1))*sin(theta0_sym).*subs(dLnu_dth_sym(nu_sym,theta0_sym),nu_sym,nu).^2./dnu_dth);
% boundary_norm_nu: boundary normalisation of the conical coefficients (set to 1 below)
boundary_norm_nu=nan(numf,N);

for f=1:numf
    % conv_factor_l: raw c_l -> Hertz-potential psi_l convention (exterior, outgoing)
    if ~converted
        conv_factor_l(f,:)  = sqrt(eta0)/(2*pi) * (1j.^(lval+1)) .* sqrt(Normp2_l) ./ sphericalRicatti_h(lval,2,ka(f));
    else 
        conv_factor_l=1;
    end
    % conv_factor_nu: raw a_nu -> Hertz-potential convention (interior).
    % Interior modes are referred to the antenna boundary r = a (the
    % radial factor below divides by ka j_nu(ka)); exterior modes to the
    % far field.
    boundary_norm_nu(f,:)=1;j_nu_1=ones(size(nu));

    conv_factor_nu(f,:) = sqrt(eta0)/(2*pi) * (j_nu_1) .* sqrt(Normp2_nu) ;
end

%% --- Compute a_nu using Eq. (26), (28), (30) of Janaswamy (2022) ---
% Tlnu(n,lindx): coupling integral T_(l,nu_n), Eq. (26), in closed form
%   2 sin(theta0) P_l(cos theta0) L_nu'(theta0) / (nu(nu+1) - l(l+1));
%   the 2 covers both halves of [theta0, pi-theta0].
% a_nu: interior conical amplitudes, Eq. (30), from the raw c_l.
% KYt: normalised admittance at r = a; Gamma_anlytc its reflection
%   coefficient; a_0_conv: TEM (nu = 0) amplitude, Normp2_0 its squared norm.

Tlnu=zeros(N,L);
Pl_th0=zeros(1,L);
a_nu = zeros(numf, N);
KYt=zeros(numf,1);
Normp2_0 = 2*pi*2*log(cot(theta0/2));

for n=1:N
   for lindx=1:L
       Pl_th0(lindx)=subs(FastLegendre_Plm(lval(lindx),0,cos(theta_sym)),theta0);
       Tlnu(n,lindx) = 2*subs( sin(theta_sym) * FastLegendre_Plm(lval(lindx),0,cos(theta_sym))...
                        *dLnu_dth_sym(nu(n),theta_sym) / (nu_factor(n)-l_factor(lindx)), theta_sym,theta0) ;  
   end
       a_nu(:,n) = ( (2*pi)/Normp2_nu(n) )  * sum(l_factor.*c_l.*Tlnu(n,:) , 2) ;
end
KYt=( sum(c_l.*Pl_th0 , 2) ) / log(cot(theta0/2));
Gamma_anlytc=(1-KYt)./(1+KYt);
a_0_conv=sqrt(eta0*Normp2_0)/(4*pi).*(1+KYt);

% Interior (a_nu_conv) and exterior (c_l_conv) coefficients in the Hertz-potential convention
a_nu_conv = a_nu .* conv_factor_nu ./ boundary_norm_nu ;
c_l_conv = c_l .* conv_factor_l ./ boundary_norm_l ;
clearvars n lindx Tlnu


%% --- Determine Conical Mode Coefficients (chiTM) ---
% Project the interior E_r onto Lnu_pnorm = L_nu_n(theta)/sqrt(Normp2_nu(n))
% to get the conical TM coefficients chiTM. The radial factor uses the
% regular (standing-wave) j_nu, referred to r = a.

Lnu_pnorm=zeros(numth_int,numph_int,N);dLnu_dth_norm=Lnu_pnorm;
chiTM  = nan(numf,N); 
prop_factor_int=chiTM;

for n = 1:N
    % Normalised conical harmonic, tiled over phi
    Lnu_pnorm(:,:,n) = repmat( double( Lnu_sym(nu(n),theta_int) ) / sqrt(Normp2_nu(n)) , [1,numph_int] ) ;

    for f = 1:numf     
   
        % Radial factor j_nu(kr) / (kr * ka j_nu(ka)) at the interior radius
        prop_factor_int(f,n) = sphericalBesselj(nu(n),kr_int(f)) / sphericalRicatti_j(nu(n),ka(f)) * boundary_norm_nu(f,n) / kr_int(f) / (j_nu_1(n));
        

        % chiTM: projection of E_r onto Lnu_pnorm, divided by -jw and the radial factor
        chiTM(f,n) = (sqrt(eta0)/mu0) * 2*pi * sym_factor_int * ...
            trapz(theta_int, Eptr_int(:,:,index_r_int,3,f).*(Lnu_pnorm(:,:,n)).*sin(theta_int), 1) ...
            / prop_factor_int(f,n) / (-1j*2*pi*Freq(f));
        
        

   end
end
%% --- TEM (nu = 0) term from the extracted exterior coefficients ---
% numerator/denominator: projections of the exterior field (from psiTM) on
% the TEM mode at r = a; their ratio is the normalised admittance Ya_norm,
% Gamma the reflection coefficient at r = a and chi0 the TEM amplitude.
numerator=ones(numf,1);denominator=ones(numf,1);
for f=1:numf
    numerator(f) = (2*pi/sqrt(Normp2_0))*...
                    sum( ...
                        psiTM(f,:) ...
                        .*( sphericalRicatti_h(lval,2,ka(f))./j_l_1 ) ...
                        *2.*( Pl_th0./sqrt(Normp2_l) )  ...
                    ,2 );
    denominator(f) = (2*pi/sqrt(Normp2_0))*...
                    sum( ...
                        psiTM(f,:) ...
                        *1j.*( sphericalRicatti_h_prime(lval,2,ka(f))./j_l_1 ) ...
                        *2.*( Pl_th0./sqrt(Normp2_l) ) ...
                    ,2 );
end

Ya_norm=numerator./denominator;
Gamma=(1-Ya_norm)./(1+Ya_norm);
chi0=sqrt(numerator.*denominator./(1-Gamma.^2));
% Scale of the analytic solution, equivalent to I0 in Janaswamy (2022)
scale_factor_ext=chi0./a_0_conv;

Gamma_input = (Y0-Yant)./(Y0+Yant);   % reflection coefficient at the FEKO port (not used below)
% TEM radial dependence (incident + reflected) at the interior radius
prop_factor_H0 = ( exp(-1j*( kr_int-ka ) ) - Gamma.*exp( 1j*( kr_int-ka ) ) )/r_int(index_r_int);

%% --- Reconstruct Exterior Near-Field Er from Coefficients ---
% E_r and H_phi at the exterior extraction radius from psiTM (FEKO) and
% c_l_conv (analytic), with the outgoing spherical-mode formulas.

Er_l  = zeros(numth_ext,numph_ext,numf,L); 
Er_l_anlytc = Er_l;

Hph_l = zeros(numth_ext,numph_ext,numf,L); Hph_l_anlytc = Hph_l;

for f = 1:numf
    for lindx = 1:L
        % From the extracted psiTM
        Er_l(:,:,f,lindx) = (-1j*2*pi*Freq(f))*l_factor(lindx)*(mu0/sqrt(eta0)) * ...
                                  psiTM(f,lindx)*prop_factor_ext(f,lindx)*(Yl0_pnorm(:,:,lindx));

        % From the analytic c_l_conv
        Er_l_anlytc(:,:,f,lindx) = (-1j*2*pi*Freq(f))*l_factor(lindx)*(mu0/sqrt(eta0)) * ...
                                         c_l_conv(f,lindx)*prop_factor_ext(f,lindx)*(Yl0_pnorm(:,:,lindx));

        % H_phi uses P_l^1(cos theta) = dP_l/dtheta (with Condon-Shortley phase)
        Pl1=FastLegendre_Plm(lval(lindx),1,cos(theta_ext));

        Hph_l(:,:,f,lindx) = (1/sqrt(eta0))*psiTM(f,lindx)*sphericalRicatti_h(lval(lindx),2,kr_ext(f))/j_l_1(lindx) ...
                                        *(-Pl1/r_ext(index_r_ext))/sqrt(Normp2_l(lindx)) / boundary_norm_l(f,lindx) ;
       
        Hph_l_anlytc(:,:,f,lindx) = (1/sqrt(eta0))*c_l_conv(f,lindx)*sphericalRicatti_h(lval(lindx),2,kr_ext(f))/j_l_1(lindx) ...
                                        *(-Pl1/r_ext(index_r_ext))/sqrt(Normp2_l(lindx)) / boundary_norm_l(f,lindx) ;
    end
    
end
clearvars lindx f



% Sum over modes
Er_exterior = sum(Er_l(:,:,:,:), 4);
Hph_exterior = sum(Hph_l, 4);
Er_exterior_anlytc = sum(Er_l_anlytc(:,:,:,:), 4);
Hph_exterior_anlytc = sum(Hph_l_anlytc(:,:,:,:), 4);

if isnan(scale_factor_ext)
    scale_factor_ext=ones(numf,1);
end
%% --- Plot Er External Reconstruction Comparison ---
% E_r vs theta at the exterior radius and f0: from psiTM, from the
% analytic c_l (times I0), and the FEKO field (markers).
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_ext/pi, abs(Er_exterior(:,index_phi,index_f)),"LineWidth",LW+3,    "LineStyle","-")
plot(theta_ext/pi, abs(Er_exterior_anlytc(:,index_phi,index_f)*scale_factor_ext(index_f)), "LineWidth",LW+1,"LineStyle","-")
plot(theta_ext/pi, abs(Eptr_ext(:,index_phi,index_r_ext,3,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*","Color",pointColor)
hold off
ylabel("magnitude (V/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_ext/pi, angle(Er_exterior(:,index_phi,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta_ext/pi, angle(Er_exterior_anlytc(:,index_phi,index_f)*scale_factor_ext(index_f))*180/pi,"LineWidth",LW,"LineStyle",":")
plot(theta_ext/pi, angle(Eptr_ext(:,1,index_r_ext,3,index_f))*180/pi,"LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",pointColor)
hold off
legend("Reconstructed $\psi^{TM}_{l0}$","Reconstructed Analytic $c_l$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$E_r$ reconstructed @r="+num2str(r_ext(index_r_ext)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $l_{max}$="+num2str(lmax),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% --- Plot Hphi External Reconstruction Comparison ---
% eta0*H_phi vs theta at the exterior radius and f0 (same three traces).
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_ext/pi, eta0*abs(Hph_exterior(:,index_phi,index_f)),"LineWidth",LW+3,    "LineStyle","-")
plot(theta_ext/pi, eta0*abs(Hph_exterior_anlytc(:,index_phi,index_f)*scale_factor_ext(index_f)), "LineWidth",LW+1,"LineStyle","-")
plot(theta_ext/pi, eta0*abs(Hptr_ext(:,index_phi,index_r_ext,1,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*","Color",pointColor)
hold off
ylabel("magnitude (V/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_ext/pi, angle(Hph_exterior(:,index_phi,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta_ext/pi, angle(Hph_exterior_anlytc(:,index_phi,index_f).*scale_factor_ext(index_f))*180/pi,"LineWidth",LW,"LineStyle",":")
plot(theta_ext/pi, angle(Hptr_ext(:,index_phi,index_r_ext,1,index_f))*180/pi,"LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",pointColor)
hold off
legend("Reconstructed $\psi^{TM}_{l0}$","Reconstructed Analytic $c_l$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$ \eta H_\phi$ reconstructed @r="+num2str(r_ext(index_r_ext)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $l_{max}$="+num2str(lmax),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% --- Reconstruct Interior Near-Field Er from Coefficients ---
% E_r and H_phi at the interior radius from chiTM (FEKO) and a_nu_conv
% (analytic, Eq. 30), with the standing-wave (j_nu) conical-mode formulas;
% H_phi also includes the TEM term (chi0 / a_0_conv, csc(theta) dependence).

Er_nu  = zeros(numth_int,numph_int,numf,N);
Er_nu_anlytc = Er_nu;
Hph_nu  = zeros(numth_int,numph_int,numf,N);
Hph_nu_anlytc = Hph_nu;

Hph_0=zeros(numth_int,numph_int,numf);Hph_0_anlytc=Hph_0;


for n = 1:N

    dLnu_dth_norm(:,:,n)=repmat( double(dLnu_dth_sym(nu(n),theta_int)) / sqrt(Normp2_nu(n)) ,[1,numph_int] ) ;

    for f = 1:numf
        % From the extracted chiTM
        Er_nu(:,:,f,n) = (-1j*2*pi*Freq(f))*nu_factor(n)*(mu0/sqrt(eta0)) * ...
                                  chiTM(f,n)*prop_factor_int(f,n)*(Lnu_pnorm(:,:,n));

        % From the analytic a_nu_conv
        Er_nu_anlytc(:,:,f,n) = (-1j*2*pi*Freq(f))*nu_factor(n)*(mu0/sqrt(eta0)) * ...
                                  a_nu_conv(f,n)*prop_factor_int(f,n)*(Lnu_pnorm(:,:,n));
        

        Hph_nu(:,:,f,n) = (1/sqrt(eta0))*chiTM(f,n)...
                            *(sphericalRicatti_j(nu(n),kr_int(f))/sphericalRicatti_j(nu(n),ka(f))) ...
                                        *(-dLnu_dth_norm(:,:,n))/r_int(index_r_int) ;
       
        Hph_nu_anlytc(:,:,f,n) = (1/sqrt(eta0))*a_nu_conv(f,n)...
                                    *(sphericalRicatti_j(nu(n),kr_int(f)) / sphericalRicatti_j(nu(n),ka(f))) ...
                                        *(-dLnu_dth_norm(:,:,n))/r_int(index_r_int) ;
        
       
    end
    
end

for f=1:numf
 Hph_0(:,:,f) = 1/sqrt(eta0)*chi0(f)*prop_factor_H0(f)*csc(theta_int)/sqrt(Normp2_0) ;
 Hph_0_anlytc(:,:,f) = 1/sqrt(eta0)*a_0_conv(f)*prop_factor_H0(f)*csc(theta_int)/sqrt(Normp2_0) ;
end
clearvars n f
% Same I0 scale for the interior analytic fields
scale_factor_int=scale_factor_ext;
if isnan(scale_factor_int)
    scale_factor_int=ones(numf,1);
end
% Sum over modes (H_phi includes the TEM term)
Er_interior  = sum(Er_nu(:,:,:,:), 4);
Hph_interior = sum(Hph_nu, 4) + Hph_0;
Er_interior_anlytc = sum(Er_nu_anlytc(:,:,:,:), 4);
Hph_interior_anlytc = sum(Hph_nu_anlytc(:,:,:,:), 4)+ Hph_0_anlytc;

%% --- Plot Er Internal Reconstruction Comparison ---
% E_r vs theta at the interior radius and f0: from chiTM, from the
% analytic a_nu (times I0), and the FEKO field (markers).
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_int/pi, abs(Er_interior(:,index_phi,index_f)),"LineWidth",LW+3,    "LineStyle","-")
plot(theta_int/pi, abs(Er_interior_anlytc(:,index_phi,index_f)*scale_factor_int(index_f)), "LineWidth",LW+1,"LineStyle","-")
plot(theta_int/pi, abs(Eptr_int(:,index_phi,index_r_int,3,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*", "Color",pointColor)
hold off
ylabel("magnitude (V/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_int/pi, angle(Er_interior(:,index_phi,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta_int/pi, angle(Er_interior_anlytc(:,index_phi,index_f)*scale_factor_int(index_f))*180/pi,"LineWidth",LW,"LineStyle",":")
plot(theta_int/pi, angle(Eptr_int(:,index_phi,index_r_int,3,index_f))*180/pi,      "LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",pointColor)
hold off
legend("Reconstructed $\chi^{TM}_{\nu_n}$","Reconstructed Analytic $a_{\nu_n}$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$E_r$ reconstructed @r="+num2str(r_int(index_r_int)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $\nu_{max}$="+num2str(max(nu)),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% --- Plot Hphi Internal Reconstruction Comparison ---
% eta0*H_phi vs theta at the interior radius and f0 (same three traces).
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_int/pi, eta0*abs(Hph_interior(:,index_phi,index_f)),"LineWidth",LW+3,    "LineStyle","-")
plot(theta_int/pi, eta0*abs(Hph_interior_anlytc(:,index_phi,index_f).*scale_factor_int(index_f)), "LineWidth",LW+1,"LineStyle","-")
plot(theta_int/pi, eta0*abs(Hptr_int(:,index_phi,index_r_int,1,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*", "Color",pointColor)
hold off
ylabel("magnitude (V/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta_int/pi, angle(Hph_interior(:,index_phi,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta_int/pi, angle(Hph_interior_anlytc(:,index_phi,index_f).*scale_factor_int(index_f))*180/pi,"LineWidth",LW,"LineStyle",":")
plot(theta_int/pi, angle(Hptr_int(:,index_phi,index_r_int,1,index_f))*180/pi,      "LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",pointColor)
hold off
legend("Reconstructed $\chi^{TM}_{\nu_n}$","Reconstructed Analytic $a_{\nu_n}$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$\eta H_\phi$ reconstructed @r="+num2str(r_int(index_r_int)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $\nu_{max}$="+num2str(max(nu)),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% Plot Coefficients (exterior)
% psiTM (left axis, solid) and I0-scaled analytic c_l (right axis, dotted).
ymax=max(abs(1.1*psiTM),[],"all");
dynamic_range=ymax; %mag
ymin=ymax-dynamic_range;

figure(fignum);fignum = fignum+1;
yyaxis left
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
plot(Freq, abs(psiTM(:,1:L)), "LineWidth",LW-0.5, "LineStyle","-")
ylabel("$|\psi^{TM}_{l0}|$","FontSize",FS)
ylim([ymin,ymax])
yyaxis right
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
set(gca,'ColorOrderIndex',1)
plot(Freq, abs(c_l_conv(:,1:L).*scale_factor_ext ), "LineWidth",LW+0.5, "LineStyle",":")
ylabel("$|c_{l}|$","FontSize",FS)
title("Magnitude Simulated and Theoretical Coefficients","FontSize",FS)
xlim([fmin,fmax]); xlabel("Frequency (Hz)","FontSize",FS)
lgd = legend(cellstr(num2str(reshape(lval,L,1))),"FontSize",FS); title(lgd,"$l$")
ylim([ymin,ymax])
grid on

figure(fignum);fignum = fignum+1;
yyaxis left
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
set(gca,'ColorOrderIndex',1)
plot(Freq, angle(psiTM(:,1:L))*180/pi, "LineWidth",LW, "LineStyle","-")
ylabel("$\angle \psi^{TM}_{l0}$","FontSize",FS)
ylim([-180 180]);
yyaxis right
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
set(gca,'ColorOrderIndex',1)
plot(Freq, angle(c_l_conv(:,1:L).*scale_factor_ext)*180/pi, "LineWidth",LW, "LineStyle",":")
ylabel("$\angle c_{l}$","FontSize",FS)
lgd = legend(cellstr(num2str(reshape(lval,L,1))),"FontSize",FS); title(lgd,"$l$")
title("Phase Simulated and Theoretical Coefficients","FontSize",FS)
xlim([fmin,fmax]); xlabel("Frequency (Hz)","FontSize",FS)
ylim([-180 180]); grid on
%% Plot Coefficients (interior)
% chi0, chiTM (left axis, solid) and I0-scaled a_0, a_nu (right axis, dotted).
ymax=max(abs(1.1*chi0),[],"all");
dynamic_range=ymax; %mag
ymin=ymax-dynamic_range;

figure(fignum);fignum = fignum+1;
yyaxis left
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
plot(Freq, abs(chi0), "LineWidth",LW-0.5, "LineStyle","-","Color",pointColor)
hold on
set(gca,'ColorOrderIndex',1)
plot(Freq, abs(chiTM(:,1:N)), "LineWidth",LW-0.5, "LineStyle","-")
hold off
ylabel("$|\chi^{TM}_{\nu_n}|$","FontSize",FS)
ylim([ymin,ymax])
yyaxis right
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
plot(Freq, abs(a_0_conv.*scale_factor_int), "LineWidth",LW-0.5, "LineStyle",":","Color",pointColor)
hold on
set(gca,'ColorOrderIndex',1)
plot(Freq, abs(a_nu_conv(:,1:N).*scale_factor_int), "LineWidth",LW+1, "LineStyle",":")
ylabel("$|a_{\nu_n}|$","FontSize",FS)
title("Magnitude Simulated and Theoretical Coefficients","FontSize",FS)
xlim([fmin,fmax]); xlabel("Frequency (Hz)","FontSize",FS)
lgd = legend(cellstr(num2str([0;reshape(nu,N,1)])),"FontSize",FS); title(lgd,"$\nu$")
ylim([ymin,ymax])
grid on

figure(fignum);fignum = fignum+1;
yyaxis left
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
set(gca,'ColorOrderIndex',1)
plot(Freq, angle(chi0)*180/pi, "LineWidth",LW, "LineStyle","-","Color",pointColor)
hold on
set(gca,'ColorOrderIndex',1)
plot(Freq, angle(chiTM(:,1:N))*180/pi, "LineWidth",LW, "LineStyle","-")
ylabel("$\angle \chi^{TM}_{\nu_n}$","FontSize",FS)
ylim([-180 180]);
yyaxis right
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
plot(Freq, angle(a_0_conv.*scale_factor_int)*180/pi, "LineWidth",LW+1, "LineStyle",":","Color",pointColor)
hold on
set(gca,'ColorOrderIndex',1)
plot(Freq, angle(a_nu_conv(:,1:N).*scale_factor_int)*180/pi, "LineWidth",LW+1, "LineStyle",":")
hold off
ylabel("$\angle a_{\nu_n}$","FontSize",FS)
lgd = legend(cellstr(num2str([0;reshape(nu,N,1)])),"FontSize",FS); title(lgd,"$\nu$")
title("Phase Simulated and Theoretical Coefficients","FontSize",FS)
xlim([fmin,fmax]); xlabel("Frequency (Hz)","FontSize",FS)
ylim([-180 180]); grid on