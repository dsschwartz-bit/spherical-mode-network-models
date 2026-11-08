%ANTENNA_MODAL_IMPEDANCES  TM_l0 modal wave impedances of a biconical antenna from FEKO near fields.
% Projects the radial near field Er on a sphere onto the m = 0, odd-l
% spherical harmonics to get the TM coefficients psiTM_l, reconstructs Er,
% Etheta and Hphi from them, and plots the modal wave impedance
% Z_l = Etheta/Hphi vs frequency, compared with the equivalent ladder
% circuits of the TM10 and TM30 modes.
%
% Inputs (in data_dir): <antenna_name>.out, <antenna_name>_NearField1.efe,
%   <antenna_name>_NearField1.hfe (spherical near-field exports).
% Output: figures only.

clear
close all
clc
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
data_dir     = "C:\path\to\FEKO\project";   % folder containing the FEKO .out/.efe/.hfe exports
antenna_name = "bicon_edgeport_symmetry";
% =========================
fileprefix   = fullfile(data_dir, antenna_name);

%% --- Plot Settings ---
LW = 2;    % LineWidth
FS = 14;   % FontSize

set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

% Marker colour for the FEKO reference points, matched to the MATLAB theme
baseColor=theme().BaseColorStyle;
if strcmp(baseColor,"dark")
    pointColor=[1,1,1];
elseif strcmp(baseColor,"light")
    pointColor=[0,0,0];
else
    pointColor=[0.5,0.5,0.5];
end
clearvars baseColor

%% --- Physical Constants ---
c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;       % Permeability of free space (H/m)
eps0 = 1/c0^2/mu0;      % Permittivity of free space (F/m)
eta0 = sqrt(mu0/eps0);  % Characteristic impedance of free space (Ohm)

%% --- Simulation Parameters ---
lmax = 9;              % Maximum spherical mode order
f0   = 5e9;             % Reference frequency for single-frequency plots (Hz)
fmin = 0.3e9;           % Minimum frequency for plot x-axis (Hz)
fmax = 5e9;             % Maximum frequency for plot x-axis (Hz)

%% --- Import FEKO Port Parameters ---
% V0: open-circuit voltage, Z0: port reference impedance,
% R0: port location, has_gnd_plane: PEC ground plane flag,
% even_odd_excitation: log2(ports) if mode symmetry detected, a: antenna apothem,
% Zant: input impedance vs Freq
[V0, Z0, R0, has_gnd_plane, even_odd_excitation, a,Zant,Freq] = ReadFEKOout(fileprefix+".out");
Y0 = 1/Z0;

% The port reference impedance must be eta0.
if abs(Z0-eta0) > 0.01
    error("antenna port impedance should be eta0 for now")
else
    Z0 = eta0;
end

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
% Eptr/Hptr: (NTheta x NPhi x NRadius x 3 x NFreq)
% Component order in dim 4: [Eph, Eth, Er]
[Eptr, Rs, theta, phi, ~] = ReadFE_spherical(fileprefix+"_NearField1.efe");
[Hptr, ~,  ~,     ~,   ~] = ReadFE_spherical(fileprefix+"_NearField1.hfe");

% Fall back apothem estimate if mesh data not available in .out file
if isnan(a)
    a = Rs(end) - c0/Freq(1)/2/pi;
end

% Frequency index for single-frequency plots; extraction on the first radius
[~,index_f] = min(abs(Freq-f0));
index_r=1;
index_phi=1;

% Determine symmetry factor: 2 if upper hemisphere only, 1 if full sphere
if max(theta) > (pi/2)
    sym_factor = 1;
else
    sym_factor = 2;
end

k = 2*pi*Freq/c0;
alpha_k=k*Rs(index_r);

numth = length(theta); numph = length(phi); numr = length(Rs);

%% --- Determine Spherical Mode Coefficients (psiTM, psiTE) ---
% Extract the TM coefficients psiTM_l by projecting the radial near field
% Er onto the normalised angular basis Yl0 (only TM is needed here).
% Only odd-l, m = 0 modes (symmetric biconical antenna).

lval = (1:2:lmax);   % Odd mode orders: l = 1, 3, 5, ...
mval = 0;              % Azimuthally symmetric (m=0 only)
L    = length(lval);

% Preallocate angular basis and normalisation
Yl0    = zeros(numth,numph,L);
psiTM  = nan(numf,L); psiTE = psiTM;boundary_norm=psiTM;prop_factor_norm=psiTM;
Norm2  = zeros(1,L);  Normp2 = Norm2;

for lindx = 1:L
    l_factor = lval(lindx)*(lval(lindx)+1);   % l(l+1)

    % Normalised angular basis: Y_l^0(theta) / sqrt(l(l+1))
    Yl0(:,:,lindx) = harmonicY(lval(lindx),0,theta,phi*ones(size(theta)), ...
        'type','real','norm',true,'phase',true) * sqrt(2^(mval~=0)) / sqrt(l_factor);

    % Norm factors (stored for reference)
    Norm2(lindx)  = (2^(mval==0))*pi * 2*factorial(lval(lindx)+mval) / ...
                    (factorial(lval(lindx)-mval)) / (2*lval(lindx)+1);
    Normp2(lindx) = Norm2(lindx)*l_factor;

    for f = 1:numf
        % Radial propagation factor h_l^(2)(kr)/kr / j^(l+1)
        boundary_norm(f,lindx) = 1;
        prop_factor_norm(f,lindx) = sphericalHankel(lval(lindx),2,alpha_k(f)) / boundary_norm(f,lindx) / alpha_k(f) / ((1j)^(lval(lindx)+1));

        % Extract psiTM by projecting Er onto Yl0 and dividing by -jw * prop factor
        psiTM(f,lindx) = (sqrt(eta0)/mu0) * 2*pi * sym_factor * ...
            trapz(theta, Eptr(:,:,index_r,3,f).*(Yl0(:,:,lindx)).*sin(theta), 1) ...
            / prop_factor_norm(f,lindx) / (-1j*2*pi*Freq(f));
    end
    
end

%% --- Reconstruct Near Fields from psiTM ---
% Er, Etheta, Hphi of each TM_l0 mode at every radius Rs(r), and the modal
% wave impedance zTM_l0 = Etheta/Hphi.

Er_modes = zeros(numth,numph,numr,numf,L);Eth_modes=zeros(numth,numph,numr,numf,L);zTM_l0=zeros(numth,numph,numr,numf,L);
Hph_modes = zeros(numth,numph,numr,numf,L);

for f = 1:numf
    for lindx = 1:L
        l_factor = lval(lindx)*(lval(lindx)+1);
        j_factor=((1j)^(lval(lindx)+1));

        for r=1:numr
            
            kr=k(f)*Rs(r);

            Er_modes(:,:,r,f,lindx) = 1j*sqrt(eta0)*l_factor*psiTM(f,lindx)*sphericalHankel(lval(lindx),2,kr)/j_factor ... 
                                        *(-Yl0(:,:,lindx)/Rs(r)) / boundary_norm(f,lindx) ;
            
            Pl1=FastLegendre_Plm(lval(lindx),1,cos(theta));

            Eth_modes(:,:,r,f,lindx) = 1j*sqrt(eta0)*psiTM(f,lindx)*sphericalRicatti_h_prime(lval(lindx),2,kr)/j_factor ...
                                        *(-Pl1/Rs(r))/sqrt(Normp2(lindx)) / boundary_norm(f,lindx) ;

            Hph_modes(:,:,r,f,lindx) = (1/sqrt(eta0))*psiTM(f,lindx)*sphericalRicatti_h(lval(lindx),2,kr)/j_factor ...
                                        *(-Pl1/Rs(r))/sqrt(Normp2(lindx)) / boundary_norm(f,lindx) ;
            
            zTM_l0(:,:,r,f,lindx)=Eth_modes(:,:,r,f,lindx)./Hph_modes(:,:,r,f,lindx) ;
        end
    end
end
clearvars alpha_k prop_factor_norm boundary_norm
% Sum over modes
Er_total = sum(Er_modes, 5);
Eth_total = sum(Eth_modes, 5);
Hph_total = sum(Hph_modes, 5);
%% --- Plot Er Reconstruction Comparison ---
% Er vs theta at the extraction radius and reference frequency:
% reconstructed from psiTM (lines) vs FEKO near field (markers).
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, abs(Er_total(:,index_phi,index_r,index_f)),"LineWidth",LW+3,    "LineStyle","-")
plot(theta/pi, abs(Eptr(:,index_phi,index_r,3,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*", "Color",pointColor)
hold off
ylabel("magnitude (V/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, angle(Er_total(:,index_phi,index_r,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta/pi, angle(Eptr(:,index_phi,index_r,3,index_f))*180/pi,      "LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",pointColor)
hold off
legend("Reconstructed $\psi^{TM}_{l0}$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$E_r$ reconstructed @r="+num2str(Rs(index_r)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $l_{max}$="+num2str(lmax),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% --- Plot Etheta Reconstruction Comparison ---
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, abs(Eth_total(:,index_phi,index_r,index_f)),"LineWidth",LW+3,"LineStyle","-")
plot(theta/pi, abs(Eptr(:,index_phi,index_r,2,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*", "Color",pointColor)
hold off
ylabel("magnitude (V/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, angle(Eth_total(:,index_phi,index_r,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta/pi, angle(Eptr(:,index_phi,index_r,2,index_f))*180/pi,      "LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",pointColor)
hold off
legend("Reconstructed $\psi^{TM}_{l0}$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$E_\theta$ reconstructed @r="+num2str(Rs(index_r)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $l_{max}$="+num2str(lmax),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% --- Plot Hphi Reconstruction Comparison ---
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, abs(Hph_total(:,index_phi,index_r,index_f)),"LineWidth",LW+3,"LineStyle","-")
plot(theta/pi, abs(Hptr(:,index_phi,index_r,1,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*", "Color",pointColor)
hold off
ylabel("magnitude (A/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, angle(Hph_total(:,index_phi,index_r,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta/pi, angle(Hptr(:,index_phi,index_r,1,index_f))*180/pi,      "LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",pointColor)
hold off
legend("Reconstructed $\psi^{TM}_{l0}$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$H_\phi$ reconstructed @r="+num2str(Rs(index_r)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $l_{max}$="+num2str(lmax),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% --- Modal Impedances vs Frequency ---
% zTM_l0 at r = Rs(index_r) and the mid-grid theta point, all l.
th_pnt=ceil(numth/2);
figure(fignum);fignum = fignum+1; 
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(Freq,real(reshape(zTM_l0(th_pnt,index_phi,index_r,:,:),numf,L)),"LineWidth",LW+2, "LineStyle","-")
hold off
ylabel("Real [$\Omega$]","FontSize",FS)

yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(Freq, imag(reshape(zTM_l0(th_pnt,index_phi,index_r,:,:),numf,L)),"LineWidth",LW+1,"LineStyle",":")
hold off
title("$Z^{TM}_{l0}$ @ $\theta/\pi=$"+num2str(theta(th_pnt)/pi,2)+", $r=$" ...
      +num2str(Rs(index_r)*1000,3)+"mm, $l_{max}$="+num2str(lmax),"FontSize",FS)
xlabel("$Freq [Hz]$","FontSize",FS)
ylabel("Imag [$\Omega$]","FontSize",FS)
ylim([-10*eta0 0])
lgd = legend(cellstr(num2str(lval.')),"FontSize",FS); title(lgd,"$l$")
grid on
%% --- Comparison with Equivalent Circuit, l = 1 ---
% TM10 ladder circuit (Chu): series C = eps0*a, shunt L = mu0*a, load eta0.
C10_TM=capacitor(eps0*a);
L10_TM=inductor(mu0*a);
R=resistor(eta0);
TM10=circuit("TM10");
add(TM10,[1,2],C10_TM)
add(TM10,[2;0],L10_TM)
add(TM10,[2;0],R)
setports(TM10,[1;0])

STM10=sparameters(TM10,Freq,eta0);
figure(fignum);fignum = fignum+1; 
plot(Freq,db(squeeze(STM10.Parameters(1,1,:))));grid

Z=s2z(STM10.Parameters,eta0);

figure(fignum);fignum = fignum+1; 
yyaxis left
set(gca,'YColor',"b")
hold on
plot(Freq,real(reshape(zTM_l0(th_pnt,index_phi,index_r,:,1),numf,1)),"LineWidth",LW-1, "LineStyle","-","Color","b")
plot(Freq,real(reshape(Z,numf,1)),"LineWidth",LW, "LineStyle","none","Marker","+","Color","b")
hold off
ylabel("Real [$\Omega$]","FontSize",FS)

yyaxis right
set(gca,'YColor',"r")
hold on
plot(Freq,imag(reshape(zTM_l0(th_pnt,index_phi,index_r,:,1),numf,1)),"LineWidth",LW-1, "LineStyle","-","Color","r")
plot(Freq,imag(reshape(Z,numf,1)),"LineWidth",LW, "LineStyle","none","Marker","x","Color","r")
hold off
ylabel("Imag [$\Omega$]","FontSize",FS)
legend("Real","Real (Circuit)","Imag","Imag (Circuit)")
grid on
%% --- Comparison with Equivalent Circuit, l = 3 ---
% TM30 ladder circuit: series C = eps0*a/3, shunt L = mu0*a/5, series
% C = eps0*a/3, shunt L = mu0*a, load eta0.
TM30=circuit("TM30");
add(TM30,[1,2],capacitor(eps0*a/3))
add(TM30,[2;0],inductor(mu0*a/5))
add(TM30,[2,3],capacitor(eps0*a/3))
add(TM30,[3;0],inductor(mu0*a))
add(TM30,[3;0],resistor(eta0))
setports(TM30,[1;0])

STM30=sparameters(TM30,Freq,eta0);
figure(fignum);fignum = fignum+1; 
plot(Freq,db(squeeze(STM30.Parameters(1,1,:))));grid

Z=s2z(STM30.Parameters,eta0);

figure(fignum);fignum = fignum+1; 
yyaxis left
set(gca,'YColor',"b")
hold on
plot(Freq,real(reshape(zTM_l0(th_pnt,index_phi,index_r,:,2),numf,1)),"LineWidth",LW-1, "LineStyle","-","Color","b")
plot(Freq,real(reshape(Z,numf,1)),"LineWidth",LW, "LineStyle","none","Marker","+","Color","b")
hold off
ylabel("Real [$\Omega$]","FontSize",FS)

yyaxis right
set(gca,'YColor',"r")
hold on
plot(Freq,imag(reshape(zTM_l0(th_pnt,index_phi,index_r,:,2),numf,1)),"LineWidth",LW-1, "LineStyle","-","Color","r")
plot(Freq,imag(reshape(Z,numf,1)),"LineWidth",LW, "LineStyle","none","Marker","x","Color","r")
hold off
ylabel("Imag [$\Omega$]","FontSize",FS)
legend("Real","Real (Circuit)","Imag","Imag (Circuit)")
grid on