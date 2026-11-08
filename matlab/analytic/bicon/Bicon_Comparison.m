%% Bicon_Comparison.m
% Compare the TM_l0 spherical-mode coefficients of a FEKO biconical
% antenna model with the analytic biconical solution (Bicon_Coeff_Symbolic).
% psiTM is extracted by projecting the FEKO near-field E_r onto the
% m = 0, odd-l harmonics; E_r is then rebuilt from psiTM and from the
% analytic coefficients (scaled to the same peak) and compared with the
% FEKO field, and both coefficient sets are plotted vs frequency.
%
% Requires, in the FEKO project folder: <antenna_name>.out and
% <antenna_name>_NearField1.efe/.hfe (port reference impedance eta0).

clear
close all
clc

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
data_dir     = "C:\path\to\FEKO\bicon";    % FEKO project folder with the exports
antenna_name = "bicon_edgeport_symmetry";    % FEKO model name (file prefix)
r_index      = 2;      % near-field sphere used (index into Rs); [] = the one nearest a + lambda_max/(2pi)
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

%% --- Simulation Parameters ---
lmax = 31;              % Maximum spherical mode order
f0   = 5e9;             % Reference frequency for single-frequency plots (Hz)
fmin = 0.3e9;           % Minimum frequency for plot x-axis (Hz)
fmax = 5e9;             % Maximum frequency for plot x-axis (Hz)

%% --- File Paths ---
fileprefix   = fullfile(data_dir, antenna_name);

%% --- Import FEKO Port Parameters ---
% V0: source voltage, Z0: port reference impedance, R0: port location,
% has_gnd_plane: PEC ground-plane flag, even_odd_excitation: symmetry flag,
% a: antenna apothem (NaN if not in the .out file), Zant: input impedance.
[V0, Z0, R0, has_gnd_plane, even_odd_excitation, a,Zant,Freq] = ReadFEKOout(fileprefix+".out");
Y0 = 1/Z0;

% The analytic comparison assumes a port reference impedance of eta0.
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
% Eptr/Hptr: (numth x numph x numr x 3 x numf); dim 4 = [E_phi, E_theta, E_r]
[Eptr, Rs, theta, phi, ~] = ReadFE_spherical(fileprefix+"_NearField1.efe");
[Hptr, ~,  ~,     ~,   ~] = ReadFE_spherical(fileprefix+"_NearField1.hfe");

% Apothem estimate if not available from the .out file
if isnan(a)
    a = Rs(end) - c0/Freq(1)/2/pi;
end

% Frequency and near-field sphere indices
[~,index_f] = min(abs(Freq-f0));
if isempty(r_index)
    [~,index_r] = min(abs(a+c0/2/pi/min(Freq)-Rs),[],'all');
else
    index_r = r_index;
end
index_phi=1;

% Symmetry factor: 2 if only the upper hemisphere was exported, else 1
if max(theta) > (pi/2)
    sym_factor = 1;
else
    sym_factor = 2;
end

k = 2*pi*Freq/c0;
alpha_k=k*Rs(index_r);

numth = length(theta); numph = length(phi); numr = length(Rs);

%% --- Determine Spherical Mode Coefficients (psiTM, psiTE) ---
% Extract the TM Hertz-potential coefficients psi_l by projecting the
% radial near field E_r onto the normalised basis Y_l0/sqrt(l(l+1)).
% Only odd-l, m = 0 modes (symmetric biconical antenna).

lval = (1:2:lmax);   % l = 1, 3, 5, ...
mval = 0;            % m = 0 only
L    = length(lval);
l_factor=lval.*(lval+1);

Yl0    = zeros(numth,numph,L);
psiTM  = nan(numf,L); psiTE = psiTM;boundary_norm=psiTM;prop_factor_norm=psiTM;
Norm2  = zeros(1,L);  Normp2 = Norm2;

for lindx = 1:L

    % Y_l^0(theta) / sqrt(l(l+1)) (harmonicY, third_party)
    Yl0(:,:,lindx) = harmonicY(lval(lindx),0,theta,phi*ones(size(theta)), ...
        'type','real','norm',true,'phase',true) * sqrt(2^(mval~=0)) / sqrt(l_factor(lindx));

    % Norm factors; Normp2 is used in the analytic conversion below
    Norm2(lindx)  = (2^(mval==0))*pi * 2*factorial(lval(lindx)+mval) / ...
                    (factorial(lval(lindx)-mval)) / (2*lval(lindx)+1);
    Normp2(lindx) = Norm2(lindx)*l_factor(lindx);

    for f = 1:numf
        % Radial factor h_l^(2)(kr)/(kr j^(l+1)) at the extraction radius
        boundary_norm(f,lindx)=1;

        prop_factor_norm(f,lindx) = sphericalHankel(lval(lindx),2,alpha_k(f)) / boundary_norm(f,lindx) / alpha_k(f) / ((1j)^(lval(lindx)+1));

        % psiTM: projection of E_r onto Yl0, divided by -jw and the radial factor
        psiTM(f,lindx) = (sqrt(eta0)/mu0) * 2*pi * sym_factor * ...
            trapz(theta, Eptr(:,:,index_r,3,f).*(Yl0(:,:,lindx)).*sin(theta), 1) ...
            / prop_factor_norm(f,lindx) / (-1j*2*pi*Freq(f));
    end
    
end

%% --- Analytic Biconical Antenna Coefficients ---
% Analytic TM coefficients of an ideal biconical antenna of half-angle
% theta0 (Bicon_Coeff_Symbolic, eigenvalue search up to nu = 200, 1 %
% switch-over threshold, raw coefficients). c_l: (numf x L);
% nu_n, dnu_n: conical eigenvalues and d(nu_n)/d(theta0).
theta0 = pi/4;   % cone half-angle (rad)
converted=false;
[c_l, nu_n, dnu_n, nu_approx_idx, dnu_approx_idx] = Bicon_Coeff_Symbolic(a, theta0, Freq, lmax,200,1,converted);fignum=fignum+1;

conv_factor=zeros(numf,L);

for f=1:numf
    if ~converted
        conv_factor(f,:)  = sqrt(eta0)/(2*pi) * (1j.^(lval+1)) .* sqrt(Normp2) ./ sphericalRicatti_h(lval, 2, k(f)*a);
    else 
        conv_factor=1;
    end
end

c_l_conv = c_l .* conv_factor .* boundary_norm ;

%% --- Reconstruct Exterior Near-Field Er from Coefficients ---
% E_r at the extraction radius from psiTM (FEKO) and from c_l_conv
% (analytic), with the same Hertz-potential formula.

Er_l  = zeros(numth,numph,numf,L); Hr_modes = Er_l;
Er_l_anlytc = Er_l;

for f = 1:numf
    for lindx = 1:L
        % From the extracted psiTM
        Er_l(:,:,f,lindx) = (-1j*2*pi*Freq(f))*l_factor(lindx)*(mu0/sqrt(eta0)) * ...
                                  psiTM(f,lindx)*prop_factor_norm(f,lindx)*(Yl0(:,:,lindx));

        % From the analytic coefficients
        Er_l_anlytc(:,:,f,lindx) = (-1j*2*pi*Freq(f))*l_factor(lindx)*(mu0/sqrt(eta0)) * ...
                                         c_l_conv(f,lindx)*prop_factor_norm(f,lindx)*(Yl0(:,:,lindx));

    end
end
clearvars alpha_k prop_factor_norm

% Sum over modes
Er_exterior        = sum(Er_l(:,:,:,:), 4);
Er_exterior_anlytc = sum(Er_l_anlytc(:,:,:,:), 4);

% Scale the analytic result to the extracted peak (shape comparison)
scale_factor = reshape(max(Er_exterior,[],[1,2]),numf,1)./reshape(max(Er_exterior_anlytc,[],[1,2]),numf,1) ;
if isnan(scale_factor)
    scale_factor=ones(numf,1);
end
%% --- Plot Er Reconstruction Comparison ---
% E_r vs theta at the extraction radius and f0: from psiTM, from the
% analytic coefficients, and the FEKO field (white markers).
figure(fignum);fignum = fignum+1; 

% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, abs(Er_exterior(:,index_phi,index_f)),"LineWidth",LW+3,    "LineStyle","-")
plot(theta/pi, abs(Er_exterior_anlytc(:,index_phi,index_f).*scale_factor(index_f)), "LineWidth",LW+1,"LineStyle","-")
plot(theta/pi, abs(Eptr(:,index_phi,index_r,3,index_f)),"LineWidth",LW-1.75, "LineStyle","none", "Marker","*", "Color",[1 1 1])
hold off
ylabel("magnitude (V/m)","FontSize",FS)

% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
hold on
plot(theta/pi, angle(Er_exterior(:,index_phi,index_f))*180/pi,"LineWidth",LW+2,"LineStyle",":")
plot(theta/pi, angle(Er_exterior_anlytc(:,index_phi,index_f)*scale_factor(index_f))*180/pi,"LineWidth",LW,"LineStyle",":")
plot(theta/pi, angle(Eptr(:,1,index_r,3,index_f))*180/pi,      "LineWidth",LW-1.75, "LineStyle","none", "Marker","^", "Color",[1 1 1])
hold off
legend("Reconstructed $\psi^{TM}_{l0}$","Reconstructed Analytic $c_l$","Computed (FEKO)","FontSize",FS)
ylim([-180 180])
title("$E_r$ reconstructed @r="+num2str(Rs(index_r)*1000,3)+"mm, "+ ...
      num2str(Freq(index_f)/1e9)+"GHz, $l_{max}$="+num2str(lmax),"FontSize",FS)
xlabel("$\theta/\pi$","FontSize",FS)
ylabel("phase (deg)","FontSize",FS)
grid on
%% Plot Coefficients
% Magnitudes: psiTM (left axis, solid) and scaled analytic (right, dotted).
ymax=max(abs(2*psiTM),[],"all");
dynamic_range=ymax; %mag
ymin=ymax-dynamic_range;

figure(fignum);fignum = fignum+1;
yyaxis left
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
plot(Freq, abs(psiTM(:,1:L)), "LineWidth",LW, "LineStyle","-")
ylabel("$dB|\psi^{TM}_{l0}|$","FontSize",FS)
ylim([ymin,ymax])
yyaxis right
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
set(gca,'ColorOrderIndex',1)
plot(Freq, abs(c_l_conv(:,1:L).*scale_factor), "LineWidth",LW, "LineStyle",":")
ylabel("$dB|C^{TM}_{l0}|$","FontSize",FS)
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
plot(Freq, angle(psiTM(:,1:L)./scale_factor)*180/pi, "LineWidth",LW, "LineStyle","-")
ylabel("$phase \psi^{TM}_{l0}$","FontSize",FS)
ylim([-180 180]);
yyaxis right
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'YColor',get(groot,'defaultAxesXColor'))
set(gca,'ColorOrderIndex',1)
plot(Freq, angle(c_l_conv(:,1:L))*180/pi, "LineWidth",LW, "LineStyle",":")
ylabel("$phase C^{TM}_{l0}$","FontSize",FS)
lgd = legend(cellstr(num2str(reshape(lval,L,1))),"FontSize",FS); title(lgd,"$l$")
title("Phase Simulated and Theoretical Coefficients","FontSize",FS)
xlim([fmin,fmax]); xlabel("Frequency (Hz)","FontSize",FS)
ylim([-180 180]); grid on