%% conical_test.m
% Read a FEKO spherical-mode (.sph) export of a dipole, rebuild its far
% field at 1 GHz with FFfromSph and plot |E_theta| vs theta, then plot
% the odd-l TM_l0 coefficients vs frequency and split Q into TM/TE (l,m) arrays.
% Used to inspect the m = 0 TM mode content that the analytic biconical
% (conical harmonic) model in this folder describes.
%
% Requires: <filepath>/<flnm>, a FEKO .sph export.
% FFfromSph and ReadSphModes_freqs are from third_party/matlab/davidson.

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
filepath = "C:\path\to\FEKO\dipole\";      % folder containing the FEKO .sph export (keep trailing \)
flnm     = "dipole_cylinder_FarField1.sph"; % .sph file name
% =========================

[Q, j_, j3D, mmax, lmax, freqs] = ReadSphModes_freqs(fullfile(filepath, flnm)) ;

angle_incr = 1 % angular step (deg) for theta and phi
theta_max=180;  
thetad = 0:angle_incr:theta_max;
theta = deg2rad(thetad);
Ntheta= length(theta);
phi_max=360; 
phid =  0:angle_incr:phi_max;
phi = deg2rad(phid);
Nphi= length(phi);

num_modes=2 *(lmax*(lmax+1)+ mmax-1) + 2;

f1=1e9;
[~,index]=min(abs(freqs-f1))
fval=freqs(index)

tstart=tic;
[E_th,E_ph,mode_counter]= FFfromSph(Q(:,index),mmax,lmax,num_modes,theta,phi,Ntheta,Nphi);
tstop=toc(tstart)

% |E_theta| vs theta in the first phi cut (E_th is Nphi x Ntheta).
figure
plot(thetad,abs(E_th(1,:)))


[Q_TM_m0, j_TM_m0, j3D_TM_m0] = ExtractTM_M0Modes(Q, lmax, j3D);

% |Q| of the odd-l TM_l0 modes (l = 1,3,...,11; row l of Q_TM_m0) vs frequency.
figure
plot(freqs/1e6, abs(Q_TM_m0(1:2:11,:)))
xlabel('Frequency (MHz)'); ylabel('|Q_{TM,l,0}|')

[Q_TM, Q_TE] = ExtractLMmodes(Q, lmax);
