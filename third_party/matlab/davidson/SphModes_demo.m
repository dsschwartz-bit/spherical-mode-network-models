% This script demonstrates the use of a spherical mode reconstruction function to
% compute far fields from a TICRA/GRASP-standard .sph file exported by
% FEKO.
%
% It firstly calls a function to read in a spherical mode file for the
% fields stored in file flnm, and then calls FFFROMSPH to compute
% the far fields. See that function for references and documentation.
%
% The script offers several standard tests using canonical radiators.
% An option is provided to limit the number of modes used - by default,
% all the stored modes in the .sph file are used.
%
% Angle theta is elevation angle (measured from z-axis downwards).
% Angle phi is azimuth angle (measured from x-axis towards y-axis).
% These are the usual antenna engineering angle conventions.
%
% The .sph files and FEKO pattern scripts used by the test cases below are
% not included in this repository.

% Author: DB Davidson, June 2024. 

clearvars
close all

angle_incr = 1 % This is the angular increment, in degrees, for both theta and phi.
theta_max=90;  
thetad = 0:angle_incr:theta_max;
theta = deg2rad(thetad);
Ntheta= length(theta);
phi_max=360; 
phid =  0:angle_incr:phi_max;
phi = deg2rad(phid);
Nphi= length(phi);



test_case = input(['Enter which test case,  1/2 wave dipole(1-default)',...
    'Hertian z 299 MHz (2) ditto x (3) ditto y (4) ditto xy (5) ditto VED array (6) ditto HED array (7)'])
if (isempty(test_case))
    test_case = 1 %
end

switch test_case
    case 1
        flnm=('dipole_FarField1_299MHz.sph')
        freq=299.792e6 % freq in Hz;
        dip_th =  8.311E-01*exp(1i*deg2rad(98.01));% FEKO
        % 8.317E-01 *exp(1i*deg2rad(97.48)); % FEKO 300 MHz
        dip_th_CIRA = -0.1083 + 0.8240i; % From D Ung codes.

    case 2
        flnm=('hertzian_dipole_FarField1_299MHz.sph')
        freq=299.792e6 % freq in Hz;
        Hertz_dip_th = 1.884E+02*exp(1i*deg2rad(90.00)) % 299.792 MHz
    case 3
        flnm=('hertzian_x_dipole_FarField1_299MHz.sph')
        freq=299.792e6 % freq in Hz;
        Hertz_x_dip_th = 1.884E+02*exp(1i*deg2rad(-90.00)) %  MHz
        % this is at phi=0, theta=0
        Hertz_x_dip_ph = 1.884E+02*exp(1i*deg2rad(90.00)) %  MHz
        % this is at phi=90, theta=90
    case 4
        flnm=('hertzian_y_dipole_FarField1_299MHz.sph')
        freq=299.792e6 % freq in Hz;
        Hertz_y_dip_ph = 1.884E+02*exp(1i*deg2rad(-90.00)) %  MHz
        % this is at phi=0, , theta=90
    case 5
        flnm=('hertzian_xy_dipole_FarField1_299MHz.sph')
        freq=299.792e6 % freq in Hz;
        Hertz_y_dip_ph = 1.884E+02*exp(1i*deg2rad(90.00)) %  MHz
        % this is at phi=135 deg, , theta=90
    case 6
        flnm=('hertzian_z_dip_array_FarField1_299MHz.sph')
        freq=299.792e6 % freq in Hz;
        hertzian_z_dip_array_xyFarField1
        Hertz_array_H = Hplane;
        hertzian_z_dip_array_yzFarField1
        Hertz_array_E = Eplane;
        % change above theta defaults for this case
        theta_max=180;% 181;
        thetad = 0:angle_incr:theta_max;
        theta = deg2rad(thetad);
        Ntheta= length(theta);
  case 7
        flnm=('hertzian_x_dip_array_FarField2_299MHz.sph')
        freq=299.792e6 % freq in Hz;
        hertzian_x_dip_array_yzFarField
        Hertz_array_H = Hplane;
        theta_max=180;% 181;
        thetad = 0:angle_incr:theta_max;
        theta = deg2rad(thetad);
        Ntheta= length(theta);
end

[Q,j,j3D,mmax,nmax] = ReadSphModes(flnm);

mode_limit = input('Limit number of modes? Enter number or hit enter to skip')
if (~isempty(mode_limit))
    nmax = mode_limit
    mmax = nmax
    num_modes=2 *(nmax*(nmax+1)+ mmax-1) + 2;
end


c = physconst('LightSpeed'); % m/s

lambda= c/freq;
k = 2*pi/lambda;



% Preallocate storage for speed. Need to know number of modes.
num_modes=2 *(nmax*(nmax+1)+ mmax-1) + 2;
E_th_mode = zeros(Nphi,Ntheta,num_modes);
E_ph_mode = zeros(Nphi,Ntheta,num_modes);


% Compute the far fields from the spherical modes.
tstart=tic;
[E_th,E_ph,mode_counter]= FFfromSph(Q,mmax,nmax,num_modes,theta,phi,Ntheta,Nphi);
tstop=toc(tstart)

% Output the results, comparing to other data - primarily fields directly
% computed by FEKO, or for SKA, to other CIRA scripts.
switch test_case
    case {1}
        % Compare value for E_theta (co-pol) at phi-0, theta=90 deg
        dip_th
        E_th(1,end)
    case {2}
        % Compare value for E_theta (co-pol) at phi-0, theta=90 deg.
        Hertz_dip_th
        E_th(1,end)
    case {3}
        % Compare value for E_theta (co-pol) at phi-0, theta=0 deg.
        Hertz_x_dip_th
        E_th(1,1)
        % Compare value for E_phi (co-pol) at phi-90, theta=90 deg.
        Hertz_x_dip_ph
        E_ph(91,end)

        for ii=1:length(phi)
            for jj=1:length(theta)
                E_th_anl(ii,jj) = -1j*188.4*cosd(thetad(jj))*cosd(phid(ii));
                E_ph_anl(ii,jj) = 1j*188.4*sind(phid(ii));
            end
        end

        figure
        plot(thetad,abs(E_th(1,:)),thetad(1:5:end),abs(E_th_anl(1,1:5:end)),'+',thetad,abs(E_ph(1,:)),thetad(1:5:end),abs(E_ph_anl(1,1,1:5:end)),'x')
        title('x-directed Hertian dipole E-plane')
        legend('E_\theta this code','E_\theta analytical','E_\phi this code','E_\phi analytical')
        figure
        plot(thetad,abs(E_th(91,:)),thetad(1:5:end),abs(E_th_anl(91,1:5:end)),'+',thetad,abs(E_ph(91,:)),thetad(1:5:end),abs(E_ph_anl(91,1,1:5:end)),'x')
        title('x-directed Hertian dipole H-plane')
        legend('E_\theta this code','E_\theta analytical','E_\phi this code','E_\phi analytical')

    case {4}
        % Compare value for E_phi (co-pol) at phi-0, theta=90 deg.
        Hertz_y_dip_ph
        E_ph(1,end)
    case {5}
        % Compare value for E_phi (co-pol) at phi-0, theta=90 deg.
        Hertz_y_dip_ph
        E_ph(136,end)
    case {6}
        figure
        plot(phid,abs(E_th(:,90)),'-',...
            Hertz_array_H(:,2),abs(Hertz_array_H(:,4)),'--')
        title('Spherical mode reconstruction test - 2 element Hertzian dipole array')
        xlabel('\phi (degrees)')
        ylabel('Field mag. excl. far field term  (V)')
        legend('H plane','H plane FEKO')
        axis([0 360 0 400])
        print -dpdf SphModeHplane
        print -dpng SphModeHplane

        figure
        plot(thetad,abs(E_th(91,:)),'-',...
            Hertz_array_E(91:end,1),abs(Hertz_array_E(91:end,4)),'--')
        title('Spherical mode reconstruction test - 2 element Hertzian dipole array')
        xlabel('\theta (degrees)')
        ylabel('Field mag. excl. far field term  (V)')
        legend('E plane','E plane FEKO')
        axis([0 180 0 400])
        print -dpdf SphModeEplane
        print -dpng SphModeEplane
    case {7}
        figure
        plot(Hertz_array_H(:,1),abs(Hertz_array_H(:,5)),'k--',...
            thetad,abs(E_ph(91,:)),'k-',flip(-thetad),abs(flip(E_ph(271,:))),'k-')
        title('Spherical mode reconstruction test - 2 element Hertzian dipole array')
        xlabel('\theta (degrees)')
        ylabel('E_\phi field mag. excl. far field term  (V)')
        legend('H plane FEKO','H plane')
        axis([-180 180 0 400])
        print -dpdf SphModeHplane_x_dip_array
        print -dpng SphModeHplane_x_dip_array 

        figure
        plot(Hertz_array_H(:,1),abs(Hertz_array_H(:,3)),'k--',...
            thetad,abs(E_th(91,:)),'k-',flip(-thetad),abs(flip(E_th(271,:))),'k-')
        title('Spherical mode reconstruction test - 2 element Hertzian dipole array')
        xlabel('\theta (degrees)')
        ylabel('E_\theta field mag. excl. far field term  (V)')
        legend('H plane FEKO','H plane')
        axis([-180 180 0 400])
        print -dpdf SphModeHplane_x_dip_array
        print -dpng SphModeHplane_x_dip_array 

end

