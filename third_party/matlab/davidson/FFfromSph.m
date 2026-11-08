function [E_th,E_ph,mode_counter]=FFfromSph(Q,mmax,nmax,num_modes,theta,phi,Ntheta,Nphi)

% This function computes the far fields E_theta and E_phi from the spherical modal coefficients Q, 
% stored in single index format,  suppressing the e{-jkr}/r term as is usual in FEKO.
% The fields are stored in phi,theta index format.
% mmax is the maximum azimuthal (phi) mode number
% nmax is the maximum elevation (theta) mode number
% theta and phi are 1D arrays  with the values of the angles in RADIANS, of
% length Ntheta and Nphi respectively. 

% The function uses the built-in MATLAB special function LEGENDRE for the associated
% Legendre functions.
%
% References:
% [1]. J.E. Hansen (ed), "Spherical Near-Field Antenna
% Measurements",IET, 2008 reprinting.
% [2]. FEKO user manual, Altair Engineering.
% [3]. GRASP user manual, TICRA.
% [4]. A. Sutinjo, "Calculating Far-Field Radiation Based on FEKO
% Spherical Wave Coefficients", Curtin University, Unpublished draft,
% August 2015.
% [5]. D.B.Davidson, A. Sutinjo, 
% "Efficient storage of embedded element patterns for low frequency radio telescopes"
% Curtin University, Submitted to PASA, June 2024.
%
% Note that the spherical modes are stored in TICRA convention,
% assuming the time convention of [1], viz. e^{-i omega t}.
% These must be multiplied by \sqrt(8*pi) [2] and conjugated.
% Note that m must also be interchanged with -m due to the time convention.
%
% The above scaling of the modes still leaves a further scaling term
% of sqrt(eta0/(2*pi)) required to fully reconstruct the FEKO field.
%
% Author: David B Davidson. 30 May 2024. 
%   Sign corrected for CmnConstant, 4 June 2024.
%   indx corrected for theta=pi (was incorrectly overwriting theta=0 value).


eta0 = 376.730313668; % post-2019 definition.

% Scale Q to align with FEKO usage.
% Also complex conjugate the coefficients to align with FEKO
% e^{j omega t} time convention.
Q = sqrt(8*pi)*conj(Q);

% Preallocate storage for speed. 
E_th_mode = zeros(Nphi,Ntheta,num_modes);
E_ph_mode = zeros(Nphi,Ntheta,num_modes);

mode_counter=0;
for n=1:nmax
    NP = legendre(n,cos(theta),'norm');
    NP(end+1,:) = 0; % to return the m+1 mode as zero.
    for m=-n:n
        m_indx = -m; % The m-mode swaps due to the GRASP time convention
        if m >= 0
            Sign = 1;
        else
            Sign = (-1)^m;
        end
        CmnConstant = sqrt((n+abs(m)+1)*(n-abs(m))); % Note - defined as positive in this code.
        %Precompute m NP(x)/sin(x) and abs(m) NP(x)/sin(x) terms, including for special cases.
        NPdsm    = NP(abs(m)+1,:)./sin(theta)*m; 
        NPdsabsm = NP(abs(m)+1,:)./sin(theta)*abs(m);
        indx0    = find(theta==0);
        if ~isempty(indx0) % special case for theta=0. See Hansen, eq A1.39
            if abs(m) ~=1
                NPdsm(indx0)    = 0;
                NPdsabsm(indx0) = 0;
            else
                Cmn = sqrt((2*n+1)/2*factorial(n-abs(m))/factorial(n+abs(m)));
                NPdsm(indx0)    = Cmn*sign(m)*n*(n+1)/2;
                NPdsabsm(indx0) = Cmn*n*(n+1)/2; % no sign function since abs(m) is always non-negative.
            end
        end
        indx1    = find(theta==pi);
        if ~isempty(indx1) % special case for theta=180 deg. See Hansen, eq A1.41
            if abs(m) ~=1
                NPdsm(indx1)    = 0;
                NPdsabsm(indx1) = 0;
            else
                Cmn = sqrt((2*n+1)/2*factorial(n-abs(m))/factorial(n+abs(m)));
                NPdsm(indx1)    = (-1)^(n+1)*Cmn*sign(m)*n*(n+1)/2;
                NPdsabsm(indx1) = (-1)^(n+1)*Cmn*n*(n+1)/2; % ditto
            end
        end
 
        % TE modes 
        s=1; 
        j_  = 2 *(n*(n+1)+ m_indx-1) + s;
        for tt = 1:Ntheta 
            % Vectorised over phi; looped over theta.
            phi_const = (exp(1i*m*phi)*Sign/sqrt(n*(n+1))).'; % Transpose for dimensional consistency below.
            E_th_mode(:,tt,j_) = -(1j^n)    *  NPdsm(tt)*Q(j_);
            E_th_mode(:,tt,j_) = phi_const.*E_th_mode(:,tt,j_);
            E_ph_mode(:,tt,j_) = -(1j^(n+1))*( NPdsabsm(tt)*Q(j_)*cos(theta(tt))...
                -CmnConstant*Q(j_)*NP(abs(m)+2,tt) );
            E_ph_mode(:,tt,j_) = phi_const.*E_ph_mode(:,tt,j_);
        end
        mode_counter = mode_counter +1;
        
        % NB! Do NOT reset fields or mode counter for s=2 case. This is a different mode set.

        % TM modes
        s=2; 
        j_  = 2 *(n*(n+1)+ m_indx-1) + s;
        for tt = 1:Ntheta
            % Vectorised over phi; looped over theta.
            phi_const = (exp(1i*m*phi)*Sign/sqrt(n*(n+1))).'; % Transpose for dimensional consistency below.
            E_th_mode(:,tt,j_) = 1i^n*    (NPdsabsm(tt)*Q(j_)*cos(theta(tt))...
                -CmnConstant*Q(j_)*NP(abs(m)+2,tt));
            E_th_mode(:,tt,j_) = phi_const.*E_th_mode(:,tt,j_) ;
            E_ph_mode(:,tt,j_) = 1i^(n+1)* NPdsm(tt) *Q(j_);
            E_ph_mode(:,tt,j_) = phi_const.*E_ph_mode(:,tt,j_);
        end
        mode_counter = mode_counter +1;
    end
end

% Now sum over all modes
E_th=sum(E_th_mode,3);
E_ph=sum(E_ph_mode,3);
% Apply scaling needed by FEKO.
sqrt_fac=sqrt(eta0/(2*pi)); 
E_th=sqrt_fac*E_th;
E_ph=sqrt_fac*E_ph;
end