function [Fieldr_mode_output, mode_counter] = NFrFromSph(Q, type, mmax, nmax, num_modes, theta, phi, Ntheta, Nphi, r, freq)
%NFRFROMSPH  Radial near field (E_r or H_r) per mode from FEKO/TICRA Q coefficients.
%   [Fieldr_mode_output, mode_counter] = NFrFromSph(Q, type, mmax, nmax,
%       num_modes, theta, phi, Ntheta, Nphi, r, freq)
%
%   type = "TM" gives E_r (TM modes only); type = "TE" gives H_r (TE modes
%   only). Fieldr_mode_output is (Ntheta x Nphi x num_modes), indexed by the
%   same single index as Q, j = 2*(n(n+1) + m_indx - 1) + s (s = 1 TE,
%   2 TM, m_indx = -m for the GRASP time convention); slots of the other
%   mode type stay zero. Sum over dimension 3 for the total field.
%
%   Same method as FFfromSph (third_party/matlab/davidson), with the
%   far-field radial factor j^(n+1) replaced by omega*h_n^(2)(kr)/kr:
%       E_r = mu0/sqrt(eta0)     * omega h_n^(2)(kr)/kr * SUM n(n+1) Y_nm Q~
%       H_r = 1j*eps0*sqrt(eta0) * omega h_n^(2)(kr)/kr * SUM n(n+1) Y_nm Q~
%   with Y_nm = sigma(m) exp(jm phi)/sqrt(2pi) * Pbar_n^|m|(cos theta)/sqrt(n(n+1))
%   (Pbar from legendre(...,'norm')) and Q~ = sqrt(8pi) conj(Q). The -1j
%   factor of the far-field expression is absorbed into Q~.
%
%   Inputs:
%     Q            (num_modes x 1) Q coefficients, TICRA/GRASP convention
%     type         "TM" or "TE"
%     mmax         maximum azimuthal index (not used)
%     nmax         maximum degree n
%     num_modes    total number of mode slots in Q
%     theta, phi   angles in radians (vectors of length Ntheta, Nphi)
%     r            radius (m);  freq  scalar frequency (Hz)
%   Outputs:
%     Fieldr_mode_output  (Ntheta x Nphi x num_modes) per-mode radial field
%     mode_counter        diagnostic counter (incremented twice per (n,m))
%
%   References:
%     J.E. Hansen (ed.), Spherical Near-Field Antenna Measurements, IET, 2008.
%     FFfromSph.m, D.B. Davidson and A. Sutinjo.

eta0  = 376.730313668; % post-2019 definition
mu0   = 4*pi*1e-7;
c0    = 299792458;
eps0=1/mu0/(c0^2);
omega = 2*pi*freq;
kr    = omega/c0 * r;

% TICRA -> FEKO convention (as in FFfromSph).
Q = sqrt(8*pi)*conj(Q);
if strcmp(type,"TE")
    s=1;
elseif strcmp(type,"TM")
    s=2;
else
    error("seconds argument should be string TE or TM")
end


Fieldr_mode = zeros(Nphi, Ntheta, num_modes);
Fieldr_mode_output = zeros(Ntheta, Nphi, num_modes);
theta=reshape(theta,1,Ntheta);
phi=reshape(phi,Nphi,1);
mode_counter = 0;
for n = 1:nmax
    NP = legendre(n, cos(theta), 'norm');
    NP(end+1,:) = 0;

    radial_factor = omega * sphericalHankel(n, 2, kr) / kr;

    for m = -n:n
        m_indx = -m; % swap due to GRASP time convention

        if m >= 0
            Sign = 1;
        else
            Sign = (-1)^m;
        end

        mode_counter = mode_counter + 1;
        j_ = 2*(n*(n+1) + m_indx - 1) + s;
        phi_const = (exp(1j*m*phi) * Sign / sqrt(2*pi) / sqrt(n*(n+1))) ;
        Fieldr_mode(:,:,j_) = n*(n+1) * radial_factor * phi_const * NP(abs(m)+1, :) * Q(j_);
        
        mode_counter = mode_counter + 1;
    end
end

% Overall scale: mu0/sqrt(eta0) for E_r (TM), 1j*eps0*sqrt(eta0) for H_r (TE).
if s==2
    conversion = mu0 / sqrt(eta0) ;
elseif s==1
    conversion = 1j * eps0 * sqrt(eta0) ;
else
    error("TE/TM type not configured correctly")
end

Fieldr_mode_output = conversion * permute(Fieldr_mode,[2 1 3 4]);
end