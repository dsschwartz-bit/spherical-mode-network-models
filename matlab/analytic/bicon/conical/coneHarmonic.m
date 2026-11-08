function Lnu = coneHarmonic(nu,theta,evenORodd)
%CONEHARMONIC  Conical harmonic L_nu(theta) for the biconical antenna region.
%   Lnu = coneHarmonic(nu,theta,evenORodd) returns, for degree nu (real,
%   non-integer, numeric or symbolic) and polar angle theta (radians):
%       "odd"  (default): L_nu = P_nu(cos theta) - P_nu(-cos theta)
%       "even":           L_nu = P_nu(cos theta) + P_nu(-cos theta)
%   with P_nu = Legendre_Pnu_m(nu,0,.). The odd combination is
%   antisymmetric about theta = pi/2 and, unlike P_nu alone, finite on
%   both cones theta0 and pi - theta0; it equals 2*M_nu(cos theta), where
%       M_nu(x) = (2/sqrt(pi)) Gamma(nu/2+1)/Gamma((nu+1)/2) sin(pi nu/2)
%                 * x * 2F1((1-nu)/2, nu/2+1; 3/2; x^2).
%   The TM modes of the region between two cones of half-angle theta0
%   (the analytic biconical antenna, Bicon_*.m) have non-integer degrees
%   nu_n given by the roots of coneHarmonic(nu,theta0,"odd") = 0.
%   Any other evenORodd value leaves Lnu unassigned (error).
%
%   See also Legendre_Pnu_m, Bicon_Coeff_Symbolic.
if nargin<3||isempty(evenORodd)
    evenORodd="odd";
end
if strcmp(evenORodd,"odd")
    Lnu = Legendre_Pnu_m(nu,0,cos(theta))-Legendre_Pnu_m(nu,0,-cos(theta)) ;
elseif strcmp(evenORodd,"even")
    Lnu = Legendre_Pnu_m(nu,0,cos(theta))+Legendre_Pnu_m(nu,0,-cos(theta)) ;
end
end
