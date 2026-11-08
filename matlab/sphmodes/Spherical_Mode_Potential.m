function Psi_lm = Spherical_Mode_Potential(l,m,prop_dir,THETA,PHI,kr)
%SPHERICAL_MODE_POTENTIAL  Scalar (Debye) potential of spherical mode (l,m).
%   Psi_lm = Spherical_Mode_Potential(l,m,prop_dir,THETA,PHI,kr) returns
%       R_l(kr) * Y_lm(THETA,PHI) / sqrt(l(l+1))
%   with Y_lm the tesseral harmonic (TesseralHarmonic).
%   prop_dir selects the radial function R_l(kr), a Riccati-Bessel
%   function scaled by j^(l+1):
%       0 : j^(l+1) * kr j_l(kr)        (standing wave, regular at origin)
%       1 : j^(l+1) * kr h_l^(1)(kr)    (incoming, exp(+jwt) convention)
%       2 : kr h_l^(2)(kr) / j^(l+1)    (outgoing)
%       3 : j^(l+1) * kr y_l(kr)
%   THETA, PHI in radians (same-size arrays); kr = k*r (scalar).
%
%   See also Spherical_Mode_Field, TesseralHarmonic.

j_l_1=1j.^(l+1);
switch prop_dir
    case 0 
        Rl_kr = sphericalRicatti_j(l,kr).*(j_l_1) ;
    case 1
        Rl_kr = sphericalRicatti_h(l,1,kr).*(j_l_1) ;
    case 2
        Rl_kr = sphericalRicatti_h(l,2,kr)./(j_l_1) ;
    case 3
        Rl_kr = sphericalRicatti_y(l,kr).*(j_l_1) ;
    otherwise
        error("propagation direction should be 0 for standing (at origin), 1 for incoming, or 2 for outgoing")
end

Psi_lm = Rl_kr.*TesseralHarmonic(l,m,THETA,PHI)/sqrt(l*(l+1));

end