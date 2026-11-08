function [jrFlm_phi,jrFlm_theta,jrFlm_r] = Spherical_Mode_Field(l,m,prop_dir,THETA,PHI,kr)
%SPHERICAL_MODE_FIELD  Vector spherical wave function of mode (l,m).
%   [jrFlm_phi,jrFlm_theta,jrFlm_r] = Spherical_Mode_Field(l,m,prop_dir,THETA,PHI,kr)
%   returns the phi, theta and r components of the mode function
%   (the field derived from Spherical_Mode_Potential, used as the
%   projection basis for near-field data):
%       r     : sqrt(l(l+1)) * R_l(kr)/kr * Y_lm
%       theta : R_l'(kr) * (grad Y_lm)_theta
%       phi   : R_l'(kr) * (grad Y_lm)_phi
%   where grad Y_lm comes from gradYlm (includes the 1/sqrt(l(l+1))
%   normalisation) and R_l' is the derivative of R_l with respect to kr.
%   prop_dir selects the radial function R_l(kr), a Riccati-Bessel
%   function scaled by j^(l+1):
%       0 : j^(l+1) * kr j_l(kr)        (standing wave, regular at origin)
%       1 : j^(l+1) * kr h_l^(1)(kr)    (incoming, exp(+jwt) convention)
%       2 : kr h_l^(2)(kr) / j^(l+1)    (outgoing)
%       3 : j^(l+1) * kr y_l(kr)
%   THETA, PHI in radians (same-size arrays); kr = k*r (scalar).
%
%   See also Spherical_Mode_Potential, gradYlm, TesseralHarmonic.

j_l_1=1j.^(l+1);
switch prop_dir
    case 0 
        Rl_kr = sphericalRicatti_j(l,kr).*(j_l_1) ;
        dRl_kr = sphericalRicatti_j_prime(l,kr).*(j_l_1) ;
    case 1
        Rl_kr = sphericalRicatti_h(l,1,kr).*(j_l_1) ;
        dRl_kr = sphericalRicatti_h_prime(l,1,kr).*(j_l_1) ;
    case 2
        Rl_kr = sphericalRicatti_h(l,2,kr)./(j_l_1) ;
        dRl_kr = sphericalRicatti_h_prime(l,2,kr)./(j_l_1) ;
    case 3
        Rl_kr = sphericalRicatti_y(l,kr).*(j_l_1) ;
        dRl_kr = sphericalRicatti_y_prime(l,kr).*(j_l_1) ;
    otherwise
        error("propagation direction should be 0 for standing (at origin), 1 for incoming, or 2 for outgoing")
end

jrFlm_r = sqrt(l*(l+1)).*(Rl_kr/kr).*TesseralHarmonic(l,m,THETA,PHI) ;
[gradYlm_th,gradYlm_ph] = gradYlm(l,m,THETA,PHI,0);
jrFlm_theta = dRl_kr.*gradYlm_th ;   
jrFlm_phi = dRl_kr.*gradYlm_ph ;

end