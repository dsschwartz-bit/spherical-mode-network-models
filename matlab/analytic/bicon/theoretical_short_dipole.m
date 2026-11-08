%% theoretical_short_dipole.m
% Near-field E_r and E_theta of a short (Hertzian) dipole of length 2a
% fed with current V0*Yant/2, and a plot of E_r (magnitude and phase) vs
% theta. Script fragment: run it inside a workspace that already defines
% V0, Yant, a, Rs, Freq, k, theta, eta0, eps0, mu0, fignum, index0_f,
% index0_r, LW and FS (e.g. after a Bicon_*.m script). Er_smalldip and
% Eth_smalldip are evaluated with implicit expansion over theta, Freq/k
% and Rs; only Er is plotted. exp(+jwt) convention.
Er_smalldip = (1/2)*(2*V0*Yant)*(2*a).*(eta0./Rs + 1./(1j*2*pi*Freq*eps0.*(Rs.^2)) ) .*exp(-1j*k.*Rs)./(4*pi*Rs).*cos(theta);
Eth_smalldip = (1/2)*(V0*Yant)*(2*a).*(1j*2*pi*Freq*mu0 + eta0./Rs + 1./(1j*2*pi*Freq*eps0.*(Rs.^2)) ) .*exp(-1j*k.*Rs)./(4*pi*Rs).*sin(theta);
fignum=fignum+1; figure_E=figure(fignum);
% Left axis: magnitude
yyaxis left
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
plot(theta, abs(Er_smalldip(:,index0_f,index0_r)), "LineWidth",LW, "LineStyle","-")
ylabel("magnitude (V/m)","FontSize",FS)
title("$E_r$ theoretical and computed @r="+num2str((Rs(index0_r)*1000),3)+"mm, "+num2str(Freq(index0_f)/1e9)+"GHz","FontSize",FS)
% Right axis: phase
yyaxis right
set(gca,'YColor',get(groot,'defaultAxesXColor'))
colororder(gca, get(groot,'defaultAxesColorOrder'))
set(gca,'ColorOrderIndex',1)
plot(theta, angle(Er_smalldip(:,index0_f,index0_r))*180/pi, "LineWidth",LW, "LineStyle",":")
ylabel("phase (deg)","FontSize",FS)
ylim([-180 180])
xlabel("$\theta/\pi$","FontSize",FS)
grid