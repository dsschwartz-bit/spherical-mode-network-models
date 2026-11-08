%TESTING_CIRCUITS  Chu equivalent circuits of the spherical-mode wave impedance.
%   Builds the RF Toolbox ladder circuit of the TM_l mode on a sphere of
%   radius a and compares its input impedance with the exact TM wave
%   impedance  Z = j*eta0*(x h_l(x))'/(x h_l(x)),  x = k*a, h_l = h_l^(2).
%   Chu's continued fraction gives the ladder, from the input:
%       series C = eps0*a/l, shunt L = mu0*a/(2l-1), series C = eps0*a/(2l-3),
%       ... down to coefficient 1, then the load eta0 (series after a C,
%       parallel after an L). l = 1 is the familiar series C, shunt L || eta0.
%   Requires RF Toolbox. Uses the workspace variable Freq (Hz) if it exists.

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'setup.m'));   % adds matlab/ and third_party/

%% --- Settings ---
l = 3;                  % mode order
a = 0.0377;             % sphere radius (m)
if ~exist('Freq', 'var'), Freq = linspace(50e6, 10e9, 400); end

%% --- Physical Constants ---
c0   = physconst('lightspeed');
mu0  = 4*pi*1e-7;       % Permeability of free space (H/m)
eps0 = 1/c0^2/mu0;      % Permittivity of free space (F/m)
eta0 = sqrt(mu0/eps0);  % Characteristic impedance of free space (Ohm)

%% --- Ladder circuit (nodes 1..N along the series branch, 0 = ground) ---
coeff = [l, 2*l-1:-2:1];   % impedance, admittance, impedance, ... levels
ckt = circuit("TM" + l);
node = 1;
for k = 1:numel(coeff)
    if mod(k, 2) == 1      % series capacitor
        add(ckt, [node, node+1], capacitor(eps0*a/coeff(k)));
        node = node + 1;
    else                   % shunt inductor
        add(ckt, [node, 0], inductor(mu0*a/coeff(k)));
    end
end
add(ckt, [node, 0], resistor(eta0));   % load: in series after a C, in parallel with a final L
setports(ckt, [1 0]);

S = sparameters(ckt, Freq, eta0);
Zc = squeeze(s2z(S.Parameters, eta0));

%% --- Exact wave impedance ---
x  = 2*pi*Freq(:)/c0*a;
h2 = @(n) sqrt(pi./(2*x)).*(besselj(n+0.5, x) - 1j*bessely(n+0.5, x));
Zex = 1j*eta0*(x.*h2(l-1) - l*h2(l))./(x.*h2(l));
fprintf('TM%d: max |Z_circuit - Z_exact|/|Z_exact| = %.2e\n', l, max(abs(Zc(:) - Zex)./abs(Zex)));

figure
plot(Freq/1e9, real(Zc), Freq/1e9, imag(Zc), Freq/1e9, real(Zex), '--', Freq/1e9, imag(Zex), '--')
xlabel('f (GHz)'); ylabel('Z (\Omega)'); grid on
legend('Re Z circuit', 'Im Z circuit', 'Re Z exact', 'Im Z exact')
title(sprintf('TM_{%d} wave impedance, a = %g m', l, a))
