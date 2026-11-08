%% generate_impedance_csv  Write lumped-load impedance tables for FEKO/ADS.
% Writes CSV files with columns [frequency (Hz), Re(Z) (Ohm), Im(Z) (Ohm)] for
% f = 10 MHz to 100.01 GHz in 100 MHz steps:
%   Zcap<C>pF.csv        - ideal capacitors, Z = 1/(j*w*C)
%   Zsw_on_<R>Ohm.csv,
%   Zsw_off_<R>kOhm.csv  - switch model: R_on or R_off in parallel with C_iso
% Output folder: <repo>/output/impedance (created if missing).

out_dir = fullfile(repo_root, 'output', 'impedance');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

f=(1:10:10001).'*1e7;
cap=[0.1,1.1,1.75,2.4]*1e-12; %Farad
Zcap=1./(1j*2*pi*f.*cap);
plot(f,db(imag(Zcap)))

for n=1:length(cap)
writematrix([f,real(Zcap(:,n)),imag(Zcap(:,n))],fullfile(out_dir,"Zcap"+num2str(cap(n)*1e12)+"pF"+".csv"))
end

Ron=0;Roff=25000; %Ohms
Ciso = 0.2e-12; %Farad
Zsw_on = 1./ ( 1/Ron + (1j*2*pi*f*Ciso) ) ;
Zsw_off = 1./ ( 1/Roff + (1j*2*pi*f*Ciso) ) ;
writematrix([f,real(Zsw_on),imag(Zsw_on)],fullfile(out_dir,"Zsw_on_"+num2str(Ron)+"Ohm"+".csv"))
writematrix([f,real(Zsw_off),imag(Zsw_off)],fullfile(out_dir,"Zsw_off_"+num2str(Roff/1e3)+"kOhm"+".csv"))
