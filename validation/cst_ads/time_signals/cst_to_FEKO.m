%CST_TO_FEKO  Convert a CST excitation time signal to a FEKO time-signal CSV.
% Reads a CST ASCII time-signal export and writes it as a two-column CSV
% (time in ns, voltage in V) with cst_to_csv, then plots the signal and its
% spectrum (plotCSTSignal).
% Input : <filepath><filename>.txt, a CST time-signal export.
% Output: output/cst_ads/csv_out/<filename>_4FEKO.csv (cst_to_csv default folder).
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(fullfile(here, 'utils')); % plotCSTSignal

%% ===== USER SETTINGS =====
filepath = "C:\path\to\CST\project\Export\time_signals\";   % folder of the CST time-signal export
% =========================
filename="excitation_sinestep";
cst_to_csv(fullfile(filepath, filename+".txt"), filename+"_4FEKO.csv", "Precision",8)
plotCSTSignal(fullfile(filepath, filename+".txt"))
