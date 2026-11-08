%% import_ADS_csv.m
% Imports the ADS export of the dipole load current under plane-wave excitation,
% Z-aligned and rotated by 45 deg (PW_exciation_dipole_comparison.csv).
%
% Input : adsCsvFile - path of the CSV, set by the calling script
%         (default: PW_exciation_dipole_comparison.csv next to this script).
%         Columns: freq [Hz], mag_rotated [A], phase_rotated [deg], mag [A], phase [deg].
% Output (workspace variables): freq, mag_rotated, phase_rotated, mag, phase.

if ~exist('adsCsvFile', 'var')
    adsCsvFile = fullfile(fileparts(mfilename('fullpath')), "PW_exciation_dipole_comparison.csv");
end

%% Set up the import options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

opts.DataLines = [2, Inf];
opts.Delimiter = ",";

opts.VariableNames = ["freq", "mag_rotated", "phase_rotated", "mag", "phase"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

PW_exciation_dipole_comparison = readtable(adsCsvFile, opts);

%% Convert to output type
freq = PW_exciation_dipole_comparison.freq;
mag_rotated = PW_exciation_dipole_comparison.mag_rotated;
phase_rotated = PW_exciation_dipole_comparison.phase_rotated;
mag = PW_exciation_dipole_comparison.mag;
phase = PW_exciation_dipole_comparison.phase;

%% Clear temporary variables
clear PW_exciation_dipole_comparison opts