function [Frequency_Hz, Zin_Ohm] = ReadFEKOimpedance(filename, dataLines)
%READFEKOIMPEDANCE  Read an impedance-vs-frequency export (.dat) from FEKO.
%  [Frequency_Hz, Zin_Ohm] = READFEKOIMPEDANCE(FILENAME) reads a
%  tab-delimited file with columns Frequency (Hz), Re(Z) (Ohm), Im(Z) (Ohm),
%  starting at line 3. Returns column vectors sorted by ascending frequency,
%  with Zin_Ohm complex.
%
%  [...] = READFEKOIMPEDANCE(FILENAME, DATALINES) reads the row interval(s)
%  DATALINES (positive scalar or N-by-2 array of [first last] rows).
%
%  Example:
%  [Frequency_Hz, Zin_Ohm] = ReadFEKOimpedance("antenna_impedance.afs.dat", [3, Inf]);
%
%  See also READTABLE, ReadFEKOreflection.

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [3, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Frequency_Hz", "Real_Ohm", "Imaginary_Ohm"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
common_mode_gndplane_noloads_afs = readtable(filename, opts);

%% Convert to output type
Frequency_Hz = common_mode_gndplane_noloads_afs.Frequency_Hz;
Zin_Ohm = common_mode_gndplane_noloads_afs.Real_Ohm + 1j*common_mode_gndplane_noloads_afs.Imaginary_Ohm;

%% Sort by frequency ascending
[Frequency_Hz, indices] = sort(Frequency_Hz, 'ascend');
Zin_Ohm = Zin_Ohm(indices);
end