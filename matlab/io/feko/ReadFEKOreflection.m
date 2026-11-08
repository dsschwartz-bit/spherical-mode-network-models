function [Frequency_Hz, Gamma] = ReadFEKOreflection(filename, dataLines)
%READFEKOREFLECTION  Read a reflection-coefficient export (.dat) from FEKO.
%  [Frequency_Hz, Gamma] = READFEKOREFLECTION(FILENAME) reads a
%  tab-delimited file with columns Frequency (Hz), |Gamma|, phase (degrees),
%  starting at line 3. Returns column vectors sorted by ascending frequency,
%  with Gamma = |Gamma|*exp(j*phase) complex.
%
%  [...] = READFEKOREFLECTION(FILENAME, DATALINES) reads the row interval(s)
%  DATALINES (positive scalar or N-by-2 array of [first last] rows).
%
%  Example:
%  [Frequency_Hz, Gamma] = ReadFEKOreflection("antenna_S11.afs.dat", [3, Inf]);
%
%  See also READTABLE, ReadFEKOimpedance.

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
opts.VariableNames = ["Frequency_Hz", "magnitude", "phase_deg"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
common_mode_gndplane_noloads_afs = readtable(filename, opts);

%% Convert to output type
Frequency_Hz = common_mode_gndplane_noloads_afs.Frequency_Hz;
Gamma = common_mode_gndplane_noloads_afs.magnitude .* exp(1j*(pi/180)*common_mode_gndplane_noloads_afs.phase_deg);

%% Sort by frequency ascending
[Frequency_Hz, indices] = sort(Frequency_Hz, 'ascend');
Gamma = Gamma(indices);
end