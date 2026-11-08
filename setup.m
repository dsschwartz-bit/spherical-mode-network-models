function setup()
%SETUP  Add the toolbox's MATLAB folders to the MATLAB path.
%   Run once per MATLAB session from the repository root, or let the
%   workflow scripts call it for you. Examples and tests add their own
%   folders when run, so they are not put on the path here.

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'matlab')));
addpath(genpath(fullfile(root, 'third_party', 'matlab')));
end
