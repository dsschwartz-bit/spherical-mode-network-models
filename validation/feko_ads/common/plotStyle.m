function S = plotStyle(varargin)
%PLOTSTYLE  Line, marker and colour settings shared by the validation/feko_ads plot scripts.
%   S = PLOTSTYLE() returns the defaults; S = PLOTSTYLE('name', value, ...) overrides fields.
%   ADS = solid line, FEKO = filled markers only (downsampled), same colour for the same case.
%   Also sets LaTeX as the default text interpreter (the axis labels and legends use LaTeX).

S.fontSize    = 12;            % axes / tick labels
S.labelSize   = 13;            % axis labels, corner labels
S.lineWidth   = 1.5;           % pt, ADS (and other) lines
S.markerSize  = 6;             % FEKO markers
S.markerEdge  = 0.6;           % pt, marker edge (darker shade of the face colour)
S.cornerLabel = true;          % e.g. (theta, phi) in the top-right corner of a panel
% Wong (2011, Nature Methods) colour order, one colour per case (angle, load, ...)
S.colors = [0      0      0;       % black
            0.902  0.624  0;       % orange
            0.337  0.706  0.914;   % sky blue
            0      0.620  0.451;   % bluish green
            0.941  0.894  0.259;   % yellow
            0      0.447  0.698;   % blue
            0.835  0.369  0;       % vermillion
            0.800  0.475  0.655];  % reddish purple
% single-case comparisons (time domain): FEKO black markers, ADS red line
S.colFEKO = [0 0 0];
S.colADS  = [0.85 0 0];
S.colAlt  = [0 0.447 0.698];   % a second ADS variant / excitation
% excitation / generator signal: thick, semi-transparent light blue
S.colExc   = [0.2 0.5 0.88];
S.lwExc    = 3;                % pt
S.alphaExc = 0.6;              % opacity

for k = 1:2:numel(varargin), S.(varargin{k}) = varargin{k + 1}; end
set(groot, 'defaultTextInterpreter', 'latex', 'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', 'defaultPolaraxesTickLabelInterpreter', 'latex');
end
