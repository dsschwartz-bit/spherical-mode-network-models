function generate_mdf(data1, data2, title, varargin)
%GENERATE_MDF  Write time and one or more voltages as a Keysight ADS MDF file.
%
%   generate_mdf(time, V, title)
%   generate_mdf(time, V, title, 'BlockName', 'BLK', ...
%                'Names', {'time','Vinc_1','Vinc_2'}, 'Comments', ["a"; "b"])
%
%   File layout:
%       ! Created Tue Oct 6 23:13:10 2026
%       ! <optional comment lines>
%       BEGIN BLK
%       % time(real) Vinc_1(real) Vinc_2(real)
%       0	-1.23e-08	-1.23e-08
%       7.25e-13	3.4e-08	1.1e-06
%       ...
%       END
%
%   Inputs
%     data1 : N x 1 independent variable (time in seconds)
%     data2 : N x K dependent variables (voltages in volts), one column
%             per variable; a single row/column vector is also accepted
%     title : output file name; '.mdf' is added if it has no
%             .mdf/.mdif/.txt extension
%
%   Options
%     'BlockName' : name after BEGIN (default 'BLK')
%     'Names'     : K+1 column names, independent variable first
%                   (default {'time','voltage'} for one column,
%                   {'time','voltage1',...} for several). These are the
%                   variable names ADS sees, so they must match the
%                   expressions used in the schematic.
%     'Precision' : significant digits (default 12)
%     'Comments'  : string array / cell array of extra comment lines,
%                   each written as "! <text>" before BEGIN

p = inputParser;
addParameter(p, 'BlockName', 'BLK');
addParameter(p, 'Names',     {});
addParameter(p, 'Precision', 12);
addParameter(p, 'Comments',  {});
parse(p, varargin{:});
o = p.Results;

% ---- shapes -------------------------------------------------------------
data1 = data1(:);
N = numel(data1);
if isvector(data2)
    data2 = data2(:);
end
if size(data2, 1) ~= N && size(data2, 2) == N
    data2 = data2.';
end
if size(data2, 1) ~= N
    error('generate_mdf:size', ...
        'data2 must have %d rows (one per time point); it is %d x %d.', ...
        N, size(data2, 1), size(data2, 2));
end
K = size(data2, 2);
if N < 2
    warning('generate_mdf:short', 'Only %d point(s) to write.', N);
end

% ---- column names -------------------------------------------------------
names = cellstr(o.Names);
if isempty(o.Names)
    if K == 1
        names = {'time', 'voltage'};
    else
        names = [{'time'}, arrayfun(@(k) sprintf('voltage%d', k), 1:K, 'UniformOutput', false)];
    end
end
if numel(names) ~= K + 1
    error('generate_mdf:names', ...
        'Names needs %d entries (time + %d variable(s)); got %d.', K + 1, K, numel(names));
end
if numel(unique(names)) ~= numel(names)
    error('generate_mdf:names', 'Column names must be unique.');
end

% ---- file name ----------------------------------------------------------
title = char(title);
[~, ~, ext] = fileparts(title);
if ~any(strcmpi(ext, {'.mdf', '.mdif', '.txt'}))
    title = [title '.mdf'];
end

fid = fopen(title, 'w');
if fid < 0
    error('generate_mdf:open', 'Could not open %s for writing.', title);
end
closer = onCleanup(@() fclose(fid));

% ---- write --------------------------------------------------------------
num    = sprintf('%%.%dg', o.Precision);
rowFmt = [strjoin(repmat({num}, 1, K + 1), '\t') '\n'];

fprintf(fid, '! Created %s\n', char(datetime('now', 'Format', 'eee MMM d HH:mm:ss yyyy')));
comments = cellstr(o.Comments);
for k = 1:numel(comments)
    fprintf(fid, '! %s\n', comments{k});
end
fprintf(fid, 'BEGIN %s\n', char(o.BlockName));
fprintf(fid, '%%%s\n', sprintf(' %s(real)', names{:}));
fprintf(fid, rowFmt, [data1, data2].');
fprintf(fid, 'END\n');
end