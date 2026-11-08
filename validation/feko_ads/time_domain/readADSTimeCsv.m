function [time, cols] = readADSTimeCsv(csvFile, expectedRefs)
%READADSTIMECSV  Read an ADS data-display CSV export (multi-line quoted header, then numeric rows).
%   [time, cols] = READADSTIMECSV(csvFile, expectedRefs)
%   csvFile      : path to the .csv
%   expectedRefs : string array of the expected trace names in file order (columns 2..end),
%                  as in the "References : [...]" line of the header.
%   time is column 1, cols the remaining columns. Stops with an error if the references in the
%   file differ from expectedRefs, so a re-export in a different order is not mislabelled.
L = readlines(csvFile);
ref = regexp(L(2), "References\s*:\s*\[(.*)\]", "tokens", "once");
if isempty(ref), error("readADSTimeCsv: no 'References' line in %s", csvFile); end
fileRefs = string(regexp(ref{1}, "'([^']*)'", "tokens"));
fileRefs = reshape(fileRefs, 1, []);
if ~isequal(fileRefs, reshape(string(expectedRefs), 1, []))
    error("readADSTimeCsv: column references in %s are\n  %s\nexpected\n  %s\nUpdate the import mapping.", ...
        csvFile, strjoin(fileRefs, " | "), strjoin(expectedRefs, " | "));
end
% Parse the numeric rows directly: readmatrix/readtable count the quoted multi-line header as
% fewer lines than it occupies, and NumHeaderLines then skips real data rows.
isData = ~cellfun(@isempty, regexp(cellstr(L), '^\s*[-+]?\d', 'once'));
nCol = numel(expectedRefs) + 1;
M = str2double(split(L(isData), ",", 2));    % one row per data line, one column per field
M = M(all(isfinite(M), 2), :);
if size(M, 2) ~= nCol || isempty(M)
    error("readADSTimeCsv: %d columns in %s, expected %d", size(M, 2), csvFile, nCol);
end
time = M(:, 1);  cols = M(:, 2:end);
end
