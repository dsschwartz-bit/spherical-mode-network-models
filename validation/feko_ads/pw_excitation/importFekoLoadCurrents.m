function feko = importFekoLoadCurrents(filepath)
%IMPORTFEKOLOADCURRENTS  Import a FEKO load-current .dat export (one column per load).
%   feko = IMPORTFEKOLOADCURRENTS(filepath) returns a struct with fields:
%     freq  - frequency vector [Hz]
%     data  - matrix, one column per load (magnitude [A] or phase [deg])
%     tags  - string array, one tag per data column, parsed from the
%             "Load_<tag>" column headers (e.g. "ph30th60")
%
%   Expected file format:
%     Line 1  : "Load"
%     Line 2  : tab-separated quoted header, e.g.
%               "Frequency"[Hz]  "Load_ph30th60 [A]"  "Load_ph60th60 [A]" ...
%     Lines 3+: tab-separated numeric data (the trailing tab produces an
%               empty column, which is dropped)

    rawLines = readlines(filepath);
    headerFields = split(rawLines(2), sprintf("\t"));
    headerFields = strtrim(headerFields);
    headerFields(headerFields == "") = [];   % drop trailing empty field

    % First header field is Frequency; the rest are the load columns
    tags = strings(numel(headerFields)-1, 1);
    for i = 2:numel(headerFields)
        tags(i-1) = regexp(headerFields(i), "ph\d+th\d+", "match", "once");
    end

    M = readmatrix(filepath, "FileType","text", "Delimiter","\t", ...
        "NumHeaderLines", 2);
    M = M(:, all(~isnan(M),1));   % drop any all-NaN trailing column

    feko.freq = M(:,1);
    feko.data = M(:,2:end);
    feko.tags = tags;
end
