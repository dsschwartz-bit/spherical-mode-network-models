function rows = parse_csv_rows(filename)
%PARSE_CSV_ROWS  Minimal RFC 4180 CSV reader.
%   rows = PARSE_CSV_ROWS(filename) splits a file into rows and fields,
%   keeping quoted fields with embedded commas or newlines (CR, LF, CRLF)
%   intact and reading "" inside quotes as a literal quote. Needed for ADS
%   CSV exports, whose quoted header cells span several lines.
%
% Returns a cell array of cell arrays: rows{i}{j} is row i, field j, as
% a raw (unquoted, "" already collapsed to ") char string. Rows may have
% differing numbers of fields if the file does.

txt = fileread(filename);
n = numel(txt);

rows = {};
row = {};
field = '';
in_quotes = false;
i = 1;
while i <= n
    c = txt(i);
    if in_quotes
        if c == '"'
            if i < n && txt(i+1) == '"'
                field = [field '"']; %#ok<AGROW>
                i = i + 2;
                continue
            else
                in_quotes = false;
                i = i + 1;
                continue
            end
        else
            field = [field c]; %#ok<AGROW>
            i = i + 1;
            continue
        end
    else
        if c == '"'
            in_quotes = true;
            i = i + 1;
            continue
        elseif c == ','
            row{end+1} = field; %#ok<AGROW>
            field = '';
            i = i + 1;
            continue
        elseif c == char(13) || c == char(10)
            row{end+1} = field; %#ok<AGROW>
            field = '';
            rows{end+1} = row; %#ok<AGROW>
            row = {};
            if c == char(13) && i < n && txt(i+1) == char(10)
                i = i + 2;
            else
                i = i + 1;
            end
            continue
        else
            field = [field c]; %#ok<AGROW>
            i = i + 1;
            continue
        end
    end
end
% flush a final field/row if the file doesn't end with a newline
if ~isempty(field) || ~isempty(row)
    row{end+1} = field; %#ok<AGROW>
    rows{end+1} = row; %#ok<AGROW>
end

% drop any trailing fully-empty row (e.g. from a final blank line)
if ~isempty(rows) && numel(rows{end}) == 1 && isempty(rows{end}{1})
    rows(end) = [];
end

end
