function data = read_ads_farfield_csv(filename)
%READ_ADS_FARFIELD_CSV  Read an ADS dataset CSV export into column vectors.
%   data = READ_ADS_FARFIELD_CSV(filename) returns a struct with one (n x 1)
%   field per CSV column, raw values in file row order (no reshaping, no
%   mag/phase combination; see compare_cst_ads_farfields.m for that).
%
% ADS puts a multi-line metadata block inside each quoted header cell, so
% the file is split with parse_csv_rows (RFC 4180 quoting). The field name
% is the first line of the header cell with non-alphanumeric characters
% removed ("mag(Vtheta)" -> magVtheta; a leading digit gets an 'f' prefix).
% Columns whose names collide raise an error. Fields ending in "prime"
% (e.g. thetaprime) are also stored without that suffix unless that name
% already exists.
%
% data.column_names / data.column_names_raw list the field names and the
% original first-line headers, in column order.

rows = parse_csv_rows(filename);
if numel(rows) < 2
    error('read_ads_farfield_csv:empty', 'File "%s" has no data rows.', filename);
end

header = rows{1};
ncols = numel(header);

column_names_raw = cell(1, ncols);
column_names     = cell(1, ncols);
for c = 1:ncols
    first_line = strtok(header{c}, sprintf('\r\n'));
    first_line = strtrim(first_line);
    column_names_raw{c} = first_line;
    sanitized = first_line(isstrprop(first_line, 'alphanum'));
    if isempty(sanitized)
        error('read_ads_farfield_csv:badHeader', ...
            'Column %d''s header ("%s") has no letters or digits to build a variable name from, in "%s".', ...
            c, first_line, filename);
    end
    if isstrprop(sanitized(1), 'digit')
        sanitized = ['f' sanitized]; %#ok<AGROW>  % a leading digit isn't a valid MATLAB identifier
    end
    column_names{c} = sanitized;
end
dup = column_names(cellfun(@(n) sum(strcmp(column_names, n)), column_names) > 1);
if ~isempty(dup)
    error('read_ads_farfield_csv:duplicateFieldName', ...
        'Columns [%s] in "%s" sanitize to duplicate field name(s) [%s] -- cannot tell them apart.', ...
        strjoin(column_names_raw, ', '), filename, strjoin(unique(dup), ', '));
end

n = numel(rows) - 1;
cols = zeros(n, ncols);
for r = 1:n
    data_row = rows{r+1};
    if numel(data_row) ~= ncols
        error('read_ads_farfield_csv:rowLength', 'Data row %d of "%s" has %d fields, expected %d.', r, filename, numel(data_row), ncols);
    end
    for c = 1:ncols
        v = str2double(data_row{c});
        if isnan(v)
            error('read_ads_farfield_csv:parseError', 'Could not parse "%s" as a number (row %d, column %d of "%s").', data_row{c}, r, c, filename);
        end
        cols(r,c) = v;
    end
end

data = struct();
for c = 1:ncols
    data.(column_names{c}) = cols(:,c);
end

% Alias "...prime" fields (e.g. theta_prime -> thetaprime) to the plain name
% (theta), without overwriting an existing column.
for c = 1:ncols
    nm = column_names{c};
    if numel(nm) > 5 && endsWith(nm, 'prime')
        stem = nm(1:end-5);
        if ~isfield(data, stem)
            data.(stem) = cols(:,c);
        end
    end
end

data.column_names     = column_names;
data.column_names_raw = column_names_raw;

end
