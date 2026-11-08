function farFieldData = importffarfield_HFSS(baseName, folder)
%IMPORTFFARFIELD_HFSS  Load HFSS far-field CSV exports.
%   farFieldData = IMPORTFFARFIELD_HFSS(baseName, folder) scans folder for
%   csv files named <baseName>_<port>_<mode>.csv (one per excitation), or a
%   single <baseName>.csv (no split by source), as written by the HFSS
%   export script python/hfss/ExportFarFields.py (run inside HFSS), and
%   returns a struct array with one entry per file.
%
%   baseName - the output_file_name base used in ExportFarFields.py.
%              By default that is "<project>_<design>_farfield_rE", e.g.
%              "rectangular_patch_load_port_farfield_rE" for a project
%              named "rectangular_patch" and a design named "load_port".
%   folder   - directory containing the csv files. Defaults to the
%              current directory if omitted.
%
%   Each element of the returned struct array has fields:
%       port     - port number (double), or [] for the no-split case
%       mode     - mode number (double), or [] for the no-split case
%       label    - "<port>_<mode>" string as it appeared in the
%                  filename, or "combined" for the no-split case
%       filename - full path to the source csv file
%       headers  - cell array of column header strings, exactly as
%                  written in the csv (e.g. "mag(rEPhi) (V)") - read as
%                  raw text, never passed through MATLAB's table/struct
%                  variable-name validation, so nothing gets sanitized,
%                  no parentheses stripped, no spaces collapsed
%       data     - numeric matrix of the csv body (rows = data points,
%                  columns match headers 1:1)
%
%   Since column names are never turned into MATLAB identifiers, get a
%   column with the helper function below rather than table/struct field
%   access:
%       col = farFieldColumn(entry, 'mag(rEPhi) (V)')
%
%   Example:
%       ff = importffarfield_HFSS('rectangular_patch_load_port_farfield_rE');
%       entry = ff([ff.port]==1 & [ff.mode]==1);
%       theta = farFieldColumn(entry, 'Theta(deg)');
%       magPhi = farFieldColumn(entry, 'mag(rEPhi) (V)');
%       plot(theta, magPhi);

    if nargin < 2 || isempty(folder)
        folder = pwd;
    end

    % Use char arrays: [a b] with a string operand would build a string
    % array instead of concatenating characters.
    baseName = char(baseName);
    folder = char(folder);

    allCsvFiles = dir(fullfile(folder, [baseName '*.csv']));
    if isempty(allCsvFiles)
        error('importffarfield_HFSS:noFilesFound', ...
            'No csv files found matching "%s" in %s', ...
            [baseName '*.csv'], folder);
    end

    farFieldData = struct('port', {}, 'mode', {}, 'label', {}, ...
                           'filename', {}, 'headers', {}, 'data', {});

    for k = 1:numel(allCsvFiles)
        fname = allCsvFiles(k).name;
        fullPath = fullfile(folder, fname);

        [~, nameNoExt, ~] = fileparts(fname);
        suffix = nameNoExt(length(baseName)+1:end);
        suffix = regexprep(suffix, '^_+', ''); %strip leading underscore

        entry.filename = fullPath;

        %Read the header line as raw text - never handed to MATLAB's
        %variable-name validator, so it can't be sanitized in any way.
        fid = fopen(fullPath, 'r');
        headerLine = fgetl(fid);
        fclose(fid);
        entry.headers = strtrim(strsplit(headerLine, ','));

        %Read the numeric body, skipping the header row.
        entry.data = readmatrix(fullPath, 'NumHeaderLines', 1);

        if k == 1
            fprintf('Detected columns in %s:\n', fname);
            fprintf('  %s\n', entry.headers{:});
        end

        %suffix is either empty ("combined", no source split) or
        %"<port>_<mode>", e.g. "1_1"
        tokens = regexp(suffix, '^(\d+)_(\d+)$', 'tokens', 'once');
        if isempty(suffix)
            entry.port = [];
            entry.mode = [];
            entry.label = 'combined';
        elseif ~isempty(tokens)
            entry.port = str2double(tokens{1});
            entry.mode = str2double(tokens{2});
            entry.label = suffix;
        else
            %doesn't match the expected <port>_<mode> pattern - keep the
            %raw suffix as the label rather than silently dropping the file
            entry.port = [];
            entry.mode = [];
            entry.label = suffix;
            warning('importffarfield_HFSS:unrecognizedSuffix', ...
                ['File "%s" has a suffix ("%s") that does not match ' ...
                 'the expected <port>_<mode> pattern - loaded anyway ' ...
                 'with an unparsed label.'], fname, suffix);
        end

        farFieldData(end+1) = entry; %#ok<AGROW>
    end
end