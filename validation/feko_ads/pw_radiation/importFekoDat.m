function s = importFekoDat(filepath)
%IMPORTFEKODAT  Import a FEKO rE far-field voltage .dat (one column per phi cut).
%   s = IMPORTFEKODAT(filepath) reads a tab-delimited FEKO export of the
%   form
%
%     Far field
%     "Frequency"[Hz]  "Phi = 60 deg [V]"  "Phi = 90 deg [V]"  "Phi = 30 deg [V]"
%     1.00000000E+08   3.48931796E-03      4.03135247E-03      2.01232712E-03
%     ...
%
%   The cut of each column is taken from its "Phi = ..." label, not its
%   position (FEKO writes the cuts in varying order). Header field i is
%   paired with data column i; the import fails if the counts disagree.
%   Columns stay in file order; downstream code looks them up by phi value
%   (findPhiColumn). s.header keeps each column's original label.
%
%   A trailing tab on each line produces an empty column, which is dropped.
%
%   Output struct (same shape as importAdsCsv, so the plotting functions
%   treat both sources identically):
%     s.source    - "FEKO"
%     s.quantity  - "mag" or "phase"
%     s.component - "Vtheta" or "Vphi" (from the filename)
%     s.unit      - "V" or "deg" (from the column header)
%     s.phi       - 1 x nPhi, phi value of each column [deg], in file order
%     s.header    - 1 x nPhi, the original column label each phi came from
%     s.freq      - nFreq x 1 frequency vector [Hz]
%     s.data      - nFreq x nPhi data matrix (column k <-> s.phi(k))
%     s.file      - the file it came from
%
%   See also IMPORTADSCSV, DOWNSAMPLEFEKO.

    arguments
        filepath (1,1) string
    end

    if ~isfile(filepath)
        error("importFekoDat:fileNotFound", "File not found: %s", filepath);
    end

    lines = readlines(filepath);

    % --- Find the header line (the one naming the phi cuts) ---------
    hdrIdx = find(contains(lines, "Phi") & contains(lines, "Frequency"), 1);
    if isempty(hdrIdx)
        hdrIdx = find(contains(lines, "Phi"), 1);
    end
    if isempty(hdrIdx)
        preview = join(lines(1:min(3,numel(lines))), " | ");
        error("importFekoDat:noHeader", ...
            ['No column header containing "Phi = ..." found in %s\n' ...
             'First lines are: %s\n' ...
             'If those look like an ADS export (a quoted "phi_prime" header ' ...
             'and comma-separated rows), the file is not a FEKO result -- ' ...
             're-export it from FEKO, or import it with importAdsCsv.'], ...
            filepath, preview);
    end

    % Header field i lines up with data column i. Only trailing empty fields
    % (from FEKO's trailing tab) are dropped; removing an inner one would
    % shift every later header onto the wrong column.
    fields = strtrim(split(lines(hdrIdx), sprintf("\t")));
    while ~isempty(fields) && fields(end) == ""
        fields(end) = [];
    end

    % --- phi value of each data column -----------------------------
    nCol = numel(fields) - 1;                % first field is Frequency
    phi  = nan(1, nCol);
    for i = 1:nCol
        tok = regexp(fields(i+1), 'Phi\s*=\s*(-?[\d.eE+-]+)', "tokens", "once");
        if ~isempty(tok)
            phi(i) = str2double(tok(1));
        end
    end
    if any(isnan(phi))
        bad = find(isnan(phi));
        error("importFekoDat:badPhiHeader", ...
            ['%s: could not read a phi value from the header of data ' ...
             'column(s) %s.\nThe header fields are:\n  %s\n' ...
             'Every data column must be labelled "Phi = <number>" -- the ' ...
             'label is what assigns a trace to a cut, not the column order.'], ...
            filepath, join(string(bad), ", "), join(fields, sprintf("\n  ")));
    end

    % --- Unit, from the first data column header, e.g. "[V]"/"[deg]"
    unit = string(regexp(fields(2), '\[([^\]]+)\]', "tokens", "once"));
    if isempty(unit) || unit == "", unit = "V"; end

    % --- Numeric block ---------------------------------------------
    M = readmatrix(filepath, "FileType","text", "Delimiter","\t", ...
        "NumHeaderLines", hdrIdx);
    M = M(~all(isnan(M), 2), :);             % drop blank trailing rows

    % Drop only trailing all-NaN columns (from the trailing tab), for the same
    % reason.
    while size(M,2) > 1 && all(isnan(M(:,end)))
        M(:,end) = [];
    end

    if size(M,2) ~= nCol + 1
        error("importFekoDat:columnMismatch", ...
            ['%s: the header labels %d data column(s) but %d were read, so ' ...
             'labels and columns cannot be matched up.\nHeader fields:\n  %s'], ...
            filepath, nCol, size(M,2)-1, join(fields, sprintf("\n  ")));
    end

    freq = M(:,1);
    data = M(:,2:end);

    header = fields(2:end).';

    % --- Quantity / component from the filename ---------------------
    [~, fname] = fileparts(filepath);
    quantity  = string(regexp(lower(fname), '(mag|phase)', "match", "once"));
    component = string(regexp(fname, 'V(theta|phi)', "match", "once"));
    if quantity == ""
        if lower(unit) == "deg"
            quantity = "phase";
        else
            quantity = "mag";
        end
    end

    s = struct( ...
        "source",    "FEKO", ...
        "quantity",  quantity, ...
        "component", component, ...
        "unit",      unit, ...
        "phi",       phi, ...
        "header",    header, ...
        "freq",      freq(:), ...
        "data",      data, ...
        "file",      string(filepath));
end