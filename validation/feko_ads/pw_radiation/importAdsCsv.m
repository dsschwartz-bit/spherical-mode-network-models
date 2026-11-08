function s = importAdsCsv(filepath)
%IMPORTADSCSV  Import an ADS rE far-field voltage CSV (long / stacked format).
%   s = IMPORTADSCSV(filepath) reads an ADS CSV of the form
%
%     "phi_prime                 ","freq                      ","mag(Vphi)  ...
%     References : ['mag(Vphi)']  References : ['mag(Vphi)']    References : []
%     Dependency : [phi_prime]    Dependency : [phi_prime,freq] ...
%     Num. Points : [3]           Num. Points : [3,97]          ...
%     Matrix Size : scalar        ...
%     Type : Real                 ..."
%     30.000,1.000E8,0.004
%     30.000,1.250E8,0.005
%     ...
%
%   i.e. a quoted multi-line header block followed by numeric rows in
%   LONG format: one row per (phi, frequency) pair. The sweep variable
%   "phi_prime" is treated as phi (the "prime" is ignored, per the
%   convention used for the rotated-loop geometry).
%
%   The data are reshaped into a wide matrix, one column per phi value.
%
%   Output struct:
%     s.source    - "ADS"
%     s.quantity  - "mag" or "phase" (parsed from the column name)
%     s.component - "Vtheta" or "Vphi"
%     s.unit      - "V" for magnitude, "deg" for phase
%     s.phi       - 1 x nPhi, phi values, ascending
%     s.header    - 1 x nPhi, label of each column
%     s.freq      - nFreq x 1 frequency vector [Hz]
%     s.data      - nFreq x nPhi data matrix (column k <-> s.phi(k))
%     s.file      - the file it came from
%
%   See also IMPORTFEKODAT.

    arguments
        filepath (1,1) string
    end

    if ~isfile(filepath)
        error("importAdsCsv:fileNotFound", "File not found: %s", filepath);
    end

    raw = fileread(filepath);

    % --- Locate the numeric rows ------------------------------------
    % Three comma-separated numbers on their own line. The header block
    % contains quotes/letters, so it can never match this pattern -- which
    % makes the parser insensitive to how many header lines ADS writes.
    numPat = ['^\s*(-?[\d.]+(?:[eE][+-]?\d+)?)\s*,' ...
              '\s*(-?[\d.]+(?:[eE][+-]?\d+)?)\s*,' ...
              '\s*(-?[\d.]+(?:[eE][+-]?\d+)?)\s*$'];

    tok = regexp(raw, numPat, "tokens", "lineanchors");
    if isempty(tok)
        error("importAdsCsv:noData", ...
            "No 3-column numeric data found in %s", filepath);
    end
    V = str2double(vertcat(tok{:}));          % nRows x 3

    % --- Parse the header block to identify the columns -------------
    hdrEnd = regexp(raw, numPat, "start", "once", "lineanchors");
    hdrTxt = raw(1:max(hdrEnd-1,1));

    % Each quoted block starts with the variable name on its own line.
    nameTok = regexp(hdrTxt, '"([^"\r\n]+)[\r\n]', "tokens");
    names   = strings(1, numel(nameTok));
    for i = 1:numel(nameTok)
        names(i) = strtrim(string(nameTok{i}{1}));
    end

    iPhi = []; iFreq = []; iVal = [];
    if numel(names) == 3
        lname = lower(names);
        iPhi  = find(startsWith(lname, "phi"), 1);
        iFreq = find(contains(lname, "freq"), 1);
        iVal  = find(contains(lname, ["mag(", "phase(", "db(", "real(", "imag("]), 1);
    end
    if isempty(iPhi) || isempty(iFreq) || isempty(iVal)
        % Fallback: ADS writes the columns in sweep order (phi, freq, value)
        warning("importAdsCsv:headerFallback", ...
            "Could not parse the header of %s; assuming column order [phi, freq, value].", ...
            filepath);
        iPhi = 1; iFreq = 2; iVal = 3;
        names = ["phi_prime", "freq", ""];
    end

    % --- Quantity / component --------------------------------------
    valName = names(iVal);
    quantity  = string(regexp(lower(valName), '^(mag|phase|db)', "match", "once"));
    component = string(regexp(valName, 'V(theta|phi)', "match", "once"));

    [~, fname] = fileparts(filepath);   % fall back to the filename
    if quantity == ""
        quantity = string(regexp(lower(fname), '(mag|phase)', "match", "once"));
    end
    if component == ""
        component = string(regexp(fname, 'V(theta|phi)', "match", "once"));
    end

    % --- Reshape long -> wide (one column per phi) ------------------
    phiVals = unique(V(:,iPhi));          % ascending
    nPhi    = numel(phiVals);

    rows0 = V(:,iPhi) == phiVals(1);
    freq0 = sort(V(rows0, iFreq));
    nFreq = numel(freq0);

    data = nan(nFreq, nPhi);
    for k = 1:nPhi
        rows = V(:,iPhi) == phiVals(k);
        [fk, ord] = sort(V(rows, iFreq));
        dk = V(rows, iVal);
        if numel(fk) ~= nFreq
            error("importAdsCsv:raggedSweep", ...
                "phi = %g has %d points but phi = %g has %d in %s.", ...
                phiVals(k), numel(fk), phiVals(1), nFreq, filepath);
        end
        if max(abs(fk - freq0)) > 1e-6 * max(freq0)
            error("importAdsCsv:freqMismatch", ...
                "Frequency grids differ between phi sweeps in %s.", filepath);
        end
        data(:,k) = dk(ord);
    end

    if quantity == "phase"
        unit = "deg";
    else
        unit = "V";
    end

    % Label of each column, so the phi <-> column mapping can be inspected
    % the same way as for FEKO (see importFekoDat).
    header = names(iPhi) + " = " + string(phiVals(:).');

    s = struct( ...
        "source",    "ADS", ...
        "quantity",  quantity, ...
        "component", component, ...
        "unit",      unit, ...
        "phi",       phiVals(:).', ...
        "header",    header, ...
        "freq",      freq0(:), ...
        "data",      data, ...
        "file",      string(filepath));
end