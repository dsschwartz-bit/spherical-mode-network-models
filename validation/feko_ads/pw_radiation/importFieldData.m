function s = importFieldData(filepath, expectedSource, phiScale)
%IMPORTFIELDDATA  Import an rE far-field voltage file, dispatching on CONTENT.
%   s = IMPORTFIELDDATA(filepath, expectedSource, phiScale) looks at what is
%   actually inside the file rather than trusting its name or extension, and
%   calls importAdsCsv or importFekoDat accordingly.
%
%   expectedSource ("ADS" or "FEKO", optional) is what the filename claims.
%   If the content says otherwise the file is still imported -- using the
%   parser that fits it -- but a warning is raised, because a mislabeled
%   export usually means the wrong file was written to that name, and an
%   "ADS vs FEKO" plot fed two ADS exports is not a comparison.
%
%   phiScale (optional, default 1) multiplies the phi values after import.
%   Use it when the sweep variable was exported in the wrong angular unit,
%   e.g. phiScale = pi/180 when phi comes back as 1719 / 3438 / 5157
%   instead of 30 / 60 / 90 (radians reported as degrees).
%
%   Recognised formats:
%     ADS  - quoted multi-line header, comma-separated "phi, freq, value"
%     FEKO - "Far field" line, tab-delimited, columns headed "Phi = NN deg"
%
%   See also IMPORTADSCSV, IMPORTFEKODAT, CHECKPHICOVERAGE.

    arguments
        filepath       (1,1) string
        expectedSource (1,1) string = ""
        phiScale       (1,1) double = 1
    end

    fmt = detectFormat(filepath);

    switch fmt
        case "FEKO", s = importFekoDat(filepath);
        case "ADS",  s = importAdsCsv(filepath);
        otherwise
            error("importFieldData:unknownFormat", ...
                ['Cannot tell what format %s is.\n' ...
                 'Expected either a FEKO export ("Far field", tab-delimited, ' ...
                 'columns headed "Phi = NN deg") or an ADS export (quoted ' ...
                 'multi-line header, comma-separated phi/freq/value rows).'], ...
                filepath);
    end

    if expectedSource ~= "" && fmt ~= upper(expectedSource)
        [~, nm, ext] = fileparts(filepath);
        warning("importFieldData:sourceMismatch", ...
            ['%s%s is named as a %s export but contains %s-format data.\n' ...
             'It has been imported as %s. Re-export it from %s before ' ...
             'trusting any comparison made with it.'], ...
            nm, ext, upper(expectedSource), fmt, fmt, upper(expectedSource));
    end

    if phiScale ~= 1
        s.phi = s.phi * phiScale;
    end
end

% ----------------------------------------------------------------------
function fmt = detectFormat(filepath)
%DETECTFORMAT  Sniff the first few lines of a file.

    if ~isfile(filepath)
        error("importFieldData:fileNotFound", "File not found: %s", filepath);
    end

    fid = fopen(filepath, "r");
    head = string(fread(fid, 4096, "*char").');
    fclose(fid);

    tab = sprintf("\t");

    if contains(head, "Phi") && contains(head, "=") && contains(head, tab)
        fmt = "FEKO";
    elseif contains(head, "Far field") && contains(head, tab)
        fmt = "FEKO";
    elseif contains(head, "Dependency :") || startsWith(strtrim(head), """")
        fmt = "ADS";
    else
        fmt = "";
    end
end