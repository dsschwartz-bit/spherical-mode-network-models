function s = importCircuitModelCsv(filepath, optsCsv)
%IMPORTCIRCUITMODELCSV  Import one ADS load-current CSV export.
%   s = IMPORTCIRCUITMODELCSV(filepath, optsCsv) reads a CSV with columns
%   freq [Hz], mag(I_Probe1.i) [A], phase(I_Probe1.i) [deg] below a quoted,
%   multi-line ADS header. optsCsv is a detectImportOptions object, reused
%   across files of the same layout.
%
%   The incidence tag (e.g. "ph30th60") is parsed from the file name and
%   returned in s.tag so the curve can be matched to the FEKO column.
%
%   Output struct fields: tag, freq, mag, phase.

    T = readtable(filepath, optsCsv);

    [~, name, ~] = fileparts(filepath);
    tg = regexp(name, "ph\d+th\d+", "match", "once");

    s.tag   = tg;
    s.freq  = T{:,1};
    s.mag   = T{:,2};
    s.phase = T{:,3};
end
