function circuit = importCircuitModelSet(dataFolder, csvFiles)
%IMPORTCIRCUITMODELSET  Import a set of ADS load-current CSV exports.
%   circuit = IMPORTCIRCUITMODELSET(dataFolder, csvFiles) returns a struct
%   array, one element per file, with fields tag, freq, mag, phase.
%   csvFiles is a cell array of file names relative to dataFolder. See
%   importCircuitModelCsv for the per-file format.

    nLoads = numel(csvFiles);
    optsCsv = detectImportOptions(fullfile(dataFolder, csvFiles{1}), "FileType","text");
    optsCsv.VariableNamingRule = "preserve";

    circuit(nLoads) = struct("tag","", "freq",[], "mag",[], "phase",[]);
    for k = 1:nLoads
        circuit(k) = importCircuitModelCsv(fullfile(dataFolder, csvFiles{k}), optsCsv);
    end
end
