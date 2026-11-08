function fp = buildDataFilePath(scriptFolder, filename)
%BUILDDATAFILEPATH  Full path of a data file stored by antenna type.
%   fp = BUILDDATAFILEPATH(scriptFolder, filename) returns
%     <scriptFolder>/<antennaType>/<antennaSubtype>/<filename>
%   where antennaType and antennaSubtype are the first two underscore-separated
%   tokens of filename, e.g.
%     "loop_bowtie_ph30th60r5.s4p" -> <scriptFolder>/loop/bowtie/loop_bowtie_ph30th60r5.s4p

    parts = split(filename, "_");
    antennaType    = parts(1);
    antennaSubtype = parts(2);
    fp = fullfile(scriptFolder, antennaType, antennaSubtype, filename);
end
