%% check_import_integrity.m
% Import check for the rotated-loop rE voltage exports read by plot_rE_voltage_comparison.m
% (loopCase = "rotated"). For each ADS .csv / FEKO .dat it prints the first raw numeric line of
% the file next to the imported values, the column -> phi mapping and the phi values in degrees.
% Output: console only.

clear; clc;

scriptFolder = fileparts(mfilename("fullpath"));
if isempty(scriptFolder), scriptFolder = pwd; end
run(fullfile(scriptFolder, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(scriptFolder, fullfile(scriptFolder, '..', 'common'));

%% ===== USER SETTINGS =====
dataFolder = scriptFolder;   % folder containing the looprotated_noloads_* exports
% =========================
components = ["Vtheta", "Vphi"];
adsFile  = @(q,c) "looprotated_noloads_ADS_rE_"  + q + "_" + c + ".csv";
fekoFile = @(q,c) "looprotated_noloads_FEKO_rE_" + q + "_" + c + ".dat";

adsPhiUnit  = "rad";
fekoPhiUnit = "deg";

for comp = components
    for q = ["mag", "phase"]
        report(fullfile(dataFolder, adsFile(q, comp)),  "ADS",  adsPhiUnit);
        report(fullfile(dataFolder, fekoFile(q, comp)), "FEKO", fekoPhiUnit);
    end
end

% ----------------------------------------------------------------------
function report(filepath, expectedSource, phiUnit)

    [~, nm, ext] = fileparts(filepath);
    fprintf("\n=== %s%s (named %s) ===\n", nm, ext, expectedSource);

    if ~isfile(filepath)
        fprintf("  MISSING\n");
        return;
    end

    % ---- Raw: first line of the file that is numeric data --------------
    lines = readlines(filepath);
    isNum = ~cellfun(@isempty, regexp(lines, '^\s*[-\d]', "once"));
    firstRaw = lines(find(isNum, 1));
    fprintf("  raw first data line : %s\n", strtrim(firstRaw));

    % ---- Imported ------------------------------------------------------
    s = importFieldData(filepath, expectedSource);
    fprintf("  parsed as           : %s %s %s\n", s.source, s.quantity, s.component);
    fprintf("  phi as exported     : [%s]\n", join(string(s.phi), ", "));
    if isfield(s, "header")
        fprintf("  column -> phi mapping (file order):\n");
        for i = 1:numel(s.phi)
            fprintf("    col %d : %-28s -> phi = %g\n", i, s.header(i), s.phi(i));
        end
    end

    sDeg = normalizePhi(s, [], 1, phiUnit);   % unit conversion only, no auto
    fprintf("  phi in degrees      : [%s]\n", join(string(round(sDeg.phi,3)), ", "));

    fprintf("  first row  (f = %.4g Hz) : [%s]\n", s.freq(1), ...
        join(string(s.data(1,:)), ", "));
    fprintf("  last row   (f = %.4g Hz) : [%s]\n", s.freq(end), ...
        join(string(s.data(end,:)), ", "));
    fprintf("  data min / max      : %g / %g   (columns are sorted by phi)\n", ...
        min(s.data(:)), max(s.data(:)));
end