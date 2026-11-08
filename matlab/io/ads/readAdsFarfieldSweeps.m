function ads = readAdsFarfieldSweeps(file, components, ads_phi)
%READADSFARFIELDSWEEPS  Far-field time signals from an ADS Data Display CSV.
%
%   ads = readAdsFarfieldSweeps(file, components, ads_phi)
%
%   file       : CSV export with a time column, V<comp>_prime columns (r*E
%                in V, e.g. Vtheta_prime, Vphi_prime) and the sweep
%                angles: theta_deg / phi_deg in degrees, or theta / phi in
%                radians. Without a theta column, theta is unknown (NaN)
%                and the file holds one sweep per phi (or, without a phi
%                column either, one direction at ads_phi).
%   components : string array, e.g. ["Theta" "Phi"]
%   ads_phi    : phi (deg) of the ADS sweep when the file has no phi
%                column (the far field was computed at this phi only)
%
%   Sweeps outside theta = 0..180 (e.g. theta_deg = -30, i.e. the far side
%   of the phi = ads_phi plane) are mapped to (|theta|, phi + 180, wrapped
%   to 0..360), with E_theta and E_phi sign-flipped, the same convention as
%   readCstFarfieldProbes (theta-hat and phi-hat reverse across the z axis).
%
%   ads : struct array, one entry per direction, with fields theta, phi
%         (deg, NaN theta = no theta column), phiInFile (true if phi came
%         from a file column), t (s, repeated transient breakpoints
%         removed) and e (N x numel(components), NaN = column missing).

[names, A] = readAdsCsv(file);
iT = find(strcmpi(names, "time"), 1);
if isempty(iT), error('ADS file needs a "time" column; found: %s', strjoin(names, ', ')); end
nC = numel(components);
iV = zeros(1, nC);
for c = 1:nC
    m = find(strcmpi(names, "V" + lower(components(c)) + "_prime"), 1);
    if isempty(m)
        fprintf('ADS: no V%s_prime column; E_%s not available\n', lower(components(c)), lower(components(c)));
    else
        iV(c) = m;
    end
end
[thCol, thScale] = angleColumn(names, "theta");
[phCol, phScale] = angleColumn(names, "phi");

if isempty(thCol) && isempty(phCol)
    key = [NaN ads_phi];  grp = ones(size(A, 1), 1);
elseif isempty(thCol)                              % phi sweep, theta unknown
    [kp, ~, grp] = unique(round(A(:, phCol), 6));
    key = [nan(numel(kp), 1), kp * phScale];
else
    cols = thCol;  if ~isempty(phCol), cols = [thCol phCol]; end
    [key, ~, grp] = unique(round(A(:, cols), 6), 'rows');
    key(:, 1) = key(:, 1) * thScale;
    if isempty(phCol)
        key(:, 2) = ads_phi;
    else
        key(:, 2) = key(:, 2) * phScale;
    end
end
ads = struct('theta', {}, 'phi', {}, 'phiInFile', {}, 't', {}, 'e', {});
for k = 1:size(key, 1)
    Ak = A(grp == k, :);
    [t, iu] = unique(Ak(:, iT), 'last');           % drop repeated breakpoints
    th = key(k, 1);  ph = key(k, 2);  sgn = 1;
    if ~isnan(th)
        th = mod(th + 180, 360) - 180;             % -180..180
        if th < 0, th = -th; ph = ph + 180; sgn = -1; end
    end
    ph = mod(ph, 360);
    ads(k).theta = th;
    ads(k).phi   = ph;
    ads(k).phiInFile = ~isempty(phCol);
    ads(k).t = t;
    ads(k).e = nan(numel(t), nC);
    for c = find(iV)
        s = 1;  if any(strcmpi(components(c), ["Theta" "Phi"])), s = sgn; end
        ads(k).e(:, c) = s * Ak(iu, iV(c));
    end
end
if isempty(phCol)
    phTxt = sprintf('%g deg (setting ads_phi; no phi column)', ads_phi);
else
    phTxt = strjoin(compose('%.4g', unique(key(:, 2)).'), ', ') + " deg";
end
if isempty(thCol)
    thTxt = "not in file";
else
    thTxt = strjoin(compose('%.4g', unique(key(:, 1)).'), ', ') + " deg";
end
fprintf('ADS: theta = %s, phi = %s\n', thTxt, phTxt);
if any(key(:, 1) < 0 | key(:, 1) > 180)
    fprintf('ADS: sweeps outside theta 0..180 mapped to (|theta|, phi+180), E_theta/E_phi sign-flipped\n');
end
end

% ------------------------------------------------------------------------
function [col, scale] = angleColumn(names, ang)
% column of an angle: <ang>_deg in degrees, else <ang> in radians
col = find(strcmpi(names, ang + "_deg"), 1);  scale = 1;
if isempty(col)
    col = find(strcmpi(names, ang), 1);  scale = 180/pi;
end
end

function [names, A] = readAdsCsv(file)
% numeric rows into A, column names from the header ("name","name",...)
L = strtrim(readlines(file));
L = L(strlength(L) > 0);
isnum = ~cellfun(@isempty, regexp(cellstr(L), '^"?\s*[-+0-9.]', 'once'));
rows = cellfun(@(s) sscanf(strrep(strrep(s, '"', ''), ',', ' '), '%f').', ...
    cellstr(L(isnum)), 'UniformOutput', false);
A = vertcat(rows{:});
tok = regexp(strjoin(L(~isnum), newline), '(?:^"|",")(\w+)', 'tokens');
names = string(cellfun(@(c) c{1}, tok, 'UniformOutput', false));
if numel(names) ~= size(A, 2)
    if size(A, 2) == 3                              % 3-column export without readable names
        names = ["time" "Vtheta_prime" "Vphi_prime"];
    else
        error('Could not read %d column names from %s.', size(A, 2), file);
    end
end
end
