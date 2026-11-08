function [Eptr, r, theta, phi, freqs] = ReadFE_spherical(filename)
% ReadFE_spherical  Read a FEKO near-field file (.efe or .hfe) on a spherical
% grid, with one "#Frequency" block per frequency.
%
% File columns: Radius Theta Phi Re/Im(Fr) Re/Im(Fth) Re/Im(Fph), sampled
% with Phi slowest, Theta next, Radius fastest; angles in degrees.
%
% Output Eptr is sized (NTheta x NPhi x NRadius x 3 x NFreq):
%   dim 4 = field component: 1 = phi, 2 = theta, 3 = r
%   dim 5 = frequency, ascending
% r     - (NRadius x 1) radii (m)
% theta - (NTheta  x 1) theta (rad)
% phi   - (NPhi    x 1) phi (rad)
% freqs - (NFreq   x 1) frequencies (Hz), ascending
%
% If filename does not exist, the name with '_fr_1' toggled before
% '_NearField' is tried (see resolve_feko_filename).
%
% Usage:
%   [Eptr, r, theta, phi, freqs] = ReadFE_spherical('dipole_cylinder_NearField1.efe')

% Try given filename first; if not found, toggle fr_1 in filename
fid = fopen(filename, 'r');
if fid < 0
    if contains(filename, '_fr_1_NearField')
        % filename has fr_1 -> try without it
        alt_filename = strrep(filename, '_fr_1_NearField', '_NearField');
    else
        % filename does not have fr_1 -> try with it
        alt_filename = regexprep(filename, '(_NearField[0-9]*)', '_fr_1$1');
    end
    fid = fopen(alt_filename, 'r');
    if fid < 0
        error('Cannot open file: %s (also tried: %s)', filename, alt_filename);
    else
        fprintf('ReadFEKO: Could not open %s, using %s instead.\n', filename, alt_filename);
    end
end

% Read all lines
alllines = {};
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    alllines{end+1} = strtrim(line);
end
fclose(fid);

nlines = length(alllines);

% --- Find all solution block start lines (lines starting with #Frequency) ---
freq_lineidx = [];
for i = 1:nlines
    if startsWith(alllines{i}, '#Frequency')
        freq_lineidx(end+1) = i;
    end
end
NFreq = length(freq_lineidx);
freqs = zeros(NFreq, 1);

% sentinel
freq_lineidx(end+1) = nlines + 1;

% --- Parse first solution block header to get dimensions ---
NRadius = 1; NTheta = 1; NPhi = 1;

for i = freq_lineidx(1):freq_lineidx(2)-1
    l = alllines{i};
    if isempty(l) || l(1) ~= '#', break; end
    content = strrep(strtrim(l(2:end)), "'", '');
    if startsWith(content,'Frequency:')
        freqs(1) = str2double(strtrim(extractAfter(content,'Frequency:')));
    elseif startsWith(content,'No. of Radius')
        tok = regexp(content, 'No\. of Radius''?\s*Samples:\s*([\d]+)', 'tokens');
        if ~isempty(tok), NRadius = str2double(tok{1}{1}); end
    elseif startsWith(content,'No. of Theta')
        tok = regexp(content, 'No\. of Theta''?\s*Samples:\s*([\d]+)', 'tokens');
        if ~isempty(tok), NTheta = str2double(tok{1}{1}); end
    elseif startsWith(content,'No. of Phi')
        tok = regexp(content, 'No\. of Phi''?\s*Samples:\s*([\d]+)', 'tokens');
        if ~isempty(tok), NPhi = str2double(tok{1}{1}); end
    end
end

NpointsPerFreq = NRadius * NTheta * NPhi;

% Preallocate outputs
Eptr      = zeros(NTheta, NPhi, NRadius, 3, NFreq);
r         = zeros(NRadius, 1);
theta_deg = zeros(NTheta,  1);  % degrees, converted to radians at end
phi_deg   = zeros(NPhi,    1);  % degrees, converted to radians at end
coords_set = false;

% --- Parse each frequency block ---
for fi = 1:NFreq

    % get frequency from block header
    for i = freq_lineidx(fi):freq_lineidx(fi+1)-1
        l = alllines{i};
        if isempty(l) || l(1) ~= '#', break; end
        content = strrep(strtrim(l(2:end)), "'", '');
        if startsWith(content,'Frequency:')
            freqs(fi) = str2double(strtrim(extractAfter(content,'Frequency:')));
        end
    end

    % collect numeric data lines for this block
    datalines = {};
    for i = freq_lineidx(fi):freq_lineidx(fi+1)-1
        l = alllines{i};
        if isempty(l) || l(1) == '#' || l(1) == '*', continue; end
        datalines{end+1} = l;
    end

    if length(datalines) ~= NpointsPerFreq
        warning('Block %d: expected %d data lines, got %d', fi, NpointsPerFreq, length(datalines));
    end

    % Parse numeric data
    % Columns: Radius Theta Phi Re(Er) Im(Er) Re(Eth) Im(Eth) Re(Eph) Im(Eph)
    data = zeros(length(datalines), 9);
    for k = 1:length(datalines)
        data(k,:) = sscanf(datalines{k}, '%f', [1,9]);
    end

    % Extract coordinate arrays from first frequency block only
    if ~coords_set
        % Radius varies fastest: first NRadius rows all share Theta(1), Phi(1)
        r         = data(1:NRadius, 1);
        % Theta: one value per NRadius-row group
        theta_deg = data(1:NRadius:NRadius*NTheta, 2);
        % Phi: one value per NRadius*NTheta-row group
        if NPhi > 1
            phi_deg = data(1:NRadius*NTheta:end, 3);
        else
            phi_deg = data(1, 3);
        end
        coords_set = true;
    end

    % Build complex field vectors
    Er  = complex(data(:,4), data(:,5));
    Eth = complex(data(:,6), data(:,7));
    Eph = complex(data(:,8), data(:,9));

    % File order: Phi slowest, Theta next, Radius fastest
    % reshape to (NRadius x NTheta x NPhi) then permute to (NTheta x NPhi x NRadius)
    Er  = permute(reshape(Er,  NRadius, NTheta, NPhi), [2,3,1]);
    Eth = permute(reshape(Eth, NRadius, NTheta, NPhi), [2,3,1]);
    Eph = permute(reshape(Eph, NRadius, NTheta, NPhi), [2,3,1]);

    % Store all radius slices: dim4 = [Eph, Eth, Er]
    Eptr(:,:,:,1,fi) = Eph;
    Eptr(:,:,:,2,fi) = Eth;
    Eptr(:,:,:,3,fi) = Er;

end


% Sort frequency dimension ascending
[freqs, sort_idx] = sort(freqs, 'ascend');
Eptr = Eptr(:,:,:,:,sort_idx);

% Convert angles from degrees to radians
theta = theta_deg * pi/180;
phi   = phi_deg   * pi/180;
end