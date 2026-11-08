function [rEpt_FF, G, theta, phi, freqs, resultType] = ReadFF_spherical(filename)
% ReadFF_spherical  Read FEKO far field file (.ffe) with multiple frequency blocks.
%
% Sampling order in file: Phi slowest, Theta fastest (no Radius in far field).
%
% Output rEpt_FF is sized (NTheta x NPhi x 2 x NFreq), where:
%   dim 1 = theta index
%   dim 2 = phi index
%   dim 3 = field vector component: 1=Eph, 2=Eth
%   dim 4 = frequency index, ascending
%   FEKO's far field is r*E*exp(+jkr) (V): the 1/r decay and exp(-jkr)
%   propagation phase are removed, hence the "r" prefix.
%
% Output G is sized (NTheta x NPhi x 3 x NFreq):
%   dim 3 = gain or directivity: 1 = phi, 2 = theta, 3 = total (dB, as
%   exported; same component order as rEpt_FF)
%
% theta      - (NTheta x 1) theta (rad)
% phi        - (NPhi   x 1) phi (rad)
% freqs      - (NFreq  x 1) frequencies (Hz), ascending
% resultType - 'Directivity' or 'Gain', from the file's "Result Type" line
%
% If filename does not exist, the name with '_fr_1' toggled before
% '_FarField' is tried (see resolve_feko_filename).
%
% Usage:
%   [rEpt_FF, G, theta, phi, freqs, resultType] = ReadFF_spherical('dipole_cylinder_FarField1.ffe')

% Try given filename first; if not found, toggle fr_1 in filename
fid = fopen(filename, 'r');
if fid < 0
    if contains(filename, '_fr_1_FarField')
        % filename has fr_1 -> try without it
        alt_filename = strrep(filename, '_fr_1_FarField', '_FarField');
    else
        % filename does not have fr_1 -> try with it
        alt_filename = regexprep(filename, '(_FarField[0-9]*)', '_fr_1$1');
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

% --- Parse first solution block header to get dimensions / result type ---
NTheta = 1; NPhi = 1;
resultType = '';

for i = freq_lineidx(1):freq_lineidx(2)-1
    l = alllines{i};
    if isempty(l) || l(1) ~= '#', break; end
    content = strrep(strtrim(l(2:end)), "'", '');
    if startsWith(content,'Frequency:')
        freqs(1) = str2double(strtrim(extractAfter(content,'Frequency:')));
    elseif startsWith(content,'No. of Theta')
        tok = regexp(content, 'No\. of Theta''?\s*Samples:\s*([\d]+)', 'tokens');
        if ~isempty(tok), NTheta = str2double(tok{1}{1}); end
    elseif startsWith(content,'No. of Phi')
        tok = regexp(content, 'No\. of Phi''?\s*Samples:\s*([\d]+)', 'tokens');
        if ~isempty(tok), NPhi = str2double(tok{1}{1}); end
    elseif startsWith(content,'Result Type:')
        resultType = strtrim(extractAfter(content,'Result Type:'));
    end
end

NpointsPerFreq = NTheta * NPhi;

% Preallocate outputs
rEpt_FF      = zeros(NTheta, NPhi, 2, NFreq);
G         = zeros(NTheta, NPhi, 3, NFreq);
theta_deg = zeros(NTheta, 1);  % degrees, converted to radians at end
phi_deg   = zeros(NPhi,   1);  % degrees, converted to radians at end
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
    % Columns: Theta Phi Re(Eth) Im(Eth) Re(Eph) Im(Eph) Dir/Gain(Theta) Dir/Gain(Phi) Dir/Gain(Total)
    data = zeros(length(datalines), 9);
    for k = 1:length(datalines)
        data(k,:) = sscanf(datalines{k}, '%f', [1,9]);
    end

    % Extract coordinate arrays from first frequency block only
    if ~coords_set
        % Theta varies fastest: first NTheta rows all share Phi(1)
        theta_deg = data(1:NTheta, 1);
        % Phi: one value per NTheta-row group
        if NPhi > 1
            phi_deg = data(1:NTheta:end, 2);
        else
            phi_deg = data(1, 2);
        end
        coords_set = true;
    end

    % Build complex field vectors
    Eth = complex(data(:,3), data(:,4));
    Eph = complex(data(:,5), data(:,6));

    % Gain/Directivity columns (real-valued, in dB as exported by FEKO)
    Gth  = data(:,7);
    Gph  = data(:,8);
    Gtot = data(:,9);

    % File order: Phi slowest, Theta fastest
    % reshape to (NTheta x NPhi)
    Eth  = reshape(Eth,  NTheta, NPhi);
    Eph  = reshape(Eph,  NTheta, NPhi);
    Gth  = reshape(Gth,  NTheta, NPhi);
    Gph  = reshape(Gph,  NTheta, NPhi);
    Gtot = reshape(Gtot, NTheta, NPhi);

    % Store field: dim3 = [Eph, Eth]
    rEpt_FF(:,:,1,fi) = Eph;
    rEpt_FF(:,:,2,fi) = Eth;

    % Store gain/directivity: dim3 = [Phi, Theta, Total]
    G(:,:,1,fi) = Gph;
    G(:,:,2,fi) = Gth;
    G(:,:,3,fi) = Gtot;

end

% Sort frequency dimension ascending
[freqs, sort_idx] = sort(freqs, 'ascend');
rEpt_FF = rEpt_FF(:,:,:,sort_idx);
G    = G(:,:,:,sort_idx);

% Convert angles from degrees to radians
theta = theta_deg * pi/180;
phi   = phi_deg   * pi/180;
end