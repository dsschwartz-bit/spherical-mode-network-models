function [Q, j_, j3D, mmax, nmax, freqs] = ReadSphModes_freqs(flnm)
%READSPHMODES_FREQS  Read spherical mode Q coefficients from a TICRA/GRASP
% .sph file exported by FEKO, supporting single or multi-frequency files.
%
% The modes are stored in the single-index compressed mode index scheme
% (see FEKO documentation). Multiple frequency blocks are detected
% automatically and all are returned.
%
% OUTPUTS:
%   Q     - Complex spherical mode coefficients.
%             Size: [num_j_indices x Nfreq]
%             Q(:, f) contains all mode coefficients for frequency freqs(f).
%             Indexed along the first dimension by the compressed mode
%             index j_; unused j_ indices are zero.
%   j_     - Vector of compressed mode indices populated in the file.
%             (Same for all frequencies; determined by nmax/mmax.)
%   j3D   - String array mapping compressed index j_ -> '(s m n)' label.
%   mmax  - Maximum azimuthal mode order (same across all frequency blocks).
%   nmax  - Maximum elevation mode order (same across all frequency blocks).
%   freqs - Vector of frequencies [Hz], length Nfreq.
%
% USAGE:
%   [Q, j_, j3D, mmax, nmax, freqs] = ReadSphModes_freqs('myantenna.sph');
%   % Access coefficients at the f-th frequency:
%   Q_f = Q(:, f);
%   % Pass to FFfromSph for a single frequency slice:
%   [E_th, E_ph, ~] = FFfromSph(Q_f, mmax, nmax, num_modes, theta, phi, Ntheta, Nphi);
%
% NOTES:
%   - Each frequency block in the .sph file must have the same nmax/mmax.
%     An error is raised if this is violated.
%   - The j_ and j3D outputs are identical across frequency blocks and are
%     only computed once (from the first block).
%
% Authors: DB Davidson, Feb 2024 (original single-frequency version).
%          Updated for multi-frequency support, Mar 2026.

%% --- Read the raw file into a table ---
% Try the given filename first; if not found, toggle "_fr_1" in the name
if exist(flnm, 'file')
    active_flnm = flnm;
else
    if contains(flnm, '_fr_1_FarField')
        % filename has fr_1 -> try without it
        alt_flnm = strrep(flnm, '_fr_1_FarField', '_FarField');
    else
        % filename does not have fr_1 -> try with it
        alt_flnm = regexprep(flnm, '(_FarField[0-9]*)', '_fr_1$1');
    end
    if exist(alt_flnm, 'file')
        fprintf('ReadSphModes_freqs: Could not find %s, using %s instead.\n', flnm, alt_flnm);
        active_flnm = alt_flnm;
    else
        error('Cannot open file: %s (also tried: %s)', flnm, alt_flnm);
    end
end
modedata = importfile_sph(active_flnm);
total_lines = size(modedata, 1);

%% --- Locate the start line of every frequency block ---
% Every block contains a "Frequency = ..." line. The raw file is re-read
% as text to find these lines and parse the frequency values.

fid = fopen(active_flnm, 'r');
raw_lines = {};
while ~feof(fid)
    raw_lines{end+1} = fgetl(fid); %#ok<AGROW>
end
fclose(fid);

% Find line numbers (1-based) of "Frequency =" lines in the raw file
freq_line_indices = [];
freqs = [];
for li = 1:length(raw_lines)
    ln = strtrim(raw_lines{li});
    if startsWith(ln, 'Frequency')
        % Parse the frequency value
        tok = regexp(ln, 'Frequency\s*=\s*([\d.E+\-]+)', 'tokens');
        if ~isempty(tok)
            freqs(end+1) = str2double(tok{1}{1}); %#ok<AGROW>
            freq_line_indices(end+1) = li; %#ok<AGROW>
        end
    end
end
Nfreq = length(freqs);
fprintf('ReadSphModes: found %d frequency block(s).\n', Nfreq);

%% --- Determine the line offset of "line 1 of modedata table" vs raw file ---
% importfile_sph skips the first 2 raw lines (TICRA header + Filename),
% so modedata row 1 corresponds to raw line 3.
raw_offset = 2; % modedata{k,..} corresponds to raw_lines{k + raw_offset}

%% --- Parse nmax/mmax from the first block ---
% In the modedata table, line 3 of the raw file = modedata row 1.
% Columns 3 and 4 hold nmax and mmax.
nmax_global = modedata{1, 3};
mmax_global = modedata{1, 4};
nmax = nmax_global;
mmax = mmax_global;

%% --- Pre-compute total number of compressed mode indices ---
num_modes = 2 * (nmax*(nmax+1) + mmax - 1) + 2;

%% --- Preallocate output arrays ---
% Q is indexed as Q(j_index, freq_index). We use the max possible j_ value
% to size the first dimension.
Q = zeros(num_modes, Nfreq);
j_    = [];
j3D  = strings(1, num_modes);

%% --- Loop over frequency blocks ---
for fi = 1:Nfreq

    % The "Frequency =" line is at raw_lines{freq_line_indices(fi)}.
    % The block header in the raw file is:
    %   freq_line_indices(fi) - 3 : TICRA line        (raw)
    %   freq_line_indices(fi) - 2 : Filename line     (raw)
    %   freq_line_indices(fi) - 1 : params line       (raw)  -> modedata row
    %   freq_line_indices(fi)     : Frequency line    (raw)
    %   freq_line_indices(fi) + 1 : zeros line        (raw)
    %   freq_line_indices(fi) + 2 : zeros line        (raw)
    %   freq_line_indices(fi) + 3 : blank line        (raw)
    %   freq_line_indices(fi) + 4 : blank line        (raw)
    %   freq_line_indices(fi) + 5 : first "m" header  (raw)  <- mode data starts

    % Convert: modedata row = raw_line - raw_offset
    % params line is at raw_line = freq_line_indices(fi) - 1
    params_modedata_row = freq_line_indices(fi) - 1 - raw_offset;

    % Verify nmax/mmax consistency across blocks
    blk_nmax = modedata{params_modedata_row, 3};
    blk_mmax = modedata{params_modedata_row, 4};
    if blk_nmax ~= nmax_global || blk_mmax ~= mmax_global
        error(['ReadSphModes: nmax/mmax mismatch at frequency block ', ...
               num2str(fi), ' (f = ', num2str(freqs(fi)/1e6), ' MHz). ', ...
               'All blocks must have the same modal truncation order.']);
    end

    % First mode data line in modedata: 5 lines after the Frequency line in raw
    % = freq_line_indices(fi) + 5, converted to modedata row index
    first_mode_raw = freq_line_indices(fi) + 5;
    lines = first_mode_raw - raw_offset;  % modedata row of first "m=0" header

    mode_counter = 1;

    %% --- m = 0 block ---
    modecheck = modedata{lines, 1};
    lines = lines + 1;
    m = 0;
    if modecheck ~= m
        error(['ReadSphModes: m=0 header not found where expected in block ', ...
               num2str(fi)]);
    end

    for n = 1:nmax
        for s = 1:2  % s=1: TE,  s=2: TM
            jidx = 2*(n*(n+1) + m - 1) + s;  % compressed mode index
            if fi == 1
                j_(mode_counter) = jidx; %#ok<AGROW>
                j3D(jidx) = string([num2str(s), ' ', num2str(m), ' ', num2str(n)]);
            end
            if s == 1
                Q(jidx, fi) = complex(modedata{lines, 1}, modedata{lines, 2});
            else
                Q(jidx, fi) = complex(modedata{lines, 3}, modedata{lines, 4});
            end
            mode_counter = mode_counter + 1;
        end
        lines = lines + 1;
    end

    %% --- m = 1 .. mmax blocks ---
    for mmode = 1:mmax
        modecheck = modedata{lines, 1};
        if modecheck ~= mmode
            error(['ReadSphModes: expected m=', num2str(mmode), ...
                   ' header in frequency block ', num2str(fi), ...
                   ', got ', num2str(modecheck)]);
        end
        lines = lines + 1;  % advance past the m-header line

        for n = abs(mmode):nmax
            for m = -mmode:2*mmode:mmode  % iterates m = -mmode then +mmode
                for s = 1:2  % s=1: TE,  s=2: TM
                    jidx = 2*(n*(n+1) + m - 1) + s;
                    if fi == 1
                        j_(mode_counter) = jidx; %#ok<AGROW>
                        j3D(jidx) = string([num2str(s), ' ', num2str(m), ' ', num2str(n)]);
                    end
                    if s == 1
                        Q(jidx, fi) = complex(modedata{lines, 1}, modedata{lines, 2});
                    else
                        Q(jidx, fi) = complex(modedata{lines, 3}, modedata{lines, 4});
                    end
                    mode_counter = mode_counter + 1;
                end
                lines = lines + 1;
            end
        end
    end

end  % frequency loop

%% --- Sort frequency dimension ascending ---
[freqs, sort_idx] = sort(freqs, 'ascend');
Q = Q(:, sort_idx);

fprintf('ReadSphModes: successfully read %d modes x %d frequencies.\n', ...
        length(j_), Nfreq);

end