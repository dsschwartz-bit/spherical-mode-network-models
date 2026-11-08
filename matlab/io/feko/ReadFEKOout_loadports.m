function [V0, Z0_ref, R0, symmetry_planes, even_odd_excitation, antenna_apothem, Freq, port_pol, secondary_ports, S, fields_config_name, near_field_name, far_field_name, primary_sparam_index, near_field_occurrence, far_field_occurrence, near_field_type, far_field_type, sparam_diagnostic] = ReadFEKOout_loadports(filename)
% ReadFEKOout_loadports  Parse a multi-configuration FEKO .out file: source/port
% data of the driven ("fields") configuration, the S-parameter matrix, and the
% secondary (loaded) ports.
%
% Extends ReadFEKOout to projects with several CONFIGURATION blocks:
%   - fields configuration: one port driven by a voltage source, the other
%     ports terminated in their physical loads; near/far fields computed.
%   - S-parameter configuration: every port excited in turn (others terminated
%     in the reference impedance); prints "S sink source re im ..." rows.
% Configurations are classified by content, not by name (classify_configurations):
% exactly one distinct excited port -> fields; >= 2 excited ports and S-matrix
% rows -> S-parameter (the one with most ports is used). If there are several
% fields configurations (one per driven port), the one driving the port with the
% lowest trailing number (EdgePort1 before EdgePort2) is the PRIMARY; the others
% are parsed by parse_one_field_config and merged into secondary_ports.
%
% Outputs (V0 ... port_pol as in ReadFEKOout, taken from the PRIMARY fields
% configuration):
%   V0                  - (Nports x 1) complex source voltage |U0|*exp(j*ARG(U0));
%                         doubled if an infinite PEC/PMC ground plane coincides
%                         with the port
%   Z0_ref              - reference impedance (Ohm); x2 for a coincident PEC
%                         ground plane, /2 for a coincident PMC ground plane
%   R0                  - 1x3 port location (m): mean of the wire-port locations
%                         or of the edge-port mesh-edge midpoints
%   symmetry_planes     - struct .XY/.XZ/.YZ (+1 PMC wall, 0 none, -1 PEC wall)
%                         and .physical_gnd (1 if an infinite ground plane
%                         coincides with the port)
%   even_odd_excitation - +log2(Nports) for ports mirrored across a magnetic wall
%                         (even mode), -log2(Nports) across an electric wall
%                         (odd mode), 0 otherwise
%   antenna_apothem     - max |r| over the mesh vertices (m), about the origin;
%                         requested interactively if the file has no mesh
%   Freq                - (F x 1) Hz: the S-parameter sweep (indexes S), or, if
%                         there is no S-parameter configuration, every FREQ value
%                         found in the file
%   port_pol            - (Nports x 3) port polarity unit vectors
%   secondary_ports     - struct array, one entry per port other than the
%                         primary port (fields below); 0x0 if none
%   S                   - (N x N x F) complex, S(sink, source, f) as printed by
%                         FEKO, normalised to the S-parameter reference
%                         impedances; [] if there is no S-parameter configuration.
%                         Input impedance: Zin_n = Z0*(1+S(n,n,:))./(1-S(n,n,:))
%   fields_config_name  - name of the primary fields configuration
%   near_field_name     - raw near-field request name of the primary configuration
%   far_field_name      - raw far-field / spherical-mode request name ('' if none)
%   near_field_occurrence, far_field_occurrence - number of earlier
%                         configurations that used the same request name
%                         (0 = first use)
%   near_field_type, far_field_type - '' (first use), 'same_freq' or
%                         'diff_freq': whether this configuration's sweep equals
%                         that of the first configuration using the name
%   primary_sparam_index - primary port's row/column in S ([] if not found)
%   sparam_diagnostic   - reason why S / primary_sparam_index are empty ('' if not)
%
% Exported field file names: the first configuration using a request name gets
% <prefix>_<name>.<ext> (.efe/.hfe/.ffe/.sph). Later reuses get
% <prefix>_<name>(k).<ext> (same sweep) or <prefix>_fr_<k>_<name>.<ext>
% (different sweep). resolve_field_file tries these candidates on disk.
%
% secondary_ports fields:
%   .name, .edge_indices (signed), .face_names {negative, positive}
%   .R2, .pol2          - 1x3 location / polarity (NaN if unresolved); from the
%                         port's own fields configuration if it has one, else
%                         from the mesh using the S-parameter edge definition
%   .Z_sparam_ref       - S-parameter reference impedance (Ohm)
%   .Z_field_load       - impedance of the load on this port in the primary
%                         fields configuration (NaN if none found)
%   .edges_match        - S-parameter edge set equals that load's edge set
%   .impedance_match    - Z_field_load equals Z_sparam_ref within 0.1%
%   .has_own_field_config, .field_config_name, .V0_own, .Z0_ref_own
%                       - from the port's own fields configuration, if any
%   .geometry_match     - own-configuration geometry agrees with the
%                         S-parameter/mesh geometry (NaN if not compared)
%   .edges_match_field_config - S-parameter edge set equals the own fields
%                         configuration's edge set (NaN if not available)
%   .impedance_match_own - Z_sparam_ref, Z_field_load and Z0_ref_own agree within
%                         0.1% (NaN if fewer than two are available)
%   .reflection_match   - S(i,i,:) agrees with (Zin-Z0_ref_own)/(Zin+Z0_ref_own)
%                         within 1% at every frequency (NaN if not checked)
%   .near_field_name, .near_field_occurrence, .near_field_type,
%   .far_field_name, .far_field_occurrence, .far_field_type - as for the primary
%   .sparam_port_index  - row/column in S (NaN if not in the S-parameter config)
%
% Printed checks for the primary port (no outputs): edge set and reference
% impedance against the S-parameter configuration, and
% Gamma = (Zin-Z0_ref)/(Zin+Z0_ref) against S(i,i,:) within 1% (complex
% difference). Zin is parsed from the "DATA OF THE VOLTAGE SOURCE" blocks for
% these checks only and is not returned.
%
% Usage:
%   [V0, Z0_ref, R0, symmetry_planes, even_odd_excitation, antenna_apothem, ...
%    Freq, port_pol, secondary_ports, S] = ReadFEKOout_loadports('loop_loadport.out')
%
% See also ReadFEKOout, resolve_field_file.

% --- Pre-scan: classify configurations (fields vs S-parameter) ---
[cfg_info, cfg_order] = classify_configurations(filename);

% --- Field-request name reuse across configurations ---
% For each configuration record the raw near/far-field request name, its 0-based
% occurrence index and the reuse type ('same_freq' / 'diff_freq', comparing its
% sweep with that of the first configuration using the name). These decide the
% exported file name (see the header); resolve_field_file finds the file.
near_name_raw   = containers.Map('KeyType','char','ValueType','any');    % config name -> raw near-field request name
far_name_raw    = containers.Map('KeyType','char','ValueType','any');    % config name -> raw far-field request name
near_name_occ   = containers.Map('KeyType','char','ValueType','double'); % config name -> 0-based occurrence index
far_name_occ    = containers.Map('KeyType','char','ValueType','double');
near_name_type  = containers.Map('KeyType','char','ValueType','any');    % config name -> '' | 'same_freq' | 'diff_freq'
far_name_type   = containers.Map('KeyType','char','ValueType','any');
near_first_cfg  = containers.Map('KeyType','char','ValueType','char');   % raw name -> name of the FIRST configuration to use it
far_first_cfg   = containers.Map('KeyType','char','ValueType','char');
near_seen_count = containers.Map('KeyType','char','ValueType','double');
far_seen_count  = containers.Map('KeyType','char','ValueType','double');

disamb_freq_tol = 1e-6;   % relative tolerance for comparing two frequency sweeps
    function tf = disamb_freqs_match(fa, fb)
        tf = (numel(fa) == numel(fb)) && ...
             (isempty(fa) || all(abs(fa(:) - fb(:)) <= disamb_freq_tol * max(abs(fb(:)), 1)));
    end

for ci = 1:numel(cfg_order)
    cname = cfg_order{ci};
    info  = cfg_info(cname);
    if ~isempty(info.near_field_name)
        raw = info.near_field_name;
        near_name_raw(cname) = raw;
        if isKey(near_seen_count, raw)
            cnt = near_seen_count(raw);
            near_name_occ(cname) = cnt;
            near_seen_count(raw) = cnt + 1;
            first_cname = near_first_cfg(raw);
            if disamb_freqs_match(info.freqs, cfg_info(first_cname).freqs)
                near_name_type(cname) = 'same_freq';
            else
                near_name_type(cname) = 'diff_freq';
            end
        else
            near_name_occ(cname) = 0;
            near_seen_count(raw) = 1;
            near_first_cfg(raw) = cname;
            near_name_type(cname) = '';
        end
    else
        near_name_raw(cname) = '';
        near_name_occ(cname) = 0;
        near_name_type(cname) = '';
    end
    if ~isempty(info.far_field_name)
        raw = info.far_field_name;
        far_name_raw(cname) = raw;
        if isKey(far_seen_count, raw)
            cnt = far_seen_count(raw);
            far_name_occ(cname) = cnt;
            far_seen_count(raw) = cnt + 1;
            first_cname = far_first_cfg(raw);
            if disamb_freqs_match(info.freqs, cfg_info(first_cname).freqs)
                far_name_type(cname) = 'same_freq';
            else
                far_name_type(cname) = 'diff_freq';
            end
        else
            far_name_occ(cname) = 0;
            far_seen_count(raw) = 1;
            far_first_cfg(raw) = cname;
            far_name_type(cname) = '';
        end
    else
        far_name_raw(cname) = '';
        far_name_occ(cname) = 0;
        far_name_type(cname) = '';
    end
end

fields_config_names = {};   % all single-driven-port configs, in file order
fields_config_ports = {};   % the port each one drives, same order
sparam_config_name  = '';
best_sparam_ports    = 1;   % need >= 2 to qualify
single_port_sparam_name = '';   % first config with S-matrix rows but only 1 port (diagnostic only)
% The S-parameter configuration needs >= 2 excited ports AND S-matrix rows:
% some projects run their per-port field solves as sequential excitations in
% one configuration, which has the same port count but no coupling data.
for ci = 1:numel(cfg_order)
    cname = cfg_order{ci};
    info = cfg_info(cname);
    if info.n_excited_ports == 1
        fields_config_names{end+1} = cname; %#ok<AGROW>
        if ~isempty(info.port_names)
            fields_config_ports{end+1} = info.port_names{1}; %#ok<AGROW>
        else
            fields_config_ports{end+1} = ''; %#ok<AGROW>
        end
        if info.has_sparam_data && isempty(single_port_sparam_name)
            single_port_sparam_name = cname;
        end
    elseif info.n_excited_ports >= 2 && info.has_sparam_data && info.n_excited_ports > best_sparam_ports
        best_sparam_ports = info.n_excited_ports;
        sparam_config_name = cname;
    end
end

% Sort the fields configurations by their driven port's trailing number
% (e.g. 'EdgePort1' -> 1, 'EdgePort2' -> 2); ties or non-numeric names fall
% back to file order (stable sort).
if ~isempty(fields_config_names)
    port_nums = inf(1, numel(fields_config_names));
    for pi_ = 1:numel(fields_config_names)
        num_tok = regexp(fields_config_ports{pi_}, '(\d+)$', 'tokens');
        if ~isempty(num_tok)
            port_nums(pi_) = str2double(num_tok{1}{1});
        end
    end
    [~, sort_idx] = sort(port_nums);   % sort() is stable in MATLAB/Octave
    fields_config_names = fields_config_names(sort_idx);
    fields_config_ports = fields_config_ports(sort_idx);
end

fields_config_name = '';
near_field_name = '';
far_field_name  = '';
near_field_occurrence = 0;
far_field_occurrence  = 0;
near_field_type = '';
far_field_type  = '';
if isempty(fields_config_names)
    % Fallback: no single-excited-port configuration found; use the first
    % configuration in the file (as ReadFEKOout does).
    if ~isempty(cfg_order)
        fields_config_name = cfg_order{1};
    end
    fprintf('ReadFEKOout_loadports: WARNING - could not uniquely identify a fields configuration (single driven port). Falling back to "%s".\n', fields_config_name);
else
    fields_config_name = fields_config_names{1};
    info = cfg_info(fields_config_name);
    near_field_name       = near_name_raw(fields_config_name);
    far_field_name        = far_name_raw(fields_config_name);
    near_field_occurrence = near_name_occ(fields_config_name);
    far_field_occurrence  = far_name_occ(fields_config_name);
    near_field_type       = near_name_type(fields_config_name);
    far_field_type        = far_name_type(fields_config_name);
    fprintf('ReadFEKOout_loadports: Configuration "%s" identified as the PRIMARY fields configuration (assigns sources + computes fields): %d driven port(s) [%s], field-calculation request found = %d.\n', ...
        fields_config_name, info.n_excited_ports, strjoin(info.port_names, ', '), info.has_field_calc);
    if ~isempty(near_field_name) || ~isempty(far_field_name)
        near_type_str = near_field_type; if isempty(near_type_str), near_type_str = 'first use'; end
        far_type_str  = far_field_type;  if isempty(far_type_str),  far_type_str  = 'first use'; end
        fprintf('ReadFEKOout_loadports: Primary port field request name(s): near-field "%s" (occurrence %d, %s), far-field/spherical-mode "%s" (occurrence %d, %s). See near_field_occurrence/far_field_occurrence/near_field_type/far_field_type output doc for how to resolve the actual file on disk -- occurrence 0 is the plain <prefix>_%s.efe/.hfe/.ffe/.sph name; occurrence >=1 needs the caller to try candidate file names based on type.\n', ...
            near_field_name, near_field_occurrence, near_type_str, far_field_name, far_field_occurrence, far_type_str, near_field_name);
    end
    if numel(fields_config_names) > 1
        other_names = fields_config_names(2:end);
        fprintf('ReadFEKOout_loadports: %d additional single-port fields configuration(s) found: %s. Each will be fully parsed for its own V0/geometry.\n', ...
            numel(other_names), strjoin(other_names, ', '));
    end
end
sparam_diagnostic = '';   % set below whenever S / primary_sparam_index will be empty
if ~isempty(sparam_config_name)
    info = cfg_info(sparam_config_name);
    fprintf('ReadFEKOout_loadports: Configuration "%s" identified as the S-PARAMETER configuration: %d sequentially-excited port(s) [%s].\n', ...
        sparam_config_name, info.n_excited_ports, strjoin(info.port_names, ', '));
elseif ~isempty(single_port_sparam_name)
    info = cfg_info(single_port_sparam_name);
    port_str = 'unknown port';
    if ~isempty(info.port_names), port_str = info.port_names{1}; end
    sparam_diagnostic = sprintf(['Configuration "%s" in this .out file is an S-parameter configuration, ' ...
        'but it only has 1 port (%s) -- there is no second port for it to couple to, so this file has no ' ...
        'multi-port S-parameter data for this script to extract. Use general_vector_spherical_modes_FEKO.m ' ...
        '(the non-loadports, single-port script) instead.'], single_port_sparam_name, port_str);
    fprintf('ReadFEKOout_loadports: %s\n', sparam_diagnostic);
else
    all_names = strjoin(cfg_order, ', ');
    if isempty(all_names), all_names = '(none)'; end
    sparam_diagnostic = sprintf(['No configuration in this .out file printed an S-parameter matrix (a ' ...
        '"S  sink  source  re  im ..." row), so there is no port-to-port coupling data to build the ' ...
        'internal-port S-matrix from. Configuration(s) found in the file: %s. Check that the FEKO project ' ...
        'includes an S-parameter request (2 or more ports, sequentially excited) -- a plain field/impedance ' ...
        'solve does not produce one.'], all_names);
    fprintf('ReadFEKOout_loadports: %s\n', sparam_diagnostic);
end

% --- Full S-parameter matrix (all frequencies) ---
S = [];
Sfreq_internal = [];
if ~isempty(sparam_config_name)
    N_sparam_ports = cfg_info(sparam_config_name).n_excited_ports;
    [S, Sfreq_internal] = parse_sparam_matrix(filename, sparam_config_name, N_sparam_ports);
    fprintf('ReadFEKOout_loadports: Parsed %d x %d S-parameter matrix across %d frequency point(s) from "%s".\n', ...
        N_sparam_ports, N_sparam_ports, numel(Sfreq_internal), sparam_config_name);
else
    fprintf('ReadFEKOout_loadports: S output will be empty (no S-parameter configuration found).\n');
end

secondary_ports = struct('name', {}, 'edge_indices', {}, 'face_names', {}, ...
    'R2', {}, 'pol2', {}, 'Z_sparam_ref', {}, 'Z_field_load', {}, ...
    'edges_match', {}, 'impedance_match', {}, 'has_own_field_config', {}, ...
    'field_config_name', {}, 'V0_own', {}, 'Z0_ref_own', {}, ...
    'geometry_match', {}, 'edges_match_field_config', {}, ...
    'impedance_match_own', {}, 'near_field_name', {}, 'near_field_occurrence', {}, 'near_field_type', {}, ...
    'far_field_name', {}, 'far_field_occurrence', {}, 'far_field_type', {}, 'sparam_port_index', {}, 'reflection_match', {});

fid = fopen(filename, 'r');
if fid < 0
    error('Cannot open file: %s', filename);
end

V0_mag_all       = [];    % accumulate per-port |U0|
V0_phase_all     = [];    % accumulate per-port ARG(U0) in degrees
cur_V0_mag       = NaN;
cur_V0_phase     = NaN;
Z0_ref       = NaN;
port_pol     = [0 0 1];   % default polarity: +z; becomes (Nports x 3) after parsing
R0           = zeros(1,3);

has_pec_plane    = false;
even_odd_excitation = 0;      % 0 = no symmetry mode, else log2(num_ports)
n_mag_walls    = 0;      % number of magnetic conducting walls
n_elec_walls   = 0;      % number of electric conducting walls
sym_planes     = {};     % list of detected symmetry planes: struct with .plane and .type
gnd_z         = NaN;
has_pmc_plane = false;
pmc_z         = NaN;
gnd_next_line = false;
pmc_double_V0 = false;  % set true if PMC ground plane coincides with in-plane wire port

port_type        = '';    % 'wire' or 'edge'
in_fields_config = false; % true only while inside the identified fields configuration
excitation_type = '';    % e.g. 'Voltage Source at a Segment', 'Voltage Source at a Wire-Node'
port_locations  = [];    % accumulate all port locations (Nx3)
edge_midpoints  = zeros(0,3);  % edge-port locations from the mesh (stay empty without a mesh)
edge_indices    = {};    % cell array of edge index groups, one cell per port
edge_face_names = {};   % cell array of {neg_face, pos_face} strings, one cell per port
edge_port_names = {};   % cell array of port names (e.g. 'EdgePort1'), one per group
face_line_pending = false;
first_port_done  = false; % tracks whether at least one port has been seen
all_ports_done   = false;
current_R0       = NaN(1,3);
in_source_block  = false;
first_freq_val   = NaN;   % frequency value of first source block
seen_src_index1  = false; % tracks if source index 1 has been seen once already

% --- Load blocks (DATA FOR EDGE LOADS), fields configuration only ---
in_load_section  = false;
load_face_pending = false;
field_loads = struct('name', {}, 'index', {}, 'circuit_type', {}, 'Z', {}, ...
    'attached_port', {}, 'edge_indices', {}, 'face_names', {});

% --- Frequency tracking: every FREQ value in the file (Freq fallback, Zin index) ---
freq_list        = [];    % distinct FREQ values, in order of appearance
last_freq        = NaN;   % most recently parsed FREQ = value

% --- Primary-port Zin, used only to validate S (check_reflection); not returned ---
zin_data          = {};    % cell: {freq_idx, complex Zin} (primary port)
in_vsource_block  = false; % true when inside a DATA OF THE VOLTAGE SOURCE block
cur_freq_idx      = 0;     % current frequency index for zin_data

% --- Single pass through file ---
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);

    % --- Detect symmetry planes (collect all; mode determined after port count known) ---
    if contains(line, 'magnetic conducting wall') || contains(line, 'electric conducting wall')
        wall_type_str = 'magnetic';
        if contains(line, 'electric'), wall_type_str = 'electric'; end
        pm = regexp(line, '(\S+-plane)', 'tokens');
        plane_name = 'unknown-plane';
        if ~isempty(pm), plane_name = pm{1}{1}; end
        % Determine plane normal: YZ->x, XZ->y, XY->z
        if contains(plane_name, 'YZ'), plane_normal = [1 0 0];
        elseif contains(plane_name, 'XZ'), plane_normal = [0 1 0];
        elseif contains(plane_name, 'XY'), plane_normal = [0 0 1];
        else, plane_normal = [0 0 0];
        end
        sym_planes{end+1} = struct('plane', plane_name, 'type', wall_type_str, 'normal', plane_normal);
        if strcmp(wall_type_str, 'magnetic'), n_mag_walls = n_mag_walls + 1;
        else, n_elec_walls = n_elec_walls + 1; end
    end

    % --- Detect ground plane ---
    if contains(line, 'AN INFINITE GROUND PLANE IS PRESENT')
        gnd_next_line = true;
    elseif gnd_next_line
        if isempty(line), continue; end
        gnd_next_line = false;
        if contains(line, 'Perfect electric conductor plane')
            tok = regexp(line, 'z\s*=\s*([\dE+\-\.]+)', 'tokens');
            if ~isempty(tok)
                has_pec_plane = true;
                gnd_z = str2double(tok{1}{1});
            end
        elseif contains(line, 'Perfect magnetic conductor plane')
            tok = regexp(line, 'z\s*=\s*([\dE+\-\.]+)', 'tokens');
            if ~isempty(tok)
                has_pmc_plane = true;
                pmc_z = str2double(tok{1}{1});
            end
        end
    end

    % --- Detect configuration block: track whether we're inside the fields configuration ---
    if contains(line, 'CONFIGURATION:')
        tok = regexp(line, 'CONFIGURATION:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            this_config = tok{1}{1};
            if strcmp(this_config, fields_config_name)
                in_fields_config = true;
            else
                if in_fields_config && first_port_done
                    % We've already collected the fields-config port(s) and are
                    % now leaving that block for good (its first occurrence is done).
                    all_ports_done = true;
                end
                in_fields_config = false;
            end
        end
    end

    % --- Stop if load data section starts (S-param solve pass) ---
    if in_fields_config && (contains(line, 'DATA FOR VERTEX LOADS') || contains(line, 'READING SOLUTION COEFFICIENTS'))
        all_ports_done = true;
    end

% --- Current frequency. A new frequency also ends port collection (for files
% without a 'READING SOLUTION COEFFICIENTS' line). ---
    freq_m = regexp(line, 'FREQ\s*=\s*([\dE+\-.]+)', 'tokens');
    if ~isempty(freq_m)
        last_freq = str2double(freq_m{1}{1});
        if in_fields_config && ~all_ports_done
            if isnan(first_freq_val)
                first_freq_val = last_freq;
            elseif abs(last_freq - first_freq_val) > first_freq_val * 1e-6
                all_ports_done = true;
            end
        end
        fidx = find(abs(freq_list - last_freq) < last_freq*1e-6, 1);
        if isempty(fidx)
            freq_list(end+1) = last_freq;
            cur_freq_idx = length(freq_list);
        else
            cur_freq_idx = fidx;
        end
    end

% --- Primary-port Zin from DATA OF THE VOLTAGE SOURCE blocks, at every
% frequency of the fields configuration ---
    if in_fields_config
        src_m = regexp(line, 'DATA OF THE VOLTAGE SOURCE NO\.\s*(\d+)', 'tokens');
        if ~isempty(src_m)
            in_vsource_block = true;
        elseif in_vsource_block
            imp_m = regexp(line, 'Impedance\s+in\s+Ohm\s+([+-]?[\dE+\-.]+)\s+([+-]?[\dE+\-.]+)', 'tokens');
            if ~isempty(imp_m) && cur_freq_idx > 0
                re_z = str2double(imp_m{1}{1});
                im_z = str2double(imp_m{1}{2});
                zin_data{end+1} = {cur_freq_idx, complex(re_z, im_z)};
                in_vsource_block = false;
            end
        end
    end

% --- Reference impedance (fields configuration). Not gated by all_ports_done:
% these lines come after 'READING SOLUTION COEFFICIENTS', which sets it. ---
    if in_fields_config
    if contains(line, 'Reference impedance in Ohm:') && isnan(Z0_ref)
        tok = regexp(line, 'Reference impedance in Ohm:\s*([\dE+\-\.]+)', 'tokens');
        if ~isempty(tok)
            Z0_ref = str2double(tok{1}{1});
        end
    end

    % Alternative format: S-PARAMETER REFERENCE IMPEDANCES AT PORTS
    % followed by lines: '   port     impedance in Ohm'
    % then '      1        1.88365E+02'
    if contains(line, 'S-PARAMETER REFERENCE IMPEDANCES AT PORTS')
        sp_ref_next = true;
    end
    if exist('sp_ref_next','var') && sp_ref_next && isnan(Z0_ref)
        tok = regexp(line, '^\s*1\s+([\dE+\-\.]+)', 'tokens');
        if ~isempty(tok)
            Z0_ref = str2double(tok{1}{1});
            sp_ref_next = false;
        end
    end
    end

    % --- Port collection (first frequency of the fields configuration only) ---
    if ~all_ports_done && in_fields_config

    % --- Load blocks (DATA FOR EDGE LOADS): each load's impedance, for the
    % cross-check against the S-parameter reference impedance. Must precede the
    % "~in_source_block -> continue" below. ---
    if contains(line, 'DATA FOR EDGE LOADS')
        in_load_section = true;
    elseif contains(line, 'READING SOLUTION COEFFICIENTS') || contains(line, 'DATA FOR THE GREEN''S FUNCTION')
        in_load_section = false;
    end

    if in_load_section
        if contains(line, 'Name:')
            tok = regexp(line, 'Name:\s*(.+)$', 'tokens');
            field_loads(end+1) = struct('name', '', 'index', NaN, 'circuit_type', '', ...
                'Z', NaN, 'attached_port', '', 'edge_indices', [], 'face_names', {{'',''}});
            if ~isempty(tok)
                field_loads(end).name = strtrim(tok{1}{1});
            end
        elseif ~isempty(field_loads) && contains(line, 'Load index:')
            tok = regexp(line, 'Load index:\s*(\d+)', 'tokens');
            if ~isempty(tok), field_loads(end).index = str2double(tok{1}{1}); end
        elseif ~isempty(field_loads) && contains(line, 'Circuit type:')
            tok = regexp(line, 'Circuit type:\s*(.+)$', 'tokens');
            if ~isempty(tok), field_loads(end).circuit_type = strtrim(tok{1}{1}); end
        elseif ~isempty(field_loads) && contains(line, 'Complex impedance:')
            tok = regexp(line, '\(\s*([+-]?[\dE+\-\.]+)\s*\+j\s*([+-]?[\dE+\-\.]+)\)', 'tokens');
            if ~isempty(tok)
                field_loads(end).Z = complex(str2double(tok{1}{1}), str2double(tok{1}{2}));
            end
        elseif ~isempty(field_loads) && contains(line, 'Attached to port:')
            tok = regexp(line, 'Attached to port:\s*(\S+)', 'tokens');
            if ~isempty(tok), field_loads(end).attached_port = tok{1}{1}; end
        elseif ~isempty(field_loads) && contains(line, '(negative faces)')
            load_face_pending = true;
        elseif ~isempty(field_loads) && load_face_pending
            load_face_pending = false;
            face_tok = regexp(strtrim(line), '(\S+)\s+(\S+)', 'tokens');
            if ~isempty(face_tok)
                field_loads(end).face_names = {face_tok{1}{1}, face_tok{1}{2}};
            end
        elseif ~isempty(field_loads) && contains(line, 'Indices of the edges:')
            tok = regexp(line, '(-?\d+)', 'tokens');
            eidx_list = [];
            for ti = 1:length(tok)
                eidx_list(end+1) = str2double(tok{ti}{1});
            end
            field_loads(end).edge_indices = eidx_list;
        end
    end

    % --- Detect source block start ---
    if contains(line, 'EXCITATION BY VOLTAGE SOURCE AT A SEGMENT')
        in_source_block = true;
        if isempty(port_type), port_type = 'wire'; end
        if isempty(excitation_type), excitation_type = 'Voltage Source at a Segment'; end
        current_R0 = NaN(1,3);
        cur_V0_mag = NaN; cur_V0_phase = NaN;
        continue
    elseif contains(line, 'EXCITATION BY VOLTAGE SOURCE AT A NODE')
        in_source_block = true;
        if isempty(port_type), port_type = 'wire'; end
        if isempty(excitation_type), excitation_type = 'Voltage Source at a Wire-Node'; end
        current_R0 = NaN(1,3);
        cur_V0_mag = NaN; cur_V0_phase = NaN;
        continue
    elseif contains(line, 'EXCITATION BY VOLTAGE SOURCE AT AN EDGE')
        in_source_block = true;
        if isempty(port_type), port_type = 'edge'; end
        if isempty(excitation_type), excitation_type = 'Voltage Source at an Edge'; end
        edge_indices{end+1} = [];      % new group for this port
        edge_face_names{end+1} = {'',''};  % placeholder {neg, pos}
        edge_port_names{end+1} = '';   % placeholder port name (e.g. 'EdgePort1')
        cur_V0_mag = NaN; cur_V0_phase = NaN;
        continue
    end

    if ~in_source_block, continue; end

    % --- Detect second pass via Source index: 1 repeating ---
    if contains(line, 'Source index:')
        tok = regexp(line, 'Source index:\s*(\d+)', 'tokens');
        if ~isempty(tok) && str2double(tok{1}{1}) == 1
            if seen_src_index1
                all_ports_done = true;
                in_source_block = false;
                continue
            else
                seen_src_index1 = true;
            end
        end
    end

    % --- Parse source fields ---
    if contains(line, '|U0| =')
        tok = regexp(line, '\|U0\|\s*=\s*([\dE+\-\.]+)', 'tokens');
        if ~isempty(tok), cur_V0_mag = str2double(tok{1}{1}); end

    elseif contains(line, 'ARG(U0) =')
        tok = regexp(line, 'ARG\(U0\)\s*=\s*([+-]?[\dE+\-\.]+)', 'tokens');
        if ~isempty(tok)
            cur_V0_phase = str2double(tok{1}{1});
            % ARG(U0) completes the source voltage: store it now
            if strcmp(port_type, 'wire') && ~isnan(cur_V0_mag)
                V0_mag_all(end+1)   = cur_V0_mag;
                V0_phase_all(end+1) = cur_V0_phase;
                if ~first_port_done, first_port_done = true; end
            end
            if strcmp(port_type, 'edge')
                V0_mag_all(end+1)   = cur_V0_mag;
                V0_phase_all(end+1) = cur_V0_phase;
                if ~first_port_done, first_port_done = true; end
            end
        end

    elseif strcmp(port_type, 'edge') && contains(line, 'Attached to port:')
        tok = regexp(line, 'Attached to port:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            edge_port_names{end} = tok{1}{1};
        end

    elseif strcmp(port_type, 'edge') && contains(line, 'Indices of the edges:')
        tok = regexp(line, '(-?\d+)', 'tokens');
        for ti = 1:length(tok)
            % Keep the sign: it selects the positive triangle (see mesh section).
            edge_indices{end}(end+1) = str2double(tok{ti}{1});
        end
        in_source_block = false;

    elseif strcmp(port_type, 'edge') && contains(line, '(negative faces)')
        face_line_pending = true;

    elseif strcmp(port_type, 'edge') && face_line_pending
        face_line_pending = false;
        face_tok = regexp(strtrim(line), '(\S+)\s+(\S+)', 'tokens');
        if ~isempty(face_tok)
            neg_name = face_tok{1}{1};
            pos_name = face_tok{1}{2};
            edge_face_names{end} = {neg_name, pos_name};
        end

    elseif strcmp(port_type, 'wire') && contains(line, 'Location of the port in m:')
        tok = regexp(line, 'x\s*=\s*([+-]?[\dE+\-\.]+),\s*y\s*=\s*([+-]?[\dE+\-\.]+),\s*z\s*=\s*([+-]?[\dE+\-\.]+)', 'tokens');
        if ~isempty(tok)
            current_R0 = [str2double(tok{1}{1}), str2double(tok{1}{2}), str2double(tok{1}{3})];
        end

    elseif strcmp(port_type, 'wire') && contains(line, 'Positive direction:')
        tok = regexp(line, 'x\s*=\s*([+-]?[\dE+\-\.]+),\s*y\s*=\s*([+-]?[\dE+\-\.]+),\s*z\s*=\s*([+-]?[\dE+\-\.]+)', 'tokens');
        if ~isempty(tok)
            this_dir = [str2double(tok{1}{1}), str2double(tok{1}{2}), str2double(tok{1}{3})];
            if size(port_pol,1) == 1 && all(port_pol == [0 0 1])
                port_pol = this_dir;   % first port: initialise
            else
                port_pol(end+1,:) = this_dir;
            end
        end
        if ~any(isnan(current_R0))
            port_locations(end+1,:) = current_R0;
        end
        in_source_block = false;
    end

    end % if ~all_ports_done

end % while

fclose(fid);

% --- Freq: the S-parameter sweep if present (it indexes S), else every
% frequency found while scanning the file ---
if ~isempty(Sfreq_internal)
    Freq = Sfreq_internal;
else
    Freq = sort(freq_list(:), 'ascend');
    if ~isempty(Freq)
        fprintf('ReadFEKOout_loadports: No S-parameter configuration found; Freq falls back to the %d frequency point(s) seen in the fields configuration "%s".\n', ...
            numel(Freq), fields_config_name);
    end
end

% --- Primary-port Zin (internal), sorted by frequency. check_reflection
% compares it with S by frequency value, not by position. ---
Zin = []; Zin_freq = [];
if ~isempty(zin_data)
    Zin_raw = NaN(length(freq_list), 1);
    for zi = 1:length(zin_data)
        Zin_raw(zin_data{zi}{1}) = zin_data{zi}{2};
    end
    [Zin_freq, freq_sort_idx] = sort(freq_list(:), 'ascend');
    Zin = Zin_raw(freq_sort_idx);
end

% --- Compute R0 as average of all port locations ---
if strcmp(port_type, 'wire') && ~isempty(port_locations)
    R0 = mean(port_locations, 1);
    fprintf('ReadFEKOout: %d %s port(s) found. R0 set to average location [%g, %g, %g] m.\n', ...
        size(port_locations,1), excitation_type, R0(1), R0(2), R0(3));
end

% --- PMC ground plane: wire port in-plane orientation ---
if has_pmc_plane
    if strcmp(port_type, 'wire') && abs(R0(3) - pmc_z) < 1e-6
        if abs(port_pol(1,1)) > 0.01 || abs(port_pol(1,2)) > 0.01
            fprintf('ReadFEKOout: Infinite PMC ground plane detected at z = %g m. Wire port is at ground plane with in-plane orientation. Applying V0 = V0*2.\n', pmc_z);
            pmc_double_V0 = true;
        else
            error('ReadFEKOout: Infinite PMC ground plane detected at z = %g m. Wire port orientation is not in-plane. Cannot determine correct V0 scaling.', pmc_z);
        end
    else
        warning('ReadFEKOout: Infinite PMC ground plane detected at z = %g m. V0 unchanged.', pmc_z);
        pmc_double_V0 = false;
    end
else
    pmc_double_V0 = false;
end

% --- Build complex voltage vector (Nports x 1) ---
if isempty(V0_mag_all)
    V0 = NaN;
else
    V0 = V0_mag_all(:) .* exp(1j * V0_phase_all(:) * pi/180);
end

if pmc_double_V0
    V0 = V0 * 2;
end

% --- Wire-port polarity correction for a separating PEC symmetry wall ---
% Wire ports separated by a PEC wall are electric images, so the tangential part
% of the positive direction flips across the wall. FEKO reports the same
% positive direction for all such ports (e.g. [0,0,-1] for ports at +x and -x
% across the YZ electric wall, where the image of port 1 is [0,0,+1]), so
% port_pol of ports 2..N is corrected here. V0 is not changed.
if strcmp(port_type, 'wire') && size(port_pol, 1) > 1
    % Find the electric wall(s) that separate the ports
    sep_elec_normals = zeros(0,3);
    if ~isempty(port_locations)
        port_spread = max(port_locations,[],1) - min(port_locations,[],1);
        for sp_i = 1:length(sym_planes)
            if ~strcmp(sym_planes{sp_i}.type, 'electric'), continue; end
            n = sym_planes{sp_i}.normal;
            if norm(n) < 0.5, continue; end
            if abs(dot(port_spread, n)) > 1e-6
                sep_elec_normals(end+1,:) = n / norm(n);
            end
        end
    end
    if ~isempty(sep_elec_normals)
        for pi_ = 2:size(port_pol, 1)
            p = port_pol(pi_,:);
            p_old = p;
            % Flip the tangential components (remove normal, flip, restore)
            for ni = 1:size(sep_elec_normals,1)
                n = sep_elec_normals(ni,:);
                p_normal = dot(p, n) * n;
                p_tangential = p - p_normal;
                p = p_normal - p_tangential;  % flip tangential
            end
            port_pol(pi_,:) = p;
            fprintf('ReadFEKOout: Wire port %d pol corrected for PEC image: [%g %g %g] -> [%g %g %g]\n', ...
                pi_, p_old(1), p_old(2), p_old(3), p(1), p(2), p(3));
        end
    end
end


% --- Parse mesh data: triangles + edges (for apothem and edge port locations) ---
antenna_apothem = NaN;
fid2 = fopen(filename, 'r');
if fid2 >= 0
    % --- Triangles (vertices, body label) and edges (KORP, KORM, POIP, POIM) ---
    in_triangles = false;
    in_edges     = false;
    triangles    = containers.Map('KeyType','int32','ValueType','any');
    tri_labels   = containers.Map('KeyType','int32','ValueType','char');  % tri_no -> face label
    edges_map    = containers.Map('KeyType','int32','ValueType','any');
    max_r2       = 0;
    found_mesh   = false;
    coord_pat    = '([+-]?[\d]+\.[\d]*E[+-][\d]+)\s+([+-]?[\d]+\.[\d]*E[+-][\d]+)\s+([+-]?[\d]+\.[\d]*E[+-][\d]+)';
    normal_pat   = '^[\s\d.E+\-]+$';
    cur_tri_no   = -1;
    cur_verts    = {};
    cur_label    = '';

    while true
        line = fgetl(fid2);
        if ~ischar(line), break; end
        line_raw = strtrim(line);

        if contains(line_raw, 'DATA OF THE METALLIC TRIANGLES')
            in_triangles = true; in_edges = false; continue
        end
        if contains(line_raw, 'DATA OF THE METALLIC EDGES')
            % Save last triangle
            if cur_tri_no > 0 && length(cur_verts) == 3
                triangles(int32(cur_tri_no)) = cur_verts;
                tri_labels(int32(cur_tri_no)) = cur_label;
            end
            in_triangles = false; in_edges = true; continue
        end
        if ~in_triangles && ~in_edges, continue; end

        if in_triangles
            % Skip normal/area lines (4 floats only)
            if ~isempty(regexp(line_raw, normal_pat, 'once')), continue; end
            tok = regexp(line_raw, coord_pat, 'tokens');
            if ~isempty(tok)
                % Check if this is a new triangle (has integer label before coords)
                tri_m = regexp(line_raw, '^(\d+)\s+(\S+)\s+', 'tokens');
                if ~isempty(tri_m)
                    % Save previous triangle
                    if cur_tri_no > 0 && length(cur_verts) == 3
                        triangles(int32(cur_tri_no)) = cur_verts;
                        tri_labels(int32(cur_tri_no)) = cur_label;
                    end
                    cur_tri_no = str2double(tri_m{1}{1});
                    cur_label  = strtrim(tri_m{1}{2});
                    cur_verts  = {[str2double(tok{1}{1}), str2double(tok{1}{2}), str2double(tok{1}{3})]};
                else
                    cur_verts{end+1} = [str2double(tok{1}{1}), str2double(tok{1}{2}), str2double(tok{1}{3})];
                end
                % Also compute apothem (relative to the origin [0 0 0], not R0)
                x = str2double(tok{1}{1});
                y = str2double(tok{1}{2});
                z = str2double(tok{1}{3});
                r2 = x^2 + y^2 + z^2;
                if r2 > max_r2, max_r2 = r2; found_mesh = true; end
            end
        end

        if in_edges
            % Edge table: no. type length Free s -1 KORP KORM POIP POIM ...
            edge_m = regexp(line_raw, '^(\d+)\s+\d+\s+[\dE+\-\.]+\s+\S+\s+\S+\s+-?\d+\s+(-?\d+)\s+(-?\d+)\s+(\d+)\s+(\d+)', 'tokens');
            if ~isempty(edge_m)
                eno  = int32(str2double(edge_m{1}{1}));
                korp = abs(str2double(edge_m{1}{2}));
                korm = abs(str2double(edge_m{1}{3}));
                poip = str2double(edge_m{1}{4});
                poim = str2double(edge_m{1}{5});
                edges_map(eno) = {korp, korm, poip, poim};
            end
        end
    end
    fclose(fid2);

    if found_mesh
        antenna_apothem = sqrt(max_r2);
        fprintf('ReadFEKOout: Antenna apothem (max |r - origin|) = %g m.\n', antenna_apothem);
    end

    % --- Compute edge port locations and polarity from mesh ---
    if strcmp(port_type, 'edge') && ~isempty(edge_indices) && triangles.Count > 0 && edges_map.Count > 0
        edge_midpoints = zeros(0,3);
        pol_vecs       = zeros(0,3);

        for pg = 1:length(edge_indices)
            group_mps = zeros(0,3);
            pg_pol    = zeros(0,3);

            for ei = 1:length(edge_indices{pg})
                eidx_signed = edge_indices{pg}(ei);
                eidx = int32(abs(eidx_signed));
                if ~isKey(edges_map, eidx), continue; end
                edata = edges_map(eidx);
                korp = int32(edata{1}); korm = int32(edata{2});
                poip = edata{3};        poim = edata{4};
                if ~isKey(triangles, korp) || ~isKey(triangles, korm), continue; end

                % The sign of the edge index selects the port's positive
                % triangle (KORP/KORM numbering is arbitrary, and the triangle
                % label does not distinguish sub-faces):
                %   negative index -> KORP (free vertex POIP) is positive
                %   positive index -> KORM (free vertex POIM) is positive
                % Checked on a wire-fed dipole (negative indices) and an
                % edge-fed bowtie (positive indices).
                if eidx_signed < 0
                    pos_tri = korp; pos_poi = poip;
                else
                    pos_tri = korm; pos_poi = poim;
                end
                verts_pos = triangles(pos_tri);
                if length(verts_pos) < 3 || pos_poi > length(verts_pos), continue; end

                % POI is the FREE vertex index of the positive triangle (1-indexed).
                % The shared edge = the two remaining vertices.
                % POI=1 -> shared = v2,v3; POI=2 -> shared = v1,v3; POI=3 -> shared = v1,v2
                all_vi = [1 2 3];
                shared_vi = all_vi(all_vi ~= pos_poi);   % two vertex indices of shared edge
                sv1 = verts_pos{shared_vi(1)};
                sv2 = verts_pos{shared_vi(2)};
                group_mps(end+1,:) = (sv1 + sv2) / 2;

                % Polarity: from the edge midpoint toward the positive
                % triangle's free vertex, perpendicular to the edge.
                edge_vec = sv2 - sv1;
                edge_len = norm(edge_vec);
                if edge_len < 1e-15, continue; end
                edge_hat = edge_vec / edge_len;
                free_pos_v = verts_pos{pos_poi};        % free vertex of positive face
                edge_mid   = (sv1 + sv2) / 2;
                vec_to_free = free_pos_v - edge_mid;
                delta_perp  = vec_to_free - dot(vec_to_free, edge_hat) * edge_hat;
                if norm(delta_perp) > 1e-15
                    pg_pol(end+1,:) = delta_perp / norm(delta_perp);
                end
            end

            if ~isempty(group_mps)
                edge_midpoints(end+1,:) = mean(group_mps, 1);
            end
            if ~isempty(pg_pol)
                mean_pol = mean(pg_pol, 1);
                pol_vecs(end+1,:) = mean_pol / norm(mean_pol);
            end
        end

        if ~isempty(edge_midpoints)
            port_locations = edge_midpoints;
            R0 = mean(edge_midpoints, 1);
            fprintf('ReadFEKOout: Edge port locations computed from mesh:\n');
            for ep_i = 1:size(edge_midpoints,1)
                fprintf('  Port %d: [%g, %g, %g] m\n', ep_i, edge_midpoints(ep_i,1), edge_midpoints(ep_i,2), edge_midpoints(ep_i,3));
            end
            fprintf('ReadFEKOout: R0 (mean of port locations) = [%g, %g, %g] m.\n', R0(1), R0(2), R0(3));
        else
            fprintf('ReadFEKOout: Edge port detected. Could not resolve location from mesh. R0 = [0 0 0].\n');
        end

        if ~isempty(pol_vecs)
            port_pol = zeros(size(pol_vecs, 1), 3);
            for pg = 1:size(pol_vecs, 1)
                port_pol(pg,:) = pol_vecs(pg,:) / norm(pol_vecs(pg,:));
            end

            % --- Resolve the polarity sign from face labels ---
            % pol points toward the sign-selected triangle. If that triangle's
            % body label differs from the body of the stated positive face,
            % flip port_pol. V0 is left as parsed.
            for pg = 1:length(edge_indices)
                if pg > length(edge_face_names), break; end
                pos_face_label = edge_face_names{pg}{2};   % stated positive face name
                if isempty(pos_face_label), continue; end
                % Body label of the sign-selected positive triangle
                korp_label = '';
                for ei = 1:length(edge_indices{pg})
                    eidx_signed = edge_indices{pg}(ei);
                    eidx = int32(abs(eidx_signed));
                    if ~isKey(edges_map, eidx), continue; end
                    edata = edges_map(eidx);
                    if eidx_signed < 0
                        pos_tri = int32(edata{1});   % KORP
                    else
                        pos_tri = int32(edata{2});   % KORM
                    end
                    if isKey(tri_labels, pos_tri)
                        korp_label = tri_labels(pos_tri);
                        break
                    end
                end
                if isempty(korp_label), continue; end
                % Positive face label 'Union1.Face1_1' -> body 'Union1'
                dot_idx = strfind(pos_face_label, '.');
                if ~isempty(dot_idx)
                    pos_body = pos_face_label(1:dot_idx(1)-1);
                else
                    pos_body = pos_face_label;
                end
                % The triangle table truncates body names to 6 characters
                % (e.g. 'right_ring' -> 'right_'); compare on that width.
                tri_label_width = 6;
                pos_body_cmp = pos_body(1:min(tri_label_width, length(pos_body)));
                if ~strcmp(pos_body_cmp, korp_label)
                    port_pol(pg,:) = -port_pol(pg,:);
                end
            end

            for pg = 1:size(port_pol, 1)
                fprintf('ReadFEKOout: Edge port %d polarity direction: [%g, %g, %g],  V0 = %.4g∠%.1f°\n', ...
                    pg, port_pol(pg,1), port_pol(pg,2), port_pol(pg,3), abs(V0(pg)), angle(V0(pg))*180/pi);
            end
        else
            fprintf('ReadFEKOout: Could not determine edge port polarity direction from mesh.\n');
            port_pol = repmat([0 0 1], length(edge_indices), 1);   % fallback to default per port
        end
    elseif strcmp(port_type, 'edge')
        fprintf('ReadFEKOout: Edge port detected. No mesh data found. R0 set to [0 0 0]. Please verify in CADFEKO.\n');
    end
end

% --- No mesh in the .out file (e.g. some patch models): ask for the apothem
% rather than propagate NaN (it sets the exp(-jka) reference plane downstream) ---
if isnan(antenna_apothem)
    fprintf('ReadFEKOout_loadports: WARNING - antenna apothem could not be determined from the mesh data in this .out file (no "DATA OF THE METALLIC TRIANGLES" section found, or it contained no usable vertex data).\n');
    antenna_apothem = [];
    while isempty(antenna_apothem) || ~isnumeric(antenna_apothem) || ~isscalar(antenna_apothem) || ~isreal(antenna_apothem) || antenna_apothem <= 0
        antenna_apothem = input('Enter antenna apothem (meters): ');
    end
end

% --- Reflection-coefficient consistency check ---
% Gamma = (Zin-Z0)/(Zin+Z0) from FEKO's source impedance must match S(i,i,:):
% |Gamma - S| <= 1% of max(|S|, |Gamma|, 0.01) at every frequency. Used for the
% primary port and for each secondary port with its own fields configuration.
freq_tol = 1e-6;   % relative tolerance for matching two frequency vectors
    function ok = freqs_match(fa, fb)
        ok = (numel(fa) == numel(fb)) && all(abs(fa(:) - fb(:)) <= freq_tol * max(abs(fb(:)), 1));
    end

refl_tol = 0.01;   % 1% tolerance, applied to the full complex difference
    function match = check_reflection(port_label, sp_idx, Zin_v, Zin_freq_v, Z0_v)
        match = NaN;
        if isempty(S) || isempty(sp_idx) || isnan(Z0_v) || isempty(Zin_v)
            return
        end
        if ~freqs_match(Zin_freq_v, Freq)
            fprintf('ReadFEKOout_loadports: WARNING - cannot check reflection coefficient consistency for port "%s": its frequency sweep does not match Freq.\n', port_label);
            return
        end
        Gamma = (Zin_v(:) - Z0_v) ./ (Zin_v(:) + Z0_v);
        Sdiag = squeeze(S(sp_idx, sp_idx, :));
        Sdiag = Sdiag(:);
        err = abs(Gamma - Sdiag);   % complex distance: both real and imaginary parts together
        tol = refl_tol * max([abs(Sdiag), abs(Gamma), 0.01*ones(size(Sdiag))], [], 2);
        bad = err > tol;
        match = ~any(bad);
        if ~match
            [~, worst_i] = max(err ./ max(tol, eps));
            fprintf('ReadFEKOout_loadports: WARNING - port "%s" reflection coefficient disagrees between S(%d,%d,:) and (Zin-Z0)/(Zin+Z0) at %d/%d frequency point(s) (worst at %.6g Hz: S=%.4g%+.4gj, Gamma=%.4g%+.4gj).\n', ...
                port_label, sp_idx, sp_idx, sum(bad), numel(bad), Freq(worst_i), ...
                real(Sdiag(worst_i)), imag(Sdiag(worst_i)), real(Gamma(worst_i)), imag(Gamma(worst_i)));
        else
            fprintf('ReadFEKOout_loadports: Port "%s" reflection coefficient (Zin-Z0)/(Zin+Z0) agrees with S(%d,%d,:) within 1%% at all %d frequency point(s).\n', ...
                port_label, sp_idx, sp_idx, numel(Sdiag));
        end
    end

% --- Locate secondary port(s) from the S-parameter configuration and
% cross-check them against the load(s) defined in the fields configuration ---
primary_sp_idx = [];    % this port's index within sparam_ports/S, if resolved
primary_port_name = '';
if strcmp(port_type, 'edge') && ~isempty(edge_port_names)
    primary_port_name = edge_port_names{1};
end
if ~isempty(sparam_config_name)
    [sparam_ports, sparam_ref_Z] = parse_sparam_config(filename, sparam_config_name);

    % --- Primary port: same edge set in the fields and S-parameter
    % configurations? ---
    if ~isempty(primary_port_name)
        primary_sp_idx = find(strcmp({sparam_ports.name}, primary_port_name), 1);
        if isempty(primary_sp_idx)
            sparam_diagnostic = sprintf(['A multi-port S-parameter configuration ("%s") was found, but the ' ...
                'primary port ("%s", from the fields configuration "%s") does not appear in its port list ' ...
                '[%s]. The primary port''s own reflection coefficient cannot be read from it. This usually ' ...
                'means the port was renamed, remeshed, or is not actually included in the S-parameter sweep.'], ...
                sparam_config_name, primary_port_name, fields_config_name, strjoin({sparam_ports.name}, ', '));
            fprintf('ReadFEKOout_loadports: WARNING - %s\n', sparam_diagnostic);
        else
            own_edges = [];
            if ~isempty(edge_indices), own_edges = edge_indices{1}; end
            sp_edges = sparam_ports(primary_sp_idx).edge_indices;
            primary_edges_match = isequal(sort(abs(own_edges)), sort(abs(sp_edges)));
            if primary_edges_match
                fprintf('ReadFEKOout_loadports: Primary port "%s" definition (edges: %s) agrees between "%s" and "%s".\n', ...
                    primary_port_name, mat2str(sort(abs(own_edges))), fields_config_name, sparam_config_name);
            else
                fprintf('ReadFEKOout_loadports: WARNING - primary port "%s" edge definition differs between "%s" (edges: %s) and "%s" (edges: %s). Port may not be the same across configurations.\n', ...
                    primary_port_name, fields_config_name, mat2str(sort(abs(own_edges))), ...
                    sparam_config_name, mat2str(sort(abs(sp_edges))));
            end

            % --- Primary port: same reference impedance in both configurations? ---
            if isKey(sparam_ref_Z, int32(primary_sp_idx))
                primary_sp_Z0 = sparam_ref_Z(int32(primary_sp_idx));
                if ~isnan(Z0_ref)
                    z0_tol = 1e-3 * max(abs(primary_sp_Z0), 1);
                    if abs(Z0_ref - primary_sp_Z0) < z0_tol
                        fprintf('ReadFEKOout_loadports: Primary port "%s" reference impedance (%.6g Ohm) agrees between "%s" and "%s".\n', ...
                            primary_port_name, Z0_ref, fields_config_name, sparam_config_name);
                    else
                        fprintf('ReadFEKOout_loadports: WARNING - primary port "%s" reference impedance differs: %.6g Ohm in "%s" vs %.6g Ohm in "%s".\n', ...
                            primary_port_name, Z0_ref, fields_config_name, primary_sp_Z0, sparam_config_name);
                    end
                end
            end
        end
    end

    for sp_i = 1:numel(sparam_ports)
        this_name = sparam_ports(sp_i).name;
        if isempty(this_name) || strcmp(this_name, primary_port_name)
            continue   % this is the primary/fields-driven port itself, skip
        end

        entry = struct('name', this_name, 'edge_indices', sparam_ports(sp_i).edge_indices, ...
            'face_names', {sparam_ports(sp_i).face_names}, 'R2', [NaN NaN NaN], ...
            'pol2', [NaN NaN NaN], 'Z_sparam_ref', NaN, 'Z_field_load', NaN, ...
            'edges_match', false, 'impedance_match', false, 'has_own_field_config', false, ...
            'field_config_name', '', 'V0_own', NaN, 'Z0_ref_own', NaN, ...
            'geometry_match', NaN, 'edges_match_field_config', NaN, ...
            'impedance_match_own', NaN, 'near_field_name', '', 'near_field_occurrence', 0, 'near_field_type', '', ...
            'far_field_name', '', 'far_field_occurrence', 0, 'far_field_type', '', 'sparam_port_index', NaN, 'reflection_match', NaN);

        entry.sparam_port_index = sp_i;

        % Reference impedance for this port, from the S-PARAMETER REFERENCE
        % IMPEDANCES AT PORTS table (indexed by port/source order).
        if isKey(sparam_ref_Z, int32(sp_i))
            entry.Z_sparam_ref = sparam_ref_Z(int32(sp_i));
        end

        % Geometry (location + polarity) from the mesh
        if exist('triangles','var') && exist('edges_map','var') && triangles.Count > 0 && edges_map.Count > 0
            [R2, pol2] = compute_edge_group_geometry(entry.edge_indices, triangles, edges_map);
            entry.R2   = R2;
            entry.pol2 = pol2;

            % --- Resolve the polarity sign from face labels (as for the
            % primary port, 6-character body-name width) ---
            if exist('tri_labels','var') && ~isempty(entry.face_names) && ~isempty(entry.face_names{2})
                pos_face_label = entry.face_names{2};
                korp_label = '';
                for ei = 1:length(entry.edge_indices)
                    eidx_signed = entry.edge_indices(ei);
                    eidx = int32(abs(eidx_signed));
                    if ~isKey(edges_map, eidx), continue; end
                    edata = edges_map(eidx);
                    if eidx_signed < 0
                        pos_tri = int32(edata{1});   % KORP
                    else
                        pos_tri = int32(edata{2});   % KORM
                    end
                    if isKey(tri_labels, pos_tri)
                        korp_label = tri_labels(pos_tri);
                        break
                    end
                end
                if ~isempty(korp_label)
                    dot_idx = strfind(pos_face_label, '.');
                    if ~isempty(dot_idx)
                        pos_body = pos_face_label(1:dot_idx(1)-1);
                    else
                        pos_body = pos_face_label;
                    end
                    tri_label_width = 6;
                    pos_body_cmp = pos_body(1:min(tri_label_width, length(pos_body)));
                    if ~strcmp(pos_body_cmp, korp_label) && all(~isnan(entry.pol2))
                        entry.pol2 = -entry.pol2;
                    end
                end
            end
        end

        % Matching load in the primary fields configuration: by attached port
        % name, else by unsigned edge-index set.
        match_idx = [];
        if ~isempty(field_loads)
            match_idx = find(strcmp({field_loads.attached_port}, this_name), 1);
        end
        if isempty(match_idx) && ~isempty(entry.edge_indices)
            want = sort(abs(entry.edge_indices));
            for fl_i = 1:numel(field_loads)
                if isequal(sort(abs(field_loads(fl_i).edge_indices)), want)
                    match_idx = fl_i;
                    break
                end
            end
        end
        if ~isempty(match_idx)
            entry.Z_field_load = field_loads(match_idx).Z;
            entry.edges_match = isequal(sort(abs(entry.edge_indices)), ...
                sort(abs(field_loads(match_idx).edge_indices)));
            if ~isnan(entry.Z_sparam_ref) && ~isnan(entry.Z_field_load)
                tol = 1e-3 * max(abs(entry.Z_sparam_ref), 1);
                entry.impedance_match = abs(entry.Z_field_load - entry.Z_sparam_ref) < tol;
            end
        end

        fprintf('ReadFEKOout_loadports: Secondary port "%s" (edges: %s):\n', this_name, mat2str(entry.edge_indices));
        fprintf('  Location from mesh: [%g, %g, %g] m,  polarity: [%g, %g, %g]\n', ...
            entry.R2(1), entry.R2(2), entry.R2(3), entry.pol2(1), entry.pol2(2), entry.pol2(3));
        if isempty(match_idx)
            fprintf('  WARNING: no matching load found for this port in the fields configuration ("%s").\n', fields_config_name);
        else
            fprintf('  S-parameter reference impedance: %.4g Ohm.  Fields-config load impedance: %.4g%+.4gj Ohm.\n', ...
                entry.Z_sparam_ref, real(entry.Z_field_load), imag(entry.Z_field_load));
            if entry.edges_match
                fprintf('  Edge indices agree with the fields-config load -- same physical port.\n');
            else
                fprintf('  WARNING: edge indices do NOT match the fields-config load edge indices.\n');
            end
            if entry.impedance_match
                fprintf('  Load impedance matches the S-parameter reference impedance.\n');
            else
                fprintf('  WARNING: load impedance does NOT match the S-parameter reference impedance.\n');
            end
        end

        secondary_ports(end+1) = entry; %#ok<AGROW>
    end
end

% --- Parse the other single-port fields configurations ---
% Each gives a directly solved V0, Zin, location and polarity for its port.
% Merged into the matching secondary_ports entry (by port name), or added as a
% new entry if the port is not in an S-parameter configuration.
if numel(fields_config_names) > 1 && exist('triangles','var') && exist('edges_map','var') && exist('tri_labels','var')
    for fc_i = 2:numel(fields_config_names)
        this_field_config = fields_config_names{fc_i};
        pinfo = parse_one_field_config(filename, this_field_config, triangles, edges_map, tri_labels);
        if isempty(pinfo.port_name)
            fprintf('ReadFEKOout_loadports: WARNING - could not identify the driven port in configuration "%s". Skipping.\n', this_field_config);
            continue
        end

        match_idx = find(strcmp({secondary_ports.name}, pinfo.port_name), 1);
        if isempty(match_idx)
            % Port not in an S-parameter configuration: new entry
            new_entry = struct('name', pinfo.port_name, 'edge_indices', pinfo.edge_indices, ...
                'face_names', {pinfo.face_names}, 'R2', [NaN NaN NaN], 'pol2', [NaN NaN NaN], ...
                'Z_sparam_ref', NaN, 'Z_field_load', NaN, 'edges_match', false, ...
                'impedance_match', false, 'has_own_field_config', false, 'field_config_name', '', ...
                'V0_own', NaN, 'Z0_ref_own', NaN, 'geometry_match', NaN, ...
                'edges_match_field_config', NaN, 'impedance_match_own', NaN, ...
                'near_field_name', '', 'near_field_occurrence', 0, 'near_field_type', '', 'far_field_name', '', ...
                'far_field_occurrence', 0, 'far_field_type', '', 'sparam_port_index', NaN, ...
                'reflection_match', NaN);
            secondary_ports(end+1) = new_entry; %#ok<AGROW>
            match_idx = numel(secondary_ports);
        end

        secondary_ports(match_idx).has_own_field_config = true;
        secondary_ports(match_idx).field_config_name    = this_field_config;
        secondary_ports(match_idx).near_field_name       = near_name_raw(this_field_config);
        secondary_ports(match_idx).far_field_name        = far_name_raw(this_field_config);
        secondary_ports(match_idx).near_field_occurrence = near_name_occ(this_field_config);
        secondary_ports(match_idx).far_field_occurrence  = far_name_occ(this_field_config);
        secondary_ports(match_idx).near_field_type       = near_name_type(this_field_config);
        secondary_ports(match_idx).far_field_type        = far_name_type(this_field_config);
        secondary_ports(match_idx).V0_own    = pinfo.V0;
        secondary_ports(match_idx).Z0_ref_own = pinfo.Z0_ref;

        % Reflection-coefficient check for this port (S(i,i,:) vs own Zin)
        if ~isempty(sparam_config_name) && exist('sparam_ports','var')
            sidx = find(strcmp({sparam_ports.name}, pinfo.port_name), 1);
            if ~isempty(sidx)
                secondary_ports(match_idx).reflection_match = check_reflection( ...
                    pinfo.port_name, sidx, pinfo.Zin, pinfo.Freq, pinfo.Z0_ref);
            end
        end

        fprintf('ReadFEKOout_loadports: Port "%s" also has its own dedicated fields configuration "%s":\n', pinfo.port_name, this_field_config);
        if ~isempty(secondary_ports(match_idx).near_field_name) || ~isempty(secondary_ports(match_idx).far_field_name)
            fprintf('  Field request name(s): near-field "%s", far-field/spherical-mode "%s". Expected files: <prefix>_%s.efe/.hfe, <prefix>_%s.ffe/.sph.\n', ...
                secondary_ports(match_idx).near_field_name, secondary_ports(match_idx).far_field_name, ...
                secondary_ports(match_idx).near_field_name, secondary_ports(match_idx).far_field_name);
        end

        % --- Reference impedance: Z_sparam_ref, Z_field_load and Z0_ref_own
        % describe the same termination and should agree. ---
        z0_vals = [secondary_ports(match_idx).Z_sparam_ref, secondary_ports(match_idx).Z_field_load, pinfo.Z0_ref];
        z0_labels = {sparam_config_name, 'primary fields-config load', this_field_config};
        z0_present = ~isnan(z0_vals);
        if sum(z0_present) >= 2
            present_vals = z0_vals(z0_present);
            present_labels = z0_labels(z0_present);
            z0_tol = 1e-3 * max(max(abs(present_vals)), 1);
            secondary_ports(match_idx).impedance_match_own = (max(present_vals) - min(present_vals)) < z0_tol;
            if secondary_ports(match_idx).impedance_match_own
                fprintf('  Reference/load impedance (%.6g Ohm) agrees across: %s.\n', ...
                    present_vals(1), strjoin(present_labels, ', '));
            else
                detail = cell(1, numel(present_vals));
                for zi = 1:numel(present_vals)
                    detail{zi} = sprintf('%s = %.6g Ohm', present_labels{zi}, present_vals(zi));
                end
                fprintf('  WARNING: reference/load impedance for port "%s" disagrees across configurations: %s.\n', ...
                    pinfo.port_name, strjoin(detail, '; '));
            end
        end

        % --- Edge set: S-parameter configuration vs this port's own fields
        % configuration ---
        if ~isempty(secondary_ports(match_idx).edge_indices) && ~isempty(pinfo.edge_indices)
            secondary_ports(match_idx).edges_match_field_config = isequal( ...
                sort(abs(secondary_ports(match_idx).edge_indices)), sort(abs(pinfo.edge_indices)));
            if secondary_ports(match_idx).edges_match_field_config
                fprintf('  Edge definition (edges: %s) agrees between "%s" and "%s".\n', ...
                    mat2str(sort(abs(pinfo.edge_indices))), sparam_config_name, this_field_config);
            else
                fprintf('  WARNING: edge definition differs between "%s" (edges: %s) and "%s" (edges: %s).\n', ...
                    sparam_config_name, mat2str(sort(abs(secondary_ports(match_idx).edge_indices))), ...
                    this_field_config, mat2str(sort(abs(pinfo.edge_indices))));
            end
        end

        % Compare with the geometry obtained so far, then use the own fields
        % configuration's geometry (as for the primary port's R0/port_pol).
        had_prior_geom = all(~isnan(secondary_ports(match_idx).R2));
        if had_prior_geom && all(~isnan(pinfo.R0))
            loc_err = norm(pinfo.R0 - secondary_ports(match_idx).R2);
            pol_err = norm(pinfo.pol - secondary_ports(match_idx).pol2);
            secondary_ports(match_idx).geometry_match = (loc_err < 1e-6) && (pol_err < 1e-6);
            if ~secondary_ports(match_idx).geometry_match
                fprintf('  WARNING: geometry from "%s" disagrees with the S-parameter-derived geometry (location err %.3g m, polarity err %.3g).\n', ...
                    this_field_config, loc_err, pol_err);
            end
        end
        secondary_ports(match_idx).R2   = pinfo.R0;
        secondary_ports(match_idx).pol2 = pinfo.pol;

        fprintf('  Location: [%g, %g, %g] m,  polarity: [%g, %g, %g]\n', ...
            pinfo.R0(1), pinfo.R0(2), pinfo.R0(3), pinfo.pol(1), pinfo.pol(2), pinfo.pol(3));
        if isnan(pinfo.V0)
            fprintf('  WARNING: could not parse V0 for this port from "%s".\n', this_field_config);
        else
            fprintf('  V0 = %.4g<%.1f deg,  %d Zin point(s) parsed (%.4g to %.4g Hz).\n', ...
                abs(pinfo.V0), angle(pinfo.V0)*180/pi, numel(pinfo.Zin), ...
                min([pinfo.Freq; NaN]), max([pinfo.Freq; NaN]));
        end
    end
end

% --- Build symmetry_planes struct ---
% symmetry_planes.XY / .XZ / .YZ: +1 = magnetic wall, 0 = none, -1 = electric wall
% symmetry_planes.physical_gnd:   1 = infinite ground plane coincident with port, 0 = none
symmetry_planes = struct('XY', 0, 'XZ', 0, 'YZ', 0, 'physical_gnd', 0);
for sp_i = 1:length(sym_planes)
    wval = 0;
    if strcmp(sym_planes{sp_i}.type, 'magnetic'), wval = +1;
    elseif strcmp(sym_planes{sp_i}.type, 'electric'), wval = -1;
    end
    if contains(sym_planes{sp_i}.plane, 'XY'), symmetry_planes.XY = wval;
    elseif contains(sym_planes{sp_i}.plane, 'XZ'), symmetry_planes.XZ = wval;
    elseif contains(sym_planes{sp_i}.plane, 'YZ'), symmetry_planes.YZ = wval;
    end
end

if has_pec_plane
    if abs(R0(3) - gnd_z) < 1e-6
        symmetry_planes.physical_gnd = 1;
        fprintf('ReadFEKOout: Infinite PEC ground plane at z = %g m. Port at z = %g m coincides. Applying V0 = V0*2.\n', gnd_z, R0(3));
        V0 = V0 * 2;
    else
        fprintf('ReadFEKOout: Infinite PEC ground plane at z = %g m. Port z = %g m does not coincide. V0 unchanged.\n', gnd_z, R0(3));
    end
end
if pmc_double_V0
    symmetry_planes.physical_gnd = 1;
end

% --- Set even_odd_excitation: symmetry wall + multiple ports ---
% Print all detected symmetry planes
if ~isempty(sym_planes)
    fprintf('ReadFEKOout: Symmetry planes detected:\n');
    for sp_i = 1:length(sym_planes)
        fprintf('  %s: %s conducting wall\n', sym_planes{sp_i}.plane, sym_planes{sp_i}.type);
    end
end

% A wall separates the ports if the port spread has a component along its normal.
num_ports = size(port_locations, 1);
n_sym_walls = n_mag_walls + n_elec_walls;
if n_sym_walls > 0 && num_ports > 1
    port_spread = max(port_locations,[],1) - min(port_locations,[],1);
    if max(port_spread) > 1e-6
        % Find which symmetry plane(s) separate the ports:
        % a plane separates the ports if its normal aligns with the port spread direction
        sep_wall_type = '';
        sep_plane_names = {};
        for sp_i = 1:length(sym_planes)
            n = sym_planes{sp_i}.normal;
            if norm(n) < 0.5, continue; end  % unknown plane, skip
            % Project port spread onto this plane's normal
            spread_along_normal = abs(dot(port_spread, n));
            if spread_along_normal > 1e-6
                sep_plane_names{end+1} = sym_planes{sp_i}.plane;
                if isempty(sep_wall_type)
                    sep_wall_type = sym_planes{sp_i}.type;
                elseif ~strcmp(sep_wall_type, sym_planes{sp_i}.type)
                    sep_wall_type = 'mixed';
                end
            end
        end
        sep_planes_str = strjoin(sep_plane_names, ', ');
        if strcmp(sep_wall_type, 'magnetic')
            even_odd_excitation = log2(num_ports);
            fprintf('ReadFEKOout: %d port(s) mirrored across %s (magnetic) -> even (common) mode. even_odd_excitation = +%g.\n', num_ports, sep_planes_str, even_odd_excitation);
        elseif strcmp(sep_wall_type, 'electric')
            even_odd_excitation = -log2(num_ports);
            fprintf('ReadFEKOout: %d port(s) mirrored across %s (electric) -> odd (differential) mode. even_odd_excitation = %g.\n', num_ports, sep_planes_str, even_odd_excitation);
        elseif strcmp(sep_wall_type, 'mixed')
            even_odd_excitation = log2(num_ports);
            fprintf('ReadFEKOout: %d port(s) separated by mixed walls (%s). even_odd_excitation = %g.\n', num_ports, sep_planes_str, even_odd_excitation);
        else
            fprintf('ReadFEKOout: %d port(s) spatially separated but no matching symmetry wall found.\n', num_ports);
        end
    else
        fprintf('ReadFEKOout: %d port(s) detected but not spatially separated - no common/differential mode assignment.\n', num_ports);
    end
elseif n_sym_walls > 0 && num_ports <= 1
    fprintf('ReadFEKOout: Symmetry planes present but single port - geometric symmetry only, no mode assignment.\n');
end

% --- Symmetry consistency check: V0 * port_pol vs separating wall type ---
% Computes the effective excitation vector per port and checks whether pairs
% of ports satisfy the required symmetry across each separating wall.
% Electric (PEC) wall: tangential components must be antisymmetric (t1 + t2 ~ 0)
% Magnetic (PMC) wall: tangential components must be symmetric   (t1 - t2 ~ 0)
if size(port_pol,1) > 1 && ~isempty(V0) && ~any(isnan(V0(:)))
    locs = port_locations;
    if isempty(locs) && ~isempty(edge_midpoints)
        locs = edge_midpoints;
    end
    eff = real(V0(:)) .* port_pol;   % effective excitation direction per port
    any_checked = false;
    for sp_i = 1:length(sym_planes)
        n     = sym_planes{sp_i}.normal;
        wtype = sym_planes{sp_i}.type;
        if norm(n) < 0.5, continue; end
        n = n / norm(n);
        for p1 = 1:size(port_pol,1)
            for p2 = p1+1:size(port_pol,1)
                if ~isempty(locs) && size(locs,1) >= max(p1,p2)
                    sep = locs(p2,:) - locs(p1,:);
                    if abs(dot(sep,n)) < 1e-9 * norm(sep), continue; end
                end
                any_checked = true;
                t1 = eff(p1,:) - dot(eff(p1,:),n)*n;
                t2 = eff(p2,:) - dot(eff(p2,:),n)*n;
                if strcmp(wtype,'electric')
                    residual = norm(t1 + t2);
                    sym_str  = 'antisymmetric';
                else
                    residual = norm(t1 - t2);
                    sym_str  = 'symmetric';
                end
                tol = 1e-4 * max(norm(t1), norm(t2));
                if residual <= tol || (norm(t1) < 1e-12 && norm(t2) < 1e-12)
                    fprintf('ReadFEKOout: Port excitation+polarity agrees with %s (%s) symmetry condition.\n', ...
                        sym_planes{sp_i}.plane, wtype);
                else
                    fprintf('ReadFEKOout: WARNING — port excitation+polarity disagrees with %s (%s) symmetry condition (residual = %.2e).\n', ...
                        sym_planes{sp_i}.plane, wtype, residual);
                end
            end
        end
    end
end

% --- Scale reference impedance based on ground plane type ---
if symmetry_planes.physical_gnd
    if has_pec_plane
        Z0_ref = Z0_ref * 2;
        fprintf('ReadFEKOout: PEC ground plane coincides with port. Applying Z0_ref = Z0_ref*2 = %g Ohm.\n', Z0_ref);
    elseif pmc_double_V0
        Z0_ref = Z0_ref / 2;
        fprintf('ReadFEKOout: PMC ground plane coincides with port. Applying Z0_ref = Z0_ref/2 = %g Ohm.\n', Z0_ref);
    end
end

% --- Reflection-coefficient consistency check (primary port) ---
if ~isempty(sparam_config_name) && exist('sparam_ports','var')
    if ~isempty(primary_sp_idx) && ~isempty(Zin)
        check_reflection(primary_port_name, primary_sp_idx, Zin, Zin_freq, Z0_ref);
    end
end

primary_sparam_index = primary_sp_idx;

fprintf('\n')

end


% =========================================================================
function [cfg_info, cfg_order] = classify_configurations(filename)
% classify_configurations  Pre-scan a FEKO .out file and summarise each
% uniquely named CONFIGURATION block from its first occurrence.
%
% Returns:
%   cfg_info  - containers.Map, configuration name -> struct:
%     .n_excited_ports - number of distinct ports after an "EXCITATION BY
%                        VOLTAGE SOURCE" header ("Attached to port:" lines in
%                        load blocks are not counted)
%     .port_names      - those port names, in order of appearance
%     .has_field_calc  - a far-field, near-field or spherical-mode request exists
%     .near_field_name - near-field request name ('' if none)
%     .far_field_name  - far-field or spherical-mode request name ('' if none)
%     .has_sparam_data - an "S sink source re im ..." row was printed
%     .freqs           - distinct frequencies (Hz) of this configuration over
%                        the WHOLE file, in order of appearance
%   cfg_order - configuration names in order of first appearance
cfg_info = containers.Map('KeyType','char','ValueType','any');
fid = fopen(filename, 'r');
if fid < 0
    cfg_order = {};
    return
end

cur_config      = '';
configs_done    = {};   % names whose first occurrence has been fully collected, in that order
in_excitation   = false;
in_load         = false;
cur_port_names  = {};
cur_has_field   = false;
cur_near_name   = '';
cur_far_name    = '';
cur_has_sparam  = false;   % config printed an "S sink source re im" row
freq_map        = containers.Map('KeyType','char','ValueType','any');   % config name -> distinct frequencies (Hz), whole file

    function flush_current()
        if ~isempty(cur_config) && ~any(strcmp(configs_done, cur_config))
            uniq_names = unique(cur_port_names, 'stable');
            cfg_info(cur_config) = struct('n_excited_ports', numel(uniq_names), ...
                'port_names', {uniq_names}, 'has_field_calc', cur_has_field, ...
                'near_field_name', cur_near_name, 'far_field_name', cur_far_name, ...
                'has_sparam_data', cur_has_sparam, 'freqs', []);
            configs_done{end+1} = cur_config;
        end
    end

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);

    if contains(line, 'CONFIGURATION:')
        tok = regexp(line, 'CONFIGURATION:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            new_config = tok{1}{1};
            if ~strcmp(new_config, cur_config)
                flush_current();
                cur_config = new_config;
                cur_port_names = {};
                cur_has_field  = false;
                cur_near_name  = '';
                cur_far_name   = '';
                cur_has_sparam = false;
                in_excitation  = false;
                in_load        = false;
            end
        end
        continue
    end

    % --- Frequency sweep of every occurrence of this config (needed to
    % classify field-request name reuse) ---
    freq_tok = regexp(line, 'FREQ\s*=\s*([\dE+\-.]+)', 'tokens');
    if ~isempty(freq_tok) && ~isempty(cur_config)
        fval = str2double(freq_tok{1}{1});
        if isKey(freq_map, cur_config)
            flist = freq_map(cur_config);
            if isempty(flist) || abs(flist(end) - fval) > abs(fval)*1e-6
                flist(end+1) = fval; %#ok<AGROW>
                freq_map(cur_config) = flist;
            end
        else
            freq_map(cur_config) = fval;
        end
    end

    if any(strcmp(configs_done, cur_config)), continue; end  % already collected this config

    % Presence of an S-matrix row (values are read by parse_sparam_matrix)
    if ~cur_has_sparam && ~isempty(regexp(line, '^S\s+\d+\s+\d+\s+[+-]?[\dE+\-\.]+\s+[+-]?[\dE+\-\.]+', 'once'))
        cur_has_sparam = true;
    end

    if contains(line, 'EXCITATION BY VOLTAGE SOURCE')
        in_excitation = true; in_load = false; continue
    end
    if contains(line, 'DATA FOR EDGE LOADS') || contains(line, 'DATA FOR LOADS') || contains(line, 'DATA FOR VERTEX LOADS')
        in_excitation = false; in_load = true; continue
    end
    if contains(line, 'READING SOLUTION COEFFICIENTS')
        in_excitation = false; in_load = false;
    end
    if contains(line, 'Far field request') || contains(line, 'Near field request') || contains(line, 'Spherical mode request')
        cur_has_field = true;
    end
    % Field-request names (FEKO exports <prefix>_<name>.<ext>); keep the first
    % of each. Far-field and spherical-mode requests usually share one name.
    if isempty(cur_near_name)
        tok = regexp(line, 'Near field request with name:\s*(\S+)', 'tokens');
        if ~isempty(tok), cur_near_name = tok{1}{1}; end
    end
    if isempty(cur_far_name)
        tok = regexp(line, 'Far field request with name:\s*(\S+)', 'tokens');
        if ~isempty(tok), cur_far_name = tok{1}{1}; end
        tok = regexp(line, 'Spherical mode request with name:\s*(\S+)', 'tokens');
        if ~isempty(tok), cur_far_name = tok{1}{1}; end
    end
    if in_excitation && contains(line, 'Attached to port:')
        tok = regexp(line, 'Attached to port:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            cur_port_names{end+1} = tok{1}{1};
        end
    end
end
flush_current();
fclose(fid);
cfg_order = configs_done;

% Attach each configuration's whole-file frequency sweep
for ci = 1:numel(cfg_order)
    cname = cfg_order{ci};
    info = cfg_info(cname);
    if isKey(freq_map, cname)
        info.freqs = freq_map(cname);
    else
        info.freqs = [];
    end
    cfg_info(cname) = info;
end
end


% =========================================================================
function [S, Sfreq] = parse_sparam_matrix(filename, sparam_config_name, N)
% parse_sparam_matrix  Read the complex S-matrix from every frequency repeat of
% the named S-parameter configuration.
%
% Within one repeat FEKO excites ports 1..N in turn and prints, for each
% source, one "S  sink  source  re  im  mag  dB  phase" line per sink.
%
% Inputs:
%   filename           - path to the .out file
%   sparam_config_name - S-parameter configuration name ('' if none)
%   N                  - number of ports
% Returns:
%   S     - (N x N x F) complex, S(sink,source,f), ascending frequency;
%           entries not reported are NaN+NaNi
%   Sfreq - (F x 1) frequencies (Hz), ascending
% Both are empty if sparam_config_name is empty or N == 0.
S = [];
Sfreq = [];
if isempty(sparam_config_name) || N == 0
    return
end

fid = fopen(filename, 'r');
if fid < 0
    return
end

in_block   = false;
cur_freq   = NaN;
cur_S      = complex(NaN(N,N), NaN(N,N));
have_data  = false;
freq_list  = [];
S_cell     = {};

    function flush_pass()
        if have_data
            freq_list(end+1) = cur_freq; %#ok<AGROW>
            S_cell{end+1} = cur_S; %#ok<AGROW>
        end
    end

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);

    if contains(line, 'CONFIGURATION:')
        tok = regexp(line, 'CONFIGURATION:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            this_config = tok{1}{1};
            if strcmp(this_config, sparam_config_name)
                flush_pass();   % save the pass we're leaving (if any)
                in_block  = true;
                cur_freq  = NaN;
                cur_S     = complex(NaN(N,N), NaN(N,N));
                have_data = false;
            else
                flush_pass();
                in_block  = false;
                have_data = false;
            end
        end
        continue
    end
    if ~in_block, continue; end

    if isnan(cur_freq)
        fm = regexp(line, 'Frequency in Hz:\s*FREQ\s*=\s*([\dE+\-\.]+)', 'tokens');
        if ~isempty(fm)
            cur_freq = str2double(fm{1}{1});
        end
    end

    sm = regexp(line, '^S\s+(\d+)\s+(\d+)\s+([+-]?[\dE+\-\.]+)\s+([+-]?[\dE+\-\.]+)', 'tokens');
    if ~isempty(sm)
        sink = round(str2double(sm{1}{1}));
        src  = round(str2double(sm{1}{2}));
        re_v = str2double(sm{1}{3});
        im_v = str2double(sm{1}{4});
        if sink >= 1 && sink <= N && src >= 1 && src <= N
            cur_S(sink, src) = complex(re_v, im_v);
            have_data = true;
        end
    end
end
flush_pass();   % last pass in the file, if the loop ended while still in_block
fclose(fid);

F = numel(freq_list);
if F == 0
    return
end
S_unsorted = complex(NaN(N,N,F), NaN(N,N,F));
for fi = 1:F
    S_unsorted(:,:,fi) = S_cell{fi};
end
[Sfreq, sidx] = sort(freq_list(:), 'ascend');
S = S_unsorted(:,:,sidx);
end


% =========================================================================
function [sparam_ports, sparam_ref_Z] = parse_sparam_config(filename, sparam_config_name)
% parse_sparam_config  From the first occurrence of the S-parameter
% configuration, read each excited port's name, signed edge indices and
% {negative, positive} face labels, and the per-port reference impedances from
% the "S-PARAMETER REFERENCE IMPEDANCES AT PORTS" table (port order = source
% order).
sparam_ports = struct('name', {}, 'edge_indices', {}, 'face_names', {});
sparam_ref_Z = containers.Map('KeyType','int32','ValueType','double');
if isempty(sparam_config_name)
    return
end

fid = fopen(filename, 'r');
if fid < 0
    return
end

in_block      = false;
started       = false;
in_ref_table  = false;
in_excitation = false;
face_pending  = false;

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);

    if contains(line, 'CONFIGURATION:')
        tok = regexp(line, 'CONFIGURATION:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            this_config = tok{1}{1};
            if strcmp(this_config, sparam_config_name)
                if started, break; end   % already collected the first occurrence
                in_block = true;
                started  = true;
            else
                if started, break; end   % left the block after collecting it
                in_block = false;
            end
        end
        continue
    end
    if ~in_block, continue; end

    if contains(line, 'S-PARAMETER REFERENCE IMPEDANCES AT PORTS')
        in_ref_table = true; continue
    end
    if in_ref_table
        tok = regexp(line, '^(\d+)\s+([\dE+\-\.]+)', 'tokens');
        if ~isempty(tok)
            sparam_ref_Z(int32(str2double(tok{1}{1}))) = str2double(tok{1}{2});
        elseif contains(line, 'Each active source')
            in_ref_table = false;
        end
    end

    if contains(line, 'EXCITATION BY VOLTAGE SOURCE')
        in_excitation = true;
        sparam_ports(end+1) = struct('name', '', 'edge_indices', [], 'face_names', {{'',''}});
        continue
    end
    if ~in_excitation, continue; end

    if contains(line, 'Attached to port:')
        tok = regexp(line, 'Attached to port:\s*(\S+)', 'tokens');
        if ~isempty(tok), sparam_ports(end).name = tok{1}{1}; end
    elseif contains(line, '(negative faces)')
        face_pending = true;
    elseif face_pending
        face_pending = false;
        face_tok = regexp(line, '(\S+)\s+(\S+)', 'tokens');
        if ~isempty(face_tok)
            sparam_ports(end).face_names = {face_tok{1}{1}, face_tok{1}{2}};
        end
    elseif contains(line, 'Indices of the edges:')
        tok = regexp(line, '(-?\d+)', 'tokens');
        eidx = [];
        for ti = 1:length(tok)
            eidx(end+1) = str2double(tok{ti}{1}); %#ok<AGROW>
        end
        sparam_ports(end).edge_indices = eidx;
        in_excitation = false;
    end
end
fclose(fid);
end


% =========================================================================
function out = parse_one_field_config(filename, config_name, triangles, edges_map, tri_labels)
% parse_one_field_config  Parse one single-port fields configuration (one
% driven port, other ports loaded): the driven port's V0, reference impedance,
% edges, location and polarity, its Zin over all frequency repeats (for the
% S-parameter consistency check only), and the first load's data.
%
% Inputs:
%   filename    - path to the .out file
%   config_name - configuration name
%   triangles, edges_map, tri_labels - mesh maps parsed by the caller
%
% Returns struct:
%   .config_name, .port_name
%   .V0           - complex source voltage (NaN if not found)
%   .Z0_ref       - reference impedance, Ohm (NaN if not found)
%   .edge_indices - signed edge indices; .face_names {negative, positive}
%   .R0, .pol     - 1x3 location / polarity from the mesh, polarity sign
%                   resolved from face labels (NaN if unresolved)
%   .Zin, .Freq   - (NFreq x 1) input impedance and frequency (Hz), ascending
%   .load_port_name, .load_edge_indices, .load_Z - the load block's port name,
%                   signed edges and complex impedance (Ohm)
out = struct('config_name', config_name, 'port_name', '', 'V0', NaN, ...
    'Z0_ref', NaN, 'edge_indices', [], 'face_names', {{'',''}}, ...
    'R0', [NaN NaN NaN], 'pol', [NaN NaN NaN], 'Zin', [], 'Freq', [], ...
    'load_port_name', '', 'load_edge_indices', [], 'load_Z', NaN);

fid = fopen(filename, 'r');
if fid < 0
    return
end

in_target         = false;
in_source_block   = false;
face_pending      = false;
got_port_info     = false;   % excitation info captured (first occurrence only)
cur_V0_mag = NaN; cur_V0_phase = NaN;

in_load_section    = false;
load_face_pending  = false; %#ok<NASGU>
got_load_info      = false;

in_vsource_block  = false;
zin_freq_list     = [];
zin_vals          = [];
last_freq         = NaN;

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);

    if contains(line, 'CONFIGURATION:')
        tok = regexp(line, 'CONFIGURATION:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            in_target = strcmp(tok{1}{1}, config_name);
        end
        continue
    end
    if ~in_target, continue; end

    fm = regexp(line, 'FREQ\s*=\s*([\dE+\-\.]+)', 'tokens');
    if ~isempty(fm)
        last_freq = str2double(fm{1}{1});
    end

    % --- Excitation block: name, V0, edges, face labels (first occurrence only) ---
    if ~got_port_info
        if contains(line, 'EXCITATION BY VOLTAGE SOURCE AT AN EDGE')
            in_source_block = true;
            cur_V0_mag = NaN; cur_V0_phase = NaN;
        elseif in_source_block && contains(line, '|U0| =')
            tok = regexp(line, '\|U0\|\s*=\s*([\dE+\-\.]+)', 'tokens');
            if ~isempty(tok), cur_V0_mag = str2double(tok{1}{1}); end
        elseif in_source_block && contains(line, 'ARG(U0) =')
            tok = regexp(line, 'ARG\(U0\)\s*=\s*([+-]?[\dE+\-\.]+)', 'tokens');
            if ~isempty(tok), cur_V0_phase = str2double(tok{1}{1}); end
        elseif in_source_block && contains(line, 'Attached to port:')
            tok = regexp(line, 'Attached to port:\s*(\S+)', 'tokens');
            if ~isempty(tok), out.port_name = tok{1}{1}; end
        elseif in_source_block && contains(line, '(negative faces)')
            face_pending = true;
        elseif in_source_block && face_pending
            face_pending = false;
            face_tok = regexp(line, '(\S+)\s+(\S+)', 'tokens');
            if ~isempty(face_tok)
                out.face_names = {face_tok{1}{1}, face_tok{1}{2}};
            end
        elseif in_source_block && contains(line, 'Indices of the edges:')
            tok = regexp(line, '(-?\d+)', 'tokens');
            eidx = [];
            for ti = 1:length(tok)
                eidx(end+1) = str2double(tok{ti}{1}); %#ok<AGROW>
            end
            out.edge_indices = eidx;
            in_source_block = false;
            if ~isnan(cur_V0_mag) && ~isnan(cur_V0_phase)
                out.V0 = cur_V0_mag * exp(1j * cur_V0_phase * pi/180);
            end
            got_port_info = true;
        end
    end

    if isnan(out.Z0_ref)
        tok = regexp(line, 'Reference impedance in Ohm:\s*([\dE+\-\.]+)', 'tokens');
        if ~isempty(tok)
            out.Z0_ref = str2double(tok{1}{1});
        end
    end

    % --- Load block: the OTHER port's termination (first occurrence only) ---
    if ~got_load_info
        if contains(line, 'DATA FOR EDGE LOADS')
            in_load_section = true;
        elseif contains(line, 'READING SOLUTION COEFFICIENTS') || contains(line, 'DATA FOR THE GREEN''S FUNCTION')
            in_load_section = false;
        end
        if in_load_section
            if contains(line, 'Complex impedance:')
                tok = regexp(line, '\(\s*([+-]?[\dE+\-\.]+)\s*\+j\s*([+-]?[\dE+\-\.]+)\)', 'tokens');
                if ~isempty(tok)
                    out.load_Z = complex(str2double(tok{1}{1}), str2double(tok{1}{2}));
                end
            elseif contains(line, 'Attached to port:')
                tok = regexp(line, 'Attached to port:\s*(\S+)', 'tokens');
                if ~isempty(tok), out.load_port_name = tok{1}{1}; end
            elseif contains(line, 'Indices of the edges:')
                tok = regexp(line, '(-?\d+)', 'tokens');
                eidx = [];
                for ti = 1:length(tok)
                    eidx(end+1) = str2double(tok{ti}{1}); %#ok<AGROW>
                end
                out.load_edge_indices = eidx;
                got_load_info = true;
            end
        end
    end

    % --- Zin over all frequency repeats of this configuration ---
    src_m = regexp(line, 'DATA OF THE VOLTAGE SOURCE NO\.\s*(\d+)', 'tokens');
    if ~isempty(src_m)
        in_vsource_block = true;
        continue
    end
    if in_vsource_block
        imp_m = regexp(line, 'Impedance\s+in\s+Ohm\s+([+-]?[\dE+\-\.]+)\s+([+-]?[\dE+\-\.]+)', 'tokens');
        if ~isempty(imp_m)
            re_z = str2double(imp_m{1}{1});
            im_z = str2double(imp_m{1}{2});
            zin_freq_list(end+1) = last_freq; %#ok<AGROW>
            zin_vals(end+1) = complex(re_z, im_z); %#ok<AGROW>
            in_vsource_block = false;
        end
    end
end
fclose(fid);

if ~isempty(zin_freq_list)
    [out.Freq, sidx] = sort(zin_freq_list(:), 'ascend');
    zv = zin_vals(:);
    out.Zin = zv(sidx);
end

% --- Geometry from the mesh, polarity sign resolved from face labels ---
if ~isempty(out.edge_indices) && ~isempty(triangles) && ~isempty(edges_map) && triangles.Count > 0 && edges_map.Count > 0
    [R0v, polv] = compute_edge_group_geometry(out.edge_indices, triangles, edges_map);
    out.R0  = R0v;
    out.pol = polv;
    if ~isempty(out.face_names) && ~isempty(out.face_names{2})
        pos_face_label = out.face_names{2};
        korp_label = '';
        for ei = 1:length(out.edge_indices)
            eidx_signed = out.edge_indices(ei);
            eidx = int32(abs(eidx_signed));
            if ~isKey(edges_map, eidx), continue; end
            edata = edges_map(eidx);
            if eidx_signed < 0
                pos_tri = int32(edata{1});
            else
                pos_tri = int32(edata{2});
            end
            if isKey(tri_labels, pos_tri)
                korp_label = tri_labels(pos_tri);
                break
            end
        end
        if ~isempty(korp_label)
            dot_idx = strfind(pos_face_label, '.');
            if ~isempty(dot_idx)
                pos_body = pos_face_label(1:dot_idx(1)-1);
            else
                pos_body = pos_face_label;
            end
            tri_label_width = 6;
            pos_body_cmp = pos_body(1:min(tri_label_width, length(pos_body)));
            if ~strcmp(pos_body_cmp, korp_label) && all(~isnan(out.pol))
                out.pol = -out.pol;
            end
        end
    end
end
end


% =========================================================================
function [loc, pol] = compute_edge_group_geometry(edge_idx_signed_list, triangles, edges_map)
% compute_edge_group_geometry  Location and polarity of one edge port from its
% signed edge indices: location = mean edge midpoint; polarity = normalised mean
% unit vector from each edge midpoint toward the free vertex of the positive
% triangle (negative index -> KORP, positive -> KORM), perpendicular to the
% edge. Outputs are [NaN NaN NaN] if they cannot be resolved.
loc = [NaN NaN NaN];
pol = [NaN NaN NaN];
group_mps = zeros(0,3);
pg_pol    = zeros(0,3);

for ei = 1:length(edge_idx_signed_list)
    eidx_signed = edge_idx_signed_list(ei);
    eidx = int32(abs(eidx_signed));
    if ~isKey(edges_map, eidx), continue; end
    edata = edges_map(eidx);
    korp = int32(edata{1}); korm = int32(edata{2});
    poip = edata{3};        poim = edata{4};
    if ~isKey(triangles, korp) || ~isKey(triangles, korm), continue; end

    if eidx_signed < 0
        pos_tri = korp; pos_poi = poip;
    else
        pos_tri = korm; pos_poi = poim;
    end
    verts_pos = triangles(pos_tri);
    if length(verts_pos) < 3 || pos_poi > length(verts_pos), continue; end

    all_vi = [1 2 3];
    shared_vi = all_vi(all_vi ~= pos_poi);
    sv1 = verts_pos{shared_vi(1)};
    sv2 = verts_pos{shared_vi(2)};
    group_mps(end+1,:) = (sv1 + sv2) / 2; %#ok<AGROW>

    edge_vec = sv2 - sv1;
    edge_len = norm(edge_vec);
    if edge_len < 1e-15, continue; end
    edge_hat = edge_vec / edge_len;
    free_pos_v = verts_pos{pos_poi};
    edge_mid   = (sv1 + sv2) / 2;
    vec_to_free = free_pos_v - edge_mid;
    delta_perp  = vec_to_free - dot(vec_to_free, edge_hat) * edge_hat;
    if norm(delta_perp) > 1e-15
        pg_pol(end+1,:) = delta_perp / norm(delta_perp); %#ok<AGROW>
    end
end

if ~isempty(group_mps)
    loc = mean(group_mps, 1);
end
if ~isempty(pg_pol)
    mean_pol = mean(pg_pol, 1);
    if norm(mean_pol) > 1e-12
        pol = mean_pol / norm(mean_pol);
    end
end
end