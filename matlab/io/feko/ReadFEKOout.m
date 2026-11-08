function [V0, Z0_ref, R0, symmetry_planes, even_odd_excitation, antenna_apothem, Zin, Freq, port_pol, internal_loads] = ReadFEKOout(filename)
% ReadFEKOout  Parse a FEKO .out file for source, port, symmetry and input
% impedance data.
%
% Port data are read from the first frequency of the first configuration in
% the file. Wire ports (voltage source at a segment / wire node) and edge ports
% (voltage source at an edge) are supported.
%
% Outputs:
%   V0                  - (Nports x 1) complex source voltage |U0|*exp(j*ARG(U0)),
%                         one entry per source in file order; doubled if an
%                         infinite PEC ground plane (or a PMC ground plane with an
%                         in-plane wire port) coincides with the port
%   Z0_ref              - reference impedance (Ohm); x2 for a coincident PEC
%                         ground plane, /2 for a coincident PMC ground plane
%   R0                  - 1x3 port location (m): mean of the wire-port locations,
%                         or of the edge-port mesh-edge midpoints ([0 0 0] if
%                         unresolved)
%   symmetry_planes     - struct .XY/.XZ/.YZ (+1 magnetic/PMC wall, 0 none,
%                         -1 electric/PEC wall) and .physical_gnd (1 if an
%                         infinite ground plane coincides with the port; use
%                         2^physical_gnd as the ground-plane factor)
%   even_odd_excitation - +log2(Nports) if the ports are mirrored across a
%                         magnetic wall (even/common mode), -log2(Nports) across
%                         an electric wall (odd/differential mode), 0 otherwise
%   antenna_apothem     - max |r| over all mesh vertices (m), about the origin
%                         (not R0); requested interactively if the file has no
%                         triangle mesh
%   Zin                 - (NFreq x NSource) complex input impedance (Ohm) from the
%                         "DATA OF THE VOLTAGE SOURCE NO. X" blocks
%   Freq                - (NFreq x 1) frequencies (Hz), ascending, rows of Zin
%   port_pol            - (Nports x 3) port polarity unit vectors:
%                         wire port: stated positive direction (corrected for a
%                         separating PEC wall, see below);
%                         edge port: from the mesh, perpendicular to the port
%                         edge toward the free vertex of the positive triangle,
%                         averaged over the port's edges; [0 0 1] if unresolved
%   internal_loads      - struct array of the loads in the "DATA FOR EDGE LOADS"
%                         section (fields .name, .index, .circuit_type, .Z
%                         (complex, Ohm), .attached_port, .edge_indices (signed),
%                         .face_names {negative, positive}); 0x0 if none. A
%                         loaded but never-excited port gives a DC path that the
%                         FEKO workflows use for the 0 Hz S-parameter point.
%
% Edge-port polarity: the sign of each entry on the "Indices of the edges:" line
% selects the positive mesh triangle (negative -> KORP, positive -> KORM); the
% result is flipped if that triangle's body differs from the stated positive face.
%
% Usage:
%   [V0, Z0_ref, R0, symmetry_planes, even_odd_excitation, antenna_apothem, Zin, Freq] = ReadFEKOout('dipole_cylinder.out')
%
% See also ReadFEKOout_loadports.

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

Zin              = [];    % input impedance (NFreq x NSource)
port_type        = '';    % 'wire' or 'edge'
first_config     = '';    % name of first configuration block
excitation_type = '';    % e.g. 'Voltage Source at a Segment', 'Voltage Source at a Wire-Node'
port_locations  = [];    % accumulate all port locations (Nx3)
edge_midpoints  = zeros(0,3);  % edge-port locations from the mesh (stay empty without a mesh)
edge_indices    = {};    % cell array of edge index groups, one cell per port
edge_face_names = {};   % cell array of {neg_face, pos_face} strings, one cell per port
face_line_pending = false;
first_port_done  = false; % tracks whether at least one port has been seen
all_ports_done   = false;
current_R0       = NaN(1,3);
in_source_block  = false;
first_freq_val   = NaN;   % frequency value of first source block
seen_src_index1  = false; % tracks if source index 1 has been seen once already

% --- Load blocks (DATA FOR EDGE LOADS), first frequency only -> internal_loads ---
in_load_section   = false;
load_face_pending = false;
internal_loads = struct('name', {}, 'index', {}, 'circuit_type', {}, 'Z', {}, ...
    'attached_port', {}, 'edge_indices', {}, 'face_names', {});

% --- Zin collection ---
Zin              = [];    % will be (NFreq x NSource) after parsing
zin_data         = {};    % cell: {freq_idx, src_no, Re, Im}
in_vsource_block = false; % true when inside a DATA OF THE VOLTAGE SOURCE block
cur_vsrc_no      = 0;     % current source number
cur_vsrc_freq_idx = 0;    % current frequency index for Zin
zin_freq_list    = [];    % ordered list of frequencies seen in vsource blocks
last_freq        = NaN;   % most recently parsed FREQ = value

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

    % --- Detect configuration block: stop if a new config starts ---
    if contains(line, 'CONFIGURATION:')
        tok = regexp(line, 'CONFIGURATION:\s*(\S+)', 'tokens');
        if ~isempty(tok)
            this_config = tok{1}{1};
            if isempty(first_config)
                first_config = this_config;
            elseif ~strcmp(this_config, first_config)
                all_ports_done = true;
                continue
            end
        end
    end

    % --- Stop if load data section starts (S-param solve pass) ---
    if contains(line, 'DATA FOR VERTEX LOADS') || contains(line, 'READING SOLUTION COEFFICIENTS')
        all_ports_done = true;
    end

    % --- Track last seen frequency (used for Zin block attribution) ---
    % Also stop port collection when a new frequency appears (fallback for files
    % that do not contain a 'READING SOLUTION COEFFICIENTS' line).
    freq_m = regexp(line, 'FREQ\s*=\s*([\dE+\-.]+)', 'tokens');
    if ~isempty(freq_m)
        last_freq = str2double(freq_m{1}{1});
        if ~all_ports_done
            if isnan(first_freq_val)
                first_freq_val = last_freq;
            elseif abs(last_freq - first_freq_val) > first_freq_val * 1e-6
                all_ports_done = true;
            end
        end
        fidx = find(abs(zin_freq_list - last_freq) < last_freq*1e-6, 1);
        if isempty(fidx)
            zin_freq_list(end+1) = last_freq;
            cur_vsrc_freq_idx = length(zin_freq_list);
        else
            cur_vsrc_freq_idx = fidx;
        end
    end

    % --- Parse Zin from DATA OF THE VOLTAGE SOURCE blocks ---
    src_m = regexp(line, 'DATA OF THE VOLTAGE SOURCE NO\.\s*(\d+)', 'tokens');
    if ~isempty(src_m)
        in_vsource_block = true;
        cur_vsrc_no = str2double(src_m{1}{1});
        continue
    end
    if in_vsource_block
        imp_m = regexp(line, 'Impedance\s+in\s+Ohm\s+([+-]?[\dE+\-.]+)\s+([+-]?[\dE+\-.]+)', 'tokens');
        if ~isempty(imp_m) && cur_vsrc_freq_idx > 0
            re_z = str2double(imp_m{1}{1});
            im_z = str2double(imp_m{1}{2});
            zin_data{end+1} = {cur_vsrc_freq_idx, cur_vsrc_no, re_z, im_z};
            in_vsource_block = false;
        end
    end

    % --- Port collection (skipped after first frequency) ---
    if ~all_ports_done

    % --- Load blocks (DATA FOR EDGE LOADS). Must precede the
    % "~in_source_block -> continue" below. ---
    if contains(line, 'DATA FOR EDGE LOADS')
        in_load_section = true;
    elseif contains(line, 'READING SOLUTION COEFFICIENTS') || contains(line, 'DATA FOR THE GREEN''S FUNCTION')
        in_load_section = false;
    end

    if in_load_section
        if contains(line, 'Name:')
            tok = regexp(line, 'Name:\s*(.+)$', 'tokens');
            internal_loads(end+1) = struct('name', '', 'index', NaN, 'circuit_type', '', ...
                'Z', NaN, 'attached_port', '', 'edge_indices', [], 'face_names', {{'',''}});
            if ~isempty(tok)
                internal_loads(end).name = strtrim(tok{1}{1});
            end
        elseif ~isempty(internal_loads) && contains(line, 'Load index:')
            tok = regexp(line, 'Load index:\s*(\d+)', 'tokens');
            if ~isempty(tok), internal_loads(end).index = str2double(tok{1}{1}); end
        elseif ~isempty(internal_loads) && contains(line, 'Circuit type:')
            tok = regexp(line, 'Circuit type:\s*(.+)$', 'tokens');
            if ~isempty(tok), internal_loads(end).circuit_type = strtrim(tok{1}{1}); end
        elseif ~isempty(internal_loads) && contains(line, 'Complex impedance:')
            tok = regexp(line, '\(\s*([+-]?[\dE+\-\.]+)\s*\+j\s*([+-]?[\dE+\-\.]+)\)', 'tokens');
            if ~isempty(tok)
                internal_loads(end).Z = complex(str2double(tok{1}{1}), str2double(tok{1}{2}));
            end
        elseif ~isempty(internal_loads) && contains(line, 'Attached to port:')
            tok = regexp(line, 'Attached to port:\s*(\S+)', 'tokens');
            if ~isempty(tok), internal_loads(end).attached_port = tok{1}{1}; end
        elseif ~isempty(internal_loads) && contains(line, '(negative faces)')
            load_face_pending = true;
        elseif ~isempty(internal_loads) && load_face_pending
            load_face_pending = false;
            face_tok = regexp(strtrim(line), '(\S+)\s+(\S+)', 'tokens');
            if ~isempty(face_tok)
                internal_loads(end).face_names = {face_tok{1}{1}, face_tok{1}{2}};
            end
        elseif ~isempty(internal_loads) && contains(line, 'Indices of the edges:')
            tok = regexp(line, '(-?\d+)', 'tokens');
            eidx_list = [];
            for ti = 1:length(tok)
                eidx_list(end+1) = str2double(tok{ti}{1});
            end
            internal_loads(end).edge_indices = eidx_list;
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


    % --- Reference impedance ---
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

    end % if ~all_ports_done

end % while

fclose(fid);

% --- Assemble Zin matrix (NFreq x NSource) and Freq vector ---
Zin = []; Freq = [];
if ~isempty(zin_data)
    nfreqs_zin = length(zin_freq_list);
    nsrcs_zin  = max(cellfun(@(x) x{2}, zin_data));
    Zin = NaN(nfreqs_zin, nsrcs_zin);
    Freq = zin_freq_list(:);  % column vector
    for zi = 1:length(zin_data)
        fi  = zin_data{zi}{1};
        si  = zin_data{zi}{2};
        re  = zin_data{zi}{3};
        im  = zin_data{zi}{4};
        Zin(fi, si) = complex(re, im);
    end
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
    fprintf('ReadFEKOout: WARNING - antenna apothem could not be determined from the mesh data in this .out file (no "DATA OF THE METALLIC TRIANGLES" section found, or it contained no usable vertex data).\n');
    antenna_apothem = [];
    while isempty(antenna_apothem) || ~isnumeric(antenna_apothem) || ~isscalar(antenna_apothem) || ~isreal(antenna_apothem) || antenna_apothem <= 0
        antenna_apothem = input('Enter antenna apothem (meters): ');
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
[Freq,indicies]=sort(Freq,'ascend');
Zin=Zin(indicies,:);

fprintf('\n')

end