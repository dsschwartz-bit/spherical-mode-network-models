function add_gmin_to_rational_fit_netlist(ckt_path, output_path, gmin_ohms)
%ADD_GMIN_TO_RATIONAL_FIT_NETLIST  Add gmin resistors from internal netlist nodes to ground.
%   ADD_GMIN_TO_RATIONAL_FIT_NETLIST(CKT_PATH) adds a 1e12 Ohm resistor
%   "Rgmin_<node> <node> 0 <gmin_ohms>" from every internal node of a
%   rational-fit SPICE netlist to ground, just before .ENDS, to help
%   convergence, and overwrites CKT_PATH.
%   ADD_GMIN_TO_RATIONAL_FIT_NETLIST(CKT_PATH, OUTPUT_PATH) writes to
%   OUTPUT_PATH instead ([] or '' = overwrite CKT_PATH).
%   ADD_GMIN_TO_RATIONAL_FIT_NETLIST(CKT_PATH, OUTPUT_PATH, GMIN_OHMS) sets
%   the resistance.
%
%   Ground ("0") and the external nodes on the .SUBCKT line are skipped:
%   they are connected in the parent schematic, and a gmin resistor there
%   only triggers an ADS "Resistor is shorted" warning.
%   Nodes are identified by element type (first letter of the name):
%       V, R, C, L : 2 nodes (tokens 2-3)
%       G, E       : 4 nodes (tokens 2-5; output +/-, control +/-)
%       F, H       : 2 nodes (tokens 2-3; token 4 is the controlling
%                    voltage source name, not a node)
%   generateSPICE netlists use V, R, C, G and F; other types raise an error.
%
%   Example:
%     add_gmin_to_rational_fit_netlist( ...
%         fullfile(repo_root,'output','SnP','rationalfit','S_dom_Q_bowtie_edgeport_45degrees_PW_excitation.s6p.ckt'), ...
%         fullfile(repo_root,'output','SnP','rationalfit','S_dom_Q_bowtie_edgeport_45degrees_PW_excitation_gmin.s6p.ckt'))

    if nargin < 2 || isempty(output_path)
        output_path = ckt_path;
    end
    if nargin < 3 || isempty(gmin_ohms)
        gmin_ohms = 1e12;
    end

    % Node-argument count per element type (see help above)
    two_node_types = 'VRCL';
    four_node_types = 'GE';
    two_node_skip_one_types = 'FH';

    % ---- Read the raw .ckt file ----
    content = fileread(ckt_path);
    lines = splitlines(content);

    % ---- Find the .SUBCKT and .ENDS line indices ----
    subckt_idx = 0;
    ends_idx = 0;
    for k = 1:numel(lines)
        line_k = strtrim(lines{k});
        if numel(line_k) >= 7 && strcmpi(line_k(1:7), '.subckt')
            subckt_idx = k;
        elseif numel(line_k) >= 5 && strcmpi(line_k(1:5), '.ends')
            ends_idx = k;
            break
        end
    end
    if subckt_idx == 0
        error('No ".SUBCKT" line found in %s', ckt_path);
    end
    if ends_idx == 0
        error('No ".ENDS" line found in %s', ckt_path);
    end

    % ---- External nodes (.SUBCKT line) -- excluded from gmin below ----
    subckt_parts = strsplit(strtrim(lines{subckt_idx}));
    % subckt_parts{1} = ".SUBCKT", subckt_parts{2} = subckt name,
    % subckt_parts{3:end} = external nodes.
    external_nodes = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for t = 3:numel(subckt_parts)
        external_nodes(subckt_parts{t}) = true;
    end
    fprintf('Found %d external node(s) on the .SUBCKT line -- excluded from gmin: %s\n', ...
            numel(subckt_parts) - 2, strjoin(subckt_parts(3:end), ', '));

    % ---- Collect every distinct internal node name ----
    % containers.Map used as a set (values unused); keys() returns them sorted.
    internal_nodes = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    for k = (subckt_idx + 1):(ends_idx - 1)
        line_k = strtrim(lines{k});
        if isempty(line_k) || line_k(1) == '*'
            continue
        end
        tokens = strsplit(line_k);
        element_letter = upper(tokens{1}(1));
        if any(element_letter == two_node_types)
            node_idxs = [2, 3];
        elseif any(element_letter == four_node_types)
            node_idxs = [2, 3, 4, 5];
        elseif any(element_letter == two_node_skip_one_types)
            node_idxs = [2, 3];
        else
            error(['Unrecognized element type ''%s'' on line %d: %s -- ' ...
                   'add its own node-argument count to one of ' ...
                   'two_node_types/four_node_types/two_node_skip_one_types ' ...
                   'above before re-running.'], element_letter, k, line_k);
        end
        for idx = node_idxs
            if idx > numel(tokens)
                error('Line %d has fewer tokens than expected for element type ''%s'': %s', k, element_letter, line_k);
            end
            node_name = tokens{idx};
            if ~strcmp(node_name, '0') && ~isKey(external_nodes, node_name)
                internal_nodes(node_name) = true;
            end
        end
    end

    node_list = keys(internal_nodes);
    fprintf('Found %d distinct internal node(s) (excluding ground and external pins) across %d element line(s).\n', ...
            numel(node_list), ends_idx - subckt_idx - 1);

    % ---- Build the new gmin resistor lines ----
    gmin_lines = cell(numel(node_list), 1);
    for i = 1:numel(node_list)
        gmin_lines{i} = sprintf('Rgmin_%s %s 0 %g', node_list{i}, node_list{i}, gmin_ohms);
    end

    % ---- Insert them immediately before .ENDS ----
    lines = [lines(1:ends_idx-1); gmin_lines; lines(ends_idx:end)];
    content = strjoin(lines, sprintf('\r\n'));

    % ---- Write the result ----
    fid = fopen(output_path, 'w');
    if fid == -1
        error('Could not open %s for writing.', output_path);
    end
    fwrite(fid, content);
    fclose(fid);

    fprintf('Added %d gmin resistor(s) of %g Ohm each, from every internal node to ground (external pins skipped).\n', numel(node_list), gmin_ohms);
    fprintf('Wrote result to %s\n', output_path);
end