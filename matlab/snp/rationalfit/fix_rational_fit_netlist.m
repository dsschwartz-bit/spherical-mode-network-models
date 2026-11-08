function fix_rational_fit_netlist(ckt_path, snp_path, output_path)
%FIX_RATIONAL_FIT_NETLIST  Rename generateSPICE netlist nodes to the SnP port names.
%   FIX_RATIONAL_FIT_NETLIST(CKT_PATH, SNP_PATH) renames the generic node
%   names of a .ckt file written by generateSPICE(fit, path, zref,
%   "SeparateReferenceNodes", true) to the port names in SNP_PATH's
%   "! Port[N] Name" comments, and overwrites CKT_PATH.
%   FIX_RATIONAL_FIT_NETLIST(CKT_PATH, SNP_PATH, OUTPUT_PATH) writes to
%   OUTPUT_PATH instead.
%
%   The 2N nodes on the .SUBCKT line are taken as N (signal, reference)
%   pairs in SnP port order and renamed to "<name>" and "<name>_ref"
%   throughout the file. generateSPICE keeps the port order of the
%   sparameters object, but check the result once for a new file type.
%   A "* SNP_SOURCE: <snp_path>" comment is prepended (see
%   GET_SNP_SOURCE_FROM_NETLIST_COMMENT).
%
%   Example:
%     fix_rational_fit_netlist( ...
%         fullfile(repo_root,'output','SnP','rationalfit','S_dom_Q_loop_loadport.s8p.ckt'), ...
%         fullfile(repo_root,'output','SnP','S_dom_Q_loop_loadport.s8p'))

    if nargin < 3
        output_path = ckt_path;
    end

    % ---- Parse port names from the SnP file's "! Port[N] Name" comments ----
    port_names = get_port_names_from_snp_comments(snp_path);
    fprintf('Found %d ports: %s\n', numel(port_names), strjoin(port_names, ', '));

    % ---- Read the raw .ckt file ----
    content = fileread(ckt_path);

    % ---- Find the .SUBCKT line and extract its node list ----
    lines = splitlines(content);
    subckt_line = '';
    for k = 1:numel(lines)
        line_k = strtrim(lines{k});
        if numel(line_k) >= 7 && strcmpi(line_k(1:7), '.subckt')
            subckt_line = line_k;
            break
        end
    end
    if isempty(subckt_line)
        error('No ".SUBCKT" line found in %s', ckt_path);
    end

    parts = strsplit(subckt_line);
    % parts{1} = ".SUBCKT", parts{2} = subckt name, parts{3:end} = nodes
    old_nodes = parts(3:end);
    fprintf('Found %d nodes on the .SUBCKT line: %s\n', numel(old_nodes), strjoin(old_nodes, ', '));

    if numel(old_nodes) ~= 2*numel(port_names)
        error(['Found %d nodes on the .SUBCKT line but expected %d ' ...
               '(2 per port for %d ports). Check that generateSPICE was ' ...
               'called with "SeparateReferenceNodes", true, and that ' ...
               'snp_path is the correct corresponding SnP file.'], ...
              numel(old_nodes), 2*numel(port_names), numel(port_names));
    end

    % ---- Map old -> new names: nodes (2i-1, 2i) -> port i signal, reference ----
    ref_suffix = '_ref';
    old_to_new = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for i = 1:numel(port_names)
        sig_idx = 2*i - 1;
        ref_idx = 2*i;
        old_to_new(old_nodes{sig_idx}) = port_names{i};
        old_to_new(old_nodes{ref_idx}) = [port_names{i} ref_suffix];
        fprintf('  %s -> %s,  %s -> %s\n', old_nodes{sig_idx}, port_names{i}, old_nodes{ref_idx}, [port_names{i} ref_suffix]);
    end

    % ---- Apply the renaming throughout the file ----
    % Exact match of whitespace-delimited tokens (SPICE netlists are
    % token based), so partial matches such as "p1" in "p10" cannot occur.
    total_count = 0;
    for k = 1:numel(lines)
        line_k = lines{k};
        tokens = strsplit(line_k);
        changed = false;
        for t = 1:numel(tokens)
            if isKey(old_to_new, tokens{t})
                tokens{t} = old_to_new(tokens{t});
                total_count = total_count + 1;
                changed = true;
            end
        end
        if changed
            lines{k} = strjoin(tokens, ' ');
        end
    end
    content = strjoin(lines, sprintf('\r\n'));

    % ---- Record the source SnP file in a comment line ----
    content = [sprintf('* SNP_SOURCE: %s\r\n', snp_path) content];

    % ---- Write the corrected file ----
    fid = fopen(output_path, 'w');
    if fid == -1
        error('Could not open %s for writing.', output_path);
    end
    fwrite(fid, content);
    fclose(fid);

    fprintf('Renamed %d node token(s) throughout the file.\n', total_count);
    fprintf('Recorded source SnP file in a ''* SNP_SOURCE:'' comment line.\n');
    fprintf('Wrote corrected netlist to %s\n', output_path);
end
