function snp_source = get_snp_source_from_netlist_comment(ckt_path)
%GET_SNP_SOURCE_FROM_NETLIST_COMMENT  SnP file recorded in a rational-fit netlist.
%   SNP_SOURCE = GET_SNP_SOURCE_FROM_NETLIST_COMMENT(CKT_PATH) returns the
%   path from the "* SNP_SOURCE: <path>" comment written by
%   FIX_RATIONAL_FIT_NETLIST (char), or '' if the netlist has none.
%   Same parsing as the ADS Python helper.
    content = fileread(ckt_path);
    lines = splitlines(content);

    snp_source = '';
    pattern = '^\*\s*SNP_SOURCE:\s*(.+)$';
    for k = 1:numel(lines)
        line_k = strtrim(lines{k});
        tokens = regexpi(line_k, pattern, 'tokens');
        if ~isempty(tokens)
            snp_source = strtrim(tokens{1}{1});
            return
        end
    end
end
