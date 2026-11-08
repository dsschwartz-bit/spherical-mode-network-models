function port_names = get_port_names_from_snp_comments(filepath)
%GET_PORT_NAMES_FROM_SNP_COMMENTS  Port names from the "! Port[N] Name" comments of a Touchstone file.
%   PORT_NAMES = GET_PORT_NAMES_FROM_SNP_COMMENTS(FILEPATH) reads the
%   comment lines before the first data line and returns a 1 x N cell array
%   of names in port order; ports without a comment get 'Port<N>'. Errors
%   if no such comments are found. Same parsing as the ADS Python helpers.
    fid = fopen(filepath, 'r');
    if fid == -1
        error('Could not open %s for reading.', filepath);
    end

    port_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
    pattern = '^!\s*Port\[(\d+)\]\s+(\S+)';

    line = fgetl(fid);
    while ischar(line)
        line = strtrim(line);
        if ~startsWith(line, '!')
            if ~isempty(line) && ~startsWith(line, '#')
                break
            end
            line = fgetl(fid);
            continue
        end
        tokens = regexpi(line, pattern, 'tokens');
        if ~isempty(tokens)
            idx = str2double(tokens{1}{1});
            name = tokens{1}{2};
            port_map(idx) = name; %#ok<NASGU>
        end
        line = fgetl(fid);
    end
    fclose(fid);

    if port_map.Count == 0
        error('No "! Port[N] Name" comments found in %s', filepath);
    end

    max_port = max(cell2mat(keys(port_map)));
    port_names = cell(1, max_port);
    for i = 1:max_port
        if isKey(port_map, i)
            port_names{i} = port_map(i);
        else
            port_names{i} = sprintf('Port%d', i);
        end
    end
end
