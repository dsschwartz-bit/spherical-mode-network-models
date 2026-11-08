function [files_by_port, freqs] = discover_cst_farfield_files(farfield_directory)
%DISCOVER_CST_FARFIELD_FILES  Find CST far-field exports and group them by
% port and frequency.
%   Matches CST's default export name "farfield (f=<freq>) [<port>].txt"
%   (e.g. "farfield (f=6) [1].txt": f = 6 in the project frequency unit,
%   port 1 excited; case-insensitive). Other files are ignored. Every port
%   must have the same set of frequencies, otherwise an error is raised.
%
% Returns:
%   files_by_port - 1 x num_ports cell; files_by_port{p} is a string array of
%                   file names (no folder) for the p-th port in ascending
%                   port number, sorted by frequency
%   freqs         - 1 x numf sorted unique frequencies as written in the
%                   file names (not scaled to Hz; apply cst_unit_scale)
listing = dir(fullfile(farfield_directory, 'farfield*.txt'));
if isempty(listing)
    error('discover_cst_farfield_files: no "farfield*.txt" files found in "%s".', farfield_directory);
end

n = numel(listing);
port_of = nan(1,n);
freq_of = nan(1,n);
name_of = cell(1,n);
for k = 1:n
    tok = regexp(listing(k).name, '^farfield\s*\(f\s*=\s*([\d.eE+-]+)\)\s*\[(\d+)\]\.txt$', 'tokens', 'ignorecase');
    if isempty(tok)
        continue
    end
    freq_of(k) = str2double(tok{1}{1});
    port_of(k) = str2double(tok{1}{2});
    name_of{k} = listing(k).name;
end
valid = ~isnan(port_of);
if ~any(valid)
    error('discover_cst_farfield_files: found %d "farfield*.txt" file(s) in "%s", but none matched the expected "farfield (f=<freq>) [<port>].txt" naming.', n, farfield_directory);
end
port_of = port_of(valid); freq_of = freq_of(valid); name_of = name_of(valid);

ports = sort(unique(port_of));
freqs = sort(unique(freq_of));
num_ports = numel(ports);
files_by_port = cell(1, num_ports);
for p = 1:num_ports
    mask = port_of == ports(p);
    port_freqs = sort(freq_of(mask));
    if numel(port_freqs) ~= numel(freqs) || max(abs(port_freqs(:) - freqs(:))) > 0
        error('discover_cst_farfield_files: port %d has a different set of frequencies than the others (found %s, expected %s) -- every port must share the same frequency list.', ...
            ports(p), mat2str(port_freqs), mat2str(freqs));
    end
    [~, order] = sort(freq_of(mask));
    names_p = name_of(mask);
    files_by_port{p} = string(names_p(order));
end
end
