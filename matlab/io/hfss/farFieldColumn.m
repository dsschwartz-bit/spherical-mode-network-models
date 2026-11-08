function col = farFieldColumn(entry, headerName)
%FARFIELDCOLUMN Look up a column by its exact csv header string.
%   col = FARFIELDCOLUMN(entry, headerName) returns entry.data(:,idx)
%   where idx is the position of headerName in entry.headers (as
%   produced by IMPORTFFARFIELD_HFSS). Errors with the list of available
%   headers if headerName isn't found, so a typo doesn't fail silently.
%
%   Example:
%       ff = importffarfield_HFSS('rectangular_patch_load_port_farfield_rE');
%       entry = ff([ff.port]==1 & [ff.mode]==1);
%       theta = farFieldColumn(entry, 'Theta(deg)');
%       magPhi = farFieldColumn(entry, 'mag(rEPhi) (V)');

    idx = find(strcmp(entry.headers, headerName), 1);
    if isempty(idx)
        error('farFieldColumn:notFound', ...
            'Column "%s" not found. Available columns:\n  %s', ...
            headerName, strjoin(entry.headers, '\n  '));
    end
    col = entry.data(:, idx);
end