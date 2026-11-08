function vars = ReadFEKOcfx(filename)
% ReadFEKOcfx  Read named variable definitions from a CADFEKO .cfx file.
%
% The .cfx file is HDF5. It is read with MATLAB's built-in h5info, h5read and
% h5readatt; no extra toolbox is needed. The first group named Variables,
% GlobalVariables, DesignVariables, ModelVariables, UserVariables or
% Expressions (breadth-first search) is read: its datasets and its scalar
% numeric attributes become struct fields.
%
% Output:
%   vars - struct, one field per variable (name made a valid MATLAB
%          identifier), e.g. vars.a = 0.075, vars.rho = 0.002.
%          Empty struct if no variables are found.
%
% Usage:
%   vars = ReadFEKOcfx('bicon_edgeport_symmetry.cfx')
%   a   = vars.a

vars = struct();

% --- Get file structure ---
try
    info = h5info(filename);
catch ME
    error('ReadFEKOcfx: Cannot read HDF5 file %s.\nError: %s', filename, ME.message);
end

% --- Search for the variable group (e.g. /Variables, /Model/Variables) ---
var_group = find_variable_group(info, filename);

if isempty(var_group)
    fprintf('ReadFEKOcfx: No variable definitions found in %s.\n', filename);
    fprintf('ReadFEKOcfx: Use h5disp(''%s'') to inspect the file structure.\n', filename);
    return
end

fprintf('ReadFEKOcfx: Found variables at path: %s\n', var_group);

% --- Read all datasets in the variable group ---
try
    ginfo = h5info(filename, var_group);
catch
    fprintf('ReadFEKOcfx: Could not read group %s.\n', var_group);
    return
end

% Read datasets (each variable is typically a dataset)
for i = 1:length(ginfo.Datasets)
    name  = ginfo.Datasets(i).Name;
    try
        val = h5read(filename, [var_group '/' name]);
        % Convert to scalar if possible
        if numel(val) == 1
            val = double(val);
        end
        % Make valid MATLAB field name
        safe_name = matlab.lang.makeValidName(name);
        vars.(safe_name) = val;
        fprintf('  %s = %g\n', name, val);
    catch
        fprintf('  %s : could not read value\n', name);
    end
end

% Read attributes (some CADFEKO versions store variables as attributes)
for i = 1:length(ginfo.Attributes)
    name = ginfo.Attributes(i).Name;
    try
        val = h5readatt(filename, var_group, name);
        if isnumeric(val) && numel(val) == 1
            safe_name = matlab.lang.makeValidName(name);
            vars.(safe_name) = double(val);
            fprintf('  %s = %g (attribute)\n', name, val);
        end
    catch
    end
end

if isempty(fieldnames(vars))
    fprintf('ReadFEKOcfx: Variable group found but no readable values extracted.\n');
    fprintf('ReadFEKOcfx: Use h5disp(''%s'', ''%s'') to inspect manually.\n', filename, var_group);
end

end


function group_path = find_variable_group(info, filename)
% Recursively search for a group whose name suggests it holds variables.

group_path = '';
candidates = {'Variables', 'GlobalVariables', 'DesignVariables', ...
              'ModelVariables', 'UserVariables', 'Expressions'};

% BFS over groups
queue = {{'/', info}};
while ~isempty(queue)
    item      = queue{1};
    queue(1)  = [];
    path      = item{1};
    node      = item{2};

    if ~isfield(node, 'Groups'), continue; end
    for g = 1:length(node.Groups)
        grp      = node.Groups(g);
        % grp.Name is the full path
        parts    = strsplit(grp.Name, '/');
        leafname = parts{end};
        for c = 1:length(candidates)
            if strcmpi(leafname, candidates{c})
                group_path = grp.Name;
                return
            end
        end
        % recurse
        try
            child = h5info(filename, grp.Name);
            queue{end+1} = {grp.Name, child}; %#ok<AGROW>
        catch
        end
    end
end

end