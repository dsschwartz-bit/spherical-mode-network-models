function root = repo_root()
%REPO_ROOT  Absolute path of the repository root folder.
%   Use it to build paths to repository files, e.g.
%   fullfile(repo_root, 'mathematica', 'tables').

root = fileparts(fileparts(mfilename('fullpath')));
end
