%% ImportAndDisplayLegendreTables.m
% Load FTheta_table.txt and FPhi_table.txt (written by
% mathematica/LegendreFThetaFPhi.wl into mathematica/tables), convert the
% entries to symbolic expressions and print each with pretty().
% Requires the Symbolic Math Toolbox.
%
% Table format (tab-delimited text, so that Excel does not reinterpret
% the expressions):
%   Row 1   : header  ("l/m", -lMax, ..., -1, 0, 1, ..., lMax)
%   Col 1   : l values
%   Cols 2+ : expression strings; entries with |m| > l are "0".

clear; clc;

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

dataDir = fullfile(repo_root, 'mathematica', 'tables');

thetaFile = fullfile(dataDir, 'FTheta_table.txt');
phiFile   = fullfile(dataDir, 'FPhi_table.txt');

%% ---- Load the tables as strings ----
% 'preserve' keeps headers such as "-2" instead of sanitising them.
opts = detectImportOptions(thetaFile, 'FileType', 'text', ...
    'TextType', 'string', 'Delimiter', '\t');
opts.VariableNamingRule = 'preserve';
opts.VariableNamesLine = 1;
opts.DataLines = [2, Inf];
% Import every column as string; otherwise a column that is mostly "0"
% placeholders is detected as numeric and its expressions become NaN.
opts.VariableTypes(:) = {'string'};
Table_theta = readtable(thetaFile, opts);

opts2 = detectImportOptions(phiFile, 'FileType', 'text', ...
    'TextType', 'string', 'Delimiter', '\t');
opts2.VariableNamingRule = 'preserve';
opts2.VariableNamesLine = 1;
opts2.DataLines = [2, Inf];
opts2.VariableTypes(:) = {'string'};
Table_phi = readtable(phiFile, opts2);

% Column 1: l values (read as strings, converted to numbers).
lVals = double(Table_theta{:, 1});

rawTheta = Table_theta{:, 2:end};
rawPhi   = Table_phi{:, 2:end};

% m values from the column headers.
mHeaders = string(Table_theta.Properties.VariableNames(2:end));
mVals = double(mHeaders);

fprintf('Raw headers found: %s\n', strjoin(mHeaders, ', '));
fprintf('Parsed m values:   %s\n\n', strjoin(string(mVals), ', '));

%% ---- Parse to symbolic (variable th) and display ----
syms th

dP_Theta_syms = parseTable(rawTheta, lVals, mVals);
P_Phi_syms   = parseTable(rawPhi, lVals, mVals);

fprintf('\n===== F_theta(l,m) =====\n');
displayTable(dP_Theta_syms, lVals, mVals, 'F_theta');

fprintf('\n===== F_phi(l,m) =====\n');
displayTable(P_Phi_syms, lVals, mVals, 'F_phi');

%% ================= Local functions ================= %%

function out = parseTable(rawStrArray, lVals, mVals)
    % str2sym on every cell; empty cells -> 0, unparsable cells -> []
    % (reported on stderr).
    [nRows, nCols] = size(rawStrArray);
    out = cell(nRows, nCols);
    for i = 1:nRows
        for j = 1:nCols
            s = rawStrArray(i, j);

            if ismissing(s) || strlength(s) == 0
                out{i, j} = sym(0);
                continue
            end

            s = strtrim(s);

            try
                out{i, j} = str2sym(s);
            catch ME
                fprintf(2, ['Failed to parse entry (l=%s, m=%s).\n' ...
                    'Raw string was: "%s"\nError: %s\n\n'], ...
                    safeStr(lVals(i)), safeStr(mVals(j)), s, ME.message);
                out{i, j} = [];
            end
        end
    end
end

function displayTable(symCellArray, lVals, mVals, label)
    % Print every non-empty entry with pretty().
    [nRows, nCols] = size(symCellArray);
    for i = 1:nRows
        for j = 1:nCols
            entry = symCellArray{i, j};
            if isempty(entry)
                continue
            end
            fprintf('%s(l=%d, m=%d) =\n', label, lVals(i), mVals(j));
            pretty(entry)
            fprintf('\n');
        end
    end
end

function str = safeStr(x)
    % Display string for a scalar, including NaN/<missing>.
    if ismissing(x)
        str = "<missing/NaN>";
    else
        str = string(x);
    end
end