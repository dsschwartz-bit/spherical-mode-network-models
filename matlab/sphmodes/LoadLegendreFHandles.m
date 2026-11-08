function [lList, mList, FThetaFun, FPhiFun] = LoadLegendreFHandles(dataDir)
%LOADLEGENDREFHANDLES  Load the Mathematica F_theta/F_phi tables as function handles.
%   [lList, mList, FThetaFun, FPhiFun] = LoadLegendreFHandles(dataDir)
%   reads FTheta_table.txt and FPhi_table.txt (tab-delimited, written by
%   mathematica/LegendreFThetaFPhi.wl) and converts every entry into a
%   vectorised function handle of theta. Requires the Symbolic Math Toolbox.
%
%   INPUTS:
%     dataDir - folder containing the two tables
%               (default: fullfile(repo_root,'mathematica','tables'))
%
%   OUTPUTS:
%     lList     - column vector of l values (rows of the table)
%     mList     - row vector of m values, -lMax..lMax (columns of the table)
%     FThetaFun - L x M cell array of function handles.
%                 FThetaFun{i,j}(THETA) evaluates
%                     -sin(theta) * d(P_bar_l^m)/dx |_{x=cos(theta)}
%                 at l = lList(i), m = mList(j).
%     FPhiFun   - L x M cell array of function handles.
%                 FPhiFun{i,j}(THETA) evaluates
%                     P_bar_l^m(cos theta) / sin(theta)
%                 at l = lList(i), m = mList(j).
%
%   The table normalisation cancels the Condon-Shortley phase in a way
%   that is not symmetric in m -> -m (entries for m and -m differ by
%   (-1)^m). Look entries up with abs(m), as gradYlm and
%   legendre_derivative_modified do; EvalLegendreF does this.
%
%   See also EvalLegendreF, gradYlm.

    if nargin < 1 || isempty(dataDir)
        dataDir = fullfile(repo_root, 'mathematica', 'tables');
    end

    thetaFile = fullfile(dataDir, 'FTheta_table.txt');
    phiFile   = fullfile(dataDir, 'FPhi_table.txt');

    opts = detectImportOptions(thetaFile, 'FileType', 'text', ...
        'TextType', 'string', 'Delimiter', '\t');
    opts.VariableNamingRule = 'preserve';
    opts.VariableNamesLine = 1;
    opts.DataLines = [2, Inf];
    opts.VariableTypes(:) = {'string'};
    Table_theta = readtable(thetaFile, opts);

    opts2 = detectImportOptions(phiFile, 'FileType', 'text', ...
        'TextType', 'string', 'Delimiter', '\t');
    opts2.VariableNamingRule = 'preserve';
    opts2.VariableNamesLine = 1;
    opts2.DataLines = [2, Inf];
    opts2.VariableTypes(:) = {'string'};
    Table_phi = readtable(phiFile, opts2);

    lList = double(Table_theta{:, 1});
    mList = double(string(Table_theta.Properties.VariableNames(2:end)));

    rawTheta = Table_theta{:, 2:end};
    rawPhi   = Table_phi{:, 2:end};

    syms th

    [nRows, nCols] = size(rawTheta);
    FThetaFun = cell(nRows, nCols);
    FPhiFun   = cell(nRows, nCols);

    for i = 1:nRows
        for j = 1:nCols
            [FThetaFun{i, j}, exprTh] = parseAndConvert( ...
                rawTheta(i, j), th, lList(i), mList(j));
            [FPhiFun{i, j}, exprPh] = parseAndConvert( ...
                rawPhi(i, j), th, lList(i), mList(j));
            %#ok<*NASGU>
        end
    end
end

function [fun, expr] = parseAndConvert(s, th, lVal, mVal)
    % Parse one table cell. Empty cells and cells that str2sym cannot
    % parse (reported on stderr) become a constant-zero function.
    % matlabFunction generates elementwise code, so the handle accepts
    % a full THETA array.
    if ismissing(s) || strlength(s) == 0
        expr = sym(0);
        fun = matlabFunction(expr, 'Vars', th);
        return
    end

    s = strtrim(s);

    try
        expr = str2sym(s);
        fun = matlabFunction(expr, 'Vars', th);
    catch ME
        fprintf(2, ['Failed to parse entry (l=%s, m=%s).\n' ...
            'Raw string was: "%s"\nError: %s\n\n'], ...
            safeStr(lVal), safeStr(mVal), s, ME.message);
        expr = sym(0);
        fun = matlabFunction(expr, 'Vars', th);
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