function y = EvalLegendreF(FFun, lList, mList, lQuery, mQuery, THETA)
%EVALLEGENDREF  Evaluate one entry of a LoadLegendreFHandles table.
%   y = EvalLegendreF(FFun, lList, mList, lQuery, mQuery, THETA) finds
%   the entry of the cell array FFun (FThetaFun or FPhiFun) for
%   l = lQuery, m = abs(mQuery) and evaluates it at THETA (radians).
%   Only abs(mQuery) is used, matching gradYlm and
%   legendre_derivative_modified, which evaluate the Legendre function at
%   |m| and apply the sign of m through the cos/sin azimuthal factor.
%
%   See also LoadLegendreFHandles.

    i = find(lList == lQuery, 1);
    j = find(mList == abs(mQuery), 1);

    if isempty(i) || isempty(j)
        error('EvalLegendreF:NotFound', ...
            ['No table entry found for l=%d, m=%d (abs used for ' ...
             'lookup). Check that the Wolfram script''s lValues ' ...
             'range covers this l.'], lQuery, mQuery);
    end

    y = FFun{i, j}(THETA);
end