function dPnm = legendre_derivative_modified(n, m, x, Pnm, normalization)
%LEGENDRE_DERIVATIVE_MODIFIED  dPbar_n^m/dx for a single (n,m) at all x.
%   dPnm = legendre_derivative_modified(n, m, x)
%   dPnm = legendre_derivative_modified(n, m, x, Pnm)
%   dPnm = legendre_derivative_modified(n, m, x, Pnm, normalization)
%
%   Single-(n,m) adaptation of legendre_derivative (Oldenhuis, in
%   third_party/) that evaluates only the orders it needs with
%   FastLegendre_Plm instead of the full m = 0..n sweep of legendre().
%
%   With normalization 'norm' (default)
%       Pbar_n^m(x) = (-1)^m sqrt((n+1/2)(n-m)!/(n+m)!) P_n^m(x),
%   where P_n^m from FastLegendre_Plm already carries the Condon-Shortley
%   phase, so Pbar has none. 'unnorm' uses P_n^m itself. Recurrence used:
%       m > 0:  dPbar_n^m/dx = [m x Pbar_n^m - sqrt((n+m)(n-m+1)) sqrt(1-x^2) Pbar_n^(m-1)] / (1-x^2)
%       m = 0:  dPbar_n^0/dx = sqrt(n(n+1)) Pbar_n^1 sqrt(1-x^2) / (1-x^2)
%   At x = +-1 the known limits are used (0 for m >= 2, -x*Inf for m = 1).
%
%   Inputs:
%     n             non-negative integer degree
%     m             integer order, |m| <= n; |m| is used throughout
%     x             real values in [-1,1], any shape (output has the same shape)
%     Pnm           optional precomputed Pbar_n^|m|(x) in the same
%                   normalisation: one row (1 x p), or an (n+1) x p matrix
%                   from which row |m|+1 is taken
%     normalization 'norm' (default) or 'unnorm'
%
%   Example:
%     th = linspace(0, pi, 101);
%     dP = legendre_derivative_modified(l, m, cos(th));
%
%   See also FastLegendre_Plm, gradYlm, legendre.

    %% Input validation
    narginchk(3, 5);

    assert(isnumeric(n) && isscalar(n) && isfinite(n) && isreal(n) && ...
           round(n)==n && n >= 0, ...
           'legendre_derivative_modified:invalid_n', ...
           'n must be a non-negative integer scalar.');

    assert(isnumeric(m) && isscalar(m) && isfinite(m) && isreal(m) && ...
           round(m)==m && abs(m) <= n, ...
           'legendre_derivative_modified:invalid_m', ...
           'm must be an integer scalar with |m| <= n.');

    m = abs(m);

    assert(~isempty(x) && isnumeric(x) && all(isreal(x(:))) && ...
           all(abs(x(~isnan(x(:)))) <= 1), ...
           'legendre_derivative_modified:invalid_x', ...
           'x must be real, numeric, nonempty, and in [-1, 1].');

    if nargin < 5 || isempty(normalization)
        normalization = 'norm';
    end
    assert(any(strcmpi(normalization, {'norm','unnorm'})), ...
           'legendre_derivative_modified:invalid_normalization', ...
           'normalization must be ''norm'' or ''unnorm''.');

    % Normalisation factor N_n^m = (-1)^m sqrt((n+1/2)(n-m)!/(n+m)!).
    % Pnm and the adjacent-order Pnm_adj must use the same normalisation.
    switch normalization
        case 'norm'
            normfac = @(nn,mm) (-1)^mm * sqrt((nn+0.5) * factorial(nn-mm) / factorial(nn+mm));
        otherwise  % unnorm
            normfac = @(nn,mm) 1;
    end

    xrow = x(:).';

    if nargin < 4 || isempty(Pnm)
        Pnm = normfac(n,m) * FastLegendre_Plm(n, m, xrow);
        Pnm = Pnm(:).';
    else
        if size(Pnm, 1) == n+1
            Pnm = Pnm(m+1, :);
        end
        Pnm = Pnm(:).';
    end

    sqx = 1 - xrow.^2;   % sin^2(theta)

    if m == 0
        % Uses Pbar_n^1 in place of Pbar_n^(-1).
        Pnm_adj = normfac(n,1) * FastLegendre_Plm(n, 1, xrow);
        Pnm_adj = Pnm_adj(:).';
        dPnm    = sqrt(n*(n+1)) .* Pnm_adj .* sqrt(sqx) ./ sqx;
    else
        Pnm_adj  = normfac(n,m-1) * FastLegendre_Plm(n, m-1, xrow);
        Pnm_adj  = Pnm_adj(:).';
        coupling = sqrt((n+m) * (n-m+1));
        dPnm     = (m .* xrow .* Pnm - coupling .* sqrt(sqx) .* Pnm_adj) ./ sqx;
    end

    % |x| = 1: the recurrence is 0/0, use the known limits.
    ispolar = abs(xrow) == 1;
    if any(ispolar)
        xp  = xrow(ispolar);   % 1 x npolar
        xpc = xp(:);           % npolar x 1 for column-safe logical indexing

        dPnm(ispolar) = 0;     % default: limit is 0 for m >= 2

        if m == 0
            % x^(n+1) * N_n^0 * n(n+1)/2
            dPnm(ispolar) = (xpc .* xp(:).^n * normfac(n,0) * n*(n+1)/2).';
        elseif m == 1
            dPnm(ispolar) = (-xpc .* inf).';
        end
    end

    dPnm = reshape(dPnm, size(x));

end