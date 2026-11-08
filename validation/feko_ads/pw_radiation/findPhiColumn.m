function idx = findPhiColumn(phiVec, phiWanted, tolDeg)
%FINDPHICOLUMN  Index of the column whose phi matches phiWanted.
%   idx = FINDPHICOLUMN(phiVec, phiWanted, tolDeg) returns the index of
%   the entry of phiVec within tolDeg (default 0.25 deg) of phiWanted, or
%   [] if there is none. This is how ADS curves and FEKO curves are
%   paired: by the phi value they represent, never by column order.
%
%   The default tolerance is loose on purpose. ADS sweeps phi_prime in
%   radians and exports it with three decimals, so 0.524 / 1.047 / 1.571
%   convert to 30.023 / 59.989 / 90.012 deg -- correct cuts, carrying up to
%   ~0.02 deg of print rounding. A tolerance of a few hundredths of a degree
%   would reject them; 0.25 deg absorbs the rounding while staying far below
%   any realistic spacing between cuts.

    if nargin < 3 || isempty(tolDeg), tolDeg = 0.25; end

    idx = find(abs(phiVec(:).' - phiWanted) <= tolDeg, 1);
end