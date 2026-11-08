function s = restrictFreqRange(s, fLimHz)
%RESTRICTFREQRANGE  Keep only the samples inside a frequency band.
%   s = RESTRICTFREQRANGE(s, fLimHz) keeps the rows of s.freq/s.data with
%   fLimHz(1) <= freq <= fLimHz(2). Pass [] to leave the data untouched.
%
%   Used to clip the FEKO sweep (0.1-6 GHz) to the band actually covered
%   by the ADS sweep (0.1-2.5 GHz), so the two are compared only where
%   both exist -- which matters most on the polar/phasor plot, where a
%   longer FEKO trace would otherwise look like a disagreement.

    if nargin < 2 || isempty(fLimHz)
        return;
    end

    tol  = 1e-6 * max(abs(fLimHz));
    keep = s.freq >= fLimHz(1) - tol & s.freq <= fLimHz(2) + tol;

    s.freq = s.freq(keep);
    s.data = s.data(keep, :);
end
