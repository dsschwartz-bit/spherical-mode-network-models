function s = downsampleFeko(s, stepHz)
%DOWNSAMPLEFEKO  Thin a data set down to a coarser frequency step.
%   s = DOWNSAMPLEFEKO(s, stepHz) keeps, for each point of the grid
%   s.freq(1) : stepHz : s.freq(end), the nearest sample of s.freq/s.data
%   (rows of s.data follow s.freq). Used to turn a dense FEKO sweep into a
%   readable set of markers.
%
%   If stepHz is empty, <= 0, or finer than the native step, s is returned
%   unchanged. Data sets on the same frequency grid decimated with the same
%   stepHz stay aligned point-for-point.

    if nargin < 2 || isempty(stepHz) || stepHz <= 0
        return;
    end

    f = s.freq;
    if numel(f) < 2 || stepHz <= median(diff(f))
        return;
    end

    targets = f(1) : stepHz : f(end);
    idx = zeros(1, numel(targets));
    for i = 1:numel(targets)
        [~, idx(i)] = min(abs(f - targets(i)));
    end
    idx = unique(idx, "stable");

    s.freq = f(idx);
    s.data = s.data(idx, :);
end
