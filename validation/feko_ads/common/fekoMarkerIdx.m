function idx = fekoMarkerIdx(freq, stepHz, k, nTraces, x, y)
%FEKOMARKERIDX Indices of the FEKO samples to show as markers (light downsampling).
%   idx = FEKOMARKERIDX(freq, stepHz, k, nTraces): one sample every stepHz (nearest real sample),
%   with the grid of trace k (of nTraces) shifted by (k-1)/nTraces*stepHz so the markers of
%   overlapping traces are staggered instead of stacked.
%   idx = FEKOMARKERIDX(freq, nMarkers, k, nTraces, x, y): (polar / parametric curves) nMarkers
%   samples evenly spaced along the arc length of the curve (x, y), also staggered per trace.
%   stepHz <= 0 or empty keeps every sample.
freq = freq(:);
if isempty(stepHz) || stepHz <= 0, idx = (1:numel(freq)).'; return; end
if nargin < 3, k = 1; end
if nargin < 4, nTraces = 1; end
shift = (k-1)/nTraces;
if nargin >= 6 && ~isempty(x)
    s = [0; cumsum(hypot(diff(x(:)), diff(y(:))))];     % arc length along the curve
    target = ((0:stepHz-1).' + shift) * s(end)/stepHz;
    ref = s;
else
    target = (freq(1) + shift*stepHz : stepHz : freq(end)).';
    ref = freq;
end
idx = zeros(numel(target), 1);
for i = 1:numel(target), [~, idx(i)] = min(abs(ref - target(i))); end
idx = unique(idx, 'stable');
end
