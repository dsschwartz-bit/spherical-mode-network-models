function mindx = getmindx(m, lmax)
% getmindx  Convert a physics azimuthal order m to its MATLAB index into mval.
%
% mval is defined as (-lmax:lmax).', so the index of m is simply m + lmax + 1.
%
% Usage:
%   mindx = getmindx(m, lmax)
%
% Example (lmax = 3, mval = -3:3):
%   getmindx(0,  3) -> 4
%   getmindx(1,  3) -> 5
%   getmindx(-3, 3) -> 1
%   mval(getmindx(1, 3)) -> 1   (recovers physics m value)
%
% Inputs:
%   m    - scalar or array of physics azimuthal orders (integers, |m| <= lmax)
%   lmax - maximum spherical mode order (defines the mval array)
%
% Output:
%   mindx - MATLAB index/indices into mval such that mval(mindx) == m

if any(abs(m) > lmax)
    error('getmindx: |m| must not exceed lmax (%d)', lmax)
end

mindx = m + lmax + 1;

end