function DecodeSphIndex(j_)
% DecodeSphIndex  Decode a compressed spherical wave index j_ into l, m, s.
%
% The single index is defined as:
%   j_ = 2*(l*(l+1) + m - 1) + s
% where s=1 (TE), s=2 (TM), m = -l..l, l = 1..lmax
%
% Usage:
%   DecodeSphIndex(4)   % j_=4 -> l=1, m=0, TM
%   DecodeSphIndex([4, 24, 60])

for k = 1:length(j_)
    jval = j_(k);
    
    % s is determined by odd/even
    s = mod(jval, 2);
    if s == 0
        s = 2;
    end
    type = "TM";
    if s == 1
        type = "TE";
    end
    
    % l*(l+1) + m = (jval - s) / 2 + 1
    lm = (jval - s) / 2 + 1;  % = l*(l+1) + m

    % find l: l*(l+1) bracket contains lm for m in -l..l
    % i.e. l*(l+1) - l <= lm <= l*(l+1) + l
    l = 1;
    while l*(l+1) + l < lm
        l = l + 1;
    end
    
    % m is the remainder
    m = lm - l*(l+1);
    
    fprintf("j_ = %d  ->  l = %d, m = %d, %s\n", jval, l, m, type);
end

end