function j_ = SphModeIndex(l, m, type)
% SphModeIndex  Compute single spherical wave index j_ from l, m, s.
%
% j_ = 2*(l*(l+1) + m - 1) + s,  s=1 TE, s=2 TM
%
% Usage:
%   j_ = SphModeIndex(1, 0, 'TM')   % -> 4
%   j_ = SphModeIndex(3, 0, 'TM')   % -> 24

if abs(m)>l || l==0
    j_=[];
else
    if strcmpi(type, 'TE')
        s = 1;
    elseif strcmpi(type, 'TM')
        s = 2;
    else
        error("type must be 'TE' or 'TM'");
    end
    
    j_ = 2*(l.*(l+1) + m - 1) + s;
end
end