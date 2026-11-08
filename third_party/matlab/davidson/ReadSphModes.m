function [Q,j_,j3D,mmax,nmax] = ReadSphModes(flnm)
%READSPHMODES  Read spherical-mode Q coefficients from a TICRA .sph file.
%   [Q, j_, j3D, mmax, nmax] = READSPHMODES(flnm) reads one frequency of
%   spherical-wave coefficients in TICRA (GRASP) .sph format.
%   Q    - complex coefficients indexed by the compressed single index
%          j_ = 2*(n*(n+1) + m - 1) + s  (s = 1 TE, s = 2 TM; see FEKO manual)
%   j_   - compressed indices in file order
%   j3D  - string labels "s m n" indexed by j_
%   mmax, nmax - maximum azimuthal and polar mode orders from the file header
%   See also ReadSphModes_freqs (multi-frequency version).
%
% DB Davidson, Feb 2024.

modedata=importfile_sph(flnm); % first two file lines skipped
headerlines=6; 
nmax=modedata{1,3};
mmax=modedata{1,4};

mode_counter=1;

lines = headerlines+1;
modecheck=modedata{lines,1};
lines = lines+1;
m = 0;
if modecheck ~= m 
    error(['Mode ',num2str(m),' incorrectly read'])
end
for n=1:nmax
    for s=1:2 % TE then TM modes
        j_(mode_counter)  = 2 *(n*(n+1)+ m-1) + s; % compressed mode index
        j3D(j_(mode_counter)) = string([num2str(s),' ',num2str(m),' ',num2str(n)]);
        if s == 1
            Q(j_(mode_counter)) = complex(modedata{lines,1},modedata{lines,2});
        else
            Q(j_(mode_counter)) = complex(modedata{lines,3},modedata{lines,4});
        end
        mode_counter = mode_counter +1;
    end
    lines = lines+1;
end

for mmode = 1:mmax
    modecheck=modedata{lines,1};
    if modecheck ~= mmode
        error(['Mode ',num2str(mmode),' incorrectly read'])
    end
    lines = lines+1; %Advance to next mode.
    for n=abs(mmode):nmax
        for m = -mmode:2*mmode:mmode
            for s=1:2 % TE then TM modes
                j_(mode_counter) = 2 *(n*(n+1)+ m-1) + s; % compressed mode index
                j3D(j_(mode_counter)) = string([num2str(s),' ',num2str(m),' ',num2str(n)]);
                if s == 1
                    Q(j_(mode_counter)) = complex(modedata{lines,1},modedata{lines,2});
                else
                    Q(j_(mode_counter)) = complex(modedata{lines,3},modedata{lines,4});
                end
                mode_counter = mode_counter +1;
            end
            lines = lines+1;
        end
    end
end