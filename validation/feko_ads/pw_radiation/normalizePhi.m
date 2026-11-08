function s = normalizePhi(s, phiList, mode, phiUnit)
%NORMALIZEPHI  Convert a data set's phi values to degrees and match them to the requested cuts.
%   s = NORMALIZEPHI(s, phiList, mode, phiUnit)
%
%   phiUnit ("deg" default, or "rad") is the angular unit of the sweep
%   variable in the file; "rad" applies 180/pi (ADS phi_prime set in
%   radians). phiList and the returned s.phi are in degrees. Only s.phi is
%   changed; s.data (including phase, in degrees in both sources) is not.
%
%   mode = "auto"  : if every phi in phiList already has a column, nothing
%                    changes. Otherwise the factors pi/180, 180/pi, 1e-3, 1e3
%                    are tried (after the phiUnit conversion) and the first
%                    that makes every requested cut line up is applied, with
%                    a message naming the file. If none works, s is returned
%                    unchanged and checkPhiCoverage reports what is missing.
%   mode = 1       : off, use the phi values as exported.
%   mode = <number>: apply that factor unconditionally.
%
%   The "auto" rescaling covers exports written in an unexpected angular
%   unit (e.g. radians already converted to degrees, reading 1719/3438/5157);
%   correcting the export is preferable.
%
%   See also CHECKPHICOVERAGE, FINDPHICOLUMN.

    if nargin < 3 || isempty(mode)
        mode = "auto";
    end
    if nargin < 4 || isempty(phiUnit)
        phiUnit = "deg";
    end

    % --- Declared unit of the source ----------------------------------
    if lower(string(phiUnit)) == "rad"
        s.phi = s.phi * 180/pi;
    end

    % --- Explicit numeric factor --------------------------------------
    if isnumeric(mode)
        if mode ~= 1
            s.phi = s.phi * mode;
        end
        return;
    end

    if mode ~= "auto"
        return;
    end

    % --- Already correct? ---------------------------------------------
    if all(arrayfun(@(p) ~isempty(findPhiColumn(s.phi, p)), phiList))
        return;
    end

    % --- Try the known unit slips -------------------------------------
    candidates = [pi/180, 180/pi, 1e-3, 1e3];
    labels     = ["pi/180 (radians exported as degrees)", ...
                  "180/pi (degrees exported as radians)", ...
                  "1e-3 (millidegrees)", ...
                  "1e3"];

    for i = 1:numel(candidates)
        scaled = s.phi * candidates(i);
        if all(arrayfun(@(p) ~isempty(findPhiColumn(scaled, p)), phiList))
            [~, nm, ext] = fileparts(s.file);
            fprintf(['  phi rescaled by %s in %s%s: [%s] -> [%s] deg.\n' ...
                     '  Fix the export when you can; this is a workaround.\n'], ...
                labels(i), nm, ext, ...
                join(string(s.phi), ", "), join(string(round(scaled,3)), ", "));
            s.phi = scaled;
            return;
        end
    end
end