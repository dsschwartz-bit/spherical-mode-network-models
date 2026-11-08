function checkPhiCoverage(s, phiList)
%CHECKPHICOVERAGE  Warn if a data set is missing any of the requested phi cuts.
%   CHECKPHICOVERAGE(s, phiList) checks that every phi in phiList has a
%   matching column in s. If some are missing, it warns and -- when the
%   file's phi values are a constant multiple of the requested ones --
%   names the phiScale that would line them up.
%
%   The common case is a sweep entered in radians but exported in degrees,
%   which turns 30/60/90 into 1719/3438/5157; phiScale = pi/180 fixes it.
%   Nothing is changed here: this only tells you what to set.
%
%   See also IMPORTFIELDDATA, FINDPHICOLUMN.

    missing = phiList(arrayfun(@(p) isempty(findPhiColumn(s.phi, p)), phiList));
    if isempty(missing)
        return;
    end

    [~, nm, ext] = fileparts(s.file);
    msg = sprintf("%s%s (%s %s): no column for phi = [%s] deg. File contains phi = [%s].", ...
        nm, ext, s.source, s.quantity, ...
        join(string(missing), ", "), join(string(s.phi), ", "));

    % Does a constant factor reconcile the two lists?
    candidates = [pi/180, 180/pi, 1/1000, 1000];
    names      = ["pi/180 (the sweep looks like radians exported as degrees)", ...
                  "180/pi (the sweep looks like degrees exported as radians)", ...
                  "1/1000 (the sweep looks like millidegrees)", ...
                  "1000 (the sweep looks like kilo-units)"];

    for i = 1:numel(candidates)
        scaled = s.phi * candidates(i);
        if all(arrayfun(@(p) ~isempty(findPhiColumn(scaled, p, 1e-2)), phiList))
            msg = msg + sprintf("\nSet the phiScale for this source to %s.", names(i));
            break;
        end
    end

    warning("checkPhiCoverage:missingPhi", "%s", msg);
end
