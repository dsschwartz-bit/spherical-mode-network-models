function [resolved_path, resolved_name] = resolve_field_file(fileprefix, raw_name, occurrence, reuse_type, ext)
% RESOLVE_FIELD_FILE  Find the exported FEKO field file for a (possibly
% reused) field-request name.
%
% Inputs (from ReadFEKOout_loadports):
%   fileprefix - path prefix of the exports (string), e.g. "<dir>/antenna"
%   raw_name   - field-request name, e.g. 'NearField1'
%   occurrence - number of earlier configurations using raw_name (0 = first)
%   reuse_type - 'same_freq' or 'diff_freq' (ignored when occurrence == 0)
%   ext        - extension without the dot, e.g. "efe"
%
% occurrence == 0: only <fileprefix>_<raw_name>.<ext> is tried.
% occurrence  > 0: FEKO renames reused exports ("fr" = frequency range):
%   'same_freq': <fileprefix>_<raw_name>(<occurrence>).<ext> tried first
%   'diff_freq': <fileprefix>_fr_<occurrence>_<raw_name>.<ext> tried first
% then the other pattern, then the plain name. The first existing file wins.
%
% Returns:
%   resolved_path - path of the file found, or "" if no candidate exists
%                   (the caller then treats the request as unavailable)
%   resolved_name - the <name> part used (without prefix/extension), so a
%                   second extension (e.g. .hfe after .efe) can be matched
%                   to the same candidate; "" if nothing was found
%
% Used by compute_port_mode_coeffs_nearfield, _farfield and _Q.
if occurrence == 0
    candidates = {char(raw_name)};
else
    same_freq_candidate = sprintf('%s(%d)', raw_name, occurrence);
    diff_freq_candidate = sprintf('fr_%d_%s', occurrence, raw_name);
    if strcmpi(reuse_type, 'diff_freq')
        candidates = {diff_freq_candidate, same_freq_candidate, char(raw_name)};
    else   % 'same_freq' or unclassified
        candidates = {same_freq_candidate, diff_freq_candidate, char(raw_name)};
    end
end
resolved_path = ""; resolved_name = "";
for ci = 1:numel(candidates)
    trial_path = fileprefix + "_" + candidates{ci} + "." + ext;
    if isfile(trial_path)
        resolved_path = trial_path;
        resolved_name = string(candidates{ci});
        return
    end
end
end
