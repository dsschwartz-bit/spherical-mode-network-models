function scale = cst_unit_scale(unit_str)
%CST_UNIT_SCALE  Scale factor from a CST unit string to the SI base unit.
%   scale = CST_UNIT_SCALE(unit_str) returns the factor converting a value
%   in unit_str (length, frequency or voltage/field) to metres, Hz or volts,
%   e.g. cst_unit_scale('mm') = 1e-3. Case-insensitive; errors on an
%   unrecognised unit. The unit strings are set in the "Unit Configuration"
%   section of general_vector_spherical_modes_loadports_CST.m.
switch lower(strtrim(unit_str))
    % Length -> metres
    case 'm',   scale = 1;
    case 'cm',  scale = 1e-2;
    case 'mm',  scale = 1e-3;
    case 'um',  scale = 1e-6;
    case 'nm',  scale = 1e-9;
    % Frequency -> Hz
    case 'hz',  scale = 1;
    case 'khz', scale = 1e3;
    case 'mhz', scale = 1e6;
    case 'ghz', scale = 1e9;
    % Voltage / field -> V (field would use the same prefixes over V/m)
    case {'v', 'v/m'},   scale = 1;
    case {'mv', 'mv/m'}, scale = 1e-3;
    case {'kv', 'kv/m'}, scale = 1e3;
    otherwise
        error('cst_unit_scale: unrecognized unit "%s".', unit_str);
end
