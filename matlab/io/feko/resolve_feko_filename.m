function active_file = resolve_feko_filename(filename, request_keyword)
% resolve_feko_filename  Try filename as-is; if not found, toggle '_fr_1'
% before the request keyword (e.g. '_NearField1', '_FarField2', '_SParameter').
%
% If the given filename contains '_fr_1_<keyword>', tries removing '_fr_1'.
% If it does not contain '_fr_1_<keyword>', tries adding '_fr_1' before keyword.
% Number after keyword is optional.
%
% Usage:
%   active_file = resolve_feko_filename('antenna_NearField1.efe', '_NearField')
%   active_file = resolve_feko_filename('antenna_fr_1_FarField2.sph', '_FarField')
%   active_file = resolve_feko_filename('antenna_SParameter.s1p', '_SParameter')

if exist(filename, 'file')
    active_file = filename;
else
    fr1_keyword = ['_fr_1' request_keyword];
    if contains(filename, fr1_keyword)
        % filename has fr_1 -> try without it
        alt_file = strrep(filename, fr1_keyword, request_keyword);
    else
        % filename does not have fr_1 -> try with it
        alt_file = regexprep(filename, ['(' request_keyword '[0-9]*)'], '_fr_1$1');
    end
    if exist(alt_file, 'file')
        fprintf('resolve_feko_filename: Could not find %s, using %s instead.\n', filename, alt_file);
        active_file = alt_file;
    else
        error('Cannot open file: %s (also tried: %s)', filename, alt_file);
    end
end
end