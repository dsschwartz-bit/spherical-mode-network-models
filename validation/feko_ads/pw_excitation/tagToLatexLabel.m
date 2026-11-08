function label = tagToLatexLabel(tag, source)
%TAGTOLATEXLABEL  Legend label from a load tag (e.g. "ph30th60", incidence phi / theta) and a source
%   name (e.g. "ADS" or "FEKO"):  "ADS, $(\phi,\theta) = (30^\circ,60^\circ)$".
%   The source is plain text outside the math (MATLAB's LaTeX interpreter drops spaces in math mode).
%   Single-quoted format: sprintf resolves the backslashes itself.

    tokens = regexp(tag, 'ph(\d+)th(\d+)', 'tokens', 'once');
    label = sprintf('%s, $(\\phi,\\theta) = (%s^\\circ,%s^\\circ)$', source, tokens{1}, tokens{2});
end
