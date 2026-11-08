function h = polarLabel(fig, S, txt)
%POLARLABEL  LaTeX label in the top-left corner of a polar figure.

h = annotation(fig, 'textbox', [0.005 0.925 0.6 0.07], 'String', txt, 'Interpreter', 'latex', ...
    'FontSize', S.labelSize, 'EdgeColor', 'none', 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', 'FitBoxToText', 'off');
end
