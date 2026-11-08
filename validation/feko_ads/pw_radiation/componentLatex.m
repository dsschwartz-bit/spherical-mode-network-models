function [symLatex, nameLatex] = componentLatex(component)
%COMPONENTLATEX  LaTeX strings for an rE far-field voltage component name.
%   [symLatex, nameLatex] = COMPONENTLATEX("Vtheta") returns
%     symLatex  = 'rE_{\theta}'   (no $ signs -- for embedding)
%     nameLatex = '$rE_{\theta}$'
%   Accepts "Vtheta" / "Vphi" (and plain "theta"/"phi").

    component = lower(string(component));

    if contains(component, "theta")
        symLatex = 'rE_{\theta}';
    elseif contains(component, "phi")
        symLatex = 'rE_{\phi}';
    else
        symLatex = char("rE_{" + component + "}");
    end

    nameLatex = ['$' symLatex '$'];
end
