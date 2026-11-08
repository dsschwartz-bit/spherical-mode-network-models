function timeDomainLegend(S)
%TIMEDOMAINLEGEND  Shared single-row legend for the time-domain figures.
%   TIMEDOMAINLEGEND(S) opens a legend figure (legendFigure, S from plotStyle) with
%   entries FEKO (.S2P) markers, ADS SnP modes markers,
%   ADS rational-fit modes line and V_gen (thick semi-transparent line).
%   Called by plot_impulse_comparison.m and plot_RFpulse_comparison.m.
fig = figure('Visible', 'off');  ax = axes(fig);  hold(ax, 'on');
h(1) = plot(ax, NaN, NaN, "LineStyle", "none", "Marker", "d", "Color", S.colFEKO);
h(2) = plot(ax, NaN, NaN, "LineStyle", "none", "Marker", "v", "Color", S.colors(4,:));   % Wong bluish green
h(3) = plot(ax, NaN, NaN, "LineStyle", "-", "Color", S.colADS);
h(4) = plot(ax, NaN, NaN, "LineStyle", "-", "Color", [S.colExc S.alphaExc], "LineWidth", S.lwExc, "Tag", "keepWidth");
styleTraces(h, S);                           % same line/marker styling as the panels
legendFigure("TimeDomain_legend", h, ...
    ["FEKO (.S2P)", "ADS (SnP modes)", "ADS (rational-fit modes)", "$V_{\mathrm{gen}}$"], 4);
close(fig);
end
