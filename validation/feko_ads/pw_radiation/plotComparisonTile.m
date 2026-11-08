function [legH, legLabels] = plotComparisonTile(ax, ads, feko, phiList, ...
    colors, freqScale, marker, markerSize, valueScale, axisLabelFontSize, markerStepHz)
%PLOTCOMPARISONTILE  Plot FEKO vs ADS rE voltage for one quantity and component.
%   [legH, legLabels] = PLOTCOMPARISONTILE(ax, ads, feko, phiList, colors,
%   freqScale, marker, markerSize, valueScale, axisLabelFontSize,
%   markerStepHz) draws one curve per entry of phiList onto axes ax: ADS as
%   a solid line, FEKO as markers only, in the same colour.
%
%   ads and feko are structs from importAdsCsv / importFekoDat. Columns are
%   located by phi value (findPhiColumn), so column order does not matter.
%   A phi missing from one source is skipped for that source with a warning.
%
%   valueScale multiplies the data (1e3 for magnitude V -> mV, 1 for phase);
%   freqScale multiplies the frequency axis (e.g. 1e-9 for GHz).
%   markerStepHz: FEKO marker step [Hz], staggered per phi (fekoMarkerIdx);
%   [] plots every sample. axisLabelFontSize is accepted but not used; the
%   caller styles the axes.
%
%   legH/legLabels are interleaved (ADS phi1, FEKO phi1, ADS phi2, ...).

    if nargin < 11, markerStepHz = []; end
    legH = gobjects(0);
    legLabels = strings(0);

    hold(ax, "on");

    for k = 1:numel(phiList)
        p   = phiList(k);
        col = colors(mod(k-1, size(colors,1)) + 1, :);

        % --- ADS: solid line ---------------------------------------
        if ~isempty(ads)
            ia = findPhiColumn(ads.phi, p);
            if isempty(ia)
                warning("plotComparisonTile:missingAdsPhi", ...
                    "No ADS column for phi = %g deg.", p);
            else
                h = plot(ax, ads.freq * freqScale, ads.data(:,ia) * valueScale, ...
                    "LineStyle","-", "LineWidth",1.8, "Color",col);
                legH(end+1)      = h;                                  %#ok<AGROW>
                legLabels(end+1) = "ADS, " + phiLatexLabel(p);          %#ok<AGROW>
            end
        end

        % --- FEKO: markers only ------------------------------------
        if ~isempty(feko)
            iF = findPhiColumn(feko.phi, p);
            if isempty(iF)
                warning("plotComparisonTile:missingFekoPhi", ...
                    "No FEKO column for phi = %g deg.", p);
            else
                m = fekoMarkerIdx(feko.freq, markerStepHz, k, numel(phiList));
                h = plot(ax, feko.freq(m) * freqScale, feko.data(m,iF) * valueScale, ...
                    "LineStyle","none", "Marker",marker, "MarkerSize",markerSize, ...
                    "MarkerEdgeColor",col, "MarkerFaceColor",col);
                legH(end+1)      = h;                                  %#ok<AGROW>
                legLabels(end+1) = "FEKO, " + phiLatexLabel(p);         %#ok<AGROW>
            end
        end
    end

    hold(ax, "off");
end