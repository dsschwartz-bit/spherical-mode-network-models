function plotTransmissionPolarTile(pax, dataset, colors, freqScale, ...
    marker, markerSize, axisLabelFontSize)
%PLOTTRANSMISSIONPOLARTILE  Plot transmission phase as a polar curve (radius = frequency).
%   PLOTTRANSMISSIONPOLARTILE(pax, dataset, colors, freqScale, marker,
%   markerSize, axisLabelFontSize) draws, onto the polar axes pax, one
%   curve per element of the dataset struct array: angular position is
%   the transmission phase (deg, converted to radians), and radius is
%   frequency. This avoids the +/-180 deg wraparound of a Cartesian
%   phase plot and makes it easy to compare several traces at once.
%
%   ADS-model data is drawn as a dense solid line; FEKO data is
%   drawn as downsampled markers only, in the same color as the
%   matching ADS curve.
%
%   dataset(k) fields: see plotTransmissionTile / buildTransmissionDataset.

    hold(pax,"on");

    for k = 1:numel(dataset)
        d = dataset(k);
        col = colors(k,:);

        thetaAds = deg2rad(d.ads.phase);
        rhoAds = d.freq * freqScale;
        polarplot(pax, thetaAds, rhoAds, "LineWidth", 1.8, "Color", col, ...
            "LineStyle","-");

        thetaFeko = deg2rad(d.fekoDS.phase.data);
        rhoFeko = d.fekoDS.phase.freq * freqScale;
        polarplot(pax, thetaFeko, rhoFeko, "LineStyle","none", ...
            "Marker", marker, "MarkerSize", markerSize, ...
            "MarkerEdgeColor", col, "MarkerFaceColor", col);
    end

    hold(pax,"off");
    grid(pax,"on");
    set(pax, "FontSize", axisLabelFontSize);
end
