function [rEpt, theta, phi] = read_cst_farfield(filename, cst_voltage_unit)
%READ_CST_FARFIELD  Read one CST far-field export (one port, one frequency).
%   [rEpt, theta, phi] = READ_CST_FARFIELD(filename, cst_voltage_unit)
%
% Format: whitespace-delimited text, 2 header lines (column names, dashed
% separator), then columns Theta [deg], Phi [deg], Abs(E), Abs(Theta),
% Phase(Theta) [deg], Abs(Phi), Phase(Phi) [deg], Ax.Ratio; Abs(E) and
% Ax.Ratio are ignored. Despite the "[V/m]" header, the values are r*E (a
% voltage, in cst_voltage_unit, V by default) and are scaled to volts.
%
% CST's grid has no closing sample at phi = 360 deg (e.g. 46 theta x 60 phi
% points, phi = 0..354). The phi = 0 slice is appended at phi + 2*pi so that
% trapz integration covers a full period. The grid must be complete and
% regular; otherwise an error is raised.
%
% Returns:
%   rEpt  - (numth x numph+1 x 2) complex r*E (V); dim 3: 1 = phi, 2 = theta
%   theta - (numth x 1) polar angle (rad), 0..pi
%   phi   - (numph+1 x 1) azimuth (rad), including the closing sample
data = readmatrix(filename, 'NumHeaderLines', 2);

theta_raw = data(:,1);
phi_raw   = data(:,2);
tol_deg = 1;
if min(theta_raw) > tol_deg || max(theta_raw) < (180 - tol_deg)
    error('read_cst_farfield: "%s" theta spans %.3g to %.3g degrees, not the expected 0 to 180.', filename, min(theta_raw), max(theta_raw));
end
if min(phi_raw) > tol_deg
    error('read_cst_farfield: "%s" phi does not start near 0 degrees (got %.3g).', filename, min(phi_raw));
end

theta_col = theta_raw * pi/180;
phi_col   = phi_raw * pi/180;

v_scale = cst_unit_scale(cst_voltage_unit);
mag_Etheta = data(:,4) * v_scale;
ang_Etheta = data(:,5) * pi/180;
mag_Ephi   = data(:,6) * v_scale;
ang_Ephi   = data(:,7) * pi/180;

theta = sort(unique(theta_col));
phi   = sort(unique(phi_col));
numth = numel(theta); numph = numel(phi);

if numth*numph ~= size(data,1)
    error('read_cst_farfield: "%s" data (%d rows) is not a complete regular grid (expected %d = %d theta x %d phi points) -- check the export for missing/duplicate points.', ...
        filename, size(data,1), numth*numph, numth, numph);
end

[~, theta_idx] = ismember(theta_col, theta);
[~, phi_idx]   = ismember(phi_col, phi);
if any(theta_idx==0) || any(phi_idx==0)
    error('read_cst_farfield: "%s" has values that did not match the unique grid (unexpected floating-point mismatch) -- cannot index safely.', filename);
end
lin_idx = sub2ind([numth, numph], theta_idx, phi_idx);

Ephi_grid   = nan(numth, numph);
Etheta_grid = nan(numth, numph);
Ephi_grid(lin_idx)   = mag_Ephi   .* exp(1j*ang_Ephi);
Etheta_grid(lin_idx) = mag_Etheta .* exp(1j*ang_Etheta);

if any(isnan(Ephi_grid(:))) || any(isnan(Etheta_grid(:)))
    error('read_cst_farfield: "%s" grid has unfilled (NaN) entries after indexing -- the (theta,phi) grid may not be as regular as it appeared.', filename);
end

% Reconstruct the closing sample trapz needs: duplicate the phi(1) slice
% at phi(1)+2*pi (same physical direction, one full period later).
phi = [phi; phi(1) + 2*pi];
Ephi_grid   = [Ephi_grid,   Ephi_grid(:,1)];
Etheta_grid = [Etheta_grid, Etheta_grid(:,1)];

rEpt = zeros(numth, numph+1, 2);
rEpt(:,:,1) = Ephi_grid;
rEpt(:,:,2) = Etheta_grid;
end
