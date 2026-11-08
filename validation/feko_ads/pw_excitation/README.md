# FEKO vs ADS: load current under plane-wave excitation

Compares the current induced in the load of a receiving antenna illuminated by a
plane wave, computed by the ADS spherical-mode circuit model and by FEKO.

Requires MATLAB R2020b or newer. Each script runs `setup.m` from the repository
root and adds `../common` itself. Results are shown as docked MATLAB figures.

## Scripts

### `plot_current_comparison.m` — bowtie, ADS vs FEKO

Bowtie antenna with a 45° edge port, plane wave incident from several directions
(phi, theta), e.g. `ph30th60`, `ph60th60`, `ph90th60`.

Inputs (in `dataFolder`, set at the top of the script):

| Source | Files | Content |
| --- | --- | --- |
| ADS  | one `.csv` per direction, tag in the file name (`..._ph30th60_...csv`) | `freq [Hz], mag(I_Probe1.i) [A], phase(I_Probe1.i) [deg]` below a quoted ADS header |
| FEKO | `..._FEKO_loadcurrents_magnitude.dat`, `..._FEKO_loadcurrents_phase.dat` | tab-delimited, one column per direction headed `"Load_<tag> [A]"` / `"Load_<tag> [deg]"` |

ADS and FEKO curves are paired by the `ph<phi>th<theta>` tag.

Outputs:

- printed tables of the mean and largest complex-current error
  `|I_ADS - I_FEKO|` (uA) per direction, ADS spline-interpolated onto the FEKO
  frequencies, from the first common frequency up to 1.8 GHz and up to the end
  of the plotted band (`errUpperHz`);
- figures `load_current_magnitude`, `load_current_phase`, `load_current_phasor`
  (polar plot of I over frequency) and a separate `load_current_legend`.
  ADS is drawn as solid lines, FEKO as thinned filled markers in the same colour.

### `plot_current_comparison_2.m` — dipole, Z-aligned vs rotated (ADS only)

Input: one ADS `.csv` (`adsCsvFile` at the top of the script) with columns
`freq, mag_rotated, phase_rotated, mag, phase`, read by `import_ADS_csv.m`.
Output: figures `rotated_load_current_magnitude|phase|phasor` and
`rotated_load_current_legend`.

## Helper functions

- `importCircuitModelSet.m`, `importCircuitModelCsv.m` — ADS load-current CSVs → struct array (`tag, freq, mag, phase`)
- `importFekoLoadCurrents.m` — FEKO load-current `.dat` → `freq`, `data` (one column per load), `tags`
- `import_ADS_csv.m` — script, reads the dipole CSV used by `plot_current_comparison_2.m`
- `plotCurrentTile.m` — magnitude or phase vs frequency (ADS line + FEKO markers)
- `plotCurrentPolarTile.m` — current phasor on polar axes
- `tagToLatexLabel.m` — legend label from a `ph..th..` tag

Shared plotting style and helpers (`plotStyle`, `newFigure`, `legendFigure`, `fekoMarkerIdx`, ...)
are in `../common/`.

Phase is plotted as exported; neither source is wrapped, unwrapped or offset.
