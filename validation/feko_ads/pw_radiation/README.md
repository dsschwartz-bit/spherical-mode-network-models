# FEKO vs ADS: rE far-field voltage (Vtheta / Vphi)

Compares the far-field voltage rE (components Vtheta and Vphi) radiated by a loop
antenna without loads, computed by the ADS spherical-mode circuit model and by
FEKO, on three phi cuts (30°, 60°, 90°).

Requires MATLAB R2020b or newer (`arguments`, `readlines`). Each script runs
`setup.m` from the repository root and adds `../common` itself.

## Scripts

- `plot_rE_voltage_comparison.m` — main script. `loopCase` selects the geometry:
  - `"not_rotated"`: loop axis along z, files `loop_noloads_*`, read from
    `loop_not_rotated/`; ADS sweep variable `phi_deg` in degrees.
  - `"rotated"`: loop rotated 90° about x, files `looprotated_noloads_*`, read from
    this folder; ADS sweep variable `phi_prime` in radians.
- `loop_not_rotated/plot_loop_rE_voltage_comparison.m` — the not-rotated case as a
  standalone script with its own marker settings and `loop_rE_voltage_*` figure names.
- `check_import_integrity.m` — prints the first raw data line of each rotated-loop
  export next to the imported values and the column → phi mapping.

## Inputs

One file per quantity (`mag`, `phase`) and component (`Vtheta`, `Vphi`):

| Source | File | Format |
| --- | --- | --- |
| ADS  | `<prefix>_ADS_rE_<q>_<c>.csv`  | quoted multi-line header, then `phi, freq, value` rows; the phi cuts are stacked (long format) |
| FEKO | `<prefix>_FEKO_rE_<q>_<c>.dat` | `Far field` line, tab-delimited, one column per cut headed `"Phi = 60 deg [V]"` (column order varies between files) |

Curves are paired by phi value (`findPhiColumn`, 0.25° tolerance), never by
column position.

## Outputs

- Printed tables of the mean and largest complex-voltage error
  `|E_ADS - E_FEKO|` (mV) per component and cut, ADS spline-interpolated onto
  the FEKO frequencies, up to 1.8 GHz and up to `freqMaxHz`.
- Docked figures:
  `<name>_<comp>_magnitude`, `_phase`, `_phasor` (polar plot of the voltage
  phasor over frequency) and a shared `<name>_legend`. ADS is a solid line,
  FEKO markers in the same colour, one colour per cut.

## Helper functions

- `importFieldData.m` — detects the format from the file content and calls the matching importer
- `importAdsCsv.m` — ADS CSV → wide matrix, one column per phi
- `importFekoDat.m` — FEKO `.dat` → matrix, one column per phi (file order, labelled)
- `normalizePhi.m` — converts the sweep variable to degrees; optional automatic rescaling
- `checkPhiCoverage.m` — warns if a requested cut is missing
- `findPhiColumn.m`, `restrictFreqRange.m`, `phiLatexLabel.m`, `componentLatex.m`
- `plotComparisonTile.m` — magnitude or phase vs frequency
- `plotComparisonPolarTile.m` — voltage phasor on polar axes

`importAdsCsv` and `importFekoDat` return the same struct
(`source, quantity, component, unit, phi, header, freq, data, file`).

## Settings

- `phiList = [30 60 90]` — cuts to plot.
- `restrictToCommonBand = true` — clip both sources to their common band (and
  `freqMaxHz`), so the polar traces cover the same sweep.
- `adsPhiUnit` / `fekoPhiUnit` — unit of the phi sweep variable in each source.
- `phiScaleMode = "auto"` — if a file's phi values do not match `phiList`
  (e.g. 1719 / 3438 / 5157 for a radian sweep converted twice), try the factors
  pi/180, 180/pi, 1e-3, 1e3 and apply the one that reproduces `phiList`; `1` disables this.

Phase is plotted as exported; neither source is wrapped, unwrapped or offset. The
phase axis is fixed to [-180, 180] only when the data already lie in that range.

## Troubleshooting

- **`No column header containing "Phi = ..." found`** — the file named as a FEKO
  export contains ADS-format data. Re-export it from FEKO.
- **`no column for phi = [...] deg`** — the requested cut is not in the file after
  the unit conversion; check `adsPhiUnit` / `phiScaleMode`.
- Before attributing a mismatch to the model, check the number of digits in the
  ADS CSV export (few decimals quantize small values) and that both tools use the
  same phase reference and time convention.
