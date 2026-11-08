# CST / ADS validation

Scripts that compare CST Studio Suite simulations with the Keysight ADS circuit
built from the spherical-mode network model (and, optionally, the Simulink link).
Each script bootstraps the path with `setup.m` and adds its own helper folders.
Set the CST/ADS data folders in the `USER SETTINGS` block at the top of each script.
Generated files go to `output/cst_ads/<subfolder>`; Touchstone models are read
from `output/SnP`.

## farfield/
| Script | Purpose |
|---|---|
| `plot_cst_farfields.m` | 3D vector plot and pattern cuts of CST far-field exports |
| `compare_cst_ads_farfields.m` | CST vs ADS far-field patterns (polar cuts) at selected frequencies |
| `compare_S11_CST_ADS.m` | S11 with port 2 shorted/open: CST, ADS and the spherical-mode model |

The CST far-field readers and plotters live in `matlab/io/cst`.

## time_signals/
| Script | Purpose |
|---|---|
| `cst_to_ads_mdf.m` | CST port time signals -> ADS MDF excitation (`mdf_out`) |
| `cst_to_FEKO.m` | CST time signal -> FEKO CSV (`csv_out`) |
| `generate_user_excitations.m` | User-defined two-port excitations for CST, ADS and Simulink (`user_excitations`) |
| `filter_excitation_gaussian.m` | Causality analysis of the S-parameter model and an excitation, Gaussian filter (`filtered_excitations`) |
| `run_causality_port2_OC_2p5GHz.m` | One causality check example |
| `compare_rE_time_at_probe_CST_ADS.m` | CST vs ADS radiated r*E(t) at probe directions |
| `plot_rE_radial_CST_ADS.m` | Radial (spoke) plots of CST and ADS r*E(t) |
| `compare_radiated_fields.m` | CST / ADS / Simulink r*E(t) comparison (`plots`, `simulink_out`) |
| `plot_group_delay_radiated.m` | Group delay of the radiated field |

The CST signal readers/writers are in `matlab/io/cst` and the ADS MDF writer and CSV
readers in `matlab/io/ads`. `utils/` holds the validation helpers (causality workflow,
`plotCSTSignal`, `radialTracePlot`).

### time_signals/ook/
On-off-keying direct antenna modulation: `ook_cst_to_ads_mdf.m` (ADS MDF from
CST port V/I), `ook_ideal_vs_radiated_spectrum.m`, `ook_plot_rE_radial.m`
(quarter/half polar plots).

## Requirements
MATLAB with RF Toolbox (S-parameters, causality checks) and Signal Processing
Toolbox; Simulink and RF Blockset for the Simulink comparisons (`validation/simulink`).
