# FEKO vs ADS circuit-model validation

These scripts compare the ADS circuit models built from spherical-mode network
models with FEKO full-wave results and plot the comparisons. The
simulation data are not included in the repository. Export them from FEKO and
ADS, then set the data paths in the `USER SETTINGS` block at the top of each script.

Each script runs `setup.m` from the repository root, so `matlab/` and
`third_party/` are on the path. It also adds its own folder and `common/` (shared
plot style and helpers) to the path. Results are shown as docked MATLAB figures;
error metrics are printed to the console.

| Example | Script(s) | Inputs | Outputs |
| --- | --- | --- | --- |
| `pw_excitation/` — load current under plane-wave excitation | `plot_current_comparison.m` (bowtie, ADS vs FEKO), `plot_current_comparison_2.m` (dipole, Z-aligned vs rotated, ADS only) | ADS `.csv` load-current exports (one per incidence direction); FEKO load-current `.dat` exports (magnitude, phase) | \|I\|, phase and phasor figures; current-error tables |
| `pw_radiation/` — rE far-field voltage (Vtheta, Vphi) of a loop | `plot_rE_voltage_comparison.m`, `loop_not_rotated/plot_loop_rE_voltage_comparison.m`, `check_import_integrity.m` | ADS `.csv` rE exports (long format, phi cuts stacked); FEKO far-field `.dat` exports (one column per phi cut) | magnitude, phase and phasor figures per component; voltage-error tables |
| `time_domain/` — Tx-loop / Rx-bowtie transient | `plot_impulse_comparison.m`, `plot_RFpulse_comparison.m` | ADS data-display `.csv` exports (`Tx_loop_Rx_bowtie_ts.csv`, `Tx_loop_Rx_bowtie_RFpulse.csv`). The FEKO result enters as an `.s2p` block in the ADS schematic | impulse / RF-pulse response figures, shared legend; relative MSE vs FEKO |
| `txrx/` — Tx-Rx link | `plot_transmission_comparison.m` | `.s4p` files that hold the ADS link (ports 1-2) and the FEKO link (ports 3-4) for each Rx direction | \|S21\| and phase figures with an \|S21\| error summary |
| `common/` | — | — | shared helpers: plot style (`plotStyle.m`, `newFigure.m`, `styleAxes.m`, `legendFigure.m`, ...), `downsampleFeko.m`, `fekoMarkerIdx.m`, `thetaTagToLatexLabel.m` |

See `pw_excitation/README.md` and `pw_radiation/README.md` for file formats and
settings. The `txrx` script needs the RF Toolbox.
