# Simulink validation

Compares the Simulink link models (built by `matlab/simulink/build_tx_rx_link_simulink.m`)
with CST and ADS. Requires Simulink, RF Toolbox, RF Blockset and Signal Processing Toolbox
(`hilbert`, `xcorr`, `hann`).

| File | Purpose |
| --- | --- |
| `compare_simulink_ADS_CST.m` | Radiated r*E(t) of the Tx-only model compared with CST and ADS |
| `plot_rE_radial_simulink.m` | Radial (spoke) plots of Simulink / CST / ADS far-field signals |
| `ook/` | On-off keying (direct antenna modulation) with a switch at the internal port |

The scripts need your own CST exports and ADS CSV exports: set their folders in the
`USER SETTINGS` block at the top of each script. Plots, caches and figures are written to
`output/simulink/` (and `output/simulink/ook/`).

## OOK (`ook/`)

- `compare_ook_CST_ADS_simulink.m` builds the RF Blockset Tx-only model with a voltage-controlled
  switch at the internal port (as `matlab/workflows/ook/run_ook_simulink_rf.m`), simulates it and
  compares it with CST and ADS (time signals, spectra, radial cuts).
- The user-defined excitations come from `validation/cst_ads/time_signals/generate_user_excitations.m`.
