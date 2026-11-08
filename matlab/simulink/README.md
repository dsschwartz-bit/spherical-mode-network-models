# Simulink link models

Simulink counterpart of the ADS Tx-to-Rx spherical-mode link
(`python/ads/create_tx_rx_link_internal_ports_TRANSIENT.py`). The `.slx` models are not
stored in the repository; `build_tx_rx_link_simulink.m` generates them from the
spherical-mode Touchstone files (`output/SnP/S_dom_*.sNp`, produced by the workflows in
`matlab/workflows/`).

## Requirements

- MATLAB, Simulink
- RF Toolbox (`sparameters`, `rational`, `rationalfit`)
- RF Blockset (default, `'Platform','rf'`): the models used for the paper results
- Simscape and Simscape Electrical (only for the alternative `'Platform','simscape'`)

## Files

| File | Purpose |
| --- | --- |
| `build_tx_rx_link_simulink.m` | Builds the model (RF Blockset by default, or Simscape), saves `<model>.slx` and `<model>_vars.mat` |
| `link_update_gains.m` | Recomputes the gains and delay after changing r, theta, phi, rotation (no rebuild) |
| `simulinkRadiatedField.m` | Radiated r*E(t) of the Tx-only model at a list of directions, with a cache (builds and simulates the missing ones) |

Workflow scripts in `matlab/workflows/`:

- `run_tx_rx_link_simulink.m`: builds and simulates a Tx-Rx link and plots the Rx feed voltage.
- `run_tx_radiated_simulink.m`: radiated r*E(t) of the Tx-only model driven by a pulse (e.g. a
  Gaussian-modulated sine from a CST port export).
- `ook/run_ook_simulink_rf.m`: on-off keying (direct antenna modulation) with a switch at the
  internal port; helpers `ook_load_excitation.m` and `ook_observation_loop.m` are in the same folder.

The comparisons with CST and ADS are in `validation/simulink/` and
`validation/cst_ads/time_signals/compare_radiated_fields.m`.

## Build and run

```matlab
run('<repo>/setup.m');                       % or let the scripts do it
snp = fullfile(repo_root, 'output', 'SnP', 'S_dom_Q_dipole_cylinder_edgeport_full.s6p');
mdl = build_tx_rx_link_simulink(snp, snp, 'r', 10, 'Waveform', 'sine', 'SourceFreq', 1e9);
out = sim(mdl);
v = out.yout.getElement('FeedRx').Values;  plot(v.Time*1e9, v.Data)
```

Or edit the settings at the top of `run_tx_rx_link_simulink.m` and run it.

- Models and their `_vars.mat` go to `output/simulink/` by default (`'OutputDir'` changes
  it). The model's PreLoadFcn loads `<model>_vars.mat` from the model's own folder, so keep
  the two files together.
- To change the geometry without rebuilding: set `r`, `theta`, `phi`, `rot_angle`, ... in the
  base workspace, run `link_update_gains`, then `sim(mdl)`.
- The spherical-harmonic gradients use the ADS tables in `mathematica/tables`
  (`FTheta_ADS_Discrete.txt`, `FPhi_ADS_Discrete.txt`); set the workspace variable
  `Ylm_data_dir` to use another folder. Without the tables `gradYlm.m` is used.
- `'TxOnly', true` builds only the radiating side (Vtheta, Vphi = r*E_theta, r*E_phi), used by
  the CST / ADS comparisons.
