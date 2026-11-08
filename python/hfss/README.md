# HFSS far-field export

`ExportFarFields.py` exports the far fields that
`matlab/workflows/general_vector_spherical_modes_loadports_HFSS.m` needs:
rE_phi and rE_theta (magnitude and phase) at every theta, phi and frequency,
one CSV per port, read by `matlab/io/hfss/importffarfield_HFSS.m`.

## Before running

- A discrete frequency sweep with a uniform step, solved with **Save Fields**.
- A far-field (Infinite Sphere) setup covering theta 0–180°, phi 0–360°.
- The antenna's Touchstone file exported from the same design (it carries
  the apothem `a` in its "Variables" comments).

## Running

In HFSS: **Tools > Run Script**, choose `ExportFarFields.py` and enter the
Script Arguments (comma separated; blank fields take the defaults):

```
output_file_name,field_type,project,design,setup,sweep,FF_setup,sources
```

`,,,,,,,` uses the active project and design, `Setup1`, `Sweep` and
`Infinite Sphere1`, and writes one file per excitation:
`<project>_<design>_farfield_rE_<port>_<mode>.csv`, next to the script.
The header of the script lists every argument.

In the USER SETTINGS of the HFSS workflow, set `filedirectory` (Touchstone
file), `farfield_folder` (CSV files; defaults to `filedirectory`),
`antenna_name` and `farfield_basename` (the output_file_name above).

The script uses the AEDT scripting interface (`ScriptEnv`, IronPython) that
HFSS runs scripts with; it is not meant to be run from a separate Python.
