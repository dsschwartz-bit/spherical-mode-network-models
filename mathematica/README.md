# Mathematica equation generators

These Wolfram Language scripts build the closed-form expressions the toolbox
needs and export them as text tables. Each script sets its working directory
to this folder (`SetDirectory[NotebookDirectory[]]`) and writes its output to
`tables/`. The generated tables are committed, so you only need to re-run a
script after changing it. Open the `.wl` file in Mathematica and evaluate it.

`Legendre_fullnorm_script.nb` is the interactive notebook version of the
Legendre derivation.

| Script | Generates (`tables/`) | Content | Used by |
| --- | --- | --- | --- |
| `LegendreFThetaFPhi.wl` | `FTheta_table.txt`, `FPhi_table.txt` | l × m tab-delimited tables of MATLAB `str2sym` strings for F_θ = −sin θ dP̄ₗₘ/dx and F_φ = P̄ₗₘ/sin θ (variable `th`, l = 0..11) | `matlab/sphmodes/LoadLegendreFHandles.m`, `matlab/sphmodes/ImportAndDisplayLegendreTables.m` (and through them `EvalLegendreF.m`, `tests/testing_gradYlm*.m`) |
| | `FTheta_ADS_Discrete.txt`, `FPhi_ADS_Discrete.txt` | ADS DataAccessComponent "Discrete" files, index i = l(l+1)+m, ADS-syntax formula in `theta` | `python/ads/create_Fdata_vars.py`, `create_Fdata_prime_vars.py`, `create_gradYlm_vars.py`, `create_gradYlm_vars_prime.py` (rebuilt together by `main_rebuild_gradYlm_data.py`); `matlab/simulink/link_update_gains.m` |
| | `FTheta_table.tex`, `FPhi_table.tex` | LaTeX longtables | documentation only |
| `rotation3D.wl` | `Rotation_ADS_Equations.txt` | ADS Eqn text for the rotation matrix R(rot_angle, Ux, Uy, Uz) (Rodrigues), with `magU` and `R11..R33` | `python/ads/create_rotationmat_vars.py` |
| | `Rotation_table.txt` | MATLAB `str2sym` strings of R | not read by the repository code |
| | `Rotation_table.tex` | LaTeX matrix | documentation only |
| `vector_field_spherical_coordinate_transform.wl` | `Cart2Sph_ADS_Equations.txt` | ADS Eqn text: `Cart2Sph_ij(theta, phi)` and the generic Cartesian template `M_ij`, `r_of_xyz`, `rho_of_xyz`, `theta_of_xyz`, `phi_of_xyz`, used to rotate a field as E′(r) = R·E(R⁻¹r) | `python/ads/create_vector_field_rotation_vars.py` |
| `spherical_coordinate_transform.wl` | `Sph2Cart_Cart2Sph_ADS_Equations.txt` | ADS Eqn text: `Sph2Cart_ij(theta, phi)` and `Cart2Sph_ij` at the rotated sample direction `x_rot, y_rot, z_rot` | not read by the repository code (alternative formulation) |

The Python ADS scripts locate the tables through `DATA_DIR = <repo>/mathematica/tables`,
which `python/ads/start_console.py` sets. The MATLAB code uses
`fullfile(repo_root, 'mathematica', 'tables')`.
