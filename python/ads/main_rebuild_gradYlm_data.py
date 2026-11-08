r"""
main_rebuild_gradYlm_data.py -- build all Legendre/gradYlm VAR cells in order.

Runs, in dependency order:
    1. create_Fdata_vars.py         -> SphericalHarmonics_lib:Ftheta_Fphi_vars
    2. create_Fdata_prime_vars.py   -> SphericalHarmonics_lib:Ftheta_Fphi_vars_prime
    3. create_gradYlm_vars.py       -> SphericalHarmonics_lib:GradYlm_vars       (needs 1)
    4. create_gradYlm_vars_prime.py -> SphericalHarmonics_lib:GradYlm_vars_prime (needs 2)

Use after regenerating FTheta_ADS_Discrete.txt / FPhi_ADS_Discrete.txt with
mathematica/LegendreFThetaFPhi.wl. Existing target cells are deleted and
rebuilt.

Usage (ADS Python console, after start_console.py):
    run_script("main_rebuild_gradYlm_data.py")
"""

from pathlib import Path

# Repository folders (defined by start_console.py; derived here if this file is exec()'d directly).
try:
    SCRIPT_DIR
except NameError:
    SCRIPT_DIR = Path(THIS_SCRIPT_PATH).resolve().parent

SCRIPTS_IN_ORDER = [
    "create_Fdata_vars.py",
    "create_Fdata_prime_vars.py",
    "create_gradYlm_vars.py",
    "create_gradYlm_vars_prime.py",
]

for script_name in SCRIPTS_IN_ORDER:
    print(f"\n===== Running {script_name} =====")
    exec(open(SCRIPT_DIR / script_name).read())

print("\n===== All four scripts completed. =====")
