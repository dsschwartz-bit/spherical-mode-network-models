r"""
create_rotationmat_vars.py -- 3D rotation matrix R(rot_angle, Ux, Uy, Uz) in ADS.

Reads mathematica/tables/Rotation_ADS_Equations.txt (written by
BuildADSEquationFile[] in rotation3D.wl: "//" comment lines followed by
"name=formula;" lines) and creates Rotation_lib:RotationMatrix_vars with one
Global-scope VAR, VAR_RotationMatrix, holding magU = sqrt(Ux^2+Uy^2+Uz^2)
and the nine entries R11..R33 of the rotation by rot_angle about the axis
(Ux, Uy, Uz). rot_angle, Ux, Uy, Uz are supplied by the parent design.

Usage (ADS Python console, after start_console.py):
    run_script("create_rotationmat_vars.py")

Requires the library Rotation_lib (create it in the ADS GUI first). The
target cell is deleted and rebuilt if it exists. Afterwards, check that the VAR has
Global scope and create a symbol view before placing the cell.
"""

from pathlib import Path
import re
import keysight.ads.de as de

# ==================== Configuration -- edit as needed ====================

# Repository folders (defined by start_console.py; derived here if this file is exec()'d directly).
try:
    SCRIPT_DIR
except NameError:
    SCRIPT_DIR = Path(THIS_SCRIPT_PATH).resolve().parent
DATA_DIR = SCRIPT_DIR.parents[1] / "mathematica" / "tables"
ROTATION_FILE = DATA_DIR / "Rotation_ADS_Equations.txt"

LIBRARY_NAME = "Rotation_lib"
NEW_CELL_NAME = "RotationMatrix_vars"   # new schematic created by this script

# ============================================================================


def parse_ads_equations(filepath):
    """Parse a flat "name=formula;" equation file (BuildADSEquationFile[]
    in rotation3D.wl) into an ordered {name: formula} dict. Blank lines
    and lines starting with "//" are skipped; the trailing ';' is
    optional."""
    entries = {}
    with open(filepath, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("//"):
                continue
            if line.endswith(";"):
                line = line[:-1]
            if "=" not in line:
                raise ValueError(f"Could not parse equation line: {raw_line!r}")
            name, formula = line.split("=", 1)
            entries[name.strip()] = formula.strip()
    return entries


# ---- Read magU plus the 9 R11..R33 formulas straight from the Wolfram export ----
rotation_vars = parse_ads_equations(ROTATION_FILE)
print(f"Parsed {len(rotation_vars)} entries from {ROTATION_FILE.name}.")

expected_names = ["magU"] + [f"R{i}{j}" for i in (1, 2, 3) for j in (1, 2, 3)]
missing = [n for n in expected_names if n not in rotation_vars]
if missing:
    raise ValueError(f"Missing expected entries in {ROTATION_FILE.name}: {missing}")

# ---- Get the active workspace and library (must already exist) ----
ws = de.active_workspace()
print(f"Active workspace: {ws}")

if LIBRARY_NAME not in ws.library_names:
    raise RuntimeError(
        f"Library '{LIBRARY_NAME}' not found -- create it in the ADS GUI first."
    )
lib = ws.open_library(LIBRARY_NAME)

# ---- Create the new schematic ----
# Rebuild: delete the cell first if it exists (create_schematic cannot
# overwrite a view, and delete_cell needs the cell's views closed).
if lib.cell_exists(NEW_CELL_NAME):
    print(f"Cell '{NEW_CELL_NAME}' already exists in '{LIBRARY_NAME}' -- deleting it first.")
    for view_name in ("schematic", "symbol"):
        try:
            de.db_uu.open_design((LIBRARY_NAME, NEW_CELL_NAME, view_name)).close_design()
        except Exception:
            pass
    lib.delete_cell(NEW_CELL_NAME)
design = de.db_uu.create_schematic((LIBRARY_NAME, NEW_CELL_NAME, "schematic"))

var_rotation = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 0), name="VAR_RotationMatrix", angle=90)
assert var_rotation.is_var_instance

with de.db_uu.Transaction(design) as transaction:
    # rot_angle, Ux, Uy, Uz are not defined here; they come from the parent
    # design (nested scope).
    n_vars = 0
    for name in expected_names:
        var_rotation.vars[name] = rotation_vars[name]
        n_vars += 1

    del var_rotation.vars["X"]

    # Global scope so the parent design can reference R11..R33 directly.
    var_rotation.set_global_scope()

    transaction.commit()

design.save_design()
print(f"Created {n_vars} variables (magU + 9 rotation matrix entries) in "
      f"VAR_RotationMatrix, in {LIBRARY_NAME}:{NEW_CELL_NAME}")
print("Reminder: manually verify Global scope is set on VAR_RotationMatrix "
      "via Component Options in the ADS GUI, and create a symbol view for "
      "RotationMatrix_vars before placing it elsewhere. The parent design "
      "must supply rot_angle, Ux, Uy, Uz as local variables.")

# ---- Lightweight spot check (string only, no evaluation needed) ----
print(f"Spot check: magU = {var_rotation.vars['magU']}")
sample_name = "R12"
print(f"Spot check: {sample_name} = {var_rotation.vars[sample_name]}")