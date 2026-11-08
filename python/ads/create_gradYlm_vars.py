r"""
create_gradYlm_vars.py -- normalized gradient of each tesseral harmonic in ADS.

Creates SphericalHarmonics_lib:GradYlm_vars, which instantiates
Ftheta_Fphi_vars as a subcircuit and adds Global-scope VARs
VAR_gradYlm_theta / VAR_gradYlm_phi holding, for l = 1..lMax, |m| <= l:

    Nl = 1/sqrt(l(l+1)),  PhiNorm = 1/sqrt(2*pi) (m = 0), 1/sqrt(pi) (m != 0)
    m = 0:  gradYlm_theta =  Nl*PhiNorm*Ftheta_l{l}_m0
            gradYlm_phi   =  0
    m < 0:  gradYlm_theta =  Nl*PhiNorm*Ftheta_l{l}_m{|m|}*sin(|m|*phi)
            gradYlm_phi   =  Nl*PhiNorm*|m|*cos(|m|*phi)*Fphi_l{l}_m{|m|}
    m > 0:  gradYlm_theta =  Nl*PhiNorm*Ftheta_l{l}_m{m}*cos(m*phi)
            gradYlm_phi   = -Nl*PhiNorm*m*sin(m*phi)*Fphi_l{l}_m{m}

The F tables are always referenced with positive m. l = 0 is skipped (Nl is
undefined and the monopole term is not radiated). This matches the
analytic gradient in tests/testing_gradYlm.m. theta and phi are supplied
by the parent design. lMax is read from FTheta_ADS_Discrete.txt.

Usage (ADS Python console, after start_console.py):
    run_script("create_gradYlm_vars.py")

Requires SphericalHarmonics_lib:Ftheta_Fphi_vars (create_Fdata_vars.py).
The target cell is deleted and rebuilt if it exists. Afterwards, check that both VARs
have Global scope and create a symbol view before placing the cell.
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
FTHETA_FILE = DATA_DIR / "FTheta_ADS_Discrete.txt"
FPHI_FILE = DATA_DIR / "FPhi_ADS_Discrete.txt"

LIBRARY_NAME = "SphericalHarmonics_lib"
FDATA_CELL_NAME = "Ftheta_Fphi_vars"   # existing cell with Ftheta_l*/Fphi_l*
NEW_CELL_NAME = "GradYlm_vars"         # new schematic created by this script

# ============================================================================


def parse_dscrdata(filepath):
    """Parse a BEGIN DSCRDATA / END DSCRDATA file into {index: formula_str},
    also returning lMax if a 'REM lMax=N' comment line is present (None
    otherwise, in which case the caller should infer it -- see
    infer_lmax)."""
    data = {}
    lmax_header = None
    in_block = False
    with open(filepath, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.upper().startswith("REM"):
                match = re.search(r"lMax\s*=\s*(\d+)", line, re.IGNORECASE)
                if match:
                    lmax_header = int(match.group(1))
                continue
            if line.upper().startswith("BEGIN"):
                in_block = True
                continue
            if line.upper().startswith("END"):
                in_block = False
                continue
            if line.startswith("%"):
                continue
            if in_block:
                parts = line.split()
                if len(parts) == 2:
                    idx_str, formula = parts
                    data[int(idx_str)] = formula
    return data, lmax_header


def infer_lmax(data, header_lmax):
    """Use the REM-header lMax if present; otherwise derive it from the
    largest index in the data (i_max = lMax(lMax+2))."""
    if header_lmax is not None:
        return header_lmax
    if not data:
        raise ValueError("No data parsed; cannot infer lMax.")
    i_max = max(data.keys())
    l_max = round(-1 + (1 + i_max) ** 0.5)
    if l_max * (l_max + 2) != i_max:
        raise ValueError(
            f"Could not cleanly infer lMax from max index {i_max} -- "
            f"the data may not cover a full 0..lMax range."
        )
    return l_max


def var_name(prefix, l, m):
    """ADS-identifier-safe variable name, e.g. gradYlm_theta_l1_mneg1."""
    m_part = f"mneg{abs(m)}" if m < 0 else f"m{m}"
    return f"{prefix}_l{l}_{m_part}"


def build_gradYlm_formulas(l, m):
    """Returns (theta_formula, phi_formula) ADS expression strings for a
    given (l,m), referencing the library's Ftheta_l{l}_m{|m|}/
    Fphi_l{l}_m{|m|} (always positive m)."""
    mm = abs(m)
    nl_expr = f"1/sqrt({l*(l+1)})"
    phinorm_expr = "1/sqrt(2*pi)" if m == 0 else "1/sqrt(pi)"
    ftheta_ref = f"Ftheta_l{l}_m{mm}"
    fphi_ref = f"Fphi_l{l}_m{mm}"

    if m == 0:
        theta_formula = f"{nl_expr}*{phinorm_expr}*{ftheta_ref}"
        phi_formula = "0"
    elif m < 0:
        theta_formula = f"{nl_expr}*{phinorm_expr}*{ftheta_ref}*sin({mm}*phi)"
        phi_formula = f"{nl_expr}*{phinorm_expr}*{mm}*cos({mm}*phi)*{fphi_ref}"
    else:
        theta_formula = f"{nl_expr}*{phinorm_expr}*{ftheta_ref}*cos({mm}*phi)"
        phi_formula = f"-{nl_expr}*{phinorm_expr}*{mm}*sin({mm}*phi)*{fphi_ref}"

    return theta_formula, phi_formula


# ---- Determine lMax from the same Wolfram-exported files ----
ftheta_data, ftheta_lmax_header = parse_dscrdata(FTHETA_FILE)
L_MAX = infer_lmax(ftheta_data, ftheta_lmax_header)
print(f"Using lMax = {L_MAX}.")

# ---- Get the active workspace and library (must already exist) ----
ws = de.active_workspace()
print(f"Active workspace: {ws}")

if LIBRARY_NAME not in ws.library_names:
    raise RuntimeError(
        f"Library '{LIBRARY_NAME}' not found -- run create_Fdata_vars.py first."
    )
lib = ws.open_library(LIBRARY_NAME)

# ---- Confirm the existing Ftheta_Fphi_vars design exists (read-only) ----
fdata_design = de.db_uu.open_design((LIBRARY_NAME, FDATA_CELL_NAME, "schematic"))
print(f"Found existing design: {fdata_design.design_name}")

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

var_theta = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 0), name="VAR_gradYlm_theta", angle=90)
assert var_theta.is_var_instance

var_phi = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 3), name="VAR_gradYlm_phi", angle=90)
assert var_phi.is_var_instance

with de.db_uu.Transaction(design) as transaction:
    # Instantiate Ftheta_Fphi_vars as a subcircuit so its global-scope
    # Ftheta_l*/Fphi_l* variables become visible here.
    design.add_instance(fdata_design.design_name, (0, 6), name="Fdata_sub", angle=0)

    # theta and phi are not defined here; they come from the parent design
    # (nested scope).
    n_theta = 0
    n_phi = 0
    for l in range(1, L_MAX + 1):   # l=0 skipped: Nl undefined there
        for m in range(-l, l + 1):
            theta_formula, phi_formula = build_gradYlm_formulas(l, m)
            var_theta.vars[var_name("gradYlm_theta", l, m)] = theta_formula
            var_phi.vars[var_name("gradYlm_phi", l, m)] = phi_formula
            n_theta += 1
            n_phi += 1

    del var_theta.vars["X"]
    del var_phi.vars["X"]

    # Global scope so the parent design can reference
    # gradYlm_theta_*/gradYlm_phi_* directly.
    var_theta.set_global_scope()
    var_phi.set_global_scope()

    transaction.commit()

design.save_design()
print(f"Created {n_theta} gradYlm_theta variables in VAR_gradYlm_theta and "
      f"{n_phi} gradYlm_phi variables in VAR_gradYlm_phi, in "
      f"{LIBRARY_NAME}:{NEW_CELL_NAME}")
print("Reminder: manually verify Global scope is set on both "
      "VAR_gradYlm_theta and VAR_gradYlm_phi via Component Options in "
      "the ADS GUI, and create a symbol view for GradYlm_vars before "
      "placing it elsewhere.")

# ---- Lightweight spot check (string only, no evaluation needed) ----
sample_name_theta = var_name("gradYlm_theta", 1, -1)
sample_name_phi = var_name("gradYlm_phi", 1, -1)
print(f"Spot check: {sample_name_theta} = {var_theta.vars[sample_name_theta]}")
print(f"Spot check: {sample_name_phi} = {var_phi.vars[sample_name_phi]}")
