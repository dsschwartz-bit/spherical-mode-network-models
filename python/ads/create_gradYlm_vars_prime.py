r"""
create_gradYlm_vars_prime.py -- tesseral-harmonic gradients in the body frame.

Creates SphericalHarmonics_lib:GradYlm_vars_prime, which instantiates
Ftheta_Fphi_vars_prime as a subcircuit and adds Global-scope VARs
VAR_gradYlm_theta_prime / VAR_gradYlm_phi_prime. The formulas are those of
create_gradYlm_vars.py with Ftheta_l*/Fphi_l* -> Ftheta_prime_l*/Fphi_prime_l*
and phi -> phi_prime:

    Nl = 1/sqrt(l(l+1)),  PhiNorm = 1/sqrt(2*pi) (m = 0), 1/sqrt(pi) (m != 0)
    m = 0:  gradYlm_theta_prime =  Nl*PhiNorm*Ftheta_prime_l{l}_m0
            gradYlm_phi_prime   =  0
    m < 0:  gradYlm_theta_prime =  Nl*PhiNorm*Ftheta_prime_l{l}_m{|m|}*sin(|m|*phi_prime)
            gradYlm_phi_prime   =  Nl*PhiNorm*|m|*cos(|m|*phi_prime)*Fphi_prime_l{l}_m{|m|}
    m > 0:  gradYlm_theta_prime =  Nl*PhiNorm*Ftheta_prime_l{l}_m{m}*cos(m*phi_prime)
            gradYlm_phi_prime   = -Nl*PhiNorm*m*sin(m*phi_prime)*Fphi_prime_l{l}_m{m}

(theta_prime, phi_prime) is the body-frame direction that the antenna
rotation maps onto the observation direction; it is produced by
VectorFieldRotate_vars (create_vector_field_rotation_vars.py) in the parent
design. These gradients drive the Tx gains in create_vcvs_tx_gain.py. The
names differ from GradYlm_vars, so both cells can be used in one design.
l = 0 is skipped; lMax is read from FTheta_ADS_Discrete.txt.

Usage (ADS Python console, after start_console.py):
    run_script("create_gradYlm_vars_prime.py")

Requires SphericalHarmonics_lib:Ftheta_Fphi_vars_prime
(create_Fdata_prime_vars.py). The target cell is deleted and rebuilt if it exists.
Afterwards, check that both VARs have Global scope and create a symbol view.
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
FDATA_PRIME_CELL_NAME = "Ftheta_Fphi_vars_prime"   # existing cell (create_Fdata_prime_vars.py)
NEW_CELL_NAME = "GradYlm_vars_prime"               # new schematic created by this script

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
    """ADS-identifier-safe variable name, e.g. gradYlm_theta_prime_l1_mneg1."""
    m_part = f"mneg{abs(m)}" if m < 0 else f"m{m}"
    return f"{prefix}_l{l}_{m_part}"


def build_gradYlm_prime_formulas(l, m):
    """Returns (theta_prime_formula, phi_prime_formula) ADS expression
    strings for a given (l,m), referencing the body-frame library's
    Ftheta_prime_l{l}_m{|m|}/Fphi_prime_l{l}_m{|m|} (always positive m)
    and using phi_prime (not phi) in the trig terms."""
    mm = abs(m)
    nl_expr = f"1/sqrt({l*(l+1)})"
    phinorm_expr = "1/sqrt(2*pi)" if m == 0 else "1/sqrt(pi)"
    ftheta_prime_ref = f"Ftheta_prime_l{l}_m{mm}"
    fphi_prime_ref = f"Fphi_prime_l{l}_m{mm}"

    if m == 0:
        theta_prime_formula = f"{nl_expr}*{phinorm_expr}*{ftheta_prime_ref}"
        phi_prime_formula = "0"
    elif m < 0:
        theta_prime_formula = f"{nl_expr}*{phinorm_expr}*{ftheta_prime_ref}*sin({mm}*phi_prime)"
        phi_prime_formula = f"{nl_expr}*{phinorm_expr}*{mm}*cos({mm}*phi_prime)*{fphi_prime_ref}"
    else:
        theta_prime_formula = f"{nl_expr}*{phinorm_expr}*{ftheta_prime_ref}*cos({mm}*phi_prime)"
        phi_prime_formula = f"-{nl_expr}*{phinorm_expr}*{mm}*sin({mm}*phi_prime)*{fphi_prime_ref}"

    return theta_prime_formula, phi_prime_formula


# ---- Determine lMax from the same Wolfram-exported files ----
ftheta_data, ftheta_lmax_header = parse_dscrdata(FTHETA_FILE)
L_MAX = infer_lmax(ftheta_data, ftheta_lmax_header)
print(f"Using lMax = {L_MAX}.")

# ---- Get the active workspace and library (must already exist) ----
ws = de.active_workspace()
print(f"Active workspace: {ws}")

if LIBRARY_NAME not in ws.library_names:
    raise RuntimeError(
        f"Library '{LIBRARY_NAME}' not found -- run create_Fdata_vars.py / "
        f"create_Fdata_prime_vars.py first."
    )
lib = ws.open_library(LIBRARY_NAME)

# ---- Confirm the existing Ftheta_Fphi_vars_prime design exists (read-only) ----
fdata_prime_design = de.db_uu.open_design((LIBRARY_NAME, FDATA_PRIME_CELL_NAME, "schematic"))
print(f"Found existing design: {fdata_prime_design.design_name}")

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

var_theta_prime = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 0), name="VAR_gradYlm_theta_prime", angle=90)
assert var_theta_prime.is_var_instance

var_phi_prime = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 3), name="VAR_gradYlm_phi_prime", angle=90)
assert var_phi_prime.is_var_instance

with de.db_uu.Transaction(design) as transaction:
    # Instantiate Ftheta_Fphi_vars_prime so its Global-scope
    # Ftheta_prime_l*/Fphi_prime_l* variables are visible here.
    design.add_instance(fdata_prime_design.design_name, (0, 6), name="Fdata_prime_sub", angle=0)

    # theta_prime and phi_prime are not defined here; they come from the
    # parent design (nested scope), i.e. from VectorFieldRotate_vars.
    n_theta_prime = 0
    n_phi_prime = 0
    for l in range(1, L_MAX + 1):   # l=0 skipped: Nl undefined there
        for m in range(-l, l + 1):
            theta_prime_formula, phi_prime_formula = build_gradYlm_prime_formulas(l, m)
            var_theta_prime.vars[var_name("gradYlm_theta_prime", l, m)] = theta_prime_formula
            var_phi_prime.vars[var_name("gradYlm_phi_prime", l, m)] = phi_prime_formula
            n_theta_prime += 1
            n_phi_prime += 1

    del var_theta_prime.vars["X"]
    del var_phi_prime.vars["X"]

    # Global scope so the parent design can reference
    # gradYlm_theta_prime_*/gradYlm_phi_prime_* directly.
    var_theta_prime.set_global_scope()
    var_phi_prime.set_global_scope()

    transaction.commit()

design.save_design()
print(f"Created {n_theta_prime} gradYlm_theta_prime variables in "
      f"VAR_gradYlm_theta_prime and {n_phi_prime} gradYlm_phi_prime variables in "
      f"VAR_gradYlm_phi_prime, in {LIBRARY_NAME}:{NEW_CELL_NAME}")
print("Reminder: manually verify Global scope is set on both "
      "VAR_gradYlm_theta_prime and VAR_gradYlm_phi_prime via Component "
      "Options in the ADS GUI, and create a symbol view for "
      "GradYlm_vars_prime before placing it elsewhere. The parent design "
      "must supply theta_prime,phi_prime as local variables -- typically "
      "from instantiating VectorFieldRotate_vars alongside this cell.")

# ---- Lightweight spot check (string only, no evaluation needed) ----
sample_name_theta_prime = var_name("gradYlm_theta_prime", 1, -1)
sample_name_phi_prime = var_name("gradYlm_phi_prime", 1, -1)
print(f"Spot check: {sample_name_theta_prime} = {var_theta_prime.vars[sample_name_theta_prime]}")
print(f"Spot check: {sample_name_phi_prime} = {var_phi_prime.vars[sample_name_phi_prime]}")
