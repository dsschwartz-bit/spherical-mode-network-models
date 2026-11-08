r"""
create_Fdata_prime_vars.py -- body-frame copy of the F_theta/F_phi functions.

Reads the same FTheta_ADS_Discrete.txt / FPhi_ADS_Discrete.txt tables as
create_Fdata_vars.py and creates SphericalHarmonics_lib:Ftheta_Fphi_vars_prime
with Global-scope VARs VAR_Ftheta_prime / VAR_Fphi_prime holding
Ftheta_prime_l{l}_m{m} / Fphi_prime_l{l}_m{m}, where every whole-word
"theta" -> "theta_prime" and "phi" -> "phi_prime".

These are used by the vector-field-rotation chain (see
create_vector_field_rotation_vars.py): the antenna pattern is evaluated at
the body-frame direction (theta_prime, phi_prime) = R^T applied to the
observation direction. A separately named copy is used, rather than
instantiating Ftheta_Fphi_vars twice with a shadowed theta, so that the
lab-frame and body-frame Global-scope names never collide.

Usage (ADS Python console, after start_console.py):
    run_script("create_Fdata_prime_vars.py")

Requires SphericalHarmonics_lib (create_Fdata_vars.py). The target cell is deleted
and rebuilt if it exists. Afterwards, check that both VARs have Global scope.
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
CELL_NAME = "Ftheta_Fphi_vars_prime"   # distinct from Ftheta_Fphi_vars

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
    largest index in the data. The largest index occurs at l=lMax,
    m=lMax: i_max = lMax(lMax+1)+lMax = lMax(lMax+2), so
    lMax = -1 + sqrt(1+i_max)."""
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


def index_to_lm_map(l_max):
    """Rebuild the i = l(l+1)+m -> (l,m) mapping used by the Wolfram export."""
    mapping = {}
    for l in range(0, l_max + 1):
        for m in range(-l, l + 1):
            i = l * (l + 1) + m
            mapping[i] = (l, m)
    return mapping


def var_name(prefix, l, m):
    """ADS-identifier-safe variable name, e.g. Ftheta_prime_l1_mneg1,
    Ftheta_prime_l1_m0."""
    m_part = f"mneg{abs(m)}" if m < 0 else f"m{m}"
    return f"{prefix}_l{l}_{m_part}"


# Whole-word substitution: "theta" -> "theta_prime", "phi" -> "phi_prime".
THETA_TO_THETA_PRIME = re.compile(r"\btheta\b")
PHI_TO_PHI_PRIME = re.compile(r"\bphi\b")


def prime_formula(formula):
    """Substitute theta -> theta_prime and phi -> phi_prime (both
    whole-word) in a raw formula string."""
    formula = THETA_TO_THETA_PRIME.sub("theta_prime", formula)
    formula = PHI_TO_PHI_PRIME.sub("phi_prime", formula)
    return formula


def populate_var_prime(var_inst, data, prefix, idx_to_lm):
    """Set var_inst.vars[name] = prime_formula(formula) for every entry,
    then delete the default 'X' placeholder that a new VAR instance
    always starts with."""
    for i, formula in data.items():
        l, m = idx_to_lm[i]
        name = var_name(prefix, l, m)
        var_inst.vars[name] = prime_formula(formula)
    del var_inst.vars["X"]


# ---- Parse the exported Discrete-format files (same files as create_Fdata_vars.py) ----
ftheta_data, ftheta_lmax_header = parse_dscrdata(FTHETA_FILE)
fphi_data, fphi_lmax_header = parse_dscrdata(FPHI_FILE)

L_MAX_ftheta = infer_lmax(ftheta_data, ftheta_lmax_header)
L_MAX_fphi = infer_lmax(fphi_data, fphi_lmax_header)

if L_MAX_ftheta != L_MAX_fphi:
    print(f"WARNING: inferred lMax differs between the F_theta file "
          f"({L_MAX_ftheta}) and F_phi file ({L_MAX_fphi}) -- using "
          f"{L_MAX_ftheta}. Check that both files came from the same "
          f"Wolfram export run.")

L_MAX = L_MAX_ftheta
idx_to_lm = index_to_lm_map(L_MAX)

print(f"Parsed {len(ftheta_data)} F_theta entries, {len(fphi_data)} F_phi "
      f"entries. Inferred lMax = {L_MAX} "
      f"({'from REM header' if ftheta_lmax_header is not None else 'from max index'}).")

# ---- Get the active workspace (library should already exist from create_Fdata_vars.py) ----
ws = de.active_workspace()
print(f"Active workspace: {ws}")

if LIBRARY_NAME not in ws.library_names:
    raise RuntimeError(
        f"Library '{LIBRARY_NAME}' not found -- run create_Fdata_vars.py first."
    )
lib = ws.open_library(LIBRARY_NAME)

# ---- Create the schematic and the two VAR instances ----
# Rebuild: delete the cell first if it exists (create_schematic cannot
# overwrite a view, and delete_cell needs the cell's views closed).
if lib.cell_exists(CELL_NAME):
    print(f"Cell '{CELL_NAME}' already exists in '{LIBRARY_NAME}' -- deleting it first.")
    for view_name in ("schematic", "symbol"):
        try:
            de.db_uu.open_design((LIBRARY_NAME, CELL_NAME, view_name)).close_design()
        except Exception:
            pass
    lib.delete_cell(CELL_NAME)
design = de.db_uu.create_schematic((LIBRARY_NAME, CELL_NAME, "schematic"))

var_ftheta_prime = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 0), name="VAR_Ftheta_prime", angle=90)
assert var_ftheta_prime.is_var_instance

var_fphi_prime = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 3), name="VAR_Fphi_prime", angle=90)
assert var_fphi_prime.is_var_instance

with de.db_uu.Transaction(design) as transaction:
    # theta_prime/phi_prime are not defined here; they come from the parent
    # design (nested scope), typically from VectorFieldRotate_vars.
    populate_var_prime(var_ftheta_prime, ftheta_data, "Ftheta_prime", idx_to_lm)
    populate_var_prime(var_fphi_prime, fphi_data, "Fphi_prime", idx_to_lm)

    # Global scope so the parent can reference Ftheta_prime_*/Fphi_prime_*.
    # The names differ from Ftheta_l*/Fphi_l*, so both cells can be live
    # in the same design.
    var_ftheta_prime.set_global_scope()
    var_fphi_prime.set_global_scope()

    transaction.commit()

design.save_design()
print(f"Created {len(ftheta_data)} Ftheta_prime variables in VAR_Ftheta_prime and "
      f"{len(fphi_data)} Fphi_prime variables in VAR_Fphi_prime, in {LIBRARY_NAME}:{CELL_NAME}")
print("Reminder: manually verify Global scope is set on both VAR_Ftheta_prime "
      "and VAR_Fphi_prime via Component Options in the ADS GUI.")

# ---- Lightweight spot check ----
sample_i = next(iter(ftheta_data))
sample_l, sample_m = idx_to_lm[sample_i]
sample_ftheta_prime_name = var_name("Ftheta_prime", sample_l, sample_m)
print(f"Spot check: {sample_ftheta_prime_name} = {var_ftheta_prime.vars[sample_ftheta_prime_name]}")

sample_i2 = next(iter(fphi_data))
sample_l2, sample_m2 = idx_to_lm[sample_i2]
sample_fphi_prime_name = var_name("Fphi_prime", sample_l2, sample_m2)
print(f"Spot check: {sample_fphi_prime_name} = {var_fphi_prime.vars[sample_fphi_prime_name]}")
