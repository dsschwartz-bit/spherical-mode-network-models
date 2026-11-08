r"""
create_vector_field_rotation_vars.py -- rotate the radiated field of an antenna.

Creates Rotation_lib:VectorFieldRotate_vars, which evaluates the field of an
antenna rotated by R(rot_angle, Ux, Uy, Uz) at a fixed lab-frame observation
point (r, theta, phi):

    E'(r) = R . E(R^T . r)

The observation direction is mapped back to the body-frame direction
(theta_prime, phi_prime) = R^T applied to (theta, phi); the antenna's native
field is sampled there and the field vector is rotated forward by R:

    (x, y, z)        = r * rhat(theta, phi)                observation point
    (x', y', z')     = R^T . (x, y, z)                      body-frame point
    (Vx, Vy, Vz)     = Sph2Cart(theta', phi') . (Vr', Vtheta', Vphi')
    (Vx', Vy', Vz')  = R . (Vx, Vy, Vz)
    (Vr, Vtheta, Vphi)_out = Cart2Sph(theta, phi) . (Vx', Vy', Vz')

Primes are written as a _prime suffix in ADS names. Cart2Sph_ij at the
observation point uses the angle form from Cart2Sph_ADS_Equations.txt.
Sph2Cart at the body-frame point uses the generic Cartesian template M_ij
(and r_of_xyz, rho_of_xyz, ...) from the same file, transposed (M_ji) and
evaluated at x_prime, y_prime, z_prime, so no inverse trig is needed inside
the matrix entries. The two matrices are at different points and are not
transposes of each other.

Pin names: the input field (native pattern, from VCVS_Tx_Gain) is on
Vr_prime, Vtheta_prime, Vphi_prime (west); the rotated output field is on
Vr, Vtheta, Vphi (east); ref is south. This deliberately reverses the
prime convention for these pins so that names chain directly
VCVS_Tx_Gain -> VectorFieldRotate_vars -> receive side.

Implementation: two Global-scope VARs (VAR_ObservationFrame,
VAR_BodyFrame), a RotationMatrix_vars subcircuit, and 27 VCVS (3 stages x
3 rows x 3 terms). Each matrix row is a series chain of three VCVS
(ideal sources in series add; their outputs cannot be paralleled since
R2 = 0), connected by net labels. VCVS terminals: 1 = control+,
2 = control-, 3 = output+, 4 = output-; defaults (R1 = 1e100 Ohm, R2 = 0,
T = 0) are kept and only G is set. theta, phi, r, rot_angle, Ux, Uy, Uz
come from the parent design; theta_prime, phi_prime are exported (Global)
for GradYlm_vars_prime.

Usage (ADS Python console, after start_console.py):
    run_script("create_vector_field_rotation_vars.py")

Requires Rotation_lib:RotationMatrix_vars (create_rotationmat_vars.py). The
target cell is deleted and rebuilt if it exists. Afterwards, check that both VARs have
Global scope.
"""

from pathlib import Path
import re
import keysight.ads.de as de
from keysight.ads.de import db_uu
from keysight.ads.de.experimental import generate_symbol as gs

# ==================== Configuration -- edit as needed ====================

# Repository folders (defined by start_console.py; derived here if this file is exec()'d directly).
try:
    SCRIPT_DIR
except NameError:
    SCRIPT_DIR = Path(THIS_SCRIPT_PATH).resolve().parent
DATA_DIR = SCRIPT_DIR.parents[1] / "mathematica" / "tables"
SPHCART_FILE = DATA_DIR / "Cart2Sph_ADS_Equations.txt"

LIBRARY_NAME = "Rotation_lib"
ROTATIONMAT_CELL_NAME = "RotationMatrix_vars"   # existing, from create_rotationmat_vars.py
NEW_CELL_NAME = "VectorFieldRotate_vars"        # new schematic created by this script

REF_NET = "ref"   # common ground/reference net; must match the symbol's "ref" pin

# ============================================================================


def parse_ads_equations(filepath):
    """Parse a flat "name=formula;" equation file (as written by
    BuildADSEquationFile[] in vector_field_spherical_coordinate_transform.wl) into an
    ordered {name: formula} dict. Lines starting with "//" (or blank
    lines) are treated as comments/whitespace and skipped; every other
    non-blank line must be 'name=formula;' (trailing ';' optional)."""
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


def matrix_var_name(prefix, i, j):
    """e.g. matrix_var_name('Cart2Sph', 1, 2) -> 'Cart2Sph_12'."""
    return f"{prefix}_{i}{j}"


def substitute_vars(formula, mapping):
    """Whole-word substitute multiple variable names simultaneously, in a
    SINGLE regex pass -- avoids any risk of one substitution's output
    being re-matched by a later pattern. Keys are tried longest-first
    (e.g. "rho" before "r") so a longer name is never shadowed by a
    shorter one matching a prefix of it at the same position."""
    keys = sorted(mapping, key=len, reverse=True)
    pattern = re.compile(r"\b(" + "|".join(re.escape(k) for k in keys) + r")\b")
    return pattern.sub(lambda m: mapping[m.group(1)], formula)


# ---- Read Cart2Sph_ij (angle form) and the generic M_ij template straight from the Wolfram export ----
sphcart_vars = parse_ads_equations(SPHCART_FILE)
print(f"Parsed {len(sphcart_vars)} entries from {SPHCART_FILE.name}.")

expected_cart2sph_angle = [matrix_var_name("Cart2Sph", i, j) for i in (1, 2, 3) for j in (1, 2, 3)]
expected_generic_m = [matrix_var_name("M", i, j) for i in (1, 2, 3) for j in (1, 2, 3)]
expected_aux = ["r_of_xyz", "rho_of_xyz", "theta_of_xyz", "phi_of_xyz_raw", "phi_of_xyz"]
missing = [n for n in expected_cart2sph_angle + expected_generic_m + expected_aux if n not in sphcart_vars]
if missing:
    raise ValueError(f"Missing expected entries in {SPHCART_FILE.name}: {missing}")

# Mapping used to instantiate the generic M(x,y,z,r,rho) template at the
# BODY-FRAME point: x,y,z,r,rho -> x_prime,y_prime,z_prime,r_prime,rho_prime, and
# (for phi_of_xyz's internal reference) phi_of_xyz_raw -> phi_prime_raw.
BODY_FRAME_MAPPING = {
    "x": "x_prime", "y": "y_prime", "z": "z_prime", "r": "r_prime", "rho": "rho_prime",
    "phi_of_xyz_raw": "phi_prime_raw",
}

# ---- Get the active workspace and library (must already exist) ----
ws = de.active_workspace()
print(f"Active workspace: {ws}")

if LIBRARY_NAME not in ws.library_names:
    raise RuntimeError(
        f"Library '{LIBRARY_NAME}' not found -- run create_rotationmat_vars.py first."
    )
lib = ws.open_library(LIBRARY_NAME)

# ---- Confirm the existing RotationMatrix_vars design exists (read-only) ----
rotmat_design = de.db_uu.open_design((LIBRARY_NAME, ROTATIONMAT_CELL_NAME, "schematic"))
print(f"Found existing design: {rotmat_design.design_name}")

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

var_obs = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 0), name="VAR_ObservationFrame", angle=90)
assert var_obs.is_var_instance

var_body = design.add_instance(("ads_datacmps", "VAR", "symbol"), (0, 3), name="VAR_BodyFrame", angle=90)
assert var_body.is_var_instance

with de.db_uu.Transaction(design) as transaction:
    # Instantiate RotationMatrix_vars as a subcircuit so its global-scope
    # R11..R33 (and magU) become visible here. rot_angle,Ux,Uy,Uz are not
    # defined by this design either -- like theta,phi,r, they're expected
    # to come from whatever design instantiates VectorFieldRotate_vars as a
    # subcircuit.
    design.add_instance(rotmat_design.design_name, (0, 6), name="RotationMatrix_sub", angle=0)

    # ---- VAR_ObservationFrame: (x,y,z) = r*rhat(theta,phi), the FIXED
    # observation position (needed to compute x_prime,y_prime,z_prime below);
    # Cart2Sph_ij read DIRECTLY from the angle-form export -- plain
    # trig in theta,phi, no substitution needed at all. ----
    var_obs.vars["x"] = "r*sin(theta)*cos(phi)"
    var_obs.vars["y"] = "r*sin(theta)*sin(phi)"
    var_obs.vars["z"] = "r*cos(theta)"
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            name = matrix_var_name("Cart2Sph", i, j)
            var_obs.vars[name] = sphcart_vars[name]
    del var_obs.vars["X"]
    var_obs.set_global_scope()

    # ---- VAR_BodyFrame: (x_prime,y_prime,z_prime) = R^T.(x,y,z) (R_ji),
    # the body-frame sampling position. r_prime, rho_prime, theta_prime,
    # phi_prime and Sph2Cart_ij (= M_ji) come from the generic
    # M(x,y,z,r,rho) template evaluated at that point via BODY_FRAME_MAPPING. ----
    for i in (1, 2, 3):
        out_name = ["x_prime", "y_prime", "z_prime"][i - 1]
        var_body.vars[out_name] = " + ".join(
            f"R{j}{i}*{['x', 'y', 'z'][j - 1]}" for j in (1, 2, 3)
        )
    var_body.vars["r_prime"] = substitute_vars(sphcart_vars["r_of_xyz"], BODY_FRAME_MAPPING)
    var_body.vars["rho_prime"] = substitute_vars(sphcart_vars["rho_of_xyz"], BODY_FRAME_MAPPING)
    var_body.vars["theta_prime"] = substitute_vars(sphcart_vars["theta_of_xyz"], BODY_FRAME_MAPPING)
    var_body.vars["phi_prime_raw"] = substitute_vars(sphcart_vars["phi_of_xyz_raw"], BODY_FRAME_MAPPING)
    var_body.vars["phi_prime"] = substitute_vars(sphcart_vars["phi_of_xyz"], BODY_FRAME_MAPPING)
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            # Sph2Cart_ij = M_ji (index-transposed), then renamed to the
            # body-frame point -- both operations happening together, at
            # the SAME (body-frame) point.
            name = matrix_var_name("Sph2Cart", i, j)
            var_body.vars[name] = substitute_vars(sphcart_vars[matrix_var_name("M", j, i)], BODY_FRAME_MAPPING)
    del var_body.vars["X"]
    var_body.set_global_scope()

    transaction.commit()

# ============================================================================
# FIELD VECTOR PIPELINE: 27 VCVS (3 rows x 9 columns), series-chained per
# row and connected purely by net labels.
#
# Row i (= component index 1..3) holds, left to right, the three
# Sph2Cart terms, the three Rotate terms and the three Cart2Sph terms that
# produce Vx/Vx_prime/Vr (row 1), Vy/Vy_prime/Vtheta (row 2) and
# Vz/Vz_prime/Vphi (row 3). Pitch is twice the VCVS size in each direction.
# The Sph2Cart inputs and Cart2Sph outputs are where the swapped external
# pin names (see module docstring) apply.
# ============================================================================

# Each stage: (stage_label, output_nets, input_nets, gain_fn(i,j))
# gain_fn returns the VAR-scope name of the (i,j) matrix entry to use as
# that VCVS's G. Row i produces output_nets[i-1] from a weighted sum over
# input_nets, using gain_fn(i, 1..3) for the three terms.
stages = [
    ("Sph2Cart", ["Vx", "Vy", "Vz"], ["Vr_prime", "Vtheta_prime", "Vphi_prime"],
     lambda i, j: matrix_var_name("Sph2Cart", i, j)),
    ("Rotate", ["Vx_prime", "Vy_prime", "Vz_prime"], ["Vx", "Vy", "Vz"],
     lambda i, j: f"R{i}{j}"),
    ("Cart2Sph", ["Vr", "Vtheta", "Vphi"], ["Vx_prime", "Vy_prime", "Vz_prime"],
     lambda i, j: matrix_var_name("Cart2Sph", i, j)),
]

# Bounding box of a placed ads_sources:VCVS instance.
VCVS_WIDTH = 1.75
VCVS_HEIGHT = 1.06875
PITCH_X = 2 * VCVS_WIDTH    # gap between adjacent VCVS (horizontal) = one VCVS width
PITCH_Y = 2 * VCVS_HEIGHT   # gap between adjacent VCVS (vertical) = one VCVS height

n_vcvs = 0

with de.db_uu.Transaction(design) as transaction:
    for row_idx in range(3):  # OUTER loop: groups by component index i=row_idx+1
        i = row_idx + 1
        y = -(row_idx * PITCH_Y)
        col = 0
        for stage_idx, (stage_label, output_nets, input_nets, gain_fn) in enumerate(stages):
            output_net = output_nets[row_idx]
            prev_net = REF_NET
            for term_idx in range(3):  # term j = term_idx+1
                j = term_idx + 1
                is_last_term = (term_idx == 2)
                out_net = output_net if is_last_term else f"{stage_label}_{output_net}_n{term_idx + 1}"
                gain_expr = gain_fn(i, j)
                input_net = input_nets[term_idx]

                inst_name = f"VCVS_{stage_label}_{output_net}_{j}"
                x = col * PITCH_X
                inst = design.add_instance(("ads_sources", "VCVS", "symbol"), (x, y), name=inst_name, angle=0)

                inst.parameters["G"].value = gain_expr

                # term_number 1=control+, 2=control-, 3=output+, 4=output-
                pins_by_term = {p.inst_term.term_number: p for p in inst.inst_pins}
                pins_by_term[3].add_label(out_net, pins_by_term[3].snap_point)
                pins_by_term[4].add_label(prev_net, pins_by_term[4].snap_point)
                pins_by_term[1].add_label(input_net, pins_by_term[1].snap_point)
                pins_by_term[2].add_label(REF_NET, pins_by_term[2].snap_point)

                n_vcvs += 1
                prev_net = out_net
                col += 1

    transaction.commit()

design.save_design()
print(f"Created x,y,z + 9 Cart2Sph_ij vars in VAR_ObservationFrame; "
      f"x_prime,y_prime,z_prime,r_prime,rho_prime,theta_prime,phi_prime_raw,phi_prime + 9 "
      f"Sph2Cart_ij vars in VAR_BodyFrame; and {n_vcvs} VCVS instances "
      f"(field-vector pipeline) in {LIBRARY_NAME}:{NEW_CELL_NAME}")
print("Reminder: manually verify Global scope is set on both VAR "
      "components via Component Options in the ADS GUI.")

# ============================================================================
# EXTERNAL PINS + SYMBOL
#
# Pins sit level with each row's first term1 (inputs, west) and last term3
# (outputs, east); ref is centred below the grid (south). find_or_add_net()
# attaches each pin to the existing net of the same name created by the
# labels above, so no wiring is needed.
# ============================================================================

PIN_CLEARANCE = 2.0
LEFT_X = 0 - PIN_CLEARANCE
RIGHT_X = (8 * PITCH_X + VCVS_WIDTH) + PIN_CLEARANCE   # rightmost column's term3 x, plus clearance
REF_X = ((8 * PITCH_X + VCVS_WIDTH) + 0) / 2            # horizontally centered across the whole grid
REF_Y = -(2 * PITCH_Y) - VCVS_HEIGHT - PIN_CLEARANCE    # below the bottom row's body

PIN_LAYER_ID = db_uu.LayerId(231)  # schematic pin layer

pin_defs = [
    # (name, position, angle, term_type)
    ("Vr_prime",     (LEFT_X, 0 * -PITCH_Y),      180.0, db_uu.TermType.INPUT),
    ("Vtheta_prime", (LEFT_X, 1 * -PITCH_Y),       180.0, db_uu.TermType.INPUT),
    ("Vphi_prime",   (LEFT_X, 2 * -PITCH_Y),       180.0, db_uu.TermType.INPUT),
    ("Vr",           (RIGHT_X, 0 * -PITCH_Y),        0.0, db_uu.TermType.OUTPUT),
    ("Vtheta",       (RIGHT_X, 1 * -PITCH_Y),        0.0, db_uu.TermType.OUTPUT),
    ("Vphi",         (RIGHT_X, 2 * -PITCH_Y),        0.0, db_uu.TermType.OUTPUT),
    ("ref",          (REF_X, REF_Y),                -90.0, db_uu.TermType.INPUT_OUTPUT),
]

with de.db_uu.Transaction(design) as transaction:
    for name, pos, angle, ttype in pin_defs:
        net = design.find_or_add_net(name)
        term = design.add_term(net, name, ttype)
        dot = design.add_dot(PIN_LAYER_ID, pos)
        pin = design.add_pin(term, dot, angle=angle)
    transaction.commit()

design.save_design()
print(f"Added {len(pin_defs)} external pins to {LIBRARY_NAME}:{NEW_CELL_NAME} "
      f"(Vr_prime/Vtheta_prime/Vphi_prime west, Vr/Vtheta/Vphi east, ref south).")

# ---- Auto-generate the symbol from these pins ----
symbol = db_uu.create_symbol(f"{LIBRARY_NAME}:{NEW_CELL_NAME}:symbol")
symbol_generator = gs.SymbolGenerator(symbol, design, 0.25, 0.25)
symbol_generator.is_dual_symbol_type = False   # False = Quad
symbol_generator.should_replace = True
symbol_generator.pin_shape = "dot"
symbol_generator.order = gs.OrderType.ORDER_LOCATION   # keeps the schematic pin geometry
symbol_generator.generate_symbol()
symbol.save_design()
print(f"Generated Quad symbol for {LIBRARY_NAME}:{NEW_CELL_NAME}.")
print("Usage reminder: the parent design (wherever VectorFieldRotate_vars "
      "is instantiated as a subcircuit) must supply theta, phi (the fixed "
      "observation direction), r, rot_angle, Ux, Uy, Uz as VAR-scope local "
      "variables, and must feed THIS design's own theta_prime,phi_prime "
      "output to GradYlm_vars_prime's theta,phi input to correctly source "
      "the Vr_prime,Vtheta_prime,Vphi_prime voltage nodes this design's "
      "west-side pins expect.")

# ---- Lightweight spot checks ----
print(f"Spot check: x_prime = {var_body.vars['x_prime']}")
print(f"Spot check: theta_prime = {var_body.vars['theta_prime']}")
print(f"Spot check: phi_prime = {var_body.vars['phi_prime']}")
print(f"Spot check: Sph2Cart_12 = {var_body.vars['Sph2Cart_12']}")
print(f"Spot check: Cart2Sph_12 = {var_obs.vars['Cart2Sph_12']}")
