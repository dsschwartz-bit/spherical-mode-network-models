r"""
create_vcvs_tx_gain.py -- build the transmit-side mode-to-field cell.

Defines create_vcvs_tx_gain(snp_path, ...), which creates the schematic and
symbol "{antenna}_VCVS_TxGain" from a Touchstone file's mode ports
(TM_L{l}_M{m} / TE_L{l}_M{m}; negative m written "Mneg{|m|}"; other ports
such as "Feed" are ignored). Each mode gets an input pin (west) and two VCVS
whose outputs are series-chained into the theta and phi field voltages,
exposed as Vtheta_prime and Vphi_prime (east), with a common ref pin
(south). Modes are placed in SnP header order. GradYlm_vars_prime is
instantiated inside the cell as GradYlmPrime_sub.

Gains (antenna body frame, see create_gradYlm_vars_prime.py):
    TM: G_theta = -gradYlm_theta_prime,  G_phi = -gradYlm_phi_prime
    TE: G_theta = -gradYlm_phi_prime,    G_phi = +gradYlm_theta_prime
E_TM is proportional to -grad Y; by duality E_TE ~ rhat x grad Y, and
rhat x thetahat = phihat, rhat x phihat = -thetahat give the TE row. There
is no parity(l) factor on the transmit side (cf. create_cccs_rx_aperture.py).

VCVS terminals: 1 = control+, 2 = control-, 3 = output+, 4 = output-.
Control sides tap each mode's input net in parallel (R1 = 1e100 Ohm, no
loading); output sides are series-chained by wires (term4 of one row to
term3 of the next), so the chain voltage is the sum of all terms.

Usage (ADS Python console, after start_console.py):
    run_script("create_vcvs_tx_gain.py")
    create_vcvs_tx_gain(r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p")
    create_vcvs_tx_gain(r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p",
                        new_cell_name="MyAntenna_VCVS_TxGain",
                        library_name="spice_import_lib")

Requires SphericalHarmonics_lib:GradYlm_vars_prime
(create_gradYlm_vars_prime.py). The parent design must define r, a_Tx, c0.
"""

import re
from pathlib import Path
import keysight.ads.de as de
from keysight.ads.de import db_uu
from keysight.ads.de.experimental import generate_symbol as gs

# ==================== Fixed configuration -- edit as needed ====================

REF_NET = "ref"
MODE_PATTERN = re.compile(r"^(TM|TE)_L(\d+)_M(neg)?(\d+)$")
# Matches the antenna-name convention S_dom_psi_{name}_... or
# S_dom_Q_{name}_... -- the name is the single token immediately after
# 'psi_' or 'Q_', up to the next underscore. E.g.
# 'S_dom_psi_dipole_cylinder_edgeport_symmetry_full.s6p' -> 'dipole'.
ANTENNA_NAME_PATTERN = re.compile(r"S_dom_(?:psi|Q)_([A-Za-z0-9]+)_")

# VCVS footprint (bounding box of a placed ads_sources:VCVS instance).
VCVS_WIDTH = 1.75
VCVS_HEIGHT = 1.06875
ROW_PITCH = VCVS_HEIGHT * 2         # component height + full-height gap
COL_PITCH = 2 * VCVS_WIDTH          # theta column / phi column horizontal separation
PIN_CLEARANCE = 2.0
GRID_START_Y = 0.0

# ============================================================================


def get_port_names_from_snp_comments(filepath):
    """Parses '! Port[N] Name' style header comments directly from the
    Touchstone file, in port order."""
    port_names = {}
    pattern = re.compile(r"^!\s*Port\[(\d+)\]\s+(\S+)", re.IGNORECASE)
    with open(filepath, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line.startswith("!"):
                if line and not line.startswith("#"):
                    break
                continue
            match = pattern.match(line)
            if match:
                port_names[int(match.group(1))] = match.group(2)
    if not port_names:
        raise ValueError(f"No '! Port[N] Name' comments found in {filepath}")
    max_port = max(port_names.keys())
    return [port_names.get(i, f"Port{i}") for i in range(1, max_port + 1)]


def var_name(prefix, l, m):
    """ADS-safe variable name, e.g. gradYlm_theta_prime_l1_mneg1 (signed
    m, same naming as create_gradYlm_vars_prime.py)."""
    m_part = f"mneg{abs(m)}" if m < 0 else f"m{m}"
    return f"{prefix}_l{l}_{m_part}"


def pt(point):
    """Coerce a PointF-like object (or plain tuple) to a plain (x,y) tuple."""
    return (point.x, point.y) if hasattr(point, "x") else tuple(point)


def extract_antenna_name(snp_path):
    """Extracts the antenna name from a Touchstone filename following
    the S_dom_psi_{name}_... / S_dom_Q_{name}_... convention -- see
    ANTENNA_NAME_PATTERN above."""
    filename = Path(snp_path).name
    match = ANTENNA_NAME_PATTERN.search(filename)
    if not match:
        raise ValueError(
            f"Could not extract antenna name from filename '{filename}' -- "
            f"expected the pattern S_dom_psi_<name>_... or S_dom_Q_<name>_..."
        )
    return match.group(1)


def build_vcvs_gain_formulas(mode_type, l, m):
    """Return (theta_gain, phi_gain) ADS expressions for mode_type "TM" or
    "TE" and (l, m). See the module docstring for the formulas."""
    if mode_type == "TM":
        theta_gain = f"(-{var_name('gradYlm_theta_prime', l, m)})"
        phi_gain = f"(-{var_name('gradYlm_phi_prime', l, m)})"
    elif mode_type == "TE":
        theta_gain = f"(-{var_name('gradYlm_phi_prime', l, m)})"
        phi_gain = f"({var_name('gradYlm_theta_prime', l, m)})"
    else:
        raise ValueError(f"Unknown mode_type '{mode_type}' (expected 'TM' or 'TE').")
    return theta_gain, phi_gain


def get_ka_factor_is_one(filepath):
    """True if the SnP header contains 'ka_factor=1' (data already
    referenced to the antenna's own radius, so no a_Tx/a_Rx correction in
    the delay); False otherwise, including when no ka_factor comment is
    present."""
    pattern = re.compile(r"ka_factor\s*=\s*([^\s,]+)", re.IGNORECASE)
    with open(filepath, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line.startswith("!"):
                if line and not line.startswith("#"):
                    break
                continue
            match = pattern.search(line)
            if match:
                return match.group(1).strip() == "1"
    return False


def create_vcvs_tx_gain(
    snp_path,
    library_name="SphericalHarmonics_lib",
    gradylm_library_name="SphericalHarmonics_lib",
    gradylm_prime_cell_name="GradYlm_vars_prime",
    new_cell_name=None,
    delay_formula=None,
):
    """Build the {antenna}_VCVS_TxGain schematic and symbol.

    Args:
        snp_path: Touchstone file whose TM_/TE_ mode ports define the
            input pins and VCVS gains.
        library_name: library the new cell is created in.
        gradylm_library_name: library holding GradYlm_vars_prime
            (independent of library_name).
        gradylm_prime_cell_name: cell providing gradYlm_theta_prime_*/
            gradYlm_phi_prime_* (create_gradYlm_vars_prime.py).
        new_cell_name: defaults to "{antenna_name}_VCVS_TxGain".
        delay_formula: delay for every VCVS, with " sec" suffix; stored
            (without suffix) in the cell-local VAR "prop_t". Default:
            "r/c0 sec" if the SnP header says ka_factor=1, else
            "((r-a_Tx)/c0) sec". create_tx_rx_link_internal_ports_TRANSIENT.py
            passes its own combined (causal) Tx+Rx delay here.

    Returns:
        The created schematic design.
    """
    if delay_formula is None:
        if get_ka_factor_is_one(snp_path):
            delay_formula = "r/c0 sec"
            print("ka_factor=1 detected in snp_path's own header comment "
                  "-- using delay_formula='r/c0 sec' (no a_Tx de-embedding "
                  "correction).")
        else:
            delay_formula = "((r-a_Tx)/c0) sec"

    if new_cell_name is None:
        antenna_name = extract_antenna_name(snp_path)
        new_cell_name = f"{antenna_name}_VCVS_TxGain"
        print(f"Antenna name '{antenna_name}' extracted from filename -- "
              f"using cell name '{new_cell_name}'.")

    # ---- Parse the Touchstone file, filter to TM/TE modes ----
    # Placement order is exactly the order modes are written in the
    # SnP file header.
    all_port_names = get_port_names_from_snp_comments(snp_path)
    modes = []
    for name in all_port_names:
        match = MODE_PATTERN.match(name)
        if match:
            mode_type = match.group(1)
            l = int(match.group(2))
            m = -int(match.group(4)) if match.group(3) else int(match.group(4))
            modes.append((mode_type, l, m, name))
    print(f"Found {len(all_port_names)} ports total, {len(modes)} modes "
          f"(non-mode ports like 'Feed' ignored), in header order: "
          f"{[mo[3] for mo in modes]}")

    # ---- Get the active workspace and library (must already exist) ----
    ws = de.active_workspace()
    print(f"Active workspace: {ws}")

    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    lib = ws.open_library(library_name)

    if gradylm_library_name not in ws.library_names:
        raise RuntimeError(f"Library '{gradylm_library_name}' not found.")
    ws.open_library(gradylm_library_name)

    # ---- Confirm the existing GradYlm_vars_prime design exists (read-only) ----
    gradylm_prime_design = de.db_uu.open_design((gradylm_library_name, gradylm_prime_cell_name, "schematic"))
    print(f"Found existing design: {gradylm_prime_design.design_name}")

    # ---- Create the new schematic ----
    # Delete first if it already exists -- re-running for the same
    # antenna otherwise fails (create_schematic can't create a view
    # that's already there). Views must be closed first, or
    # delete_cell fails with "cell has open cellviews".
    if lib.cell_exists(new_cell_name):
        print(f"Cell '{new_cell_name}' already exists in '{library_name}' "
              f"-- deleting it first.")
        for view_name in ("schematic", "symbol"):
            try:
                de.db_uu.open_design((library_name, new_cell_name, view_name)).close_design()
            except Exception:
                pass
        lib.delete_cell(new_cell_name)
    design = de.db_uu.create_schematic((library_name, new_cell_name, "schematic"))

    # Pin/subcircuit geometry (depends only on len(modes)).
    west_x = 0 - PIN_CLEARANCE
    east_x = COL_PITCH + VCVS_WIDTH + PIN_CLEARANCE
    ref_x = (0 + COL_PITCH + VCVS_WIDTH) / 2
    ref_y = GRID_START_Y - (len(modes) * ROW_PITCH) - PIN_CLEARANCE
    last_row_y = GRID_START_Y - ((len(modes) - 1) * ROW_PITCH)
    vtheta_prime_y = last_row_y
    vphi_prime_y = last_row_y - ROW_PITCH / 2

    n_vcvs = 0
    with de.db_uu.Transaction(design) as transaction:
        # prop_t holds the delay expression once (without the " sec"
        # suffix); every VCVS's "T" is set to "prop_t sec".
        prop_t_value = delay_formula[:-len(" sec")].strip() if delay_formula.endswith(" sec") else delay_formula
        var_inst = design.add_instance(
            ("ads_datacmps", "VAR", "symbol"), (west_x, GRID_START_Y + ROW_PITCH), name="VAR_prop_t", angle=90
        )
        var_inst.vars["prop_t"] = prop_t_value
        del var_inst.vars["X"]

        # Instantiate GradYlm_vars_prime as a subcircuit so its
        # global-scope gradYlm_theta_prime_*/gradYlm_phi_prime_*
        # variables become visible here, placed to the RIGHT of the
        # Vtheta_prime/Vphi_prime output pins, vertically between the two.
        design.add_instance(gradylm_prime_design.design_name,
                             (east_x + PIN_CLEARANCE, (vtheta_prime_y + vphi_prime_y) / 2),
                             name="GradYlmPrime_sub", angle=0)

        prev_inst_theta = None
        prev_inst_phi = None
        for row_idx, (mode_type, l, m, port_name) in enumerate(modes):
            y = GRID_START_Y - (row_idx * ROW_PITCH)
            is_last = (row_idx == len(modes) - 1)

            theta_gain, phi_gain = build_vcvs_gain_formulas(mode_type, l, m)

            # theta-chain VCVS
            inst_t = design.add_instance(("ads_sources", "VCVS", "symbol"),
                                          (0, y), name=f"VCVS_TxTheta_{port_name}", angle=0)
            inst_t.parameters["G"].value = theta_gain
            inst_t.parameters["T"].value = "prop_t sec"
            pins_t = {p.inst_term.term_number: p for p in inst_t.inst_pins}
            pins_t[1].add_label(port_name, pt(pins_t[1].snap_point))       # control+ <- mode input
            pins_t[2].add_label(REF_NET, pt(pins_t[2].snap_point))         # control- -> ref (parallel tap)
            # Free ends of the output chain: first instance's term3 is the
            # Vtheta_prime output, last instance's term4 is ref. All other
            # term3/term4 pins are joined by the wire drawn on the next row.
            if row_idx == 0:
                pins_t[3].add_label("Vtheta_prime", pt(pins_t[3].snap_point))  # first: output+ -> final output
            else:
                prev_pins_t = {p.inst_term.term_number: p for p in prev_inst_theta.inst_pins}
                design.add_wire([pt(prev_pins_t[4].snap_point), pt(pins_t[3].snap_point)])
            if is_last:
                pins_t[4].add_label(REF_NET, pt(pins_t[4].snap_point))     # last: output- -> ref
            n_vcvs += 1
            prev_inst_theta = inst_t

            # phi-chain VCVS
            inst_p = design.add_instance(("ads_sources", "VCVS", "symbol"),
                                          (COL_PITCH, y), name=f"VCVS_TxPhi_{port_name}", angle=0)
            inst_p.parameters["G"].value = phi_gain
            inst_p.parameters["T"].value = "prop_t sec"
            pins_p = {p.inst_term.term_number: p for p in inst_p.inst_pins}
            pins_p[1].add_label(port_name, pt(pins_p[1].snap_point))
            pins_p[2].add_label(REF_NET, pt(pins_p[2].snap_point))
            if row_idx == 0:
                pins_p[3].add_label("Vphi_prime", pt(pins_p[3].snap_point))
            else:
                prev_pins_p = {p.inst_term.term_number: p for p in prev_inst_phi.inst_pins}
                design.add_wire([pt(prev_pins_p[4].snap_point), pt(pins_p[3].snap_point)])
            if is_last:
                pins_p[4].add_label(REF_NET, pt(pins_p[4].snap_point))
            n_vcvs += 1
            prev_inst_phi = inst_p

        transaction.commit()

    design.save_design()
    print(f"Created {n_vcvs} VCVS instances ({len(modes)} modes x 2 chains) in "
          f"{library_name}:{new_cell_name}")

    # ---- External pins ----
    pin_layer_id = db_uu.LayerId(231)

    pin_defs = []
    for row_idx, (mode_type, l, m, port_name) in enumerate(modes):
        y = GRID_START_Y - (row_idx * ROW_PITCH)
        pin_defs.append((port_name, (west_x, y), 180.0, db_uu.TermType.INPUT))
    pin_defs.append(("Vtheta_prime", (east_x, vtheta_prime_y), 0.0, db_uu.TermType.OUTPUT))
    pin_defs.append(("Vphi_prime", (east_x, vphi_prime_y), 0.0, db_uu.TermType.OUTPUT))
    pin_defs.append(("ref", (ref_x, ref_y), -90.0, db_uu.TermType.INPUT_OUTPUT))

    with de.db_uu.Transaction(design) as transaction:
        for name, pos, angle, ttype in pin_defs:
            net = design.find_or_add_net(name)
            term = design.add_term(net, name, ttype)
            dot = design.add_dot(pin_layer_id, pos)
            design.add_pin(term, dot, angle=angle)
        transaction.commit()

    design.save_design()
    print(f"Added {len(pin_defs)} external pins to {library_name}:{new_cell_name}.")

    # ---- Symbol generation ----
    symbol = db_uu.create_symbol(f"{library_name}:{new_cell_name}:symbol")
    symbol_generator = gs.SymbolGenerator(symbol, design, 0.25, 0.25)
    symbol_generator.is_dual_symbol_type = False
    symbol_generator.should_replace = True
    symbol_generator.pin_shape = "dot"
    symbol_generator.order = gs.OrderType.ORDER_LOCATION
    symbol_generator.generate_symbol()
    symbol.save_design()
    print(f"Generated Quad symbol for {library_name}:{new_cell_name}.")
    print("Usage reminder: the parent design must supply r, a_Tx, c0 as "
          "VAR-scope local variables (for the T delay parameter), and each "
          "mode's input pin should be driven by that mode's excitation.")

    return design


print("Ready: create_vcvs_tx_gain(snp_path, ...) is defined. Call it with "
      "your Touchstone file path to build the cell -- see the module "
      "docstring's Usage section.")
