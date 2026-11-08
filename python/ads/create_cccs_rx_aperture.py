r"""
create_cccs_rx_aperture.py -- build the receive-side mode-projection cell.

Defines create_cccs_rx_aperture(snp_path, ...), which creates the schematic
and symbol "{antenna}_CCCS_RxAperture" from a Touchstone file's mode ports
(TM_L{l}_M{m} / TE_L{l}_M{m}; negative m written "Mneg{|m|}"; other ports
such as "Feed" are ignored). Inputs Itheta, Iphi (west) are the field
currents; the mu0*r inductors that produce them from Vtheta/Vphi live in the
parent link schematic. Each mode gets two CCCS (theta- and phi-sensing),
whose control sides form two series chains and whose outputs are paralleled
onto that mode's output pin (east). Modes are placed in SnP header order.

Gains, with parity(l) = (-1)^(l-1) (receive side only) and 2*(2*pi) the
receive-aperture factor:
    TM: G_theta = parity*2*(2*pi)*(-gradYlm_theta),  G_phi = parity*2*(2*pi)*(-gradYlm_phi)
    TE: G_theta = parity*2*(2*pi)*(-gradYlm_phi),    G_phi = parity*2*(2*pi)*(+gradYlm_theta)
(E_TM ~ -grad Y, E_TE ~ rhat x grad Y; see create_vcvs_tx_gain.py.) The
unprimed GradYlm_vars are used: the receiver projects onto its own mode
basis. CCCS terminals: 1 = control+, 2 = control-, 3/4 = output -/+ for TM;
TE swaps 3/4, which applies an extra sign flip. R2 keeps its default
(1e100 Ohm, ideal current source).

Usage (ADS Python console, after start_console.py):
    run_script("create_cccs_rx_aperture.py")
    create_cccs_rx_aperture(r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p")
    create_cccs_rx_aperture(r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p",
                            new_cell_name="MyAntenna_CCCS_RxAperture",
                            library_name="spice_import_lib")

Requires SphericalHarmonics_lib:GradYlm_vars (create_gradYlm_vars.py). The
parent design must define a_Rx and c0 (used by the delay).
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

# CCCS footprint (same symbol size as VCVS).
CCCS_WIDTH = 1.75
CCCS_HEIGHT = 1.06875
ROW_PITCH = CCCS_HEIGHT * 2   # component height + full-height gap
COL_PITCH = 2 * CCCS_WIDTH
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
    """ADS-safe variable name, e.g. gradYlm_theta_l1_mneg1 (signed m,
    same naming as create_gradYlm_vars.py)."""
    m_part = f"mneg{abs(m)}" if m < 0 else f"m{m}"
    return f"{prefix}_l{l}_{m_part}"


def parity_prefix(l):
    """(-1)^(l-1) as an expression prefix: '' for odd l, '-' for even l.
    Applies to both TM and TE; the far-field negative-gradient sign is a
    separate factor in build_cccs_gain_formulas."""
    return "" if l % 2 == 1 else "-"


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


def build_cccs_gain_formulas(mode_type, l, m):
    """Return (theta_gain, phi_gain) ADS expressions for mode_type "TM" or
    "TE" and (l, m). See the module docstring for the formulas."""
    sign = parity_prefix(l)
    if mode_type == "TM":
        theta_gain = f"{sign}2*(2*pi)*(-{var_name('gradYlm_theta', l, m)})"
        phi_gain = f"{sign}2*(2*pi)*(-{var_name('gradYlm_phi', l, m)})"
    elif mode_type == "TE":
        theta_gain = f"{sign}2*(2*pi)*(-{var_name('gradYlm_phi', l, m)})"
        phi_gain = f"{sign}2*(2*pi)*({var_name('gradYlm_theta', l, m)})"
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


def create_cccs_rx_aperture(
    snp_path,
    library_name="SphericalHarmonics_lib",
    gradylm_library_name="SphericalHarmonics_lib",
    gradylm_cell_name="GradYlm_vars",
    new_cell_name=None,
    delay_formula=None,
):
    """Build the {antenna}_CCCS_RxAperture schematic and symbol.

    Args:
        snp_path: Touchstone file whose TM_/TE_ mode ports define the
            output pins and CCCS gains.
        library_name: library the new cell is created in.
        gradylm_library_name: library holding GradYlm_vars (independent
            of library_name).
        gradylm_cell_name: cell providing the unprimed gradYlm_theta_*/
            gradYlm_phi_* variables (create_gradYlm_vars.py).
        new_cell_name: defaults to "{antenna_name}_CCCS_RxAperture".
        delay_formula: value of the cell-local VAR "prop_t" used as every
            CCCS's delay T. Default: "0" if the SnP header says
            ka_factor=1, else "-a_Rx/c0". A negative delay is fine for
            S-parameter analysis but non-causal in transient analysis;
            create_tx_rx_link_internal_ports_TRANSIENT.py passes "0" and
            moves the a_Rx term into the Tx delay.

    Returns:
        The created schematic design.
    """
    if delay_formula is None:
        if get_ka_factor_is_one(snp_path):
            delay_formula = "0"
            print("ka_factor=1 detected in snp_path's own header comment "
                  "-- using delay_formula='0' (no a_Rx de-embedding "
                  "correction).")
        else:
            delay_formula = "-a_Rx/c0"

    if new_cell_name is None:
        antenna_name = extract_antenna_name(snp_path)
        new_cell_name = f"{antenna_name}_CCCS_RxAperture"
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

    # ---- Confirm the existing GradYlm_vars design exists (read-only) ----
    gradylm_design = de.db_uu.open_design((gradylm_library_name, gradylm_cell_name, "schematic"))
    print(f"Found existing design: {gradylm_design.design_name}")

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

    west_x = 0 - PIN_CLEARANCE
    east_x = COL_PITCH + CCCS_WIDTH + PIN_CLEARANCE
    ref_x = (0 + COL_PITCH + CCCS_WIDTH) / 2
    ref_y = GRID_START_Y - (len(modes) * ROW_PITCH) - PIN_CLEARANCE
    last_row_y = GRID_START_Y - ((len(modes) - 1) * ROW_PITCH)
    output_y_center = (GRID_START_Y + last_row_y) / 2

    n_cccs = 0
    with de.db_uu.Transaction(design) as transaction:
        # prop_t holds the delay expression once, referenced by every
        # CCCS's "T" parameter. delay_formula has no " sec" suffix here.
        var_inst = design.add_instance(
            ("ads_datacmps", "VAR", "symbol"), (west_x, GRID_START_Y + ROW_PITCH), name="VAR_prop_t", angle=90
        )
        var_inst.vars["prop_t"] = delay_formula
        del var_inst.vars["X"]

        design.add_instance(gradylm_design.design_name,
                             (east_x + PIN_CLEARANCE, output_y_center),
                             name="GradYlm_sub", angle=0)

        prev_inst_theta = None
        prev_inst_phi = None
        for row_idx, (mode_type, l, m, port_name) in enumerate(modes):
            y = GRID_START_Y - (row_idx * ROW_PITCH)
            is_last = (row_idx == len(modes) - 1)

            theta_gain, phi_gain = build_cccs_gain_formulas(mode_type, l, m)

            output_net = port_name

            inst_t = design.add_instance(("ads_sources", "CCCS", "symbol"),
                                          (0, y), name=f"CCCS_RxTheta_{port_name}", angle=0)
            inst_t.parameters["G"].value = theta_gain
            inst_t.parameters["T"].value = "prop_t"
            inst_t.parameters["R1"].value = "0 Ohm"
            pins_t = {p.inst_term.term_number: p for p in inst_t.inst_pins}
            if row_idx == 0:
                pins_t[1].add_label("Itheta", pt(pins_t[1].snap_point))
            else:
                prev_pins_t = {p.inst_term.term_number: p for p in prev_inst_theta.inst_pins}
                design.add_wire([pt(prev_pins_t[2].snap_point), pt(pins_t[1].snap_point)])
            if is_last:
                pins_t[2].add_label(REF_NET, pt(pins_t[2].snap_point))
            # Output pins: TM uses term3=ref, term4=mode output; TE swaps
            # them, which negates the current source (the TE sign flip is
            # applied in the wiring, not in G).
            if mode_type == "TE":
                pins_t[3].add_label(output_net, pt(pins_t[3].snap_point))
                pins_t[4].add_label(REF_NET, pt(pins_t[4].snap_point))
            else:
                pins_t[3].add_label(REF_NET, pt(pins_t[3].snap_point))
                pins_t[4].add_label(output_net, pt(pins_t[4].snap_point))
            n_cccs += 1
            prev_inst_theta = inst_t

            inst_p = design.add_instance(("ads_sources", "CCCS", "symbol"),
                                          (COL_PITCH, y), name=f"CCCS_RxPhi_{port_name}", angle=0)
            inst_p.parameters["G"].value = phi_gain
            inst_p.parameters["T"].value = "prop_t"
            inst_p.parameters["R1"].value = "0 Ohm"
            pins_p = {p.inst_term.term_number: p for p in inst_p.inst_pins}
            if row_idx == 0:
                pins_p[1].add_label("Iphi", pt(pins_p[1].snap_point))
            else:
                prev_pins_p = {p.inst_term.term_number: p for p in prev_inst_phi.inst_pins}
                design.add_wire([pt(prev_pins_p[2].snap_point), pt(pins_p[1].snap_point)])
            if is_last:
                pins_p[2].add_label(REF_NET, pt(pins_p[2].snap_point))
            if mode_type == "TE":
                pins_p[3].add_label(output_net, pt(pins_p[3].snap_point))
                pins_p[4].add_label(REF_NET, pt(pins_p[4].snap_point))
            else:
                pins_p[3].add_label(REF_NET, pt(pins_p[3].snap_point))
                pins_p[4].add_label(output_net, pt(pins_p[4].snap_point))
            n_cccs += 1
            prev_inst_phi = inst_p

        transaction.commit()

    design.save_design()
    print(f"Created {n_cccs} CCCS instances ({len(modes)} modes x 2 chains) in "
          f"{library_name}:{new_cell_name}")

    pin_layer_id = db_uu.LayerId(231)

    pin_defs = [
        ("Itheta", (west_x, GRID_START_Y), 180.0, db_uu.TermType.INPUT),
        ("Iphi", (west_x, GRID_START_Y - ROW_PITCH / 2), 180.0, db_uu.TermType.INPUT),
    ]
    for row_idx, (mode_type, l, m, port_name) in enumerate(modes):
        y = GRID_START_Y - (row_idx * ROW_PITCH)
        pin_defs.append((port_name, (east_x, y), 0.0, db_uu.TermType.OUTPUT))
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

    symbol = db_uu.create_symbol(f"{library_name}:{new_cell_name}:symbol")
    symbol_generator = gs.SymbolGenerator(symbol, design, 0.25, 0.25)
    symbol_generator.is_dual_symbol_type = False
    symbol_generator.should_replace = True
    symbol_generator.pin_shape = "dot"
    symbol_generator.order = gs.OrderType.ORDER_LOCATION
    symbol_generator.generate_symbol()
    symbol.save_design()
    print(f"Generated Quad symbol for {library_name}:{new_cell_name}.")
    print("Usage reminder: the parent design must supply a_Rx, c0 as "
          "VAR-scope local variables (for the T delay parameter). "
          "Itheta/Iphi must be driven by currents through EXTERNAL "
          "mu0*r inductors fed from Vtheta/Vphi -- not built inside "
          "this block.")

    return design


print("Ready: create_cccs_rx_aperture(snp_path, ...) is defined. Call it "
      "with your Touchstone file path to build the cell -- see the "
      "module docstring's Usage section.")
