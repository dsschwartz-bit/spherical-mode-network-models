r"""
create_tx_rx_link.py -- build a Tx-to-Rx spherical-mode link schematic in ADS (S-parameters).

create_tx_rx_link(tx_snp_path, rx_snp_path, ...) creates cell "Tx<tx>_Rx<rx>":
    Tx antenna -> eta0 mode resistors -> <tx>_VCVS_TxGain -> VectorFieldRotate_vars
    -> L = mu0*r (theta and phi paths) -> <rx>_CCCS_RxAperture -> eta0 mode resistors
    -> Rx antenna, with Term1/Term2 (Z0) on the Feed ports, a VAR block (DEFAULT_VARS,
    a/Z0 from the Touchstone headers, then var_overrides) and an S_Param controller
    over the overlapping frequency range at the coarser of the two steps.
Each antenna is its GUI-imported rational-fit cell if present, else "<base>_TOUCHESTONE"
(see add_dual_view_snp_instance); all reference pins are grounded. Antennas with
internal ports (e.g. "load1") need create_tx_rx_link_internal_ports.py.

ADS API notes: built-in components bind instance pins by number, the custom subcircuits
by name (label_instance_pins_in_order / _by_name). Every ground uses its own GROUND
symbol, since ADS does not alias a user-named "ref" net to "gnd!".

Requires: VectorFieldRotate_vars / RotationMatrix_vars in Rotation_lib (rotation_library_name).
create_vcvs_tx_gain.py, create_cccs_rx_aperture.py and create_snp_dual_view_cells.py
are loaded automatically from SCRIPT_DIR.

Usage (ADS Python console, after start_console.py):
    run_script("create_tx_rx_link.py")
    create_tx_rx_link(r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p",
                      r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p",
                      library_name="SphericalHarmonics_lib")
"""

import re
from pathlib import Path
import keysight.ads.de as de
from keysight.ads.de import db_uu
from keysight.ads.de.experimental import generate_symbol as gs

# ============================ Fixed configuration ============================

# Repository folders (defined by start_console.py; derived here if this file is exec()'d directly).
try:
    SCRIPT_DIR
except NameError:
    SCRIPT_DIR = Path(THIS_SCRIPT_PATH).resolve().parent  # noqa: F821

# Load the helper scripts (safe to repeat). create_snp_dual_view_cells.py provides
# create_snp_component_cell, redraw_symbol_for_imported_cell and get_gui_import_base_name;
# its extract_antenna_name / get_port_names_from_snp_comments differ, so the vcvs/cccs
# scripts are loaded again afterwards to restore theirs.
exec(open(SCRIPT_DIR / "create_vcvs_tx_gain.py").read())
exec(open(SCRIPT_DIR / "create_cccs_rx_aperture.py").read())
exec(open(SCRIPT_DIR / "create_snp_dual_view_cells.py").read())
exec(open(SCRIPT_DIR / "create_vcvs_tx_gain.py").read())
exec(open(SCRIPT_DIR / "create_cccs_rx_aperture.py").read())

ANTENNA_NAME_PATTERN = re.compile(r"S_dom_(?:psi|Q)_([A-Za-z0-9]+)_")
# Same as the vcvs/cccs scripts' MODE_PATTERN; redefined so it does not depend on
# whichever copy is currently in the shared namespace.
MODE_PATTERN = re.compile(r"^(TM|TE)_L(\d+)_M(neg)?(\d+)$")

# VAR-block defaults (eps0 in pF/m, mu0 in nH/m, r and a in m, angles in rad).
# Override per call with var_overrides.
DEFAULT_VARS = {
    "eps0": "8.8542",
    "mu0": "1256.6",
    "eta0": "376.73",
    "Z0": "eta0",
    "r": "10",
    "theta": "90*(pi/180)",
    "phi": "0*(pi/180)",
    "a_Rx": "0.075",
    "a_Tx": "a_Rx",
    "rot_angle": "0*(pi/180)",
    "Ux": "1",
    "Uy": "1",
    "Uz": "1",
}

PIN_LAYER_ID_NUM = 231
PIN_CLEARANCE = 1.5
BLOCK_SPACING_X = VCVS_WIDTH   # horizontal gap between major blocks (VCVS_WIDTH from create_vcvs_tx_gain.py)
MODE_PIN_PITCH = 2.0     # vertical spacing between stacked per-mode resistors

# ============================================================================


def get_port_names_from_snp_comments(filepath):
    """Return port names from the Touchstone '! Port[N] Name' header comments, in port order."""
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


def extract_antenna_name(snp_path):
    """Return <name> from a filename S_dom_psi_<name>_... or S_dom_Q_<name>_...
    (first alphanumeric word only)."""
    filename = Path(snp_path).name
    match = ANTENNA_NAME_PATTERN.search(filename)
    if not match:
        raise ValueError(
            f"Could not extract antenna name from filename '{filename}' -- "
            f"expected the pattern S_dom_psi_<name>_... or S_dom_Q_<name>_..."
        )
    return match.group(1)


def get_port_count_from_file_extension(snp_path):
    """Get the port count from the .sNNp file extension, e.g. '.s6p' -> 6."""
    extension = Path(snp_path).suffix[1:]
    digits = re.findall(r"\d+", extension)
    if digits:
        return int(digits[0])
    raise ValueError(f"No port number found in file extension: {Path(snp_path).suffix}")


def parse_snp_frequency_range(filepath):
    """Return (start_hz, stop_hz, step_hz) of a Touchstone file, assuming a
    uniform linear sweep.

    Touchstone v1 wraps large S-matrices over continuation lines without a
    frequency, so all numbers are read as one stream and split into blocks
    of 1 + 2*N*N values (N from the file extension). A leading 0 Hz point
    is skipped.
    """
    port_count = get_port_count_from_file_extension(filepath)
    values_per_block = 1 + 2 * port_count * port_count

    unit_multipliers = {"hz": 1.0, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}
    unit_mult = None
    flat_values = []
    with open(filepath, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("!"):
                continue
            if line.startswith("#"):
                tokens = line[1:].split()
                if tokens:
                    unit_mult = unit_multipliers.get(tokens[0].lower())
                continue
            for token in line.split():
                try:
                    flat_values.append(float(token))
                except ValueError:
                    continue

    if len(flat_values) < 2 * values_per_block:
        raise ValueError(
            f"Not enough data in {filepath} for even 2 complete frequency "
            f"blocks (port_count={port_count}, expected {values_per_block} "
            f"values per block, found {len(flat_values)} values total)."
        )
    if len(flat_values) % values_per_block != 0:
        print(f"WARNING: {filepath} has {len(flat_values)} numeric values, "
              f"not an exact multiple of {values_per_block} (expected per "
              f"frequency block for {port_count} ports) -- frequency "
              f"extraction below may be misaligned.")

    unit_mult = unit_mult if unit_mult else 1e9
    freqs = [flat_values[i] * unit_mult for i in range(0, len(flat_values), values_per_block)]

    # A leading DC point (for transient/DC use) would skew the step estimate.
    if freqs and freqs[0] == 0.0:
        print(f"Skipping leading DC (0 Hz) point in {filepath}.")
        freqs = freqs[1:]

    if len(freqs) < 2:
        raise ValueError(f"Fewer than 2 (non-DC) frequency points found in {filepath}")
    start_hz = freqs[0]
    stop_hz = freqs[-1]
    step_hz = (stop_hz - start_hz) / (len(freqs) - 1)
    return start_hz, stop_hz, step_hz


A_Z0_HEADER_PATTERN = re.compile(r"a=([\d.eE+-]+)m.*?Z0=([\d.eE+-]+)ohms", re.IGNORECASE)


def parse_a_and_z0_from_header(filepath):
    """Return (a [m], Z0 [Ohm]) from the 'a=<value>m ... Z0=<value>ohms' header
    comment written by general_vector_spherical_modes.m, or (None, None)."""
    with open(filepath, "r") as f:
        for raw_line in f:
            if not raw_line.startswith("!"):
                break
            match = A_Z0_HEADER_PATTERN.search(raw_line)
            if match:
                return float(match.group(1)), float(match.group(2))
    return None, None


def pt(point):
    """Coerce a PointF-like object (or plain tuple) to a plain (x,y) tuple."""
    return (point.x, point.y) if hasattr(point, "x") else tuple(point)


def label_instance_pins_in_order(instance, net_names):
    """Label pins by term number (net_names[0] -> term 1, ...).

    For built-in components (number-bound terms); use
    label_instance_pins_by_name for the custom subcircuits. Extra pins are
    left unlabeled.
    """
    pins_by_term = {p.inst_term.term_number: p for p in instance.inst_pins}
    for i, net_name in enumerate(net_names, start=1):
        pins_by_term[i].add_label(net_name, pt(pins_by_term[i].snap_point))


def ground_pin(design, pin, name, mirror=None):
    """Place a GROUND symbol on a pin (coincident pins connect without a wire).

    Each grounded pin gets its own GROUND symbol: all GROUND symbols share
    "gnd!", whereas a user-named "ref" net cannot be aliased to it.
    mirror: optional db_uu.MirrorType, to match a mirrored parent (e.g. Term2).
    """
    p = pt(pin.snap_point)
    gnd = design.add_instance(("ads_rflib", "GROUND", "symbol"), p, name=name, angle=-90, mirror=mirror)
    return gnd


def label_instance_pins_by_name(design, instance, name_to_net=None, ground_names=None):
    """Label pins of a custom subcircuit by term name.

    The custom subcircuits (VCVS_Tx_Gain, CCCS_Rx_Aperture,
    VectorFieldRotate_vars) have name-bound terms; accessing them by
    term_number raises "This InstTerm is referenced by name".

    Args:
        name_to_net: pin name -> net label.
        ground_names: pin names to ground with their own GROUND symbol.
    Pins in neither are left unlabeled with a warning.
    """
    name_to_net = name_to_net or {}
    ground_names = ground_names or []
    for p in instance.inst_pins:
        term_name = p.inst_term.term_name
        if term_name in ground_names:
            ground_pin(design, p, f"GND_{instance.name}_{term_name}")
        elif term_name in name_to_net:
            p.add_label(name_to_net[term_name], pt(p.snap_point))
        else:
            print(f"WARNING: pin '{term_name}' on instance '{instance.name}' "
                  f"not found in name_to_net or ground_names -- left unlabeled.")


def add_dual_view_snp_instance(design, filepath, port_names, origin, name, side_tag, library_name, already_redrawn=None, mirror=None):
    """Place the antenna cell: the GUI-imported rational-fit cell
    (get_gui_import_base_name) if it exists in library_name, else
    "<base>_TOUCHESTONE".

    The TOUCHESTONE cell is rebuilt and an existing rational-fit cell's
    symbol redrawn once per run (tracked in already_redrawn, shared by the
    Tx and Rx calls). All reference pins are grounded (N+1 or 2N pin
    scheme). Signal pins are labelled "Feed<side_tag>" or
    "<side_tag>_<port>"; non-mode ports such as "load1" are not terminated.

    Returns (snp, sig_pins_by_port).
    """
    base_name = get_gui_import_base_name(filepath)
    rational_fit_cell_name = base_name
    snp_cell_name = f"{base_name}_TOUCHESTONE"
    ws = de.active_workspace()
    library = ws.open_library(library_name)
    valid_pin_counts = (len(port_names) + 1, 2 * len(port_names))

    snp_redraw_key = (library_name, snp_cell_name)
    if already_redrawn is not None and snp_redraw_key in already_redrawn:
        print(f"'{snp_cell_name}' already rebuilt earlier in this run, not rebuilding again.")
    else:
        print(f"Building/rebuilding '{snp_cell_name}' in '{library_name}' to keep it current.")
        create_snp_component_cell(filepath, library_name=library_name, new_cell_name=snp_cell_name)
        if already_redrawn is not None:
            already_redrawn.add(snp_redraw_key)

    rational_fit_exists = library.cell_exists(rational_fit_cell_name)
    if rational_fit_exists:
        check_design = de.db_uu.open_design((library_name, rational_fit_cell_name, "schematic"))
        existing_pin_count = len(check_design.terms)
        check_design.close_design()
        if existing_pin_count not in valid_pin_counts:
            raise RuntimeError(
                f"'{rational_fit_cell_name}' in '{library_name}' has {existing_pin_count} pins, but "
                f"{valid_pin_counts[0]} (SingRefPin) or {valid_pin_counts[1]} (OneRefPinPerPort) are "
                f"expected for {len(port_names)} ports -- rebuild it via the hybrid workflow "
                f"(GUI import + redraw_symbol_for_imported_cell), or delete it to fall back to "
                f"'{snp_cell_name}'."
            )
        rf_redraw_key = (library_name, rational_fit_cell_name)
        if already_redrawn is not None and rf_redraw_key in already_redrawn:
            print(f"Found existing '{rational_fit_cell_name}' -- already redrawn this run.")
        else:
            print(f"Found existing '{rational_fit_cell_name}' -- redrawing its symbol to keep it current.")
            redraw_symbol_for_imported_cell(library_name, rational_fit_cell_name, snp_path=filepath)
            if already_redrawn is not None:
                already_redrawn.add(rf_redraw_key)

    cell_name_to_use = rational_fit_cell_name if rational_fit_exists else snp_cell_name
    print(f"Using '{cell_name_to_use}' for the {side_tag} side.")

    snp = design.add_instance((library_name, cell_name_to_use, "symbol"), origin, name=name, mirror=mirror)

    # Pins are bound by NUMBER: SingRefPin -> signals 1..N, shared ref N+1;
    # OneRefPinPerPort -> signal 2i-1, ref 2i for port i.
    pins_by_term = {p.inst_term.term_number: p for p in snp.inst_pins}
    sig_pins_by_port = {}
    n = len(port_names)
    if len(snp.inst_pins) == n + 1:
        for i, port_name in enumerate(port_names, start=1):
            sig_pins_by_port[port_name] = pins_by_term[i]
        ground_pin(design, pins_by_term[n + 1], f"GND_{name}_ref")
    elif len(snp.inst_pins) == 2 * n:
        for i, port_name in enumerate(port_names, start=1):
            sig_pins_by_port[port_name] = pins_by_term[2 * i - 1]
            ground_pin(design, pins_by_term[2 * i], f"GND_{name}_ref_{port_name}")
    else:
        raise ValueError(
            f"'{cell_name_to_use}' has {len(snp.inst_pins)} pins, but expected {n + 1} "
            f"(SingRefPin) or {2 * n} (OneRefPinPerPort) for {n} ports."
        )
    for port_name, sig_pin in sig_pins_by_port.items():
        net_label = f"{port_name}{side_tag}" if port_name == "Feed" else f"{side_tag}_{port_name}"
        sig_pin.add_label(net_label, pt(sig_pin.snap_point))
    return snp, sig_pins_by_port


def create_tx_rx_link(
    tx_snp_path,
    rx_snp_path,
    library_name="SphericalHarmonics_lib",
    rotation_library_name="Rotation_lib",
    vector_field_rotate_cell_name="VectorFieldRotate_vars",
    new_cell_name=None,
    var_overrides=None,
    freq_units="GHz",
):
    """Build the Tx-to-Rx link schematic (topology in the module docstring).

    Args:
        tx_snp_path, rx_snp_path: Tx and Rx Touchstone files.
        library_name: library for the new cell and the gain/aperture/antenna cells.
        rotation_library_name: library holding VectorFieldRotate_vars (default Rotation_lib).
        vector_field_rotate_cell_name: name of that cell.
        new_cell_name: default "Tx<tx_name>_Rx<rx_name>".
        var_overrides: dict overriding DEFAULT_VARS / header values by name.
        freq_units: units for the S_Param Start/Stop/Step fields.
    Returns:
        The schematic Design.
    """
    tx_name = extract_antenna_name(tx_snp_path)
    rx_name = extract_antenna_name(rx_snp_path)
    if new_cell_name is None:
        new_cell_name = f"Tx{tx_name}_Rx{rx_name}"
        print(f"Tx antenna '{tx_name}', Rx antenna '{rx_name}' -- using "
              f"cell name '{new_cell_name}'.")

    # ---- Build (or rebuild) the Tx gain and Rx aperture cells ----
    create_vcvs_tx_gain(tx_snp_path, library_name=library_name)
    create_cccs_rx_aperture(rx_snp_path, library_name=library_name)
    vcvs_tx_gain_cell = f"{tx_name}_VCVS_TxGain"
    cccs_rx_aperture_cell = f"{rx_name}_CCCS_RxAperture"

    # ---- Parse both Touchstone files (ports + frequency range) ----
    tx_port_names = get_port_names_from_snp_comments(tx_snp_path)
    rx_port_names = get_port_names_from_snp_comments(rx_snp_path)
    tx_start, tx_stop, tx_step = parse_snp_frequency_range(tx_snp_path)
    rx_start, rx_stop, rx_step = parse_snp_frequency_range(rx_snp_path)
    combined_start = max(tx_start, rx_start)
    combined_stop = min(tx_stop, rx_stop)
    combined_step = max(tx_step, rx_step)
    print(f"Tx sweep: {tx_start/1e9:.4g}-{tx_stop/1e9:.4g} GHz, step {tx_step/1e9:.4g} GHz")
    print(f"Rx sweep: {rx_start/1e9:.4g}-{rx_stop/1e9:.4g} GHz, step {rx_step/1e9:.4g} GHz")
    print(f"Combined (overlap range, coarser step): "
          f"{combined_start/1e9:.4g}-{combined_stop/1e9:.4g} GHz, step {combined_step/1e9:.4g} GHz")
    if combined_start >= combined_stop:
        raise ValueError("Tx and Rx frequency ranges do not overlap.")

    # ---- a/Z0 from the file headers override DEFAULT_VARS; var_overrides wins ----
    var_values = dict(DEFAULT_VARS)

    tx_a, tx_z0 = parse_a_and_z0_from_header(tx_snp_path)
    if tx_a is not None:
        var_values["a_Tx"] = str(tx_a)
        print(f"Tx file header: a_Tx = {tx_a} m")
    rx_a, rx_z0 = parse_a_and_z0_from_header(rx_snp_path)
    if rx_a is not None:
        var_values["a_Rx"] = str(rx_a)
        print(f"Rx file header: a_Rx = {rx_a} m")

    # Z0: Tx value takes precedence (warn on mismatch); written as "eta0" if equal to it.
    parsed_z0 = tx_z0 if tx_z0 is not None else rx_z0
    if tx_z0 is not None and rx_z0 is not None and abs(tx_z0 - rx_z0) > 0.01:
        print(f"WARNING: Tx file Z0={tx_z0} and Rx file Z0={rx_z0} disagree -- "
              f"using Tx's value ({tx_z0}).")
    if parsed_z0 is not None:
        if abs(parsed_z0 - 376.73) < 0.01:
            var_values["Z0"] = "eta0"
        else:
            var_values["Z0"] = str(parsed_z0)
        print(f"Header Z0 = {parsed_z0} Ohm -> Z0 = {var_values['Z0']}")

    if var_overrides:
        var_values.update(var_overrides)

    # ---- Get workspace/libraries, open required existing designs ----
    ws = de.active_workspace()
    print(f"Active workspace: {ws}")
    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    ws.open_library(library_name)
    if rotation_library_name not in ws.library_names:
        raise RuntimeError(f"Library '{rotation_library_name}' not found.")
    ws.open_library(rotation_library_name)

    vfr_design = de.db_uu.open_design(
        (rotation_library_name, vector_field_rotate_cell_name, "schematic")
    )
    vcvs_design = de.db_uu.open_design((library_name, vcvs_tx_gain_cell, "schematic"))
    cccs_design = de.db_uu.open_design((library_name, cccs_rx_aperture_cell, "schematic"))

    # ---- Create the new schematic ----
    design = de.db_uu.create_schematic((library_name, new_cell_name, "schematic"))

    tx_n = len(tx_port_names)
    rx_n = len(rx_port_names)

    x_tx_snp = 0.0
    x_vcvs = BLOCK_SPACING_X + 1
    x_vfr = 7.5
    x_l = 8.0
    x_cccs = x_l + BLOCK_SPACING_X + 0.5
    x_rx_snp = 15.0
    y_top = 0.0

    # VAR1/SP1 row above the main chain (approximate; bounding boxes not queried).
    y_above = y_top + 10.0
    x_var = x_tx_snp
    x_sp = x_vcvs

    # Top of the resistor columns below the antenna blocks (approximate).
    y_below_snp = -1.0

    with de.db_uu.Transaction(design) as transaction:
        # ---- VAR block (above the main chain) ----
        var_inst = design.add_instance(
            ("ads_datacmps", "VAR", "symbol"), (x_var, y_above), name="VAR1", angle=90
        )
        for name, value in var_values.items():
            var_inst.vars[name] = value
        del var_inst.vars["X"]

        # ---- Tx / Rx antenna cells (Rational Fit if present, else TOUCHESTONE; refs grounded) ----
        already_redrawn_snp_cells = set()
        tx_snp, tx_sig_pins = add_dual_view_snp_instance(design, tx_snp_path, tx_port_names, (x_tx_snp, y_top), "SnP1", side_tag="Tx", library_name=library_name, already_redrawn=already_redrawn_snp_cells)
        rx_snp, rx_sig_pins = add_dual_view_snp_instance(design, rx_snp_path, rx_port_names, (x_rx_snp, y_top), "SnP2", side_tag="Rx", library_name=library_name, already_redrawn=already_redrawn_snp_cells, mirror=db_uu.MirrorType.MIRROR_Y)

        # ---- Per-mode eta0 resistors below each antenna block ----
        # Only TM_/TE_ mode ports (MODE_PATTERN); Feed and internal ports are skipped.
        # R has term1 at local (0,0), term2 at (1,0); angle=-90 puts the grounded
        # term2 directly below term1. Tx_/Rx_ prefixes keep the two sides' nets
        # separate (e.g. in a self-link).
        for i, name in enumerate(tx_port_names):
            if not MODE_PATTERN.match(name):
                continue
            r_inst = design.add_instance(
                ("ads_rflib", "R", "symbol"), (x_tx_snp, y_below_snp - i * MODE_PIN_PITCH), name=f"R_Tx_{name}", angle=-90
            )
            r_inst.parameters["R"].value = "eta0 Ohm"
            r_pins = {p.inst_term.term_number: p for p in r_inst.inst_pins}
            r_pins[1].add_label(f"Tx_{name}", pt(r_pins[1].snap_point))
            ground_pin(design, r_pins[2], f"GND_R_Tx_{name}")
        for i, name in enumerate(rx_port_names):
            if not MODE_PATTERN.match(name):
                continue
            r_inst = design.add_instance(
                ("ads_rflib", "R", "symbol"), (x_rx_snp, y_below_snp - i * MODE_PIN_PITCH), name=f"R_Rx_{name}", angle=-90
            )
            r_inst.parameters["R"].value = "eta0 Ohm"
            r_pins = {p.inst_term.term_number: p for p in r_inst.inst_pins}
            r_pins[1].add_label(f"Rx_{name}", pt(r_pins[1].snap_point))
            ground_pin(design, r_pins[2], f"GND_R_Rx_{name}")

        # ---- Term1 (Tx feed), Term2 (Rx feed), placed on the antenna Feed pins ----
        # Term's signal pin is at its origin, so it coincides with the Feed pin;
        # its reference pin gets its own GROUND.
        tx_feed_pos = pt(tx_sig_pins["Feed"].snap_point)
        term1 = design.add_instance(("ads_simulation", "Term", "symbol"), tx_feed_pos, name="Term1", angle=180)
        term1.parameters["Num"].value = "1"
        term1.parameters["Z"].value = "Z0 Ohm"
        term1_pins = {p.inst_term.term_number: p for p in term1.inst_pins}
        term1_pins[1].add_label("FeedTx", pt(term1_pins[1].snap_point))
        ground_pin(design, term1_pins[2], "GND_Term1")

        rx_feed_pos = pt(rx_sig_pins["Feed"].snap_point)
        # Term pins: pin1 at local (0,0), pin2 at (1,0). MIRROR_Y plus angle=180
        # keeps pin2 (ground) pointing away from the circuit; pin1 stays on Feed.
        term2 = design.add_instance(("ads_simulation", "Term", "symbol"), rx_feed_pos, name="Term2", angle=180, mirror=db_uu.MirrorType.MIRROR_Y)
        term2.parameters["Num"].value = "2"
        term2.parameters["Z"].value = "Z0 Ohm"
        term2_pins = {p.inst_term.term_number: p for p in term2.inst_pins}
        term2_pins[1].add_label("FeedRx", pt(term2_pins[1].snap_point))
        ground_pin(design, term2_pins[2], "GND_Term2", mirror=db_uu.MirrorType.MIRROR_Y)

        # ---- VCVS_Tx_Gain instance: mode inputs (matching Tx SnP order), Vtheta_prime, Vphi_prime, ref (grounded) ----
        vcvs_inst = design.add_instance((library_name, vcvs_tx_gain_cell, "symbol"), (x_vcvs, y_top), name="VCVS_TxGain_inst", angle=0)
        vcvs_name_to_net = {name: f"Tx_{name}" for name in tx_port_names}
        vcvs_name_to_net.update({"Vtheta_prime": "Vtheta_prime", "Vphi_prime": "Vphi_prime"})
        label_instance_pins_by_name(design, vcvs_inst, vcvs_name_to_net, ground_names=["ref"])

        # ---- VectorFieldRotate_vars instance ----
        vfr_inst = design.add_instance((rotation_library_name, vector_field_rotate_cell_name, "symbol"), (x_vfr, y_top), name="VectorFieldRotate_inst", angle=0)
        # Vr_prime/Vr are grounded (gradYlm has no radial component), as is ref.
        vfr_name_to_net = {"Vtheta_prime": "Vtheta_prime", "Vphi_prime": "Vphi_prime", "Vtheta": "Vtheta", "Vphi": "Vphi"}
        label_instance_pins_by_name(design, vfr_inst, vfr_name_to_net, ground_names=["Vr_prime", "Vr", "ref"])

        # ---- Two inductors: Vtheta->Itheta, Vphi->Iphi ----
        l_theta = design.add_instance(("ads_rflib", "L", "symbol"), (x_l, y_top + 1.5), name="L1", angle=0)
        l_theta.parameters["L"].value = "(mu0*r) nH"
        label_instance_pins_in_order(l_theta, ["Vtheta", "Itheta"])

        l_phi = design.add_instance(("ads_rflib", "L", "symbol"), (x_l, y_top - MODE_PIN_PITCH + 1.5), name="L2", angle=0)
        l_phi.parameters["L"].value = "(mu0*r) nH"
        label_instance_pins_in_order(l_phi, ["Vphi", "Iphi"])

        # ---- CCCS_Rx_Aperture instance: Itheta, Iphi, mode outputs (matching Rx SnP order), ref (grounded) ----
        cccs_inst = design.add_instance((library_name, cccs_rx_aperture_cell, "symbol"), (x_cccs, y_top), name="CCCS_RxAperture_inst", angle=0)
        cccs_name_to_net = {name: f"Rx_{name}" for name in rx_port_names}
        cccs_name_to_net.update({"Itheta": "Itheta", "Iphi": "Iphi"})
        label_instance_pins_by_name(design, cccs_inst, cccs_name_to_net, ground_names=["ref"])

        # ---- S_Param simulation controller (above the main chain) ----
        freq_divisors = {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}
        divide_by = freq_divisors[freq_units]
        sp = design.add_instance(("ads_simulation", "S_Param", "symbol"), (x_sp, y_above), name="SP1", angle=0)
        sp.parameters["Start"].value = f"{combined_start / divide_by} {freq_units}"
        sp.parameters["Stop"].value = f"{combined_stop / divide_by} {freq_units}"
        sp.parameters["Step"].value = f"{combined_step / divide_by} {freq_units}"
        # SPSS_WSP is the "Perform Stability Analysis" option; off to save time.
        sp.parameters["SPSS_WSP"].value = "no"

        transaction.commit()

    design.save_design()
    print(f"Created {new_cell_name} in {library_name}: Tx SnP ({tx_n} ports), "
          f"Rx SnP ({rx_n} ports), {vcvs_tx_gain_cell}, VectorFieldRotate_vars, "
          f"2 inductors, {cccs_rx_aperture_cell}, {tx_n - 1 + rx_n - 1} eta0 "
          f"resistors, Term1/Term2, S_Param controller.")
    print("Usage reminder: all of eps0, mu0, eta0, Z0, r, theta, phi, a_Tx, "
          "a_Rx, rot_angle, Ux, Uy, Uz are set in this cell's own VAR1 block "
          "(global scope); override them with var_overrides.")

    return design


print("Ready: create_tx_rx_link(tx_snp_path, rx_snp_path, ...) is defined "
      "(create_vcvs_tx_gain/create_cccs_rx_aperture auto-loaded above). "
      "See the module docstring's Usage section.")
