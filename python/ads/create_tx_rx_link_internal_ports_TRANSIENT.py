r"""
create_tx_rx_link_internal_ports_TRANSIENT.py -- Tx-to-Rx link with a transient setup.

As create_tx_rx_link_internal_ports.py (see its docstring), plus a Tran controller
next to the S_Param controller. Differences:
  - Causal delays: the whole propagation delay r/c0 (reduced by a_Tx and/or a_Rx for
    antennas whose Touchstone "ka_factor=" is not 1) is put on the Tx VCVS block and
    the Rx CCCS block gets delay 0, so no block needs a negative delay.
  - Excitation: a VtSine source in series with Term1 (Term1 "-" -> VtSine -> Tx Feed
    reference). Its waveform parameters are left at ADS defaults.
  - Tran1: StopTime 10 ns, MaxTimeStep = 1/(2*f_stop)/6 (12 points per period at
    f_stop), MaxIters 20, ImpMaxFreq = f_stop, no passivity enforcement; saved nodes
    Vtheta, Vphi, FeedRx and the reference nets; branch currents off.
The default cell name is "Tx<tx>_Rx<rx>_internal_ports_transient"; an existing cell of
that name is replaced.

Requires: VectorFieldRotate_vars / RotationMatrix_vars in Rotation_lib (rotation_library_name).
create_vcvs_tx_gain.py, create_cccs_rx_aperture.py and create_snp_dual_view_cells.py
are loaded automatically from SCRIPT_DIR.

Usage (ADS Python console, after start_console.py):
    run_script("create_tx_rx_link_internal_ports_TRANSIENT.py")
    create_tx_rx_link_internal_ports_TRANSIENT(
        r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p",
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
# its extract_antenna_name (full stem) differs from the vcvs/cccs scripts' version
# (first word only, used for tx_name/rx_name and the VCVS/CCCS cell names), so those
# two scripts are loaded again afterwards to restore theirs. get_ka_factor_is_one
# also comes from the vcvs/cccs scripts.
exec(open(SCRIPT_DIR / "create_vcvs_tx_gain.py").read())
exec(open(SCRIPT_DIR / "create_cccs_rx_aperture.py").read())
exec(open(SCRIPT_DIR / "create_snp_dual_view_cells.py").read())
exec(open(SCRIPT_DIR / "create_vcvs_tx_gain.py").read())
exec(open(SCRIPT_DIR / "create_cccs_rx_aperture.py").read())

ANTENNA_NAME_PATTERN = re.compile(r"S_dom_(?:psi|Q)_([A-Za-z0-9]+)_")
MODE_PATTERN = re.compile(r"^(TM|TE)_L(\d+)_M(neg)?(\d+)$")

# VAR-block defaults (eps0 in pF/m, mu0 in nH/m, r and a in m, angles in rad).
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
BLOCK_SPACING_X = VCVS_WIDTH
MODE_PIN_PITCH = 2.0

# ============================================================================


def parse_snp_frequency_range(filepath):
    """Return (start_hz, stop_hz, step_hz); same as in create_tx_rx_link.py."""
    port_count = get_port_count_from_file_extension_internal(filepath)
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

    if freqs and freqs[0] == 0.0:
        print(f"Skipping leading DC (0 Hz) point in {filepath}.")
        freqs = freqs[1:]

    if len(freqs) < 2:
        raise ValueError(f"Fewer than 2 (non-DC) frequency points found in {filepath}")
    start_hz = freqs[0]
    stop_hz = freqs[-1]
    step_hz = (stop_hz - start_hz) / (len(freqs) - 1)
    return start_hz, stop_hz, step_hz


def get_port_count_from_file_extension_internal(snp_path):
    extension = Path(snp_path).suffix[1:]
    digits = re.findall(r"\d+", extension)
    if digits:
        return int(digits[0])
    raise ValueError(f"No port number found in file extension: {Path(snp_path).suffix}")


A_Z0_HEADER_PATTERN = re.compile(r"a=([\d.eE+-]+)m.*?Z0=([\d.eE+-]+)ohms", re.IGNORECASE)


def parse_a_and_z0_from_header(filepath):
    """Return (a [m], Z0 [Ohm]) from the header, or (None, None); same as in create_tx_rx_link.py."""
    with open(filepath, "r") as f:
        for raw_line in f:
            if not raw_line.startswith("!"):
                break
            match = A_Z0_HEADER_PATTERN.search(raw_line)
            if match:
                return float(match.group(1)), float(match.group(2))
    return None, None


def pt(point):
    return (point.x, point.y) if hasattr(point, "x") else tuple(point)


def label_instance_pins_in_order(instance, net_names):
    pins_by_term = {p.inst_term.term_number: p for p in instance.inst_pins}
    for i, net_name in enumerate(net_names, start=1):
        pins_by_term[i].add_label(net_name, pt(pins_by_term[i].snap_point))


def ground_pin(design, pin, name, mirror=None):
    p = pt(pin.snap_point)
    gnd = design.add_instance(("ads_rflib", "GROUND", "symbol"), p, name=name, angle=-90, mirror=mirror)
    return gnd


def label_instance_pins_by_name(design, instance, name_to_net=None, ground_names=None):
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


def add_snp_instance_internal_ports(design, filepath, port_names, max_freq_hz, origin, name, side_tag, mirror=None):
    """Place a raw native SnP instance (OneRefPinPerPort, ImpMaxFreq = max_freq_hz).

    Not used by create_tx_rx_link_internal_ports_TRANSIENT, which uses
    add_dual_view_snp_instance_internal_ports; kept as a standalone
    alternative. Port i has signal term 2i-1 and reference term 2i; mode
    references are grounded, Feed and internal-port references left for
    the caller. Returns (snp, ref_pins_by_port).
    """
    snp = design.add_instance(
        db_uu.LCVName("ads_datacmps", "SnP", "symbol"), origin, name=name, mirror=mirror
    )
    port_count = len(port_names)
    longest = max(port_names, key=len)
    port_names_string = f"1 {len(longest)} {' '.join(port_names)}"
    max_freq_ghz = max_freq_hz / 1e9
    with de.db.Transaction(design, "Edit params") as transaction:
        snp.parameters["File"].value = f'"{filepath}"'
        snp.parameters["NumPorts"].value = str(port_count)
        snp.parameters["Type"].value = '"touchstone"'
        snp.parameters["PinSpacing"].value = "Loose"
        snp.parameters["HasReferencePin"].value = "OneRefPinPerPort"
        snp.parameters["PinLayout"].value = "DualInLine"
        snp.parameters["InterpMode"].value = '"linear"'
        snp.parameters["ImpMaxFreq"].value = f"{max_freq_ghz} GHz"
        snp.parameters["EnforcePassivity"].value = "no"
        snp.parameters["port_name_list"].value = port_names_string
        snp.invoke_item_parameter_changed_callback(
            ["File", "NumPorts", "Type", "PinSpacing", "HasReferencePin",
             "PinLayout", "InterpMode", "ImpMaxFreq", "EnforcePassivity", "port_name_list"]
        )
        transaction.commit()

    pins_by_term = {p.inst_term.term_number: p for p in snp.inst_pins}
    ref_pins_by_port = {}

    for i, port_name in enumerate(port_names, start=1):
        sig_term = 2 * i - 1
        ref_term = 2 * i
        sig_pin = pins_by_term[sig_term]
        ref_pin = pins_by_term[ref_term]
        ref_pins_by_port[port_name] = ref_pin

        net_label = f"{port_name}{side_tag}" if port_name == "Feed" else f"{side_tag}_{port_name}"
        sig_pin.add_label(net_label, pt(sig_pin.snap_point))

        if port_name == "Feed":
            # Wired by the caller (to VtSine on Tx, to Term2 on Rx).
            pass
        elif MODE_PATTERN.match(port_name):
            ground_pin(design, ref_pin, f"GND_{name}_ref_{port_name}")
        # else: internal port (e.g. "load1") -- left for the caller.

    return snp, ref_pins_by_port


def add_dual_view_snp_instance_internal_ports(design, filepath, port_names, origin, name, side_tag, library_name, already_redrawn=None, mirror=None):
    """Place the antenna cell and label/ground its pins per port type.

    Same as in create_tx_rx_link_internal_ports.py: uses the GUI-imported
    rational-fit cell named get_gui_import_base_name(filepath) if it exists
    in library_name, else "<base>_TOUCHESTONE" (always rebuilt); an existing
    rational-fit cell has its pin count checked and its symbol redrawn.
    The SnP parameters (file, ports, ...) are set inside the cell; the
    TOUCHESTONE cell does not set ImpMaxFreq (Tran1 sets it instead).

    Reference scheme, detected from the pin count:
      - N+1 pins (SingRefPin): the shared reference is returned for every
        port and not grounded (the caller wires it, as Feed's reference,
        to VtSine). Not allowed when the antenna has internal ports: each
        internal port needs its own reference, so a ValueError is raised
        (use a OneRefPinPerPort cell, e.g. "<base>_TOUCHESTONE").
      - 2N pins (OneRefPinPerPort): mode references grounded; Feed and
        internal-port references left for the caller.

    Args:
        already_redrawn: set shared by the Tx and Rx calls so each cell is
            rebuilt/redrawn only once per run; redrawing a symbol that is
            already instantiated leaves the earlier instance's pins misplaced.
    Returns:
        (snp, ref_pins_by_port, sig_pins_by_port)
    """
    base_name = get_gui_import_base_name(filepath)
    rational_fit_cell_name = base_name
    snp_cell_name = f"{base_name}_TOUCHESTONE"
    ws = de.active_workspace()
    library = ws.open_library(library_name)
    # A rational-fit cell may use either reference scheme.
    valid_pin_counts = (len(port_names) + 1, 2 * len(port_names))

    # ---- <base>_TOUCHESTONE: rebuilt via create_snp_component_cell ----
    snp_redraw_key = (library_name, snp_cell_name)
    if already_redrawn is not None and snp_redraw_key in already_redrawn:
        print(f"'{snp_cell_name}' in '{library_name}' -- already "
              f"rebuilt earlier in this run (e.g. Tx and Rx sharing "
              f"the same antenna), not rebuilding again.")
    else:
        print(f"Building/rebuilding '{snp_cell_name}' in "
              f"'{library_name}' via create_snp_component_cell to "
              f"keep it current.")
        create_snp_component_cell(filepath, library_name=library_name, new_cell_name=snp_cell_name)
        if already_redrawn is not None:
            already_redrawn.add(snp_redraw_key)

    # ---- Rational-fit cell (name == base_name): symbol redrawn if it exists ----
    rational_fit_exists = library.cell_exists(rational_fit_cell_name)
    if rational_fit_exists:
        check_design = de.db_uu.open_design((library_name, rational_fit_cell_name, "schematic"))
        existing_pin_count = len(check_design.terms)
        check_design.close_design()
        if existing_pin_count not in valid_pin_counts:
            raise RuntimeError(
                f"'{rational_fit_cell_name}' in '{library_name}' has "
                f"{existing_pin_count} pins, but {valid_pin_counts[0]} "
                f"(SingRefPin) or {valid_pin_counts[1]} (OneRefPinPerPort) "
                f"are expected for {len(port_names)} ports -- its SnP "
                f"file must have changed since it was last built. Rebuild "
                f"it via the hybrid workflow (GUI import + "
                f"redraw_symbol_for_imported_cell) before re-running "
                f"this, or delete it to fall back to "
                f"'{snp_cell_name}' instead."
            )
        rf_redraw_key = (library_name, rational_fit_cell_name)
        if already_redrawn is not None and rf_redraw_key in already_redrawn:
            print(f"Found existing '{rational_fit_cell_name}' in "
                  f"'{library_name}' -- already redrawn earlier in this "
                  f"run, not redrawing again.")
        else:
            print(f"Found existing '{rational_fit_cell_name}' in "
                  f"'{library_name}' -- redrawing its symbol to keep it "
                  f"current.")
            redraw_symbol_for_imported_cell(library_name, rational_fit_cell_name, snp_path=filepath)
            if already_redrawn is not None:
                already_redrawn.add(rf_redraw_key)

    # ---- Place the rational-fit cell if present ----
    cell_name_to_use = rational_fit_cell_name if rational_fit_exists else snp_cell_name
    print(f"Using '{cell_name_to_use}' for this instance.")

    snp = design.add_instance(
        (library_name, cell_name_to_use, "symbol"), origin, name=name, mirror=mirror
    )

    # Terms are number-bound (term_name raises). SingRefPin: signals 1..N,
    # shared ref N+1. OneRefPinPerPort: port i has signal 2i-1, ref 2i.
    pins_by_term = {p.inst_term.term_number: p for p in snp.inst_pins}
    ref_pins_by_port = {}
    sig_pins_by_port = {}

    if len(snp.inst_pins) == len(port_names) + 1:
        # SingRefPin: one shared reference, not grounded; the caller wires it
        # (as Feed's reference) to VtSine. With internal ports this would tie
        # every port reference to one node, so it is refused.
        internal = [p for p in port_names if p != "Feed" and not MODE_PATTERN.match(p)]
        if internal:
            raise ValueError(
                f"'{cell_name_to_use}' has a single shared reference pin, but "
                f"{side_tag} has internal ports {internal}: each port needs its "
                f"own reference. Use a OneRefPinPerPort cell (delete the "
                f"rational-fit cell to fall back to '{snp_cell_name}', or "
                f"re-import the fit with one reference pin per port).")
        shared_ref_pin = pins_by_term[len(port_names) + 1]
        for i, port_name in enumerate(port_names, start=1):
            sig_pin = pins_by_term[i]
            sig_pins_by_port[port_name] = sig_pin
            ref_pins_by_port[port_name] = shared_ref_pin
            net_label = f"{port_name}{side_tag}" if port_name == "Feed" else f"{side_tag}_{port_name}"
            sig_pin.add_label(net_label, pt(sig_pin.snap_point))
    elif len(snp.inst_pins) == 2 * len(port_names):
        # OneRefPinPerPort: one reference pin per port.
        for i, port_name in enumerate(port_names, start=1):
            sig_pin = pins_by_term[2 * i - 1]
            ref_pin = pins_by_term[2 * i]
            sig_pins_by_port[port_name] = sig_pin
            ref_pins_by_port[port_name] = ref_pin

            net_label = f"{port_name}{side_tag}" if port_name == "Feed" else f"{side_tag}_{port_name}"
            sig_pin.add_label(net_label, pt(sig_pin.snap_point))

            if port_name == "Feed":
                # Wired by the caller (to VtSine on Tx, to Term2 on Rx).
                pass
            elif MODE_PATTERN.match(port_name):
                ground_pin(design, ref_pin, f"GND_{name}_ref_{port_name}")
            # else: internal port (e.g. "load1") -- the caller adds a Term.
    else:
        raise ValueError(
            f"'{cell_name_to_use}' has {len(snp.inst_pins)} pins, but "
            f"expected {len(port_names) + 1} (SingRefPin) or "
            f"{2 * len(port_names)} (OneRefPinPerPort) for "
            f"{len(port_names)} ports."
        )

    return snp, ref_pins_by_port, sig_pins_by_port


def create_tx_rx_link_internal_ports_TRANSIENT(
    tx_snp_path,
    rx_snp_path,
    library_name="SphericalHarmonics_lib",
    rotation_library_name="Rotation_lib",
    vector_field_rotate_cell_name="VectorFieldRotate_vars",
    new_cell_name=None,
    var_overrides=None,
    freq_units="GHz",
):
    """Build the Tx-to-Rx link schematic with VtSine excitation, causal
    delays and a Tran controller (see module docstring).

    Args: as for create_tx_rx_link(); new_cell_name defaults to
    "Tx<tx_name>_Rx<rx_name>_internal_ports_transient".
    Returns:
        The schematic Design.
    """
    tx_name = extract_antenna_name(tx_snp_path)
    rx_name = extract_antenna_name(rx_snp_path)
    if new_cell_name is None:
        new_cell_name = f"Tx{tx_name}_Rx{rx_name}_internal_ports_transient"
        print(f"Tx antenna '{tx_name}', Rx antenna '{rx_name}' -- using "
              f"cell name '{new_cell_name}'.")

    # Put the whole delay on the Tx VCVS (Rx CCCS delay = 0) so no block has a
    # negative delay. The a_Tx / a_Rx term is subtracted only for an antenna
    # whose "ka_factor=" header is not 1 (get_ka_factor_is_one).
    tx_needs_correction = not get_ka_factor_is_one(tx_snp_path)
    rx_needs_correction = not get_ka_factor_is_one(rx_snp_path)
    if tx_needs_correction and rx_needs_correction:
        combined_vcvs_delay = "((r-a_Tx-a_Rx)/c0) sec"
    elif tx_needs_correction:
        combined_vcvs_delay = "((r-a_Tx)/c0) sec"
    elif rx_needs_correction:
        combined_vcvs_delay = "((r-a_Rx)/c0) sec"
    else:
        combined_vcvs_delay = "r/c0 sec"
    print(f"Tx ka_factor=1: {not tx_needs_correction}, Rx ka_factor=1: "
          f"{not rx_needs_correction} -- combined VCVS delay_formula="
          f"'{combined_vcvs_delay}'.")

    create_vcvs_tx_gain(tx_snp_path, library_name=library_name,
                         delay_formula=combined_vcvs_delay)
    create_cccs_rx_aperture(rx_snp_path, library_name=library_name,
                             delay_formula="0")
    vcvs_tx_gain_cell = f"{tx_name}_VCVS_TxGain"
    cccs_rx_aperture_cell = f"{rx_name}_CCCS_RxAperture"

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

    # a/Z0 from the file headers override DEFAULT_VARS; var_overrides wins.
    var_values = dict(DEFAULT_VARS)

    tx_a, tx_z0 = parse_a_and_z0_from_header(tx_snp_path)
    if tx_a is not None:
        var_values["a_Tx"] = str(tx_a)
        print(f"Tx file header: a_Tx = {tx_a} m")
    rx_a, rx_z0 = parse_a_and_z0_from_header(rx_snp_path)
    if rx_a is not None:
        var_values["a_Rx"] = str(rx_a)
        print(f"Rx file header: a_Rx = {rx_a} m")

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

    ws = de.active_workspace()
    print(f"Active workspace: {ws}")
    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    library = ws.open_library(library_name)
    if rotation_library_name not in ws.library_names:
        raise RuntimeError(f"Library '{rotation_library_name}' not found.")
    ws.open_library(rotation_library_name)

    vfr_design = de.db_uu.open_design(
        (rotation_library_name, vector_field_rotate_cell_name, "schematic")
    )
    vcvs_design = de.db_uu.open_design((library_name, vcvs_tx_gain_cell, "schematic"))
    cccs_design = de.db_uu.open_design((library_name, cccs_rx_aperture_cell, "schematic"))

    # Replace an existing cell; its views must be closed first or delete_cell
    # fails with "cell has open cellviews".
    if library.cell_exists(new_cell_name):
        print(f"Cell '{new_cell_name}' already exists in '{library_name}' "
              f"-- deleting it first.")
        for view_name in ("schematic", "symbol"):
            try:
                de.db_uu.open_design((library_name, new_cell_name, view_name)).close_design()
            except Exception:
                pass
        library.delete_cell(new_cell_name)
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

    y_above = y_top + 10.0
    x_var = x_tx_snp
    x_sp = x_vcvs
    x_tran = x_vfr

    y_below_snp = -1.0

    with de.db_uu.Transaction(design) as transaction:
        var_inst = design.add_instance(
            ("ads_datacmps", "VAR", "symbol"), (x_var, y_above), name="VAR1", angle=90
        )
        for var_name_key, value in var_values.items():
            var_inst.vars[var_name_key] = value
        del var_inst.vars["X"]

        # ---- Tx antenna (mode refs grounded; Feed ref via VtSine to Term1; internal refs to Terms) ----
        already_redrawn_snp_cells = set()
        tx_snp, tx_ref_pins_by_port, tx_sig_pins_by_port = add_dual_view_snp_instance_internal_ports(design, tx_snp_path, tx_port_names, (x_tx_snp, y_top), "SnP1", side_tag="Tx", library_name=library_name, already_redrawn=already_redrawn_snp_cells)

        # ---- Rx antenna (Feed ref wired directly to Term2) ----
        rx_snp, rx_ref_pins_by_port, rx_sig_pins_by_port = add_dual_view_snp_instance_internal_ports(design, rx_snp_path, rx_port_names, (x_rx_snp, y_top), "SnP2", side_tag="Rx", library_name=library_name, already_redrawn=already_redrawn_snp_cells, mirror=db_uu.MirrorType.MIRROR_Y)

        # ---- Per-mode eta0 resistors below each antenna block (term2 grounded) ----
        for i, port_name in enumerate(tx_port_names):
            if not MODE_PATTERN.match(port_name):
                continue
            r_inst = design.add_instance(
                ("ads_rflib", "R", "symbol"), (x_tx_snp, y_below_snp - i * MODE_PIN_PITCH), name=f"R_Tx_{port_name}", angle=-90
            )
            r_inst.parameters["R"].value = "eta0 Ohm"
            r_pins = {p.inst_term.term_number: p for p in r_inst.inst_pins}
            r_pins[1].add_label(f"Tx_{port_name}", pt(r_pins[1].snap_point))
            ground_pin(design, r_pins[2], f"GND_R_Tx_{port_name}")
        for i, port_name in enumerate(rx_port_names):
            if not MODE_PATTERN.match(port_name):
                continue
            r_inst = design.add_instance(
                ("ads_rflib", "R", "symbol"), (x_rx_snp, y_below_snp - i * MODE_PIN_PITCH), name=f"R_Rx_{port_name}", angle=-90
            )
            r_inst.parameters["R"].value = "eta0 Ohm"
            r_pins = {p.inst_term.term_number: p for p in r_inst.inst_pins}
            r_pins[1].add_label(f"Rx_{port_name}", pt(r_pins[1].snap_point))
            ground_pin(design, r_pins[2], f"GND_R_Rx_{port_name}")

        tx_feed_pos = pt(tx_sig_pins_by_port["Feed"].snap_point)
        term1 = design.add_instance(("ads_simulation", "Term", "symbol"), tx_feed_pos, name="Term1", angle=180)
        term1.parameters["Num"].value = "1"
        term1.parameters["Z"].value = "Z0 Ohm"
        term1_pins = {p.inst_term.term_number: p for p in term1.inst_pins}
        term1_pins[1].add_label("FeedTx", pt(term1_pins[1].snap_point))

        # VtSine in series with Term1: Term1 "-" -> VtSine -> Tx Feed reference.
        term1_neg_pos = pt(term1_pins[2].snap_point)
        vtsine = design.add_instance(("ads_sources", "VtSine", "symbol"), term1_neg_pos, name="VtSine1", angle=180)
        vtsine_pins = {p.inst_term.term_number: p for p in vtsine.inst_pins}
        design.add_wire([pt(vtsine_pins[2].snap_point), pt(tx_ref_pins_by_port["Feed"].snap_point)])
        tx_ref_pins_by_port["Feed"].add_label("FeedTx_ref", pt(tx_ref_pins_by_port["Feed"].snap_point))

        # ---- Terms on internal (non-Feed, non-mode) ports, e.g. "load1" ----
        # Term numbering: Tx Feed = 1, Tx internal ports = 2..(1+n_tx_int),
        # then Rx Feed = 2+n_tx_int, Rx internal ports continue after it.
        # Each Term sits on its port's signal pin (same placement rule as
        # Feed's Term) with its "-" terminal wired to that port's own ref pin.
        ref_net_names = []  # reference-pin net names (added to Tran outputs)
        def _add_internal_port_terms(snp_port_names, sig_pins, ref_pins, side_tag, first_num, mirror=None):
            num = first_num
            for port_name in snp_port_names:
                if port_name == "Feed" or MODE_PATTERN.match(port_name):
                    continue
                t = design.add_instance(("ads_simulation", "Term", "symbol"), pt(sig_pins[port_name].snap_point),
                                        name=f"Term_{side_tag}_{port_name}", angle=180, mirror=mirror)
                t.parameters["Num"].value = str(num)
                t.parameters["Z"].value = "Z0 Ohm"
                t_pins = {p.inst_term.term_number: p for p in t.inst_pins}
                t_pins[1].add_label(f"{side_tag}_{port_name}", pt(t_pins[1].snap_point))
                design.add_wire([pt(t_pins[2].snap_point), pt(ref_pins[port_name].snap_point)])
                ref_net = f"{side_tag}_{port_name}_ref"
                ref_pins[port_name].add_label(ref_net, pt(ref_pins[port_name].snap_point))
                ref_net_names.append(ref_net)
                num += 1
            return num

        n_tx_int = sum(1 for p in tx_port_names if p != "Feed" and not MODE_PATTERN.match(p))
        _add_internal_port_terms(tx_port_names, tx_sig_pins_by_port, tx_ref_pins_by_port, "Tx", 2)

        rx_feed_pos = pt(rx_sig_pins_by_port["Feed"].snap_point)
        term2 = design.add_instance(("ads_simulation", "Term", "symbol"), rx_feed_pos, name="Term2", angle=180, mirror=db_uu.MirrorType.MIRROR_Y)
        term2.parameters["Num"].value = str(2 + n_tx_int)
        term2.parameters["Z"].value = "Z0 Ohm"
        term2_pins = {p.inst_term.term_number: p for p in term2.inst_pins}
        term2_pins[1].add_label("FeedRx", pt(term2_pins[1].snap_point))
        design.add_wire([pt(term2_pins[2].snap_point), pt(rx_ref_pins_by_port["Feed"].snap_point)])
        rx_ref_pins_by_port["Feed"].add_label("FeedRx_ref", pt(rx_ref_pins_by_port["Feed"].snap_point))
        _add_internal_port_terms(rx_port_names, rx_sig_pins_by_port, rx_ref_pins_by_port, "Rx", 3 + n_tx_int,
                                 mirror=db_uu.MirrorType.MIRROR_Y)

        vcvs_inst = design.add_instance((library_name, vcvs_tx_gain_cell, "symbol"), (x_vcvs, y_top), name="VCVS_TxGain_inst", angle=0)
        vcvs_name_to_net = {port_name: f"Tx_{port_name}" for port_name in tx_port_names}
        vcvs_name_to_net.update({"Vtheta_prime": "Vtheta_prime", "Vphi_prime": "Vphi_prime"})
        label_instance_pins_by_name(design, vcvs_inst, vcvs_name_to_net, ground_names=["ref"])

        vfr_inst = design.add_instance((rotation_library_name, vector_field_rotate_cell_name, "symbol"), (x_vfr, y_top), name="VectorFieldRotate_inst", angle=0)
        vfr_name_to_net = {"Vtheta_prime": "Vtheta_prime", "Vphi_prime": "Vphi_prime", "Vtheta": "Vtheta", "Vphi": "Vphi"}
        label_instance_pins_by_name(design, vfr_inst, vfr_name_to_net, ground_names=["Vr_prime", "Vr", "ref"])

        l_theta = design.add_instance(("ads_rflib", "L", "symbol"), (x_l, y_top + 1.5), name="L1", angle=0)
        l_theta.parameters["L"].value = "(mu0*r) nH"
        label_instance_pins_in_order(l_theta, ["Vtheta", "Itheta"])

        l_phi = design.add_instance(("ads_rflib", "L", "symbol"), (x_l, y_top - MODE_PIN_PITCH + 1.5), name="L2", angle=0)
        l_phi.parameters["L"].value = "(mu0*r) nH"
        label_instance_pins_in_order(l_phi, ["Vphi", "Iphi"])

        cccs_inst = design.add_instance((library_name, cccs_rx_aperture_cell, "symbol"), (x_cccs, y_top), name="CCCS_RxAperture_inst", angle=0)
        cccs_name_to_net = {port_name: f"Rx_{port_name}" for port_name in rx_port_names}
        cccs_name_to_net.update({"Itheta": "Itheta", "Iphi": "Iphi"})
        label_instance_pins_by_name(design, cccs_inst, cccs_name_to_net, ground_names=["ref"])

        freq_divisors = {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}
        divide_by = freq_divisors[freq_units]
        sp = design.add_instance(("ads_simulation", "S_Param", "symbol"), (x_sp, y_above), name="SP1", angle=0)
        sp.parameters["Start"].value = f"{combined_start / divide_by} {freq_units}"
        sp.parameters["Stop"].value = f"{combined_stop / divide_by} {freq_units}"
        sp.parameters["Step"].value = f"{combined_step / divide_by} {freq_units}"
        # SPSS_WSP is the "Perform Stability Analysis" option; off to save time.
        sp.parameters["SPSS_WSP"].value = "no"

        # ---- Tran controller, alongside S_Param (not replacing it) ----
        # MaxTimeStep = 1/(2*f_stop)/6, i.e. 12 points per period at f_stop;
        # MaxIters = 20 Newton iterations per step. TruncTol and ChargeTol stay
        # at ADS defaults.
        max_time_step_sec = 1 / (2 * combined_stop) / 6
        tran_max_iters = 20
        tran = design.add_instance(("ads_simulation", "Tran", "symbol"), (x_tran, y_above), name="Tran1", angle=0)
        tran.parameters["StopTime"].value = "10 nsec"
        tran.parameters["MaxTimeStep"].value = f"{max_time_step_sec} sec"
        tran.parameters["MaxIters"].value = str(tran_max_iters)
        tran.parameters["ImpEnforcePassivity"].value = "no"
        tran.parameters["ImpMaxFreq"].value = f"{combined_stop / divide_by} {freq_units}"
        tran.parameters["UseNodeNestLevel"].value = "no"
        tran.parameters["NodeName"].value = (
            ['"Vtheta"', '"Vphi"', '"FeedRx"', '"FeedTx_ref"', '"FeedRx_ref"']
            + [f'"{nm}"' for nm in ref_net_names]
        )
        tran.parameters["UseCurrentNestLevel"].value = "no"

        transaction.commit()

    design.save_design()
    print(f"Created {new_cell_name} in {library_name}: Tx {tx_snp.master_lcv_name} "
          f"instance ({tx_n} ports), Rx "
          f"{rx_snp.master_lcv_name} instance ({rx_n} ports), "
          f"{vcvs_tx_gain_cell} "
          f"(delay={combined_vcvs_delay}), VectorFieldRotate_vars, "
          f"2 inductors, {cccs_rx_aperture_cell} (delay=0), VtSine1 (Tx "
          f"side, in series with Term1), Term1/Term2, S_Param controller, "
          f"Tran controller (StopTime=10 nsec, "
          f"MaxTimeStep={max_time_step_sec} sec (Nyquist/6), "
          f"MaxIters={tran_max_iters}, TruncTol/ChargeTol at ADS defaults, ImpEnforcePassivity=no, "
          f"ImpMaxFreq={combined_stop / divide_by} {freq_units}, node "
          f"voltage output restricted to Vtheta/Vphi/FeedRx, branch "
          f"currents disabled).")
    print("Usage reminder: all of eps0, mu0, eta0, Z0, r, theta, phi, a_Tx, "
          "a_Rx, rot_angle, Ux, Uy, Uz are set in this cell's own VAR1 block "
          "(global scope) -- override via var_overrides if needed. Internal "
          "ports (e.g. 'load1') each get their own Term (Z0), numbered "
          "Tx Feed=1, Tx internal ports next, then Rx Feed, then Rx internal ports. "
          "VtSine1's own sine-wave parameters (amplitude, frequency, etc.) "
          "and Tran1's StartTime are left at ADS defaults -- configure "
          "these manually for your actual desired excitation and "
          "simulation start. Tran1's "
          "ImpMaxFreq/ImpEnforcePassivity/NoiseBandwidth VALUES are set "
          "correctly by this script, but showing them as visible annotations "
          "on the schematic (like StopTime/MaxTimeStep already are by "
          "default) has no confirmed programmatic path -- toggle 'Display "
          "Parameter on Schematic' manually via the GUI for these three, "
          "same as any other parameter you want shown.")

    return design


print("Ready: create_tx_rx_link_internal_ports_TRANSIENT(tx_snp_path, "
      "rx_snp_path, ...) is defined (create_vcvs_tx_gain/"
      "create_cccs_rx_aperture/create_snp_component_cell all auto-loaded "
      "above -- Tx/Rx antenna blocks are now this antenna's own "
      "\"{antenna}_SnP\" cell, auto-built if needed, not a raw SnP "
      "instance). See the module docstring's Usage section.")
