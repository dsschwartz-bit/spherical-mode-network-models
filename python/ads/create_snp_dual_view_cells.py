r"""
create_snp_dual_view_cells.py -- build swappable SnP / rational-fit antenna cells in ADS.

Main functions:
  create_snp_component_cell(snp_path, ...)    cell wrapping the native SnP block
      (OneRefPinPerPort), named "<base>_TOUCHESTONE" to pair with the GUI-imported
      rational-fit cell of the same antenna (see get_gui_import_base_name).
  redraw_symbol_for_imported_cell(lib, cell)  redraw only the symbol of a rational-fit
      cell created by File > Import > Netlist File, so its pins match the SnP layout.
  create_snp_rational_fit_cell[_internal_ports](netlist_path, ...)  NetlistInclude-based
      cells. Known not to work: the resulting cell has no model_def (and the API has no
      setter), so it simulates with S11 = 1. Use the GUI import + redraw instead.
  close_cell_views(lib, cell)                 close a cell's views opened in this session.

All symbols place pins at the native SnP component's own coordinates, so the cells and a
plain SnP instance are interchangeable without re-wiring. NetlistInclude connects by net
name only, so netlist node names must equal the SnP port names.

Usage (ADS Python console, after start_console.py):
    run_script("create_snp_dual_view_cells.py")
    create_snp_component_cell(r"<repo>\output\SnP\S_dom_psi_<antenna>.s6p",
                              library_name="spice_import_lib")
    # In the GUI: File > Import > Netlist File -> the rational-fit .ckt -> Options:
    # Dialect "Spice 3", Output Format "ADS Netlist" -> choose the subcircuit with the
    # right pin count. Do not rename the imported cell (e.g. "s_dom_q_<antenna>_s12p").
    redraw_symbol_for_imported_cell("spice_import_lib", "s_dom_q_<antenna>_s12p")
"""

import re
import subprocess
from pathlib import Path
import keysight.ads.de as de
from keysight.ads.de import db_uu

# ============================ Fixed configuration ============================

ANTENNA_PREFIX_PATTERN = re.compile(r"^S_dom_(?:psi|Q)_")
MODE_PATTERN = re.compile(r"^(TM|TE)_L(\d+)_M(neg)?(\d+)$")
PIN_LAYER_ID_NUM = 231
PIN_SPACING_Y = 2.0

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
    """Return the filename stem with any S_dom_psi_/S_dom_Q_ prefix removed."""
    stem = Path(snp_path).stem
    return ANTENNA_PREFIX_PATTERN.sub("", stem)


def get_gui_import_base_name(filepath):
    """Return the cell name ADS's GUI netlist import would give this file.

    Drops a trailing ".ckt", lowercases, and replaces "." with "_", so the
    netlist and its Touchstone file give the same name (e.g.
    "S_dom_Q_x.s12p.ckt" and "S_dom_Q_x.s12p" -> "s_dom_q_x_s12p"). The
    imported cell must keep this name: the netlist content refers to it.
    """
    name = Path(filepath).name
    if name.lower().endswith(".ckt"):
        name = name[: -len(".ckt")]
    return name.lower().replace(".", "_")


def pt(point):
    """Coerce a PointF-like object (or plain tuple) to a plain (x,y) tuple."""
    return (point.x, point.y) if hasattr(point, "x") else tuple(point)


def ground_pin(design, pin, name):
    """Place a GROUND symbol on a pin (coincident pins connect without a wire)."""
    p = pt(pin.snap_point)
    return design.add_instance(("ads_rflib", "GROUND", "symbol"), p, name=name, angle=-90)


def draw_symbol_directly(pin_specs, library_name, cell_name, stub_length=0.25, title_text=None):
    """Draw a symbol view (rectangle body, pin stubs, pin labels) at exact coordinates.

    Args:
        pin_specs: list of (net_name, position, angle) -- the same list passed to
            add_external_pins_at_positions, so symbol and schematic agree. angle is
            0.0 (right), 180.0 (left) or 90.0 (bottom, e.g. a shared "ref" pin).
        stub_length: body inset from the pins on each populated side.
        title_text: optional label drawn above the body (e.g. "SnP", "Rational Fit").
    Pin names are written inside the body, like the native SnP symbol.
    """
    symbol = db_uu.create_symbol(f"{library_name}:{cell_name}:symbol")
    layer_id = db_uu.LayerId(PIN_LAYER_ID_NUM)
    font_name = "Roboto Condensed Black"
    font_height = 0.08
    text_margin = 0.05
    away_from_center_offset = 0.1

    left_xs = [pos[0] for _, pos, angle in pin_specs if angle == 180.0]
    right_xs = [pos[0] for _, pos, angle in pin_specs if angle == 0.0]
    bottom_ys = [pos[1] for _, pos, angle in pin_specs if angle == 90.0]
    all_x = [pos[0] for _, pos, _ in pin_specs]
    all_y = [pos[1] for _, pos, _ in pin_specs]

    body_left = (min(left_xs) + stub_length) if left_xs else min(all_x)
    body_right = (max(right_xs) - stub_length) if right_xs else max(all_x)
    body_bottom = (min(bottom_ys) + stub_length) if bottom_ys else min(all_y)
    body_top = max(all_y)

    with de.db_uu.Transaction(symbol) as transaction:
        symbol.add_rectangle(layer_id, (body_left, body_bottom), (body_right, body_top))

        for term_number, (net_name, pos, angle) in enumerate(pin_specs, start=1):
            if angle == 180.0:
                body_edge = (body_left, pos[1])
                text_origin = (body_left + text_margin - away_from_center_offset, pos[1])
            elif angle == 0.0:
                body_edge = (body_right, pos[1])
                text_width = symbol.calculate_text_width(net_name, font_height, font_name)
                text_origin = (body_right - text_margin - text_width + away_from_center_offset, pos[1])
            elif angle == 90.0:
                body_edge = (pos[0], body_bottom)
                text_width = symbol.calculate_text_width(net_name, font_height, font_name)
                text_origin = (pos[0] - text_width / 2, body_bottom + text_margin)
            else:
                raise ValueError(
                    f"draw_symbol_directly: unsupported angle {angle} for pin "
                    f"'{net_name}' -- only 0.0, 180.0, or 90.0 are handled."
                )
            symbol.add_line(layer_id, [pos, body_edge])
            net = symbol.find_or_add_net(net_name)
            # add_term creates name-bound terms; add_numbered_term gives the
            # number-bound terms the native SnP uses. A design cannot mix the two.
            term = symbol.add_numbered_term(net, net_name, term_number, db_uu.TermType.INPUT_OUTPUT)
            dot = symbol.add_dot(layer_id, pos)
            symbol.add_pin(term, dot, angle=angle)
            symbol.add_text(layer_id, net_name, text_origin, font_name, font_height)

        if title_text is not None:
            title_font_height = 0.16
            title_width = symbol.calculate_text_width(title_text, title_font_height, font_name)
            title_center_x = (body_left + body_right) / 2
            title_origin = (title_center_x - title_width / 2, body_top + text_margin)
            symbol.add_text(layer_id, title_text, title_origin, font_name, title_font_height)
        transaction.commit()

    symbol.save_design()
    title_note = f" Titled '{title_text}' above the body." if title_text else ""
    print(f"Drew symbol directly for {library_name}:{cell_name} "
          f"({len(pin_specs)} pins), matching the schematic's own exact "
          f"native pin coordinates -- no SymbolGenerator involved. Each "
          f"pin's name is labeled inside the body, matching the native "
          f"SnP block's own convention.{title_note}")
    return symbol


def create_snp_component_cell(
    snp_path,
    library_name="SphericalHarmonics_lib",
    new_cell_name=None,
):
    """Create a cell wrapping a native SnP block.

    Uses OneRefPinPerPort (N signal + N reference pins, all external and
    unconnected; ground the Feed/mode references after instantiating), the
    same scheme as the rational-fit netlists from fix_rational_fit_netlist,
    so both cells are pin-for-pin swappable. An existing cell of the same
    name is deleted first.

    Args:
        snp_path: Touchstone (.sNp) file.
        library_name: target library.
        new_cell_name: default "<get_gui_import_base_name(snp_path)>_TOUCHESTONE".
    Returns:
        The schematic Design (both views closed).
    """
    if new_cell_name is None:
        base_name = get_gui_import_base_name(snp_path)
        new_cell_name = f"{base_name}_TOUCHESTONE"
        print(f"Base name '{base_name}' extracted -- using cell name "
              f"'{new_cell_name}'.")

    port_names = get_port_names_from_snp_comments(snp_path)
    port_count = len(port_names)
    print(f"Found {port_count} ports: {port_names}")

    ws = de.active_workspace()
    print(f"Active workspace: {ws}")
    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    library = ws.open_library(library_name)

    # Creating a cell that already exists fails, so delete it first.
    if library.cell_exists(new_cell_name):
        print(f"Cell '{new_cell_name}' already exists in '{library_name}' "
              f"-- deleting it first.")
        library.delete_cell(new_cell_name)

    design = de.db_uu.create_schematic((library_name, new_cell_name, "schematic"))

    # Place the (kept) SnP instance and read its pin positions, so the external
    # pins below coincide with its pins and match the native SnP layout.
    snp, positions = get_oneref_snp_positions(design, snp_path, port_names, (0, 0), name="SnP1")

    pins_by_term = {p.inst_term.term_number: p for p in snp.inst_pins}
    ref_suffix = "_ref"
    with de.db_uu.Transaction(design) as transaction:
        for i, name in enumerate(port_names, start=1):
            sig_pin = pins_by_term[2 * i - 1]
            ref_pin = pins_by_term[2 * i]
            sig_pin.add_label(name, pt(sig_pin.snap_point))
            ref_pin.add_label(f"{name}{ref_suffix}", pt(ref_pin.snap_point))
        transaction.commit()

    pin_specs = []
    for name in port_names:
        sig_pos, sig_angle, ref_pos, ref_angle = positions[name]
        pin_specs.append((name, sig_pos, sig_angle))
        pin_specs.append((f"{name}{ref_suffix}", ref_pos, ref_angle))
    add_external_pins_at_positions(design, pin_specs)

    design.save_design()
    print(f"Created {library_name}:{new_cell_name} (SnP-based, "
          f"OneRefPinPerPort), {port_count} signal pins + {port_count} "
          f"reference pins at the SnP's own native positions, all "
          f"external and UNCONNECTED -- ground Feed/mode reference "
          f"pins after instantiating.")

    symbol = draw_symbol_directly(pin_specs, library_name, new_cell_name, title_text="SnP")

    design.close_design()
    symbol.close_design()
    return design


def place_and_configure_singref_snp(design, snp_path, port_names, origin, name="SnP1"):
    """Place an SnP instance (SingRefPin, DualInLine) at origin.

    Returns (instance, pins_by_term).
    """
    port_count = len(port_names)
    longest = max(port_names, key=len)
    port_names_string = f"1 {len(longest)} {' '.join(port_names)}"

    with de.db_uu.Transaction(design) as transaction:
        snp = design.add_instance(
            db_uu.LCVName("ads_datacmps", "SnP", "symbol"), origin, name=name
        )
        transaction.commit()

    with de.db.Transaction(design, "Edit params") as transaction:
        snp.parameters["File"].value = f'"{snp_path}"'
        snp.parameters["NumPorts"].value = str(port_count)
        snp.parameters["Type"].value = '"touchstone"'
        snp.parameters["PinSpacing"].value = "Loose"
        snp.parameters["HasReferencePin"].value = "SingRefPin"
        snp.parameters["PinLayout"].value = "DualInLine"
        snp.parameters["InterpMode"].value = '"linear"'
        snp.parameters["EnforcePassivity"].value = "no"
        snp.parameters["port_name_list"].value = port_names_string
        snp.invoke_item_parameter_changed_callback(
            ["File", "NumPorts", "Type", "PinSpacing", "HasReferencePin",
             "PinLayout", "InterpMode", "EnforcePassivity", "port_name_list"]
        )
        transaction.commit()

    pins_by_term = {p.inst_term.term_number: p for p in snp.inst_pins}
    return snp, pins_by_term


def get_singref_snp_positions(design, snp_path, port_names, origin, name="SnP1"):
    """Place a SingRefPin SnP with its bottom-left pin at origin and read its pins.

    A probe instance is placed first to measure the bottom-left pin offset,
    deleted, and the returned instance placed at the corrected position.
    Signal pin angle is 180 (left column) or 0 (right column), from the
    midpoint of the observed x range; the reference pin is 90 (bottom).

    Returns (instance, {name: (position, angle)}, (ref_position, ref_angle)).
    """
    origin_x, origin_y = origin
    port_count = len(port_names)

    # ---- Probe pass: measure the bottom-left offset, then discard ----
    probe, probe_pins_by_term = place_and_configure_singref_snp(
        design, snp_path, port_names, origin, name=f"{name}_probe"
    )
    probe_x = []
    probe_y = []
    for i in range(1, port_count + 1):
        pos = pt(probe_pins_by_term[i].snap_point)
        probe_x.append(pos[0])
        probe_y.append(pos[1])
    probe_ref_pos = pt(probe_pins_by_term[port_count + 1].snap_point)
    probe_x.append(probe_ref_pos[0])
    probe_y.append(probe_ref_pos[1])
    bl_offset_x = min(probe_x) - origin_x
    bl_offset_y = min(probe_y) - origin_y
    probe.delete_object()

    # ---- Final pass: place at the corrected position, so the bottom-left pin lands exactly at origin ----
    corrected_origin = (origin_x - bl_offset_x, origin_y - bl_offset_y)
    snp, pins_by_term = place_and_configure_singref_snp(design, snp_path, port_names, corrected_origin, name)

    raw_positions = {}
    all_x = []
    for i, name_i in enumerate(port_names, start=1):
        pos = pt(pins_by_term[i].snap_point)
        raw_positions[name_i] = pos
        all_x.append(pos[0])

    ref_pos_raw = pt(pins_by_term[port_count + 1].snap_point)
    all_x.append(ref_pos_raw[0])

    mid_x = (min(all_x) + max(all_x)) / 2

    positions = {}
    for name_i in port_names:
        pos = raw_positions[name_i]
        angle = 180.0 if pos[0] < mid_x else 0.0
        positions[name_i] = (pos, angle)

    ref_angle = 90.0
    ref_pos = ref_pos_raw

    return snp, positions, (ref_pos, ref_angle)


def place_and_configure_oneref_snp(design, snp_path, port_names, origin, name="SnP1"):
    """Place an SnP instance (OneRefPinPerPort, DualInLine) at origin.

    The full parameter set is applied so the pin geometry matches a real
    instance exactly. Returns (instance, pins_by_term).
    """
    port_count = len(port_names)
    longest = max(port_names, key=len)
    port_names_string = f"1 {len(longest)} {' '.join(port_names)}"

    with de.db_uu.Transaction(design) as transaction:
        snp = design.add_instance(
            db_uu.LCVName("ads_datacmps", "SnP", "symbol"), origin, name=name
        )
        transaction.commit()

    with de.db.Transaction(design, "Edit params") as transaction:
        snp.parameters["File"].value = f'"{snp_path}"'
        snp.parameters["NumPorts"].value = str(port_count)
        snp.parameters["Type"].value = '"touchstone"'
        snp.parameters["PinSpacing"].value = "Loose"
        snp.parameters["HasReferencePin"].value = "OneRefPinPerPort"
        snp.parameters["PinLayout"].value = "DualInLine"
        snp.parameters["InterpMode"].value = '"linear"'
        snp.parameters["EnforcePassivity"].value = "no"
        snp.parameters["port_name_list"].value = port_names_string
        snp.invoke_item_parameter_changed_callback(
            ["File", "NumPorts", "Type", "PinSpacing", "HasReferencePin",
             "PinLayout", "InterpMode", "EnforcePassivity", "port_name_list"]
        )
        transaction.commit()

    pins_by_term = {p.inst_term.term_number: p for p in snp.inst_pins}
    return snp, pins_by_term


def get_oneref_snp_positions(design, snp_path, port_names, origin, name="SnP1"):
    """Place a OneRefPinPerPort SnP with its bottom-left pin at origin and read its pins.

    Same probe-then-place approach as get_singref_snp_positions. Angles
    (180 left / 0 right) are taken from the actual pin x positions, since
    the DualInLine column assignment depends on the port count. The caller
    may delete the returned instance (instance.delete_object()).

    Returns (instance, {name: (sig_pos, sig_angle, ref_pos, ref_angle)}).
    """
    origin_x, origin_y = origin

    # ---- Probe pass: measure the bottom-left offset, then discard ----
    probe, probe_pins_by_term = place_and_configure_oneref_snp(
        design, snp_path, port_names, origin, name=f"{name}_probe"
    )
    probe_x = []
    probe_y = []
    for i in range(1, len(port_names) + 1):
        sig_pos = pt(probe_pins_by_term[2 * i - 1].snap_point)
        ref_pos = pt(probe_pins_by_term[2 * i].snap_point)
        probe_x += [sig_pos[0], ref_pos[0]]
        probe_y += [sig_pos[1], ref_pos[1]]
    bl_offset_x = min(probe_x) - origin_x
    bl_offset_y = min(probe_y) - origin_y
    probe.delete_object()

    # ---- Final pass: place at the corrected position, so the bottom-left pin lands exactly at origin ----
    corrected_origin = (origin_x - bl_offset_x, origin_y - bl_offset_y)
    snp, pins_by_term = place_and_configure_oneref_snp(design, snp_path, port_names, corrected_origin, name)

    raw_positions = {}
    all_x = []
    for i, name_i in enumerate(port_names, start=1):
        sig_pos = pt(pins_by_term[2 * i - 1].snap_point)
        ref_pos = pt(pins_by_term[2 * i].snap_point)
        raw_positions[name_i] = (sig_pos, ref_pos)
        all_x += [sig_pos[0], ref_pos[0]]

    mid_x = (min(all_x) + max(all_x)) / 2

    positions = {}
    for name_i in port_names:
        sig_pos, ref_pos = raw_positions[name_i]
        sig_angle = 180.0 if sig_pos[0] < mid_x else 0.0
        ref_angle = 180.0 if ref_pos[0] < mid_x else 0.0
        positions[name_i] = (sig_pos, sig_angle, ref_pos, ref_angle)
    return snp, positions


def add_external_pins_at_positions(design, pin_specs):
    """Add one external pin per (net_name, position, angle) in pin_specs.

    Uses number-bound terms (add_numbered_term, as the native SnP does; a
    design cannot mix name- and number-bound terms), numbered 1.. in list
    order -- the same numbering draw_symbol_directly uses.
    """
    pin_layer_id = db_uu.LayerId(PIN_LAYER_ID_NUM)
    with de.db_uu.Transaction(design) as transaction:
        for term_number, (net_name, position, angle) in enumerate(pin_specs, start=1):
            net = design.find_or_add_net(net_name)
            term = design.add_numbered_term(net, net_name, term_number, db_uu.TermType.INPUT_OUTPUT)
            dot = design.add_dot(pin_layer_id, position)
            design.add_pin(term, dot, angle=angle)
        transaction.commit()



def get_snp_source_from_netlist_comment(netlist_path):
    """Return the path in the 'SNP_SOURCE: <path>' comment written by
    fix_rational_fit_netlist, or None. Accepts '*' (SPICE) and ';' (ADS
    Netlist, as produced by nettrans) comment characters."""
    pattern = re.compile(r"^[*;]\s*SNP_SOURCE:\s*(.+)$", re.IGNORECASE)
    with open(netlist_path, "r") as f:
        for raw_line in f:
            match = pattern.match(raw_line.strip())
            if match:
                return match.group(1).strip()
    return None


def translate_netlist_to_ads(input_path, output_path=None, dialect="-s3"):
    """Translate a SPICE netlist to ADS Netlist format with nettrans.

    NetlistInclude does not read raw SPICE correctly. Flags: "-g" (ADS
    Netlist output) and "-n" (keep node-name case, needed for net-name
    matching). stdin=DEVNULL avoids "OSError: [WinError 6]" inside the ADS
    console.

    Args:
        input_path: SPICE netlist (.net/.ckt).
        output_path: default "<stem with '.'->'_'>_ads.net" next to the input
            (e.g. "foo.s8p.ckt" -> "foo_s8p_ads.net").
        dialect: nettrans input flag ("-s2", "-s3", "-p", "-lt"); "-s3"
            matches MATLAB generateSPICE output.
    Returns:
        Path of the translated file.
    """
    input_path = Path(input_path)
    if output_path is None:
        sanitized_stem = input_path.stem.replace(".", "_")
        output_path = input_path.parent / f"{sanitized_stem}_ads.net"
    else:
        output_path = Path(output_path)

    result = subprocess.run(
        ["nettrans", str(input_path), str(output_path), dialect, "-g", "-n"],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
    )
    print(f"nettrans stdout: {result.stdout}")
    if result.stderr:
        print(f"nettrans stderr: {result.stderr}")
    if result.returncode != 0:
        raise RuntimeError(
            f"nettrans failed (return code {result.returncode}) translating "
            f"{input_path} -- see output above. Check nettrans.log (likely "
            f"in the current working directory) for details."
        )
    if not output_path.exists():
        raise RuntimeError(
            f"nettrans reported success but {output_path} was not created -- "
            f"check nettrans.log for details."
        )
    print(f"Translated {input_path} -> {output_path}")
    return output_path


def parse_ads_netlist_define_blocks(netlist_path):
    """Return [(name, node_list), ...] for every "define <name> ( nodes )"
    block in an ADS Netlist file, in file order.

    Handles backslash line continuations. Translated files can contain more
    than one block (including a 0-node one), so callers pick by node count.
    """
    with open(netlist_path, "r") as f:
        content = f.read()
    blocks = []
    for match in re.finditer(r"define\s+(\S+)\s*\((.*?)\)", content, re.DOTALL):
        name = match.group(1)
        node_block = match.group(2).replace("\\", " ")
        nodes = node_block.split()
        blocks.append((name, nodes))
    if not blocks:
        raise ValueError(f"No 'define(...)' block found in {netlist_path}")
    return blocks


def parse_netlist_subckt_nodes(netlist_path):
    """Return the node list of the first ".SUBCKT" line in a raw SPICE netlist.

    Diagnostic helper; translated ADS netlists use "define" blocks instead
    (see parse_ads_netlist_define_blocks).
    """
    with open(netlist_path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if line[:7].lower() == ".subckt":
                parts = line.split()
                return parts[2:]
    raise ValueError(f"No '.SUBCKT' line found in {netlist_path}")


def close_cell_views(library_name, cell_name):
    """Close a cell's schematic and symbol views opened in this session, so
    GUI actions such as delete do not report the cell as in use."""
    ws = de.active_workspace()
    for view_name in ("schematic", "symbol"):
        try:
            design = de.db_uu.open_design((library_name, cell_name, view_name))
            design.close_design()
        except Exception:
            pass


def redraw_symbol_for_imported_cell(library_name, cell_name, snp_path=None, title_text="Rational Fit"):
    r"""Replace the symbol of a GUI-imported rational-fit cell with one matching
    the native SnP pin layout. The schematic (and its model_def) is untouched.

    Pin positions come from a temporary "_pin_layout_scratch" cell, deleted
    before returning. The reference scheme (N, N+1 or 2N terms) is inferred
    from the existing term count.

    Args:
        library_name, cell_name: the imported cell, with the name the GUI gave
            it (see get_gui_import_base_name); do not rename it.
        snp_path: Touchstone file for port names; default read from the
            SNP_SOURCE comment of the cell's NetlistInclude file.
        title_text: label drawn above the symbol body.
    Returns:
        The redrawn symbol Design (both views closed).

    Example:
        redraw_symbol_for_imported_cell("spice_import_lib", "s_dom_q_<antenna>_s12p")
    """
    ws = de.active_workspace()
    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    library = ws.open_library(library_name)
    if not library.cell_exists(cell_name):
        raise RuntimeError(
            f"Cell '{cell_name}' not found in '{library_name}' -- this "
            f"function expects an ALREADY-EXISTING cell, created via "
            f"ADS's own GUI import (File > Import > Netlist File), not "
            f"one built by this file's own create_snp_rational_fit_cell/"
            f"internal_ports."
        )

    existing_design = de.db_uu.open_design((library_name, cell_name, "schematic"))
    existing_names = [t.name for t in existing_design.terms]
    print(f"Found {len(existing_names)} existing external term(s) in "
          f"{library_name}:{cell_name}: {existing_names}")

    if snp_path is None:
        netlist_inc = None
        for inst in existing_design.instances:
            if "NetlistInclude" in str(inst.master_lcv_name):
                netlist_inc = inst
                break
        if netlist_inc is None:
            raise RuntimeError(
                f"No NetlistInclude instance found in {library_name}:"
                f"{cell_name} -- cannot auto-detect snp_path. Pass it "
                f"explicitly instead."
            )
        include_dirs = netlist_inc.parameters["IncludePath"].value
        include_files = netlist_inc.parameters["IncludeFiles"].value
        netlist_full_path = Path(include_dirs[0]) / include_files[0][0]
        recorded_source = get_snp_source_from_netlist_comment(netlist_full_path)
        if recorded_source is None:
            raise ValueError(
                f"No snp_path provided, and no 'SNP_SOURCE:' comment "
                f"found in {netlist_full_path} to read it from "
                f"automatically. Pass snp_path explicitly instead."
            )
        snp_path = recorded_source
        print(f"snp_path not provided -- auto-detected from "
              f"{netlist_inc.name}'s own included file's 'SNP_SOURCE:' "
              f"comment: {snp_path}")

    port_names = get_port_names_from_snp_comments(snp_path)
    port_count = len(port_names)
    print(f"Found {port_count} ports (from the SnP file): {port_names}")

    if len(existing_names) == port_count:
        scheme = "N"
    elif len(existing_names) == port_count + 1:
        scheme = "N+1"
    elif len(existing_names) == 2 * port_count:
        scheme = "2N"
    else:
        raise ValueError(
            f"Existing cell has {len(existing_names)} external term(s), "
            f"but expected {port_count} (implicit ground), "
            f"{port_count + 1} (single shared reference), or "
            f"{2 * port_count} (per-port reference) for {port_count} "
            f"ports -- cannot determine which reference scheme this "
            f"cell actually uses."
        )
    print(f"Existing term count ({len(existing_names)}) matches the "
          f"'{scheme}' reference scheme.")

    scratch_cell_name = "_pin_layout_scratch"
    if library.cell_exists(scratch_cell_name):
        # Leftover from an interrupted run: remove it so create_schematic succeeds.
        stale_design = de.db_uu.open_design((library_name, scratch_cell_name, "schematic"))
        stale_design.close_design()
        library.delete_cell(scratch_cell_name)
    scratch_design = de.db_uu.create_schematic((library_name, scratch_cell_name, "schematic"))
    try:
        if scheme == "2N":
            scratch_snp, positions = get_oneref_snp_positions(scratch_design, snp_path, port_names, (0, 0), name="scratch_SnP")
            scratch_snp.delete_object()
            pin_specs = []
            for name in existing_names:
                if name in positions:
                    pos, angle = positions[name][0], positions[name][1]
                elif name.endswith("_ref") and name[:-len("_ref")] in positions:
                    base = name[:-len("_ref")]
                    pos, angle = positions[base][2], positions[base][3]
                else:
                    raise ValueError(
                        f"Existing term '{name}' does not match any port "
                        f"name or '{{port}}_ref' pattern from {port_names} "
                        f"-- cannot map it to a native pin position."
                    )
                pin_specs.append((name, pos, angle))
        else:
            scratch_snp, positions, (ref_pos, ref_angle) = get_singref_snp_positions(scratch_design, snp_path, port_names, (0, 0), name="scratch_SnP")
            scratch_snp.delete_object()
            pin_specs = [(name, pos, angle) for name, (pos, angle) in positions.items() if name in existing_names]
            if scheme == "N+1":
                extra_names = [n for n in existing_names if n not in positions]
                if len(extra_names) != 1:
                    raise ValueError(
                        f"Expected exactly 1 existing term not matching a "
                        f"port name (the shared reference), found "
                        f"{len(extra_names)}: {extra_names}"
                    )
                pin_specs.append((extra_names[0], ref_pos, ref_angle))
    finally:
        # Library.delete_cell fails with "cell has open cellviews" unless the
        # design is closed first.
        scratch_design.close_design()
        library.delete_cell(scratch_cell_name)

    existing_design.cell.delete_view("symbol")
    symbol = draw_symbol_directly(pin_specs, library_name, cell_name, title_text=title_text)
    existing_design.close_design()
    symbol.close_design()
    print(f"Redrew ONLY the symbol view of {library_name}:{cell_name} -- "
          f"its schematic (and model_def) are untouched. Both views "
          f"closed -- ready for a GUI rename/delete if desired.")
    return symbol


def create_snp_rational_fit_cell(
    netlist_path,
    snp_path=None,
    library_name="SphericalHarmonics_lib",
    new_cell_name=None,
    ref_suffix="_ref",
):
    r"""Create a cell wrapping a rational-fit netlist via NetlistInclude.

    Known not to simulate correctly (S11 = 1; see module docstring). Use the
    GUI import + redraw_symbol_for_imported_cell instead.

    The netlist is translated with nettrans, and its node count selects the
    reference scheme: N (implicit ground, no "ref" pin), N+1 (shared "ref"
    pin, tied by a 0-Ohm resistor if the node has another name), or 2N
    (delegates to create_snp_rational_fit_cell_internal_ports).

    Args:
        netlist_path: SPICE rational-fit netlist whose node names equal the
            SnP port names.
        snp_path: Touchstone file for the port names; default from the
            netlist's SNP_SOURCE comment (a mismatch only prints a warning).
        library_name: target library.
        new_cell_name: default "<antenna>_RationalFit".
        ref_suffix: reference-node suffix (2N case only).
    Returns:
        The schematic Design (both views closed).
    """
    recorded_source = get_snp_source_from_netlist_comment(netlist_path)
    if snp_path is None:
        if recorded_source is None:
            raise ValueError(
                f"No snp_path provided, and no 'SNP_SOURCE:' comment "
                f"found in {netlist_path} to read it from automatically. "
                f"Either pass snp_path explicitly, or regenerate the "
                f"netlist via fix_rational_fit_netlist.py/.m, which "
                f"records this comment automatically."
            )
        snp_path = recorded_source
        print(f"snp_path not provided -- read from netlist's own "
              f"'SNP_SOURCE:' comment: {snp_path}")
    elif recorded_source is not None and Path(recorded_source) != Path(str(snp_path)):
        print(f"WARNING: netlist's recorded source ('{recorded_source}') "
              f"does not match the snp_path explicitly passed in "
              f"('{snp_path}') -- proceeding with the explicitly-provided "
              f"snp_path, but double-check these are meant to be paired.")

    port_names = get_port_names_from_snp_comments(snp_path)
    port_count = len(port_names)
    print(f"Found {port_count} ports (from the SnP file, for pin "
          f"matching only): {port_names}")

    # NetlistInclude needs ADS Netlist format, not raw SPICE.
    translated_path = translate_netlist_to_ads(netlist_path)

    # Classify by node count before setting new_cell_name, so the 2N delegate
    # still applies its own default name. The translated file may hold several
    # define blocks; use the first whose node count fits a known scheme.
    blocks = parse_ads_netlist_define_blocks(translated_path)
    netlist_nodes = None
    for block_name, block_nodes in blocks:
        if len(block_nodes) in (port_count, port_count + 1, 2 * port_count):
            netlist_nodes = block_nodes
            print(f"Using define block '{block_name}' ({len(block_nodes)} "
                  f"nodes) from {translated_path}.")
            break
    if netlist_nodes is None:
        block_summary = ", ".join(f"'{n}' ({len(ns)} nodes)" for n, ns in blocks)
        raise ValueError(
            f"None of the define(...) block(s) found in {translated_path} "
            f"have a node count matching {port_count} (implicit ground), "
            f"{port_count + 1} (single shared reference), or "
            f"{2 * port_count} (per-port reference) for {port_count} ports. "
            f"Found: {block_summary}"
        )

    if len(netlist_nodes) == 2 * port_count:
        print(f"Netlist has {len(netlist_nodes)} nodes (2x the {port_count} "
              f"ports) -- a per-port reference scheme. Automatically "
              f"delegating to create_snp_rational_fit_cell_internal_ports.")
        return create_snp_rational_fit_cell_internal_ports(
            netlist_path,
            snp_path=snp_path,
            library_name=library_name,
            new_cell_name=new_cell_name,
            ref_suffix=ref_suffix,
        )
    elif len(netlist_nodes) == port_count:
        has_shared_ref = False
        extra_ref_name = None
    elif len(netlist_nodes) == port_count + 1:
        extra_nodes = [n for n in netlist_nodes if n not in port_names]
        if len(extra_nodes) != 1:
            raise ValueError(
                f"Netlist has {port_count + 1} nodes but its signal node "
                f"names don't cleanly match the SnP file's {port_count} "
                f"port names (found {len(extra_nodes)} unmatched node(s): "
                f"{extra_nodes}) -- expected exactly 1 (the shared "
                f"reference node)."
            )
        has_shared_ref = True
        extra_ref_name = extra_nodes[0]
        print(f"Netlist uses a single shared reference node: '{extra_ref_name}'.")
    else:
        raise ValueError(
            f"Netlist has {len(netlist_nodes)} nodes, but expected {port_count} "
            f"(implicit ground) or {port_count + 1} (single shared reference) "
            f"or {2 * port_count} (per-port reference) for {port_count} ports."
        )

    # N / N+1 cases only (2N returned above).
    if new_cell_name is None:
        antenna_name = extract_antenna_name(snp_path)
        new_cell_name = f"{antenna_name}_RationalFit"
        print(f"Antenna name '{antenna_name}' extracted -- using cell name "
              f"'{new_cell_name}'.")

    translated_path = Path(translated_path)
    directory = str(translated_path.parent)
    filename = translated_path.name

    ws = de.active_workspace()
    print(f"Active workspace: {ws}")
    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    library = ws.open_library(library_name)

    # Creating a cell that already exists fails, so delete it first.
    if library.cell_exists(new_cell_name):
        print(f"Cell '{new_cell_name}' already exists in '{library_name}' "
              f"-- deleting it first.")
        library.delete_cell(new_cell_name)

    design = de.db_uu.create_schematic((library_name, new_cell_name, "schematic"))

    with de.db_uu.Transaction(design) as transaction:
        netlist_inc = design.add_instance(
            ("ads_simulation", "NetlistInclude", "symbol"), (5, 0), name="NetlistInclude1"
        )
        netlist_inc.parameters["IncludePath"].value = [directory]
        netlist_inc.parameters["IncludeFiles"].value = [[filename, ""]]
        transaction.commit()

    # Scratch SnP instance, used only to read native pin positions, then deleted.
    scratch_snp, positions, (ref_pos, ref_angle) = get_singref_snp_positions(design, snp_path, port_names, (0, 0), name="scratch_SnP")
    scratch_snp.delete_object()

    pin_specs = [(name, pos, angle) for name, (pos, angle) in positions.items()]
    if has_shared_ref:
        pin_specs.append(("ref", ref_pos, ref_angle))
    add_external_pins_at_positions(design, pin_specs)

    if has_shared_ref and extra_ref_name.lower() != "ref":
        # Tie the netlist's reference node to the "ref" pin with a 0-Ohm resistor.
        with de.db_uu.Transaction(design) as transaction:
            tie_pos = (ref_pos[0], ref_pos[1] - 1.0)
            tie_r = design.add_instance(("ads_rflib", "R", "symbol"), tie_pos, name="R_ref_tie")
            tie_r.parameters["R"].value = "0 Ohm"
            tie_pins = {p.inst_term.term_number: p for p in tie_r.inst_pins}
            tie_pins[1].add_label(extra_ref_name, pt(tie_pins[1].snap_point))
            tie_pins[2].add_label("ref", pt(tie_pins[2].snap_point))
            transaction.commit()

    design.save_design()
    if has_shared_ref:
        print(f"Created {library_name}:{new_cell_name} (rational-fit based), "
              f"NetlistInclude1 pointing at the translated netlist "
              f"{translated_path} (translated from {netlist_path} via "
              f"nettrans), "
              f"{port_count} signal pins + 1 shared reference pin ('ref') "
              f"at the SnP's own native positions, matching "
              f"create_snp_component_cell's symbol exactly. 'ref' is left "
              f"UNCONNECTED (not grounded) -- ground/wire it after "
              f"instantiating, same as create_snp_component_cell's own "
              f"'ref' pin.")
    else:
        print(f"Created {library_name}:{new_cell_name} (rational-fit based), "
              f"NetlistInclude1 pointing at the translated netlist "
              f"{translated_path} (translated from {netlist_path} via "
              f"nettrans), "
              f"{port_count} signal pins at the SnP's own native positions. "
              f"No 'ref' pin -- the netlist uses SPICE's implicit global "
              f"ground ('0') internally.")

    symbol = draw_symbol_directly(pin_specs, library_name, new_cell_name, title_text="Rational Fit")

    design.close_design()
    symbol.close_design()

    return design


def create_snp_rational_fit_cell_internal_ports(
    netlist_path,
    snp_path=None,
    library_name="SphericalHarmonics_lib",
    new_cell_name=None,
    ref_suffix="_ref",
):
    r"""Create a NetlistInclude cell with one external reference pin per port.

    Known not to simulate correctly (S11 = 1; see module docstring). Use the
    GUI import + redraw_symbol_for_imported_cell instead.

    Expects a netlist from MATLAB generateSPICE(..., "SeparateReferenceNodes",
    true) processed by fix_rational_fit_netlist, i.e. node pairs
    "<name>"/"<name><ref_suffix>". All pins are left unconnected; ground the
    Feed/mode references after instantiating.

    Args:
        netlist_path: fixed SPICE rational-fit netlist.
        snp_path: Touchstone file for the port names; default from the
            netlist's SNP_SOURCE comment.
        library_name: target library.
        new_cell_name: default "<antenna>_RationalFit".
        ref_suffix: reference-node suffix used by fix_rational_fit_netlist.
    Returns:
        The schematic Design (both views closed).
    """
    recorded_source = get_snp_source_from_netlist_comment(netlist_path)
    if snp_path is None:
        if recorded_source is None:
            raise ValueError(
                f"No snp_path provided, and no 'SNP_SOURCE:' comment "
                f"found in {netlist_path} to read it from automatically. "
                f"Either pass snp_path explicitly, or regenerate the "
                f"netlist via fix_rational_fit_netlist.py/.m, which "
                f"records this comment automatically."
            )
        snp_path = recorded_source
        print(f"snp_path not provided -- read from netlist's own "
              f"'SNP_SOURCE:' comment: {snp_path}")
    elif recorded_source is not None and Path(recorded_source) != Path(str(snp_path)):
        print(f"WARNING: netlist's recorded source ('{recorded_source}') "
              f"does not match the snp_path explicitly passed in "
              f"('{snp_path}') -- proceeding with the explicitly-provided "
              f"snp_path, but double-check these are meant to be paired.")

    if new_cell_name is None:
        antenna_name = extract_antenna_name(snp_path)
        new_cell_name = f"{antenna_name}_RationalFit"
        print(f"Antenna name '{antenna_name}' extracted -- using cell name "
              f"'{new_cell_name}'.")

    port_names = get_port_names_from_snp_comments(snp_path)
    print(f"Found {len(port_names)} ports (from the SnP file, for pin "
          f"matching only): {port_names}")

    # NetlistInclude needs ADS Netlist format, not raw SPICE.
    translated_path = translate_netlist_to_ads(netlist_path)
    directory = str(translated_path.parent)
    filename = translated_path.name

    ws = de.active_workspace()
    print(f"Active workspace: {ws}")
    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    library = ws.open_library(library_name)

    # Creating a cell that already exists fails, so delete it first.
    if library.cell_exists(new_cell_name):
        print(f"Cell '{new_cell_name}' already exists in '{library_name}' "
              f"-- deleting it first.")
        library.delete_cell(new_cell_name)

    design = de.db_uu.create_schematic((library_name, new_cell_name, "schematic"))

    with de.db_uu.Transaction(design) as transaction:
        netlist_inc = design.add_instance(
            ("ads_simulation", "NetlistInclude", "symbol"), (5, 0), name="NetlistInclude1"
        )
        netlist_inc.parameters["IncludePath"].value = [directory]
        netlist_inc.parameters["IncludeFiles"].value = [[filename, ""]]
        transaction.commit()

    # Scratch SnP instance, used only to read native pin positions, then deleted.
    scratch_snp, positions = get_oneref_snp_positions(design, snp_path, port_names, (0, 0), name="scratch_SnP")
    scratch_snp.delete_object()

    pin_specs = []
    for name in port_names:
        sig_pos, sig_angle, ref_pos, ref_angle = positions[name]
        pin_specs.append((name, sig_pos, sig_angle))
        pin_specs.append((f"{name}{ref_suffix}", ref_pos, ref_angle))
    add_external_pins_at_positions(design, pin_specs)

    design.save_design()
    print(f"Created {library_name}:{new_cell_name} (rational-fit based, "
          f"internal-ports-aware), NetlistInclude1 pointing at the "
          f"translated netlist {translated_path} (translated from "
          f"{netlist_path} via nettrans), {len(port_names)} signal pins + "
          f"{len(port_names)} reference pins at the native SnP's own "
          f"positions, all external and unconnected -- ground Feed/mode "
          f"reference pins after instantiating.")

    symbol = draw_symbol_directly(pin_specs, library_name, new_cell_name, title_text="Rational Fit")

    design.close_design()
    symbol.close_design()

    return design


print("Ready: create_snp_component_cell, redraw_symbol_for_imported_cell, "
      "create_snp_rational_fit_cell[_internal_ports] (not recommended: gives "
      "S11 = 1; use the GUI import + redraw_symbol_for_imported_cell instead) "
      "and translate_netlist_to_ads are defined. See the module docstring.")