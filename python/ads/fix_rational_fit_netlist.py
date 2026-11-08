r"""
fix_rational_fit_netlist.py -- give a rational-fit SPICE netlist real port names.

MATLAB's generateSPICE(Rfit, path, S.Impedance, "SeparateReferenceNodes",
true) names each port's signal node "po{i}" and its own reference node
"neg{i}" (used throughout that port's sub-circuit). This script renames,
throughout the file,

    po{i}  -> the i-th port name of the SnP file (e.g. "Feed", "TM_L1_M0")
    neg{i} -> that port name + ref_suffix (e.g. "Feed_ref")

so that the netlist connects by name to the cells built by
create_snp_dual_view_cells.py. Port order follows the SnP file, which
sparameters/rationalfit/generateSPICE preserve. A "* SNP_SOURCE: <path>"
line is prepended so the netlist can be matched to its SnP file later
(get_snp_source_from_netlist_comment). The MATLAB equivalent is
matlab/snp/rationalfit/fix_rational_fit_netlist.m.

Pure Python (no ADS calls); can be used in the ADS console or any Python.

Usage (ADS Python console, after start_console.py):
    run_script("fix_rational_fit_netlist.py")
    fix_rational_fit_netlist(r"<repo>\output\SnP\rationalfit\S_dom_Q_<antenna>.s8p.ckt",
                             r"<repo>\output\SnP\S_dom_Q_<antenna>.s8p")
"""

import re
from pathlib import Path

# ==================== Fixed configuration -- edit as needed ====================

PO_TOKEN_PATTERN = re.compile(r"\bpo(\d+)\b")
NEG_TOKEN_PATTERN = re.compile(r"\bneg(\d+)\b")

# ============================================================================


def get_port_names_from_snp_comments(filepath):
    """Parse '! Port[N] Name' header comments of a Touchstone file and
    return the port names in port order."""
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


def fix_rational_fit_netlist(ckt_path, snp_path, output_path=None, ref_suffix="_ref"):
    """Rename po{i}/neg{i} nodes in a generateSPICE netlist to SnP port names.

    Args:
        ckt_path: raw .ckt file from MATLAB generateSPICE with
            SeparateReferenceNodes=true.
        snp_path: matching Touchstone file (port names, in order).
        output_path: output file; default is ckt_path with "_fixed"
            inserted before the extension (the original is kept).
        ref_suffix: appended to each port name for its reference node,
            e.g. "Feed" -> "Feed_ref".

    Returns:
        Path of the written file.
    """
    ckt_path = Path(ckt_path)
    if output_path is None:
        output_path = ckt_path.with_name(ckt_path.stem + "_fixed" + ckt_path.suffix)
    else:
        output_path = Path(output_path)

    port_names = get_port_names_from_snp_comments(snp_path)
    print(f"Found {len(port_names)} ports in {snp_path}: {port_names}")

    # po{i}/neg{i} (1-indexed) -> port_names[i-1] / port_names[i-1]+ref_suffix
    po_to_port_name = {i + 1: name for i, name in enumerate(port_names)}
    neg_to_ref_name = {i + 1: f"{name}{ref_suffix}" for i, name in enumerate(port_names)}

    with open(ckt_path, "r") as f:
        content = f.read()

    po_count = 0
    neg_count = 0

    def replace_po(match):
        nonlocal po_count
        po_count += 1
        return po_to_port_name[int(match.group(1))]

    def replace_neg(match):
        nonlocal neg_count
        neg_count += 1
        return neg_to_ref_name[int(match.group(1))]

    content = PO_TOKEN_PATTERN.sub(replace_po, content)
    content = NEG_TOKEN_PATTERN.sub(replace_neg, content)

    # Record the source SnP file in a comment line, so later code (see
    # get_snp_source_from_netlist_comment) or a human can verify a given
    # rational-fit netlist actually corresponds to the SnP file it's
    # about to be paired with, rather than assuming.
    snp_path_str = str(Path(snp_path))
    content = f"* SNP_SOURCE: {snp_path_str}\r\n" + content

    with open(output_path, "w", newline="") as f:
        f.write(content)

    print(f"Renamed {po_count} po{{i}} token(s) and {neg_count} neg{{i}} token(s) "
          f"throughout the file.")
    print(f"Recorded source SnP file in a '* SNP_SOURCE:' comment line.")
    print(f"Wrote corrected netlist to {output_path}")
    return output_path


def get_snp_source_from_netlist_comment(ckt_path):
    """Return the SnP path recorded in a '* SNP_SOURCE: <path>' line by
    fix_rational_fit_netlist, or None if the netlist has no such line."""
    pattern = re.compile(r"^\*\s*SNP_SOURCE:\s*(.+)$", re.IGNORECASE)
    with open(ckt_path, "r") as f:
        for raw_line in f:
            match = pattern.match(raw_line.strip())
            if match:
                return match.group(1).strip()
    return None


print("Ready: fix_rational_fit_netlist(ckt_path, snp_path, ...) is defined. "
      "See the module docstring's Usage section.")
