r"""
inspect_design.py -- dump the contents of an existing ADS design.

Defines inspect_design(cellview_ref, to_file=False), which prints every
instance with its parameters and pin-to-net connections, then the design's
external (symbol) pins. Read-only. The external-pin attribute name is not
documented, so several candidates (terms, pins, ...) are tried and the
relevant dir(design) entries are printed if none exists.

Usage (ADS Python console, after start_console.py):
    run_script("inspect_design.py")
    inspect_design(("SphericalHarmonics_lib", "Ftheta_Fphi_vars", "schematic"))
    inspect_design(("<library>", "<cell>", "schematic"), to_file=True)
        # -> <repo>/output/ads/<library>_<cell>_schematic.txt
"""

import contextlib
from pathlib import Path
import keysight.ads.de as de
from keysight.ads.de import db_uu

# ==================== Configuration -- edit as needed ====================

# Repository folders (defined by start_console.py; derived here if this file is exec()'d directly).
try:
    SCRIPT_DIR
except NameError:
    SCRIPT_DIR = Path(THIS_SCRIPT_PATH).resolve().parent
OUTPUT_DIR = SCRIPT_DIR.parents[1] / "output" / "ads"   # <repo>/output/ads

# ============================================================================


def inspect_design(cellview_ref, to_file=False, output_dir=OUTPUT_DIR):
    """Print instances, parameters, pin-to-net connections and external
    pins of an existing design.

    Args:
        cellview_ref: (library, cell, view), e.g.
            ("SphericalHarmonics_lib", "Ftheta_Fphi_vars", "schematic").
        to_file: if True, write the report to
            "{library}_{cell}_{view}.txt" in output_dir (overwritten)
            instead of printing it.
        output_dir: folder for to_file output (created if missing).

    Returns:
        The opened design.
    """
    if to_file:
        library, cell, view = cellview_ref
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        output_path = Path(output_dir) / f"{library}_{cell}_{view}.txt"
        with open(output_path, "w") as f:
            with contextlib.redirect_stdout(f):
                design = _inspect_design_impl(cellview_ref)
        print(f"Output written to {output_path}")
        return design
    else:
        return _inspect_design_impl(cellview_ref)


def _inspect_design_impl(cellview_ref):
    """Print the inspection report to stdout (console or redirected file)."""
    design = de.db_uu.open_design(cellview_ref)
    print(f"=== Design: {design.design_name} ===\n")

    # ---- Instances: parameters + pins/nets ----
    print("--- Instances ---")
    for inst in design.instances:
        print(f"\nInstance: {inst.name}")

        try:
            param_names = list(inst.parameters)
        except Exception as e:
            print(f"  (could not list parameters: {e})")
            param_names = []
        for pname in param_names:
            try:
                pval = inst.parameters[pname].value
                print(f"  param {pname} = {pval}")
            except Exception as e:
                print(f"  param {pname} = <error reading value: {e}>")

        try:
            for pin in inst.inst_pins:
                try:
                    term_num = pin.inst_term.term_number
                except Exception:
                    term_num = "?"
                try:
                    net_obj = pin.net
                    net_name = getattr(net_obj, "name", net_obj)
                except Exception as e:
                    net_name = f"<error reading net: {e}>"
                print(f"  pin term{term_num} -> net {net_name}")
        except Exception as e:
            print(f"  (could not list pins: {e})")

    # ---- External (design-level) pins; attribute name varies, so try several ----
    print("\n--- External pins (design-level) ---")
    found = False
    for attr_name in ["terms", "pins", "external_pins", "design_pins", "design_terms"]:
        if hasattr(design, attr_name):
            found = True
            print(f"(found via design.{attr_name})")
            try:
                for item in getattr(design, attr_name):
                    print(f"  {item}")
            except Exception as e:
                print(f"  (attribute exists but could not iterate: {e})")
            break
    if not found:
        print("No confirmed external-pin attribute found among the "
              "candidates tried. Pin/term-related attributes on this "
              "design object:")
        print([x for x in dir(design) if 'pin' in x.lower() or 'term' in x.lower()])

    return design


print("Ready. Call inspect_design((library, cell, view)) to print to the console, or")
print("inspect_design((library, cell, view), to_file=True) to write to an")
print("auto-named text file (library_cell_view.txt) in OUTPUT_DIR, e.g.:")
print('  inspect_design(("SphericalHarmonics_lib", "Ftheta_Fphi_vars", "schematic"))')
