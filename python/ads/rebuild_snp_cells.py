r"""
rebuild_snp_cells.py -- refresh the antenna cells of one or more Touchstone files.

Defines rebuild_snp_cells(snp_paths, ...), which, for each SnP file, without
touching any Tx/Rx link cell:
  - rebuilds "{base_name}_TOUCHESTONE" with create_snp_component_cell
    (overwriting manual edits);
  - if the Rational Fit cell "{base_name}" (created by File > Import >
    Netlist File, see create_snp_dual_view_cells.py) exists, redraws its
    symbol with redraw_symbol_for_imported_cell, and optionally replaces
    its netlist: pass netlist_paths={snp_path: raw_spice_netlist}; the
    netlist is translated with nettrans, checked against the cell
    (subcircuit name == base_name, node count == number of terminals) and
    copied over the NetlistInclude target (old file kept as .bak). The
    netlist refresh has not been tested in ADS. A missing Rational Fit
    cell is only reported; it cannot be created from Python.
base_name is get_gui_import_base_name(snp_path) (lowercase filename with
"." -> "_"), as used by the link scripts.

workspace_path=None uses the open workspace; otherwise any other open
workspace is closed (unsaved work is not saved) and workspace_path opened.
A redrawn symbol that is already placed may keep stale pin positions:
reopen the cell windows and re-place the instance if the port list changed.

Usage (ADS Python console, after start_console.py):
    run_script("rebuild_snp_cells.py")
    rebuild_snp_cells(r"<repo>\output\SnP\S_dom_psi_<antenna>.s14p")
    rebuild_snp_cells([path1, path2], workspace_path=r"C:\path\to\your\ADS_workspace_wrk",
                      library_name="SphericalHarmonics_lib")
    rebuild_snp_cells(path1, netlist_paths={path1: r"<repo>\output\SnP\S_dom_psi_<antenna>.s14p.ckt"})
"""

from pathlib import Path
import keysight.ads.de as de

# Repository folders (defined by start_console.py; derived here if this file is exec()'d directly).
try:
    SCRIPT_DIR
except NameError:
    SCRIPT_DIR = Path(THIS_SCRIPT_PATH).resolve().parent

exec(open(SCRIPT_DIR / "create_snp_dual_view_cells.py").read())


def _open_workspace(workspace_path):
    """Returns the workspace to use, opening/switching to workspace_path
    if given (None = whatever is already open)."""
    if workspace_path is None:
        if not de.workspace_is_open():
            raise RuntimeError("No workspace is open -- pass workspace_path.")
        return de.active_workspace()
    wanted = Path(workspace_path)
    if not wanted.exists():
        raise FileNotFoundError(f"Workspace path does not exist: {wanted}")
    if de.workspace_is_open():
        ws = de.active_workspace()
        if wanted.name.lower() in str(ws).lower():
            print(f"Workspace already open: {ws}")
            return ws
        print(f"A different workspace is open ({ws}) -- closing it to open {wanted}.")
        ws.close()
    ws = de.open_workspace(wanted)
    print(f"Opened workspace: {ws}")
    return ws


def _refresh_rational_fit_netlist(library_name, base_name, snp_path, raw_netlist_path):
    """Translates raw_netlist_path (nettrans) and overwrites the file the
    existing Rational Fit cell's NetlistInclude points at. Checks the new
    netlist's define-block name/node count against the cell first."""
    import shutil
    design = de.db_uu.open_design((library_name, base_name, "schematic"))
    try:
        n_terms = len(design.terms)
        inc = next((i for i in design.instances if "NetlistInclude" in str(i.master_lcv_name)), None)
        if inc is None:
            raise RuntimeError(f"No NetlistInclude in {library_name}:{base_name}.")
        target = Path(inc.parameters["IncludePath"].value[0]) / inc.parameters["IncludeFiles"].value[0][0]
    finally:
        design.close_design()
    tmp = Path(raw_netlist_path).with_name(f"{Path(raw_netlist_path).stem.replace('.', '_')}_rebuild_ads.net")
    translated = translate_netlist_to_ads(raw_netlist_path, tmp)
    blocks = parse_ads_netlist_define_blocks(translated)
    good = [b for b in blocks if len(b[1]) == n_terms]
    if not good:
        raise ValueError(f"New netlist has define blocks with {[len(b[1]) for b in blocks]} nodes; "
                         f"cell '{base_name}' has {n_terms} terminals -- port list changed? Re-import via GUI.")
    if good[0][0].lower() != base_name.lower():
        raise ValueError(f"New netlist's subcircuit is '{good[0][0]}', but the cell is '{base_name}' -- "
                         f"the name must match; regenerate with the matching filename.")
    if target.exists():
        shutil.copy2(target, str(target) + ".bak")
    shutil.copy2(translated, target)
    print(f"Refreshed netlist {target} (old copy: {target}.bak).")


def rebuild_snp_cells(snp_paths, library_name="SphericalHarmonics_lib", workspace_path=None,
                      rebuild_touchstone=True, redraw_rational_fit=True, netlist_paths=None):
    """Rebuilds {base}_TOUCHESTONE and/or redraws the {base} Rational Fit
    symbol for each Touchstone file in snp_paths (a path or a list of
    paths). netlist_paths: optional {snp_path: raw_netlist_path} to also
    refresh the Rational Fit cell's netlist. Returns {base_name:
    {"touchstone": bool, "rational_fit": str}}, where rational_fit is
    "redrawn", "netlist+redrawn", "missing" or "skipped"."""
    if isinstance(snp_paths, (str, Path)):
        snp_paths = [snp_paths]
    netlist_paths = {str(k): v for k, v in (netlist_paths or {}).items()}
    ws = _open_workspace(workspace_path)
    if library_name not in ws.library_names:
        raise RuntimeError(f"Library '{library_name}' not found.")
    library = ws.open_library(library_name)

    results = {}
    for snp_path in snp_paths:
        base_name = get_gui_import_base_name(snp_path)
        if base_name in results:
            print(f"'{base_name}' already handled in this run -- skipping duplicate.")
            continue
        snp_cell = f"{base_name}_TOUCHESTONE"
        res = {"touchstone": False, "rational_fit": "skipped"}
        print(f"\n===== {Path(snp_path).name} -> base name '{base_name}' =====")

        if rebuild_touchstone:
            if library.cell_exists(snp_cell):
                close_cell_views(library_name, snp_cell)
            create_snp_component_cell(snp_path, library_name=library_name, new_cell_name=snp_cell)
            res["touchstone"] = True
            print(f"Rebuilt '{snp_cell}'.")

        if redraw_rational_fit:
            if library.cell_exists(base_name):
                raw_netlist = netlist_paths.get(str(snp_path))
                if raw_netlist:
                    _refresh_rational_fit_netlist(library_name, base_name, snp_path, raw_netlist)
                redraw_symbol_for_imported_cell(library_name, base_name, snp_path=snp_path)
                res["rational_fit"] = "netlist+redrawn" if raw_netlist else "redrawn"
                print(f"Redrew symbol of Rational Fit cell '{base_name}'.")
            else:
                res["rational_fit"] = "missing"
                print(f"No Rational Fit cell '{base_name}' in '{library_name}' -- create it with "
                      f"File > Import > Netlist File (name must stay '{base_name}'), then re-run.")
        results[base_name] = res

    print("\nSummary:")
    for base_name, res in results.items():
        print(f"  {base_name}: TOUCHESTONE rebuilt={res['touchstone']}, Rational Fit {res['rational_fit']}")
    return results


print("Ready: rebuild_snp_cells(snp_paths, library_name='SphericalHarmonics_lib', workspace_path=None, "
      "rebuild_touchstone=True, redraw_rational_fit=True, netlist_paths=None). See the module docstring.")
