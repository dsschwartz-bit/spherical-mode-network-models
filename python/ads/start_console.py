r"""
start_console.py -- run first in every ADS Python console session.

Imports the ADS design-environment modules used by the other scripts,
records where this repository lives, and opens the ADS workspace set in
USER SETTINGS below.

Usage (ADS Python console, Ctrl+Shift+P), replacing <repo> with your clone:

    exec(open(THIS_SCRIPT_PATH := r"<repo>\python\ads\start_console.py").read())

After that, load any other script in this folder with run_script(), e.g.

    run_script("create_tx_rx_link.py")
"""

from pathlib import Path
import re
import keysight.ads.de as de
from keysight.ads.de import db_uu
from keysight.ads.de.experimental import generate_symbol as gs

# ============================== USER SETTINGS ==============================
# ADS workspace to open (created beforehand in ADS).
WORKSPACE_PATH = Path(r"C:\path\to\your\ADS_workspace_wrk")
# ===========================================================================

# Repository locations, derived from this file's own path (exec() does not
# set __file__, so the exec line above passes it in as THIS_SCRIPT_PATH).
SCRIPT_DIR = Path(THIS_SCRIPT_PATH).resolve().parent  # noqa: F821
REPO_DIR = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_DIR / "mathematica" / "tables"


def run_script(name):
    """Execute another script from python/ads/ in the console's namespace."""
    global THIS_SCRIPT_PATH
    THIS_SCRIPT_PATH = str(SCRIPT_DIR / name)
    exec(open(THIS_SCRIPT_PATH).read(), globals())


if de.workspace_is_open():
    ws = de.active_workspace()
    print(f"Workspace already open: {ws}")
elif WORKSPACE_PATH.exists():
    ws = de.open_workspace(WORKSPACE_PATH)
    print(f"Opened workspace: {ws}")
else:
    ws = None
    print(f"WORKSPACE_PATH not found: {WORKSPACE_PATH}. Edit USER SETTINGS "
          f"in start_console.py or open a workspace with de.open_workspace().")

print(f"Ready. SCRIPT_DIR = {SCRIPT_DIR}; DATA_DIR = {DATA_DIR}")
