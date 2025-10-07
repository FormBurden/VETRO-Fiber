# modules/nap/io.py
from __future__ import annotations

import json
import os
from typing import Any, Dict, Tuple

from modules import config
from modules.nap.siting import propose_naps

def propose_naps_from_files(
    conduit_filename: str,
    vaults_filename: str,
    service_locations_filename: str,
    t3_vault_filename: str,
    nap_out_name: str = "nap_proposed.geojson",
    drop_out_name: str = "fiber-drop_proposed.geojson",
) -> Tuple[str, str]:
    """
    File I/O wrapper — requires a T3 vault file to anchor traversal.
    All paths are relative to modules.config.DATA_DIR.
    Returns absolute paths to the two output files under DATA_DIR.
    """
    base = config.DATA_DIR

    with open(os.path.join(base, conduit_filename), "r", encoding="utf-8") as f:
        conduit_gj = json.load(f)
    with open(os.path.join(base, vaults_filename), "r", encoding="utf-8") as f:
        vaults_gj = json.load(f)
    with open(os.path.join(base, service_locations_filename), "r", encoding="utf-8") as f:
        sl_gj = json.load(f)
    with open(os.path.join(base, t3_vault_filename), "r", encoding="utf-8") as f:
        t3_gj = json.load(f)

    nap_fc, drop_fc = propose_naps(conduit_gj, vaults_gj, sl_gj, t3_gj)

    nap_path = os.path.join(base, nap_out_name)
    drop_path = os.path.join(base, drop_out_name)
    os.makedirs(os.path.dirname(nap_path) or base, exist_ok=True)
    os.makedirs(os.path.dirname(drop_path) or base, exist_ok=True)

    with open(nap_path, "w", encoding="utf-8") as f:
        json.dump(nap_fc, f, indent=2)
    with open(drop_path, "w", encoding="utf-8") as f:
        json.dump(drop_fc, f, indent=2)

    return nap_path, drop_path
