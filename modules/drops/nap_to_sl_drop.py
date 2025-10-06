# modules/drops/nap_to_sl_drop.py
import os, re, json, datetime
from typing import Optional, Tuple, Dict, Any, List

from modules import config  # <- required: always use modules.config.DATA_DIR
from modules.fiber.colors import fiber_num_to_color_label


# ---------- helpers ----------

_RANGE_RE = re.compile(r"\((?:\d+\s*ct,?\s*)?(\d+)\s*-\s*(\d+)\)", re.IGNORECASE)

# Extract the FIRST numeric selection inside parentheses and convert it to color positions (1..12)
# We favor the first number or range that appears right after the 'ct,' anchor if present.
def _parse_primary_color_positions(nap_id: str) -> set[int]:
    """
    Returns a set of 1..12 positions derived from the primary numeric(s) in the NAP ID.
    Rules:
      • Only read the LEFT side of the parentheses content BEFORE 'Tie Point' or 'Tie-Point'
      • Accept a range X-Y or a single number X
      • Convert absolute numbers to 1..12 positions via ((n-1) % 12) + 1

    Examples:
      "(24ct, 5-7)"                                   -> {5,6,7}
      "(48ct, 30-31 / Tie Point 48ct ...)"            -> {6,7}   (reads 30-31 ONLY, ignores right side)
      "(48ct, 28 / Tie Point 48ct 29-33 to 24ct 5-9)" -> {4}     (reads 28 ONLY)
      "(24ct, 18)"                                    -> {6}

    If nothing parseable is found: returns empty set().
    """
    import re

    if not isinstance(nap_id, str):
        return set()

    # Grab inside the first parentheses: "...( <inner> )..."
    m = re.search(r"\(([^)]*)\)", nap_id)
    if not m:
        return set()
    inner = m.group(1)

    # Prefer text after the first comma (skip the 'xxct,' prefix if present)
    sub = inner.split(",", 1)[1] if "," in inner else inner

    # CRITICAL: ignore everything to the RIGHT of 'Tie Point' or 'Tie-Point'
    sub = re.split(r"\bTie[- ]?Point\b", sub, flags=re.IGNORECASE)[0]

    # First try a range on the LEFT side
    m_range = re.search(r"\b(\d+)\s*-\s*(\d+)\b", sub)
    if m_range:
        a, b = int(m_range.group(1)), int(m_range.group(2))
        if a > b:
            a, b = b, a
        return { ((n - 1) % 12) + 1 for n in range(a, b + 1) }

    # Otherwise, a single number on the LEFT side
    m_single = re.search(r"\b(\d+)\b", sub)
    if m_single:
        n = int(m_single.group(1))
        return { ((n - 1) % 12) + 1 }

    return set()

def _placement_from_nap(props: Dict[str, Any]) -> str:
    """
    Placement priority:
      1) 'Nap Location' if set (e.g., 'Underground', 'Aerial')
      2) infer from 'Type' (contains 'Aerial' or 'UG')
      3) default 'Underground'
    """
    val = props.get("Nap Location")
    if isinstance(val, str) and val.strip():
        return val.strip()
    t = (props.get("Type") or "").lower()
    if "aerial" in t:
        return "Aerial"
    if "ug" in t or "under" in t:
        return "Underground"
    return "Underground"

def _copy_meta(props: Dict[str, Any]) -> Dict[str, Any]:
    """
    Carry through project/plan fields when present so the new features
    drop cleanly into the same context; we do not invent values.
    """
    keep = ["v_project","v_project_id","v_plan","v_plan_id"]
    return {k: props[k] for k in keep if k in props}


# ---------- core builder ----------

def build_drops_from_files(
    nap_filename: str,
    service_locations_filename: str,
    out_filename: str = "drops_phase1.geojson",
) -> str:
    """
    Phase 1 engine:
      - Load NAP + Service Location GeoJSONs from modules.config.DATA_DIR
      - For each NAP, parse (X-Y) range from 'ID'
      - Match SLs with same 'NAP #' and 'Fiber #' in [X..Y]
      - Emit LineString from NAP point -> SL point
      - Props mimic 'fiber-drop' template minimally:
          v_layer='Fiber - Drop', Placement, Color
        (we also include 'NAP #' and 'Fiber #' for traceability)

    Returns the absolute output path written under DATA_DIR.
    """
    base = config.DATA_DIR

    nap_path = os.path.join(base, nap_filename)
    sl_path  = os.path.join(base, service_locations_filename)
    out_path = os.path.join(base, out_filename)

    with open(nap_path, "r", encoding="utf-8") as f:
        nap = json.load(f)
    with open(sl_path, "r", encoding="utf-8") as f:
        sl = json.load(f)

    out = _build_drops_geojson(nap, sl)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return out_path


def _build_drops_geojson(
    nap_geojson: Dict[str, Any],
    sl_geojson: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Internal: constructs the FeatureCollection (no I/O).
    - Always includes "Fiber Capacity": 72
    - Excludes v_project/v_project_id/v_plan/v_plan_id and any v_* timestamps
    - MATCH RULE: Compare by 1..12 color position. For SL Fiber # F:
        SL_pos = ((F - 1) % 12) + 1
      For a NAP, parse the FIRST number or X-Y range inside the parentheses
      and convert those absolute numbers to color positions; if SL_pos is in
      that set, it’s a match.
    """
    # Index SLs by NAP #
    sl_by_nap: Dict[int, List[Dict[str, Any]]] = {}
    for s in sl_geojson.get("features", []):
        sp = s.get("properties", {})
        nap_num = sp.get("NAP #")
        try:
            if nap_num is not None:
                nap_num = int(nap_num)
                sl_by_nap.setdefault(nap_num, []).append(s)
        except Exception:
            continue

    out_features: List[Dict[str, Any]] = []

    for nap in nap_geojson.get("features", []):
        np = nap.get("properties", {})
        nap_num = np.get("NAP #")
        try:
            if nap_num is None:
                continue
            nap_num = int(nap_num)
        except Exception:
            continue

        # NEW: parse primary color positions from the NAP ID
        nap_positions = _parse_primary_color_positions(np.get("ID", ""))
        if not nap_positions:
            # Phase 1: skip NAPs without a parseable primary number/range
            continue

        # Require NAP geometry as Point
        ngeom = nap.get("geometry", {})
        if (ngeom.get("type") != "Point") or not ngeom.get("coordinates"):
            continue
        nap_xy = ngeom["coordinates"]

        placement = _placement_from_nap(np)

        for sl in sl_by_nap.get(nap_num, []):
            sp = sl.get("properties", {})
            fnum = sp.get("Fiber #")
            try:
                if fnum is None:
                    continue
                fnum = int(fnum)
            except Exception:
                continue

            # Compare by color position (1..12) instead of absolute fiber number.
            sl_pos = ((fnum - 1) % 12) + 1
            if sl_pos not in nap_positions:
                continue

            sgeom = sl.get("geometry", {})
            if (sgeom.get("type") != "Point") or not sgeom.get("coordinates"):
                continue
            sl_xy = sgeom["coordinates"]

            props = {
                "v_layer": "Fiber - Drop",
                "Fiber Capacity": 72,
                "Placement": placement,
                "Color": fiber_num_to_color_label(fnum),  # keep the label users expect
                # Optional traceability during development (remove later if you want):
                "NAP #": nap_num,
                "Fiber #": fnum,
            }

            out_features.append({
                "type": "Feature",
                "properties": props,
                "geometry": {
                    "type": "LineString",
                    "coordinates": [nap_xy, sl_xy],
                },
            })

    return {"type": "FeatureCollection", "features": out_features}
