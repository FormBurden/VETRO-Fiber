# modules/conduit/validator.py
# Validate conduit feature attributes and emit a simple QA report + optional annotated GeoJSON.
# ALWAYS read/write under modules.config.DATA_DIR (project rule).

from __future__ import annotations
import os, json, re, datetime
from typing import Dict, Any, List, Optional, Tuple

from modules import config  # <- use modules.config.DATA_DIR


_ALLOWED_DT = {"24ct", "48ct", "72ct"}
_ALLOWED_SECTION = {f"Section {i}" for i in range(1, 7)}

_CONDUIT_TYPE_RE = re.compile(
    r'^\s*(?P<count>\d+)\s*x\s*1\.25\s*"?\s*(?P<roadshot>Road\s*Shot)?\s*$', re.IGNORECASE
)


def _norm(x: Any) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip()
    return s if s and s.lower() != "null" else None


def _parse_conduit_type(s: Optional[str]) -> Dict[str, Any]:
    """
    Returns a dict:
      {recognized: bool, n: int|None, is_2x: bool, is_road_shot: bool, raw: str|None}
    Only recognizes the project-approved:
      "1 x 1.25\"", "2 x 1.25\"", "2 x 1.25\" Road Shot"
    """
    raw = s
    s = _norm(s) or ""
    m = _CONDUIT_TYPE_RE.match(s)
    if not m:
        return {"recognized": False, "n": None, "is_2x": False, "is_road_shot": False, "raw": raw}
    n = int(m.group("count"))
    return {
        "recognized": True,
        "n": n,
        "is_2x": (n == 2),
        "is_road_shot": bool(m.group("roadshot")),
        "raw": raw,
    }


def validate_conduit_feature(props: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Validate ONE conduit feature's attributes using the (corrected) rules:

      • If Distribution Type 2 is set → Conduit Type must be 2× OR 2× Road Shot.
      • If Conduit Type is 2× and Distribution Type 2 is null → Backhaul must be null.
      • If Backhaul Distribution Type is set → Backhaul must be “Yes” AND Conduit Type must be 2× (Road Shot allowed).

    Notes:
      - 'null' (string) and JSON null are treated as None.
      - Section and basic DT value checks are still validated as before.
    """
    issues: List[Dict[str, Any]] = []

    cid = _norm(props.get("ID")) or "(no ID)"

    ct_info = _parse_conduit_type(props.get("Conduit Type"))
    dt1 = _norm(props.get("Distribution Type"))
    dt2 = _norm(props.get("Distribution Type 2"))
    bh  = _norm(props.get("Backhaul"))
    bhdt = _norm(props.get("Backhaul Distribution Type"))
    section = _norm(props.get("Section"))

    # --- Basic value checks (unchanged) ---
    if not ct_info["recognized"]:
        issues.append({
            "severity": "warn",
            "code": "conduit_type_unrecognized",
            "msg": f"{cid}: Conduit Type '{props.get('Conduit Type')}' not in recognized set "
                   f'["1 x 1.25\"", "2 x 1.25\"", "2 x 1.25\" Road Shot"]. This feature may be ignored by automation.'
        })

    if dt1 and dt1 not in _ALLOWED_DT:
        issues.append({"severity": "error", "code": "dt1_invalid",
                       "msg": f"{cid}: Distribution Type '{dt1}' must be one of {_ALLOWED_DT}."})

    if dt2 and dt2 not in _ALLOWED_DT:
        issues.append({"severity": "error", "code": "dt2_invalid",
                       "msg": f"{cid}: Distribution Type 2 '{dt2}' must be one of {_ALLOWED_DT}."})

    if bhdt and bhdt not in _ALLOWED_DT:
        issues.append({"severity": "error", "code": "backhaul_dt_invalid",
                       "msg": f"{cid}: Backhaul Distribution Type '{bhdt}' must be one of {_ALLOWED_DT}."})

    if section and section not in _ALLOWED_SECTION:
        issues.append({"severity": "warn", "code": "section_unexpected",
                       "msg": f"{cid}: Section '{section}' is outside the expected 'Section 1'..'Section 6'."})

    # --- Corrected core rules ---

    # A) If Distribution Type 2 is set → Conduit Type must be 2× (Road Shot allowed).
    if dt2 and not ct_info["is_2x"]:
        issues.append({"severity": "error", "code": "dt2_requires_2x",
                       "msg": f"{cid}: Distribution Type 2 present → Conduit Type must be '2 x 1.25\"' "
                              f"(Road Shot allowed). Found '{props.get('Conduit Type')}'."})

    # B) If Conduit Type is 2× and Distribution Type 2 is null → Backhaul must be null.
    if ct_info["is_2x"] and (dt2 is None):
        if bh is not None:
            issues.append({"severity": "error", "code": "2x_dt2_null_backhaul_must_null",
                           "msg": f"{cid}: 2× conduit with no 'Distribution Type 2' requires Backhaul to be null (not set)."})

    # C) If Backhaul Distribution Type is set → Backhaul must be “Yes” AND Conduit Type must be 2× (Road Shot allowed).
    if bhdt:
        if (bh or "").lower() != "yes":
            issues.append({"severity": "error", "code": "backhaul_dt_requires_yes",
                           "msg": f"{cid}: 'Backhaul Distribution Type' present → Backhaul must be 'Yes'."})
        if not ct_info["is_2x"]:
            issues.append({"severity": "error", "code": "backhaul_dt_requires_2x",
                           "msg": f"{cid}: 'Backhaul Distribution Type' present → Conduit Type must be 2× "
                                  f"(Road Shot allowed). Found '{props.get('Conduit Type')}'."})

    # Advisory: recommend DT1 is filled for any conduit considered by automation.
    if not dt1:
        issues.append({"severity": "warn", "code": "dt1_missing",
                       "msg": f"{cid}: 'Distribution Type' is empty; recommend setting it for clarity."})

    return issues



def validate_conduit_file(
    conduit_filename: str,
    *,
    write_annotated_geojson: bool = True,
    annotated_out_name: str = "conduit_validated.geojson",
    report_out_name: Optional[str] = None
) -> Dict[str, Any]:
    """
    Load conduit GeoJSON from modules.config.DATA_DIR and validate all features.
    If write_annotated_geojson=True, writes a new GeoJSON that adds:
      props.qa_ok (bool) and props.qa_issues (list of strings).
    Also writes a JSON QA report (defaults to 'conduit_qa-YYYYmmdd-HHMMSS.json').

    Returns a dict:
      {
        "total": int,
        "considered": int,    # features with recognized conduit type
        "ok": int,
        "issues": int,
        "report_path": str,
        "annotated_path": Optional[str]
      }
    """
    base = config.DATA_DIR
    in_path = os.path.join(base, conduit_filename)

    with open(in_path, "r", encoding="utf-8") as f:
        gj = json.load(f)

    features = gj.get("features", [])
    out_features = []
    report: List[Dict[str, Any]] = []

    total = len(features)
    considered = 0
    ok = 0
    issues_count = 0

    for feat in features:
        props = feat.get("properties", {}) or {}
        iss = validate_conduit_feature(props)

        ct_info = _parse_conduit_type(props.get("Conduit Type"))
        if ct_info["recognized"]:
            considered += 1

        qa_ok = len([i for i in iss if i["severity"] == "error"]) == 0
        if qa_ok:
            ok += 1
        issues_count += len(iss)

        # For annotated GeoJSON, embed compact messages (strings)
        if write_annotated_geojson:
            ap = dict(props)
            ap["qa_ok"] = qa_ok
            ap["qa_issues"] = [f"[{i['severity']}] {i['code']}: {i['msg']}" for i in iss]
            out_features.append({"type": "Feature", "properties": ap, "geometry": feat.get("geometry")})

        # For the report, keep structured rows
        report.append({
            "id": props.get("ID"),
            "Conduit Type": props.get("Conduit Type"),
            "Distribution Type": props.get("Distribution Type"),
            "Distribution Type 2": props.get("Distribution Type 2"),
            "Backhaul": props.get("Backhaul"),
            "Backhaul Distribution Type": props.get("Backhaul Distribution Type"),
            "Section": props.get("Section"),
            "qa_ok": qa_ok,
            "issues": iss,
        })

    # Write annotated GeoJSON (optional)
    annotated_path = None
    if write_annotated_geojson:
        annotated = {"type": "FeatureCollection", "features": out_features}
        annotated_path = os.path.join(base, annotated_out_name)
        os.makedirs(os.path.dirname(annotated_path), exist_ok=True)
        with open(annotated_path, "w", encoding="utf-8") as f:
            json.dump(annotated, f, indent=2)

    # Write QA report
    ts = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    if not report_out_name:
        report_out_name = f"conduit_qa-{ts}.json"
    report_path = os.path.join(base, report_out_name)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump({
            "summary": {
                "total": total,
                "considered": considered,
                "ok": ok,
                "issues": issues_count
            },
            "rows": report
        }, f, indent=2)

    return {
        "total": total, "considered": considered, "ok": ok, "issues": issues_count,
        "report_path": report_path, "annotated_path": annotated_path
    }
