# modules/conduit/nap_siting.py
from __future__ import annotations
import os, json
from typing import Dict, Any, List, Tuple, Sequence, Optional

from modules import config
from modules.geo.geometry import (
    distance_feet, l_distance_feet,
    point_segment_distance_feet, project_point_onto_polyline_feet
)

# Tunables (project-style constants)
MAX_L_DROP_FT = 250.0
VAULT_ON_LINE_FT = 15.0
NAP_MAX_SL = 4


def _iter_lines(conduit_gj: Dict[str, Any]) -> List[Tuple[int, List[List[float]], Dict[str, Any]]]:
    out = []
    for i, f in enumerate(conduit_gj.get("features", [])):
        geom = f.get("geometry") or {}
        if geom.get("type") != "LineString":
            continue
        coords = geom.get("coordinates") or []
        if len(coords) < 2:
            continue
        out.append((i, coords, f.get("properties") or {}))
    return out


def _iter_points(gj: Dict[str, Any]) -> List[Tuple[int, List[float], Dict[str, Any]]]:
    out = []
    for i, f in enumerate(gj.get("features", [])):
        geom = f.get("geometry") or {}
        if geom.get("type") != "Point":
            continue
        coords = geom.get("coordinates")
        if not coords or len(coords) != 2:
            continue
        out.append((i, coords, f.get("properties") or {}))
    return out


def _vaults_on_conduit(conduit_gj: Dict[str, Any], vaults_gj: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Return vaults that lie on a conduit segment (within VAULT_ON_LINE_FT),
    annotated with which conduit feature they sit on and their order along the line.
    """
    lines = _iter_lines(conduit_gj)
    vaults = _iter_points(vaults_gj)

    results: List[Dict[str, Any]] = []
    for vi, vpt, vprops in vaults:
        best = None
        for li, line_coords, lprops in lines:
            # Distance from point to the polyline (min over segments)
            min_d = float("inf")
            for j in range(1, len(line_coords)):
                a, b = line_coords[j-1], line_coords[j]
                d = point_segment_distance_feet(vpt, a, b)
                if d < min_d:
                    min_d = d
            if min_d <= VAULT_ON_LINE_FT:
                along_ft, _ = project_point_onto_polyline_feet(vpt, line_coords)
                best = (li, along_ft, line_coords, lprops)
                break
        if best:
            li, along_ft, line_coords, lprops = best
            results.append({
                "vault_index": vi,
                "coords": vpt,
                "properties": vprops,
                "conduit_index": li,
                "conduit_props": lprops,
                "along_ft": along_ft
            })
    # Order by line and along distance
    results.sort(key=lambda r: (r["conduit_index"], r["along_ft"]))
    return results


def _build_l_linestring(a: Sequence[float], b: Sequence[float]) -> List[List[float]]:
    """
    Build an L-shaped polyline [a -> intermediate -> b] choosing the
    elbow that yields the same total L-distance but looks tidy:
    first horizontal, then vertical (lon first, then lat).
    """
    ax, ay = a
    bx, by = b
    elbow = [bx, ay]  # move in X, then Y
    return [list(a), elbow, list(b)]


def _propose_naps_and_drops(
    conduit_gj: Dict[str, Any],
    vaults_gj: Dict[str, Any],
    sl_gj: Dict[str, Any]
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Very first-pass siting:
      1) Find vaults on conduit lines.
      2) For each such vault, attach up to 4 SLs within 250 ft (L distance).
      3) Emit a NAP point per qualified vault and draft Fiber-Drop lines.

    We treat these as Underground by default (conduit implies UG).
    """
    vaults_on = _vaults_on_conduit(conduit_gj, vaults_gj)
    sl_points = _iter_points(sl_gj)

    # Build quick index of SLs that are within MAX_L_DROP_FT of each vault
    nap_features: List[Dict[str, Any]] = []
    drop_features: List[Dict[str, Any]] = []

    # Simple counter for NAP # proposal (scoped to this run)
    nap_counter = 0

    for v in vaults_on:
        v_xy = v["coords"]
        nearby: List[Tuple[float, int, List[float], Dict[str, Any]]] = []
        for si, s_xy, s_props in sl_points:
            dL = l_distance_feet(v_xy, s_xy)
            if dL <= MAX_L_DROP_FT:
                nearby.append((dL, si, s_xy, s_props))
        if not nearby:
            continue

        # Sort by L-distance, take up to NAP_MAX_SL
        nearby.sort(key=lambda x: x[0])
        chosen = nearby[:NAP_MAX_SL]

        # Propose a NAP at this vault
        nap_counter += 1
        nap_props = {
            "v_layer": "NAP (proposed)",
            "Nap Location": "Underground",
            "NAP # (proposed)": nap_counter,
            # Carry through some project/plan context if present on conduit:
            "conduit_id": (v["conduit_props"] or {}).get("ID"),
        }
        nap_features.append({
            "type": "Feature",
            "properties": nap_props,
            "geometry": {"type": "Point", "coordinates": list(v_xy)}
        })

        # Draft drops for the chosen SLs
        for dL, si, s_xy, s_props in chosen:
            drop_props = {
                "v_layer": "Fiber - Drop (draft)",
                "Placement": "Underground",
                "Distance (ft)": round(dL, 1),
                "NAP # (proposed)": nap_counter,
                "Service Location": s_props.get("Address") or s_props.get("Name") or s_props.get("ID"),
            }
            drop_features.append({
                "type": "Feature",
                "properties": drop_props,
                "geometry": {
                    "type": "LineString",
                    "coordinates": _build_l_linestring(v_xy, s_xy)
                }
            })

    nap_fc = {"type": "FeatureCollection", "features": nap_features}
    drop_fc = {"type": "FeatureCollection", "features": drop_features}
    return nap_fc, drop_fc


def propose_naps_from_files(
    conduit_filename: str,
    vaults_filename: str,
    service_locations_filename: str,
    nap_out_name: str = "nap_proposed.geojson",
    drop_out_name: str = "fiber-drop_proposed.geojson"
) -> Tuple[str, str]:
    """
    I/O wrapper:
      - Reads conduit*, vaults*, and service-location* from modules.config.DATA_DIR
      - Writes proposed 'nap_proposed.geojson' and 'fiber-drop_proposed.geojson'
    Returns (nap_path, drop_path)
    """
    base = config.DATA_DIR
    with open(os.path.join(base, conduit_filename), "r", encoding="utf-8") as f:
        conduit_gj = json.load(f)
    with open(os.path.join(base, vaults_filename), "r", encoding="utf-8") as f:
        vaults_gj = json.load(f)
    with open(os.path.join(base, service_locations_filename), "r", encoding="utf-8") as f:
        sl_gj = json.load(f)

    nap_fc, drop_fc = _propose_naps_and_drops(conduit_gj, vaults_gj, sl_gj)

    nap_path = os.path.join(base, nap_out_name)
    drop_path = os.path.join(base, drop_out_name)
    with open(nap_path, "w", encoding="utf-8") as f:
        json.dump(nap_fc, f, indent=2)
    with open(drop_path, "w", encoding="utf-8") as f:
        json.dump(drop_fc, f, indent=2)

    return nap_path, drop_path
