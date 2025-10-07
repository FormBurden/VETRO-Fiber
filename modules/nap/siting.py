# modules/nap/siting.py
from __future__ import annotations

from typing import Any, Dict, List, Tuple

from modules.geo.geometry import l_distance_feet
from modules.drops.lshape import build_l_drop_linestring
from modules.conduit.traversal import dfs_ordered_vaults

# Tunables preserved from original
MAX_L_DROP_FT = 250.0
NAP_MAX_SL = 4

def _iter_points(gj: Dict[str, Any]) -> List[Tuple[int, List[float], Dict[str, Any]]]:
    out: List[Tuple[int, List[float], Dict[str, Any]]] = []
    for i, f in enumerate(gj.get("features", [])):
        geom = f.get("geometry") or {}
        if geom.get("type") != "Point":
            continue
        coords = geom.get("coordinates")
        if not coords or len(coords) != 2:
            continue
        out.append((i, coords, f.get("properties") or {}))
    return out

def propose_naps(
    conduit_gj: Dict[str, Any],
    vaults_gj: Dict[str, Any],
    sl_gj: Dict[str, Any],
    t3_gj: Dict[str, Any],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Place NAPs at vaults in traversal order and attach up to NAP_MAX_SL nearby SLs
    with Underground L-shaped drops.
    Returns: (nap_feature_collection, drop_feature_collection)
    """
    ordered_vaults = dfs_ordered_vaults(conduit_gj, vaults_gj, t3_gj)
    sl_points = _iter_points(sl_gj)

    nap_features: List[Dict[str, Any]] = []
    drop_features: List[Dict[str, Any]] = []
    nap_counter = 0

    for _, v_xy, v_props in ordered_vaults:
        nearby: List[Tuple[float, int, List[float], Dict[str, Any]]] = []
        for si, s_xy, s_props in sl_points:
            dL = l_distance_feet(v_xy, s_xy)
            if dL <= MAX_L_DROP_FT:
                nearby.append((dL, si, s_xy, s_props))
        if not nearby:
            continue

        nearby.sort(key=lambda x: x[0])
        chosen = nearby[:NAP_MAX_SL]

        nap_counter += 1
        nap_props = {
            "v_layer": "NAP (proposed)",
            "Nap Location": "Underground",
            "NAP #": nap_counter,
            "vault_id": (v_props.get("ID") or v_props.get("Name")),
        }
        nap_features.append(
            {"type": "Feature", "properties": nap_props, "geometry": {"type": "Point", "coordinates": list(v_xy)}}
        )

        for dL, si, s_xy, s_props in chosen:
            drop_props = {
                "v_layer": "Fiber - Drop (draft)",
                "Placement": "Underground",
                "Distance (ft)": round(dL, 1),
                "NAP #": nap_counter,
                "Service Location": s_props.get("Address") or s_props.get("Name") or s_props.get("ID"),
            }
            drop_features.append(
                {
                    "type": "Feature",
                    "properties": drop_props,
                    "geometry": {"type": "LineString", "coordinates": build_l_drop_linestring(v_xy, s_xy)},
                }
            )

    return (
        {"type": "FeatureCollection", "features": nap_features},
        {"type": "FeatureCollection", "features": drop_features},
    )
