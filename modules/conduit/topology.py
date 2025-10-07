# modules/conduit/topology.py
from __future__ import annotations

from typing import Any, DefaultDict, Dict, List, Optional, Sequence, Tuple
from collections import defaultdict

from modules.geo.geometry import distance_feet, point_segment_distance_feet

# Snap tolerance for vault-to-line association (feet)
VAULT_ON_LINE_FT = 15.0


def iter_lines(conduit_gj: Dict[str, Any]) -> List[Tuple[int, List[List[float]], Dict[str, Any]]]:
    out: List[Tuple[int, List[List[float]], Dict[str, Any]]] = []
    for i, f in enumerate(conduit_gj.get("features", [])):
        geom = f.get("geometry") or {}
        if geom.get("type") != "LineString":
            continue
        coords = geom.get("coordinates") or []
        if len(coords) < 2:
            continue
        out.append((i, coords, f.get("properties") or {}))
    return out


def iter_points(gj: Dict[str, Any]) -> List[Tuple[int, List[float], Dict[str, Any]]]:
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


def _norm_val(x: Any) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip()
    return s if s and s.lower() != "null" else None


def node_key(coord: Sequence[float]) -> Tuple[int, int]:
    """Quantize lon/lat to ~1e-7 deg for robust vertex matching."""
    lon, lat = coord
    return (int(round(lon * 1e7)), int(round(lat * 1e7)))


def build_conduit_topology(conduit_gj: Dict[str, Any]):
    """
    Returns:
      feature_nodes[fid] -> [node_keys ...] in geometry order
      node_to_features[node_key] -> {fid, ...}
      feature_props[fid] -> selected props for traversal rules
    """
    lines = iter_lines(conduit_gj)

    feature_nodes: Dict[int, List[Tuple[int, int]]] = {}
    node_to_features: DefaultDict[Tuple[int, int], set] = defaultdict(set)
    feature_props: Dict[int, Dict[str, Any]] = {}

    for fid, coords, props in lines:
        nodes = [node_key(c) for c in coords]
        feature_nodes[fid] = nodes
        for nk in nodes:
            node_to_features[nk].add(fid)

        feature_props[fid] = {
            "Distribution Type": _norm_val(props.get("Distribution Type")),
            "Section": _norm_val(props.get("Section")),
            "Fiber Letter #1": _norm_val(props.get("Fiber Letter #1")),
            "Fiber Letter #2": _norm_val(props.get("Fiber Letter #2")),
            "Conduit Type": _norm_val(props.get("Conduit Type")),
            "ID": _norm_val(props.get("ID")),
        }

    return feature_nodes, node_to_features, feature_props


def assign_vaults_to_features(
    conduit_gj: Dict[str, Any],
    feature_nodes: Dict[int, List[Tuple[int, int]]],
    vaults_gj: Dict[str, Any],
) -> DefaultDict[int, List[Tuple[float, int, List[float], Dict[str, Any]]]]:
    """
    Associate each vault point to the nearest feature polyline (if within VAULT_ON_LINE_FT).
    Returns v_by_feature[fid] -> list of (along_ft, vault_index, vault_coord, vault_props), sorted by along_ft.
    """
    lines = {fid: (coords, props) for fid, coords, props in iter_lines(conduit_gj)}

    def _nearest_on_feature(fid: int, p: List[float]) -> Tuple[float, float]:
        coords, _ = lines[fid]
        best = (float("inf"), 0.0)
        along = 0.0
        for i in range(1, len(coords)):
            a, b = coords[i - 1], coords[i]
            d = point_segment_distance_feet(p, a, b)
            if d < best[0]:
                # crude along by closer endpoint
                da = distance_feet(a, p)
                db = distance_feet(b, p)
                seg_len = distance_feet(a, b)
                best = (d, along + (0.0 if da <= db else seg_len))
            along += distance_feet(a, b)
        return best  # (min_d_ft, along_ft)

    v_by_feature: DefaultDict[int, List[Tuple[float, int, List[float], Dict[str, Any]]]] = defaultdict(list)
    for vi, vpt, vprops in iter_points(vaults_gj):
        best = (float("inf"), None, 0.0)  # (dist_ft, fid, along_ft)
        for fid in feature_nodes.keys():
            d_ft, along_ft = _nearest_on_feature(fid, vpt)
            if d_ft < best[0]:
                best = (d_ft, fid, along_ft)
        if best[1] is not None and best[0] <= VAULT_ON_LINE_FT:
            v_by_feature[best[1]].append((best[2], vi, vpt, vprops))

    for fid, items in v_by_feature.items():
        items.sort(key=lambda t: t[0])
    return v_by_feature
