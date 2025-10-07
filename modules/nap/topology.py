# modules/nap/topology.py
from __future__ import annotations
from typing import Dict, Any, List, Tuple, Sequence, Optional, DefaultDict
from collections import defaultdict
import math

from modules.geo.geometry import (
    distance_feet, l_distance_feet,
    point_segment_distance_feet
)
from .constants import VAULT_ON_LINE_FT

def iter_lines(conduit_gj: Dict[str, Any]) -> List[Tuple[int, List[List[float]], Dict[str, Any]]]:
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

def iter_points(gj: Dict[str, Any]) -> List[Tuple[int, List[float], Dict[str, Any]]]:
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

def norm_val(x: Any) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip()
    return s if s and s.lower() != "null" else None

def node_key(coord: Sequence[float]) -> Tuple[int, int]:
    lon, lat = coord
    return (int(round(lon * 1e7)), int(round(lat * 1e7)))

def build_conduit_topology(conduit_gj: Dict[str, Any]):
    """
    Return:
      - feature_nodes[fid] -> [node_keys...] in geometry order
      - node_to_features[node_key] -> {fid, ...}
      - feature_props[fid] -> dict of branch-relevant props
    """
    lines = iter_lines(conduit_gj)
    feature_nodes: Dict[int, List[Tuple[int,int]]] = {}
    node_to_features: DefaultDict[Tuple[int,int], set] = defaultdict(set)
    feature_props: Dict[int, Dict[str, Any]] = {}

    for fid, coords, props in lines:
        nodes = [node_key(c) for c in coords]
        feature_nodes[fid] = nodes
        for nk in nodes:
            node_to_features[nk].add(fid)
        feature_props[fid] = {
            "Distribution Type": norm_val(props.get("Distribution Type")),
            "Section":           norm_val(props.get("Section")),
            "Fiber Letter #1":   norm_val(props.get("Fiber Letter #1")),
            "Fiber Letter #2":   norm_val(props.get("Fiber Letter #2")),
            "Conduit Type":      norm_val(props.get("Conduit Type")),
            "ID":                norm_val(props.get("ID")),
        }

    # vaults_on_feature is computed by assign_vaults_to_features (needs vaults)
    vaults_on_feature: DefaultDict[int, List[Tuple[float,int,List[float],Dict[str,Any]]]] = defaultdict(list)
    return feature_nodes, node_to_features, feature_props, vaults_on_feature

def assign_vaults_to_features(
    conduit_gj: Dict[str, Any],
    feature_nodes: Dict[int, List[Tuple[int,int]]],
    vaults_gj: Dict[str, Any]
) -> DefaultDict[int, List[Tuple[float,int,List[float],Dict[str,Any]]]]:
    """
    Map each vault to nearest conduit polyline if within VAULT_ON_LINE_FT.
    Store (along_ft, vault_index, vault_coord, vault_props) sorted by along.
    """
    lines = {fid:(coords, props) for fid, coords, props in iter_lines(conduit_gj)}

    def nearest_on_feature(fid: int, p: List[float]) -> Tuple[float, float]:
        coords, _ = lines[fid]
        best = (float("inf"), 0.0)
        along = 0.0
        for i in range(1, len(coords)):
            a, b = coords[i-1], coords[i]
            d = point_segment_distance_feet(p, a, b)
            if d < best[0]:
                da = distance_feet(a, p)
                db = distance_feet(b, p)
                seg_len = distance_feet(a, b)
                best = (d, along + (0.0 if da <= db else seg_len))
            along += distance_feet(a, b)
        return best  # (min_d_ft, along_ft)

    vaults = iter_points(vaults_gj)
    v_by_feature: DefaultDict[int, List[Tuple[float,int,List[float],Dict[str,Any]]]] = defaultdict(list)
    for vi, vpt, vprops in vaults:
        best = (float("inf"), None, 0.0)
        for fid in feature_nodes.keys():
            d_ft, along_ft = nearest_on_feature(fid, vpt)
            if d_ft < best[0]:
                best = (d_ft, fid, along_ft)
        if best[1] is not None and best[0] <= VAULT_ON_LINE_FT:
            v_by_feature[best[1]].append((best[2], vi, vpt, vprops))

    for fid, items in v_by_feature.items():
        items.sort(key=lambda t: t[0])

    return v_by_feature

def vec_meters(p1: Sequence[float], p2: Sequence[float]) -> Tuple[float, float]:
    lon1, lat1 = p1
    lon2, lat2 = p2
    lat_mid = (lat1 + lat2) / 2.0
    m_per_deg_lat = 111132.92 - 559.82 * math.cos(math.radians(2 * lat_mid)) + 1.175 * math.cos(math.radians(4 * lat_mid)) - 0.0023 * math.cos(math.radians(6 * lat_mid))
    m_per_deg_lon = (111412.84 * math.cos(math.radians(lat_mid)) - 93.5 * math.cos(math.radians(3 * lat_mid)) + 0.118 * math.cos(math.radians(5 * lat_mid)))
    return ((lon2 - lon1) * m_per_deg_lon, (lat2 - lat1) * m_per_deg_lat)

def find_start_node_from_t3(
    t3_gj: Dict[str, Any],
    feature_nodes: Dict[int, List[Tuple[int,int]]]
) -> Optional[Tuple[int,int]]:
    t3_pts: List[List[float]] = []
    for f in t3_gj.get("features", []):
        g = f.get("geometry") or {}
        if g.get("type") == "Point" and isinstance(g.get("coordinates"), (list, tuple)) and len(g["coordinates"]) == 2:
            t3_pts.append(g["coordinates"])

    if not t3_pts:
        return None

    all_nodes = []
    for nodes in feature_nodes.values():
        all_nodes.extend(nodes)
    if not all_nodes:
        return None
    all_nodes = list(dict.fromkeys(all_nodes))

    def centroid(nk: Tuple[int,int]) -> List[float]:
        return [nk[0] / 1e7, nk[1] / 1e7]

    t3 = t3_pts[0]
    best = (float("inf"), None)
    for nk in all_nodes:
        d = distance_feet(t3, centroid(nk))
        if d < best[0]:
            best = (d, nk)
    return best[1]
