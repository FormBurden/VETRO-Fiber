from __future__ import annotations
import os, json, math
from typing import Dict, Any, List, Tuple, Sequence, Optional, DefaultDict
from collections import defaultdict

from modules import config
from modules.geo.geometry import (
    distance_feet, l_distance_feet,
    point_segment_distance_feet
)

# Tunables
MAX_L_DROP_FT = 250.0
VAULT_ON_LINE_FT = 15.0
NODE_SNAP_FT = 20.0
NAP_MAX_SL = 4


# ------------------------ basic iterators ------------------------

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


# ------------------------ topology helpers ------------------------

def _norm_val(x: Any) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip()
    return s if s and s.lower() != "null" else None

def _node_key(coord: Sequence[float]) -> Tuple[int, int]:
    """Quantize lon/lat to ~1e-7 deg (~1 cm-ish) for robust vertex matching."""
    lon, lat = coord
    return (int(round(lon * 1e7)), int(round(lat * 1e7)))

def _build_conduit_topology(conduit_gj: Dict[str, Any]):
    """
    Build:
      - feature_nodes[fid] -> [node_keys...] in geometry order
      - node_to_features[node_key] -> {fid, ...}
      - feature_props[fid] -> dict with keys we branch on:
            'Distribution Type', 'Section', 'Fiber Letter #1', 'Fiber Letter #2', 'Conduit Type', 'ID'
      - vaults_on_feature[fid] -> list of (along_ft, vault_index, vault_coord, vault_props)
    Vaults are associated to the *nearest segment* of a feature if within VAULT_ON_LINE_FT.
    """
    lines = _iter_lines(conduit_gj)
    feature_nodes: Dict[int, List[Tuple[int,int]]] = {}
    node_to_features: DefaultDict[Tuple[int,int], set] = defaultdict(set)
    feature_props: Dict[int, Dict[str, Any]] = {}

    for fid, coords, props in lines:
        nodes = [_node_key(c) for c in coords]
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

    # Pre-compute vault assignment to features/along
    def _polyline_along_feet(coords: List[List[float]], p: List[float]) -> Tuple[float, float]:
        """Return (min_distance_ft, along_ft) from p to the polyline."""
        best = (float("inf"), 0.0)
        along = 0.0
        for i in range(1, len(coords)):
            a, b = coords[i-1], coords[i]
            d = point_segment_distance_feet(p, a, b)
            if d < best[0]:
                # crude along: decide by closer endpoint
                da = distance_feet(a, p)
                db = distance_feet(b, p)
                seg_len = distance_feet(a, b)
                best = (d, along + (0.0 if da <= db else seg_len))
            along += distance_feet(a, b)
        return best

    vaults_on_feature: DefaultDict[int, List[Tuple[float,int,List[float],Dict[str,Any]]]] = defaultdict(list)
    # Fill later in _assign_vaults_to_features (needs vaults)
    return feature_nodes, node_to_features, feature_props, vaults_on_feature

def _assign_vaults_to_features(
    conduit_gj: Dict[str, Any],
    feature_nodes: Dict[int, List[Tuple[int,int]]],
    vaults_gj: Dict[str, Any]
) -> DefaultDict[int, List[Tuple[float,int,List[float],Dict[str,Any]]]]:
    """
    For each vault, find the nearest feature polyline within VAULT_ON_LINE_FT and
    compute its along-distance (feet) from the feature's start vertex.
    """
    lines = {fid:(coords, props) for fid, coords, props in _iter_lines(conduit_gj)}
    def _nearest_on_feature(fid: int, p: List[float]) -> Tuple[float, float]:
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

    vaults = _iter_points(vaults_gj)
    v_by_feature: DefaultDict[int, List[Tuple[float,int,List[float],Dict[str,Any]]]] = defaultdict(list)
    for vi, vpt, vprops in vaults:
        best = (float("inf"), None, 0.0)  # (dist_ft, fid, along_ft)
        for fid in feature_nodes.keys():
            d_ft, along_ft = _nearest_on_feature(fid, vpt)
            if d_ft < best[0]:
                best = (d_ft, fid, along_ft)
        if best[1] is not None and best[0] <= VAULT_ON_LINE_FT:
            v_by_feature[best[1]].append((best[2], vi, vpt, vprops))

    # sort each list by along distance
    for fid, items in v_by_feature.items():
        items.sort(key=lambda t: t[0])

    return v_by_feature

def _is_tie_point(current_props: Dict[str, Any], neighbor_props: Dict[str, Any]) -> bool:
    """
    A neighbor is a 'Tie-Point' branch if EITHER Distribution Type OR Fiber Letter #1 differs.
    Section and Conduit Type are carried along for context but don't trigger branching alone.
    """
    return (current_props.get("Distribution Type") != neighbor_props.get("Distribution Type")) or \
           (current_props.get("Fiber Letter #1")    != neighbor_props.get("Fiber Letter #1"))


def _find_start_node_from_t3(vaults_gj: Dict[str, Any], feature_nodes: Dict[int, List[Tuple[int,int]]]) -> Optional[Tuple[int,int]]:
    """
    Try to find a 'T3' vault and snap it to the nearest node; else return None.
    """
    candidates = []
    for _, vpt, vprops in _iter_points(vaults_gj):
        blob = " ".join(str(vprops.get(k,"")) for k in ("Name","ID","Type","v_layer","Label"))
        if "t3" in blob.lower():
            candidates.append(vpt)
    if not candidates:
        return None
    # Pick the first T3 we find and snap to closest node
    all_nodes = []
    for nodes in feature_nodes.values():
        all_nodes.extend(nodes)
    # de-duplicate
    all_nodes = list(dict.fromkeys(all_nodes))
    # pick closest in feet
    def _centroid(node_key):  # reconstruct approximate lon/lat from key
        return [node_key[0]/1e7, node_key[1]/1e7]
    t3 = candidates[0]
    best = (float("inf"), None)
    for nk in all_nodes:
        d = distance_feet(t3, _centroid(nk))
        if d < best[0]:
            best = (d, nk)
    return best[1]


# ------------------------ traversal + drafting ------------------------

def _build_l_linestring(a: Sequence[float], b: Sequence[float]) -> List[List[float]]:
    ax, ay = a
    bx, by = b
    elbow = [bx, ay]  # lon-first then lat
    return [list(a), elbow, list(b)]

def _dfs_ordered_vaults(
    conduit_gj: Dict[str, Any],
    vaults_gj: Dict[str, Any]
) -> List[Tuple[int, List[float], Dict[str, Any]]]:
    """
    Depth-first traversal starting at T3 (if present). At each node, detour into
    Tie-Point branches (Distribution Type OR Fiber Letter #1 differs) *before*
    continuing along the current feature. Returns a list of vault tuples in visit order.
    """
    feature_nodes, node_to_features, feature_props, _ = _build_conduit_topology(conduit_gj)
    v_by_feature = _assign_vaults_to_features(conduit_gj, feature_nodes, vaults_gj)

    # Map feature -> its conduit coords (to know vertex ordering length)
    lines = {fid:(coords, props) for fid, coords, props in _iter_lines(conduit_gj)}

    # Choose a start node: T3 if possible, else an endpoint node of any feature
    start_node = _find_start_node_from_t3(vaults_gj, feature_nodes)
    if start_node is None:
        # pick any endpoint (node that belongs to exactly one feature or endpoint of a polyline)
        endpoints = []
        for fid, nodes in feature_nodes.items():
            if nodes:
                endpoints.extend([nodes[0], nodes[-1]])
        start_node = endpoints[0] if endpoints else None
    if start_node is None:
        # nothing to traverse
        ordered: List[Tuple[int, List[float], Dict[str, Any]]] = []
        # still return vaults in any stable order to avoid empty output
        for fid, lst in v_by_feature.items():
            for _, vi, vpt, vprops in lst:
                ordered.append((vi, vpt, vprops))
        return ordered

    # Find a starting feature touching that node (pick any)
    start_fid = next(iter(node_to_features.get(start_node, [])), None)
    if start_fid is None:
        # same fallback as above
        ordered: List[Tuple[int, List[float], Dict[str, Any]]] = []
        for fid, lst in v_by_feature.items():
            for _, vi, vpt, vprops in lst:
                ordered.append((vi, vpt, vprops))
        return ordered

    visited_fids: set[int] = set()
    result: List[Tuple[int, List[float], Dict[str, Any]]] = []

    def _emit_vaults_along(fid: int, forward: bool):
        """Append vaults assigned to this feature in geometric order."""
        lst = v_by_feature.get(fid, [])
        if not lst:
            return
        items = lst if forward else list(reversed(lst))
        for _, vi, vpt, vprops in items:
            result.append((vi, vpt, vprops))

    def _traverse(fid: int, entry_node: Tuple[int,int]):
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, _ = lines[fid]
        nodes = feature_nodes[fid]

        # direction: make nodes[0] == entry_node if possible
        forward = True
        if nodes and nodes[-1] == entry_node:
            nodes = list(reversed(nodes))
            coords = list(reversed(coords))
            forward = True
        # If entry_node not at either end (junction at interior vertex), keep forward=True and rely on vault alongs.

        # Walk vertex by vertex
        current_props = feature_props[fid]

        # Emit vaults encountered on this feature first (in this feature's direction),
        # but we will interleave DFS at each vertex prior to moving beyond it.
        # For simplicity, emit all now, since we only need NAP ordering by feature visit order.
        _emit_vaults_along(fid, forward=True)

        for nk in nodes:
            neighbor_fids = [nf for nf in node_to_features.get(nk, set()) if nf != fid and nf not in visited_fids]
            if not neighbor_fids:
                continue

            # Split neighbors into tie-point vs same-path
            tie_neighbors = []
            same_neighbors = []
            for nf in neighbor_fids:
                if _is_tie_point(current_props, feature_props[nf]):
                    tie_neighbors.append(nf)
                else:
                    same_neighbors.append(nf)

            # Detour into tie-points first
            for nf in tie_neighbors:
                _traverse(nf, nk)

            # After tie-points, you *may* also have same-path splits; traverse them next.
            for nf in same_neighbors:
                _traverse(nf, nk)

        # end traverse

    _traverse(start_fid, start_node)
    return result


def _propose_naps_and_drops(
    conduit_gj: Dict[str, Any],
    vaults_gj: Dict[str, Any],
    sl_gj: Dict[str, Any]
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Revised siting that respects path order and tie-point branching.
    We still:
      - Place a NAP at a vault IF it serves at least one SL within 250' (L).
      - Underground L-shaped drops; Aerial straight-only will be added later.
    """
    ordered_vaults = _dfs_ordered_vaults(conduit_gj, vaults_gj)
    sl_points = _iter_points(sl_gj)

    # Quick proximity index SLs→nearby vault (<=250' L)
    nap_features: List[Dict[str, Any]] = []
    drop_features: List[Dict[str, Any]] = []

    # Build list of candidate SLs per vault (recomputed per vault for clarity)
    nap_counter = 0

    for vi, v_xy, v_props in ordered_vaults:
        nearby: List[Tuple[float, int, List[float], Dict[str, Any]]] = []
        for si, s_xy, s_props in sl_points:
            dL = l_distance_feet(v_xy, s_xy)
            if dL <= MAX_L_DROP_FT:
                nearby.append((dL, si, s_xy, s_props))
        if not nearby:
            # Not every vault becomes a NAP; skip if no SLs
            continue

        nearby.sort(key=lambda x: x[0])
        chosen = nearby[:NAP_MAX_SL]

        nap_counter += 1
        nap_props = {
            "v_layer": "NAP (proposed)",
            "Nap Location": "Underground",  # conduit implies UG; aerial logic later
            "NAP # (proposed)": nap_counter,
            # carry some vault context if present
            "vault_id": _norm_val(v_props.get("ID")) or _norm_val(v_props.get("Name")),
        }
        nap_features.append({
            "type": "Feature",
            "properties": nap_props,
            "geometry": {"type": "Point", "coordinates": list(v_xy)}
        })

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
    I/O wrapper (unchanged signature) — pass the **raw** conduit file if available.
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
