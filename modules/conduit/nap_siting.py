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

def _vec_meters(p1: Sequence[float], p2: Sequence[float]) -> Tuple[float, float]:
    """Approx local-meter vector from p1->p2 using crude lat/lon scaling (good enough for angles)."""
    lon1, lat1 = p1
    lon2, lat2 = p2
    lat_mid = (lat1 + lat2) / 2.0
    # meters per degree
    m_per_deg_lat = 111132.92 - 559.82 * math.cos(math.radians(2 * lat_mid)) + 1.175 * math.cos(math.radians(4 * lat_mid)) - 0.0023 * math.cos(math.radians(6 * lat_mid))
    m_per_deg_lon = (111412.84 * math.cos(math.radians(lat_mid)) - 93.5 * math.cos(math.radians(3 * lat_mid)) + 0.118 * math.cos(math.radians(5 * lat_mid)))
    return ((lon2 - lon1) * m_per_deg_lon, (lat2 - lat1) * m_per_deg_lat)

def _angle_between(ax: float, ay: float, bx: float, by: float) -> float:
    """Unsigned angle (radians) between 2 vectors a and b in meters-space."""
    da = math.hypot(ax, ay)
    db = math.hypot(bx, by)
    if da == 0 or db == 0:
        return math.pi  # worst case
    dot = (ax * bx + ay * by) / (da * db)
    dot = max(-1.0, min(1.0, dot))
    return math.acos(dot)


def _is_tie_point(current_props: Dict[str, Any], neighbor_props: Dict[str, Any]) -> bool:
    """
    A neighbor is a 'Tie-Point' branch when it is in the SAME Section as the
    current feature, but represents a different cable:
      - Distribution Type differs OR
      - Fiber Letter #1 differs OR
      - Conduit Type differs

    We only treat OTHER Sections as out-of-scope for the current pass; those are
    traversed later when we move to their Section.
    """
    if (current_props.get("Section") or "").lower() != (neighbor_props.get("Section") or "").lower():
        return False  # different Section is not a tie-point detour now; will be walked later

    return (
        (current_props.get("Distribution Type") or "") != (neighbor_props.get("Distribution Type") or "")
        or (current_props.get("Fiber Letter #1") or "") != (neighbor_props.get("Fiber Letter #1") or "")
        or (current_props.get("Conduit Type") or "") != (neighbor_props.get("Conduit Type") or "")
    )


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
    Traversal that:
      • Walks Sections in order: Section 1..6, then unlabeled (None).
      • Within the active Section, picks a chain start (endpoint if available).
      • Continues along *same cable* first (same Section, Distribution Type, Conduit Type, Fiber Letter #1),
        choosing the next feature that keeps the heading as straight as possible (min turn angle).
      • At each junction, detours FIRST into tie-points (same Section but different cable),
        then resumes the main cable.
      • Emits vaults in the exact order we flow along features (not all-at-once).
    """
    # --- Build topology and feature meta ---
    feature_nodes, node_to_features, feature_props, _ = _build_conduit_topology(conduit_gj)
    v_by_feature = _assign_vaults_to_features(conduit_gj, feature_nodes, vaults_gj)
    lines = {fid:(coords, props) for fid, coords, props in _iter_lines(conduit_gj)}

    # quick helpers on nodes <-> lon/lat
    def _nk_to_xy(nk: Tuple[int,int]) -> List[float]:
        return [nk[0] / 1e7, nk[1] / 1e7]

    def _feature_along_vertices(fid: int) -> Tuple[List[List[float]], List[Tuple[int,int]], List[float]]:
        """Return (coords, nodes, along_ft at each vertex)."""
        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        along = [0.0]
        for i in range(1, len(coords)):
            seg = distance_feet(coords[i-1], coords[i])
            along.append(along[-1] + seg)
        return coords, nodes, along

    # group features by Section label for ordered processing
    def _sec_norm(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        s = s.strip()
        return s if s and s.lower() != "null" else None

    section_order = [f"Section {i}" for i in range(1,7)]
    # Build: section -> set(fid)
    by_section: Dict[Optional[str], List[int]] = {}
    for fid in feature_nodes.keys():
        sec = _sec_norm(feature_props[fid].get("Section"))
        by_section.setdefault(sec, []).append(fid)

    # ensure stable ordering inside lists
    for lst in by_section.values():
        lst.sort(key=lambda f: (str(feature_props[f].get("ID") or ""), f))

    visited_fids: set[int] = set()
    result: List[Tuple[int, List[float], Dict[str, Any]]] = []

    def _emit_vaults_until(fid: int, up_to_along_ft: float, v_state: Dict[int, int]):
        """Emit vaults on this fid while their along<=up_to_along_ft (direction aware via v_state pointer)."""
        lst = v_by_feature.get(fid, [])
        if not lst:
            return
        # Ensure list sorted by along
        lst = sorted(lst, key=lambda t: t[0])
        # cursor state
        i = v_state.get(fid, 0)
        while i < len(lst) and lst[i][0] <= up_to_along_ft + 1e-6:
            _, vi, vpt, vprops = lst[i]
            result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i

    def _emit_remaining_vaults(fid: int, v_state: Dict[int,int]):
        lst = v_by_feature.get(fid, [])
        if not lst:
            return
        lst = sorted(lst, key=lambda t: t[0])
        i = v_state.get(fid, 0)
        while i < len(lst):
            _, vi, vpt, vprops = lst[i]
            result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i

    def _same_cable(a: int, b: int) -> bool:
        pa = feature_props[a]; pb = feature_props[b]
        return (
            (pa.get("Section") or "").lower() == (pb.get("Section") or "").lower() and
            (pa.get("Distribution Type") or "") == (pb.get("Distribution Type") or "") and
            (pa.get("Conduit Type") or "") == (pb.get("Conduit Type") or "") and
            (pa.get("Fiber Letter #1") or "") == (pb.get("Fiber Letter #1") or "")
        )

    # choose starts inside a section: prefer endpoints (nodes that touch only one feature from this section)
    def _section_endpoints(section_key: Optional[str]) -> List[Tuple[int, Tuple[int,int]]]:
        fids = [f for f in by_section.get(section_key, []) if f not in visited_fids]
        endpoints: List[Tuple[int, Tuple[int,int]]] = []
        section_set = set(fids)
        for fid in fids:
            nodes = feature_nodes[fid]
            if not nodes:
                continue
            for nk in (nodes[0], nodes[-1]):
                # how many features of THIS section touch this node?
                touching = [nf for nf in node_to_features.get(nk, set()) if nf in section_set]
                if len(touching) == 1:
                    endpoints.append((fid, nk))
        # stable ordering: leftmost->rightmost then bottom->top
        endpoints.sort(key=lambda t: ( _nk_to_xy(t[1])[0], _nk_to_xy(t[1])[1] ))
        return endpoints

    # find a predecessor heading vector for angle comparisons
    def _last_segment_vec(coords: List[List[float]], forward: bool) -> Tuple[float,float]:
        if forward:
            p1, p2 = coords[-2], coords[-1]
        else:
            p1, p2 = coords[1], coords[0]
        return _vec_meters(p1, p2)

    def _pick_next_same_cable(fid: int, tail_node: Tuple[int,int]) -> Optional[int]:
        """Among neighbors at tail_node, choose SAME-CABLE feature that yields the straightest continuation."""
        candidates = [nf for nf in node_to_features.get(tail_node, set())
                      if nf != fid and nf not in visited_fids and _same_cable(fid, nf)]
        if not candidates:
            return None
        # current direction vector at the tail of this feature:
        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        forward = (nodes[-1] == tail_node)  # we are leaving from the end we just arrived at
        ax, ay = _last_segment_vec(coords, forward=True) if forward else _last_segment_vec(coords, forward=False)

        best = (float("inf"), None)
        tail_xy = _nk_to_xy(tail_node)
        for nf in sorted(candidates, key=lambda f: (str(feature_props[f].get("ID") or ""), f)):
            ncoords, _ = lines[nf]
            nnodes = feature_nodes[nf]
            # pick the first segment vector of nf leaving this node
            if nnodes[0] == tail_node and len(ncoords) >= 2:
                bx, by = _vec_meters(ncoords[0], ncoords[1])
            elif nnodes[-1] == tail_node and len(ncoords) >= 2:
                bx, by = _vec_meters(ncoords[-1], ncoords[-2])
            else:
                # node is interior; pick whichever adjacent vertex exists
                k = nnodes.index(tail_node)
                if 0 < k < len(nnodes)-1:
                    # prefer the forward part
                    bx, by = _vec_meters(_nk_to_xy(nnodes[k]), _nk_to_xy(nnodes[k+1]))
                else:
                    continue
            ang = _angle_between(ax, ay, bx, by)
            if ang < best[0]:
                best = (ang, nf)
        return best[1]

    def _walk_feature(fid: int, entry_node: Tuple[int,int], v_state: Dict[int,int]):
        """Walk along one feature, vertex by vertex, interleaving tie-point detours and vault emission."""
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        props = feature_props[fid]

        # align direction so nodes[0] == entry_node if possible
        if nodes and nodes[-1] == entry_node:
            nodes = list(reversed(nodes))
            coords = list(reversed(coords))

        # precompute along per-vertex, and vaults along
        coords0, nodes0, along = _feature_along_vertices(fid)

        # if reversed, recompute along on the reversed coords
        if nodes0[0] != nodes[0]:
            coords0 = list(reversed(coords0))
            nodes0  = list(reversed(nodes0))
            # recompute along for reversed orientation
            along = [0.0]
            for i in range(1, len(coords0)):
                along.append(along[-1] + distance_feet(coords0[i-1], coords0[i]))

        # map: node_key -> along_ft at that vertex
        along_at_node = {nodes0[i]: along[i] for i in range(len(nodes0))}

        # emit vaults as we pass each vertex, detouring into tie-points first
        for i, nk in enumerate(nodes):
            # 1) emit vaults up to this vertex
            _emit_vaults_until(fid, along_at_node[nk], v_state)

            # 2) handle tie-point branches at this junction (same Section, different cable)
            neighbor_fids = [nf for nf in node_to_features.get(nk, set()) if nf != fid and nf not in visited_fids]
            tie_neighbors = [nf for nf in neighbor_fids if _is_tie_point(props, feature_props[nf])]
            for nf in sorted(tie_neighbors, key=lambda f: (str(feature_props[f].get("ID") or ""), f)):
                _walk_feature(nf, nk, v_state)

        # at the tail, try to continue on SAME-CABLE first (straightest)
        tail = nodes[-1]
        next_same = _pick_next_same_cable(fid, tail)
        if next_same is not None:
            _walk_feature(next_same, tail, v_state)

        # finally, flush any remaining vaults on this feature
        _emit_remaining_vaults(fid, v_state)

    # --- master order over sections: 1..6, then unlabeled ---
    sections_to_run: List[Optional[str]] = [s for s in section_order if s in by_section] + ([None] if None in by_section else [])
    v_state: Dict[int,int] = {}  # per-feature vault cursor

    for sec in sections_to_run:
        # walk all components in this section
        while True:
            # pick an unvisited starting feature in this section
            candidates = [f for f in by_section.get(sec, []) if f not in visited_fids]
            if not candidates:
                break

            # prefer endpoints (true ends in this section); else fall back to any node
            endpoints = _section_endpoints(sec)
            start_fid, start_nk = None, None
            for (fid, nk) in endpoints:
                if fid in candidates:
                    start_fid, start_nk = fid, nk
                    break
            if start_fid is None:
                # no endpoints; pick a stable feature and start at its first node
                start_fid = sorted(candidates, key=lambda f: (str(feature_props[f].get("ID") or ""), f))[0]
                start_nk = feature_nodes[start_fid][0]

            _walk_feature(start_fid, start_nk, v_state)

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
