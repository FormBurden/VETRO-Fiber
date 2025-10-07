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


def _find_start_node_from_t3(
    t3_gj: Dict[str, Any],
    feature_nodes: Dict[int, List[Tuple[int,int]]]
) -> Optional[Tuple[int,int]]:
    """
    Use a dedicated T3 GeoJSON (Point feature(s)) to pick the start node.
    We snap the FIRST T3 point we find to the nearest network node.

    If no valid T3 point exists, return None (the caller will fall back).
    """
    # collect all T3 points
    t3_pts: List[List[float]] = []
    for f in t3_gj.get("features", []):
        g = f.get("geometry") or {}
        if g.get("type") == "Point" and isinstance(g.get("coordinates"), (list, tuple)) and len(g["coordinates"]) == 2:
            t3_pts.append(g["coordinates"])

    if not t3_pts:
        return None

    # all nodes in network
    all_nodes = []
    for nodes in feature_nodes.values():
        all_nodes.extend(nodes)
    if not all_nodes:
        return None
    # de-dup
    all_nodes = list(dict.fromkeys(all_nodes))

    # centroid helper
    def _centroid(nk: Tuple[int,int]) -> List[float]:
        return [nk[0] / 1e7, nk[1] / 1e7]

    # “You can pick any one” — just use the first T3 point
    t3 = t3_pts[0]
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
    vaults_gj: Dict[str, Any],
    t3_gj: Dict[str, Any]
) -> List[Tuple[int, List[float], Dict[str, Any]]]:
    """
    Traversal that implements the user's field rules:

    Start at T3:
      • Pick the touching conduit with the HIGHEST Distribution Type (72 > 48 > 24).
      • If tied, prefer Fiber Letter #1 = 'A' first, then alphabetical.

    At each visited node (depth-first):
      1) Emit any NAP located at this node (junction first).
      2) Detour into LOWER Distribution Type branches first (e.g., 48 -> 24).
         • If multiple lowers exist, visit smaller fiber-count first (e.g., 24 before 48).
      3) If CURRENT conduit is 24ct, also detour equal-24 branches now with Fiber Letter order:
            'B' first, then 'A', then others alphabetically.
      4) Continue along SAME CABLE (same Section, Distribution Type, Conduit Type, Fiber Letter #1)
         using the straightest-turn rule (min angle).
      5) Finally, process any remaining same-Section ties (equal/higher) that were deferred.
    """
    import re

    # --- topology & assignments ---
    feature_nodes, node_to_features, feature_props, _ = _build_conduit_topology(conduit_gj)
    v_by_feature = _assign_vaults_to_features(conduit_gj, feature_nodes, vaults_gj)
    lines = {fid: (coords, props) for fid, coords, props in _iter_lines(conduit_gj)}

    def _nk_to_xy(nk: Tuple[int,int]) -> List[float]:
        return [nk[0] / 1e7, nk[1] / 1e7]

    def _feature_along_vertices(fid: int) -> Tuple[List[List[float]], List[Tuple[int,int]], List[float]]:
        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        along = [0.0]
        for i in range(1, len(coords)):
            seg = distance_feet(coords[i-1], coords[i])
            along.append(along[-1] + seg)
        return coords, nodes, along

    def _sec_norm(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        s = s.strip()
        return s if s and s.lower() != "null" else None

    def _dist_num(s: Optional[str]) -> int:
        if not s: return 0
        m = re.search(r"(\d+)", s)
        return int(m.group(1)) if m else 0

    section_order = [f"Section {i}" for i in range(1,7)]
    by_section: Dict[Optional[str], List[int]] = {}
    for fid in feature_nodes.keys():
        sec = _sec_norm(feature_props[fid].get("Section"))
        by_section.setdefault(sec, []).append(fid)
    for lst in by_section.values():
        lst.sort(key=lambda f: (str(feature_props[f].get("ID") or ""), f))

    visited_fids: set[int] = set()
    result: List[Tuple[int, List[float], Dict[str, Any]]] = []

    # --- shared-node NAP de-dup (by rounded XY) ---
    seen_vault_xy: set[Tuple[int,int]] = set()
    def _round_xy(pt: List[float], scale: int = 8) -> Tuple[int,int]:
        return (int(round(pt[0]*(10**scale))), int(round(pt[1]*(10**scale))))

    # --- emit helpers ---
    def _emit_vaults_until(fid: int, up_to_along_ft: float, v_state: Dict[int,int]):
        lst = v_by_feature.get(fid, [])
        if not lst:
            return
        lst = sorted(lst, key=lambda t: t[0])  # (along_ft, nap_idx, [x,y], props)
        i = v_state.get(fid, 0)
        while i < len(lst) and lst[i][0] <= up_to_along_ft + 1e-6:
            _, vi, vpt, vprops = lst[i]
            key = _round_xy(vpt)
            if key not in seen_vault_xy:
                seen_vault_xy.add(key)
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
            key = _round_xy(vpt)
            if key not in seen_vault_xy:
                seen_vault_xy.add(key)
                result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i

    # --- identity + angles ---
    def _same_cable(a: int, b: int) -> bool:
        pa = feature_props[a]; pb = feature_props[b]
        sec_eq  = (pa.get("Section") or "").strip().lower() == (pb.get("Section") or "").strip().lower()
        # numeric compare (handles '24ct', '48 ct', etc.)
        dist_eq = _dist_num(pa.get("Distribution Type")) == _dist_num(pb.get("Distribution Type"))
        let_a   = (pa.get("Fiber Letter #1") or "").strip().upper()
        let_b   = (pb.get("Fiber Letter #1") or "").strip().upper()
        # SAME line if Section + Distribution match.
        # If letters exist for both and match, great; if letters missing on either side, still same line.
        # (We deliberately IGNORE Conduit Type differences like Road Shot vs 1x1.25".)
        return sec_eq and dist_eq and (let_a == let_b or not (let_a and let_b))

    def _last_seg_vec(coords: List[List[float]], forward: bool) -> Tuple[float,float]:
        if forward:
            p1, p2 = coords[-2], coords[-1]
        else:
            p1, p2 = coords[1], coords[0]
        return _vec_meters(p1, p2)

    def _angle_between(ax: float, ay: float, bx: float, by: float) -> float:
        da = math.hypot(ax, ay) or 1e-9
        db = math.hypot(bx, by) or 1e-9
        cosv = max(-1.0, min(1.0, (ax*bx + ay*by) / (da*db)))
        return math.degrees(math.acos(cosv))

    def _pick_next_same_cable(fid: int, tail_node: Tuple[int,int]) -> Optional[int]:
        candidates = [
            nf for nf in node_to_features.get(tail_node, set())
            if nf != fid and nf not in visited_fids and _same_cable(fid, nf)
        ]
        if not candidates:
            return None
        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        forward = (nodes[-1] == tail_node)
        ax, ay = _last_seg_vec(coords, forward=True) if forward else _last_seg_vec(coords, forward=False)

        best = (float("inf"), None)
        for nf in sorted(candidates, key=lambda f: (str(feature_props[f].get("ID") or ""), f)):
            ncoords, _ = lines[nf]
            nnodes = feature_nodes[nf]
            if nnodes[0] == tail_node and len(ncoords) >= 2:
                bx, by = _vec_meters(ncoords[0], ncoords[1])
            elif nnodes[-1] == tail_node and len(ncoords) >= 2:
                bx, by = _vec_meters(ncoords[-1], ncoords[-2])
            else:
                k = nnodes.index(tail_node)
                if 0 < k < len(nnodes)-1:
                    bx, by = _vec_meters(_nk_to_xy(nnodes[k]), _nk_to_xy(nnodes[k+1]))
                else:
                    continue
            ang = _angle_between(ax, ay, bx, by)
            if ang < best[0]:
                best = (ang, nf)
        return best[1]

    # Is neighbor a same-Section tie (not same cable) we might detour into?
    def _is_tie_neighbor(base_fid: int, neighbor_fid: int) -> bool:
        pa = feature_props[base_fid]; pb = feature_props[neighbor_fid]
        same_section = _sec_norm(pa.get("Section")) == _sec_norm(pb.get("Section"))
        return same_section and not _same_cable(base_fid, neighbor_fid)

    # Start choice at T3: highest Distribution Type; tie -> A first, then alphabetical, then ID
    def _choose_start_fid(touching: List[int]) -> int:
        def score(fid: int) -> Tuple[int, int, str, str, int]:
            p = feature_props[fid]
            dist = _dist_num(p.get("Distribution Type"))
            letter = (p.get("Fiber Letter #1") or "").strip().upper()
            return (-dist, 0 if letter == "A" else 1, letter, str(p.get("ID") or ""), fid)
        return sorted(touching, key=score)[0]

    # Sort helper for equal-24 detours: B -> A -> others (alphabetical)
    def _letter_rank_for_equal24(letter: str) -> Tuple[int, str]:
        L = (letter or "").strip().upper()
        if L == "B": return (0, L)
        if L == "A": return (1, L)
        return (2, L)  # then C, D, ...

    def _section_endpoints(section_key: Optional[str]) -> List[Tuple[int, Tuple[int,int]]]:
        fids = [f for f in by_section.get(section_key, []) if f not in visited_fids]
        endpoints: List[Tuple[int, Tuple[int,int]]] = []
        section_set = set(fids)
        for fid in fids:
            nodes = feature_nodes[fid]
            if not nodes:
                continue
            for nk in (nodes[0], nodes[-1]):
                touching = [nf for nf in node_to_features.get(nk, set()) if nf in section_set]
                if len(touching) == 1:
                    endpoints.append((fid, nk))
        endpoints.sort(key=lambda t: (_nk_to_xy(t[1])[0], _nk_to_xy(t[1])[1]))
        return endpoints

    def _walk_feature(fid: int, entry_node: Tuple[int,int], v_state: Dict[int,int]):
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        props = feature_props[fid]
        cur_dist = _dist_num(props.get("Distribution Type"))

        # orient from entry_node to the other end
        if nodes and nodes[-1] == entry_node:
            nodes = list(reversed(nodes))
            coords = list(reversed(coords))

        # along-from-start per vertex (for vault emission)
        coords0, nodes0, along = _feature_along_vertices(fid)
        if nodes0[0] != nodes[0]:
            coords0  = list(reversed(coords0))
            nodes0   = list(reversed(nodes0))
            along    = [0.0]
            for i in range(1, len(coords0)):
                along.append(along[-1] + distance_feet(coords0[i-1], coords0[i]))
        along_at_node = {nodes0[i]: along[i] for i in range(len(nodes0))}

        # helper: distance from an XY to another feature's polyline (feet)
        def _min_dist_to_feature(xy: List[float], other_fid: int) -> float:
            ocoords, _ = lines[other_fid]
            best = float("inf")
            for i in range(1, len(ocoords)):
                d = point_segment_distance_feet(xy, ocoords[i-1], ocoords[i])
                if d < best:
                    best = d
            return best

        # helper: neighbors at node nk including “snapped” mid-segment tees within NODE_SNAP_FT
        def _neighbors_at_node(nk: Tuple[int,int]) -> List[int]:
            xy = _nk_to_xy(nk)

            # 1) True shared-vertex neighbors
            shared = {nf for nf in node_to_features.get(nk, set())
                    if nf != fid and nf not in visited_fids and _is_tie_neighbor(fid, nf)}

            # 2) Snap-to-polyline neighbors (mid-segment tees in the SAME Section, not same cable)
            same_section_fids = [nf for nf in by_section.get(_sec_norm(props.get("Section")), []) if nf != fid]
            snapped = set()
            for nf in same_section_fids:
                if nf in visited_fids or _same_cable(fid, nf):
                    continue
                if _min_dist_to_feature(xy, nf) <= NODE_SNAP_FT:
                    snapped.add(nf)

            return list(shared | snapped)

        # traverse each vertex of this feature
        for nk in nodes:
            # 1) Emit NAPs located up to this node first (junction numbered before detours)
            _emit_vaults_until(fid, along_at_node[nk], v_state)

            # 2) Build neighbor lists with Distribution-based priorities
            neighbor_fids = _neighbors_at_node(nk)

            lowers, equal24, equal_non24, higher = [], [], [], []
            for nf in neighbor_fids:
                ndist = _dist_num(feature_props[nf].get("Distribution Type"))
                if ndist < cur_dist:
                    lowers.append(nf)
                elif ndist == cur_dist == 24:
                    equal24.append(nf)
                elif ndist == cur_dist:
                    equal_non24.append(nf)
                else:
                    higher.append(nf)

            # 3) Detour into LOWER first (smallest fiber count first)
            for nf in sorted(
                lowers,
                key=lambda f: (_dist_num(feature_props[f].get("Distribution Type")),
                            str(feature_props[f].get("ID") or ""), f)
            ):
                _walk_feature(nf, nk, v_state)

            # 4) If CURRENT is 24ct, then detour equal-24 now (B -> A -> others)
            for nf in sorted(
                equal24,
                key=lambda f: _letter_rank_for_equal24(feature_props[f].get("Fiber Letter #1"))
                            + (str(feature_props[f].get("ID") or ""),)
            ):
                _walk_feature(nf, nk, v_state)

            # (We purposely do NOT walk equal_non24/higher at this point; they’re deferred.)

        # 5) Continue along SAME CABLE (straightest turn from the tail node)
        tail = nodes[-1]
        next_same = _pick_next_same_cable(fid, tail)
        if next_same is not None:
            _walk_feature(next_same, tail, v_state)

    # 6) After continuing the mainline, process any remaining same-Section ties we deferred
    #    (equal_non24 first, then higher) — include snapped mid-segment neighbors at endpoints.
    for nk in (nodes[0], nodes[-1]):
        neighbor_fids = _neighbors_at_node(nk)

        en24 = [nf for nf in neighbor_fids
                if _dist_num(feature_props[nf].get("Distribution Type")) == cur_dist != 24]
        hi   = [nf for nf in neighbor_fids
                if _dist_num(feature_props[nf].get("Distribution Type"))  > cur_dist]

        for nf in sorted(
            en24,
            key=lambda f: (str(feature_props[f].get("Fiber Letter #1") or ""),
                           str(feature_props[f].get("ID") or ""), f)
        ):
            _walk_feature(nf, nk, v_state)

        for nf in sorted(
            hi,
            key=lambda f: (-_dist_num(feature_props[f].get("Distribution Type")),
                           str(feature_props[f].get("Fiber Letter #1") or ""),
                           str(feature_props[f].get("ID") or ""), f)
        ):
            _walk_feature(nf, nk, v_state)

    # 7) trailing NAPs on this feature
    _emit_remaining_vaults(fid, v_state)


    # master Section order
    sections_to_run: List[Optional[str]] = [s for s in section_order if s in by_section] + ([None] if None in by_section else [])
    v_state: Dict[int,int] = {}

    for sec in sections_to_run:
        section_fids = set(f for f in by_section.get(sec, []) if f not in visited_fids)
        if not section_fids:
            continue

        # T3-anchored start (highest Distribution Type; tie -> A first)
        start_nk = _find_start_node_from_t3(t3_gj, feature_nodes)
        start_fid = None
        if start_nk is not None:
            touching = [nf for nf in (node_to_features.get(start_nk, set()) or []) if nf in section_fids]
            if touching:
                start_fid = _choose_start_fid(touching)

        # Walk all connected components in this Section
        while True:
            candidates = [f for f in by_section.get(sec, []) if f not in visited_fids]
            if not candidates:
                break
            if start_fid is None:
                endpoints = _section_endpoints(sec)
                pick = None
                for (fid, nk) in endpoints:
                    if fid in candidates:
                        pick = (fid, nk)
                        break
                if pick is None:
                    fid = sorted(candidates, key=lambda f: (str(feature_props[f].get("ID") or ""), f))[0]
                    nk = feature_nodes[fid][0]
                    pick = (fid, nk)
                start_fid, start_nk = pick

            _walk_feature(start_fid, start_nk, v_state)
            start_fid = None

    return result



def _propose_naps_and_drops(
    conduit_gj: Dict[str, Any],
    vaults_gj: Dict[str, Any],
    sl_gj: Dict[str, Any],
    t3_gj: Dict[str, Any]
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Siting that respects path order and tie-point branching, starting from T3 anchor.
    Underground: L-shaped drops; (Aerial straight-only will be added later).
    """
    ordered_vaults = _dfs_ordered_vaults(conduit_gj, vaults_gj, t3_gj)
    sl_points = _iter_points(sl_gj)

    nap_features: List[Dict[str, Any]] = []
    drop_features: List[Dict[str, Any]] = []

    nap_counter = 0
    for vi, v_xy, v_props in ordered_vaults:
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
            "NAP #": nap_counter,  # <— renamed
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
                "NAP #": nap_counter,  # <— renamed
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
    t3_vault_filename: str,
    nap_out_name: str = "nap_proposed.geojson",
    drop_out_name: str = "fiber-drop_proposed.geojson"
) -> Tuple[str, str]:
    """
    I/O wrapper — now requires a T3 vault file to anchor traversal.
    All paths are relative to modules.config.DATA_DIR.
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

    nap_fc, drop_fc = _propose_naps_and_drops(conduit_gj, vaults_gj, sl_gj, t3_gj)

    nap_path = os.path.join(base, nap_out_name)
    drop_path = os.path.join(base, drop_out_name)
    with open(nap_path, "w", encoding="utf-8") as f:
        json.dump(nap_fc, f, indent=2)
    with open(drop_path, "w", encoding="utf-8") as f:
        json.dump(drop_fc, f, indent=2)

    return nap_path, drop_path
