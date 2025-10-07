# modules/nap/traversal.py
from __future__ import annotations
from typing import Dict, Any, List, Tuple, Optional, DefaultDict, Sequence
import math
import re
from collections import defaultdict

from modules.geo.geometry import (
    distance_feet, l_distance_feet,
    point_segment_distance_feet
)
from .constants import NODE_SNAP_FT
from .topology import (
    iter_lines, iter_points, norm_val,
    build_conduit_topology, assign_vaults_to_features,
    vec_meters, find_start_node_from_t3
)

def _sec_norm(s: Optional[str]) -> Optional[str]:
    if s is None:
        return None
    s = s.strip()
    return s if s and s.lower() != "null" else None

def _dist_num(s: Optional[str]) -> int:
    if not s: return 0
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else 0

def dfs_ordered_vaults(
    conduit_gj: Dict[str, Any],
    vaults_gj: Dict[str, Any],
    t3_gj: Dict[str, Any]
) -> List[Tuple[int, List[float], Dict[str, Any]]]:
    """
    Depth-first traversal following your field rules (T3-anchored, Distribution priority,
    equal-24 'B' before 'A', straightest-turn continuation, and deferred same-section ties).
    Returns [(vault_index, [lon,lat], vault_props), ...] in the visit order.
    """
    # --- topology & assignments ---
    feature_nodes, node_to_features, feature_props, _ = build_conduit_topology(conduit_gj)
    v_by_feature = assign_vaults_to_features(conduit_gj, feature_nodes, vaults_gj)
    lines = {fid: (coords, props) for fid, coords, props in iter_lines(conduit_gj)}

    def nk_to_xy(nk: Tuple[int,int]) -> List[float]:
        return [nk[0] / 1e7, nk[1] / 1e7]

    def feature_along_vertices(fid: int) -> Tuple[List[List[float]], List[Tuple[int,int]], List[float]]:
        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        along = [0.0]
        for i in range(1, len(coords)):
            seg = distance_feet(coords[i-1], coords[i])
            along.append(along[-1] + seg)
        return coords, nodes, along

    section_order = [f"Section {i}" for i in range(1,7)]
    by_section: Dict[Optional[str], List[int]] = {}
    for fid in feature_nodes.keys():
        sec = _sec_norm(feature_props[fid].get("Section"))
        by_section.setdefault(sec, []).append(fid)
    for lst in by_section.values():
        lst.sort(key=lambda f: (str(feature_props[f].get("ID") or ""), f))

    visited_fids: set[int] = set()
    result: List[Tuple[int, List[float], Dict[str, Any]]] = []

    # de-dup NAPs by rounded XY
    seen_vault_xy: set[Tuple[int,int]] = set()
    def round_xy(pt: List[float], scale: int = 8) -> Tuple[int,int]:
        return (int(round(pt[0]*(10**scale))), int(round(pt[1]*(10**scale))))

    def emit_vaults_until(fid: int, up_to_along_ft: float, v_state: Dict[int,int]):
        lst = v_by_feature.get(fid, [])
        if not lst: return
        lst = sorted(lst, key=lambda t: t[0])
        i = v_state.get(fid, 0)
        while i < len(lst) and lst[i][0] <= up_to_along_ft + 1e-6:
            _, vi, vpt, vprops = lst[i]
            key = round_xy(vpt)
            if key not in seen_vault_xy:
                seen_vault_xy.add(key)
                result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i

    def emit_remaining_vaults(fid: int, v_state: Dict[int,int]):
        lst = v_by_feature.get(fid, [])
        if not lst: return
        lst = sorted(lst, key=lambda t: t[0])
        i = v_state.get(fid, 0)
        while i < len(lst):
            _, vi, vpt, vprops = lst[i]
            key = round_xy(vpt)
            if key not in seen_vault_xy:
                seen_vault_xy.add(key)
                result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i

    def same_cable(a: int, b: int) -> bool:
        pa = feature_props[a]; pb = feature_props[b]
        sec_eq  = (pa.get("Section") or "").strip().lower() == (pb.get("Section") or "").strip().lower()
        dist_eq = _dist_num(pa.get("Distribution Type")) == _dist_num(pb.get("Distribution Type"))
        let_a   = (pa.get("Fiber Letter #1") or "").strip().upper()
        let_b   = (pb.get("Fiber Letter #1") or "").strip().upper()
        return sec_eq and dist_eq and (let_a == let_b or not (let_a and let_b))

    def last_seg_vec(coords: List[List[float]], forward: bool) -> Tuple[float,float]:
        if forward:
            p1, p2 = coords[-2], coords[-1]
        else:
            p1, p2 = coords[1], coords[0]
        return vec_meters(p1, p2)

    def angle_between(ax: float, ay: float, bx: float, by: float) -> float:
        da = math.hypot(ax, ay) or 1e-9
        db = math.hypot(bx, by) or 1e-9
        cosv = max(-1.0, min(1.0, (ax*bx + ay*by) / (da*db)))
        return math.degrees(math.acos(cosv))

    def pick_next_same_cable(fid: int, tail_node: Tuple[int,int]) -> Optional[int]:
        candidates = [
            nf for nf in node_to_features.get(tail_node, set())
            if nf != fid and nf not in visited_fids and same_cable(fid, nf)
        ]
        if not candidates:
            return None
        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        forward = (nodes[-1] == tail_node)
        ax, ay = last_seg_vec(coords, forward=True) if forward else last_seg_vec(coords, forward=False)

        best = (float("inf"), None)
        for nf in sorted(candidates, key=lambda f: (str(feature_props[f].get("ID") or ""), f)):
            ncoords, _ = lines[nf]
            nnodes = feature_nodes[nf]
            if nnodes[0] == tail_node and len(ncoords) >= 2:
                bx, by = vec_meters(ncoords[0], ncoords[1])
            elif nnodes[-1] == tail_node and len(ncoords) >= 2:
                bx, by = vec_meters(ncoords[-1], ncoords[-2])
            else:
                k = nnodes.index(tail_node)
                if 0 < k < len(nnodes)-1:
                    bx, by = vec_meters([nnodes[k][0]/1e7, nnodes[k][1]/1e7], [nnodes[k+1][0]/1e7, nnodes[k+1][1]/1e7])
                else:
                    continue
            ang = angle_between(ax, ay, bx, by)
            if ang < best[0]:
                best = (ang, nf)
        return best[1]

    def is_tie_neighbor(base_fid: int, neighbor_fid: int) -> bool:
        pa = feature_props[base_fid]; pb = feature_props[neighbor_fid]
        same_section = _sec_norm(pa.get("Section")) == _sec_norm(pb.get("Section"))
        return same_section and not same_cable(base_fid, neighbor_fid)

    def letter_rank_equal24(letter: str) -> Tuple[int, str]:
        L = (letter or "").strip().upper()
        if L == "B": return (0, L)
        if L == "A": return (1, L)
        return (2, L)

    def choose_start_fid(touching: List[int]) -> int:
        def score(fid: int) -> Tuple[int, int, str, str, int]:
            p = feature_props[fid]
            dist = _dist_num(p.get("Distribution Type"))
            letter = (p.get("Fiber Letter #1") or "").strip().upper()
            return (-dist, 0 if letter == "A" else 1, letter, str(p.get("ID") or ""), fid)
        return sorted(touching, key=score)[0]

    def section_endpoints(section_key: Optional[str]) -> List[Tuple[int, Tuple[int,int]]]:
        fids = [f for f in by_section.get(section_key, []) if f not in visited_fids]
        endpoints: List[Tuple[int, Tuple[int,int]]] = []
        section_set = set(fids)
        for fid in fids:
            nodes = feature_nodes[fid]
            if not nodes:
                continue
            for nk in (nodes[0], nodes[-1]):
                touching = [nf for nf in (node_to_features.get(nk, set()) or []) if nf in section_set]
                if len(touching) == 1:
                    endpoints.append((fid, nk))
        endpoints.sort(key=lambda t: (nk_to_xy(t[1])[0], nk_to_xy(t[1])[1]))
        return endpoints

    def neighbors_at_node(fid: int, nk: Tuple[int,int]) -> List[int]:
        # shared-vertex neighbors + snapped mid-segment neighbors within NODE_SNAP_FT
        xy = nk_to_xy(nk)
        pa = feature_props[fid]
        same_section_fids = [nf for nf in by_section.get(_sec_norm(pa.get("Section")), []) if nf != fid]

        def min_dist_to_feature(xy_pt: List[float], other_fid: int) -> float:
            ocoords, _ = lines[other_fid]
            best = float("inf")
            for i in range(1, len(ocoords)):
                d = point_segment_distance_feet(xy_pt, ocoords[i-1], ocoords[i])
                if d < best:
                    best = d
            return best

        shared = {nf for nf in (node_to_features.get(nk, set()) or [])
                  if nf != fid and nf not in visited_fids and is_tie_neighbor(fid, nf)}

        snapped = set()
        for nf in same_section_fids:
            if nf in visited_fids or same_cable(fid, nf):
                continue
            if min_dist_to_feature(xy, nf) <= NODE_SNAP_FT:
                snapped.add(nf)

        return list(shared | snapped)

    def walk_feature(fid: int, entry_node: Tuple[int,int], v_state: Dict[int,int]):
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        props = feature_props[fid]
        cur_dist = _dist_num(props.get("Distribution Type"))

        if nodes and nodes[-1] == entry_node:
            nodes = list(reversed(nodes))
            coords = list(reversed(coords))

        coords0, nodes0, along = feature_along_vertices(fid)
        if nodes0[0] != nodes[0]:
            coords0  = list(reversed(coords0))
            nodes0   = list(reversed(nodes0))
            along    = [0.0]
            for i in range(1, len(coords0)):
                along.append(along[-1] + distance_feet(coords0[i-1], coords0[i]))
        along_at_node = {nodes0[i]: along[i] for i in range(len(nodes0))}

        for nk in nodes:
            # 1) emit vaults up to this vertex
            emit_vaults_until(fid, along_at_node[nk], v_state)

            # 2) neighbors and priority buckets
            neighbor_fids = neighbors_at_node(fid, nk)
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

            # 3) lower detours (smallest fiber count first)
            for nf in sorted(lowers, key=lambda f: (_dist_num(feature_props[f].get("Distribution Type")),
                                                    str(feature_props[f].get("ID") or ""), f)):
                walk_feature(nf, nk, v_state)

            # 4) equal-24 detours (B -> A -> others)
            for nf in sorted(equal24, key=lambda f: (letter_rank_equal24(feature_props[f].get("Fiber Letter #1")),
                                                     str(feature_props[f].get("ID") or ""), f)):
                walk_feature(nf, nk, v_state)

            # equal_non24/higher are deferred

        # 5) continue along same cable by straightest angle
        tail = nodes[-1]
        next_same = pick_next_same_cable(fid, tail)
        if next_same is not None:
            walk_feature(next_same, tail, v_state)

        # 6) after mainline: process deferred equals (non-24) then higher, from endpoints
        for nk in (nodes[0], nodes[-1]):
            neighbor_fids = neighbors_at_node(fid, nk)
            en24 = [nf for nf in neighbor_fids
                    if _dist_num(feature_props[nf].get("Distribution Type")) == cur_dist != 24]
            hi   = [nf for nf in neighbor_fids
                    if _dist_num(feature_props[nf].get("Distribution Type"))  > cur_dist]

            for nf in sorted(en24, key=lambda f: (str(feature_props[f].get("Fiber Letter #1") or ""),
                                                  str(feature_props[f].get("ID") or ""), f)):
                walk_feature(nf, nk, v_state)
            for nf in sorted(hi, key=lambda f: (-_dist_num(feature_props[f].get("Distribution Type")),
                                                str(feature_props[f].get("Fiber Letter #1") or ""),
                                                str(feature_props[f].get("ID") or ""), f)):
                walk_feature(nf, nk, v_state)

        # 7) trailing vaults on this feature
        emit_remaining_vaults(fid, v_state)

    # master Section order
    sections_to_run: List[Optional[str]] = [s for s in section_order if s in by_section] + ([None] if None in by_section else [])
    v_state: Dict[int,int] = {}

    for sec in sections_to_run:
        section_fids = set(f for f in by_section.get(sec, []) if f not in visited_fids)
        if not section_fids:
            continue

        # T3-anchored start (highest Distribution; tie -> A first)
        start_nk = find_start_node_from_t3(t3_gj, feature_nodes)
        start_fid = None
        if start_nk is not None:
            touching = [nf for nf in (node_to_features.get(start_nk, set()) or []) if nf in section_fids]
            if touching:
                def score(fid: int) -> Tuple[int, int, str, str, int]:
                    p = feature_props[fid]
                    dist = _dist_num(p.get("Distribution Type"))
                    letter = (p.get("Fiber Letter #1") or "").strip().upper()
                    return (-dist, 0 if letter == "A" else 1, letter, str(p.get("ID") or ""), fid)
                start_fid = sorted(touching, key=score)[0]

        while True:
            candidates = [f for f in by_section.get(sec, []) if f not in visited_fids]
            if not candidates:
                break
            if start_fid is None:
                # pick an endpoint if possible; else smallest ID
                # (We recompute endpoints within this component)
                def section_endpoints_local() -> List[Tuple[int, Tuple[int,int]]]:
                    fids = [f for f in by_section.get(sec, []) if f not in visited_fids]
                    endpoints: List[Tuple[int, Tuple[int,int]]] = []
                    section_set = set(fids)
                    for fid in fids:
                        nodes = feature_nodes[fid]
                        if not nodes: continue
                        for nk in (nodes[0], nodes[-1]):
                            touching = [nf for nf in (node_to_features.get(nk, set()) or []) if nf in section_set]
                            if len(touching) == 1:
                                endpoints.append((fid, nk))
                    endpoints.sort(key=lambda t: (nk_to_xy(t[1])[0], nk_to_xy(t[1])[1]))
                    return endpoints

                eps = section_endpoints_local()
                if eps:
                    start_fid, start_nk = eps[0]
                else:
                    fid = sorted(candidates, key=lambda f: (str(feature_props[f].get("ID") or ""), f))[0]
                    start_nk = feature_nodes[fid][0]
                    start_fid = fid

            walk_feature(start_fid, start_nk, v_state)
            start_fid = None

    return result
