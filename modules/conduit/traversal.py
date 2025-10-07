# modules/conduit/traversal.py
from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from modules.geo.geometry import distance_feet, point_segment_distance_feet
from modules.conduit.topology import (
    build_conduit_topology,
    assign_vaults_to_features,
    iter_lines,
)

# Node snap tolerance (feet) for tee detection near a node
NODE_SNAP_FT = 20.0


def _vec_meters(p1: Sequence[float], p2: Sequence[float]) -> Tuple[float, float]:
    """Approx local-meter vector from p1->p2 using crude lat/lon scaling; good enough for angles."""
    lon1, lat1 = p1
    lon2, lat2 = p2
    lat_mid = (lat1 + lat2) / 2.0
    m_per_deg_lat = (
        111132.92
        - 559.82 * math.cos(math.radians(2 * lat_mid))
        + 1.175 * math.cos(math.radians(4 * lat_mid))
        - 0.0023 * math.cos(math.radians(6 * lat_mid))
    )
    m_per_deg_lon = (
        111412.84 * math.cos(math.radians(lat_mid))
        - 93.5 * math.cos(math.radians(3 * lat_mid))
        + 0.118 * math.cos(math.radians(5 * lat_mid))
    )
    return ((lon2 - lon1) * m_per_deg_lon, (lat2 - lat1) * m_per_deg_lat)


def _angle_between(ax: float, ay: float, bx: float, by: float) -> float:
    """Unsigned angle (degrees) between two meter-space vectors."""
    da = math.hypot(ax, ay) or 1e-9
    db = math.hypot(bx, by) or 1e-9
    cosv = max(-1.0, min(1.0, (ax * bx + ay * by) / (da * db)))
    return math.degrees(math.acos(cosv))


def _sec_norm(s: Optional[str]) -> Optional[str]:
    if s is None:
        return None
    s = s.strip()
    return s if s and s.lower() != "null" else None


def _dist_num(s: Optional[str]) -> int:
    if not s:
        return 0
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else 0


def dfs_ordered_vaults(
    conduit_gj: Dict[str, Any],
    vaults_gj: Dict[str, Any],
    t3_gj: Dict[str, Any],
) -> List[Tuple[int, List[float], Dict[str, Any]]]:
    """
    Depth-first traversal of the conduit network producing vaults (NAPs) in the field-rule order.

    Start at T3:
      • Pick touching conduit with HIGHEST Distribution Type (72 > 48 > 24).
      • If tied, prefer Fiber Letter #1 = 'A' first, then alphabetical.

    At each visited node (depth-first):
      1) Emit any NAP located at this node (junction first).
      2) Detour into LOWER Distribution branches first (smaller fiber count first).
      3) If CURRENT is 24ct, also detour equal-24 now in order B → A → others.
      4) Continue along SAME CABLE (same Section/Distribution/ConduitType/Letter#1) using straightest turn.
      5) Finally, process remaining same-Section ties (equal_non24, then higher).
    """
    feature_nodes, node_to_features, feature_props = build_conduit_topology(conduit_gj)
    v_by_feature = assign_vaults_to_features(conduit_gj, feature_nodes, vaults_gj)
    lines = {fid: (coords, props) for fid, coords, props in iter_lines(conduit_gj)}

    def _nk_to_xy(nk: Tuple[int, int]) -> List[float]:
        return [nk[0] / 1e7, nk[1] / 1e7]

    def _feature_along_vertices(fid: int) -> Tuple[List[List[float]], List[Tuple[int, int]], List[float]]:
        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        along = [0.0]
        for i in range(1, len(coords)):
            seg = distance_feet(coords[i - 1], coords[i])
            along.append(along[-1] + seg)
        return coords, nodes, along

    def _neighbor_rank(
        current_fid: int,
        prev_node_xy: Optional[List[float]],
        at_node_xy: List[float],
        candidate_fid: int,
        is_start: bool
    ) -> tuple:
        """
        Sort key for neighbor selection.

        Category order (lower sorts first):
          0) START PICK ONLY: highest Distribution first (reverse)  [applies when is_start=True]
          1) Lower Distribution detours (smaller fiber count first)
          2) (current is 24ct) equal-24 detours in order B -> A -> others
          3) Continue along SAME CABLE using straightest turn (smaller angle first)
          4) Remaining same-Section ties: equal_non24 first, then higher
          5) Everything else

        Angle is the unsigned angle between incoming and outgoing vectors (smaller is straighter).
        """
        pa = feature_props[current_fid]
        pb = feature_props[candidate_fid]

        dist_a = _dist_num(pa.get("Distribution Type"))
        dist_b = _dist_num(pb.get("Distribution Type"))

        sec_a = (pa.get("Section") or "").strip().lower()
        sec_b = (pb.get("Section") or "").strip().lower()

        let_b = (pb.get("Fiber Letter #1") or "").strip().upper()

        # default angle if we can't compute
        angle = 180.0
        if prev_node_xy is not None:
            coords_b, nodes_b, _ = _feature_along_vertices(candidate_fid)
            if coords_b and (coords_b[0] == at_node_xy or coords_b[-1] == at_node_xy):
                next_xy = coords_b[1] if coords_b[0] == at_node_xy and len(coords_b) > 1 else (
                          coords_b[-2] if coords_b[-1] == at_node_xy and len(coords_b) > 1 else at_node_xy)
            else:
                if coords_b:
                    d0 = distance_feet(coords_b[0], at_node_xy)
                    d1 = distance_feet(coords_b[-1], at_node_xy)
                    next_xy = coords_b[0] if d0 < d1 else coords_b[-1]
                else:
                    next_xy = at_node_xy

            ax, ay = _vec_meters(prev_node_xy, at_node_xy)
            bx, by = _vec_meters(at_node_xy, next_xy)
            angle = _angle_between(ax, ay, bx, by)

        if is_start:
            let_order = 0 if (pb.get("Fiber Letter #1") or "").strip().upper() == "A" else 1
            return (0, -dist_b, let_order, angle)

        if dist_b < dist_a:
            return (1, dist_b, angle)

        if dist_a == 24 and dist_b == 24:
            letter_rank = 2
            if let_b == "B":
                letter_rank = 0
            elif let_b == "A":
                letter_rank = 1
            return (2, letter_rank, angle)

        if _same_cable(current_fid, candidate_fid):
            return (3, angle)

        if sec_a == sec_b:
            if dist_b == dist_a and dist_a != 24:
                let_order = 0 if let_b == "A" else 1
                return (4, 0, let_order, angle)
            if dist_b > dist_a:
                return (4, 1, dist_b, angle)

        return (5, dist_b, angle)

    # ----- helpers to anchor at T3 -----
    def _find_start_node_from_t3(
        t3_gj: Dict[str, Any],
        feature_nodes: Dict[int, List[Tuple[int, int]]],
    ) -> Optional[Tuple[int, int]]:
        """Snap the first T3 point to the nearest network node."""
        t3_pts: List[List[float]] = []
        for f in t3_gj.get("features", []):
            g = f.get("geometry") or {}
            if g.get("type") == "Point" and isinstance(g.get("coordinates"), (list, tuple)) and len(g["coordinates"]) == 2:
                t3_pts.append(g["coordinates"])
        if not t3_pts:
            return None

        all_nodes: List[Tuple[int, int]] = []
        for nodes in feature_nodes.values():
            all_nodes.extend(nodes)
        if not all_nodes:
            return None

        all_nodes = list(dict.fromkeys(all_nodes))  # de-dup
        t3 = t3_pts[0]
        best = (float("inf"), None)
        for nk in all_nodes:
            d = distance_feet(t3, _nk_to_xy(nk))
            if d < best[0]:
                best = (d, nk)
        return best[1]

    def _choose_start_fid(touching: List[int]) -> int:
        # Highest Distribution first; tie → 'A' first, then alphabetical by letter, then ID for stability
        def score(fid: int) -> Tuple[int, int, str, str, int]:
            p = feature_props[fid]
            dist = _dist_num(p.get("Distribution Type"))
            letter = (p.get("Fiber Letter #1") or "").strip().upper()
            return (-dist, 0 if letter == "A" else 1, letter, str(p.get("ID") or ""), fid)
        return sorted(touching, key=score)[0]

    # ----- section buckets (stable within) -----
    section_order = [f"Section {i}" for i in range(1, 7)]
    by_section: Dict[Optional[str], List[int]] = {}
    def _sec_norm_local(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        s = s.strip()
        return s if s and s.lower() != "null" else None
    for fid in feature_nodes.keys():
        sec = _sec_norm_local(feature_props[fid].get("Section"))
        by_section.setdefault(sec, []).append(fid)
    for lst in by_section.values():
        lst.sort(key=lambda f: (str(feature_props[f].get("ID") or ""), f))

    visited_fids: set[int] = set()
    result: List[Tuple[int, List[float], Dict[str, Any]]] = []

    # shared-node NAP de-dup (by rounded XY)
    seen_vault_xy: set[Tuple[int, int]] = set()
    def _round_xy(pt: List[float], scale: int = 8) -> Tuple[int, int]:
        return (int(round(pt[0] * (10 ** scale))), int(round(pt[1] * (10 ** scale))))

    def _emit_vaults_until(fid: int, up_to_along_ft: float, v_state: Dict[int, int]):
        """
        Emit vaults (NAPs) on feature fid up to the oriented 'along' distance reached so far.
        Orientation aware (handles walking a feature in reverse).
        """
        lst = v_by_feature.get(fid, [])
        if not lst:
            return

        flip = getattr(_emit_vaults_until, "flip", {}).get(fid, False)
        total = getattr(_emit_vaults_until, "total", {}).get(fid, None)

        if total is None:
            coords, _ = lines[fid]
            tl = 0.0
            for i in range(1, len(coords)):
                tl += distance_feet(coords[i - 1], coords[i])
            total = tl
            _emit_vaults_until.total = getattr(_emit_vaults_until, "total", {})
            _emit_vaults_until.total[fid] = total

        def to_oriented(entry):
            orig_along, vi, vpt, vprops = entry
            oriented = (total - orig_along) if flip else orig_along
            return (oriented, vi, vpt, vprops)

        oriented = [to_oriented(t) for t in lst]
        oriented.sort(key=lambda t: t[0])

        i = v_state.get(fid, 0)
        if i and (i > len(oriented) or oriented[i - 1][0] > up_to_along_ft + 1e-6):
            already = oriented[:i]
            last_along = already[-1][0] if already else -1.0
            i = 0
            while i < len(oriented) and oriented[i][0] <= last_along + 1e-6:
                i += 1

        while i < len(oriented) and oriented[i][0] <= up_to_along_ft + 1e-6:
            _, vi, vpt, vprops = oriented[i]
            key = _round_xy(vpt)
            if key not in seen_vault_xy:
                seen_vault_xy.add(key)
                result.append((vi, vpt, vprops))
            i += 1

        v_state[fid] = i

    def _emit_remaining_vaults(fid: int, v_state: Dict[int, int]):
        """Emit any vaults (NAPs) left on this feature in oriented order at the end of the walk."""
        lst = v_by_feature.get(fid, [])
        if not lst:
            return

        flip = getattr(_emit_vaults_until, "flip", {}).get(fid, False)
        total = getattr(_emit_vaults_until, "total", {}).get(fid, None)
        if total is None:
            coords, _ = lines[fid]
            tl = 0.0
            for i in range(1, len(coords)):
                tl += distance_feet(coords[i - 1], coords[i])
            total = tl
            _emit_vaults_until.total = getattr(_emit_vaults_until, "total", {})
            _emit_vaults_until.total[fid] = total

        def to_oriented(entry):
            orig_along, vi, vpt, vprops = entry
            oriented = (total - orig_along) if flip else orig_along
            return (oriented, vi, vpt, vprops)

        oriented = [to_oriented(t) for t in lst]
        oriented.sort(key=lambda t: t[0])

        i = v_state.get(fid, 0)
        while i < len(oriented):
            _, vi, vpt, vprops = oriented[i]
            key = _round_xy(vpt)
            if key not in seen_vault_xy:
                seen_vault_xy.add(key)
                result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i

    def _same_cable(a: int, b: int) -> bool:
        """
        Two features are the *same cable* iff:
        - Section matches (case-insensitive, trimmed),
        - Distribution Type numeric matches,
        - Fiber Letter #1 matches when BOTH present; if either side is missing a letter, still treat as same cable.
        Conduit Type is intentionally ignored for 'same cable' detection.

        Rationale (matches your NAP 1–24 rules):
        • "Two 48ct 'A'" at a vault = same cable (continue straight unless detouring),
        • Equal-24 detours prefer 'B' first, then 'A',
        • Do NOT split a same-Section/same-Distribution line just because letter is absent on one segment
            or because Conduit Type differs in the attribute table.
        """
        pa = feature_props[a]
        pb = feature_props[b]

        # Section must match
        sec_a = (pa.get("Section") or "").strip().lower()
        sec_b = (pb.get("Section") or "").strip().lower()
        if sec_a != sec_b:
            return False

        # Numeric Distribution Type must match (72 / 48 / 24)
        if _dist_num(pa.get("Distribution Type")) != _dist_num(pb.get("Distribution Type")):
            return False

        # Letter logic: if both sides have a letter, they must match; if either is missing, still "same cable"
        la = (pa.get("Fiber Letter #1") or "").strip().upper()
        lb = (pb.get("Fiber Letter #1") or "").strip().upper()
        if la and lb:
            return la == lb
        return True


    def _last_seg_vec(coords: List[List[float]], forward: bool) -> Tuple[float, float]:
        if forward:
            p1, p2 = coords[-2], coords[-1]
        else:
            p1, p2 = coords[1], coords[0]
        return _vec_meters(p1, p2)

    def _pick_next_same_cable(fid: int, tail_node: Tuple[int, int]) -> Optional[int]:
        candidates = [
            nf
            for nf in node_to_features.get(tail_node, set())
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
                if tail_node in nnodes:
                    k = nnodes.index(tail_node)
                    if 0 < k < len(nnodes) - 1:
                        bx, by = _vec_meters(_nk_to_xy(nnodes[k]), _nk_to_xy(nnodes[k + 1]))
                    else:
                        continue
                else:
                    continue
            ang = _angle_between(ax, ay, bx, by)
            if ang < best[0]:
                best = (ang, nf)
        return best[1]

    def _min_dist_to_feature(xy: List[float], other_fid: int) -> float:
        ocoords, _ = lines[other_fid]
        best = float("inf")
        for i in range(1, len(ocoords)):
            d = point_segment_distance_feet(xy, ocoords[i - 1], ocoords[i])
            if d < best:
                best = d
        return best

    def _neighbors_at_node(fid: int, nk: Tuple[int, int]) -> List[int]:
        """Neighbors at a node in the same Section but different cable; includes snapped tees within NODE_SNAP_FT."""
        props = feature_props[fid]
        xy = _nk_to_xy(nk)
        shared = {
            nf
            for nf in (node_to_features.get(nk, set()) or set())
            if nf != fid
            and nf not in visited_fids
            and (_sec_norm_local(props.get("Section")) == _sec_norm_local(feature_props[nf].get("Section")))
            and not _same_cable(fid, nf)
        }
        snapped = set()
        for nf in (by_section.get(_sec_norm_local(props.get("Section")), []) or []):
            if nf == fid or nf in visited_fids or _same_cable(fid, nf):
                continue
            if _min_dist_to_feature(xy, nf) <= NODE_SNAP_FT:
                snapped.add(nf)
        return list(shared | snapped)

    v_state: Dict[int, int] = {}

    def _walk_feature(fid: int, entry_node: Tuple[int, int]):
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        props = feature_props[fid]
        cur_dist = _dist_num(props.get("Distribution Type"))

        flip = bool(nodes and nodes[-1] == entry_node)

        if not hasattr(_emit_vaults_until, "flip"):
            _emit_vaults_until.flip = {}
        _emit_vaults_until.flip[fid] = flip

        if not hasattr(_emit_vaults_until, "total"):
            _emit_vaults_until.total = {}
        if fid not in _emit_vaults_until.total:
            tl = 0.0
            for i in range(1, len(coords)):
                tl += distance_feet(coords[i - 1], coords[i])
            _emit_vaults_until.total[fid] = tl

        along_orig = [0.0]
        for i in range(1, len(coords)):
            along_orig.append(along_orig[-1] + distance_feet(coords[i - 1], coords[i]))
        total_len = along_orig[-1] if along_orig else 0.0

        along_at_node_orig = {nodes[i]: along_orig[i] for i in range(len(nodes))}
        if not flip:
            along_at_node = along_at_node_orig
        else:
            along_at_node = {nk: (total_len - a) for nk, a in along_at_node_orig.items()}

        prev_xy: Optional[List[float]] = None
        for nk in nodes:
            _emit_vaults_until(fid, along_at_node[nk], v_state)

            at_xy = _nk_to_xy(nk)
            is_start = bool(start_pair and fid == start_pair[0] and nk == start_pair[1])

            neighbor_fids = _neighbors_at_node(fid, nk)
            neighbors = list(neighbor_fids)
            neighbors.sort(
                key=lambda fid2: _neighbor_rank(fid, prev_xy, at_xy, fid2, is_start)
            )
            neighbor_fids = neighbors

            lowers: List[int] = []
            equal24: List[int] = []
            equal_non24: List[int] = []
            higher: List[int] = []
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

            for nf in sorted(
                lowers,
                key=lambda f: (_dist_num(feature_props[f].get("Distribution Type")), str(feature_props[f].get("ID") or ""), f),
            ):
                _walk_feature(nf, nk)

            def _letter_rank_for_equal24(letter: str) -> Tuple[int, str]:
                L = (letter or "").strip().upper()
                if L == "B":
                    return (0, L)
                if L == "A":
                    return (1, L)
                return (2, L)

            for nf in sorted(
                equal24,
                key=lambda f: _letter_rank_for_equal24(feature_props[f].get("Fiber Letter #1")) + (str(feature_props[f].get("ID") or ""),),
            ):
                _walk_feature(nf, nk)

            # Continue along SAME CABLE from the tail (straightest turn)
            tail = nodes[-1] if not flip else nodes[0]
            next_same = _pick_next_same_cable(fid, tail)
            if next_same is not None:
                _walk_feature(next_same, tail)

            # Deferred same-Section ties at endpoints
            for nk2 in (nodes[0], nodes[-1]):
                neighbor_fids2 = _neighbors_at_node(fid, nk2)
                en24 = [nf for nf in neighbor_fids2 if _dist_num(feature_props[nf].get("Distribution Type")) == cur_dist != 24]
                hi = [nf for nf in neighbor_fids2 if _dist_num(feature_props[nf].get("Distribution Type")) > cur_dist]

                for nf in sorted(
                    en24,
                    key=lambda f: (str(feature_props[f].get("Fiber Letter #1") or ""), str(feature_props[f].get("ID") or ""), f),
                ):
                    _walk_feature(nf, nk2)

                for nf in sorted(
                    hi,
                    key=lambda f: (
                        -_dist_num(feature_props[f].get("Distribution Type")),
                        str(feature_props[f].get("Fiber Letter #1") or ""),
                        str(feature_props[f].get("ID") or ""),
                        f,
                    ),
                ):
                    _walk_feature(nf, nk2)

            _emit_remaining_vaults(fid, v_state)
            prev_xy = at_xy

    # ===== Master traversal order, T3-anchored =====
    start_nk = _find_start_node_from_t3(t3_gj, feature_nodes)
    preferred_section: Optional[str] = None
    start_pair: Optional[Tuple[int, Tuple[int, int]]] = None

    if start_nk is not None:
        touching = [nf for nf in (node_to_features.get(start_nk, set()) or [])]
        if touching:
            start_fid = _choose_start_fid(touching)
            preferred_section = (feature_props[start_fid].get("Section") or "").strip() or None
            start_pair = (start_fid, start_nk)

    sections_to_run: List[Optional[str]] = []
    if preferred_section is not None and preferred_section in by_section:
        sections_to_run.append(preferred_section)
    for s in [f"Section {i}" for i in range(1, 7)]:
        if s in by_section and s not in sections_to_run:
            sections_to_run.append(s)
    if None in by_section:
        sections_to_run.append(None)

    for sec in sections_to_run:
        if start_pair and preferred_section == sec and start_pair[0] not in visited_fids:
            _walk_feature(start_pair[0], start_pair[1])

        while True:
            candidates = [f for f in (by_section.get(sec, []) or []) if f not in visited_fids]
            if not candidates:
                break
            fid0 = candidates[0]
            nk0 = feature_nodes[fid0][0]
            _walk_feature(fid0, nk0)

    return result
