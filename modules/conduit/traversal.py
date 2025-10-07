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
    Depth-first traversal of the conduit network producing vaults in the field rule order.

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

        # all nodes in the network
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

    def _dist_num(s: Optional[str]) -> int:
        if not s:
            return 0
        m = re.search(r"(\d+)", s)
        return int(m.group(1)) if m else 0

    def _sec_norm(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        s = s.strip()
        return s if s and s.lower() != "null" else None

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
    for fid in feature_nodes.keys():
        sec = _sec_norm(feature_props[fid].get("Section"))
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
        If this feature is being walked in reverse (tail→head), remap each vault's stored
        'along_ft' (which is in ORIGINAL feature order) to ORIENTED space:

            oriented_along = total_len - original_along

        We keep per-feature 'flip' and 'total_len' in attributes on this function
        (they are written by _walk_feature when the walk for that fid starts).
        """
        lst = v_by_feature.get(fid, [])
        if not lst:
            return

        # Look up orientation meta for this feature (populated by _walk_feature)
        flip = getattr(_emit_vaults_until, "flip", {}).get(fid, False)
        total = getattr(_emit_vaults_until, "total", {}).get(fid, None)

        # If we don't have total, compute once (sum of segment lengths in ORIGINAL order)
        if total is None:
            coords, _ = lines[fid]
            tl = 0.0
            for i in range(1, len(coords)):
                tl += distance_feet(coords[i - 1], coords[i])
            total = tl
            # Cache so subsequent calls are consistent
            _emit_vaults_until.total = getattr(_emit_vaults_until, "total", {})
            _emit_vaults_until.total[fid] = total

        # Build a view of this feature's vaults in the ORIENTED space
        # Each entry: (oriented_along, vi, vpt, vprops)
        def to_oriented(entry):
            orig_along, vi, vpt, vprops = entry
            oriented = (total - orig_along) if flip else orig_along
            return (oriented, vi, vpt, vprops)

        oriented = [to_oriented(t) for t in lst]
        oriented.sort(key=lambda t: t[0])  # ensure increasing along in oriented space

        # Emit in oriented order up to our current position
        i = v_state.get(fid, 0)
        # If v_state was recorded against a previous sorting, re-find our position
        if i and (i > len(oriented) or oriented[i - 1][0] > up_to_along_ft + 1e-6):
            # resync by finding first index with along > already-emitted last along
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
        """
        Emit any vaults (NAPs) left on this feature in oriented order at the end of the walk.
        Orientation handling matches _emit_vaults_until (see that docstring).
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
        while i < len(oriented):
            _, vi, vpt, vprops = oriented[i]
            key = _round_xy(vpt)
            if key not in seen_vault_xy:
                seen_vault_xy.add(key)
                result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i


    def _same_cable(a: int, b: int) -> bool:
        pa = feature_props[a]
        pb = feature_props[b]
        sec_eq = (pa.get("Section") or "").strip().lower() == (pb.get("Section") or "").strip().lower()
        dist_eq = _dist_num(pa.get("Distribution Type")) == _dist_num(pb.get("Distribution Type"))
        let_a = (pa.get("Fiber Letter #1") or "").strip().upper()
        let_b = (pb.get("Fiber Letter #1") or "").strip().upper()
        # SAME line if Section + Distribution match; if letters exist on both and match, good.
        # If letters missing on either side, still same line. Conduit Type differences are ignored here.
        return sec_eq and dist_eq and (let_a == let_b or not (let_a and let_b))

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
        # true shared-vertex neighbors
        shared = {
            nf
            for nf in (node_to_features.get(nk, set()) or set())
            if nf != fid
            and nf not in visited_fids
            and (_sec_norm(props.get("Section")) == _sec_norm(feature_props[nf].get("Section")))
            and not _same_cable(fid, nf)
        }
        # snapped tees (nearby mid-segment) in same Section
        snapped = set()
        for nf in (by_section.get(_sec_norm(props.get("Section")), []) or []):
            if nf == fid or nf in visited_fids or _same_cable(fid, nf):
                continue
            if _min_dist_to_feature(xy, nf) <= NODE_SNAP_FT:
                snapped.add(nf)
        return list(shared | snapped)

    def _letter_rank_for_equal24(letter: str) -> Tuple[int, str]:
        L = (letter or "").strip().upper()
        if L == "B":
            return (0, L)
        if L == "A":
            return (1, L)
        return (2, L)

    v_state: Dict[int, int] = {}

    def _walk_feature(fid: int, entry_node: Tuple[int, int]):
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, _ = lines[fid]
        nodes = feature_nodes[fid]
        props = feature_props[fid]
        cur_dist = _dist_num(props.get("Distribution Type"))

        # Determine if our walk direction is reversed vs ORIGINAL feature order
        # We consider "forward" to be ORIGINAL nodes[0] -> nodes[-1].
        # If the entry node is the last node, we're walking in reverse.
        flip = bool(nodes and nodes[-1] == entry_node)

        # Cache this orientation for emission helpers
        if not hasattr(_emit_vaults_until, "flip"):
            _emit_vaults_until.flip = {}
        _emit_vaults_until.flip[fid] = flip

        # Pre-compute and cache total length once (used by emission helpers)
        if not hasattr(_emit_vaults_until, "total"):
            _emit_vaults_until.total = {}
        if fid not in _emit_vaults_until.total:
            tl = 0.0
            for i in range(1, len(coords)):
                tl += distance_feet(coords[i - 1], coords[i])
            _emit_vaults_until.total[fid] = tl

        # Build along-distance per *ORIENTED* vertex position for thresholding:
        # If flip, we want along_at_node increasing from the entry end toward the far end.
        # We'll compute along in ORIGINAL space, then convert to oriented threshold:
        along_orig = [0.0]
        for i in range(1, len(coords)):
            along_orig.append(along_orig[-1] + distance_feet(coords[i - 1], coords[i]))
        total_len = along_orig[-1] if along_orig else 0.0

        # Map ORIGINAL node keys to along-from-start
        along_at_node_orig = {nodes[i]: along_orig[i] for i in range(len(nodes))}

        # Convert to ORIENTED threshold: if walking in reverse, oriented_along = total - original_along
        if not flip:
            along_at_node = along_at_node_orig
        else:
            along_at_node = {nk: (total_len - a) for nk, a in along_at_node_orig.items()}

        # Traverse vertices in ORIGINAL order (nodes as-is); thresholds are in oriented space,
        # so emission will still happen in the correct walk direction.
        for nk in nodes:
            _emit_vaults_until(fid, along_at_node[nk], v_state)

            neighbor_fids = _neighbors_at_node(fid, nk)
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

        # Continue along SAME CABLE from the tail (same as before)
        tail = nodes[-1] if not flip else nodes[0]
        next_same = _pick_next_same_cable(fid, tail)
        if next_same is not None:
            _walk_feature(next_same, tail)

        # Deferred same-Section ties at endpoints
        for nk in (nodes[0], nodes[-1]):
            neighbor_fids = _neighbors_at_node(fid, nk)
            en24 = [nf for nf in neighbor_fids if _dist_num(feature_props[nf].get("Distribution Type")) == cur_dist != 24]
            hi = [nf for nf in neighbor_fids if _dist_num(feature_props[nf].get("Distribution Type")) > cur_dist]

            for nf in sorted(
                en24,
                key=lambda f: (str(feature_props[f].get("Fiber Letter #1") or ""), str(feature_props[f].get("ID") or ""), f),
            ):
                _walk_feature(nf, nk)

            for nf in sorted(
                hi,
                key=lambda f: (
                    -_dist_num(feature_props[f].get("Distribution Type")),
                    str(feature_props[f].get("Fiber Letter #1") or ""),
                    str(feature_props[f].get("ID") or ""),
                    f,
                ),
            ):
                _walk_feature(nf, nk)

        _emit_remaining_vaults(fid, v_state)



    # ===== Master traversal order, T3-anchored =====
    # Prefer the Section that contains the T3-touching start feature first
    start_nk = _find_start_node_from_t3(t3_gj, feature_nodes)
    preferred_section: Optional[str] = None
    start_pair: Optional[Tuple[int, Tuple[int, int]]] = None

    if start_nk is not None:
        touching = [nf for nf in (node_to_features.get(start_nk, set()) or [])]
        if touching:
            start_fid = _choose_start_fid(touching)
            preferred_section = _sec_norm(feature_props[start_fid].get("Section"))
            start_pair = (start_fid, start_nk)

    # Build the section iteration list: preferred Section first (if any), then Section 1..6, then None
    sections_to_run: List[Optional[str]] = []
    if preferred_section is not None and preferred_section in by_section:
        sections_to_run.append(preferred_section)
    for s in section_order:
        if s in by_section and s not in sections_to_run:
            sections_to_run.append(s)
    if None in by_section:
        sections_to_run.append(None)

    # Walk each Section's connected components
    for sec in sections_to_run:
        # Walk preferred start in this section if it matches
        if start_pair and preferred_section == sec and start_pair[0] not in visited_fids:
            _walk_feature(start_pair[0], start_pair[1])

        # Then walk any remaining components in this Section
        while True:
            candidates = [f for f in (by_section.get(sec, []) or []) if f not in visited_fids]
            if not candidates:
                break
            # Try to pick an endpoint start (nicer ordering); else just pick the first remaining
            fid0 = candidates[0]
            nk0 = feature_nodes[fid0][0]
            _walk_feature(fid0, nk0)

    return result