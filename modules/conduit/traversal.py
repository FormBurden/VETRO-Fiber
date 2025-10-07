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

    # PLACE THIS HELPER INSIDE dfs_ordered_vaults, near other locals

    def _neighbor_rank(current_fid: int,
                    came_from_xy: Optional[List[float]],
                    at_xy: List[float],
                    nb_fid: int,
                    is_start_node: bool) -> Tuple:
        """
        Rank neighbors according to your field rules:

        1) Immediately after arriving at a node: we will emit NAPs first (handled elsewhere).
        2) Detour preference:
        • LOWER Distribution first (smaller fiber count first).
        • If CURRENT is 24ct, also detour equal-24 now in order: B → A → others.
        3) Continue SAME CABLE (Section+Distribution+Letter#1) with straightest turn.
        4) Remaining ties by section order and then ID for stability.

        We encode this as a sortable key (smaller = higher priority).
        """
        pa = feature_props[current_fid]
        pb = feature_props[nb_fid]

        cur_dist = _dist_num(pa.get("Distribution Type"))
        nb_dist  = _dist_num(pb.get("Distribution Type"))

        cur_letter = (pa.get("Fiber Letter #1") or "").strip().upper()
        nb_letter  = (pb.get("Fiber Letter #1") or "").strip().upper()

        same_line = _same_cable(current_fid, nb_fid)

        # Section bucket (stable final tie-breaker)
        sec = _sec_norm(pb.get("Section"))
        try:
            sec_bucket = section_order.index(sec) if sec in section_order else len(section_order)
        except Exception:
            sec_bucket = len(section_order)

        # Straightness (prefer straighter)
        # If we don't have a previous direction vector, set a neutral angle.
        if came_from_xy is not None:
            ax, ay = _vec_meters(came_from_xy, at_xy)
        else:
            ax, ay = (0.0, 1.0)  # arbitrary "up" to keep key consistent

        # nb direction: from node (at_xy) toward the first vertex of nb feature away from the node
        nb_coords, nb_nodes, _ = _feature_along_vertices(nb_fid)
        # pick the endpoint in nb_coords that is farther from at_xy to define the outgoing vector
        if distance_feet(at_xy, nb_coords[0]) < distance_feet(at_xy, nb_coords[-1]):
            bx, by = _vec_meters(at_xy, nb_coords[-1])
        else:
            bx, by = _vec_meters(at_xy, nb_coords[0])
        angle = _angle_between(ax, ay, bx, by)

        # Priority layers:

        # A) SAME CABLE continuation should generally be preferred after required detours
        same_line_flag = 0 if same_line else 1

        # B) Detours: lower distribution first
        # encode as: lower nb_dist → smaller key
        lower_dist_key = nb_dist

        # C) Special 24ct equal handling when current is 24:
        # if equal 24: B first, then A, then others
        equal24_order = 9  # default high
        if cur_dist == 24 and nb_dist == 24:
            if nb_letter == "B":
                equal24_order = 0
            elif nb_letter == "A":
                equal24_order = 1
            else:
                equal24_order = 2

        # D) Straightness (smaller angle is straighter)
        straight_key = angle

        # E) Final tie-breakers
        id_key = str(pb.get("ID") or "")

        # Build the rank tuple. The ordering logic is:
        #   - Prefer required detours (lower dist) before same-line continuation,
        #     but if both are same-line, straightness wins.
        #
        # We do this by making a composite that first prefers lower distribution,
        # then (only among equals) prefer equal-24 ordering, then whether it's same-line,
        # then straightness, then section bucket, then ID.
        #
        # Note: is_start_node can be used if you want to tweak behavior at the very first step
        # (e.g., prefer highest dist at start as you specified elsewhere), but our start-fid
        # choice already handles that so we keep the same ranking here.
        return (
            lower_dist_key,       # lower first
            equal24_order,        # 24ct B→A→others when equal 24
            same_line_flag,       # same line earlier
            straight_key,         # straighter earlier
            sec_bucket,           # section stability
            id_key,               # id stability
        )


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


    def _walk_feature(fid: int,
                    enter_node: Tuple[int, int],
                    came_from_xy: Optional[List[float]],
                    v_state: Dict[int, int]) -> None:
        """
        Walk a single feature from the node we just entered.
        Critical fixes:
        • Set orientation so we walk AWAY from the junction (tie-point rule).
        • Immediately EMIT vaults/NAPs at this node BEFORE any branch traversal.
        • At every subsequent node along this feature, emit the vault(s) at that node first,
            then consider neighbor features in the _neighbor_rank() order.
        """
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, nodes, along = _feature_along_vertices(fid)

        # Find which vertex/segment matches the enter_node; choose direction AWAY from it.
        # nodes is a list of node-keys per vertex; locate matching index.
        try:
            enter_idx = nodes.index(enter_node)
        except ValueError:
            # If we can't find it (shouldn't happen), pick the closer endpoint by distance.
            d0 = distance_feet(_nk_to_xy(nodes[0]), _nk_to_xy(enter_node))
            d1 = distance_feet(_nk_to_xy(nodes[-1]), _nk_to_xy(enter_node))
            enter_idx = 0 if d0 < d1 else (len(nodes) - 1)

        # Orientation: if we entered at head (index 0), we go forward; if at tail (last), we go backward.
        go_forward = (enter_idx == 0)
        # If we entered somewhere in the middle (rare), pick the side that increases distance from the entry.
        if 0 < enter_idx < (len(nodes) - 1):
            # compare remaining length forward vs backward; choose the longer (towards cable "end")
            remain_fwd = along[-1] - along[enter_idx]
            remain_back = along[enter_idx]
            go_forward = remain_fwd >= remain_back

        # Record flip meta for _emit_vaults_until/_emit_remaining_vaults
        flip_map = getattr(_emit_vaults_until, "flip", {})
        flip_map = dict(flip_map)
        flip_map[fid] = (not go_forward)  # flipped if we are walking tail→head
        _emit_vaults_until.flip = flip_map

        # Emit vault(s) AT THE ENTRY NODE immediately (up_to_along=0 in oriented space).
        _emit_vaults_until(fid, 0.0, v_state)

        # Start walking along this feature, vertex by vertex, in the chosen direction.
        if go_forward:
            vertex_iter = range(enter_idx + 1, len(nodes))
            prev_xy = coords[enter_idx]
        else:
            vertex_iter = range(enter_idx - 1, -1, -1)
            prev_xy = coords[enter_idx]

        for vi in vertex_iter:
            cur_xy = coords[vi]
            # Compute oriented distance “so far” on this feature
            oriented_along = along[vi] if go_forward else (along[-1] - along[vi])
            # Emit all vaults up to this oriented position (this includes node-coincident vaults)
            _emit_vaults_until(fid, oriented_along, v_state)

            # At this node, consider neighbor features (branches) BEFORE continuing.
            nk = nodes[vi]
            touching = [tfid for tfid in node_to_features.get(nk, []) if tfid != fid]

            if touching:
                # Order neighbors per rules.
                is_start_node = (came_from_xy is None)
                at_xy = cur_xy
                neighbors = list(touching)
                neighbors.sort(
                    key=lambda nb: _neighbor_rank(fid, came_from_xy, at_xy, nb, is_start_node)
                )

                # Depth-first into neighbors in ranked order.
                for nb in neighbors:
                    if nb in visited_fids:
                        continue
                    # WALK the neighbor, entering at this same node “nk”.
                    _walk_feature(nb, nk, at_xy, v_state)

            # update prev direction for straightness calc on subsequent steps
            came_from_xy = prev_xy
            prev_xy = cur_xy

        # At the end of this feature, emit any remaining vaults (safety).
        _emit_remaining_vaults(fid, v_state)


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
