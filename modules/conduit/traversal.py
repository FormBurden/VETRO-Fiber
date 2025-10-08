# modules/conduit/traversal.py
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from modules.geo.geometry import distance_feet
from modules.conduit.topology import (
    build_conduit_topology,
    assign_vaults_to_features,
    iter_lines,
)


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
      4) Continue along SAME CABLE (same Section/Distribution/Letter#1).
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
    seen_vault_ids: set[Tuple[int, int]] = set()
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
            key = (fid, vi)
            if key not in seen_vault_ids:
                seen_vault_ids.add(key)
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
            key = (fid, vi)
            if key not in seen_vault_ids:
                seen_vault_ids.add(key)
                result.append((vi, vpt, vprops))
            i += 1
        v_state[fid] = i

    def _emit_node_vaults(nk: Tuple[int, int], tol_ft: float = 1.0) -> None:
        """
        Emit any vault(s)/NAP(s) that lie within tol_ft of the given node coordinate,
        regardless of which feature they are assigned to.

        This guarantees the start T3 node becomes NAP 1 if a NAP exists there
        (even when the point is not exactly coincident with the snapped node).
        """
        xy = _nk_to_xy(nk)

        # check every feature that touches this node
        for fid in node_to_features.get(nk, []):
            lst = v_by_feature.get(fid, [])
            if not lst:
                continue

            # Ensure oriented metadata exists for this fid (so later calls don't recompute)
            if fid not in getattr(_emit_vaults_until, "total", {}):
                coords, _ = lines[fid]
                tl = 0.0
                for i in range(1, len(coords)):
                    tl += distance_feet(coords[i - 1], coords[i])
                _emit_vaults_until.total = getattr(_emit_vaults_until, "total", {})
                _emit_vaults_until.total[fid] = tl

            # Emit all vaults that lie within tol_ft of the node XY
            for _, vi, vpt, vprops in lst:
                if distance_feet(vpt, xy) <= tol_ft:
                    key = (fid, vi)
                    if key not in seen_vault_ids:
                        seen_vault_ids.add(key)
                        result.append((vi, vpt, vprops))


    def _same_cable(a: int, b: int) -> bool:
        """
        Two features are the *same cable* iff:
        - Section matches (case-insensitive, trimmed),
        - Distribution Type numeric matches,
        - Fiber Letter #1 matches when BOTH present; if either side is missing a letter, still treat as same cable.
        Conduit Type is intentionally ignored for 'same cable' detection.
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

    def _neighbor_rank(
        current_fid: int,
        came_from_xy: Optional[List[float]],
        at_xy: List[float],
        nb_fid: int,
        is_start_node: bool,
        prefer_same_line_first: bool,
    ) -> Tuple:
        """
        Attribute-only neighbor ranking (no geometry).

        Order logic:
          A) If prefer_same_line_first=True (current feature ENDS at this node),
             then SAME CABLE neighbors get top priority; otherwise push same-cable back.
          B) Detours: LOWER Distribution first (24 < 48 < 72).
          C) If CURRENT is 24ct and neighbor is also 24ct → prefer B, then A, then others.
          D) Among different cables within the same section, prefer 'A' first, then alphabetical.
          E) Stable fallbacks: section bucket → ID.
        """
        pa = feature_props[current_fid]
        pb = feature_props[nb_fid]

        cur_dist = _dist_num(pa.get("Distribution Type"))
        nb_dist  = _dist_num(pb.get("Distribution Type"))

        nb_letter  = (pb.get("Fiber Letter #1") or "").strip().upper()

        # A) same-cable priority depends on end-of-feature status
        same_line = _same_cable(current_fid, nb_fid)
        if prefer_same_line_first:
            same_line_key = 0 if same_line else 1
        else:
            same_line_key = 1 if same_line else 0

        # B) lower distribution first
        lower_dist_key = nb_dist

        # C) special 24/24 ordering: B → A → others
        equal24_order = 9
        if cur_dist == 24 and nb_dist == 24:
            if nb_letter == "B":
                equal24_order = 0
            elif nb_letter == "A":
                equal24_order = 1
            else:
                equal24_order = 2

        # D) within different-cable ties, A-first then alphabetical
        letter_bucket = (0 if nb_letter == "A" else 1, nb_letter)

        # E) stable section + id
        sec = _sec_norm(pb.get("Section"))
        try:
            sec_bucket = section_order.index(sec) if sec in section_order else len(section_order)
        except Exception:
            sec_bucket = len(section_order)

        id_key = str(pb.get("ID") or "")

        return (
            same_line_key,   # A
            lower_dist_key,  # B
            equal24_order,   # C
            letter_bucket,   # D
            sec_bucket,      # E
            id_key,          # E
        )

    v_state: Dict[int, int] = {}

    def _walk_feature(
        fid: int,
        enter_node: Tuple[int, int],
        came_from_xy: Optional[List[float]],
        v_state: Dict[int, int],
    ) -> None:
        """
        Walk a single feature from the node we just entered.

        Key behavior:
          • Emit NAPs at the entry node.
          • While traversing along this feature, at each node:
              - If the feature continues past this node → ONLY detour branches (exclude same-cable).
              - If the feature ENDS at this node → prefer SAME CABLE neighbor first (continue backbone),
                then other branches per attribute rules.
        """
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, nodes, along = _feature_along_vertices(fid)

        # Find which vertex/segment matches the enter_node; choose direction AWAY from it.
        try:
            enter_idx = nodes.index(enter_node)
        except ValueError:
            d0 = distance_feet(_nk_to_xy(nodes[0]), _nk_to_xy(enter_node))
            d1 = distance_feet(_nk_to_xy(nodes[-1]), _nk_to_xy(enter_node))
            enter_idx = 0 if d0 < d1 else (len(nodes) - 1)

        # Orientation: head→tail if we entered at head; tail→head if we entered at tail.
        go_forward = (enter_idx == 0)
        if 0 < enter_idx < (len(nodes) - 1):
            # middle entry: pick the side with more remaining length
            remain_fwd = along[-1] - along[enter_idx]
            remain_back = along[enter_idx]
            go_forward = remain_fwd >= remain_back

        # Record flip meta for the emitters (so along is oriented correctly)
        flip_map = getattr(_emit_vaults_until, "flip", {})
        flip_map = dict(flip_map)
        flip_map[fid] = (not go_forward)  # flipped if walking tail→head
        _emit_vaults_until.flip = flip_map

        # Emit NAP(s) at entry node
        _emit_vaults_until(fid, 0.0, v_state)

        # Iterate vertices in walking direction
        if go_forward:
            vertex_iter = range(enter_idx + 1, len(nodes))
            prev_xy = coords[enter_idx]
            tail_index = len(nodes) - 1
        else:
            vertex_iter = range(enter_idx - 1, -1, -1)
            prev_xy = coords[enter_idx]
            tail_index = 0

        for vi in vertex_iter:
            cur_xy = coords[vi]

            # Oriented-along distance at this vertex
            oriented_along = along[vi] if go_forward else (along[-1] - along[vi])
            _emit_vaults_until(fid, oriented_along, v_state)

            # Node bookkeeping
            nk = nodes[vi]
            touching = [tfid for tfid in node_to_features.get(nk, []) if tfid != fid and tfid not in visited_fids]

            if touching:
                # Determine if the current feature ENDS at this node
                ends_here = (vi == tail_index)

                # If the feature continues past this node, EXCLUDE same-cable neighbors here.
                # (We stay on this feature to continue the backbone; only take detours.)
                nbrs = list(touching)
                if not ends_here:
                    nbrs = [nb for nb in nbrs if not _same_cable(fid, nb)]

                # Rank neighbors; when at an end, same-cable first; otherwise push same-cable back.
                neighbors = sorted(
                    nbrs,
                    key=lambda nb: _neighbor_rank(
                        fid,
                        came_from_xy,
                        cur_xy,
                        nb,
                        is_start_node=(came_from_xy is None),
                        prefer_same_line_first=ends_here,
                    ),
                )

                # Depth-first into neighbors in ranked order.
                for nb in neighbors:
                    if nb in visited_fids:
                        continue
                    _walk_feature(nb, nk, cur_xy, v_state)

            # advance prior-dir anchor
            came_from_xy = prev_xy
            prev_xy = cur_xy

        # At the end of this feature, emit any remaining vaults (safety)
        _emit_remaining_vaults(fid, v_state)

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

            # Emit the NAP(s) at the T3 node BEFORE any walking
            _emit_node_vaults(start_nk)

    # Build the section run order (preferred section first, then the rest, then None)
    sections_to_run: List[Optional[str]] = []
    if preferred_section is not None and preferred_section in by_section:
        sections_to_run.append(preferred_section)
    for s in [f"Section {i}" for i in range(1, 7)]:
        if s in by_section and s not in sections_to_run:
            sections_to_run.append(s)
    if None in by_section:
        sections_to_run.append(None)

    # Walk, starting with the chosen start feature at the T3 node (if any)
    for sec in sections_to_run:
        if start_pair and preferred_section == sec and start_pair[0] not in visited_fids:
            _walk_feature(start_pair[0], start_pair[1], None, {})  # came_from_xy=None, fresh v_state

        # Then finish any remaining features in this section
        while True:
            candidates = [f for f in (by_section.get(sec, []) or []) if f not in visited_fids]
            if not candidates:
                break
            fid0 = candidates[0]
            nk0 = feature_nodes[fid0][0]
            _walk_feature(fid0, nk0, None, {})

    return result
