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
    Traverse the conduit network and emit vaults (NAPs) in the desired field order.

    START RULE (simple + explicit):
      • At the T3 node, look at the touching conduits and read their "Section".
      • Start with Section 1, walk ALL touching features that are Section 1 (ordering is deterministic).
      • Then Section 2, then Section 3, … (in strictly increasing section number).
      • No clockwise/CCW tie-breakers at T3.

    DURING TRAVERSAL (unchanged from your rules):
      • If the current feature CONTINUES through a node → only take detours (exclude same-cable) before continuing.
      • If the current feature ENDS at a node → prefer SAME-CABLE first (continue backbone), then other branches.
      • Detour preference: LOWER Distribution first (24 < 48 < 72).
      • 24/24 special: prefer B → A → others.
      • Within different-cable ties in the same Section: 'A' first, then alphabetical.
      • Emit node-coincident NAPs with 20 ft tolerance when at the T3 node, and emit along-feature NAPs as you advance.
    """
    # --- Build topology and indices ---
    feature_nodes, node_to_features, feature_props = build_conduit_topology(conduit_gj)
    v_by_feature = assign_vaults_to_features(conduit_gj, feature_nodes, vaults_gj)
    lines = {fid: (coords, props) for fid, coords, props in iter_lines(conduit_gj)}

    # ----- helpers -----
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

    # Section parsing / normalization
    def _sec_norm(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        s = str(s).strip()
        return s if s and s.lower() != "null" else None

    def _sec_num(s: Optional[str]) -> int:
        """
        Extract a numeric section index from strings like 'Section 1', '1', 'Sec 2', etc.
        Unknown/None returns a large number so it sorts last.
        """
        if not s:
            return 10**9
        m = re.search(r"(\d+)", str(s))
        return int(m.group(1)) if m else 10**9

    # Build per-section feature buckets (used after T3 start group)
    by_section: Dict[Optional[str], List[int]] = {}
    for fid in feature_nodes.keys():
        sec = _sec_norm(feature_props[fid].get("Section"))
        by_section.setdefault(sec, []).append(fid)
    for lst in by_section.values():
        # Stable fallback ordering inside a section (ID then fid)
        lst.sort(key=lambda f: (str(feature_props[f].get("ID") or ""), f))

    visited_fids: set[int] = set()
    result: List[Tuple[int, List[float], Dict[str, Any]]] = []
    seen_vault_ids: set[Tuple[int, int]] = set()  # de-dup by (feature id, vault index)

    def _dist_num(s: Optional[str]) -> int:
        if not s:
            return 0
        m = re.search(r"(\d+)", s)
        return int(m.group(1)) if m else 0

    # ---------- Emitters ----------
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
            # Caller moved "backwards" relative to last emit; restart from 0
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

    def _emit_node_vaults(nk: Tuple[int, int], tol_ft: float = 20.0) -> None:
        """Emit any vault(s)/NAP(s) within tol_ft of this node, across all touching features."""
        xy = _nk_to_xy(nk)
        for fid in node_to_features.get(nk, []):
            lst = v_by_feature.get(fid, [])
            if not lst:
                continue
            if fid not in getattr(_emit_vaults_until, "total", {}):
                coords, _ = lines[fid]
                tl = 0.0
                for i in range(1, len(coords)):
                    tl += distance_feet(coords[i - 1], coords[i])
                _emit_vaults_until.total = getattr(_emit_vaults_until, "total", {})
                _emit_vaults_until.total[fid] = tl
            for _, vi, vpt, vprops in lst:
                if distance_feet(vpt, xy) <= tol_ft:
                    key = (fid, vi)
                    if key not in seen_vault_ids:
                        seen_vault_ids.add(key)
                        result.append((vi, vpt, vprops))

    # ---------- Cable identity + neighbor ranking (kept) ----------
    def _same_cable(a: int, b: int) -> bool:
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
        # Letter logic: if both sides have a letter, they must match; if either missing, still "same cable"
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
        Attribute-driven neighbor ranking (no geometry angles).
        Order:
          A) If prefer_same_line_first=True (current ENDS here), SAME-CABLE gets top priority; else push same-cable back.
          B) LOWER Distribution first (24 < 48 < 72).
          C) If CURRENT is 24 and neighbor is also 24 → prefer B → A → others.
          D) Among different cables within the same Section, prefer 'A' first, then alphabetical.
          E) Stable fallbacks: Section bucket (1..n, then None) → ID.
        """
        pa = feature_props[current_fid]
        pb = feature_props[nb_fid]
        cur_dist = _dist_num(pa.get("Distribution Type"))
        nb_dist = _dist_num(pb.get("Distribution Type"))
        nb_letter = (pb.get("Fiber Letter #1") or "").strip().upper()

        # A) same-cable priority conditioned on end-of-feature
        same_line = _same_cable(current_fid, nb_fid)
        if prefer_same_line_first:
            same_line_key = 0 if same_line else 1
        else:
            same_line_key = 1 if same_line else 0

        # B) lower distribution first
        lower_dist_key = nb_dist

        # C) special 24/24 preference: B → A → others
        equal24_order = 9
        if cur_dist == 24 and nb_dist == 24:
            if nb_letter == "B":
                equal24_order = 0
            elif nb_letter == "A":
                equal24_order = 1
            else:
                equal24_order = 2

        # D) letter preference within different-cable ties (A first, then alphabetical)
        letter_bucket = (0 if nb_letter == "A" else 1, nb_letter)

        # E) section + id fallback (ensures stable ordering)
        sec = _sec_norm(pb.get("Section"))
        sec_bucket = _sec_num(sec)  # smaller section number first, None/unknown last (big number)
        id_key = str(pb.get("ID") or "")

        return (
            same_line_key,   # A
            lower_dist_key,  # B  (smaller = higher priority)
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
        """
        if fid in visited_fids:
            return
        visited_fids.add(fid)

        coords, nodes, along = _feature_along_vertices(fid)

        # Determine entry orientation
        try:
            enter_idx = nodes.index(enter_node)
        except ValueError:
            # snap to closest end if the node isn't exactly an endpoint index
            d0 = distance_feet(_nk_to_xy(nodes[0]), _nk_to_xy(enter_node))
            d1 = distance_feet(_nk_to_xy(nodes[-1]), _nk_to_xy(enter_node))
            enter_idx = 0 if d0 < d1 else (len(nodes) - 1)

        go_forward = (enter_idx == 0)
        if 0 < enter_idx < (len(nodes) - 1):
            remain_fwd = along[-1] - along[enter_idx]
            remain_back = along[enter_idx]
            go_forward = remain_fwd >= remain_back

        # Record orientation for emitters
        flip_map = getattr(_emit_vaults_until, "flip", {})
        flip_map = dict(flip_map)
        flip_map[fid] = (not go_forward)  # flipped if we walk tail→head
        _emit_vaults_until.flip = flip_map

        # Emit at the entry point first
        _emit_vaults_until(fid, 0.0, v_state)

        # Walk vertices in oriented order
        if go_forward:
            vertex_iter = range(enter_idx + 1, len(nodes))
            tail_index = len(nodes) - 1
        else:
            vertex_iter = range(enter_idx - 1, -1, -1)
            tail_index = 0

        prev_xy = coords[enter_idx]
        for vi in vertex_iter:
            cur_xy = coords[vi]
            oriented_along = along[vi] if go_forward else (along[-1] - along[vi])
            _emit_vaults_until(fid, oriented_along, v_state)

            nk = nodes[vi]
            touching = [tfid for tfid in node_to_features.get(nk, []) if tfid != fid and tfid not in visited_fids]
            if touching:
                ends_here = (vi == tail_index)
                nbrs = list(touching)
                if not ends_here:
                    # Continuing feature → exclude same-cable here (detours only)
                    nbrs = [nb for nb in nbrs if not _same_cable(fid, nb)]

                neighbors = sorted(
                    nbrs,
                    key=lambda nb: _neighbor_rank(
                        fid,
                        prev_xy,
                        cur_xy,
                        nb,
                        is_start_node=False,
                        prefer_same_line_first=ends_here,
                    ),
                )

                # DFS into neighbors (keeps B>A & the rest of the logic)
                for nb in neighbors:
                    _walk_feature(nb, nk, cur_xy, v_state)

            prev_xy = cur_xy

        # End of feature: flush remaining vaults
        _emit_remaining_vaults(fid, v_state)

    # ---------- T3: Section-first start selection (simple) ----------
    def _find_start_node_from_t3(
        t3_gj_in: Dict[str, Any],
        feature_nodes_in: Dict[int, List[Tuple[int, int]]],
    ) -> Optional[Tuple[int, int]]:
        t3_pts: List[List[float]] = []
        for f in t3_gj_in.get("features", []):
            g = f.get("geometry") or {}
            if g.get("type") == "Point" and isinstance(g.get("coordinates"), (list, tuple)) and len(g["coordinates"]) == 2:
                t3_pts.append(g["coordinates"])
        if not t3_pts:
            return None

        all_nodes_list: List[Tuple[int, int]] = []
        for nds in feature_nodes_in.values():
            all_nodes_list.extend(nds)
        if not all_nodes_list:
            return None
        all_nodes_list = list(dict.fromkeys(all_nodes_list))  # de-dup

        t3_xy = t3_pts[0]
        best = (float("inf"), None)
        for nk in all_nodes_list:
            d = distance_feet(t3_xy, _nk_to_xy(nk))
            if d < best[0]:
                best = (d, nk)
        return best[1]

    start_nk = _find_start_node_from_t3(t3_gj, feature_nodes)

    # Build initial starts strictly by Section number at T3: 1, then 2, then 3, ...
    start_fids_ordered: List[int] = []
    if start_nk is not None:
        touching = list(node_to_features.get(start_nk, []))
        if touching:
            # Map touching → (section_number, fid)
            touch_pairs = [(_sec_num(_sec_norm(feature_props[fid].get("Section"))), fid) for fid in touching]
            # Group by section number
            from collections import defaultdict
            sec_to_fids: Dict[int, List[int]] = defaultdict(list)
            for sn, fid in touch_pairs:
                sec_to_fids[sn].append(fid)
            # Build ordered list: smallest section number to largest
            for sn in sorted(sec_to_fids.keys()):
                # Deterministic order within a Section at T3:
                # **HIGHER Distribution first** (72 > 48 > 24), then letter 'A' first, then alphabetical; then ID.
                def _t3_key(fid: int) -> Tuple[int, Tuple[int, str], str]:
                    p = feature_props[fid]
                    dist = _dist_num(p.get("Distribution Type"))  # 24, 48, 72...
                    letter = (p.get("Fiber Letter #1") or "").strip().upper()

                    # At T3: always prefer 'A' first (no special B>A for 24-count here)
                    letter_rank = (0 if letter == "A" else 1, letter)

                    id_key = str(p.get("ID") or "")
                    # NOTE: -dist => HIGHER distribution first (72→48→24)
                    return (-dist, letter_rank, id_key)


                start_fids_ordered.extend(sorted(sec_to_fids[sn], key=_t3_key))

    # Walk all chosen starts Section-by-Section (starting set)
    v_state = {}
    for sfid in start_fids_ordered:
        if sfid in visited_fids:
            continue
        _walk_feature(sfid, start_nk, None, v_state)

    # Finish remaining features by global Section order (1..n then None)
    remaining_secs = sorted(by_section.keys(), key=_sec_num)
    for sec in remaining_secs:
        while True:
            candidates = [f for f in (by_section.get(sec, []) or []) if f not in visited_fids]
            if not candidates:
                break
            fid0 = candidates[0]
            nk0 = feature_nodes[fid0][0]
            _walk_feature(fid0, nk0, None, v_state)

    return result
