# modules/drops/lshape.py
from __future__ import annotations
from typing import List, Sequence

def build_l_drop_linestring(a: Sequence[float], b: Sequence[float]) -> List[List[float]]:
    """Return an L-shaped LineString (lon,lat) path from point a -> b: a -> (b.x, a.y) -> b"""
    ax, ay = a
    bx, by = b
    elbow = [bx, ay]
    return [list(a), elbow, list(b)]
