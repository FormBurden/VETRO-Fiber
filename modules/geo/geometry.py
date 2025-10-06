# modules/geo/geometry.py
# Small, dependency-free helpers for feet distances using lat/lon.
from __future__ import annotations
import math
from typing import Sequence, Tuple

FEET_PER_METER = 3.28084

def _meters_per_deg_lat(lat_deg: float) -> float:
    # Accurate enough for engineering-scale distances
    lat = math.radians(lat_deg)
    return (111132.92 - 559.82 * math.cos(2 * lat) + 1.175 * math.cos(4 * lat) - 0.0023 * math.cos(6 * lat))

def _meters_per_deg_lon(lat_deg: float) -> float:
    lat = math.radians(lat_deg)
    return (111412.84 * math.cos(lat) - 93.5 * math.cos(3 * lat) + 0.118 * math.cos(5 * lat))

def _dx_dy_meters(p1: Sequence[float], p2: Sequence[float]) -> Tuple[float, float]:
    lon1, lat1 = p1
    lon2, lat2 = p2
    lat_mid = (lat1 + lat2) / 2.0
    dx_m = (lon2 - lon1) * _meters_per_deg_lon(lat_mid)
    dy_m = (lat2 - lat1) * _meters_per_deg_lat(lat_mid)
    return dx_m, dy_m

def distance_feet(p1: Sequence[float], p2: Sequence[float]) -> float:
    dx_m, dy_m = _dx_dy_meters(p1, p2)
    return math.hypot(dx_m, dy_m) * FEET_PER_METER

def l_distance_feet(p1: Sequence[float], p2: Sequence[float]) -> float:
    dx_m, dy_m = _dx_dy_meters(p1, p2)
    return (abs(dx_m) + abs(dy_m)) * FEET_PER_METER

def point_segment_distance_feet(p: Sequence[float], a: Sequence[float], b: Sequence[float]) -> float:
    # Work in local meters for numerical stability
    dx_m, dy_m = _dx_dy_meters(a, b)
    px_m, py_m = _dx_dy_meters(a, p)
    seg2 = dx_m * dx_m + dy_m * dy_m
    if seg2 == 0:
        # a==b
        return math.hypot(px_m, py_m) * FEET_PER_METER
    t = max(0.0, min(1.0, (px_m * dx_m + py_m * dy_m) / seg2))
    proj_x = t * dx_m
    proj_y = t * dy_m
    dist_m = math.hypot(px_m - proj_x, py_m - proj_y)
    return dist_m * FEET_PER_METER

def polyline_length_feet(coords: Sequence[Sequence[float]]) -> float:
    total = 0.0
    for i in range(1, len(coords)):
        total += distance_feet(coords[i-1], coords[i])
    return total

def project_point_onto_polyline_feet(p: Sequence[float], line: Sequence[Sequence[float]]) -> Tuple[float, Tuple[float, float]]:
    """
    Returns (distance_along_feet, closest_point_lonlat).
    """
    best = (float("inf"), 0.0, (line[0][0], line[0][1]))  # (dist_ft, dist_along_ft, point)
    dist_along = 0.0
    for i in range(1, len(line)):
        a, b = line[i-1], line[i]
        # Distance from p to segment a-b
        d_ft = point_segment_distance_feet(p, a, b)

        # Approximate projection along length: move to the closer vertex if exact footpoint not computed
        # For our grouping/order use-case, picking the closer vertex is adequate.
        da = distance_feet(a, p)
        db = distance_feet(b, p)
        along_here = dist_along + (0.0 if da <= db else distance_feet(a, b))

        if d_ft < best[0]:
            best = (d_ft, along_here, (a[0], a[1]) if da <= db else (b[0], b[1]))
        dist_along += distance_feet(a, b)
    return best[1], best[2]
