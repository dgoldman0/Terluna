"""People leaving a fire: the hydraulic model of movement (Nelson and Mowrer, as in the SFPE Handbook), the required
safe egress time, and the throughput of an evacuation by air.

In the hydraulic model a crowd's walking speed falls with its density, S = k - a k D (a = 0.266 m2 per person; k =
1.40 m/s in corridors, aisles, ramps and doorways), with a free speed of 0.85 k below 0.54 persons per m2, and the
flow through a door or corridor peaks at 1.3 persons a second per metre of effective width (the clear width less
15 cm at each door jamb).

At lunar gravity people walk as fast as on Earth: the preferred change from walking to running comes at 1.42 m/s
in lunar gravity on parabolic flights (De Witt et al. 2014). What changes is grip: a foot can push sideways or brake
with about a sixth of Earth's force, so people take six times as long to stop, start and turn, and a dense crowd's
stop-and-go flow through doors is likely slower. There are no measurements; `flow_factor` carries an assumed range.
"""
from __future__ import annotations
import numpy as np

SPEED_CONSTANT_M_S = 1.40
DENSITY_COEFFICIENT_M2 = 0.266
FREE_DENSITY = 0.54
MAX_SPECIFIC_FLOW = 1.3        # persons per second per metre of effective width
DOOR_BOUNDARY_M = 0.15         # lost at each jamb


def speed_m_s(density_per_m2, k: float = SPEED_CONSTANT_M_S, factor: float = 1.0):
    """Walking speed in a crowd of this density: k - a k D, or the free speed 0.85 k below 0.54 persons/m2."""
    d = np.asarray(density_per_m2, float)
    s = np.where(d < FREE_DENSITY, 0.85 * k, k - DENSITY_COEFFICIENT_M2 * k * d)
    return factor * np.clip(s, 0.0, None)


def specific_flow(density_per_m2, k: float = SPEED_CONSTANT_M_S, factor: float = 1.0):
    """Persons per second per metre of effective width: speed times density."""
    return speed_m_s(density_per_m2, k, factor) * np.asarray(density_per_m2, float)


def effective_width_m(door_widths_m) -> float:
    return float(sum(max(w - 2.0 * DOOR_BOUNDARY_M, 0.0) for w in door_widths_m))


def movement_s(people: float, travel_m: float, door_widths_m, flow_factor: float = 1.0,
               free_speed_m_s: float = 0.85 * SPEED_CONSTANT_M_S) -> dict:
    """Time for everyone to leave a space: the longer of walking the longest route and the queue at the doors."""
    walk = travel_m / (free_speed_m_s * flow_factor)
    queue = people / (effective_width_m(door_widths_m) * MAX_SPECIFIC_FLOW * flow_factor)
    return dict(walk_s=walk, queue_s=queue, movement_s=max(walk, queue))


def rset_s(detection_s: float, alarm_s: float, premovement_s: float, movement: float) -> float:
    """Required safe egress time: detection, alarm, the time people take to start moving, and the movement."""
    return detection_s + alarm_s + premovement_s + movement


def aerial_evacuation(people: float, rim_m: float, berth_spacing_m: float, ship_capacity: float,
                      cycle_s: float) -> dict:
    """Clearing a floor from its rim by air: berths along the rim, each loading a ship of the given capacity once a
    cycle (docking, loading, casting off and the next ship coming in)."""
    berths = int(rim_m // berth_spacing_m)
    rate = berths * ship_capacity / cycle_s
    return dict(berths=berths, persons_per_s=rate, clear_s=people / rate, ship_trips=int(np.ceil(people / ship_capacity)))
