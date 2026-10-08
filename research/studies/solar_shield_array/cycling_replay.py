"""Replay the published return arc without ignored optimizer files."""
import json

import numpy as np

from protection.dynamics.ephemeris import Ephemeris
from protection.dynamics.cycling import CyclingEnvironment
from .cycling_search import HERE, CENTRAL, SCENARIO, ReturnProblem
from .cycling_analysis import digest, ROOT
from shared import constants as K
from shared.provenance import constants_changed


def replay():
    product = json.loads((HERE/"results/cycling.json").read_text())
    for name, expected in product["producer"]["source_hashes"].items():
        if digest(ROOT/name) != expected:
            raise ValueError(f"Stale source: {name}")
    if constants_changed(product["producer"]["constants"]):
        raise ValueError("Stale constants")
    arc = product["independently_replayed_return_arc"]
    t = np.array(arc["control_times_s"])
    eph = Ephemeris(epoch=product["epoch_tdb"])
    samples = eph.sample(np.arange(t[0]-3600, t[-1]+7200, 1800.))
    eph.close()
    env = CyclingEnvironment(samples, product["areal_mass_kg_m2"], CENTRAL)
    problem = ReturnProblem(env, arc["radius_seed_km"]*1000, t[-1]-t[0], t[0], len(t)-1)
    nodes = np.zeros((len(t), 8))
    nodes[0, :6] = np.array(arc["initial_icrf_state_m_m_s"])/problem.scale
    nodes[:, 6:] = arc["control_angles_rad"]
    trace = problem.propagate(nodes.ravel(), step=600., rtol=2e-11)
    position = float(np.linalg.norm(trace["local_closure"][:3]))
    velocity = float(np.linalg.norm(trace["local_closure"][3:]))
    if abs(position-arc["closure_position_m"]) > 1 or abs(velocity-arc["closure_velocity_m_s"]) > 1e-5:
        raise ValueError("Published return arc did not reproduce to 1 m / 10 micrometres/s")
    return dict(days=(t[-1]-t[0])/K.JULIAN_DAY, closure_position_m=position,
                closure_velocity_m_s=velocity, replay_pass=True)


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2))
