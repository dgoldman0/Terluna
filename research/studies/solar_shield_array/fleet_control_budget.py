"""Ideal lateral photon-control scales, separate from fleet demonstrations."""
import json
import math
from pathlib import Path

from shared import constants as K
from shared.provenance import constants_used
from .cycling_search import CENTRAL, HERE, ROOT
from .fleet_run import digest


def main():
    sigma = .05
    face_on = 2*CENTRAL*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*sigma)
    transverse = face_on*2/(3*math.sqrt(3))
    speed = math.sqrt(K.MOON_GM/15e6)
    sources = {str(Path(__file__).relative_to(ROOT)): digest(__file__),
        'research/studies/solar_shield_array/scenario.json': digest(HERE/'scenario.json'),
        'research/studies/solar_shield_array/cycling_search.py': digest(HERE/'cycling_search.py')}
    result = dict(schema='terluna.research.solar-shield-fleet-control-budget/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources)),
        evidence='Analytical 1-AU ideal specular filter limits; no avoidance maneuver or controlled fleet was propagated',
        areal_mass_kg_m2=sigma, diverted_bolometric_fraction=CENTRAL,
        face_on_photon_acceleration_m_s2=face_on,
        maximum_transverse_photon_acceleration_m_s2=transverse,
        maximizing_normal_cone_deg=math.degrees(math.asin(1/math.sqrt(3))),
        rest_to_rest_translation_hours={str(distance): 2*math.sqrt(distance/transverse)/3600
            for distance in [50., 100., 10000.]},
        lunar_circular_speed_at_15000km_m_s=speed,
        along_track_phase_error_for_50m_s=50/speed,
        reading_rule='Translation times are the ideal symmetric bang-bang scale for a bounded one-axis acceleration, zero initial/final lateral velocity and no gravitational coupling. Solar availability, Earth-safe beams, spectral/coverage constraints and changing orbits must enter an actual maneuver. The phase-error figure is a state-estimation tolerance, not a commanded update interval.',
        excluded_loads=['collectors', 'habitats', 'attitude and control hardware beyond the assumed base optical allocation'],
        delivered_electrical_power_W=None)
    (HERE/'results/fleet_control_budget.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
