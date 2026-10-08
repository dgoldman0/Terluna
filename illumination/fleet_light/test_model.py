"""The fleet's night light: geometry checks and the stored product's provenance."""
import hashlib
import json
from pathlib import Path

import numpy as np

from illumination.fleet_light import model

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRODUCT = json.loads((HERE / 'results' / 'night_light.json').read_text())


def test_a_tile_facing_the_moon_returns_light_at_its_own_impact_parameter():
    P = np.array([18000.0, 900.0, 3000.0])
    n = P / np.linalg.norm(P)
    ray = np.array([1.0, 0, 0]) - 2.0 * (np.array([1.0, 0, 0]) @ n) * n
    t = -P @ ray
    assert abs(np.linalg.norm(P + t * ray) - np.hypot(P[1], P[2])) < 1e-6


def test_untilted_glints_land_only_near_the_terminator():
    lux = PRODUCT['specular']['random_0']['lux_by_alpha_band']['lux']
    assert sum(lux[:6]) == 0 and sum(lux[6:]) > 0


def test_rolling_the_night_half_dims_the_glow_and_turns_the_glints_away():
    rolls = PRODUCT['night_roll']
    mid = [rolls[r]['diffuse_lux_at_equator'][0] for r in ('0', '30', '60', '80')]
    assert mid == sorted(mid, reverse=True) and mid[-1] < mid[0] / 10
    assert all(rolls[r]['specular_mean_lux_over_night_hemisphere'] == 0 for r in ('30', '60', '80'))


def test_product_binds_its_producer_and_inputs():
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    for name, h in PRODUCT['producer']['files'].items():
        assert digest(HERE.parent / name) == h, name
    for name, h in PRODUCT['producer']['inputs'].items():
        assert digest(ROOT / name) == h, name
