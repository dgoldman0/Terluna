"""Mie efficiencies against their limits and Bohren and Huffman's worked example; the haze product's provenance."""
import hashlib
import json
from pathlib import Path

import numpy as np

from illumination.aerosol import mie

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def test_rayleigh_limit():
    for m in (1.5 + 0j, 1.33 + 0.01j):
        x = 0.01
        _, qsca, g = mie.efficiencies(x, m)
        assert abs(qsca / (8 / 3 * x ** 4 * abs((m ** 2 - 1) / (m ** 2 + 2)) ** 2) - 1) < 1e-3
        assert abs(g) < 1e-3


def test_large_spheres_extinguish_twice_their_area():
    assert abs(mie.efficiencies(1000.0, 1.5 + 0j)[0] - 2.0) < 0.03


def test_bohren_and_huffman_example():
    qext, qsca, _ = mie.efficiencies(5.213, 1.55 + 0j)
    assert abs(qext - 3.1054) < 2e-3 and abs(qsca - qext) < 1e-9


def test_haze_product_binds_its_producer_and_inputs():
    path = HERE / 'results' / 'haze.json'
    d = json.loads(path.read_text())
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    for name, h in d['producer']['files'].items():
        assert digest(HERE.parent / name) == h, name
    for name, h in d['producer']['inputs'].items():
        assert digest(ROOT / name) == h, name
    assert 0.6 < d['molecular_tau_550'] < 0.8
    for region in d['regions'].values():
        for case in region['cases'].values():
            assert case['day']['visibility_km'] <= d['clear_visibility_km'] + 1
