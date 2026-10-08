"""Checks of the sea-slope product: its provenance, its parts and how they add."""
import hashlib
import json

import numpy as np

from research.studies.sea_appearance import slopes as S
from shared.provenance import constants_changed


def test_the_committed_slopes_are_current():
    product = json.loads(S.RESULT.read_text())
    assert product["schema"] == "terluna.research.sea-slopes/1"
    for path, expected in product["producer"]["files"].items():
        assert hashlib.sha256((S.ROOT / path).read_bytes()).hexdigest() == expected, path
    assert not constants_changed(product["producer"]["constants"])
    assert product["series"]["sha256"] == hashlib.sha256(S.SERIES.read_bytes()).hexdigest()
    assert product["inputs"]["research/studies/sea_appearance/results/lighting_calendar.json"] == \
        hashlib.sha256(S.CALENDAR.read_bytes()).hexdigest()


def test_each_coast_adds_resolved_and_short_waves_at_the_hours_with_spectra():
    product = json.loads(S.RESULT.read_text())
    with np.load(S.SERIES) as z:
        hours = z["hours"]
        for coast in product["coasts"]:
            key = coast["short"].replace(" ", "_").lower()
            short = z[f"{key}/short_covariance"].reshape(-1, 2, 2)
            assert short.shape[0] == len(hours)
            assert np.all(np.linalg.eigvalsh(short.astype(float)) >= -1e-9)
            glassy = z[f"{key}/short_level"] == 0
            assert glassy.sum() == coast["hours_without_short_waves"]
            assert np.all(short[glassy] == 0)
            assert np.all(z[f"{key}/u_star"][glassy] <= product["short_wave_threshold_u_star_m_s"] + 1e-6)
            resolved = z[f"{key}/resolved_covariance"].reshape(-1, 2, 2)
            assert len(resolved) == coast["resolved_samples"]
            if len(resolved):
                at = np.searchsorted(hours, z[f"{key}/resolved_hours"])
                total = np.trace(resolved + short[at], axis1=1, axis2=2)
                assert abs(np.median(total) - coast["total_mss"]["p50"]) < 1e-5


def test_cox_and_munk_reference():
    assert S.cox_munk_clean(0.0) == 0.003
    assert abs(S.cox_munk_clean(10.0) - (0.003 + 5.12e-3 * 10.0 * 1.024)) < 0.002


def test_rotation_keeps_the_principal_slopes():
    cov = S.rotated(0.02, 0.01, 30.0)
    assert np.allclose(np.linalg.eigvalsh(cov), [0.01, 0.02])
    direction = np.array([np.cos(np.radians(30)), np.sin(np.radians(30))])
    assert np.isclose(direction @ cov @ direction, 0.02)
