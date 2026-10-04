"""Checks of the shore month's pieces: the GCM's Sun clock, PlaSim's records and the month matching."""
import struct

import numpy as np

from climate.waves import shore_month as sm

CLOCK = dict(step_s=1800.0, sidereal_day_s=27.321661 * 86400, sidereal_year_s=365.25636 * 86400,
             first_snapshot_step=0.0)


def test_sun_starts_over_the_far_side_and_moves_west_once_per_synodic_month():
    synodic_s = 1 / (1 / CLOCK["sidereal_day_s"] - 1 / CLOCK["sidereal_year_s"])
    steps = np.array([0.0, 0.25, 0.5, 1.0]) * synodic_s / CLOCK["step_s"]
    assert np.allclose(sm.wrap(sm.subsolar_at_step(steps, CLOCK) - [180.0, 90.0, 0.0, 180.0]), 0, atol=1e-6)
    hours = np.array([0.0, 0.25 * synodic_s / 3600])
    assert np.allclose(sm.subsolar_at_hour(hours, CLOCK), [180.0, 90.0], atol=1e-6)


def test_plasim_records_are_read_by_name(tmp_path):
    def record(payload):
        return struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))
    path = tmp_path / "plasim_status"
    path.write_bytes(record(b"nstep".ljust(16)) + record(struct.pack("<i", 1468662))
                     + record(b"naccuout".ljust(16)) + record(struct.pack("<i", 26)))
    assert struct.unpack("<i", sm.plasim_value(path, "nstep"))[0] == 1468662
    assert struct.unpack("<i", sm.plasim_value(path, "naccuout"))[0] == 26


def test_matching_times_find_each_pass_of_the_sun():
    jd = np.arange(0.0, 90.0, 0.125)
    sun = sm.wrap(100.0 - 360.0 * jd / 29.5)                # passes 40 E at 29.5 * 60 / 360 days, once a month
    passes = sm.matching_times(jd, sun, 40.0)
    assert np.allclose(passes, 29.5 * 60 / 360 + 29.5 * np.arange(3), atol=1e-9)
