"""Checks of the rendering coasts' local terrain grids."""
import math

import numpy as np
import pytest

from geography import coast_terrain
from shared.constants import MOON_RADIUS


def test_local_grid_points_lie_where_the_projection_puts_them():
    lat, lon = coast_terrain.destination(163.2, -33.7, np.array(0.0), np.array(0.0))
    assert (float(lat), float(lon)) == pytest.approx((-33.7, 163.2))
    arc = MOON_RADIUS * math.radians(1.0)
    lat, lon = coast_terrain.destination(163.2, -33.7, np.array(0.0), np.array(arc))
    assert float(lat) == pytest.approx(-32.7, abs=1e-9) and float(lon) == pytest.approx(163.2, abs=1e-9)
    lat, lon = coast_terrain.destination(10.0, 0.0, np.array(arc), np.array(0.0))
    assert float(lat) == pytest.approx(0.0, abs=1e-9) and float(lon) == pytest.approx(11.0, abs=1e-9)


def test_the_manifest_pins_every_slice_it_lists():
    manifest = coast_terrain.manifest()
    for entry in manifest["slices"]:
        first, last = entry["selected_byte_range"]
        assert last - first + 1 == entry["selected_bytes"] == entry["row_count"] * entry["row_bytes"]
        assert first == entry["first_row"] * entry["row_bytes"]
        assert len(entry["selected_sha256"]) == 64 and entry["lines"] * entry["row_bytes"] == entry["source_size_bytes"]
        # Pixel centres: latitude (LINE_PROJECTION_OFFSET - row) / 256 spans the coast's band of rows.
        top = (entry["line_projection_offset"] - entry["first_row"]) / entry["pixels_per_degree"]
        bottom = (entry["line_projection_offset"] - entry["first_row"] - entry["row_count"] + 1) / entry["pixels_per_degree"]
        assert top - bottom == pytest.approx(2 * coast_terrain.HALF_BAND_DEG, abs=0.01)
