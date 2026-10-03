"""Check native terrain orientation and exclusion of a disconnected lake."""
import hashlib
import json

import numpy as np

from geography import coastal_grid


def test_coastal_export_preserves_pixel_centres_and_connected_water(tmp_path, monkeypatch):
    terrain = tmp_path/"ldem_16.img"
    source = np.memmap(terrain, mode="w+", dtype="<i2", shape=(2880, 5760))
    # Five rows north-to-south at 0.03125..0.28125 E and -0.03125..-0.28125 N.
    heights = np.array([[0, 0, 30, 30, 30],
                        [1, 2, 30, 30, 30],
                        [3, 4, 30, 30, 30],
                        [5, 6, 30, 30, 0],
                        [7, 8, 30, 30, 30]])
    source[1440:1445, :5] = heights*2
    source.flush()
    del source
    gravity = tmp_path/"jggrx_0420a_sha.tab"
    gravity.write_text("synthetic geoid replaced below\n")
    files = {p.name: p for p in (terrain, gravity)}
    manifest = dict(files=[dict(name=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                               bytes=p.stat().st_size, credit="synthetic test") for p in files.values()])
    monkeypatch.setattr(coastal_grid.fetch_inputs, "manifest", lambda: manifest)
    monkeypatch.setattr(coastal_grid.fetch_inputs, "path", lambda name: files[name])
    monkeypatch.setattr(coastal_grid, "geoid", lambda grid, max_degree: np.zeros((5, 5)))
    monkeypatch.setattr(coastal_grid, "ROOT", tmp_path)
    (tmp_path/"geography/products").mkdir(parents=True)
    (tmp_path/"geography/topography.py").write_text("synthetic test\n")
    (tmp_path/"geography/coastal_grid.py").write_text("synthetic test\n")
    monkeypatch.setattr(coastal_grid, "__file__", str(tmp_path/"geography/coastal_grid.py"))
    (tmp_path/"shared").mkdir()
    (tmp_path/"shared/constants.json").write_text("{}\n")
    np.savez(tmp_path/"geography/products/atlas_28pct_4ppd.npz",
             metadata=json.dumps(dict(schema="terluna.geography.atlas-grid/1", sea_level_m=10)))
    product = coastal_grid.sector(west=.03125, south=-.28125, width=.25, height=.25)
    np.testing.assert_array_equal(product["height_m"], heights[::-1])
    assert product["longitude_deg"].tolist() == [.03125, .09375, .15625, .21875, .28125]
    assert product["latitude_deg"].tolist() == [-.28125, -.21875, -.15625, -.09375, -.03125]
    assert product["wet"].sum() == 10
    np.testing.assert_array_equal(product["depth_m"][:, :2], 10-heights[::-1, :2])
    assert product["depth_m"][1, 4] == -1
    assert json.loads(product["metadata"])["excluded_wet_nodes"] == 1
