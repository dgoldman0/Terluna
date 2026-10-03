"""Native tile registration, water connectivity and pinned coastal bytes."""
import hashlib
import json

import numpy as np
import pytest

from geography import shore_terrain


def test_regional_tile_reads_north_first_rows_into_north_increasing_grid():
    source = dict(pixels_per_degree=4,first_row=344,row_count=16,samples=720)
    lon, lat, rows, cols = shore_terrain.native_indices((1.125,1.125,.5,.5),source)
    np.testing.assert_array_equal(lon,[1.125,1.375,1.625])
    np.testing.assert_array_equal(lat,[1.125,1.375,1.625])
    np.testing.assert_array_equal(rows,[11,10,9])
    np.testing.assert_array_equal(cols,[4,5,6])
    for bounds in ((1.,1.125,.5,.5),(1.125,-.125,.5,.5),(1.125,1.125,.3,.5)):
        with pytest.raises(ValueError):
            shore_terrain.native_indices(bounds,source)


def test_diagonal_contacts_and_isolated_pools_keep_their_water_identity():
    depth = np.array([[10.,5.,-2.,-2.],[10.,-1.,3.,3.],[10.,-1.,3.,3.]])
    wet, bottom = shore_terrain.flooded_component(depth)
    np.testing.assert_array_equal(wet,[[True,True,False,False],[True,False,False,False],[True,False,False,False]])
    np.testing.assert_array_equal(bottom[wet],depth[wet])
    assert (bottom[~wet] <= -1).all()


def test_changed_terrain_bytes_are_rejected_before_processing(tmp_path,monkeypatch):
    label, data = tmp_path/"label.lbl",tmp_path/"slice.img"
    label.write_bytes(b'label'); data.write_bytes(b'good')
    manifest = tmp_path/"source.json"
    manifest.write_text(json.dumps(dict(local_label="label.lbl",local_image="slice.img",
        label_sha256=hashlib.sha256(b'label').hexdigest(),selected_sha256=hashlib.sha256(b'good').hexdigest(),
        label_bytes=5,selected_bytes=4)))
    monkeypatch.setattr(shore_terrain,"ROOT",tmp_path)
    monkeypatch.setattr(shore_terrain,"MANIFEST",manifest)
    shore_terrain.restore()
    data.write_bytes(b'evil')
    with pytest.raises(ValueError,match="Pinned coastal input differs"):
        shore_terrain.restore()
