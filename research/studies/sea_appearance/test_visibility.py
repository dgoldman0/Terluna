"""Checks of angular units, conserved light and the texture counterfactual."""
import zipfile

import numpy as np
import pytest

from research.studies.sea_appearance import visibility as v


def test_grid_solid_angle_matches_rectangular_aperture():
    ppd, field = 30, 8
    _, _, omega, _ = v.grid(ppd, field)
    # Independent exact solid angle of a centred rectangular plane aperture.
    a = round(ppd * field) / (2 * ppd * 180 / np.pi)
    expected = 4 * np.arctan(a*a / np.sqrt(1 + 2*a*a))
    assert omega.sum() == pytest.approx(expected, rel=1e-7)
    _, _, finer, _ = v.grid(2*ppd, field)
    assert abs(finer.sum()-expected) < .26*abs(omega.sum()-expected)


def test_counterfactual_preserves_light_and_limb_but_removes_texture():
    x, y, omega, cosine = v.grid(100, 4)
    r = np.hypot(x, y) / np.tan(np.radians(.9))
    pattern = np.where(r<1, 1 + .35*np.sin(x*600)*np.cos(y*450), 0)
    direct = v.normalize_beam(pattern, omega, cosine, 2.56)
    reference = direct + .53
    counter, roi = v.featureless_inner(reference, r, omega, cosine)
    assert np.array_equal(counter[r>=.9], reference[r>=.9])
    assert np.sum((counter-.53)*omega*cosine) == pytest.approx(2.56, rel=1e-12)
    assert np.sum(direct*omega*cosine) == pytest.approx(2.56, rel=1e-12)
    assert np.std(counter[roi]) < .05*np.std(reference[roi])
    assert counter.min() > 0


@pytest.mark.parametrize('value', [0, -1, np.nan, np.inf])
def test_invalid_light_cannot_be_normalized(value):
    with pytest.raises(ValueError):
        v.normalize_beam(np.array([value]), np.ones(1), np.ones(1), 2.56)


def test_pds_byte_offset_and_little_endian_sample_contract(tmp_path):
    header = ('RECORD_BYTES = 512\r\n^IMAGE = 2\r\nLINES = 2\r\n'
              'LINE_SAMPLES = 3\r\nSAMPLE_TYPE = PC_REAL\r\nSAMPLE_BITS = 32\r\n'
              'MAP_RESOLUTION = 76 <PIX/DEG>\r\nLINE_PROJECTION_OFFSET = 5319.5\r\n'
              'SAMPLE_PROJECTION_OFFSET = -0.5\r\nWESTERNMOST_LONGITUDE = 0\r\n'
              'PRODUCT_VERSION_ID = 1.2\r\nEND\r\n').encode()
    values = np.array([[.125, .25, .5], [1, 2, 4]], dtype='<f4')
    path = tmp_path/'control.IMG'
    path.write_bytes(header.ljust(512, b' ')+values.tobytes())
    actual, meta = v.pds_tile(path)
    np.testing.assert_array_equal(actual, values)
    assert meta['sample_offset'] == -.5
    assert meta['ppd'] == 76


def test_external_model_hashes_reject_modified_code(tmp_path):
    archive=tmp_path/'model.zip'; directory=tmp_path/'model'; directory.mkdir()
    with zipfile.ZipFile(archive,'w') as z:
        z.writestr('model/run.m','disp(1);')
    (directory/'run.m').write_text('disp(1);')
    assert v.verify_model(archive,directory)['file_count']==1
    (directory/'run.m').write_text('disp(2);')
    with pytest.raises(ValueError,match='mismatch'):
        v.verify_model(archive,directory)
