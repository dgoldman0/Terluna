"""The selected view must keep terrain, orbital geometry and source admission consistent."""
import json
from pathlib import Path
import numpy as np
import pytest
from research.studies.sea_appearance import nobili


def test_nobili_stated_geometry_clears_the_native_rim():
    p=json.loads((Path(nobili.__file__).with_name('results')/'nobili-night.json').read_text())
    earth=p['earth'];terrain=p['terrain']['profile']
    rim=np.interp(earth['azimuth_deg'],terrain['azimuth_deg'],terrain['elevation_deg'])
    assert earth['elevation_deg']-earth['diameter_deg']/2>rim+2
    assert 1.9<earth['diameter_deg']<2.05
    assert .63<earth['lit_fraction']<.68
    assert p['light']['sun_elevation_deg'] < -85
    assert p['terrain']['water_depth_at_observer_m'] > 100
    assert p['cloud']['model_sun_longitude_deg']==pytest.approx(p['cloud']['ephemeris_sun_longitude_deg'],abs=.002)


def test_input_admission_rejects_changed_bytes(tmp_path):
    f=tmp_path/'input';f.write_bytes(b'original');expected=nobili.digest(f)
    nobili.checked(f,expected);f.write_bytes(b'changed')
    with pytest.raises(ValueError,match='hash mismatch'):nobili.checked(f,expected)
