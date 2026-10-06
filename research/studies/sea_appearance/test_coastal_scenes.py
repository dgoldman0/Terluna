"""Admission and geometry checks for the two selected cloud scenarios."""
import json
from pathlib import Path
import numpy as np
import pytest
from research.studies.sea_appearance import nobili
from research.studies.sea_appearance.cloud_candidates import matching_date
from geography.lunar_ephemeris import earth_and_sun

HERE = Path(__file__).parent


def test_native_land_rejected_before_lighting_or_weather(monkeypatch, tmp_path):
    manifest = tmp_path/'inputs.json'
    manifest.write_text(json.dumps(dict(site_lon_lat_deg=[148.3811,0.], horizons=dict(jd_tt=2463761.),
        cloud=dict(snapshot=300,observer_column=1))))
    monkeypatch.setattr(nobili, 'terrain', lambda *args, **kwargs: dict(water_depth_at_observer_m=-3151.))
    with pytest.raises(ValueError, match='viewpoint is land'):
        nobili.build(tmp_path/'build', tmp_path/'missing-climate', manifest)


@pytest.mark.parametrize('name', ['fecunditatis-earth-clouds','smythii-twilight-clouds'])
def test_selected_coasts_have_water_and_match_weather_solar_phase(name):
    p = json.loads((HERE/'results'/f'{name}.json').read_text())
    assert p['terrain']['water_depth_at_observer_m'] > 100
    assert p['cloud']['model_sun_longitude_deg'] == pytest.approx(p['cloud']['ephemeris_sun_longitude_deg'],abs=.002)
    target = p['cloud']['model_sun_longitude_deg']
    jd = matching_date(target, p['jd_tt']+.1)
    s = earth_and_sun(np.array([jd]))['sun'][0]
    lon = np.degrees(np.arctan2(s[1],s[0]))%360
    assert abs((lon-target+180)%360-180) < 1e-6
    assert jd == pytest.approx(p['jd_tt'],abs=1e-7)


def test_earth_visible_and_absent_are_distinct_orbital_cases():
    first = json.loads((HERE/'results/fecunditatis-earth-clouds.json').read_text())
    second = json.loads((HERE/'results/smythii-twilight-clouds.json').read_text())
    assert first['earth']['lit_fraction'] > .99
    assert 41 < first['earth']['elevation_deg'] < 44
    assert 1.9 < first['earth']['diameter_deg'] < 2
    assert first['light']['sun_elevation_deg'] < -50
    assert second['earth']['elevation_deg']+second['earth']['diameter_deg']/2 < -6
    assert -4 < second['light']['sun_elevation_deg'] < -3
    for manifest in ['fecunditatis_earth_clouds_inputs.json','smythii_twilight_clouds_inputs.json']:
        p = json.loads((HERE/manifest).read_text())
        assert p['selection']['overhead_tau'] < 1


def test_fetch_products_are_bound_to_the_scene_and_reduce_unlimited_waves():
    for scene, waves in [('fecunditatis-earth-clouds','fecunditatis-cloud-waves'),
                         ('smythii-twilight-clouds','smythii-twilight-waves')]:
        scene_path = HERE/'results'/f'{scene}.json'
        p = json.loads(scene_path.read_text())
        w = json.loads((HERE/'results'/f'{waves}.json').read_text())
        assert w['scene_sha256'] == nobili.digest(scene_path)
        assert w['significant_height_m'] <= p['wave_sensitivity'][0]['significant_height_m']
        assert 40000 < w['fetch_m'] < 300000
        assert w['fetch_source']['first_land_height_m'] >= 0
        assert .84 <= w['selected_inverse_wave_age'] <= 5
