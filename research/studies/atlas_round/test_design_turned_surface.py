"""Checks of the climate model's surface for the run on the design surface: the runner's experiment A28D takes only its
own surface; once built, the settled turn under the design water is lighter-loaded and about a degree from the final
water's, the record holds that turn's frame and the design routing's water kept through the turn with less water than
the final surface (the parted sea's floor left dry), files named by the routing's content, and the runner takes the
record."""
import hashlib
import json

import pytest

from climate.gcm import exoplasim_run as er
from climate.gcm import turned
from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import design_settled_turn as DS
from research.studies.atlas_round import design_turned_surface as DT
from research.studies.atlas_round import settled_turn as ST
from research.studies.atlas_round import turned_surface as TS


def test_the_design_surfaces_experiment_takes_only_its_own_surface(tmp_path):
    lakes = er.EXPERIMENTS['A28D']['lakes']
    assert lakes != er.EXPERIMENTS['A28S']['lakes']
    path = tmp_path / 'design_turned_surface_x.json'
    for surface, match in ((f'moon_28pct_water_lakes_{lakes}_turned_x_gaussian_T21.nc', 'missing'),
                           (f"moon_28pct_water_lakes_{er.EXPERIMENTS['A28S']['lakes']}_turned_x_gaussian_T21.nc",
                            'another surface')):
        path.write_text(json.dumps(dict(schema=turned.SCHEMA, surface=surface, files={surface: '0' * 64},
                                        frame=dict(state='settled', tilt_to_ecliptic_deg=1.7))))
        with pytest.raises(ValueError, match=match):
            er.frame_setting(path, er.EXPERIMENTS['A28D'])


@pytest.fixture(scope='module')
def record():
    if not (ST.OUT.exists() and DS.OUT.exists()):
        pytest.skip('run python -m research.studies.atlas_round.design_settled_turn first')
    path = DT.record_path(hashlib.sha256(DS.OUT.read_bytes()).hexdigest())
    if not path.exists():
        pytest.skip('run python -m research.studies.atlas_round.design_turned_surface first')
    final = TS.record_path(hashlib.sha256(ST.OUT.read_bytes()).hexdigest())
    return path, json.loads(path.read_text()), json.loads(final.read_text())


def test_the_design_water_settles_to_its_own_turn():
    if not DS.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.design_settled_turn first')
    d, f = json.loads(DS.OUT.read_text()), json.loads(ST.OUT.read_text())
    assert d['final_water_mass_kg']['all'] < f['final_water_mass_kg']['all']   # the parted sea holds less water
    against = d['against_the_final_water']
    assert against['settled'] == f['cases']['final_water']['settled']
    assert 0.0 < against['frame_turn_deg'] < 3.0
    assert d['series_settled']['reproduced'] == f['series_settled']['reproduced']   # the method reproduces as before
    assert DR.rel(DR.OUT) in d['inputs'] and str(DR.PRODUCT) in d['ignored_inputs']


def test_the_record_turns_the_design_water_into_its_settled_frame(record):
    _, r, final = record
    st = json.loads(DS.OUT.read_text())
    assert r['frame']['rotation_today_to_turned'] == st['frame_today_to_turned']
    assert r['frame']['tilt_to_ecliptic_deg'] == st['cases']['final_water']['settled']['tilt_to_ecliptic_deg']
    assert r['against_the_final_surface']['frame'] == final['frame']
    assert r['against_the_final_surface']['frame_turn_deg'] == pytest.approx(
        st['against_the_final_water']['frame_turn_deg'], abs=1e-3)
    today, now = r['shares']['today'], r['shares']['turned']
    assert now['water'] == pytest.approx(today['water'], abs=1e-3)
    assert now['seas'] == pytest.approx(today['seas'], abs=1e-3)
    assert today['seas'] == pytest.approx(json.loads(DR.OUT.read_text())['summary']['sea_share'], abs=2e-3)
    assert today['seas'] < final['shares']['today']['seas']                    # the parted sea leaves its floor dry
    assert r['water']['lake_cells_on_seas'] == 0
    lake_cells = json.loads(DR.OUT.read_text())['summary']['lake_cells']        # lakes' rims at their own level, depth 0
    assert r['water']['wet_cells_at_or_below_their_ground'] < 0.01 * lake_cells
    assert r['t21']['turned']['water_share'] == pytest.approx(now['water'], abs=1e-3)
    for name, sha in r['files'].items():
        local = TS.PRODUCTS / name
        if local.exists():
            assert hashlib.sha256(local.read_bytes()).hexdigest() == sha


def test_the_runner_takes_the_record_for_the_design_surfaces_experiment(record):
    path, r, _ = record
    assert r['ignored_inputs'][str(DR.PRODUCT)][:16] == er.EXPERIMENTS['A28D']['lakes']
    if not (er.PRODUCTS / r['surface']).is_file():
        pytest.skip('the surface file is not in products')
    frame = er.frame_setting(path, er.EXPERIMENTS['A28D'])
    assert (frame['surface'], frame['tilt_deg'], frame['state']) == \
        (r['surface'], r['frame']['tilt_to_ecliptic_deg'], TS.STATE)
