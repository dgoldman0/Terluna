"""Empirical units, habit axes and invalid representation checks."""
import pytest
from atmosphere.snow import takahashi_growth as lab
from atmosphere.snow import scale_crystal_growth as first


def test_extended_dendrite_fit_preserves_previous_benchmark():
    for minute in range(3,31):
        old=first.observed(minute);new=lab.observed(-14.4,minute)
        assert new['mass_kg']==old['mass_kg']
        assert new['equatorial_dimension_m']==old['diameter_m']
        assert new['polar_dimension_m']==old['thickness_m']


def test_column_span_uses_long_axis_not_plate_width():
    end=lab.observed(-5.3,30)
    assert 1.85e-3 < end['diameter_m'] < 1.87e-3
    assert .135e-3 < end['equatorial_dimension_m'] < .137e-3
    assert end['diameter_m']==end['polar_dimension_m']
    assert 7.55e-9 < end['mass_kg'] < 7.57e-9
    assert lab.fall_speed(-5.3,30)==pytest.approx(.286)


def test_visual_transcription_of_thick_plate_mass_and_dendrite_fall_speed():
    # Original row is 2.5e-8, not the easily misread 1.5e-8.
    p=lab.observed(-8.6,3)
    assert 1.32e-10 < p['mass_kg'] < 1.34e-10
    assert lab.fall_speed(-14.4,30)==pytest.approx(.297)


@pytest.mark.parametrize('temperature',[-8.6,-12.2,-16.5,-18.2,-22.])
def test_unrepresentable_initial_spheroids_are_rejected(temperature):
    # 916.8 is the explicitly pinned SCALE reference material density.
    with pytest.raises(ValueError,match='above bulk ice'):
        lab.initial(temperature,916.8)


@pytest.mark.parametrize('temperature',[-3.7,-5.3,-10.6,-14.4,-20.1])
def test_primary_initial_spheroids_preserve_measurements(temperature):
    import math
    p=lab.initial(temperature,916.8);obs=lab.observed(temperature,3)
    assert 4*math.pi/3*p['equatorial_radius_m']**2*p['polar_radius_m']*p['density_kg_m3']==pytest.approx(obs['mass_kg'],rel=2e-15)
    assert 2*p['equatorial_radius_m']==obs['equatorial_dimension_m']
    assert 2*p['polar_radius_m']==obs['polar_dimension_m']


@pytest.mark.parametrize('temperature,age',[(-4,3),(-5.3,2.9),(-10.6,12.1),(-14.4,float('nan'))])
def test_no_unmeasured_temperature_or_age_extrapolation(temperature,age):
    with pytest.raises(ValueError):lab.observed(temperature,age)


def test_measurement_reader_rejects_unidentified_bytes(tmp_path):
    from atmosphere.snow import takahashi_measurements as data
    p=tmp_path/'unidentified.xlsx';p.write_bytes(b'wrong source')
    with pytest.raises(ValueError,match='SHA256'):data.read_rows(p)


def test_recovered_measurements_retain_missing_axes_and_original_row_ids():
    from atmosphere.snow import takahashi_measurements as data
    from atmosphere.snow.scale_crystal_growth import CACHE
    p=CACHE/'growth_environment_sources/1991_takahashi.xlsx'
    if not p.exists():pytest.skip('Optional external workbook is not present')
    rows=data.read_rows(p)
    assert len(rows)==169
    assert len({r['excel_row'] for r in rows})==169
    assert sum(r['c_full_mm'] is None for r in rows)==79
    selected={r['excel_row']:r for r in rows}
    assert selected[114]['mass_kg']==pytest.approx(6.62e-8)
    assert selected[114]['a_full_mm']==pytest.approx(3.77)
    assert selected[114]['spheroid_density_kg_m3'] is None
    assert 887 < selected[51]['spheroid_density_kg_m3'] < 888
    assert selected[72]['c_full_mm']==pytest.approx(.8923)
    assert 'source_c_axis_decimal_anomaly_unresolved' in selected[72]['quality_flags']
    assert 'maximum_dimension_mm' in selected[72]['calibration_excluded_observables']


def test_growth_patch_rejects_unidentified_native_source():
    from atmosphere.snow import scale_growth_controls as controls
    with pytest.raises(ValueError,match='Unidentified native'):controls.patch('changed source')
