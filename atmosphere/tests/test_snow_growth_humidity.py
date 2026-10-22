"""Input identity and physical mapping checks for the independent humidity test."""
import hashlib
import json
import math
from pathlib import Path
import pytest
from atmosphere.snow import growth_humidity_validation as humidity
from atmosphere.snow.flakes import R_VAPOUR, ice_saturation_pressure


def row(id_='G'):
    return next(r for r in json.loads(humidity.DATA.read_text())['cases'] if r['id']==id_)


def test_density_excess_is_not_mistaken_for_relative_humidity():
    r=row();c=humidity.configuration(r,'native',10.,.05,919.)
    t=c['temperature_K']
    recovered=(c['vapor_pressure_Pa']-ice_saturation_pressure(t))/(R_VAPOUR*t)
    assert recovered==pytest.approx(.17e-3,rel=2e-14)
    assert c['pressure_Pa']==101000.
    assert c['duration_s']==600.
    assert c['growth_mode']==2


@pytest.mark.parametrize('span_um',[5.,10.,20.])
def test_reported_full_seed_span_and_dense_mass_mapping(span_um):
    c=humidity.configuration(row(),'both',span_um,.0125,919.)
    assert c['initial_equatorial_m']==pytest.approx(span_um*1e-6)
    assert c['initial_polar_m']==c['initial_equatorial_m']
    volume=math.pi/6*c['initial_equatorial_m']**2*c['initial_polar_m']
    assert c['initial_mass_kg']/volume==pytest.approx(919.)
    assert c['mean_mode']==1 and c['directional_mode']==1


def test_temperature_probe_preserves_observed_density_excess():
    base=humidity.configuration(row(),'native',10.,.05,919.)
    cold=humidity.configuration(row(),'native',10.,.05,919.,-.3)
    assert cold['temperature_K']==pytest.approx(base['temperature_K']-.3)
    assert cold['vapor_pressure_Pa']!=base['vapor_pressure_Pa']
    assert (cold['vapor_pressure_Pa']-ice_saturation_pressure(cold['temperature_K']))/(R_VAPOUR*cold['temperature_K'])==pytest.approx(.17e-3)


@pytest.mark.parametrize('kwargs',[{'dt':0.},{'dt':float('nan')},{'dt':.07},{'seed_span_um':4.},{'ice_density':0.},{'temperature_offset_C':float('nan')},{'mode':'invented'}])
def test_invalid_experimental_inputs_rejected(kwargs):
    values=dict(row=row(),mode='native',seed_span_um=10.,dt=.05,ice_density=919.,temperature_offset_C=0.)
    values.update(kwargs)
    with pytest.raises(ValueError):humidity.configuration(**values)


def test_observation_transcription_preserves_printed_ids_and_uncertainty():
    d=json.loads(humidity.DATA.read_text())
    assert [r['id'] for r in d['cases']]==list('ABCDEFGHIJKLM')
    assert row('G')['observed_diameter_mm']==1.46
    assert row('H')['vapor_density_excess_g_m3']==.02
    assert d['initial_particle_mass_kg'] is None
    assert d['initial_aspect_ratio'] is None
    assert 'nominal group' in d['temperature_status']


def test_observation_source_tampering_is_rejected(tmp_path,monkeypatch):
    pdf=tmp_path/'source.pdf';pdf.write_bytes(b'controlled source bytes')
    digest=hashlib.sha256(pdf.read_bytes()).hexdigest()
    d=json.loads(humidity.DATA.read_text());d['source']['sha256']=digest
    metadata=tmp_path/'source.json';metadata.write_text(json.dumps(d))
    monkeypatch.setattr(humidity,'SOURCE',pdf);monkeypatch.setattr(humidity,'SOURCE_SHA',digest);monkeypatch.setattr(humidity,'DATA',metadata)
    assert len(humidity.observations()['cases'])==13
    pdf.write_bytes(b'different source')
    with pytest.raises(ValueError,match='source identity'):humidity.observations()


def test_refinement_records_missing_endpoint_as_incomplete():
    assert humidity.compare_refinement({'complete':False},{'complete':True})=={'complete_pair':False}
