"""Check physical forcing invariants and population accounting in attribution runs."""
import json
import re
from pathlib import Path
import numpy as np
import pytest
from atmosphere.snow import scale_size_diagnosis as d


def test_earth_scaling_preserves_pressure_path_and_column_water_per_seed(tmp_path):
    source=d.lunar.CACHE/'lunar/manifest.json'
    if not source.exists(): pytest.skip('Optional pinned SCALE build cache')
    (tmp_path/'lunar').mkdir()
    (tmp_path/'lunar/manifest.json').write_bytes(source.read_bytes())
    case=d.configure('test_scaled',gravity=d.STANDARD_GRAVITY,depth_factor=d.DEPTH_RATIO,cache=tmp_path)
    s=json.loads((case/'settings.json').read_text())
    _,_,original,states=d.lunar.specification()
    assert s['depth_m']*s['gravity_m_s2']==pytest.approx(original['depth_m']*d.MOON_SURFACE_GRAVITY)
    assert s['profile']==states
    forcing=(case/'lunar_forcing.txt').read_text().splitlines()
    depth,dx,flux,rate,duration,start=map(float,forcing[1].split())
    assert rate*depth==pytest.approx(3/3600)
    assert rate*depth/flux==pytest.approx(8.144044315643984e-8,rel=1e-7)
    assert duration==original['total_time_s']
    for fname in ('init.conf','run.conf'):
        text=(case/fname).read_text()
        g=float(re.search(r'CONST_GRAV\s*=\s*([^,]+)',text).group(1).replace('d','e'))
        dz=float(re.search(r'\bDZ\s*=\s*([^,]+)',text).group(1).replace('d','e'))
        assert g==pytest.approx(d.STANDARD_GRAVITY)
        assert dz*s['layers']==pytest.approx(depth)
    supplied=np.array([list(map(float,line.split())) for line in forcing[3:]])
    np.testing.assert_allclose(supplied,np.array([[v[k] for k in ('t','p','e')] for v in states]),rtol=0,atol=0)


def test_arrival_mass_and_constituent_number_close_independently(tmp_path):
    # 20 one-monomer particles plus 10 three-monomer aggregates: 50 original
    # crystals in total; the warm rain packet must not enter the ice account.
    a=np.array([[70,0,20,1e-6,.001,300,1,1,10],
                [80,1,10,4e-6,.002,200,.5,3,10],
                [80,1,100,8e-6,.003,1000,1,7,1]],dtype='>f8')
    a.tofile(tmp_path/'lunar_arrivals.bin')
    settings=dict(gather_start_s=60.,duration_s=120.,horizontal_cells=2,dx_m=1.,
        original_settings=dict(injected_number_flux_m2_s=50/120,snowfall_mm_h=3.))
    p=d.population(tmp_path,dict(settings=settings))
    assert p['number_flux_m2_s']==pytest.approx(.25)
    assert p['monomer_flux_to_injected']==pytest.approx(1)
    assert p['number_survival_fraction']==pytest.approx(.6)
    assert p['mean_monomers']==pytest.approx(5/3)
    assert p['mean_mass_kg']*p['number_flux_m2_s']==pytest.approx(5e-7)
    assert p['mass_per_arriving_monomer_kg']==pytest.approx(1.2e-6)
    assert p['near_sphere_mass_fraction']==pytest.approx(1/3)


@pytest.mark.parametrize('kw',[{'gravity':0},{'depth_factor':float('nan')},{'source_factor':-1},{'water_factor':0}])
def test_invalid_factors_leave_no_run_directory(tmp_path,kw):
    with pytest.raises(ValueError):d.configure('invalid_factor',cache=tmp_path,**kw)
    assert not list(tmp_path.iterdir())


def test_independent_source_and_water_scalings(tmp_path):
    source=d.lunar.CACHE/'lunar/manifest.json'
    if not source.exists(): pytest.skip('Optional pinned SCALE build cache')
    (tmp_path/'lunar').mkdir()
    (tmp_path/'lunar/manifest.json').write_bytes(source.read_bytes())
    case=d.configure('test_supply',source_factor=.1,water_factor=10.,depth_factor=d.DEPTH_RATIO,cache=tmp_path)
    s=json.loads((case/'settings.json').read_text())
    depth,_,flux,rate,_,_=map(float,(case/'lunar_forcing.txt').read_text().splitlines()[1].split())
    assert flux==pytest.approx(s['original_settings']['injected_number_flux_m2_s']/10)
    assert rate*depth==pytest.approx(30/3600)
    assert rate*depth/flux==pytest.approx(8.144044315643984e-6,rel=1e-7)


def test_exact_cdf_finds_difference_between_sparse_queries():
    # Queries at .001 and .002 would see identical 0/1 CDFs. The exact
    # empirical distributions are completely separated between the points.
    q=d.weighted_cdf_distance([.0011],[2.],[.0012],[100.])
    assert q['distance']==pytest.approx(1.)
    assert q['diameter_m']==pytest.approx(.0011)


def test_exact_cdf_combines_tied_dimensions_and_mass_weights():
    q=d.weighted_cdf_distance([1.,1.,3.],[1.,2.,1.],[1.,2.,3.],[1.,2.,1.])
    assert q['distance']==pytest.approx(.5)
    assert q['diameter_m']==pytest.approx(1.)


def test_exact_cdf_gives_each_column_equal_weight(tmp_path,monkeypatch):
    monkeypatch.setattr(d.lunar,'CACHE',tmp_path)
    for name,diameter in [('left',1.),('right',2.)]:
        path=tmp_path/'runs'/name;path.mkdir(parents=True)
        (path/'settings.json').write_text(json.dumps(dict(gather_start_s=0.,horizontal_cells=2)))
        rows=np.array([[1,0,1,100,diameter,300,1,1,10],
                       [1,1,1,1,3,300,1,1,10]],dtype='>f8')
        rows.tofile(path/'lunar_arrivals.bin')
    assert d.exact_column_cdf('left','right')['distance']==pytest.approx(.5)
