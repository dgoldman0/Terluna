"""Directional spectrum transfer, coastal boundaries and physical energy units."""
from datetime import timedelta
import os
from pathlib import Path
import re

import numpy as np
import pytest

from climate.waves.coastal_checks import nesting_spectra
from climate.waves.model import ORIGIN
from climate.waves.pilot import execute
from climate.waves.shore import boundary_audit, coastal_files, stationary_boundary
from climate.waves.shore_checks import group_velocity, spectral_transport
from climate.waves.shore_region import remap_spectra
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


def spectral_fixture():
    return dict(locations=np.array([[0.,0.],[1.,0.]]),frequency=np.array([.01,.02,.04,.08]),
                direction=np.arange(4)*90.,spectra=np.arange(32,dtype=float).reshape(1,2,4,4)/100,
                available=np.array([[True,False]]))


def test_selected_directional_spectrum_preserves_variance_and_dry_locations(tmp_path):
    original = spectral_fixture()
    text = stationary_boundary(original,0)
    stationary = tmp_path/'stationary.spc'
    stationary.write_text(text.replace('RFREQ','AFREQ'))
    parsed=nesting_spectra(stationary)
    assert parsed['dates'] == [None] and parsed['frequency_reference'] == 'AFREQ'
    np.testing.assert_allclose(parsed['spectra'][0,0],original['spectra'][0,0],atol=1e-8,rtol=0)
    header, values = text.split('-99\n',1)
    header = header.replace('SWAN 1\n','SWAN 1\nTIME\n1\n')
    path = tmp_path/'boundary.nest'
    path.write_text(header+'-99\n'+ORIGIN.strftime('%Y%m%d.%H%M%S')+'\n'+values
                    +(ORIGIN+timedelta(hours=3)).strftime('%Y%m%d.%H%M%S')+'\n'+values)
    selected = nesting_spectra(path,{ORIGIN+timedelta(hours=3)})
    assert selected['dates'] == [ORIGIN+timedelta(hours=3)]
    np.testing.assert_allclose(selected['spectra'][0,0],original['spectra'][0,0],atol=1e-8,rtol=0)
    np.testing.assert_array_equal(selected['available'],original['available'])
    with pytest.raises(ValueError,match='absent'):
        nesting_spectra(path,{ORIGIN+timedelta(hours=6)})


def test_wet_boundary_requires_a_parent_spectrum():
    terrain = dict(longitude_deg=np.array([0.,1.]),latitude_deg=np.array([0.,1.]),
                   wet=np.ones((2,2),bool),metadata=dict(pixels_per_degree=1))
    boundary = spectral_fixture()
    with pytest.raises(ValueError,match='lacks parent spectra'):
        boundary_audit(terrain,boundary,0)
    terrain['wet'][0,1]=False
    assert boundary_audit(terrain,boundary,0)['missing_parent_spectra'] == 0


def test_group_speed_recovers_deep_and_shallow_water_limits():
    gravity=MOON_SURFACE_GRAVITY
    frequency=np.array([.1,.2,.4])
    np.testing.assert_allclose(group_velocity(frequency,10000),gravity/(4*np.pi*frequency),rtol=1e-12)
    np.testing.assert_allclose(group_velocity(frequency/10000,.5),np.sqrt(gravity*.5),rtol=1e-7)


def test_opposing_waves_retain_separate_incoming_and_outgoing_power():
    frequency=np.array([.1,.2,.4])
    direction=np.array([0.,90.,180.,270.])
    density=np.zeros((3,4))
    density[:,0]=2
    density[:,2]=1
    power=spectral_transport(frequency,direction,density,1000,1025,MOON_SURFACE_GRAVITY,[1.,0.])
    assert power['shoreward_w_m'] == pytest.approx(2*power['seaward_w_m'])
    assert power['net_shoreward_w_m'] == pytest.approx(power['seaward_w_m'])
    assert power['all_directions_w_m'] == pytest.approx(3*power['seaward_w_m'])
    np.testing.assert_allclose(power['vector_w_m'],[power['net_shoreward_w_m'],0.],atol=1e-10)


def test_coastal_spectrum_interpolation_preserves_shape_and_bounds_its_reach():
    source=spectral_fixture()
    source['locations']=np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.],[.5,.5]])
    shape=source['spectra'][0,0]
    source['spectra']=np.repeat(shape[None,None],5,axis=1)
    source['spectra'][0,-1]=1e6
    source['available']=np.array([[True,True,True,True,False]])
    xy=np.array([[.5,.5],[0.,0.],[1.,1.]])
    mapped,audit=remap_spectra(source,xy,np.array([True,True,False]),1.)
    np.testing.assert_allclose(mapped['spectra'][0,:2],np.repeat(shape[None],2,axis=0),atol=1e-14)
    assert not mapped['available'][0,2] and mapped['spectra'][0,2].sum()==0
    assert audit['maximum_nearest_distance_in_source_cells']<1
    with pytest.raises(ValueError,match='beyond'):
        remap_spectra(source,np.array([[4.,4.]]),np.array([True]),1.)


@pytest.fixture
def executable():
    value = os.environ.get('TERLUNA_SWAN_EXECUTABLE')
    if not value:
        pytest.skip('Set TERLUNA_SWAN_EXECUTABLE for the energy-transport control')
    return Path(value).resolve(strict=True)


def test_transport_matches_lunar_group_velocity_direction_and_true_energy(executable,tmp_path):
    # A small flat sea, one frequency bin and one direction; all source terms
    # are disabled except depth breaking, which stays zero in deep water.
    x = np.arange(5)/32
    terrain = dict(longitude_deg=x,latitude_deg=x,depth_m=np.full((5,5),1000.))
    wind = dict(lon=x,lat=x,u=np.zeros((5,5)),v=np.zeros((5,5)))
    locations=np.array([(v,x[0]) for v in x[1:]]+[(x[-1],v) for v in x[1:]]+
                       [(v,x[-1]) for v in x[-2::-1]]+[(x[0],v) for v in x[-2::-1]])
    ratio=MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    boundary=dict(locations=locations,frequency=np.geomspace(.03*ratio,3*ratio,49),
                  direction=np.arange(36)*10.+5,available=np.ones((1,len(locations)),bool))
    cases={}
    for direction, true_energy in ((5,False),(95,False),(5,True)):
        density=np.zeros((1,len(locations),49,36))
        density[:,:,30,(direction-5)//10]=.001
        boundary['spectra']=density
        settings=f'SET GRAV {MOON_SURFACE_GRAVITY:.12g} RHO 1000 INRHOG {int(true_energy)}'
        files=coastal_files(terrain,stationary_boundary(boundary,0),wind,settings,stride=1)
        files['INPUT']=re.sub(r'INPGRID WIND[^\n]*\nREADINP WIND[^\n]*\n','',files['INPUT'])
        files['INPUT']=files['INPUT'].replace('GEN3 KOMEN DRAG WU AGROW','GEN3 KOMEN DRAG WU\nOFF QUAD\nOFF WCAP')
        directory=tmp_path/f'd{direction}_energy{true_energy}'
        execute(executable,directory,files,['shore.tbl'])
        data=np.loadtxt(directory/'shore.tbl').reshape(5,5,11)[1:-1,1:-1].reshape(-1,11)
        assert np.all(data[:,7] == 0)
        variance=data[:,2]**2/16
        cg=MOON_SURFACE_GRAVITY*data[:,3]/(4*np.pi)
        scale=1000*MOON_SURFACE_GRAVITY if true_energy else 1.
        expected=variance*cg*scale
        heading=np.radians(direction)
        expected_vector=expected[:,None]*[np.cos(heading),np.sin(heading)]
        np.testing.assert_allclose(data[:,8:10],expected_vector,rtol=.002)
        cases[direction,true_energy]=data
    np.testing.assert_array_equal(cases[5,False][:,:8],cases[5,True][:,:8])
    np.testing.assert_allclose(cases[5,True][:,8:10],cases[5,False][:,8:10]*1000*MOON_SURFACE_GRAVITY,rtol=1e-6)
