"""Physical controls preserve the checkpoint numerics and explicit forcing."""
import numpy as np
import pytest
from atmosphere.snow import flakes as sf, aggregation as ag, sectional
from atmosphere.snow.habits import GEOMETRIES
from atmosphere.snow.column_experiments import Kinetics, QuadraticContinuation, solve_column

G = GEOMETRIES['bf95_dmax_dendrites_area']


def test_extension_reproduces_checkpoint_column_with_same_physics():
    args = dict(geometry=G, layers=3, mass_bins=28, maximum_diameter_m=.03,
                fill=.04, gather=.06, max_step_s=13.)
    reference = sectional.column(sf.earth_column(), -2., 1e4, 1., **args)
    result = solve_column(sf.earth_column(), -2., 1e4, 1., **args)
    for key in ('arriving_numbers_m2', 'final_number_m3', 'time_mean_number_m3'):
        np.testing.assert_allclose(result[key], reference[key], rtol=1e-12, atol=1e-24)
    for key, value in reference['mass_ledger_kg_m2'].items():
        assert result['mass_ledger_kg_m2'][key] == pytest.approx(value, rel=1e-12, abs=1e-14)
    assert result['diagnostics']['steps'] == reference['diagnostics']['steps']


def test_outer_encounter_area_changes_kernel_without_changing_fall_speeds():
    d = np.geomspace(2e-5, .3, 25); m = G.mass(d)
    world = sf.earth_column(); state = world.state(world.height_of(sf.T0-10))
    v0, re0, k0 = Kinetics('legacy')(m, d, G, world.gravity, state)
    v1, re1, k1 = Kinetics('outer_area')(m, d, G, world.gravity, state)
    np.testing.assert_array_equal(v0, v1)
    np.testing.assert_array_equal(re0, re1)
    np.testing.assert_allclose(k1*G.area_ratio, k0, rtol=2e-15)
    np.testing.assert_array_equal(np.diag(k1), np.zeros(d.size))


def test_quadratic_continuation_preserves_seeds_and_matches_mass_at_join():
    geometry = QuadraticContinuation(G)
    d = np.array([2e-5, 1e-4, 4e-4, .01, .1, 1.])
    np.testing.assert_array_equal(geometry.mass(d[:4]), G.mass(d[:4]))
    np.testing.assert_allclose(geometry.size(geometry.mass(d)), d, rtol=3e-15)
    assert geometry.mass(.02)/geometry.mass(.01) == pytest.approx(4.)
    assert geometry.mass(.01*(1+1e-9))/G.mass(.01) == pytest.approx(1+2e-9, rel=1e-15)


class ScaledVelocity:
    """Controlled speed perturbation with collisions disabled."""
    def __init__(self, factor): self.factor = factor
    def description(self): return {'test_only_speed_factor': self.factor}
    def __call__(self, m, d, geometry, gravity, state):
        v, re = ag.speeds(d, m, gravity, state, geometry.area_ratio)
        return self.factor*v, re, np.zeros((len(m),len(m)))


def test_fixed_source_flux_and_duration_survive_changed_drag():
    args = dict(geometry=G, layers=2, mass_bins=24, maximum_diameter_m=.03,
                fill=.02, gather=.03, injection_velocity_m_s=.1)
    a = solve_column(sf.earth_column(), -2., 1e4, 1., kinetics=ScaledVelocity(1), **args)
    b = solve_column(sf.earth_column(), -2., 1e4, 1., kinetics=ScaledVelocity(2), **args)
    assert a['seconds'] == b['seconds']
    for key in ('injected_mass_flux_kg_m2_s', 'injected_number_flux_m2_s'):
        assert a['settings'][key] == b['settings'][key]
    assert b['settings']['implied_seed_number_per_m3'] == .5*a['settings']['implied_seed_number_per_m3']
    for result in (a,b):
        assert abs(result['mass_ledger_kg_m2']['relative_residual']) < 2e-13
        assert result['diagnostics']['minimum_number_m3'] >= 0
    assert a['mass_ledger_kg_m2']['injected'] == pytest.approx(b['mass_ledger_kg_m2']['injected'], rel=1e-14)


def test_limited_growth_column_preserves_mass_ledger_and_fixed_forcing():
    args=dict(geometry=G,layers=3,mass_bins=32,maximum_diameter_m=.03,
              fill=.04,gather=.06,max_step_s=13.)
    donor=solve_column(sf.earth_column(),-2.,1e4,1.,growth_scheme='donor',**args)
    limited=solve_column(sf.earth_column(),-2.,1e4,1.,growth_scheme='limited',**args)
    assert limited['settings']['growth_scheme']=='limited'
    assert limited['seconds']==donor['seconds']
    assert limited['settings']['injected_mass_flux_kg_m2_s']==donor['settings']['injected_mass_flux_kg_m2_s']
    assert abs(limited['mass_ledger_kg_m2']['relative_residual'])<1e-12
    assert limited['diagnostics']['minimum_number_m3']>=0
    assert np.dot(limited['arriving_numbers_m2'],limited['masses'])>0
    assert not np.array_equal(limited['final_number_m3'],donor['final_number_m3'])
    # The implemented cap applies to the actual discrete face-flux deposition.
    depth=limited['layer_edges_m'][-1]-limited['layer_edges_m'][0]
    supplied=1./3600.*limited['diagnostics']['final_time_s']
    assert limited['mass_ledger_kg_m2']['deposited']<=supplied*(1+1e-12)


def test_limited_growth_boundary_loss_is_recorded_without_clipping():
    result=solve_column(sf.earth_column(),-2.,1e4,3.,geometry=G,layers=2,mass_bins=16,
                        minimum_diameter_m=5e-5,maximum_diameter_m=4e-4,
                        fill=.1,gather=.1,growth_scheme='limited')
    assert result['mass_ledger_kg_m2']['upper_domain_exit']>0
    assert not result['domain_adequacy']['passed']
    assert abs(result['mass_ledger_kg_m2']['relative_residual'])<1e-12


def test_batch_growth_preserves_full_column_population_and_ledgers_exactly():
    args=dict(geometry=G,layers=3,mass_bins=32,maximum_diameter_m=.03,
              fill=.04,gather=.06,max_step_s=13.,growth_scheme='limited')
    loop=solve_column(sf.earth_column(),-2.,1e4,1.,batch_growth=False,**args)
    batch=solve_column(sf.earth_column(),-2.,1e4,1.,batch_growth=True,**args)
    for key in ('arriving_numbers_m2','final_number_m3','time_mean_number_m3'):
        np.testing.assert_array_equal(loop[key],batch[key])
    assert loop['mass_ledger_kg_m2']==batch['mass_ledger_kg_m2']
    assert loop['diagnostics']==batch['diagnostics']
