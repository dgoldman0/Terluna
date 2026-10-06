import json

import numpy as np

from shared.provenance import constants_changed
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from . import exhaust_isolation as ex

PRODUCT = json.loads(ex.OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.exhaust-isolation/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(ex.FILES)
    for name, value in producer['files'].items():
        assert ex.digest(ex.ROOT/name) == value, name
    assert producer['inputs']['results/holding.json'] == ex.digest(ex.HOLDING)
    assert producer['inputs']['results/zoned_aperture.json'] == ex.digest(ex.ZONED)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == ex.digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_the_allowance_follows_the_september_diagnostic():
    allowance = PRODUCT['allowance']
    assert abs(allowance['binding_energy_MJ_kg']-.94) < .01
    for row in allowance['per_budget']:
        assert np.isclose(row['deposited_power_MW'], row['budget_kg_s']*allowance['binding_energy_MJ_kg']*10)
        assert np.isclose(row['share_of_jet_power'], row['deposited_power_MW']*1e6/(allowance['jet_power_TW']*1e12))


def test_a_narrow_cone_on_the_axis_takes_the_axial_intensity():
    axis = np.array([[0., 0., 1.]])
    half = np.radians(1.)
    for profile in ex.PLUMES.values():
        share = ex.cap_share(profile, axis, axis, np.array([half]))[0]
        expect = ex.intensity(profile, np.array(0.))*2*np.pi*(1-np.cos(half))/ex.total_intensity(profile)
        assert np.isclose(share, expect, rtol=.05)


def test_pair_axes_are_unit_vectors_canted_from_the_resultant_and_clear_of_the_moon():
    force = np.array([[1., .1, -.05]])
    moon = np.array([[-8e7, 3e6, 1e6]])
    for plane in ('clear', 'through_moon'):
        a, b = ex.plume_axes(force, moon, np.radians(45.), plane)
        for axis in (a, b):
            assert np.isclose(ex.length(axis), 1.).all()
            assert np.isclose(np.degrees(np.arccos(-np.sum(axis*ex.unit(force), axis=-1))), 45.).all()
    a, b = ex.plume_axes(force, moon, np.radians(45.), 'clear')
    angles = [np.degrees(np.arccos(np.sum(x*ex.unit(moon), axis=-1)))[0] for x in (a, b)]
    assert min(angles) >= 45.-1e-9


def test_canting_further_sends_less_of_every_plume_into_the_protected_sphere():
    rows = PRODUCT['direct_path']['clear']
    assert [r['cant_deg'] for r in rows] == sorted(r['cant_deg'] for r in rows)
    for name in ex.PLUMES:
        shares = [r[f'{name}_share'] for r in rows]
        assert all(later < earlier for earlier, later in zip(shares, shares[1:]))
    assert PRODUCT['direct_path']['through_moon'][0]['clearance_min_deg'] < rows[0]['clearance_min_deg']


def test_hall_plumes_send_far_more_into_the_protected_sphere_than_gridded_ion_plumes():
    row = PRODUCT['direct_path']['clear'][0]
    assert row['hall_BPT4000_share'] > 10*row['gridded_ion_NEXT_share']


def test_every_slow_gas_particle_has_one_fate():
    gas = PRODUCT['slow_gas']
    for name in ex.DESIGN['molar_mass_kg_mol']:
        g = gas[name]
        assert np.isclose(g['protected_sphere']+g['near_earth']+g['escaped']+g['still_flying'], 1.)
        cumulative = g['cumulative_by_day']
        assert all(b >= a for a, b in zip(cumulative, cumulative[1:]))
        assert np.isclose(cumulative[-1], g['protected_sphere'])
        delivered = g['delivered_by_lifetime']
        assert delivered['10_days'] <= delivered['30_days'] <= delivered['100_days'] <= g['protected_sphere']


def test_halving_the_integrator_step_keeps_the_slow_gas_fates():
    check = PRODUCT['slow_gas']['integrator_check']
    assert check['check_step_s'] == check['step_s']/2
    coarse, fine = check['protected_sphere_count']
    assert abs(coarse-fine) <= max(2, .05*fine)
    assert check['fate_changed'] <= .01*check['particles']


def test_the_slow_gas_alone_brings_more_than_the_allowance_for_one_kg_per_second():
    allowance = {row['budget_kg_s']: row['deposited_power_MW'] for row in PRODUCT['allowance']['per_budget']}
    for species in PRODUCT['slow_gas_implied']['by_species'].values():
        assert min(species['30_days']['deposited_MW']) > allowance[1.]
