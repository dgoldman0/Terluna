"""Identity and extensive-moment checks on the separate repaired populations."""
import json
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import repair
from shared.provenance import constants_changed


def test_changed_column_input_is_rejected_before_running():
    relative='research/studies/climate_targets/results/climate_targets.json'
    with pytest.raises(RuntimeError,match='Input bytes'):
        repair.identity({'inputs':{relative:'0'*64}})


def test_published_repair_preserves_model_identity_and_arriving_moments():
    if not repair.OUT.exists():
        pytest.skip('Sectional comparison product has not been generated')
    product=json.loads(repair.OUT.read_text())
    assert product['schema']==repair.SCHEMA
    ids={row['id'] for row in product['runs']}
    assert ids.isdisjoint(product['unperformed_cases'])
    assert ids | set(product['unperformed_cases']) == {c['id'] for c in repair.cases()}
    for row in product['runs']:
        for key in ('files','inputs'):
            for path,expected in row['producer'][key].items():
                assert repair.sha(repair.ROOT/path)==expected,path
        assert not constants_changed(row['producer']['constants'])
        mass=np.asarray(row['arriving_mass_pivots_kg'])
        number=np.asarray(row['arriving_numbers_m2'])
        seconds=row['gathering_seconds']
        assert (number>=0).all() and np.isfinite(number).all()
        ledger=row['mass_ledger_kg_m2']
        supplied=ledger['injected']+ledger['deposited']
        escaped=ledger['lower_domain_exit']+ledger['upper_domain_exit']
        residual=supplied-ledger['freezing_exit']-escaped-ledger['remaining']
        assert residual==pytest.approx(ledger['residual'],abs=1e-12*supplied)
        assert ledger['relative_residual']==pytest.approx(residual/supplied,abs=1e-12)
        assert abs(residual/supplied) < 1e-10
        domain=row['domain_adequacy']
        assert domain['mass_fraction_exited']==pytest.approx(escaped/supplied,rel=1e-13)
        assert domain['passed'] is (escaped/supplied <= domain['tolerance'])
        assert row['numerical_status']==('domain_passed_refinement_pending' if domain['passed'] else 'finite_domain_inadequate')
        assert row['summary']['snowfall_mm_h']==pytest.approx(np.dot(mass,number)*3600/seconds,rel=1e-13)
        assert row['freezing_level']['totals']['mass_flux_kg_m2_s']==pytest.approx(np.dot(mass,number)/seconds,rel=1e-13)
        for layer in [row['freezing_level']]+row['layers']:
            n,a,m,nf,mf=[np.asarray(layer[k]) for k in ('number_m3','projected_area_m2','mass_kg_m3','number_flux_m2_s','mass_flux_kg_m2_s')]
            assert all(np.isfinite(x).all() and (x>=0).all() for x in (n,a,m,nf,mf))
            for key,value in [('number_m3',n.sum()),('mass_kg_m3',m.sum()),('projected_area_m_inv',np.dot(n,a)),
                              ('number_flux_m2_s',nf.sum()),('mass_flux_kg_m2_s',mf.sum())]:
                assert layer['totals'][key]==pytest.approx(value,rel=2e-13)
        assert row['domain_adequacy']['physical_size_cutoff'] is False
