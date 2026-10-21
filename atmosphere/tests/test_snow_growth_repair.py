"""Admission checks for saved full-column growth-repair controls."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]


def test_saved_growth_repair_controls_keep_identity_and_water_budgets():
    product=json.loads((ROOT/'atmosphere/snow/results/growth_repair.json').read_text())
    assert product['schema']=='terluna.atmosphere.snow-growth-repair/1'
    assert {r['id'] for r in product['cases']} >= {'earth_limited_m144','moon_limited_m96'}
    for row in product['cases']+product['references']:
        for group in ('files','inputs'):
            for path,expected in row['producer'][group].items():
                assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected
        ledger=row['mass_ledger_kg_m2']
        source=ledger['injected']+ledger['deposited']
        exits=ledger['lower_domain_exit']+ledger['upper_domain_exit']
        residual=source-ledger['freezing_exit']-exits-ledger['remaining']
        assert abs(residual/source)<1e-10
        assert row['domain_adequacy']['mass_fraction_exited']==pytest.approx(exits/source,rel=1e-12)
        assert row['domain_adequacy']['passed']==(exits/source<=row['domain_adequacy']['tolerance'])
        assert not row['domain_adequacy']['physical_size_cutoff']
    for row in product['cases']:
        assert row['settings']['growth_scheme']=='limited'
        masses=np.array(row['masses_kg']); numbers=np.array(row['arriving_numbers_m2'])
        assert np.isfinite(numbers).all() and (numbers>=0).all()
        arriving_mass=np.dot(masses,numbers)
        assert row['summary']['snowfall_mm_h']==pytest.approx(arriving_mass*3600/row['gathering_seconds'],rel=1e-13)
        assert arriving_mass<=row['mass_ledger_kg_m2']['freezing_exit']*(1+1e-12)
        supplied=row['specification']['snowfall_mm_h']/3600*row['diagnostics']['final_time_s']
        assert row['mass_ledger_kg_m2']['deposited']<=supplied*(1+1e-10)
