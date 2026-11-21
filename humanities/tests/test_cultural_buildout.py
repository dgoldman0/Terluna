"""Gross funding uses, no-receipt exposure and shared build-out controls."""
import json

import pytest

from humanities.economics import accounts, buildout


def test_buildout_regenerates_and_preserves_bound_annual_products():
    saved=json.loads((accounts.HERE/'results/buildout.json').read_text())
    assert buildout.results()==saved
    for path,identity in saved['inputs']['products'].items():
        assert accounts.digest(accounts.ROOT/path)==identity['sha256']


@pytest.mark.parametrize('case',['lower','middle','higher'])
def test_combined_buildout_contains_sports_once_and_reconciles_financing(case):
    out=buildout.results()
    sport=out['programmes']['sports']['cases'][case]
    full=out['programmes']['entertainment']['cases'][case]
    increment=out['entertainment_increment_after_sports'][case]
    assert increment['installed_assets_U']==30000000
    assert full['gross_uses_before_financing_U']==pytest.approx(sport['gross_uses_before_financing_U']+increment['total_U'])
    for row in (sport,full):
        debt=row['optional_inherited_debt']
        sources=debt['principal_U']+debt['other_funding_U']
        uses=sum(row['uses'].values())+debt['construction_interest_U']+debt['issuance_fees_U']+debt['debt_reserve_U']
        assert sources==pytest.approx(uses)


def test_runway_keeps_minimum_costs_without_success_pay_or_duplicate_marketing():
    out=buildout.results()['programmes']
    assert out['sports']['base_payroll_U_year']==28800000
    assert out['sports']['cases']['lower']['uses']['no_receipt_operating_runway_U']==pytest.approx((28800000+41041600)/2)
    assert out['entertainment']['cases']['higher']['uses']['no_receipt_operating_runway_U']==pytest.approx(34240000+53413600)
    assert out['entertainment']['cases']['higher']['uses']['renewal_reserve_U']==11500000


def test_promotion_changes_funding_without_changing_the_audience_target():
    out=buildout.results()['programmes']
    for name,target in [('sports',16000000),('entertainment',28000000)]:
        assert out[name]['paid_audience_target']==target
        lo=out[name]['cases']['lower']['uses']['initial_paid_promotion_U']
        hi=out[name]['cases']['higher']['uses']['initial_paid_promotion_U']
        assert lo==pytest.approx(target*0.8)
        assert hi==pytest.approx(target*8)
