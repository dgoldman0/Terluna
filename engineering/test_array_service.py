import pytest

from engineering.array_service import ServiceCase, coverage_multiplicity, missing_ray_fraction, simulate_service


def test_optical_overlap_integrates_exact_multiplicity():
    p = coverage_multiplicity(1.25, 1.5, 1, 1)
    assert sum(p.values()) == pytest.approx(1)
    assert sum(n*w for n, w in p.items()) == pytest.approx(1.25*1.5)
    assert missing_ray_fraction(0, p) == 0
    assert missing_ray_fraction(1, p) == pytest.approx(1)
    assert missing_ray_fraction(.1, coverage_multiplicity(1, 1, 1, 1)) == pytest.approx(.1)
    assert missing_ray_fraction(.1, coverage_multiplicity(2, 2, 1, 1)) == pytest.approx(.1**4)


def test_zero_failure_make_before_break_keeps_all_tiles_operational():
    p = simulate_service(ServiceCase(10000, failure_rate_per_year=0, duration_days=100))
    assert p['peak_missing_fraction'] == 0
    assert p['planned_completed'] > 0
    assert p['peak_in_transit'] > 0
    assert p['balance_residual'] == {'failed': 0, 'planned': 0}


def test_service_shortage_accumulates_renewal_backlog_without_inventing_failures():
    p = simulate_service(ServiceCase(10000, failure_rate_per_year=0, capacity_over_planned=.5,
                                    duration_days=200))
    assert p['final_renewal_waiting'] > 150
    assert p['failed_arrivals'] == 0
    assert p['balance_residual'] == {'failed': 0, 'planned': 0}


def test_failure_replacement_is_delayed_and_conserved():
    p = simulate_service(ServiceCase(1000, failure_rate_per_year=0, duration_days=10,
                                    capacity_over_planned=10000, correlated_failure_count=100,
                                    correlated_failure_day=1, failed_slot_response_days=2, step_days=.125))
    assert p['failed_arrivals'] == p['failed_completed'] == 100
    assert p['peak_missing_fraction'] == pytest.approx(.1)
    assert p['mean_missing_fraction'] == pytest.approx(.02, abs=.002)


def test_outage_removes_capacity_and_recovers_without_lost_jobs():
    base = dict(population=10000, failure_rate_per_year=0, duration_days=100)
    a = simulate_service(ServiceCase(**base))
    b = simulate_service(ServiceCase(**base, outage_start_days=20, outage_duration_days=10))
    assert b['peak_renewal_waiting'] > a['peak_renewal_waiting']
    assert b['balance_residual'] == a['balance_residual']


def test_time_refinement_of_finite_transit_shock():
    base = dict(population=1000, failure_rate_per_year=0, duration_days=20,
                capacity_over_planned=10000, correlated_failure_count=100,
                correlated_failure_day=1, failed_slot_response_days=2)
    coarse = simulate_service(ServiceCase(**base, step_days=.25))
    fine = simulate_service(ServiceCase(**base, step_days=.0625))
    assert coarse['mean_missing_fraction'] == pytest.approx(fine['mean_missing_fraction'], abs=.0013)


@pytest.mark.parametrize('q', [-.1, 1.1, float('nan')])
def test_bad_outage_fraction_rejected(q):
    with pytest.raises(ValueError):
        missing_ray_fraction(q, {1: 1})
