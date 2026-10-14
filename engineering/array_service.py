"""Fleet renewal and failure queues; population accounting, not orbital traffic.

Every planned exchange leaves its old tile operating until its replacement arrives.
An unexpected failure opens a slot immediately. Replacement jobs have finite transit
time and share a finite dispatch capacity, with failed slots served first. The
factory/spare supply is assumed sufficient for that capacity; no material appears
as an orbital insertion or instant repair in this model.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import heapq
import math

import numpy as np

from shared import constants as K


@dataclass(frozen=True)
class ServiceCase:
    population: int
    life_years: float = 15.0
    failure_rate_per_year: float = 0.01
    capacity_over_planned: float = 1.2
    failed_slot_response_days: float = 1.0
    planned_transit_days: float = 27.0
    duration_days: float = K.JULIAN_YEAR_DAYS
    step_days: float = 0.125
    seed: int = 4917
    outage_start_days: float | None = None
    outage_duration_days: float = 0.0
    correlated_failure_count: int = 0
    correlated_failure_day: float = 90.0


def coverage_multiplicity(width_m: float, height_m: float, pitch_x_m: float,
                          pitch_y_m: float) -> dict[int, float]:
    """Exact area shares covered by 1, 2, or 4 tiles in a periodic rectangle grid.

    Each width/pitch is in [1,2]. This ideal flat service lattice is only an
    optical comparison; it is not the real curved fleet or its inclined rays.
    """
    values = (width_m, height_m, pitch_x_m, pitch_y_m)
    if not all(math.isfinite(x) and x > 0 for x in values):
        raise ValueError("Optical dimensions must be finite and positive")
    x, y = width_m / pitch_x_m, height_m / pitch_y_m
    if not (1 <= x <= 2 and 1 <= y <= 2):
        raise ValueError("This lattice requires one or two covers in each direction")
    return {1: (2-x)*(2-y), 2: (x-1)*(2-y)+(2-x)*(y-1), 4: (x-1)*(y-1)}


def missing_ray_fraction(off_fraction: float, multiplicity: dict[int, float]) -> float:
    """Expected uncovered fraction for *independent*, uniform whole-tile outages.

    Correlated failures invalidate independence. This result excludes baseline
    seams, pinholes, finite-Sun effects and spectral transmission through film.
    """
    if not math.isfinite(off_fraction) or not 0 <= off_fraction <= 1:
        raise ValueError("Outage fraction must lie in [0,1]")
    if (not multiplicity or any(n < 1 or int(n) != n or w < 0 for n, w in multiplicity.items())
            or not math.isclose(sum(multiplicity.values()), 1, abs_tol=1e-12)):
        raise ValueError("Multiplicity must be a probability distribution over positive integers")
    return float(sum(w * off_fraction**n for n, w in multiplicity.items()))


def simulate_service(case: ServiceCase) -> dict:
    """Seeded, aggregate count simulation of the entire named population.

    A uniformly staggered initial age distribution generates planned demand.
    Constant failure hazard is a scenario, not a measured component lifetime.
    Random failures do not cancel already-dispatched planned jobs: the counts
    deliberately charge both, a conservative capacity approximation. Planned
    backlog is a count of overdue renewal requests, not an invented failure rate.
    """
    for key in ("life_years", "capacity_over_planned", "failed_slot_response_days",
                "planned_transit_days", "duration_days", "step_days"):
        if not math.isfinite(getattr(case, key)) or getattr(case, key) <= 0:
            raise ValueError(f"{key} must be finite and positive")
    if case.population < 1 or int(case.population) != case.population:
        raise ValueError("population must be a positive integer")
    if not math.isfinite(case.failure_rate_per_year) or case.failure_rate_per_year < 0:
        raise ValueError("Failure hazard must be nonnegative")
    if (not math.isfinite(case.correlated_failure_count) or int(case.correlated_failure_count) != case.correlated_failure_count
            or case.correlated_failure_count < 0 or case.correlated_failure_count > case.population):
        raise ValueError("Correlated batch must fit population")
    if not math.isfinite(case.outage_duration_days) or case.outage_duration_days < 0:
        raise ValueError("Outage duration must be nonnegative")
    if not math.isfinite(case.correlated_failure_day) or case.correlated_failure_day < 0:
        raise ValueError("Correlated failure day must be finite and nonnegative")
    if case.outage_start_days is not None and (not math.isfinite(case.outage_start_days) or case.outage_start_days < 0):
        raise ValueError("Outage start must be finite and nonnegative")
    rng = np.random.default_rng(case.seed)
    planned_rate = case.population / (case.life_years*K.JULIAN_YEAR_DAYS)
    dispatch_rate = planned_rate * case.capacity_over_planned
    planned_pending = failed_pending = failed_transit = planned_transit = 0
    failed_arrivals = failed_completed = planned_arrivals = planned_completed = 0
    planned_dispatched = failed_dispatched = 0
    planned_fraction = capacity_fraction = 0.0
    events = []
    trace = []
    peak_missing = peak_planned = peak_in_transit = 0
    missing_integral = 0.0
    shock_done = False
    t = 0.0
    while t < case.duration_days-1e-12:
        dt = min(case.step_days, case.duration_days-t)
        end = t+dt
        # Complete transit jobs on this time grid; error bounded by one timestep.
        while events and events[0][0] <= end+1e-12:
            _, kind, count = heapq.heappop(events)
            if kind == "failed":
                failed_transit -= count
                failed_completed += count
            else:
                planned_transit -= count
                planned_completed += count
        operational = case.population-failed_pending-failed_transit
        chance = -math.expm1(-case.failure_rate_per_year*dt/K.JULIAN_YEAR_DAYS)
        failures = int(rng.binomial(operational, chance))
        if not shock_done and end >= case.correlated_failure_day:
            failures += min(case.correlated_failure_count, operational-failures)
            shock_done = True
        failed_pending += failures
        failed_arrivals += failures
        planned_fraction += planned_rate*dt
        arrivals = int(planned_fraction)
        planned_fraction -= arrivals
        planned_pending += arrivals
        planned_arrivals += arrivals
        # Only dispatch time inside an outage is removed; timestep cannot skip it.
        down = 0.0
        if case.outage_start_days is not None:
            down = max(0.0, min(end, case.outage_start_days+case.outage_duration_days)
                       -max(t, case.outage_start_days))
        capacity_fraction += dispatch_rate*(dt-down)
        slots = int(capacity_fraction)
        capacity_fraction -= slots
        repair = min(slots, failed_pending)
        renew = min(slots-repair, planned_pending)
        if repair:
            failed_pending -= repair
            failed_transit += repair
            failed_dispatched += repair
            heapq.heappush(events, (end+case.failed_slot_response_days, "failed", repair))
        if renew:
            planned_pending -= renew
            planned_transit += renew
            planned_dispatched += renew
            heapq.heappush(events, (end+case.planned_transit_days, "planned", renew))
        missing = failed_pending+failed_transit
        peak_missing = max(peak_missing, missing)
        peak_planned = max(peak_planned, planned_pending)
        peak_in_transit = max(peak_in_transit, failed_transit+planned_transit)
        missing_integral += missing*dt
        if not trace or math.floor(end) > math.floor(t) or end == case.duration_days:
            trace.append(dict(day=end, missing=missing, failed_waiting=failed_pending,
                              renewal_waiting=planned_pending, in_transit=failed_transit+planned_transit))
        t = end
    balance = dict(failed=failed_arrivals-failed_completed-failed_pending-failed_transit,
                   planned=planned_arrivals-planned_completed-planned_pending-planned_transit)
    if any(balance.values()):
        raise AssertionError("Population queue balance failed")
    return dict(configuration=asdict(case), planned_per_day=planned_rate,
                dispatch_capacity_per_day=dispatch_rate,
                expected_failure_per_day_at_full_population=case.population*case.failure_rate_per_year/K.JULIAN_YEAR_DAYS,
                capacity_over_all_expected_jobs=dispatch_rate/(planned_rate+case.population*case.failure_rate_per_year/K.JULIAN_YEAR_DAYS),
                failed_arrivals=failed_arrivals, failed_completed=failed_completed,
                planned_arrivals=planned_arrivals, planned_completed=planned_completed,
                failed_dispatched=failed_dispatched, planned_dispatched=planned_dispatched,
                final_failed_waiting=failed_pending, final_failed_transit=failed_transit,
                final_renewal_waiting=planned_pending, final_renewal_transit=planned_transit,
                mean_missing_fraction=missing_integral/(case.duration_days*case.population),
                peak_missing_fraction=peak_missing/case.population,
                peak_in_transit=peak_in_transit, peak_renewal_waiting=peak_planned,
                balance_residual=balance, daily_trace=trace,
                interpretation="Count simulation with finite transit and dispatch; no orbital insertion, physical debris, "
                "fault-containment success, spare manufacturing, or spatial coverage is demonstrated.")
