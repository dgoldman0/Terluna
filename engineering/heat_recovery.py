"""Work recovered from a steady heat source, and the radiator area that rejects what remains, SI throughout.

The source's heat flow and temperature are imposed and stay as they are when recovery is added, and all recovered
electricity is used and ends as heat inside the boundary. A stream that cools through a range of temperatures needs
its own available-work account; converter mass, the source's response and temperature drops are left out. Radiator
area is total emitting surface, so a panel radiating from both faces has half of it as footprint.
"""
from math import isfinite

from shared.constants import STEFAN_BOLTZMANN


def _nonnegative(**values):
    for name, value in values.items():
        if not isfinite(value) or value < 0:
            raise ValueError(f'{name} must be finite and nonnegative')


def recover_heat(heat_W, hot_K, cold_K, carnot_fraction,
                 recovered_overhead_fraction=0., external_support_W=0.):
    """Split existing heat into engine rejection and recovered work.

Overhead drawn from the recovered work comes out of the processor's share and adds no heat; only support power
brought from outside adds to the heat rejected. carnot_fraction scales the Carnot efficiency, the ideal bound.
"""
    _nonnegative(heat_W=heat_W, hot_K=hot_K, cold_K=cold_K,
                 carnot_fraction=carnot_fraction,
                 recovered_overhead_fraction=recovered_overhead_fraction,
                 external_support_W=external_support_W)
    if not 0 < cold_K <= hot_K:
        raise ValueError('Temperatures must satisfy 0 < cold_K <= hot_K')
    if carnot_fraction > 1 or recovered_overhead_fraction > 1:
        raise ValueError('Fractions must lie in [0, 1]')
    ceiling = 1-cold_K/hot_K
    work = heat_W*ceiling*carnot_fraction
    overhead = work*recovered_overhead_fraction
    processor = work-overhead
    engine_heat = heat_W-work
    rejection = engine_heat+processor+overhead+external_support_W
    return dict(carnot_efficiency=ceiling,
                assumed_conversion_efficiency=ceiling*carnot_fraction,
                recovered_electric_W=work,
                recovered_overhead_W=overhead,
                processor_electric_W=processor,
                engine_rejection_W=engine_heat,
                external_support_W=external_support_W,
                internal_heat_to_radiators_W=rejection,
                added_internal_heat_W=rejection-heat_W,
                energy_balance_residual_W=rejection-heat_W-external_support_W)


def radiator_area(heat_W, temperature_K, emissivity,
                  background_K=0., absorbed_external_W_m2=0.):
    """Emitting area that rejects internal heat to a uniform background at background_K, viewed fully.

absorbed_external_W_m2 is the external power each square metre already absorbs besides that background (sunlight,
planetary infrared). It grows with the area, so the gross radiated power exceeds the internal heat by what the area
absorbs.
"""
    _nonnegative(heat_W=heat_W, temperature_K=temperature_K,
                 emissivity=emissivity, background_K=background_K,
                 absorbed_external_W_m2=absorbed_external_W_m2)
    if temperature_K == 0 or not 0 < emissivity <= 1:
        raise ValueError('Positive temperature and emissivity in (0, 1] required')
    net = emissivity*STEFAN_BOLTZMANN*(temperature_K**4-background_K**4)-absorbed_external_W_m2
    if net <= 0:
        raise ValueError('No positive net heat rejection at the specified temperature')
    area = heat_W/net
    absorbed = area*(emissivity*STEFAN_BOLTZMANN*background_K**4+absorbed_external_W_m2)
    gross = area*emissivity*STEFAN_BOLTZMANN*temperature_K**4
    return dict(net_flux_W_m2=net, emitting_area_m2=area,
                absorbed_environment_W=absorbed, gross_radiated_W=gross,
                net_heat_rejected_W=gross-absorbed)
