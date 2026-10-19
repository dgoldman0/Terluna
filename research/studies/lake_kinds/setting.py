"""The resort at 301 K, built as the waters' oxygen study builds it at 300 K.

The author set the resort at 301 K on 18 October (the register on domain/climate, from climate_targets' best value
within 300-302 K). The waters' oxygen study (research/studies/lake_oxygen) builds its resort from the lock's weather,
the design case's 3-day means, moved by the share of the way between the design case and the hot world's test that
climate_targets gives the target on the ground's temperature, with the lapse rate that carries the air from the GCM
cell down to each lake taken from the target's own sky ladder over land. Here both come from the same product's two
runs: the share is (T - T_design) / (T_test - T_design), as climate_targets' interpolate() reckons it, and the ladder
is interpolated height by height with that share, the lapse rate read from it between 0.9 and 5 km as lake_oxygen
reads it. At 300 K this reproduces climate_targets' own share and ladder exactly, and so lake_oxygen's forcing.

Everything that runs the lake model at 301 K imports Setting301 from here.
"""
from __future__ import annotations

import json

from research.studies.lake_oxygen import run as R

TARGET_K = 301.0
CHECK_K = 300.0                          # the waters' oxygen study's resort, reproduced to prove the construction
LAPSE_LOW_KM, LAPSE_HIGH_KM = 0.9, 5     # lake_oxygen's Setting._lapse reads the ladder between these heights


def interpolation(targets, ground_k):
    """The share of the way to the hot world's test at a ground temperature, the sky ladder over land there and its
    lapse rate (K per m by which the air warms going down), from climate_targets' design and test runs."""
    a, b = targets['runs']['design'], targets['runs']['test']
    share = (ground_k - a['ground_k']) / (b['ground_k'] - a['ground_k'])
    ladder = [dict(height_km=x['height_km'], air_c=x['air_c'] + share * (y['air_c'] - x['air_c']),
                   pressure_atm=x['pressure_atm'] + share * (y['pressure_atm'] - x['pressure_atm']))
              for x, y in zip(a['ladder'], b['ladder'])]
    by = {r['height_km']: r['air_c'] for r in ladder}
    lapse = (by[LAPSE_LOW_KM] - by[LAPSE_HIGH_KM]) / ((LAPSE_HIGH_KM - LAPSE_LOW_KM) * 1000.0)
    climate = {k: a[k] + share * (b[k] - a[k]) for k in ('sea_median_c', 'land_air_2m_median_c', 'land_air_median_c',
                                                          'land_rain_mm_day')}
    return dict(ground_k=ground_k, share=share, lapse_k_per_m=lapse, ladder=ladder, climate=climate)


class Setting301(R.Setting):
    """lake_oxygen's setting with its resort moved to ground_k: the share of the way and the lapse rate that its
    resort_windows() and windows() read are replaced by the interpolation at that temperature; the lock, the hot world
    and every other input stay lake_oxygen's own."""

    def __init__(self, ground_k=TARGET_K, quick=False):
        super().__init__(quick=quick)
        self.ground_k = float(ground_k)
        self.interpolation = interpolation(self.targets, self.ground_k)
        self.share['resort'] = self.interpolation['share']
        self.lapse['resort'] = self.interpolation['lapse_k_per_m']


def check_against_climate_targets(st: R.Setting):
    """At 300 K the construction must give climate_targets' own share and ladder, and the lapse rate lake_oxygen
    reads from that target's ladder."""
    target = next(t for t in st.targets['targets'] if t['name'] == R.TARGET_NAMES['resort'])
    mine = interpolation(st.targets, target['ground_k'])
    theirs = R.Setting._lapse(st, 'resort')              # always the 300 K target's own ladder
    return dict(ground_k=target['ground_k'], share=mine['share'],
                climate_targets_share=target['interpolated_share_of_the_way'],
                share_equal=mine['share'] == target['interpolated_share_of_the_way'],
                ladder_equal=all(m['air_c'] == t['air_c'] for m, t in zip(mine['ladder'],
                                                                          target['climate']['ladder'])),
                lapse_k_per_m=mine['lapse_k_per_m'], lake_oxygen_lapse_k_per_m=theirs,
                lapse_equal=mine['lapse_k_per_m'] == theirs)


def describe(st: Setting301):
    i = st.interpolation
    return dict(ground_k=i['ground_k'], share_of_the_way=i['share'], lapse_k_per_km=i['lapse_k_per_m'] * 1000.0,
                seas_median_c=i['climate']['sea_median_c'], land_air_2m_median_c=i['climate']['land_air_2m_median_c'],
                land_air_0_9km_median_c=i['climate']['land_air_median_c'],
                land_rain_mm_day=i['climate']['land_rain_mm_day'],
                design_ground_k=st.targets['runs']['design']['ground_k'],
                test_ground_k=st.targets['runs']['test']['ground_k'])


if __name__ == '__main__':
    s = Setting301()
    print(json.dumps(describe(s), indent=1))
    print(json.dumps(check_against_climate_targets(Setting301(CHECK_K)), indent=1))
