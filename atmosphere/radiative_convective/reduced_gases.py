"""How much methane and N2O warm the Moon's deep column, line by line.

    python -m atmosphere.radiative_convective.reduced_gases      # writes results/reduced_gases.json

The ecology domain asks how far the reduced gases warm the Moon (biosphere/ecology/questions.md, item 3), and the
author made the answer the next work on 12 October 2026: the levels at which the gases are held wait for it. In the
bare column under the shield no OH is made, so methane leaves only into soils and settles near 22-29 ppm at central
values, and N2O has no sink on any human timescale and builds 0.36-0.55 ppb a year (studies/domain_synthesis). This
calculation adds HITRAN's CH4 and N2O lines (inputs_added.json) to the clear-sky columns of the inverse
climate calculation (climate.py) and records, at fixed temperature profiles, the drop in outgoing longwave radiation
(the instantaneous clear-sky forcing at the top, no stratospheric adjustment) and the change of net flux at the
tropopause. The Moon's 1.2 atm column holds about seven times Earth's air above each square metre, so each ppm is a
seven times thicker absorber and the bands saturate sooner. A gas's optical depth is linear in its amount at these
fractions, so each column's optical depth per ppm is computed once and scaled.

The Earth control at 1 atm checks the method against AR6's effective radiative forcing for the gases' rise from 1750
to 2019 (Forster et al. 2021; the concentrations from Gulev et al. 2021, IPCC AR6 WGI). The ratio of AR6's value
to this clear-sky instantaneous one converts the Moon's values to all-sky adjusted estimates, as trace_gases.py does
for the fluorinated gases; that ratio carries Earth's clouds and adjustments to the Moon.

Warming follows through two responses: the column's own clear-sky one, the change of outgoing longwave less absorbed
sunlight with surface temperature in the inverse-climate rows near the design case (fixed relative humidity, so the
water-vapour feedback is in, clouds are not), and the 1.0-1.6 W/m2/K of provisioning's heat block. Methane's
near-infrared absorption of sunlight is not included.

The two gases join the column's spectroscopy only while this calculation runs (register): their line files and
isotopologues are added to spectroscopy.py's tables, and their inputs are read from inputs_added.json. Other
products bind spectroscopy.py, optics.py and inputs.json by hash, so those stay as they are, and nothing changes for
the gases they already hold.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

from atmosphere.radiative_convective import climate as cl, thermodynamics as th, optics, longwave as lw
from atmosphere.radiative_convective import spectroscopy as sp, fetch_inputs, fetch_added_inputs

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'reduced_gases.json'
SCHEMA = 'terluna.atmosphere.reduced-gas-forcing/1'
INVERSE = HERE / 'results' / 'inverse_climate.csv'
POPULATION = ROOT / 'provisioning' / 'results' / 'population.json'
CODE = ('reduced_gases.py', 'optics.py', 'spectroscopy.py', 'longwave.py', 'climate.py', 'thermodynamics.py',
        'fetch_added_inputs.py')
# The added gases' line files and isotopologues, from https://hitran.org/docs/iso-meta/ (retrieved 2026-10-12).
GAS_LINES = {'CH4': ('hitran_ch4_0-3500.par',), 'N2O': ('hitran_n2o_0-3500.par',)}
GAS_ISOTOPOLOGUES = {
    (4, '1'): (21, 44.001062), (4, '2'): (22, 44.998096), (4, '3'): (23, 44.998096),
    (4, '4'): (24, 46.005308), (4, '5'): (25, 45.005278),
    (6, '1'): (32, 16.031300), (6, '2'): (33, 17.034655), (6, '3'): (34, 17.037475), (6, '4'): (35, 18.040830),
}
LINES = ('hitran_ch4_0-3500.par', 'hitran_n2o_0-3500.par', 'hitran_h2o_0-3500.par', 'hitran_co2_0-3500.par')
S_MIN = 1e-28                    # the 296-K intensity floor the column uses for CO2

SCENARIOS = dict(earth_1atm=cl.Scenario('earth_1.0atm', th.EARTH, 101325.0, 400.0),
                 moon_12atm=cl.Scenario('moon_1.2atm', th.MOON, 121590.0, 400.0))
TS = (290.0, 295.0, 300.0)
DESIGN_TS = 295.0                # the design case near 295 K (294.9 K, climate/gcm)
# The synthesis's levels: Earth's preindustrial and present air, the screens' steady states and build-ups.
CH4_PPM = (0.7292, 1.8663, 3.4, 10.0, 21.6, 29.0, 60.0, 100.0, 230.0)
N2O_PPM = (0.2701, 0.3321, 0.55, 1.0, 2.0, 3.1, 4.2)
BOTH_PPM = ((1.8663, 0.3321), (21.6, 0.55), (29.0, 0.55), (230.0, 4.2))
# AR6 effective radiative forcing for 1750-2019 (W/m2) and the concentrations it spans (ppm).
AR6 = dict(CH4=dict(erf_w_m2=0.54, range_w_m2=(0.43, 0.65), from_ppm=0.7292, to_ppm=1.8663),
           N2O=dict(erf_w_m2=0.21, range_w_m2=(0.18, 0.24), from_ppm=0.2701, to_ppm=0.3321))
DESIGN_DIMMING = 0.95            # the titania stack dimmed a further 5% (decisions register)
EVIDENCE = ("Line-by-line clear-sky radiative transfer in the inverse-climate columns with HITRAN CH4 and N2O lines "
            "added at fixed temperature and humidity profiles: instantaneous forcing at the top and the net change at "
            "the tropopause, an Earth control set against AR6's effective forcing, and warming through stated "
            "responses. No clouds, no stratospheric adjustment, no shortwave absorption by the gases, no chemistry.")
READING_RULE = ("Mixing ratios ppm of dry air; forcing W/m2, positive warming; toa is the drop in outgoing longwave, "
                "tropopause the rise in net downward flux there. allsky scales clear-sky by AR6's ratio on the Earth "
                "control. warming_k is forcing over a response in W/m2/K: the column's clear-sky response near the "
                "design case, and provisioning's 1.0-1.6.")


def input_records() -> dict:
    """Name -> SHA-256 of every external input, in inputs.json and inputs_added.json."""
    manifests = (fetch_inputs.INPUTS.manifest(), fetch_added_inputs.manifest())
    return {f['name']: f['sha256'] for m in manifests for f in m['files']}


def register():
    """Add CH4 and N2O to the column's spectroscopy for this calculation: their tables' entries, and their inputs read
    from inputs_added.json. Lookups of every other input go to inputs.json as before."""
    sp.LINE_FILES.update(GAS_LINES)
    sp.ISOTOPOLOGUES.update(GAS_ISOTOPOLOGUES)
    added = {f['name'] for f in fetch_added_inputs.manifest()['files']}

    def path(name):
        return fetch_added_inputs.path(name) if name in added else fetch_inputs.INPUTS.path(name)

    def manifest():
        m = fetch_inputs.INPUTS.manifest()
        return dict(m, files=m['files'] + fetch_added_inputs.manifest()['files'])
    fetch_inputs.path, fetch_inputs.manifest = path, manifest


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fluxes(col, nu, tau):
    up, down = lw.fluxes(nu, tau, col['t_k'], col['t_k'][0])
    return lw.integrate(nu, up), lw.integrate(nu, down)


def tropopause_level(col) -> int:
    return int(np.argmin(np.abs(np.log(col['p_pa'] / col['tropopause_pa']))))


def per_ppm(col, gas, cfg):
    """Optical depth (layers x wavenumbers) of 1 ppm of `gas`, well mixed in dry air: each layer's cross-section from
    the gas's (ln p, T) node table, optics.py's path for air-broadened gases, times the layer's molecules. It is
    assembled here so the optics module, which other products bind by hash, stays as it is."""
    gas_cfg = replace(cfg, molecules=(gas,), s_min=((gas, S_MIN),), continuum=False, cia=False, tabulated=(gas,))
    optics._load(gas_cfg)
    grid = gas_cfg.grid
    layers = optics.layer_inputs(col)
    needed = set()
    for layer in layers:
        i0, _, j0, _ = optics._brackets(layer['p'], layer['t'])
        needed.update((i0 + di, j0 + dj) for di in (0, 1) for dj in (0, 1))
    missing = [(grid, gas_cfg, gas, i, j) for i, j in sorted(needed)
               if not optics._node_path(gas_cfg, gas, i, j).is_file()]
    optics._map(optics._node, missing, gas_cfg)
    molecules = col['layer_column_cm2'] * (1 - col['layer_x_h2o']) * 1e-6
    tau = np.vstack([optics._tabulated(gas_cfg, gas, layer['p'], layer['t']) * molecules[k]
                     for k, layer in enumerate(layers)])
    return grid.nu, tau


def column_cases(col, cfg):
    """Forcing at every case for one column."""
    nu, tau = optics.optical_depth(col, cfg)
    _, tau_ch4 = per_ppm(col, 'CH4', cfg)
    _, tau_n2o = per_ppm(col, 'N2O', cfg)
    k = tropopause_level(col)
    up0, down0 = fluxes(col, nu, tau)

    def case(ch4, n2o):
        up, down = fluxes(col, nu, tau + ch4 * tau_ch4 + n2o * tau_n2o)
        return dict(toa=float(up0[-1] - up[-1]), tropopause=float((down[k] - up[k]) - (down0[k] - up0[k])))
    return dict(olr_w_m2=float(up0[-1]), tropopause_pa=float(col['p_pa'][k]),
                ch4={f'{x:g}': case(x, 0.0) for x in CH4_PPM},
                n2o={f'{y:g}': case(0.0, y) for y in N2O_PPM},
                both={f'{x:g}+{y:g}': case(x, y) for x, y in BOTH_PPM})


def clear_sky_response(rows, scenario, shield='titania_stack', dimming=DESIGN_DIMMING, ts=(290.0, 300.0)):
    """d(OLR - ASR)/dTs (W/m2/K) between two surface temperatures in the inverse-climate rows, with the design's
    shield and dimming on the absorbed sunlight."""
    olr, asr = {}, {}
    for r in rows:
        if r['scenario'] != scenario or float(r['ts_k']) not in ts or r['humidity'] != 'manabe_wetherald':
            continue
        if float(r['co2_ppm']) != 400.0 or float(r['stratosphere_k']) != 200.0:
            continue
        if r['olr']:
            olr[float(r['ts_k'])] = float(r['olr'])
        if r['shield'] == shield and r['asr']:
            asr[float(r['ts_k'])] = dimming * float(r['asr'])
    a, b = ts
    return ((olr[b] - asr[b]) - (olr[a] - asr[a])) / (b - a), dict(olr=olr, asr=asr)


def results(processes=6):
    register()
    cfg = replace(optics.LONGWAVE, processes=processes)
    columns, seconds = {}, {}
    for name, sc in SCENARIOS.items():
        columns[name] = {}
        for ts in TS:
            start = time.time()
            columns[name][f'{ts:g}'] = column_cases(sc.column(ts, sc.longwave_levels), cfg)
            seconds[f'{name} {ts:g}'] = round(time.time() - start, 1)
            print(name, ts, seconds[f'{name} {ts:g}'], 's', flush=True)
    rows = list(csv.DictReader(INVERSE.open()))
    response, used = clear_sky_response(rows, 'moon_1.2atm')
    heat = json.loads(POPULATION.read_text())['heat']
    lam = [tw * heat['moon_w_m2_per_tw'] for tw in heat['moon_tw_per_k']]
    # AR6 on the Earth control: its effective forcing over this clear-sky top-of-air forcing at 295 K.
    earth = columns['earth_1atm'][f'{DESIGN_TS:g}']
    ratios = {}
    for gas, ar6 in AR6.items():
        key = gas.lower()
        clear = earth[key][f"{ar6['to_ppm']:g}"]['toa'] - earth[key][f"{ar6['from_ppm']:g}"]['toa']
        ratios[gas] = dict(clear_sky_toa_w_m2=clear, ar6_erf_w_m2=ar6['erf_w_m2'],
                           ar6_range_w_m2=list(ar6['range_w_m2']),
                           ratio=ar6['erf_w_m2'] / clear)
    moon = columns['moon_12atm'][f'{DESIGN_TS:g}']
    warming = {}
    for gas, cases in (('CH4', moon['ch4']), ('N2O', moon['n2o'])):
        r = ratios[gas]['ratio']
        warming[gas] = {x: dict(clear_sky_toa_w_m2=c['toa'], allsky_w_m2=r * c['toa'],
                                warming_k_clear_sky_response=c['toa'] / response,
                                warming_k_heat_block=[r * c['toa'] / lam[1], r * c['toa'] / lam[0]])
                        for x, c in cases.items()}
    both_ratio = 0.5 * (ratios['CH4']['ratio'] + ratios['N2O']['ratio'])
    warming['both'] = {x: dict(clear_sky_toa_w_m2=c['toa'], allsky_w_m2=both_ratio * c['toa'],
                               warming_k_clear_sky_response=c['toa'] / response,
                               warming_k_heat_block=[both_ratio * c['toa'] / lam[1], both_ratio * c['toa'] / lam[0]],
                               sum_of_parts_toa_w_m2=(moon['ch4'][x.split('+')[0]]['toa']
                                                      + moon['n2o'][x.split('+')[1]]['toa']))
                       for x, c in moon['both'].items()}
    manifest = input_records()
    return dict(
        schema=SCHEMA,
        producer=dict(files={f'atmosphere/radiative_convective/{c}': digest(HERE / c) for c in CODE},
                      inputs={**{f'inputs/{n}': manifest[n] for n in LINES},
                              'atmosphere/radiative_convective/results/inverse_climate.csv': digest(INVERSE),
                              'provisioning/results/population.json': digest(POPULATION)},
                      seconds=seconds, processes=processes),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        mixing_ratios_ppm=dict(ch4=list(CH4_PPM), n2o=list(N2O_PPM), both=[list(x) for x in BOTH_PPM]),
        surface_temperatures_k=list(TS), design_ts_k=DESIGN_TS, line_intensity_floor=S_MIN,
        columns=columns, ar6_check=ratios,
        responses=dict(clear_sky_w_m2_k=response, clear_sky_basis=dict(scenario='moon_1.2atm', shield='titania_stack',
                                                                       dimming=DESIGN_DIMMING, rows=used),
                       heat_block_w_m2_k=lam),
        moon_at_design=warming)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--processes', type=int, default=6)
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args(argv)
    product = results(args.processes)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(product, indent=1) + '\n')
    print('wrote', args.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
