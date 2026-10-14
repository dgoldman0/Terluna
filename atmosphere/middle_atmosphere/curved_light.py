"""The Moon's curved air in the middle atmosphere: what its longer-lit upper air does to the exobase, and how far the
GCM's misplaced sunlight would move the ground.

    python -m atmosphere.middle_atmosphere.curved_light       # results/curved_light.json, about an hour on 4 cores

The middle-atmosphere column (equilibrium.py) takes the Sun as a global and diurnal mean through a pseudo-spherical
two-stream column, as flat a footing as the GCM's. climate/results/gcm/twilight.json gives, shell by shell, the
sunlight the curved air puts into the design column's air per unit mass against that flat column (its shells_w_m2,
sphere_shallow over pseudo_spherical): about 5% more in the troposphere, 14-15% more from 100 to 360 km and about 20%
more at the thermosphere's 0.3 Pa base, where the air stays lit for long after the ground below has lost the Sun.
Two checks:

- exobase: the three design cases behind the exobase table (cooling in local thermodynamic equilibrium; CO2's 15 um
  band out of it with all the absorbed near-infrared sunlight turned to heat; and with only its collisional share)
  solved again with their near-infrared heating (500 nm-5 um, which holds all the upper air's solar heating behind
  the titania stack) scaled by that ratio, and the thermal column (escape.py) run from the new 0.3 Pa base with the
  heats the exobase table uses;
- ground: the design case at 298 K solved again with its near-infrared heating shared among the GCM's ten layers as
  PlaSim shares its absorbed sunlight (twilight.json's PlaSim against the line-by-line column, the column total kept),
  at the same fixed surface temperature, and again as PlaSim shares it with its CO2 and O2 bands and the Moon's slant
  factors, at the water and Rayleigh multipliers refitted with them (climate/gcm/plasim_shortwave.py; the GCM's
  A28_dim5_moon_geometry_gases). The change in the top-of-air imbalance over the column's own change of imbalance
  with surface temperature is how far the ground would move, the column's convection and humidity as the model holds
  them.

Each solve is kept on the drive, so a stopped run resumes.
"""
from __future__ import annotations

import contextlib
import hashlib
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from atmosphere.radiative_convective import optics, thermodynamics as th
from atmosphere.middle_atmosphere import chemistry as ch, equilibrium as eq, escape, run as runner

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'curved_light.json'
TWILIGHT = ROOT / 'climate' / 'results' / 'gcm' / 'twilight.json'
LUNAR_FIT = ROOT / 'climate' / 'results' / 'gcm' / 'shortwave_fit_A28_dim5_lunar.json'
GASES = ROOT / 'climate' / 'results' / 'gcm' / 'shortwave_gases.json'
STORE = Path('/media/projectspace/terluna-research/gcm-runs/twilight/middle_atmosphere')
SCHEMA = 'terluna.atmosphere.curved-light/1'
EXOBASE_CASES = ('moon_1.2atm_titania_stack', 'moon_1.2atm_titania_stack_nonlte',
                 'moon_1.2atm_titania_stack_nonlte_nir_collisional')
GROUND_CASE = 'moon_1.2atm_titania_stack_ts298'
BALANCE_CASES = ('moon_1.2atm_titania_stack_ts278', 'moon_1.2atm_titania_stack', 'moon_1.2atm_titania_stack_ts298',
                 'moon_1.2atm_titania_stack_ts308')
BASE_PA = 0.3
PROCESSES = 4


def log(msg, t0=[time.time()]):
    print(f'{time.time() - t0[0]:8.1f} s  {msg}', file=sys.stderr, flush=True)


def twilight():
    return json.loads(TWILIGHT.read_text())


def curved_ratio(product):
    """The curved light's ratio of heating per unit mass, by the log of each shell's mean pressure."""
    s = product['shells_w_m2']
    p = np.array(s['pressure_edges_pa'])
    mid = np.sqrt(p[1:] * p[:-1])
    ratio = np.array(s['sphere_shallow']) / np.array(s['pseudo_spherical'])
    return np.log(mid)[::-1], ratio[::-1]                          # rising log pressure


def plasim_shares(product):
    """PlaSim's share of the column's absorbed sunlight in each GCM layer against the line-by-line column's, with the
    layers' upper edges in sigma (top first)."""
    whole = product['global_mean_w_m2']['whole_spectrum']
    plasim, lbl = np.array(whole['plasim']['air_by_layer']), np.array(whole['pseudo_spherical']['air_by_layer'])
    sigma = np.array(product['gcm']['sigma'])
    half = np.r_[0.0, [2 * s for s in sigma[:1]]]                  # PlaSim's half levels from its full levels
    for s in sigma[1:]:
        half = np.r_[half, 2 * s - half[-1]]
    return half, (plasim / plasim.sum()) / (lbl / lbl.sum())


def plasim_band_shares(product):
    """As plasim_shares, with PlaSim's CO2 and O2 bands and the Moon's slant factors (plasim_shortwave.py's lunar
    column) at the water and Rayleigh multipliers refitted with them, to three decimals as the GCM test sets them, on
    the ledger's column and solar elevations; and its global-mean absorption by layer, W/m2 of the ledger's sunlight."""
    from illumination.surface_light import model as surface
    from climate.gcm import plasim_shortwave as ps, twilight as tw
    col, air, _, _ = surface.build(surface.COLUMNS[tw.WORLD])
    refit = json.loads(LUNAR_FIT.read_text())['forms']['scale']
    k, rayleigh = round(refit['k'], 3), round(refit['rayleigh_scale'], 3)
    lunar = ps.lunar_column(tw.plasim_column(col, air.molar_mass, tw.SKY_ALBEDO), json.loads(GASES.read_text()), product)
    column = dict(lunar, rcoeff=ps.RCOEFF_PER_RAYSCALE * rayleigh)
    grid = np.radians(product['by_elevation']['elevation_deg'])
    layers = np.zeros((len(column['sigma']), grid.size))
    for i, el in enumerate(grid):
        mu = math.sin(el)
        if mu > 0.0:
            net = ps.net_flux(column, mu, 'scale', k)
            layers[:, i] = mu * (net[:-1] - net[1:])
    plasim = tw.sphere_share(layers.T, grid)
    lbl = np.array(product['global_mean_w_m2']['whole_spectrum']['pseudo_spherical']['air_by_layer'])
    return (plasim / plasim.sum()) / (lbl / lbl.sum()), plasim * product['sunlight']['covered_w_m2'] / 4.0


@contextlib.contextmanager
def heating_hook(factor, keep_total=False):
    """While inside, the solver's near-infrared heating of each layer is multiplied by factor(col) (per layer); with
    keep_total the column's total is held at what it was."""
    original = eq.near_infrared

    def hooked(col, *args, **kwargs):
        nir = original(col, *args, **kwargs)
        f = factor(col)
        layers = nir['layers'] * f
        scale = float(nir['layers'].sum() / layers.sum()) if keep_total else 1.0
        nir['layers'] = layers * scale
        nir['nodes'] = nir['nodes'] * (f * scale)[:, None]
        return nir

    eq.near_infrared = hooked
    try:
        yield
    finally:
        eq.near_infrared = original


def curved_factor(product):
    lp, ratio = curved_ratio(product)
    return lambda col: np.interp(np.log(col['layer_p_pa']), lp, ratio)


def plasim_factor(product, share=None):
    """The factor sharing the heating as PlaSim does, by its own layers' shares or by share."""
    half, own = plasim_shares(product)
    share = own if share is None else share

    def factor(col):
        sigma = col['layer_p_pa'] / col['p_pa'][0]
        layer = np.clip(np.searchsorted(half, sigma) - 1, 0, share.size - 1)
        return share[layer]
    return factor


def solved(label, case, factor=None, keep_total=False):
    """A case solved with the hook, kept on the drive by label: its layer profile and budget."""
    path = STORE / f'{label}.json'
    if path.exists():
        return json.loads(path.read_text())
    log(f'solving {label}')
    eq.SW_CONFIG = replace(optics.SHORTWAVE, processes=PROCESSES)
    start = time.time()
    with heating_hook(factor, keep_total) if factor else contextlib.nullcontext():
        res = eq.Solver().solve(case)
    m = ch.fixed_densities(res.col)['M']
    out = dict(label=label, case=case.name, seconds=time.time() - start,
               summary=runner.summarise(res, time.time() - start),
               p_pa=res.col['layer_p_pa'].tolist(), t_k=res.col['layer_t_k'].tolist(),
               O=(res.state.n[ch.INDEX['O']] / m).tolist(), sw_layers_w_m2=np.asarray(res.sw_layers).tolist(),
               budget={k: float(v) for k, v in res.budget.items()})
    STORE.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.part')
    tmp.write_text(json.dumps(out) + '\n')
    tmp.replace(path)
    log(f'{label}: T(0.3 Pa) {escape.base_temperature(out, BASE_PA):.1f} K, net at the top {out["budget"]["net_toa"]:.2f} W/m2')
    return out


def stored(case_name):
    """A case as run.py left it: its layer profile and summary."""
    prof = escape._profile(HERE / 'results', case_name)
    summary = json.loads((HERE / 'results' / 'cases' / f'{case_name}.json').read_text())
    return dict(label='stored', case=case_name, summary=summary, p_pa=prof['p_pa'].tolist(), t_k=prof['t_k'].tolist(),
                O=prof['O'].tolist())


def exobase(solution, summary):
    """The thermal column from the solution's 0.3 Pa base, with the exobase table's heats (escape.py)."""
    prof = {k: np.asarray(solution[k]) for k in ('p_pa', 't_k', 'O')}
    t_base = escape.base_temperature(prof, BASE_PA)
    dry = float(summary['dry_pressure_pa'])
    air = th.earthlike_air(dry, 400.0)
    o2 = air.fractions['O2'] / (air.fractions['O2'] + air.fractions['N2'])
    f_o = float(np.exp(np.interp(np.log(BASE_PA), np.log(prof['p_pa'])[::-1],
                                 np.log(np.maximum(prof['O'], 1e-30))[::-1])))
    mixing = ch.Mixing(scale=float(summary['kzz_scale']))
    eddy = lambda p, m=mixing, pt=float(summary['tropopause_pa']): m.profile(p, pt)
    heats = dict(swept=escape.HEATS, film_quiet=(escape.film_heat(),),
                 film_maximum=(escape.film_heat(escape.SOLAR_MAXIMUM, escape.XRAY_SOLAR_MAXIMUM),),
                 leak_quiet=(escape.leakage_heat(),),
                 leak_maximum=(escape.leakage_heat(activity=escape.SOLAR_MAXIMUM, xray_activity=escape.XRAY_SOLAR_MAXIMUM),))
    rows = []
    for source, qs in heats.items():
        for q in qs:
            cfg = escape.ColumnConfig(surface_pressure_pa=float(summary['surface_pressure_pa']),
                                      surface_temperature_k=float(summary['surface_temperature_k']),
                                      lower_temperature_k=t_base, lower_pressure_pa=BASE_PA, oxygen_mole_fraction=o2)
            row = dict(heat_source=source, deposited_heat_w_m2=float(q))
            try:
                s, column, _ = escape.solve_column(q, cfg)
                o = escape.atomic_oxygen_loss(column, f_o, escape.lower_boundary(cfg)[3], eddy)
                row.update(exobase_temperature_k=s['exobase_temperature_k'], exobase_radius_R=s['exobase_radius_R'],
                           molecular_loss_kg_s=s['molecular_loss_kg_s'], atomic_o_loss_eddy_kg_s=o['loss_eddy_kg_s'],
                           flags=' '.join(s['domain_flags']))
            except (ValueError, RuntimeError, FloatingPointError) as exc:
                row.update(error=str(exc))
            rows.append(row)
    return dict(base_temperature_k=t_base, atomic_o_fraction_at_base=f_o, rows=rows)


def main(argv=None) -> int:
    product = twilight()
    cases = runner.cases()
    lp, ratio = curved_ratio(product)
    out = dict(schema=SCHEMA,
               evidence=('One-dimensional, global- and diurnal-mean, clear-sky radiative-convective and photochemical '
                         'equilibrium (equilibrium.py) with its near-infrared solar heating scaled, layer by layer, by '
                         'the ratio of heating per unit mass the twilight ledger finds for the curved air against a '
                         'flat column of the same design air, or shared among the GCM\'s layers as PlaSim shares it, '
                         'without and with its CO2 and O2 bands; the '
                         'thermal column above the 0.3 Pa base as escape.py runs it. The ratio comes from a column with '
                         'the design case\'s blended profile and the ledger\'s limits (clear sky, local thermodynamic '
                         'equilibrium in its optics, a top at 0.1 Pa); the middle atmosphere keeps its own limits (no '
                         'waves, tides or day-night circulation).'),
               reading_rule=('exobase: for each design case, the 0.3 Pa base and the exobase at each deposited heat as '
                             'stored (run.py, recomputed here with the current escape.py) and with the curved light; '
                             'ground: the design case at 298 K with its stored net at the top of the air, solved again '
                             'and with PlaSim\'s sharing as run in the design case (and, under with_bands, with its CO2 '
                             'and O2 bands), and the shift of the ground, the change in net over minus the slope of net '
                             'with surface temperature across the balance cases.'),
               ratio=dict(log_pressure_pa=list(lp), heating_ratio=list(ratio)))
    out['exobase'] = {}
    for name in EXOBASE_CASES:
        base = stored(name)
        curved = solved(f'{name}_curved', cases[name], curved_factor(product))
        out['exobase'][name] = dict(stored=exobase(base, base['summary']), curved=exobase(curved, base['summary']),
                                    curved_solve_seconds=curved['seconds'])
        b, c = out['exobase'][name]['stored'], out['exobase'][name]['curved']
        log(f'{name}: base {b["base_temperature_k"]:.1f} -> {c["base_temperature_k"]:.1f} K')
    slopes = [(float(stored(n)['summary']['surface_temperature_k']), float(stored(n)['summary']['net_toa']))
              for n in BALANCE_CASES]
    ts, net = np.array(slopes).T
    slope = float(np.polyfit(ts, net, 1)[0])
    again = solved(f'{GROUND_CASE}_again', cases[GROUND_CASE])
    shared = solved(f'{GROUND_CASE}_plasim', cases[GROUND_CASE], plasim_factor(product), keep_total=True)
    d_net = shared['budget']['net_toa'] - again['budget']['net_toa']
    half, share = plasim_shares(product)
    out['ground'] = dict(case=GROUND_CASE, stored_net_toa_w_m2=float(stored(GROUND_CASE)['summary']['net_toa']),
                         solved_again_net_toa_w_m2=again['budget']['net_toa'],
                         plasim_shared_net_toa_w_m2=shared['budget']['net_toa'], change_w_m2=d_net,
                         balance=dict(surface_temperature_k=list(ts), net_toa_w_m2=list(net), slope_w_m2_k=slope),
                         ground_shift_k=-d_net / slope, gcm_half_levels_sigma=list(half), plasim_share_ratio=list(share),
                         t_k=dict(again=again['t_k'], plasim=shared['t_k'], p_pa=again['p_pa']))
    log(f'ground: net at the top {again["budget"]["net_toa"]:.2f} -> {shared["budget"]["net_toa"]:.2f} W/m2, '
        f'slope {slope:.2f} W/m2/K, ground {-d_net / slope:+.2f} K')
    band_share, band_w_m2 = plasim_band_shares(product)
    bands = solved(f'{GROUND_CASE}_plasim_bands', cases[GROUND_CASE], plasim_factor(product, band_share), keep_total=True)
    d_bands = bands['budget']['net_toa'] - again['budget']['net_toa']
    whole = product['global_mean_w_m2']['whole_spectrum']
    out['ground']['with_bands'] = dict(
        plasim_shared_net_toa_w_m2=bands['budget']['net_toa'], change_w_m2=d_bands, ground_shift_k=-d_bands / slope,
        plasim_share_ratio=list(band_share), t_k=bands['t_k'],
        air_by_layer_w_m2=dict(line_by_line=whole['pseudo_spherical']['air_by_layer'],
                               plasim_as_run=whole['plasim']['air_by_layer'], plasim_with_bands=list(band_w_m2)))
    log(f'ground with the bands: net at the top {bands["budget"]["net_toa"]:.2f} W/m2, ground {-d_bands / slope:+.2f} K')
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    out['producer'] = dict(domain='atmosphere', files={f'middle_atmosphere/{Path(m.__file__).name}': digest(m.__file__)
                                                         for m in (sys.modules[__name__], eq, escape, runner, ch)},
                           inputs={str(p.relative_to(ROOT)): digest(p) for p in (TWILIGHT, LUNAR_FIT, GASES)})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix('.tmp')
    tmp.write_text(json.dumps(out, indent=1, default=float) + '\n')
    tmp.replace(OUT)
    log(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
