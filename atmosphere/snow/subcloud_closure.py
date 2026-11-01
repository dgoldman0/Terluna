"""PlaSim's snow sublimation below a cloud base from McSnow's capped lunar flakes: the adopted closure's tables and its
reference implementation (mcsnow_subcloud.md, the fifth test; register, 30 October 2026).

Below a cloud base at height h1 above the freezing level the snow's mass flux F falls, going down, as

    dF/dz = s(z) k(z) [(F / F(h1)) / (F^(z) / F^(h1))]^(-q(z)) F,

with s the subsaturation over ice of the air the snow meets, k the share of the flux the flakes sublimate per metre at
unit subsaturation, q how k rises as the flakes shrink, and F^ the flux line of the matched references. k and q come
from the flake law (capacitance 0.25 D, Hall and Pruppacher's ventilation, Heymsfield and Westbrook's fall at lunar
gravity, Brown and Francis's mass law; aggregation.py) applied to the superparticles of McSnow reference columns by
0.5 km layer over each column's last half: k with McSnow's own fall speeds in the flux, q from the flakes' masses scaled
by 0.7 and their sizes by 0.7^(1/1.9). The references form two families, S from -15 C (16.28 km above the freezing
level) and D from -30 C (32.56 km), each at 0.1, 0.3, 1 and 3 mm/h at the reference crystal flux, every one on the
capped build (aggregates within 14 cm; mcsnow_subcloud.py).

A column's snow takes each family's references matched by ln F(h1) against their fluxes at h1, linear between
neighbours and held at the ends, and ln k, q and ln F^ by height. A snow top H between 16.28 and 32.56 km interpolates
the families at the same height with weight (H - 16.28 km) / (32.56 km - 16.28 km) on D, each family by depth below its
own top; a top at or below 16.28 km takes family S, and one at or above 32.56 km family D, by depth below H, held at
the deepest layer tabulated. F^ is the least-squares line of a reference's flux on height through its layers from the
freezing level up to the one holding h1 + 500 m.

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.subcloud_closure      # writes results/subcloud_closure.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from atmosphere.snow import aggregation as agg
from atmosphere.snow import flakes as sf
from atmosphere.snow import mcsnow_growth_aggregation as growagg
from atmosphere.snow import mcsnow_subcloud as sc

HERE = Path(__file__).resolve().parent
PRODUCT = HERE / 'results/subcloud_closure.json'
SCHEMA = 'terluna.atmosphere.snow-subcloud-closure/1'
RUNS = sc.CACHE / 'runs'
LAYER = 500.0                     # m, the tables' layers
STEP = 10.0                       # m, the reference implementation's step
R_SHRINK = 0.7                    # the mass ratio that gives q
TOP_S, TOP_D = 16280.0, 32560.0   # m above the freezing level
RATES = ('q0.1', 'q0.3', 'q1', 'q3')
FAMILIES = {'S': (tuple(f'cap_ref_{r}_f1' for r in RATES), TOP_S),
            'D': (tuple(f'cap_rdeep_{r}_f1' for r in RATES), TOP_D)}
SENSITIVITY_10CM = {'cap_ref_q1_f1': 'capl_ref_q1_f1', 'cap_rdeep_q0.3_f1': 'capl_rdeep_q0.3_f1',
                    'cap_rdeep_q1_f1': 'capl_rdeep_q1_f1'}


# --- the tables, from the references' snapshots --------------------------------------------------------------------

def read(name):
    """A run's snapshots over its last half, each an array of height, mass, size, monomers, multiplicity and fall speed,
    with unset fall speeds and seeds below twice the seed mass left out."""
    rec = json.loads((RUNS / name / 'run.json').read_text())
    end = rec['configuration']['general_nml']['time_end']
    raw = np.fromfile(RUNS / name / 'snapshots.bin', dtype='<f8')
    snaps, i = [], 0
    while i < raw.size:
        t, n = raw[i], int(raw[i + 1])
        a = raw[i + 2:i + 2 + 6 * n].reshape(n, 6)
        i += 2 + 6 * n
        if t >= end / 2.0 - 1.0:
            v, m = a[:, 5], a[:, 1]
            ok = np.isfinite(v) & (v > 1e-6) & (v < 50.0) & (m > 2.0 * growagg.SEED_MASS)
            snaps.append(a[ok])
    grid = rec['configuration']['grid_nml']
    return dict(name=name, snaps=snaps, top=float(grid['dom_top']), area=float(grid['box_area']),
                executable_sha256=rec['executable_sha256'])


def layer_k(a, xi, s, column, zc, r=1.0, own_speed=True):
    """The flake law's k over the superparticles s of a pooled array: their deposition at unit subsaturation over the
    flux they carry, masses scaled by r and sizes by r^(1/1.9); the flux with McSnow's fall speeds (own_speed) or the
    flake law's."""
    st = column.state(column.z_freeze + zc)
    unit = dict(st, e=2.0 * sf.ice_saturation_pressure(st['t']))
    m, d = a[s, 1] * r, a[s, 2] * r ** (1.0 / sf.MASS_EXPONENT)
    v, re = agg.speeds(d, m, column.gravity, st)
    g = agg.deposition(d, re, unit)
    return float((xi[s] * g).sum() / (xi[s] * m * (a[s, 5] if own_speed else v)).sum())


def profile(run, column):
    """k, q and the flux (kg/m2/s) by 0.5 km layer from the freezing level, over the pooled snapshots."""
    a = np.concatenate(run['snaps'])
    xi = a[:, 4] / len(run['snaps'])
    rows = []
    for lo in np.arange(0.0, run['top'], LAYER):
        zc = lo + LAYER / 2.0
        s = (a[:, 0] >= lo) & (a[:, 0] < lo + LAYER)
        if s.sum() < 5:
            continue
        q = -math.log(layer_k(a, xi, s, column, zc, R_SHRINK, False) / layer_k(a, xi, s, column, zc, 1.0, False)) \
            / math.log(R_SHRINK)
        rows.append(dict(z_m=float(zc), depth_m=float(run['top'] - zc), k_per_m=layer_k(a, xi, s, column, zc), q=q,
                         flux_kg_m2_s=float((xi[s] * a[s, 1] * a[s, 5]).sum() / (LAYER * run['area']))))
    return rows


# --- the closure, from the tables ------------------------------------------------------------------------------------

class Reference:
    """One reference column's table: k and q by depth below its top, its flux by height, its flux line."""

    def __init__(self, table, top):
        self.top = top
        rows = sorted(table, key=lambda x: x['z_m'])
        self.z = np.array([x['z_m'] for x in rows])
        self.k = np.array([x['k_per_m'] for x in rows])
        self.q = np.array([x['q'] for x in rows])
        self.f = np.array([x['flux_kg_m2_s'] for x in rows])

    def clamp(self, zr):
        """A height in the reference's own column, held within its lowest and highest layers' centres."""
        return min(max(zr, float(self.z[0])), float(self.z[-1]))

    def at(self, zr):
        """ln k and q at height zr in the reference's own column."""
        zr = self.clamp(zr)
        return math.log(float(np.interp(zr, self.z, self.k))), float(np.interp(zr, self.z, self.q))

    def flux(self, zr):
        return float(np.interp(self.clamp(zr), self.z, self.f))

    def line(self, zr1):
        """The least-squares line (a, b) of flux on height through the layers from the freezing level up to the one
        holding zr1 + 500 m (two layers at least)."""
        sel = self.z - LAYER / 2.0 <= max(self.clamp(zr1), LAYER / 2.0) + LAYER
        b, a = np.polyfit(self.z[sel], self.f[sel], 1)
        return float(a), float(b)


def load(product=None, sensitivity=None):
    """The families' references from the product (or the file); sensitivity '10cm' swaps in the 10 cm references."""
    product = product or json.loads(PRODUCT.read_text())
    out = {}
    for fam, spec in product['families'].items():
        refs = []
        for name in spec['references']:
            use = product['sensitivity_10cm']['swap'].get(name, name) if sensitivity == '10cm' else name
            table = product['tables'][use]
            refs.append(Reference(table, spec['top_m']))
        out[fam] = dict(refs=refs, top=spec['top_m'])
    return out


def match(x, nodes):
    """Weights over the references from ln F against their ln fluxes: linear between neighbours, held outside."""
    nodes = np.asarray(nodes)
    order = np.argsort(nodes)
    xs = nodes[order]
    w = np.zeros(len(nodes))
    if x <= xs[0]:
        w[order[0]] = 1.0
    elif x >= xs[-1]:
        w[order[-1]] = 1.0
    else:
        j = int(np.searchsorted(xs, x)) - 1
        f = (x - xs[j]) / (xs[j + 1] - xs[j])
        w[order[j]], w[order[j + 1]] = 1.0 - f, f
    return w


def family_terms(fam, h1, f_match, z_match, shift):
    """ln k, q and ln F^ of one family as functions of height z, its references matched by ln f_match against their
    fluxes at z_match; shift moves heights into the references' own columns (z_ref = z + shift, held within them)."""
    refs = fam['refs']
    w = match(math.log(f_match), [math.log(r.flux(z_match + shift)) for r in refs])
    used = [(wi, r, r.line(h1 + shift)) for wi, r in zip(w, refs) if wi > 0.0]
    lnk = lambda z: sum(wi * r.at(z + shift)[0] for wi, r, _ in used)
    q = lambda z: sum(wi * r.at(z + shift)[1] for wi, r, _ in used)
    lnf = lambda z: sum(wi * math.log(max(a + b * r.clamp(z + shift), 1e-30)) for wi, r, (a, b) in used)
    return (lnk, q, lnf), [float(x) for x in w]


def terms(families, h1, f_match, top, z_match=None):
    """The closure's ln k, q and ln F^ by height for a snow top at `top` (m above the freezing level), matched by
    ln f_match at z_match (h1 unless given), with the families' weights and the weight on D."""
    z_match = h1 if z_match is None else z_match
    s_top, d_top = families['S']['top'], families['D']['top']
    if top <= s_top:
        t, w = family_terms(families['S'], h1, f_match, z_match, s_top - top)
        return t, dict(S=w), 0.0
    if top >= d_top:
        t, w = family_terms(families['D'], h1, f_match, z_match, d_top - top)
        return t, dict(D=w), 1.0
    wd = (top - s_top) / (d_top - s_top)
    ts, ws = family_terms(families['S'], h1, f_match, z_match, 0.0)
    td, wdd = family_terms(families['D'], h1, f_match, z_match, 0.0)
    mix = tuple((lambda z, a=a, b=b: (1.0 - wd) * a(z) + wd * b(z)) for a, b in zip(ts, td))
    return mix, dict(S=ws, D=wdd), wd


def closure(families, h1, f_h1, top, subsaturation, z_low=0.0, step=STEP, edges=None, base_layer=False):
    """F from h1 down to z_low (m above the freezing level) for snow arriving at h1 with flux f_h1 (kg/m2/s) from a
    top at `top`, in air whose subsaturation over ice at height z is subsaturation(z). Returns D = ln(F(h1)/F(z_low)),
    F(z_low), F at each of `edges` passed on the way, and the families' weights. With base_layer, f_h1 is the flux in
    the 500 m layer above h1, as McSnow's tests measured it: the references are matched there, the flux is carried down
    their line to h1 first, and D runs from that layer."""
    z_match = h1 + LAYER / 2.0 if base_layer else h1
    (lnk, q, lnf), weights, wd = terms(families, h1, f_h1, top, z_match)
    f_start = f_h1
    if base_layer:
        f_h1 = f_start * math.exp(lnf(h1) - lnf(z_match))
    fhat1 = math.exp(lnf(h1))
    f = f_h1
    grid = list(np.arange(h1, z_low, -step)) + [z_low]
    marks = sorted(edges or [], reverse=True)
    at = {}
    for z_hi, z_lo in zip(grid[:-1], grid[1:]):
        zm = 0.5 * (z_hi + z_lo)
        if f > 0.0:
            ratio = (f / f_h1) / (math.exp(lnf(zm)) / fhat1)
            f *= math.exp(-subsaturation(zm) * math.exp(lnk(zm)) * ratio ** (-q(zm)) * (z_hi - z_lo))
        while marks and marks[0] >= z_lo:
            at[marks.pop(0)] = f
    d = math.log(f_start / f) if f > 0.0 else math.inf
    return dict(D=d, f_low=f, f_h1=f_h1, at=at, weights=weights, weight_on_D=wd)


# --- the product ------------------------------------------------------------------------------------------------------

def test_vectors(product):
    """The closure on the fifth test's resolved box runs' settings, from the tables (the base layer's flux, constant
    subsaturation down to the low layer's centre), for the climate model's copy to reproduce."""
    families = load(product)
    out = []
    for case in product['inputs']['test_cases']:
        r = closure(families, case['h1_m'], case['f_base_kg_m2_s'], case['top_m'], lambda z, s=case['subsaturation']: s,
                    case['z_low_m'], base_layer=True)
        out.append(dict(case, D=r['D'], weight_on_D=r['weight_on_D'], weights=r['weights']))
    return out


FIFTH_TEST_CASES = (   # name, box, top (m), z_low (m), F_base (mm/h), D by the fifth test's closure (from snapshots)
    ('cap_sub_sea_c', 'sea', TOP_S, 250.0, 0.6128, 0.3579), ('cap_sub_land_c', 'land', TOP_S, 250.0, 0.4274, 0.6742),
    ('cap_sub_sea_l', 'sea', TOP_S, 250.0, 0.2050, 0.9069), ('cap_sub_land_l', 'land', TOP_S, 250.0, 0.1411, 2.3014),
    ('cap_sdeep_land_q1', 'land', TOP_D, 250.0, 0.6577, 0.3846),
    ('cap_sdeep_sea_q0.1', 'sea', TOP_D, 250.0, 0.0846, 1.0067),
    ('cap_sdeep_land_q0.1', 'land', TOP_D, 250.0, 0.0730, 1.9256),
    ('cap_sub_sea_q0.1', 'sea', TOP_S, 750.0, 0.0722, 4.0659),
    ('cap_sub_land_q0.1', 'land', TOP_S, 5250.0, 0.0477, 4.6318),
    ('cap_sub_sea_dry', 'sea_dry', TOP_S, 250.0, 0.6275, 0.6715),
    ('cap_sdeep_land_dry', 'land_dry', TOP_D, 250.0, 0.6173, 1.1937))


def run():
    column, binding = sf.moon_column()
    names = [n for fam in FAMILIES.values() for n in fam[0]] + list(SENSITIVITY_10CM.values())
    tables, runs = {}, {}
    for name in names:
        r = read(name)
        tables[name] = profile(r, column)
        runs[name] = dict(snapshots=len(r['snaps']), top_m=r['top'], executable_sha256=r['executable_sha256'],
                          run_json_sha256=hashlib.sha256((RUNS / name / 'run.json').read_bytes()).hexdigest())
    cases = []
    for name, box, top, z_low, f_mm_h, d_fifth in FIFTH_TEST_CASES:
        b = sc.BOXES[box]
        cases.append(dict(name=name, h1_m=b['h1_m'], subsaturation=-b['ssat_below'], top_m=top, z_low_m=z_low,
                          f_base_kg_m2_s=f_mm_h / 3600.0, D_fifth_test=d_fifth))
    product = dict(
        schema=SCHEMA,
        evidence=('Tables from McSnow reference columns on the capped lunar build, by the flake law over their '
                  'superparticles; the closure they define followed McSnow within its bands on nine of the eleven '
                  'resolved box runs of the fifth test and overstated the loss in 0.1 mm/h snow from -15 C (1.48, '
                  '1.65 in e-foldings). No lunar observation; the cap stands in for a breakup law.'),
        reading_rule=('Below a cloud base h1 (m above the freezing level), dF/dz = s k [(F/F(h1))/(F^(z)/F^(h1))]^-q '
                      'F going down, s the subsaturation over ice. Each family matches its references by ln F(h1) '
                      'against their fluxes at h1 (linear between neighbours, held at the ends) and averages ln k, q '
                      'and ln F^ with those weights; F^ is the least-squares line of a reference flux on height '
                      'through its layers from 0 up to the one holding h1 + 500 m. A snow top H between the '
                      'families\' tops interpolates the two at the same height with weight (H - top_S)/(top_D - '
                      'top_S) on D, each by depth below its own top; H at or below top_S takes S, and at or above '
                      'top_D takes D, by depth below H, held at the deepest layer. Heights in a family are its '
                      'references\' own: z_ref = z + top_f - H for the single-family cases. k is per metre at unit '
                      'subsaturation.'),
        assumptions=['k and q carry the thermodynamics of the 301 K design ladder at each height; a column warmer or '
                     'colder than the ladder takes them as they are.',
                     'A snow top above 32.56 km takes family D by depth below the top, held at the deepest layer, '
                     'where the capped flakes have gathered under the cap; no run tested such a top.',
                     'In 0.1 mm/h snow from -15 C, where nearly all the snow sublimates, the closure overstates the '
                     'loss by about five points of the flux at the cloud base: the stated limit of its adoption.'],
        families={f: dict(references=list(refs), top_m=top) for f, (refs, top) in FAMILIES.items()},
        sensitivity_10cm=dict(swap=SENSITIVITY_10CM,
                              note='The 10 cm references in place of their 14 cm twins; the others unchanged.'),
        tables=tables, runs=runs,
        inputs=dict(layer_m=LAYER, step_m=STEP, shrink_ratio=R_SHRINK, test_cases=cases,
                    column=dict(gravity=column.gravity, z_freeze_m=column.z_freeze, lapse_k_per_m=column.lapse),
                    climate_binding=binding.get('inputs', {}),
                    files_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in (Path(agg.__file__), Path(sf.__file__))}),
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    product['test_vectors'] = test_vectors(product)
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    with sc.subcloud():
        product = run()
    for v in product['test_vectors']:
        print(f"{v['name']}: D from the tables {v['D']:.4f}, the fifth test's {v['D_fifth_test']:.4f} "
              f"({v['D'] / v['D_fifth_test']:.3f}); weight on D {v['weight_on_D']:.2f}")


if __name__ == '__main__':
    main()
