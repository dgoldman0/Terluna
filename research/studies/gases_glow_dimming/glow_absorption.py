"""Where infrared from space is absorbed in the Moon's air: the fleet's glow and the Moon's own returned infrared,
by the atmosphere domain's correlated-k model, beside the GCM's broadband treatment of them.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.gases_glow_dimming.glow_absorption   # about half a minute

The GCM carries the glow as infrared added at the top of the air and taken down by PlaSim's broadband
transmissivities, which are fitted to the air's own emission. The fleet's films are cold, 150-190 K on their faces and
up to about 300 K inside the deepest stacks (the glow study, research/array-industry), so their glow falls mostly in the
far infrared and the CO2 band. The infrared they return is the Moon's own outgoing infrared, which comes back down the
paths it left by. This check takes both through the design column of the middle atmosphere (its clear titania-stack
columns at 288 and 298 K, which the glow study also uses for the Moon's outgoing infrared), band by band and g-point
by g-point:
- the glow as grey emission at 150-300 K, isotropic from the sky;
- the returned infrared at the angle and g-point at which it left, as a specular film returns it.
It gives the share of each that reaches the ground and the share absorbed in each of PlaSim's ten layers.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'glow_absorption.json'
SCHEMA = 'terluna.research.glow-absorption/1'
MIDDLE = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results'
COLUMNS = ('moon_1.2atm_titania_stack', 'moon_1.2atm_titania_stack_ts298')
DESIGN_SURFACE_K = 294.9                     # the GCM design case, A28_dim5_moon years 20-29
GLOW_K = (150.0, 170.0, 190.0, 230.0, 270.0, 300.0)
# PlaSim's half levels at T21 with 10 layers, from its full levels 0.038 ... 0.983 (each full level the mean of the
# half levels about it), top first.
PLASIM_SIGMA_HALF = (0.0, 0.076, 0.162, 0.260, 0.374, 0.500, 0.634, 0.764, 0.882, 0.966, 1.0)
EVIDENCE = ('The atmosphere domain\'s correlated-k thermal radiation (16 bands, 16 g-points, H2O lines and MT_CKD '
            'continuum, CO2, O3 and collision-induced absorption, checked against its line-by-line model) on the middle '
            'atmosphere\'s clear titania-stack columns at 1.2 atm. Clear sky: clouds would stop more of both before the '
            'ground. The glow is taken as grey and isotropic over a range of film temperatures; the films\' own '
            'spectral emissivity and the fleet\'s place in the sky are not followed.')
READING_RULE = ('Shares are of the infrared arriving at the top of the air. ground is the share reaching the surface; '
                'by_plasim_layer[k] the share absorbed in PlaSim\'s layer k, top first, between PLASIM_SIGMA_HALF[k] and '
                '[k+1] of the surface pressure; above_tropopause the share absorbed above the column\'s tropopause. '
                'at_design interpolates the two columns linearly to the design surface temperature.')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rnd(x, sig=4):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def input_files():
    return [MIDDLE / 'cases' / f'{n}.json' for n in COLUMNS] + [MIDDLE / 'profiles' / f'{n}.csv' for n in COLUMNS]


def shares(left, sigma, trop):
    """From the share still travelling down at each level (surface first), the shares reaching the ground and
    absorbed in each PlaSim layer and above the tropopause."""
    absorbed = left[1:] - left[:-1]                      # absorbed in each layer of the column, surface first
    lower, upper = sigma[:-1], sigma[1:]                 # each layer's sigma at its bottom and top
    layers = []
    for top, bottom in zip(PLASIM_SIGMA_HALF[:-1], PLASIM_SIGMA_HALF[1:]):
        # each column layer's absorption shared among PlaSim's layers by the air they hold in common
        overlap = np.clip(np.minimum(lower, bottom) - np.maximum(upper, top), 0.0, None) / (lower - upper)
        layers.append(float((absorbed * overlap).sum()))
    return dict(ground=float(left[0]), by_plasim_layer=layers, above_tropopause=float(absorbed[trop:].sum()))


def column(name, model, plan):
    from atmosphere.radiative_convective import ck, longwave as lw
    from atmosphere.middle_atmosphere import lbl_check
    case = plan[name]
    _, t, trop, o3, col, _ = lbl_check.rebuild(case, MIDDLE)
    tau = model.optical_depth(col, {'O3': o3 * col['layer_column_cm2']})       # (layers, bands, g), surface first
    nl = tau.shape[0]
    depth = np.zeros((nl + 1, ck.NB, ck.NG))                                  # from the top down to each level
    for j in range(nl - 1, -1, -1):
        depth[j] = depth[j + 1] + tau[j]
    mu, wmu = lw.angles(4)
    # the outgoing radiance at the top by angle, band and g-point: longwave's linear-source upward sweep
    b = ck.band_planck(col['t_k'])
    up = []
    for m in mu:
        i = np.broadcast_to(ck.band_planck(np.array([col['t_k'][0]]))[0][:, None], (ck.NB, ck.NG)).copy()
        for k in range(nl):
            e, a, c = lw._linear_source_terms(tau[k] / m)
            i = i * e + b[k + 1][:, None] * a - (b[k + 1] - b[k])[:, None] * c
        up.append(i * ck.W)
    weight = [2 * math.pi * wt * m for m, wt in zip(mu, wmu)]
    olr = sum(w * u for w, u in zip(weight, up))
    returned = sum(w * u[None] * np.exp(-depth / m) for w, u, m in zip(weight, up, mu)).sum(axis=(1, 2)) / olr.sum()
    sigma = col['p_pa'] / col['p_pa'][0]
    out = dict(surface_k=float(case.surface_temperature_k), olr_w_m2=float(olr.sum()),
               returned_infrared=shares(returned, sigma, trop), glow={})
    for temp in GLOW_K:
        src = ck.band_planck(np.array([temp]))[0][:, None] * ck.W[None]
        down = sum(w * src[None] * np.exp(-depth / m) for w, m in zip(weight, mu)).sum(axis=(1, 2))
        out['glow'][f'{temp:.0f}_k'] = shares(down / (math.pi * src.sum()), sigma, trop)
    return out


def interpolate(a, b, x):
    if isinstance(a, dict):
        return {k: interpolate(a[k], b[k], x) for k in a}
    if isinstance(a, list):
        return [interpolate(p, q, x) for p, q in zip(a, b)]
    return a + (b - a) * x


def results():
    from atmosphere.radiative_convective import ck
    from atmosphere.middle_atmosphere import run as rn
    plan = rn.cases()
    model = ck.CKLongwave()
    cols = {name: column(name, model, plan) for name in COLUMNS}
    c0, c1 = (cols[n] for n in COLUMNS)
    x = (DESIGN_SURFACE_K - c0['surface_k']) / (c1['surface_k'] - c0['surface_k'])
    at_design = {k: interpolate(c0[k], c1[k], x) for k in ('returned_infrared', 'glow')}
    return dict(
        schema=SCHEMA,
        producer=dict(files={str(Path(__file__).resolve().relative_to(ROOT)): sha256(__file__)},
                      inputs={str(p.relative_to(ROOT)): sha256(p) for p in input_files()},
                      ck_table_key=ck.table_key()),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        design=dict(columns=list(COLUMNS), design_surface_k=DESIGN_SURFACE_K, glow_k=list(GLOW_K),
                    plasim_sigma_half=list(PLASIM_SIGMA_HALF)),
        columns=rnd(cols), at_design=rnd(at_design))


def dumps(product) -> str:
    return json.dumps(product, indent=2, allow_nan=False) + '\n'


def main(argv=None) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()))
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
