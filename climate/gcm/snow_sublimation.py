"""PlaSim's falling snow sublimated below a cloud base from the lunar flakes' own sizes and fall speeds, in place of the
Earth-tuned constant gamma (register, 30 October 2026): the adopted closure, read from the snow work's product, applied
to many columns at once.

The flakes' mass flux F falls below a cloud base at height h1 above the freezing level, going down, as

    dF/dz = s(z) k(z) [(F / F(h1)) / (F^(z) / F^(h1))]^(-q(z)) F,

s the subsaturation over ice of the air the snow meets, k the share of the flux the flakes sublimate per metre at unit
subsaturation, q how k rises as they shrink, F^ the matched references' flux line (atmosphere/snow/results/
subcloud_closure.json, from McSnow's capped lunar columns; its reading rule is followed here line for line). The
product's two families, from -15 C (16.28 km above the freezing level) and -30 C (32.56 km), each match their four
references by ln F against the references' fluxes at h1; a snow top H between the families' tops interpolates them at
the same height with weight (H - 16.28 km) / (32.56 km - 16.28 km) on the deep family, and a top beyond either takes
that family by depth below H, held within its references' columns. With a precipitating fraction sigma the matching
takes the flux where the snow falls, F / sigma; the share sublimated applies to the whole flux.

`closure_columns` is the closure on columns as rows; its explicit steps follow the snow work's reference implementation
(atmosphere/snow/subcloud_closure.py), which this copy reproduces on the product's test vectors.

PlaSim's snow (nsnowsub 1, in rain_drops.py's patch of mkrain; snowprep and snowsubl here, line for line in Fortran):

- **Graupel and flakes.** PlaSim forms its precipitation as rain, which mkrain freezes in a cold layer: all that a
  cold layer forms is its new snow, split by kind, 0.75 of the convective and 0.4 of the large-scale as graupel (CM1's
  shares on the design boxes), the rest as flakes; rain arriving from a warm layer above and frozen joins the
  graupel. A layer forming snow is in cloud and sublimates nothing; its lower edge is the flakes' cloud base h1, where
  the closure starts afresh with the flux leaving it.
- **Heights.** Each layer's edges by height above the column's freezing level, the first crossing of 0 C going up from
  the ground. A column colder than 0 C to the ground takes its freezing level below the ground by the references'
  ladder lapse (0.921 K/km, in the product), so that height above it keeps standing for temperature; no test reached
  such a column. The snow's top is its formation-weighted height, the mean of its forming layers' middles weighted by
  the snow they form (nsnowtop 1), or its highest forming layer's middle (nsnowtop 2, a sensitivity).
- **The precipitating fraction** (rules_fraction2.txt). At the cloud base the flakes' flux F falls at local rates
  spread lognormally around ln f* = max(ln c + beta ln F, ln F) with spread s, the two design boxes' pooled numbers;
  the closure runs at the five Gauss-Hermite rates across that spread and the share kept is their weighted mean.
  The box-mean air overstates the boxes' sublimation by 1.19 against each column's own air, which no fraction
  carries; that limit is stated and never applied (AIR_FACTOR).
- **Graupel** sublimates by CM1's own law, Morrison's EPSG (QVI - QV) / ABI with his ventilation, density and fall
  speeds at lunar gravity, its number falling with its mass so that its size holds; the size is CM1's graupel mean
  diameter at the boxes' cloud bases for its typical local rate, set where it forms.
- **Air and limits.** The deficit is over ice (Murray's form with PlaSim's constant), the closure's steps exact and
  at most snowstep metres; in a layer the flakes and the graupel together take at most what would saturate it over
  ice in the step, with PlaSim's own latent-heat correction over ice; a warm layer melts all, as before, and melted
  graupel makes rain at the lunar snow's melted size, as the rain's patch takes it.
- **Shadows** (nsnowshad 1), from the same formation and air and never fed back: no fraction, the fixed fraction
  (the boxes' snowing share), the 10 cm references, all the snow as flakes, graupel at the sea box's and at the land
  box's sizes, and PlaSim's gamma. Their sublimation and melting, the adopted snow's by kind, the snow formed, the
  freezing level, the flakes' flux at their cloud bases and the part of it matched below the lightest reference, and
  what the layers' saturation over ice held back are PlaSim's codes 501-523, means over each output interval
  (OUTPUT_PATCHES).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from shared.constants import (CM1_DRY_AIR_GAS_CONSTANT, CM1_GRAUPEL_DENSITY, CM1_GRAUPEL_FALL_CAP,
                              CM1_GRAUPEL_FALL_COEFFICIENT, CM1_GRAUPEL_FALL_EXPONENT, CM1_GRAUPEL_MEAN_DIAMETER_CAP,
                              CM1_GRAUPEL_MEAN_DIAMETER_FLOOR, CM1_GRAVITY, CM1_ICE_VENTILATION_F1,
                              CM1_ICE_VENTILATION_F2, CM1_RAIN_DENSITY_EXPONENT, CM1_RAIN_DENSITY_REFERENCE_PRESSURE,
                              CM1_RAIN_DENSITY_REFERENCE_TEMPERATURE, CM1_VAPOUR_DIFFUSIVITY_COEFFICIENT,
                              CM1_VAPOUR_DIFFUSIVITY_EXPONENT, CM1_VISCOSITY_COEFFICIENT, CM1_VISCOSITY_SUTHERLAND,
                              MURRAY_ICE_A, MURRAY_ICE_B)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRODUCT = ROOT / 'atmosphere' / 'snow' / 'results' / 'subcloud_closure.json'
SCHEMA = 'terluna.atmosphere.snow-subcloud-closure/1'
LAYER = 500.0                      # m, the product's layers
TINY = 1.0e-30


class Family:
    """One family's references as arrays: heights (rows) by reference (columns), with each reference's prefix sums for
    its flux line."""

    def __init__(self, tables, names, top):
        self.top = float(top)
        self.names = list(names)
        rows = [sorted(tables[n], key=lambda x: x['z_m']) for n in names]
        self.z = [np.array([x['z_m'] for x in r]) for r in rows]
        self.lnk = [np.log(np.array([x['k_per_m'] for x in r])) for r in rows]
        self.q = [np.array([x['q'] for x in r]) for r in rows]
        self.f = [np.array([x['flux_kg_m2_s'] for x in r]) for r in rows]
        self.sums = [np.vstack([np.zeros((1, 5)), np.cumsum(np.column_stack([np.ones_like(z), z, z * z, f, z * f]),
                                                            axis=0)]) for z, f in zip(self.z, self.f)]

    def clamp(self, r, zr):
        return np.clip(zr, self.z[r][0], self.z[r][-1])

    def flux(self, r, zr):
        return np.interp(self.clamp(r, zr), self.z[r], self.f[r])

    def line(self, r, zr1):
        """Each column's least-squares line (a, b) of reference r's flux on height through its layers from the
        freezing level up to the one holding zr1 + 500 m (two layers at least)."""
        z = self.z[r]
        upto = np.maximum(self.clamp(r, zr1), LAYER / 2.0) + LAYER
        m = np.searchsorted(z - LAYER / 2.0, upto, side='right')       # layers whose lower edge lies at or below
        m = np.clip(m, 2, z.size)
        s = self.sums[r][m]
        n, sz, szz, sf, szf = s[:, 0], s[:, 1], s[:, 2], s[:, 3], s[:, 4]
        b = (n * szf - sz * sf) / (n * szz - sz * sz)
        return (sf - b * sz) / n, b


def load(path: Path = PRODUCT, sensitivity: str | None = None) -> dict:
    """The product's families; sensitivity '10cm' swaps in its 10 cm references."""
    product = json.loads(Path(path).read_text())
    if product['schema'] != SCHEMA:
        raise ValueError(f"expected {SCHEMA}, found {product['schema']}")
    swap = product['sensitivity_10cm']['swap'] if sensitivity == '10cm' else {}
    fams = {f: Family(product['tables'], [swap.get(n, n) for n in spec['references']], spec['top_m'])
            for f, spec in product['families'].items()}
    return dict(families=fams, product=product)


def match_weights(x, nodes):
    """Weights (columns as rows) over a family's references from ln F (x) against their ln fluxes (nodes, rows by
    reference): linear between neighbours in sorted order, held at the ends."""
    order = np.argsort(nodes, axis=1)
    xs = np.take_along_axis(nodes, order, axis=1)
    n = nodes.shape[1]
    w = np.zeros_like(nodes)
    rows = np.arange(nodes.shape[0])
    below, above = x <= xs[:, 0], x >= xs[:, -1]
    j = np.clip(np.sum(xs <= x[:, None], axis=1) - 1, 0, n - 2)
    lo, hi = xs[rows, j], xs[rows, j + 1]
    f = np.where(hi > lo, (x - lo) / np.where(hi > lo, hi - lo, 1.0), 0.0)
    f = np.where(below, 0.0, np.where(above, 1.0, f))
    w[rows, order[rows, j]] = 1.0 - f
    w[rows, order[rows, j + 1]] += f
    return w


class Terms:
    """The closure's ln k, q and ln F^ by height for each column, from the families' references matched at z_match by
    ln f_match, for snow tops `top` (m above the freezing level)."""

    def __init__(self, families, h1, f_match, z_match, top):
        s_fam, d_fam = families['S'], families['D']
        top = np.asarray(top, float)
        self.wd = np.clip((top - s_fam.top) / (d_fam.top - s_fam.top), 0.0, 1.0)
        self.parts = []
        for fam, weight, shift in ((s_fam, 1.0 - self.wd, np.where(top <= s_fam.top, s_fam.top - top, 0.0)),
                                   (d_fam, self.wd, np.where(top >= d_fam.top, d_fam.top - top, 0.0))):
            if not np.any(weight > 0.0):
                continue
            nodes = np.column_stack([np.log(np.maximum(fam.flux(r, z_match + shift), TINY))
                                     for r in range(len(fam.names))])
            w = match_weights(np.log(np.maximum(f_match, TINY)), nodes)
            lines = [fam.line(r, h1 + shift) for r in range(len(fam.names))]
            self.parts.append((fam, weight, shift, w, lines))

    def at(self, z):
        lnk, q, lnf = (np.zeros_like(z) for _ in range(3))
        for fam, weight, shift, w, lines in self.parts:
            for r, (a, b) in enumerate(lines):
                wr = weight * w[:, r]
                if not np.any(wr > 0.0):
                    continue
                zr = fam.clamp(r, z + shift)
                lnk += wr * np.interp(zr, fam.z[r], fam.lnk[r])
                q += wr * np.interp(zr, fam.z[r], fam.q[r])
                lnf += wr * np.log(np.maximum(a + b * zr, TINY))
        return lnk, q, lnf


def closure_columns(families, h1, f_h1, top, edges, subsaturation, sigma=1.0, step=10.0, base_layer=False,
                    exact=False):
    """The flakes' flux down each column's layers below its cloud base, columns as rows: h1, f_h1 (kg/m2/s), top and
    sigma by column; edges (columns by layer edge, top first, m above the freezing level, the first at h1) and
    subsaturation (columns by layer). Steps of at most `step` metres: explicit, as the reference implementation, or with
    exact the equation solved over each step with its coefficients held at the step's middle, u^q falling by
    q s k R^q ds for u = F / F(h1) and R = F^(z) / F^(h1), which holds its accuracy on PlaSim's thick layers. Returns
    the flux at each edge. With base_layer, f_h1 is the flux in the 500 m layer above h1, as McSnow's tests measured it,
    carried down the matched line to h1 first."""
    h1, f_h1, top = (np.asarray(x, float) for x in (h1, f_h1, top))
    sigma = np.broadcast_to(np.asarray(sigma, float), h1.shape)
    z_match = h1 + LAYER / 2.0 if base_layer else h1
    terms = Terms(families, h1, f_h1 / sigma, z_match, top)
    f = f_h1.copy()
    if base_layer:
        f = f * np.exp(terms.at(h1)[2] - terms.at(z_match)[2])
    f1, lnf1 = f.copy(), terms.at(h1)[2]
    edges = np.asarray(edges, float)
    out = np.zeros_like(edges)
    out[:, 0] = f
    for j in range(edges.shape[1] - 1):
        thick = np.maximum(edges[:, j] - edges[:, j + 1], 0.0)
        n = max(1, int(math.ceil(float(thick.max()) / step)))
        ds = thick / n
        s = np.asarray(subsaturation)[:, j]
        for i in range(n):
            zm = edges[:, j] - (i + 0.5) * ds
            lnk, q, lnf = terms.at(zm)
            if exact:
                u = np.where(f1 > 0.0, f / np.where(f1 > 0.0, f1, 1.0), 0.0)
                rate = s * np.exp(lnk + q * (lnf - lnf1)) * ds              # s k R^q ds
                small = np.abs(q) < 1e-6
                qs = np.where(small, 1.0, q)
                uq = np.maximum(u, TINY) ** qs - qs * rate
                u = np.where(small, u * np.exp(-np.minimum(rate, 700.0)),
                             np.where(uq > 0.0, np.maximum(uq, 0.0) ** (1.0 / qs), 0.0))
                f = np.where(f > TINY, u * f1, 0.0)
                continue
            ratio = np.where(f > 0.0, (f / np.where(f1 > 0.0, f1, 1.0)) / np.exp(lnf - lnf1), 1.0)
            f = np.where(f > TINY, f * np.exp(-np.minimum(s * np.exp(lnk) * np.maximum(ratio, TINY) ** (-q) * ds,
                                                          700.0)), 0.0)
        out[:, j + 1] = f
    return out


# ===================================================================================== PlaSim's snow, nsnowsub 1

RHOSU = CM1_RAIN_DENSITY_REFERENCE_PRESSURE / (CM1_DRY_AIR_GAS_CONSTANT * CM1_RAIN_DENSITY_REFERENCE_TEMPERATURE)
EXPMAX = 50.0                     # the most a step's exponent may take, as the rain's (e^50)
UNSET = -1.0e30                   # no forming layer yet
_X5, _W5 = np.polynomial.hermite_e.hermegauss(5)
GAUSS_X = tuple(float(x) for x in _X5)                        # the five-point Gauss-Hermite rule, unit spread
GAUSS_W = tuple(float(w) for w in _W5 / _W5.sum())
# The precipitating fraction (rules_fraction2.txt on the research drive, session-tools/snow_sublimation, 30 October;
# fraction_test2.json, SHA256 fd714507d49bcc64...): at the design boxes' cloud bases the snow's flux falls at local
# rates spread lognormally, flux-weighted, around ln f* = max(ln c + beta ln F, ln F), F the box-mean flux (kg/m2/s),
# with spread s. Measured from each box's snow fields over its 236 3-hourly outputs, each box's numbers gave the other's
# sublimation 1.15 and 1.20 times the column-by-column closure (FT1 holds) and missed one land band at 1.64, where the
# box-mean air alone gives 1.55. The two boxes' pooled numbers, 448 outputs, adopted.
FRACTION_LN_C = -2.8932739162750796
FRACTION_BETA = 0.6273036030191783
FRACTION_SPREAD = 1.0467781310653141
# The fixed-fraction shadow: the mean of the boxes' snowing fractions, the share of their columns carrying at least
# 0.001 mm/h of snow at the cloud base (0.5131 sea, 0.5822 land; fraction_test.json and rules_fraction.txt).
FIXED_FRACTION = 0.5476508706302966
# The box-mean air's own overstatement, the column-by-column closure on the box-mean air over the closure on each
# column's own air (1.194 sea, 1.186 land; fraction_test2_air.json): reported beside the runs and never applied.
AIR_FACTOR = 1.19
# Graupel's shares of PlaSim's snow by kind (register, 30 October): on the design boxes CM1's graupel makes about 0.75
# of the frozen flux in the updraft classes PlaSim's convection stands for and about 0.4 elsewhere (box_graupel.py and
# box_graupel.json on the research drive, SHA256 cd1daaf2c6b72bce...).
GRAUPEL_SHARE_CONVECTIVE = 0.75
GRAUPEL_SHARE_LARGE_SCALE = 0.4
# CM1's graupel mean diameter (1/lambda, flux-weighted) at the boxes' cloud bases by its local rate (box_graupel.json):
# the classes 0-0.1, 0.1-0.3, 0.3-1, 1-3, 3-10 and above 10 mm/h, each at its geometric centre and the end classes
# half a decade beyond their edge, read linearly in ln rate between them and held beyond; the build takes the boxes'
# mean.
GRAUPEL_RATE_MM_H = (0.1 / math.sqrt(10.0), math.sqrt(0.03), math.sqrt(0.3), math.sqrt(3.0), math.sqrt(30.0),
                     10.0 * math.sqrt(10.0))
GRAUPEL_DIAMETER_MM = {
    'sea': (0.671233086424835, 0.7765059844900017, 0.8168877037039266, 0.8912509170459648, 1.2412311473101745,
            1.626313811221787),
    'land': (0.5841864167560554, 0.6794590105134561, 0.7442346026950745, 0.7779165236872287, 0.8864821774308892,
             1.1461838631196377)}
GRAUPEL_DIAMETER_MM['mean'] = tuple(0.5 * (a + b) for a, b in zip(GRAUPEL_DIAMETER_MM['sea'],
                                                                  GRAUPEL_DIAMETER_MM['land']))
GRAUPEL_SIZES = ('mean', 'sea', 'land')
GAMMA_4B = math.gamma(4.0 + CM1_GRAUPEL_FALL_EXPONENT)           # Morrison's CONS7 times 6
GAMMA_VENT = math.gamma(2.5 + CM1_GRAUPEL_FALL_EXPONENT / 2.0)  # his CONS11
# The variants the scheme carries down each column: the adopted one, which PlaSim takes, and with nsnowshad 1 its
# shadows, from the same formation and air and never fed back. Each: its name, its kind (1 the closure on the flakes
# with graupel beside them, 2 graupel alone beside the adopted flakes, 3 PlaSim's gamma), the closure's matching (1 the
# spread of local rates, 2 no fraction, 3 the fixed fraction), its table set (1 the 14 cm references, 2 the 10 cm),
# whether graupel is split off, and the graupel's sizes (1 the boxes' mean, 2 the sea box's, 3 the land box's).
VARIANTS = (('adopted', 1, 1, 1, 1, 1), ('no_fraction', 1, 2, 1, 1, 1), ('fixed_fraction', 1, 3, 1, 1, 1),
            ('cap_10cm', 1, 1, 2, 1, 1), ('no_graupel', 1, 1, 1, 0, 0), ('graupel_sea_sizes', 2, 0, 1, 1, 2),
            ('graupel_land_sizes', 2, 0, 1, 1, 3), ('gamma', 3, 0, 1, 0, 0))
NV, NCV, NNODE, NSLOT = len(VARIANTS), 5, 5, 8
# The diagnostics, means over each output interval: water fluxes in m/s summed over the column, the freezing level's
# height in m (below the ground where the column is colder than 0 C at the ground); for the adopted snow, the flakes'
# flux at their cloud bases, the part of it whose local rates lie below every active family's lightest reference there
# (where the matching holds at that reference, as for snow lighter than 0.1 mm/h), and the sublimation that a layer's
# saturation over ice held back.
DIAGNOSTICS = {
    501: ('snow_formed', 'frozen_precipitation_formed_in_the_column', 'm s-1'),
    502: ('snow_flakes_sublimated', 'flakes_sublimated_in_the_column', 'm s-1'),
    503: ('snow_graupel_sublimated', 'graupel_sublimated_in_the_column', 'm s-1'),
    504: ('snow_flakes_melted', 'flakes_melted_in_the_column', 'm s-1'),
    505: ('snow_graupel_melted', 'graupel_melted_in_the_column', 'm s-1'),
    506: ('snow_nofrac_sublimated', 'frozen_precipitation_sublimated_with_no_fraction_shadow', 'm s-1'),
    507: ('snow_nofrac_melted', 'frozen_precipitation_melted_with_no_fraction_shadow', 'm s-1'),
    508: ('snow_fixfrac_sublimated', 'frozen_precipitation_sublimated_with_a_fixed_fraction_shadow', 'm s-1'),
    509: ('snow_fixfrac_melted', 'frozen_precipitation_melted_with_a_fixed_fraction_shadow', 'm s-1'),
    510: ('snow_cap10_sublimated', 'frozen_precipitation_sublimated_with_10_cm_references_shadow', 'm s-1'),
    511: ('snow_cap10_melted', 'frozen_precipitation_melted_with_10_cm_references_shadow', 'm s-1'),
    512: ('snow_nograupel_sublimated', 'frozen_precipitation_sublimated_all_as_flakes_shadow', 'm s-1'),
    513: ('snow_nograupel_melted', 'frozen_precipitation_melted_all_as_flakes_shadow', 'm s-1'),
    514: ('snow_gsea_sublimated', 'graupel_sublimated_at_the_sea_box_sizes_shadow', 'm s-1'),
    515: ('snow_gsea_melted', 'graupel_melted_at_the_sea_box_sizes_shadow', 'm s-1'),
    516: ('snow_gland_sublimated', 'graupel_sublimated_at_the_land_box_sizes_shadow', 'm s-1'),
    517: ('snow_gland_melted', 'graupel_melted_at_the_land_box_sizes_shadow', 'm s-1'),
    518: ('snow_gamma_sublimated', 'frozen_precipitation_sublimated_with_gamma_shadow', 'm s-1'),
    519: ('snow_gamma_melted', 'frozen_precipitation_melted_with_gamma_shadow', 'm s-1'),
    520: ('freezing_level', 'height_of_the_freezing_level_above_the_ground', 'm'),
    521: ('snow_flakes_base', 'flakes_flux_at_their_cloud_bases', 'm s-1'),
    522: ('snow_flakes_base_light', 'flakes_flux_at_their_cloud_bases_matched_below_the_lightest_reference', 'm s-1'),
    523: ('snow_capped', 'sublimation_the_layers_saturation_over_ice_held_back', 'm s-1')}
NDIAG = len(DIAGNOSTICS)
FREEZING_SLOT, BASE_SLOT, LIGHT_SLOT, CAPPED_SLOT = 19, 20, 21, 22
# (sublimated, melted) slots of each variant in the diagnostics (0-based); the adopted one's are split by kind
SLOTS = {0: ((1, 2), (3, 4)), 1: (5, 6), 2: (7, 8), 3: (9, 10), 4: (11, 12), 5: (13, 14), 6: (15, 16), 7: (17, 18)}
PARAMETERS = dict(nsnowsub=1, nsnowtop=1, nsnowshad=1, snowstep=1000.0, snowgconv=GRAUPEL_SHARE_CONVECTIVE,
                  snowglarge=GRAUPEL_SHARE_LARGE_SCALE, snowlnc=FRACTION_LN_C, snowbeta=FRACTION_BETA,
                  snowspread=FRACTION_SPREAD, snowfixed=FIXED_FRACTION)
NAMELIST_NAMES = tuple(PARAMETERS)
INTEGER_NAMES = ('nsnowsub', 'nsnowtop', 'nsnowshad')


def scheme_tables(path: Path = PRODUCT) -> dict:
    """The product's references as the scheme reads them: the 14 cm references of the shallow and deep families, then
    the 10 cm ones that the sensitivity swaps in; ln k, q and the flux on their 500 m rows from 250 m, each row's
    running sums for the flux lines; the two table sets as slots (four from each family); the families' tops and the
    ladder's lapse."""
    product = json.loads(Path(path).read_text())
    if product['schema'] != SCHEMA:
        raise ValueError(f"expected {SCHEMA}, found {product['schema']}")
    fams = product['families']
    names = list(fams['S']['references']) + list(fams['D']['references'])
    swap = product['sensitivity_10cm']['swap']
    refs = names + [swap[n] for n in names if n in swap]
    family = [0] * 4 + [1] * 4 + [names.index(n) // 4 for n in names if n in swap]
    sets = (tuple(range(8)), tuple(refs.index(swap.get(n, n)) for n in names))
    rows = {n: sorted(product['tables'][n], key=lambda x: x['z_m']) for n in refs}
    nz = [len(rows[n]) for n in refs]
    nzmax = max(nz)
    z0, dz = LAYER / 2.0, LAYER
    lnk, q, f = (np.zeros((len(refs), nzmax)) for _ in range(3))
    ps = np.zeros((len(refs), nzmax + 1, 5))
    for r, n in enumerate(refs):
        for i, x in enumerate(rows[n]):
            if x['z_m'] != z0 + i * dz:
                raise ValueError(f'{n}: rows off the 500 m grid at {x["z_m"]}')
            lnk[r, i], q[r, i], f[r, i] = math.log(x['k_per_m']), x['q'], x['flux_kg_m2_s']
            z = z0 + i * dz
            ps[r, i + 1] = ps[r, i] + np.array([1.0, z, z * z, f[r, i], z * f[r, i]])
    return dict(refs=refs, family=family, sets=sets, nz=nz, z0=z0, dz=dz, lnk=lnk, q=q, f=f, ps=ps,
                tops=(float(fams['S']['top_m']), float(fams['D']['top_m'])),
                lapse=float(product['inputs']['column']['lapse_k_per_m']))


TABLES = scheme_tables()


def _clamp(t, r, zr):
    return min(max(zr, t['z0']), t['z0'] + (t['nz'][r] - 1) * t['dz'])


def _interp(t, a, r, zr):
    """a[r] at a height already held within reference r's rows, linear between them."""
    i = min(int((zr - t['z0']) / t['dz']), t['nz'][r] - 2)
    return a[r, i] + (zr - (t['z0'] + i * t['dz'])) * (a[r, i + 1] - a[r, i]) / t['dz']


def _line(t, r, zr1):
    """Reference r's least-squares line (a, b) of flux on height through its rows from the freezing level up to the one
    holding zr1 + 500 m, two at least (Family.line)."""
    upto = max(_clamp(t, r, zr1), t['dz'] / 2.0) + t['dz']
    m = 0
    for i in range(t['nz'][r]):
        if t['z0'] + i * t['dz'] - t['dz'] / 2.0 <= upto:
            m = i + 1
    m = min(max(m, 2), t['nz'][r])
    n, sz, szz, sf, szf = t['ps'][r, m]
    b = (n * szf - sz * sf) / (n * szz - sz * sz)
    return (sf - b * sz) / n, b


def _match(x, nodes):
    """match_weights for one column: linear in ln F between the sorted references' ln fluxes, held at the ends."""
    order = sorted(range(len(nodes)), key=lambda k: nodes[k])
    xs = [nodes[k] for k in order]
    j = min(max(sum(1 for v in xs if v <= x) - 1, 0), len(xs) - 2)
    lo, hi = xs[j], xs[j + 1]
    f = (x - lo) / (hi - lo) if hi > lo else 0.0
    if x <= xs[0]:
        f = 0.0
    elif x >= xs[-1]:
        f = 1.0
    w = [0.0] * len(nodes)
    w[order[j]] = 1.0 - f
    w[order[j + 1]] += f
    return w


def ice_saturation(zt, p, c):
    """PlaSim's saturation mixing ratio formed over ice, Murray's form with PlaSim's ra1, and its correction factor."""
    qsi = c['rdbrv'] * c['ra1'] * math.exp(MURRAY_ICE_A * (zt - c['TMELT']) / (zt - MURRAY_ICE_B)) / p
    qsi = min(qsi, c['rdbrv'])
    cor = 1.0 / (1.0 - (1.0 / c['rdbrv'] - 1.0) * qsi)
    return max(qsi * cor, 0.0), cor


def graupel_size(fg, size, prm):
    """CM1's graupel mean diameter (m) at graupel flux fg (kg/m2/s), read at its typical local rate."""
    lnf = math.log(max(fg, TINY))
    x = max(prm['snowlnc'] + prm['snowbeta'] * lnf, lnf) + math.log(3600.0)
    xs = [math.log(r) for r in GRAUPEL_RATE_MM_H]
    ys = GRAUPEL_DIAMETER_MM[GRAUPEL_SIZES[size - 1]]
    if x <= xs[0]:
        return ys[0] * 1e-3
    if x >= xs[-1]:
        return ys[-1] * 1e-3
    i = max(k for k in range(len(xs) - 1) if xs[k] <= x)
    return (ys[i] + (x - xs[i]) * (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])) * 1e-3


def graupel_k(dg, zt, zq, qsi, abi, rho, p, gfac):
    """The share of graupel's flux sublimated per metre (Morrison's EPSG (QVI - QV) / ABI over the flux, its number
    falling with its mass so that its size holds), for mean diameter dg in air of temperature zt and density rho."""
    lam = 1.0 / min(max(dg, CM1_GRAUPEL_MEAN_DIAMETER_FLOOR), CM1_GRAUPEL_MEAN_DIAMETER_CAP)
    dens = (RHOSU / rho) ** CM1_RAIN_DENSITY_EXPONENT                 # Morrison's DUM, as for his rain
    agn = gfac * dens * CM1_GRAUPEL_FALL_COEFFICIENT
    b = CM1_GRAUPEL_FALL_EXPONENT
    vm = min(agn * GAMMA_4B / 6.0 / lam ** b, CM1_GRAUPEL_FALL_CAP * gfac * dens)
    mu = CM1_VISCOSITY_COEFFICIENT * zt ** 1.5 / (zt + CM1_VISCOSITY_SUTHERLAND)
    dv = CM1_VAPOUR_DIFFUSIVITY_COEFFICIENT * zt ** CM1_VAPOUR_DIFFUSIVITY_EXPONENT / p
    sc = mu / (rho * dv)
    w = CM1_ICE_VENTILATION_F1 / lam ** 2 + CM1_ICE_VENTILATION_F2 * math.sqrt(agn * rho / mu) * sc ** (1.0 / 3.0) \
        * GAMMA_VENT / lam ** (2.5 + b / 2.0)
    return 2.0 * rho * lam ** 4 * dv * w * max(qsi - zq, 0.0) / (CM1_GRAUPEL_DENSITY * vm * abi)


def snowprep(s: dict, c: dict, prm: dict) -> dict:
    """The column state before mkrain's level loop: each layer's edges by height above the freezing level, from the
    present temperature and humidity (the first crossing of 0 C going up from the ground, interpolated; a column
    colder than 0 C to the ground takes its freezing level below the ground, by the ladder's lapse), and every variant
    empty."""
    n, nlev = s['dq'].shape
    t, q = s['dt'], s['dq']
    zh = np.zeros((n, nlev + 1))
    z0 = np.zeros(n)
    for j in range(n):
        h = np.zeros(nlev + 1)
        zf = np.zeros(nlev)
        for k in range(nlev - 1, -1, -1):
            tv = t[j, k] * (1.0 + (1.0 / c['rdbrv'] - 1.0) * q[j, k])
            dz = c['gascon'] * tv / c['ga'] * s['dsigma'][k] / s['sigma'][k]
            h[k] = h[k + 1] + dz
            zf[k] = h[k + 1] + 0.5 * dz
        pts = [(0.0, t[j, nlev])] + [(zf[k], t[j, k]) for k in range(nlev - 1, -1, -1)]
        zfr = None
        for (zl, tl), (zu, tu) in zip(pts[:-1], pts[1:]):
            if tl > c['TMELT'] >= tu:
                zfr = zl + (tl - c['TMELT']) / (tl - tu) * (zu - zl)
                break
        if zfr is None:
            zfr = -(c['TMELT'] - t[j, nlev]) / TABLES['lapse'] if t[j, nlev] <= c['TMELT'] else h[0]
        z0[j] = zfr
        zh[j] = h - zfr
    shape = (n, NV)
    st = dict(zh=zh, z0=z0, prev_c=np.zeros(n), prev_l=np.zeros(n),
              fl=np.zeros((n, NV, 2)), gr=np.zeros((n, NV, 2)), dg=np.zeros(shape),
              seg=np.zeros(n, bool), h1=np.zeros(n), top=np.zeros(n), wd=np.zeros(n), shift=np.zeros((n, 2)),
              rlnf=np.zeros((n, len(TABLES['refs']))), ra=np.zeros((n, len(TABLES['refs']))),
              rb=np.zeros((n, len(TABLES['refs']))), rl1=np.zeros((n, len(TABLES['refs']))),
              f1=np.zeros(shape), act=np.zeros(shape, bool), nn=np.zeros(shape, int),
              u=np.zeros((n, NV, NNODE)), w=np.zeros((n, NV, NNODE, NSLOT)), lnf1=np.zeros((n, NV, NNODE)),
              topw=np.zeros(n), topz=np.zeros(n), high=np.full(n, UNSET), diag=np.zeros((n, NDIAG)),
              gfac=(c['ga'] / CM1_GRAVITY) ** ((CM1_GRAUPEL_FALL_EXPONENT + 1.0) / 3.0))
    st['diag'][:, FREEZING_SLOT] = z0
    return st


def _variants(prm):
    return range(NV) if prm['nsnowshad'] == 1 else range(1)


def _start(st, j, edge, prm, t=TABLES):
    """A segment from layer edge `edge` (its cloud base): the families' weights and shifts for the snow's top, each
    reference's flux at the cloud base, its flux line and that line's value there; then each closure variant's nodes."""
    s_top, d_top = t['tops']
    h1 = st['zh'][j, edge]
    top = st['topz'][j] / st['topw'][j] if prm['nsnowtop'] == 1 else st['high'][j]
    wd = min(max((top - s_top) / (d_top - s_top), 0.0), 1.0)
    shift = (s_top - top if top <= s_top else 0.0, d_top - top if top >= d_top else 0.0)
    st['seg'][j], st['h1'][j], st['top'][j], st['wd'][j] = True, h1, top, wd
    st['shift'][j] = shift
    for r in range(len(t['refs'])):
        zr1 = _clamp(t, r, h1 + shift[t['family'][r]])
        st['rlnf'][j, r] = math.log(max(_interp(t, t['f'], r, zr1), TINY))
        a, b = _line(t, r, h1 + shift[t['family'][r]])
        st['ra'][j, r], st['rb'][j, r] = a, b
        st['rl1'][j, r] = math.log(max(a + b * zr1, TINY))
    mass = st['mass'][j]
    for v in _variants(prm):
        if VARIANTS[v][1] != 1:
            continue
        f1 = (st['fl'][j, v, 0] + st['fl'][j, v, 1]) * mass
        st['f1'][j, v] = f1
        st['act'][j, v] = f1 > TINY
        if not st['act'][j, v]:
            continue
        lnf1 = math.log(f1)
        kind = VARIANTS[v][2]
        if kind == 1:
            mu = max(prm['snowlnc'] + prm['snowbeta'] * lnf1, lnf1)
            rates = [mu + prm['snowspread'] * x for x in GAUSS_X]
        else:
            rates = [lnf1 if kind == 2 else lnf1 - math.log(prm['snowfixed'])]
        slots = t['sets'][VARIANTS[v][3] - 1]
        st['nn'][j, v] = len(rates)
        if v == 0:
            st['diag'][j, BASE_SLOT] += f1 / 1000.0
        for node, x in enumerate(rates):
            w = [0.0] * NSLOT
            light = True
            for fam, fw in ((0, 1.0 - wd), (1, wd)):
                if fw > 0.0:
                    nodes = [st['rlnf'][j, slots[k]] for k in range(4 * fam, 4 * fam + 4)]
                    light = light and x <= min(nodes)
                    m = _match(x, nodes)
                    for k in range(4):
                        w[4 * fam + k] = fw * m[k]
            if v == 0 and light:
                st['diag'][j, LIGHT_SLOT] += (GAUSS_W[node] if len(rates) == NNODE else 1.0) * f1 / 1000.0
            l1 = 0.0
            for k in range(NSLOT):
                if w[k] > 0.0:
                    l1 += w[k] * st['rl1'][j, slots[k]]
            st['w'][j, v, node], st['lnf1'][j, v, node], st['u'][j, v, node] = w, l1, 1.0


def _integrate(st, j, zt_, zb_, sub, prm, t=TABLES):
    """The closure variants' flakes through one layer, from zt_ down to zb_ (m above the freezing level), in exact steps
    of at most snowstep metres: each step's reference values once, then every node of every variant."""
    vs = [v for v in _variants(prm) if VARIANTS[v][1] == 1 and st['act'][j, v]]
    if not vs:
        return
    nref = 11 if any(VARIANTS[v][3] == 2 for v in vs) else 8
    thick = zt_ - zb_
    if thick <= 0.0:
        return
    n = max(1, math.ceil(thick / prm['snowstep']))
    ds = thick / n
    lk, qq, ll = [0.0] * nref, [0.0] * nref, [0.0] * nref
    for i in range(n):
        zm = zt_ - (i + 0.5) * ds
        for r in range(nref):
            zr = _clamp(t, r, zm + st['shift'][j, t['family'][r]])
            lk[r], qq[r] = _interp(t, t['lnk'], r, zr), _interp(t, t['q'], r, zr)
            ll[r] = math.log(max(st['ra'][j, r] + st['rb'][j, r] * zr, TINY))
        for v in vs:
            slots = t['sets'][VARIANTS[v][3] - 1]
            f1 = st['f1'][j, v]
            for node in range(st['nn'][j, v]):
                u = st['u'][j, v, node]
                if u * f1 <= TINY:
                    st['u'][j, v, node] = 0.0
                    continue
                lnk = qn = lnf = 0.0
                for k in range(NSLOT):
                    wk = st['w'][j, v, node, k]
                    if wk > 0.0:
                        r = slots[k]
                        lnk += wk * lk[r]
                        qn += wk * qq[r]
                        lnf += wk * ll[r]
                rate = sub * math.exp(min(lnk + qn * (lnf - st['lnf1'][j, v, node]), EXPMAX)) * ds
                if abs(qn) < 1e-6:
                    u = u * math.exp(-min(rate, 700.0))
                else:
                    uq = max(u, TINY) ** qn - qn * rate
                    u = uq ** (1.0 / qn) if uq > 0.0 else 0.0
                st['u'][j, v, node] = u


def _kept(st, j, v):
    """The share of a closure variant's flakes at its cloud base still falling (its nodes' weighted mean)."""
    if st['nn'][j, v] == NNODE:
        return sum(w * u for w, u in zip(GAUSS_W, st['u'][j, v]))
    return float(st['u'][j, v, 0])


def _take(a, loss):
    """Take loss (in mkrain's units) from a pair (convective, large-scale) by their shares."""
    tot = a[0] + a[1]
    if tot <= 0.0 or loss <= 0.0:
        return
    fc = a[0] / tot
    a[0] = max(a[0] - loss * fc, 0.0)
    a[1] = max(a[1] - loss * (1.0 - fc), 0.0)


def snowsubl(jlev, zq, zt, zqsat, zcor, zprsc, zprsl, zdqdts, zdtdt, s, c, prm, st):
    """One layer of every column (snowsubl in Fortran, line for line), after mkrain's conversions: the snow formed
    there split into graupel and flakes; a layer forming snow starts the flakes' segment at its lower edge, without
    sublimation; a cold layer below sublimates the flakes by the closure and the graupel by Morrison's law in air
    unsaturated over ice, together at most what would saturate the layer over ice in the step; a warm layer melts all.
    zprsc, zprsl (the adopted variant's snow) and the layer's tendency columns zdqdts, zdtdt change in place."""
    nlev = s['dq'].shape[1]
    st['mass'] = s['dp'] / c['ga']
    vs = _variants(prm)
    gsh = (prm['snowgconv'], prm['snowglarge'])
    for j in range(zprsc.size):
        mass = st['mass'][j]
        if s['dt'][j, jlev] <= c['TMELT']:                # in a cold layer all it forms is snow, PlaSim's rain frozen
            form = (s['dprcl'][j, jlev], s['dprll'][j, jlev])
        else:
            form = (s['dprscl'][j, jlev], s['dprsll'][j, jlev])
        exp_ = (st['prev_c'][j] + form[0], st['prev_l'][j] + form[1])
        warm = s['dt'][j, jlev] > c['TMELT'] or (jlev == nlev - 1 and s['dt'][j, nlev] > c['TMELT'])
        st['diag'][j, 0] += (form[0] + form[1]) * mass / 1000.0
        for v in vs:                                          # the layer's new snow, split by kind
            if VARIANTS[v][4] == 1:
                for kk in range(2):
                    if VARIANTS[v][1] == 1:                    # graupel alone carries no flakes
                        st['fl'][j, v, kk] += (1.0 - gsh[kk]) * form[kk]
                    st['gr'][j, v, kk] += gsh[kk] * form[kk]
            else:
                for kk in range(2):
                    st['fl'][j, v, kk] += form[kk]
        if warm:                                              # all of it melts
            for v in vs:
                sl = SLOTS[v][1]
                fm = (st['fl'][j, v, 0] + st['fl'][j, v, 1]) * mass / 1000.0
                gm = (st['gr'][j, v, 0] + st['gr'][j, v, 1]) * mass / 1000.0
                if v == 0:
                    st['diag'][j, sl[0]] += fm
                    st['diag'][j, sl[1]] += gm
                else:
                    st['diag'][j, sl] += fm + gm
                st['fl'][j, v], st['gr'][j, v], st['dg'][j, v], st['act'][j, v] = 0.0, 0.0, 0.0, False
            st['seg'][j] = False
            st['topw'][j], st['topz'][j], st['high'][j] = 0.0, 0.0, UNSET
            refr = (zprsc[j], zprsl[j])                          # rain frozen again at a cold ground
        else:
            refr = (max(zprsc[j] - exp_[0], 0.0), max(zprsl[j] - exp_[1], 0.0))
        for v in vs:                                          # rain frozen in the layer joins the graupel
            kind, split = VARIANTS[v][1], VARIANTS[v][4]
            for kk in range(2):
                if kind == 3:
                    st['fl'][j, v, kk] += refr[kk]
                elif split == 1:
                    st['gr'][j, v, kk] += refr[kk]
        forming = form[0] + form[1] > 0.0 and not warm
        if forming:
            zmid = 0.5 * (st['zh'][j, jlev] + st['zh'][j, jlev + 1])
            st['topw'][j] += form[0] + form[1]
            st['topz'][j] += (form[0] + form[1]) * zmid
            if st['high'][j] == UNSET:
                st['high'][j] = zmid
            _start(st, j, jlev + 1, prm)
            for v in vs:
                if VARIANTS[v][4] == 1:
                    fg = (st['gr'][j, v, 0] + st['gr'][j, v, 1]) * mass
                    st['dg'][j, v] = graupel_size(fg, VARIANTS[v][5], prm) if fg > TINY else 0.0
        if not warm:
            p = s['dp'][j] * s['sigma'][jlev]
            qsi, cor = ice_saturation(zt[j], p, c)
            zlcps = c['ALS'] / (c['acpd'] * (1.0 + c['ADV'] * s['dq'][j, jlev]))
            abi = 1.0 + zlcps * MURRAY_ICE_A * (c['TMELT'] - MURRAY_ICE_B) * qsi * cor / (zt[j] - MURRAY_ICE_B) ** 2
            sub = max(1.0 - zq[j] / qsi, 0.0) if qsi > 0.0 else 0.0
            cap = max(qsi - zq[j], 0.0) * s['dsigma'][jlev] / c['deltsec2'] / abi * mass
            thick = st['zh'][j, jlev] - st['zh'][j, jlev + 1]
            tv = zt[j] * (1.0 + (1.0 / c['rdbrv'] - 1.0) * zq[j])
            rho = p / (c['gascon'] * tv)
            if not forming and sub > 0.0:
                if not st['seg'][j] and any(st['fl'][j, v, 0] + st['fl'][j, v, 1] > 0.0 for v in vs
                                            if VARIANTS[v][1] == 1):
                    st['topw'][j], st['topz'][j], st['high'][j] = 1.0, st['zh'][j, jlev], st['zh'][j, jlev]
                    _start(st, j, jlev, prm)                     # a segment from this layer's upper edge
                kept0 = [_kept(st, j, v) if VARIANTS[v][1] == 1 and st['act'][j, v] else 1.0 for v in range(NV)]
                u0 = st['u'][j].copy()
                _integrate(st, j, st['zh'][j, jlev], st['zh'][j, jlev + 1], sub, prm)
                lf_main = 0.0
                for v in vs:
                    kind = VARIANTS[v][1]
                    if kind == 3:
                        continue
                    lf = 0.0
                    if kind == 1 and st['act'][j, v]:
                        ff = (st['fl'][j, v, 0] + st['fl'][j, v, 1]) * mass
                        k1 = _kept(st, j, v)
                        lf = ff * (1.0 - k1 / kept0[v]) if kept0[v] > 0.0 else ff
                    fg = (st['gr'][j, v, 0] + st['gr'][j, v, 1]) * mass
                    lg = 0.0
                    if fg > 0.0:
                        if st['dg'][j, v] <= 0.0:
                            st['dg'][j, v] = graupel_size(fg, VARIANTS[v][5], prm)
                        kg = graupel_k(st['dg'][j, v], zt[j], zq[j], qsi, abi, rho, p, st['gfac'])
                        lg = fg * (1.0 - math.exp(-min(kg * thick, EXPMAX)))
                    room = cap if kind == 1 else max(cap - lf_main, 0.0)
                    if lf + lg > room:
                        fac = room / (lf + lg)
                        if v == 0:
                            st['diag'][j, CAPPED_SLOT] += (lf + lg - room) / 1000.0
                        lf, lg = lf * fac, lg * fac
                        if kind == 1 and st['act'][j, v]:
                            st['u'][j, v] = u0[v] - (u0[v] - st['u'][j, v]) * fac
                    _take(st['fl'][j, v], lf / mass)
                    _take(st['gr'][j, v], lg / mass)
                    if v == 0:
                        lf_main = lf
                        st['diag'][j, SLOTS[0][0][0]] += lf / 1000.0
                        st['diag'][j, SLOTS[0][0][1]] += lg / 1000.0
                        dq = (lf + lg) / mass / s['dsigma'][jlev]
                        zdqdts[j] += dq
                        zdtdt[j] -= dq * zlcps
                    else:
                        st['diag'][j, SLOTS[v][0]] += (lf + lg) / 1000.0
            if prm['nsnowshad'] == 1 and zq[j] < zqsat[j]:   # PlaSim's gamma, its own snow, never fed back
                zlcpe = c['ALS'] / (c['acpd'] * (1.0 + c['ADV'] * s['dq'][j, jlev]))
                for kk in (1, 0):
                    if st['fl'][j, 7, kk] > 0.0:
                        d = min(c['gamma'] * (zqsat[j] - zq[j]) * s['dsigma'][jlev] / c['deltsec2']
                                / (1.0 + zlcpe * c['ra2'] * (c['TMELT'] - c['ra4']) * zqsat[j] * zcor[j]
                                   / (zt[j] - c['ra4']) ** 2), st['fl'][j, 7, kk])
                        st['fl'][j, 7, kk] -= d
                        st['diag'][j, SLOTS[7][0]] += d * mass / 1000.0
        zprsc[j] = st['fl'][j, 0, 0] + st['gr'][j, 0, 0]
        zprsl[j] = st['fl'][j, 0, 1] + st['gr'][j, 0, 1]
        st['prev_c'][j], st['prev_l'][j] = zprsc[j], zprsl[j]


# ------------------------------------------------------------------------------------------------- the Fortran

MARKER = 'Terluna: snow from its flakes'


def _flist(values, per_line=4) -> str:
    rows = [', '.join(values[i:i + per_line]) for i in range(0, len(values), per_line)]
    return ' &\n     &  ' + ', &\n     &  '.join(rows) + ' &\n     & '


def declarations(prm: dict = PARAMETERS, t: dict = TABLES) -> str:
    """rainmod's new variables for the snow: the namelist's settings, the closure's tables and the graupel's law as
    parameters, and the column state mkrain's level loop carries."""
    nr, nzmax = len(t['refs']), t['lnk'].shape[1]
    tab = lambda a: _flist([repr(float(a[r, i])) for r in range(nr) for i in range(nzmax)])
    ints = lambda xs: ', '.join(str(int(x)) for x in xs)
    reals = lambda xs: ', '.join(repr(float(x)) for x in xs)
    size_tab = [GRAUPEL_DIAMETER_MM[k][i] * 1e-3 for k in GRAUPEL_SIZES for i in range(len(GRAUPEL_RATE_MM_H))]
    slots = []
    for v in range(NV):
        sub, melt = SLOTS[v]
        slots += [sub[0] + 1, melt[0] + 1, sub[1] + 1, melt[1] + 1] if v == 0 else [sub + 1, melt + 1, 0, 0]
    return (
        f"      integer :: nsnowsub = 0  ! {MARKER} (nevapprec 2): 1 the closure and Morrison's graupel\n"
        "      integer :: nsnowtop = 1  ! the snow's top: 1 formation-weighted height, 2 its highest forming layer\n"
        "      integer :: nsnowshad = 1 ! 1 carries the shadows beside the adopted snow, never fed back\n"
        f"      real :: snowstep = {prm['snowstep']!r} ! the longest exact step through a layer (m)\n"
        f"      real :: snowgconv = {prm['snowgconv']!r}, snowglarge = {prm['snowglarge']!r} ! graupel's shares\n"
        f"      real :: snowlnc = {prm['snowlnc']!r}, snowbeta = {prm['snowbeta']!r} ! the fraction's law\n"
        f"      real :: snowspread = {prm['snowspread']!r}, snowfixed = {prm['snowfixed']!r}\n"
        f"      integer, parameter :: NSNR = {nr}, NSNZ = {nzmax}, NSNV = {NV}, NSNC = {NCV}, NSNN = {NNODE}\n"
        f"      integer, parameter :: NSNS = {NSLOT}, NSNG = {len(GRAUPEL_RATE_MM_H)}\n"
        f"      integer, parameter :: NSNFRZ = {FREEZING_SLOT + 1}, NSNBASE = {BASE_SLOT + 1}, NSNLIGHT = "
        f"{LIGHT_SLOT + 1}, NSNCAP = {CAPPED_SLOT + 1}\n"
        f"      real, parameter :: snz0 = {t['z0']!r}, sndz = {t['dz']!r} ! the references' rows (m)\n"
        f"      real, parameter :: sntops = {t['tops'][0]!r}, sntopd = {t['tops'][1]!r} ! the families' tops (m)\n"
        f"      real, parameter :: snlapse = {t['lapse']!r} ! the ladder's lapse (K/m)\n"
        f"      integer, parameter :: snnz(NSNR) = (/{ints(t['nz'])}/)\n"
        f"      integer, parameter :: snfam(NSNR) = (/{ints([f + 1 for f in t['family']])}/)\n"
        f"      integer, parameter :: snset(NSNS,2) = reshape((/{ints([r + 1 for s in t['sets'] for r in s])}/), &\n"
        "     &  (/NSNS,2/))\n"
        "      real, parameter :: snlk(NSNZ,NSNR) = reshape((/" + tab(t['lnk']) + "/), (/NSNZ,NSNR/))\n"
        "      real, parameter :: snq(NSNZ,NSNR) = reshape((/" + tab(t['q']) + "/), (/NSNZ,NSNR/))\n"
        "      real, parameter :: snf(NSNZ,NSNR) = reshape((/" + tab(t['f']) + "/), (/NSNZ,NSNR/))\n"
        "      real, parameter :: sngx(NSNN) = (/" + _flist([repr(x) for x in GAUSS_X], 3) + "/)\n"
        "      real, parameter :: sngw(NSNN) = (/" + _flist([repr(x) for x in GAUSS_W], 3) + "/)\n"
        f"      integer, parameter :: snvkind(NSNV) = (/{ints(v[1] for v in VARIANTS)}/)\n"
        f"      integer, parameter :: snvmat(NSNV) = (/{ints(v[2] for v in VARIANTS)}/)\n"
        f"      integer, parameter :: snvset(NSNV) = (/{ints(v[3] for v in VARIANTS)}/)\n"
        f"      integer, parameter :: snvsplit(NSNV) = (/{ints(v[4] for v in VARIANTS)}/)\n"
        f"      integer, parameter :: snvsize(NSNV) = (/{ints(v[5] for v in VARIANTS)}/)\n"
        "      integer, parameter :: snslot(4,NSNV) = reshape((/" + _flist([str(x) for x in slots], 16)
        + "/), (/4,NSNV/))\n"
        "      real, parameter :: sngrate(NSNG) = (/" + _flist([repr(math.log(r)) for r in GRAUPEL_RATE_MM_H], 3)
        + "/) ! ln mm/h\n"
        "      real, parameter :: sngdia(NSNG,3) = reshape((/" + _flist([repr(x) for x in size_tab]) + "/), &\n"
        "     &  (/NSNG,3/))\n"
        f"      real, parameter :: sngra = {CM1_GRAUPEL_FALL_COEFFICIENT!r}, sngrb = {CM1_GRAUPEL_FALL_EXPONENT!r}"
        f", sngrho = {CM1_GRAUPEL_DENSITY!r}, sngcap = {CM1_GRAUPEL_FALL_CAP!r} ! Morrison's graupel\n"
        f"      real, parameter :: sngdmin = {CM1_GRAUPEL_MEAN_DIAMETER_FLOOR!r}, sngdmax = "
        f"{CM1_GRAUPEL_MEAN_DIAMETER_CAP!r}, sngf1 = {CM1_ICE_VENTILATION_F1!r}, sngf2 = {CM1_ICE_VENTILATION_F2!r}\n"
        f"      real, parameter :: sng4b = {GAMMA_4B!r}, sngvt = {GAMMA_VENT!r}, sngref = {CM1_GRAVITY!r}\n"
        f"      real, parameter :: snrhosu = {RHOSU!r}, snexpd = {CM1_RAIN_DENSITY_EXPONENT!r}\n"
        f"      real, parameter :: sndv0 = {CM1_VAPOUR_DIFFUSIVITY_COEFFICIENT!r}, sndvx = "
        f"{CM1_VAPOUR_DIFFUSIVITY_EXPONENT!r}, snmu0 = {CM1_VISCOSITY_COEFFICIENT!r}, snmus = "
        f"{CM1_VISCOSITY_SUTHERLAND!r}\n"
        f"      real, parameter :: snmura = {MURRAY_ICE_A!r}, snmurb = {MURRAY_ICE_B!r} ! Murray's form over ice\n"
        f"      real, parameter :: sntiny = {TINY!r}, snxmax = {EXPMAX!r}, snunset = {UNSET!r}\n"
        "      real :: snps(5,0:NSNZ,NSNR), sngfac\n"
        "      integer :: nsninit = 0\n"
        "      real :: snzh(NHOR,NLEP), snz0c(NHOR), snprevc(NHOR), snprevl(NHOR)\n"
        "      real :: snfl(2,NSNV,NHOR), sngr(2,NSNV,NHOR), sndg(NSNV,NHOR), snf1(NSNV,NHOR)\n"
        "      integer :: snseg(NHOR), snact(NSNV,NHOR), snnn(NSNV,NHOR)\n"
        "      real :: snh1(NHOR), sntop(NHOR), snwd(NHOR), snshift(2,NHOR)\n"
        "      real :: snrlnf(NSNR,NHOR), snra(NSNR,NHOR), snrb(NSNR,NHOR), snrl1(NSNR,NHOR)\n"
        "      real :: snu(NSNN,NSNV,NHOR), snw(NSNS,NSNN,NSNV,NHOR), snlnf1(NSNN,NSNV,NHOR)\n"
        "      real :: sntopw(NHOR), sntopz(NHOR), snhigh(NHOR)\n")


ROUTINES = r"""
!     ===================
!     SUBROUTINE SNOWINIT
!     ===================

      subroutine snowinit
      use rainmod
!
!     Terluna: snow from its flakes: each reference's running sums of its rows for the flux lines, and the gravity
!     factor of Morrison's graupel speeds, (g/9.81)**((b+1)/3)
!
      integer :: jr,jz
      real :: zz
      do jr=1,NSNR
       snps(:,0,jr)=0.
       do jz=1,snnz(jr)
        zz=snz0+REAL(jz-1)*sndz
        snps(1,jz,jr)=snps(1,jz-1,jr)+1.
        snps(2,jz,jr)=snps(2,jz-1,jr)+zz
        snps(3,jz,jr)=snps(3,jz-1,jr)+zz*zz
        snps(4,jz,jr)=snps(4,jz-1,jr)+snf(jz,jr)
        snps(5,jz,jr)=snps(5,jz-1,jr)+zz*snf(jz,jr)
       enddo
      enddo
      sngfac=(ga/sngref)**((sngrb+1.)/3.)
      nsninit=1
      return
      end subroutine snowinit

!     ===================
!     FUNCTION SNOWCLAMP
!     ===================

      real function snowclamp(kr,pz)
      use rainmod
!     Terluna: a height held within reference kr's rows
      integer :: kr
      real :: pz
      snowclamp=MIN(MAX(pz,snz0),snz0+REAL(snnz(kr)-1)*sndz)
      return
      end function snowclamp

!     =================
!     FUNCTION SNOWTAB
!     =================

      real function snowtab(ktab,kr,pz)
      use rainmod
!     Terluna: reference kr's ln k (ktab 1), q (2) or flux (3) at a height held within its rows, linear between them
      integer :: ktab,kr,ji
      real :: pz,za,zb
      ji=MIN(INT((pz-snz0)/sndz),snnz(kr)-2)+1
      if(ktab == 1) then
       za=snlk(ji,kr)
       zb=snlk(ji+1,kr)
      elseif(ktab == 2) then
       za=snq(ji,kr)
       zb=snq(ji+1,kr)
      else
       za=snf(ji,kr)
       zb=snf(ji+1,kr)
      endif
      snowtab=za+(pz-(snz0+REAL(ji-1)*sndz))*(zb-za)/sndz
      return
      end function snowtab

!     ===================
!     SUBROUTINE SNOWLINE
!     ===================

      subroutine snowline(kr,pzr1,pa,pb)
      use rainmod
!     Terluna: reference kr's least-squares line of flux on height through its rows from the freezing level up to the
!     one holding pzr1 + 500 m, two at least
      integer :: kr,jz,im
      real :: pzr1,pa,pb,zupto,zn,zsz,zszz,zsf,zszf
      real :: snowclamp
      zupto=MAX(snowclamp(kr,pzr1),sndz/2.)+sndz
      im=0
      do jz=1,snnz(kr)
       if(snz0+REAL(jz-1)*sndz-sndz/2. <= zupto) im=jz
      enddo
      im=MIN(MAX(im,2),snnz(kr))
      zn=snps(1,im,kr)
      zsz=snps(2,im,kr)
      zszz=snps(3,im,kr)
      zsf=snps(4,im,kr)
      zszf=snps(5,im,kr)
      pb=(zn*zszf-zsz*zsf)/(zn*zszz-zsz*zsz)
      pa=(zsf-pb*zsz)/zn
      return
      end subroutine snowline

!     ====================
!     SUBROUTINE SNOWMATCH
!     ====================

      subroutine snowmatch(px,pnodes,pw)
      use rainmod
!     Terluna: weights over four references from ln F against their ln fluxes, linear between neighbours in sorted
!     order and held at the ends
      real :: px,pnodes(4),pw(4),zxs(4),zlo,zhi,zf
      integer :: iord(4),ja,jb,jt,jj
      do ja=1,4
       iord(ja)=ja
      enddo
      do ja=2,4
       jb=ja
       do while(jb > 1)
        if(pnodes(iord(jb)) < pnodes(iord(jb-1))) then
         jt=iord(jb)
         iord(jb)=iord(jb-1)
         iord(jb-1)=jt
         jb=jb-1
        else
         exit
        endif
       enddo
      enddo
      do ja=1,4
       zxs(ja)=pnodes(iord(ja))
      enddo
      jj=0
      do ja=1,4
       if(zxs(ja) <= px) jj=jj+1
      enddo
      jj=MIN(MAX(jj-1,0),2)+1
      zlo=zxs(jj)
      zhi=zxs(jj+1)
      zf=0.
      if(zhi > zlo) zf=(px-zlo)/(zhi-zlo)
      if(px <= zxs(1)) then
       zf=0.
      elseif(px >= zxs(4)) then
       zf=1.
      endif
      pw(:)=0.
      pw(iord(jj))=1.-zf
      pw(iord(jj+1))=pw(iord(jj+1))+zf
      return
      end subroutine snowmatch

!     ====================
!     SUBROUTINE SNOWPREP
!     ====================

      subroutine snowprep
      use rainmod
!
!     Terluna: snow from its flakes (climate/gcm/snow_sublimation.py, snowprep): each layer's edges by height above
!     the column's freezing level, from the present temperature and humidity: the first crossing of 0 C going up from
!     the ground, interpolated, and in a column colder than 0 C to the ground a level below the ground by the
!     references' ladder lapse; and every variant empty
!
      integer :: jhor,jlev,jp
      real :: zh(NLEP),zf(NLEV),ztv,zdz,zfr,zzl,ztl,zzu,ztu,zfound
      if(nsninit == 0) call snowinit
      do jhor=1,NHOR
       zh(NLEP)=0.
       do jlev=NLEV,1,-1
        ztv=dt(jhor,jlev)*(1.+(1./rdbrv-1.)*dq(jhor,jlev))
        zdz=gascon*ztv/ga*dsigma(jlev)/sigma(jlev)
        zh(jlev)=zh(jlev+1)+zdz
        zf(jlev)=zh(jlev+1)+0.5*zdz
       enddo
       zfound=0.
       zzl=0.
       ztl=dt(jhor,NLEP)
       do jp=NLEV,1,-1
        zzu=zf(jp)
        ztu=dt(jhor,jp)
        if(ztl > TMELT .and. TMELT >= ztu) then
         zfr=zzl+(ztl-TMELT)/(ztl-ztu)*(zzu-zzl)
         zfound=1.
         exit
        endif
        zzl=zzu
        ztl=ztu
       enddo
       if(zfound == 0.) then
        if(dt(jhor,NLEP) <= TMELT) then
         zfr=-(TMELT-dt(jhor,NLEP))/snlapse
        else
         zfr=zh(1)
        endif
       endif
       snz0c(jhor)=zfr
       snzh(jhor,:)=zh(:)-zfr
       asnow(jhor,NSNFRZ)=asnow(jhor,NSNFRZ)+zfr
      enddo
      snprevc(:)=0.
      snprevl(:)=0.
      snfl(:,:,:)=0.
      sngr(:,:,:)=0.
      sndg(:,:)=0.
      snseg(:)=0
      snact(:,:)=0
      snnn(:,:)=0
      snf1(:,:)=0.
      snu(:,:,:)=0.
      snw(:,:,:,:)=0.
      snlnf1(:,:,:)=0.
      sntopw(:)=0.
      sntopz(:)=0.
      snhigh(:)=snunset
      nsnacc=nsnacc+1
      return
      end subroutine snowprep

!     =====================
!     SUBROUTINE SNOWSTART
!     =====================

      subroutine snowstart(jhor,kedge,pmass,knv)
      use rainmod
!
!     Terluna: a segment from layer edge kedge, the snow's cloud base: the families' weights and shifts for the snow's
!     top, each reference's flux there, its flux line and that line's value there; then each closure variant's nodes
!     (snow_sublimation._start)
!
      integer :: jhor,kedge,knv,jr,jv,jn,jk,jfam,inn,islot,ilight
      real :: pmass,zh1,ztop,zwd,zs1,zs2,zr1,za,zb,zfl,zlnf1,zmu,zx,zfw,zl1
      real :: zrates(NSNN),znodes(4),zm(4),zwk(NSNS)
      real :: snowclamp,snowtab
      zh1=snzh(jhor,kedge)
      if(nsnowtop == 1) then
       ztop=sntopz(jhor)/sntopw(jhor)
      else
       ztop=snhigh(jhor)
      endif
      zwd=MIN(MAX((ztop-sntops)/(sntopd-sntops),0.),1.)
      zs1=0.
      if(ztop <= sntops) zs1=sntops-ztop
      zs2=0.
      if(ztop >= sntopd) zs2=sntopd-ztop
      snseg(jhor)=1
      snh1(jhor)=zh1
      sntop(jhor)=ztop
      snwd(jhor)=zwd
      snshift(1,jhor)=zs1
      snshift(2,jhor)=zs2
      do jr=1,NSNR
       zr1=snowclamp(jr,zh1+snshift(snfam(jr),jhor))
       snrlnf(jr,jhor)=LOG(MAX(snowtab(3,jr,zr1),sntiny))
       call snowline(jr,zh1+snshift(snfam(jr),jhor),za,zb)
       snra(jr,jhor)=za
       snrb(jr,jhor)=zb
       snrl1(jr,jhor)=LOG(MAX(za+zb*zr1,sntiny))
      enddo
      do jv=1,knv
       if(snvkind(jv) /= 1) cycle
       zfl=(snfl(1,jv,jhor)+snfl(2,jv,jhor))*pmass
       snf1(jv,jhor)=zfl
       snact(jv,jhor)=0
       if(zfl > sntiny) snact(jv,jhor)=1
       if(snact(jv,jhor) == 0) cycle
       zlnf1=LOG(zfl)
       if(snvmat(jv) == 1) then
        zmu=MAX(snowlnc+snowbeta*zlnf1,zlnf1)
        inn=NSNN
        do jn=1,NSNN
         zrates(jn)=zmu+snowspread*sngx(jn)
        enddo
       elseif(snvmat(jv) == 2) then
        inn=1
        zrates(1)=zlnf1
       else
        inn=1
        zrates(1)=zlnf1-LOG(snowfixed)
       endif
       snnn(jv,jhor)=inn
       if(jv == 1) asnow(jhor,NSNBASE)=asnow(jhor,NSNBASE)+zfl/1000.
       do jn=1,inn
        zx=zrates(jn)
        zwk(:)=0.
        ilight=1
        do jfam=1,2
         if(jfam == 1) then
          zfw=1.-zwd
         else
          zfw=zwd
         endif
         if(zfw > 0.) then
          do jk=1,4
           znodes(jk)=snrlnf(snset(4*(jfam-1)+jk,snvset(jv)),jhor)
          enddo
          if(zx > MINVAL(znodes)) ilight=0
          call snowmatch(zx,znodes,zm)
          do jk=1,4
           zwk(4*(jfam-1)+jk)=zfw*zm(jk)
          enddo
         endif
        enddo
        if(jv == 1 .and. ilight == 1) then
         if(inn == NSNN) then
          asnow(jhor,NSNLIGHT)=asnow(jhor,NSNLIGHT)+sngw(jn)*zfl/1000.
         else
          asnow(jhor,NSNLIGHT)=asnow(jhor,NSNLIGHT)+zfl/1000.
         endif
        endif
        zl1=0.
        do jk=1,NSNS
         if(zwk(jk) > 0.) then
          islot=snset(jk,snvset(jv))
          zl1=zl1+zwk(jk)*snrl1(islot,jhor)
         endif
        enddo
        snw(:,jn,jv,jhor)=zwk(:)
        snlnf1(jn,jv,jhor)=zl1
        snu(jn,jv,jhor)=1.
       enddo
      enddo
      return
      end subroutine snowstart

!     ===================
!     SUBROUTINE SNOWINT
!     ===================

      subroutine snowint(jhor,pztop,pzbot,psub,knv)
      use rainmod
!
!     Terluna: the closure variants' flakes through one layer, from pztop down to pzbot (m above the freezing level),
!     in exact steps of at most snowstep metres: each step's reference values once, then every node of every variant
!     (snow_sublimation._integrate): u**q falls by q s k R**q ds
!
      integer :: jhor,knv,jv,jn,jk,jr,jstep,isteps,nref,iany,iuse(NSNV),islot
      real :: pztop,pzbot,psub,zthick,zds,zzm,zzr,zu,zlnk,zqn,zlnf,zrate,zuq,zwk
      real :: zlk(NSNR),zqq(NSNR),zll(NSNR)
      real :: snowclamp,snowtab
      iany=0
      nref=8
      do jv=1,knv
       iuse(jv)=0
       if(snvkind(jv) == 1 .and. snact(jv,jhor) == 1) then
        iuse(jv)=1
        iany=1
        if(snvset(jv) == 2) nref=NSNR
       endif
      enddo
      if(iany == 0) return
      zthick=pztop-pzbot
      if(zthick <= 0.) return
      isteps=MAX(1,CEILING(zthick/snowstep))
      zds=zthick/REAL(isteps)
      do jstep=1,isteps
       zzm=pztop-(REAL(jstep-1)+0.5)*zds
       do jr=1,nref
        zzr=snowclamp(jr,zzm+snshift(snfam(jr),jhor))
        zlk(jr)=snowtab(1,jr,zzr)
        zqq(jr)=snowtab(2,jr,zzr)
        zll(jr)=LOG(MAX(snra(jr,jhor)+snrb(jr,jhor)*zzr,sntiny))
       enddo
       do jv=1,knv
        if(iuse(jv) == 0) cycle
        do jn=1,snnn(jv,jhor)
         zu=snu(jn,jv,jhor)
         if(zu*snf1(jv,jhor) <= sntiny) then
          snu(jn,jv,jhor)=0.
          cycle
         endif
         zlnk=0.
         zqn=0.
         zlnf=0.
         do jk=1,NSNS
          zwk=snw(jk,jn,jv,jhor)
          if(zwk > 0.) then
           islot=snset(jk,snvset(jv))
           zlnk=zlnk+zwk*zlk(islot)
           zqn=zqn+zwk*zqq(islot)
           zlnf=zlnf+zwk*zll(islot)
          endif
         enddo
         zrate=psub*EXP(MIN(zlnk+zqn*(zlnf-snlnf1(jn,jv,jhor)),snxmax))*zds
         if(ABS(zqn) < 1.E-6) then
          zu=zu*EXP(-MIN(zrate,700.))
         else
          zuq=MAX(zu,sntiny)**zqn-zqn*zrate
          if(zuq > 0.) then
           zu=zuq**(1./zqn)
          else
           zu=0.
          endif
         endif
         snu(jn,jv,jhor)=zu
        enddo
       enddo
      enddo
      return
      end subroutine snowint

!     ===================
!     FUNCTION SNOWKEPT
!     ===================

      real function snowkept(jhor,kv)
      use rainmod
!     Terluna: the share of a closure variant's flakes at its cloud base still falling, its nodes' weighted mean
      integer :: jhor,kv,jn
      if(snnn(kv,jhor) == NSNN) then
       snowkept=0.
       do jn=1,NSNN
        snowkept=snowkept+sngw(jn)*snu(jn,kv,jhor)
       enddo
      else
       snowkept=snu(1,kv,jhor)
      endif
      return
      end function snowkept

!     ====================
!     FUNCTION SNOWGSIZE
!     ====================

      real function snowgsize(pfg,ksize)
      use rainmod
!     Terluna: CM1's graupel mean diameter (m) at graupel flux pfg (kg/m2/s), read at its typical local rate
      integer :: ksize,jg,ji
      real :: pfg,zlnf,zx
      zlnf=LOG(MAX(pfg,sntiny))
      zx=MAX(snowlnc+snowbeta*zlnf,zlnf)+LOG(3600.)
      if(zx <= sngrate(1)) then
       snowgsize=sngdia(1,ksize)
      elseif(zx >= sngrate(NSNG)) then
       snowgsize=sngdia(NSNG,ksize)
      else
       ji=1
       do jg=1,NSNG-1
        if(sngrate(jg) <= zx) ji=jg
       enddo
       snowgsize=sngdia(ji,ksize)+(zx-sngrate(ji))*(sngdia(ji+1,ksize)-sngdia(ji,ksize))                         &
     &          /(sngrate(ji+1)-sngrate(ji))
      endif
      return
      end function snowgsize

!     =================
!     FUNCTION SNOWGK
!     =================

      real function snowgk(pdg,pt,pq,pqsi,pabi,prho,pp)
      use rainmod
!     Terluna: the share of graupel's flux sublimated per metre, Morrison's EPSG (QVI - QV) / ABI over its flux, its
!     number falling with its mass so that its size holds
      real :: pdg,pt,pq,pqsi,pabi,prho,pp,zlam,zdens,zagn,zvm,zmu,zdv,zsc,zw
      zlam=1./MIN(MAX(pdg,sngdmin),sngdmax)
      zdens=(snrhosu/prho)**snexpd
      zagn=sngfac*zdens*sngra
      zvm=MIN(zagn*sng4b/6./zlam**sngrb,sngcap*sngfac*zdens)
      zmu=snmu0*pt**1.5/(pt+snmus)
      zdv=sndv0*pt**sndvx/pp
      zsc=zmu/(prho*zdv)
      zw=sngf1/zlam**2+sngf2*SQRT(zagn*prho/zmu)*zsc**(1./3.)*sngvt/zlam**(2.5+sngrb/2.)
      snowgk=2.*prho*zlam**4*zdv*zw*MAX(pqsi-pq,0.)/(sngrho*zvm*pabi)
      return
      end function snowgk

!     =================
!     SUBROUTINE SNOWTAKE
!     =================

      subroutine snowtake(pa,ploss)
      use rainmod
!     Terluna: take ploss (mkrain's units) from a convective and large-scale pair by their shares
      real :: pa(2),ploss,ztot,zfc
      ztot=pa(1)+pa(2)
      if(ztot <= 0. .or. ploss <= 0.) return
      zfc=pa(1)/ztot
      pa(1)=MAX(pa(1)-ploss*zfc,0.)
      pa(2)=MAX(pa(2)-ploss*(1.-zfc),0.)
      return
      end subroutine snowtake

!     ====================
!     SUBROUTINE SNOWSUBL
!     ====================

      subroutine snowsubl(jlev,pq,pt,pqsat,pcor,pprsc,pprsl,pdqdt,pdtdt)
      use rainmod
!
!     Terluna: snow from its flakes (nsnowsub 1; climate/gcm/snow_sublimation.py, snowsubl, line for line). One layer
!     of every column after mkrain's conversions: the layer's new snow split into graupel and flakes; a layer forming
!     snow starts the flakes' segment at its lower edge, without sublimation; a cold layer below sublimates the flakes
!     by the closure and the graupel by Morrison's law in air unsaturated over ice, together at most what would
!     saturate the layer over ice in the step; a warm layer melts all. The shadows follow from the same formation and
!     air and are never fed back; pprsc, pprsl and the layer's tendencies change for the adopted snow alone.
!
      integer :: jlev,jhor,jv,jk,knv,iform,iwarm,ikind,isplit
      real :: pq(NHOR),pt(NHOR),pqsat(NHOR),pcor(NHOR),pprsc(NHOR),pprsl(NHOR),pdqdt(NHOR),pdtdt(NHOR)
      real :: zmass,zform(2),zexp(2),zrefr(2),zgsh(2),zfm,zgm,zmid,zfg,zp,zqsi,zcor,zlcps,zabi,zsub,zcap
      real :: zthick,ztv,zrho,zlfm,zff,zk1,zlf,zlg,zkg,zroom,zfac,zdq,zd,zlcpe,zfl
      real :: zkept0(NSNV),zu0(NSNN,NSNV)
      real :: snowkept,snowgsize,snowgk
      knv=1
      if(nsnowshad == 1) knv=NSNV
      zgsh(1)=snowgconv
      zgsh(2)=snowglarge
      do jhor=1,NHOR
       zmass=dp(jhor)/ga
       if(dt(jhor,jlev) <= TMELT) then   ! in a cold layer all it forms is snow, PlaSim's rain frozen
        zform(1)=dprcl(jhor,jlev)
        zform(2)=dprll(jhor,jlev)
       else
        zform(1)=dprscl(jhor,jlev)
        zform(2)=dprsll(jhor,jlev)
       endif
       zexp(1)=snprevc(jhor)+zform(1)
       zexp(2)=snprevl(jhor)+zform(2)
       iwarm=0
       if(dt(jhor,jlev) > TMELT) iwarm=1
       if(jlev == NLEV .and. dt(jhor,NLEP) > TMELT) iwarm=1
       asnow(jhor,1)=asnow(jhor,1)+(zform(1)+zform(2))*zmass/1000.
       do jv=1,knv
        do jk=1,2
         if(snvsplit(jv) == 1) then
          if(snvkind(jv) == 1) snfl(jk,jv,jhor)=snfl(jk,jv,jhor)+(1.-zgsh(jk))*zform(jk)
          sngr(jk,jv,jhor)=sngr(jk,jv,jhor)+zgsh(jk)*zform(jk)
         else
          snfl(jk,jv,jhor)=snfl(jk,jv,jhor)+zform(jk)
         endif
        enddo
       enddo
       if(iwarm == 1) then
        do jv=1,knv
         zfm=(snfl(1,jv,jhor)+snfl(2,jv,jhor))*zmass/1000.
         zgm=(sngr(1,jv,jhor)+sngr(2,jv,jhor))*zmass/1000.
         if(jv == 1) then
          asnow(jhor,snslot(2,jv))=asnow(jhor,snslot(2,jv))+zfm
          asnow(jhor,snslot(4,jv))=asnow(jhor,snslot(4,jv))+zgm
         else
          asnow(jhor,snslot(2,jv))=asnow(jhor,snslot(2,jv))+zfm+zgm
         endif
         snfl(:,jv,jhor)=0.
         sngr(:,jv,jhor)=0.
         sndg(jv,jhor)=0.
         snact(jv,jhor)=0
        enddo
        snseg(jhor)=0
        sntopw(jhor)=0.
        sntopz(jhor)=0.
        snhigh(jhor)=snunset
        zrefr(1)=pprsc(jhor)
        zrefr(2)=pprsl(jhor)
       else
        zrefr(1)=MAX(pprsc(jhor)-zexp(1),0.)
        zrefr(2)=MAX(pprsl(jhor)-zexp(2),0.)
       endif
       do jv=1,knv
        do jk=1,2
         if(snvkind(jv) == 3) then
          snfl(jk,jv,jhor)=snfl(jk,jv,jhor)+zrefr(jk)
         elseif(snvsplit(jv) == 1) then
          sngr(jk,jv,jhor)=sngr(jk,jv,jhor)+zrefr(jk)
         endif
        enddo
       enddo
       iform=0
       if(zform(1)+zform(2) > 0. .and. iwarm == 0) iform=1
       if(iform == 1) then
        zmid=0.5*(snzh(jhor,jlev)+snzh(jhor,jlev+1))
        sntopw(jhor)=sntopw(jhor)+(zform(1)+zform(2))
        sntopz(jhor)=sntopz(jhor)+(zform(1)+zform(2))*zmid
        if(snhigh(jhor) == snunset) snhigh(jhor)=zmid
        call snowstart(jhor,jlev+1,zmass,knv)
        do jv=1,knv
         if(snvsplit(jv) == 1) then
          zfg=(sngr(1,jv,jhor)+sngr(2,jv,jhor))*zmass
          sndg(jv,jhor)=0.
          if(zfg > sntiny) sndg(jv,jhor)=snowgsize(zfg,snvsize(jv))
         endif
        enddo
       endif
       if(iwarm == 0) then
        zp=dp(jhor)*sigma(jlev)
        zqsi=rdbrv*ra1*EXP(snmura*(pt(jhor)-TMELT)/(pt(jhor)-snmurb))/zp
        zqsi=MIN(zqsi,rdbrv)
        zcor=1./(1.-(1./rdbrv-1.)*zqsi)
        zqsi=MAX(zqsi*zcor,0.)
        zlcps=ALS/(acpd*(1.+ADV*dq(jhor,jlev)))
        zabi=1.+zlcps*snmura*(TMELT-snmurb)*zqsi*zcor/(pt(jhor)-snmurb)**2
        zsub=0.
        if(zqsi > 0.) zsub=MAX(1.-pq(jhor)/zqsi,0.)
        zcap=MAX(zqsi-pq(jhor),0.)*dsigma(jlev)/deltsec2/zabi*zmass
        zthick=snzh(jhor,jlev)-snzh(jhor,jlev+1)
        ztv=pt(jhor)*(1.+(1./rdbrv-1.)*pq(jhor))
        zrho=zp/(gascon*ztv)
        if(iform == 0 .and. zsub > 0.) then
         if(snseg(jhor) == 0) then
          zfl=0.
          do jv=1,knv
           if(snvkind(jv) == 1) zfl=MAX(zfl,snfl(1,jv,jhor)+snfl(2,jv,jhor))
          enddo
          if(zfl > 0.) then
           sntopw(jhor)=1.
           sntopz(jhor)=snzh(jhor,jlev)
           snhigh(jhor)=snzh(jhor,jlev)
           call snowstart(jhor,jlev,zmass,knv)
          endif
         endif
         do jv=1,knv
          zkept0(jv)=1.
          if(snvkind(jv) == 1 .and. snact(jv,jhor) == 1) zkept0(jv)=snowkept(jhor,jv)
          zu0(:,jv)=snu(:,jv,jhor)
         enddo
         call snowint(jhor,snzh(jhor,jlev),snzh(jhor,jlev+1),zsub,knv)
         zlfm=0.
         do jv=1,knv
          ikind=snvkind(jv)
          if(ikind == 3) cycle
          zlf=0.
          if(ikind == 1 .and. snact(jv,jhor) == 1) then
           zff=(snfl(1,jv,jhor)+snfl(2,jv,jhor))*zmass
           zk1=snowkept(jhor,jv)
           if(zkept0(jv) > 0.) then
            zlf=zff*(1.-zk1/zkept0(jv))
           else
            zlf=zff
           endif
          endif
          zfg=(sngr(1,jv,jhor)+sngr(2,jv,jhor))*zmass
          zlg=0.
          if(zfg > 0.) then
           if(sndg(jv,jhor) <= 0.) sndg(jv,jhor)=snowgsize(zfg,snvsize(jv))
           zkg=snowgk(sndg(jv,jhor),pt(jhor),pq(jhor),zqsi,zabi,zrho,zp)
           zlg=zfg*(1.-EXP(-MIN(zkg*zthick,snxmax)))
          endif
          if(ikind == 1) then
           zroom=zcap
          else
           zroom=MAX(zcap-zlfm,0.)
          endif
          if(zlf+zlg > zroom) then
           zfac=zroom/(zlf+zlg)
           if(jv == 1) asnow(jhor,NSNCAP)=asnow(jhor,NSNCAP)+(zlf+zlg-zroom)/1000.
           zlf=zlf*zfac
           zlg=zlg*zfac
           if(ikind == 1 .and. snact(jv,jhor) == 1) then
            snu(:,jv,jhor)=zu0(:,jv)-(zu0(:,jv)-snu(:,jv,jhor))*zfac
           endif
          endif
          call snowtake(snfl(:,jv,jhor),zlf/zmass)
          call snowtake(sngr(:,jv,jhor),zlg/zmass)
          if(jv == 1) then
           zlfm=zlf
           asnow(jhor,snslot(1,jv))=asnow(jhor,snslot(1,jv))+zlf/1000.
           asnow(jhor,snslot(3,jv))=asnow(jhor,snslot(3,jv))+zlg/1000.
           zdq=(zlf+zlg)/zmass/dsigma(jlev)
           pdqdt(jhor)=pdqdt(jhor)+zdq
           pdtdt(jhor)=pdtdt(jhor)-zdq*zlcps
          else
           asnow(jhor,snslot(1,jv))=asnow(jhor,snslot(1,jv))+(zlf+zlg)/1000.
          endif
         enddo
        endif
        if(nsnowshad == 1 .and. pq(jhor) < pqsat(jhor)) then
         zlcpe=ALS/(acpd*(1.+ADV*dq(jhor,jlev)))
         do jk=2,1,-1
          if(snfl(jk,NSNV,jhor) > 0.) then
           zd=MIN(gamma*(pqsat(jhor)-pq(jhor))*dsigma(jlev)/deltsec2                                         &
     &        /(1.0+zlcpe*ra2*(TMELT-ra4)*pqsat(jhor)*pcor(jhor)/(pt(jhor)-ra4)**2),snfl(jk,NSNV,jhor))
           snfl(jk,NSNV,jhor)=snfl(jk,NSNV,jhor)-zd
           asnow(jhor,snslot(1,NSNV))=asnow(jhor,snslot(1,NSNV))+zd*zmass/1000.
          endif
         enddo
        endif
       endif
       pprsc(jhor)=snfl(1,1,jhor)+sngr(1,1,jhor)
       pprsl(jhor)=snfl(2,1,jhor)+sngr(2,1,jhor)
       snprevc(jhor)=pprsc(jhor)
       snprevl(jhor)=pprsl(jhor)
      enddo
      return
      end subroutine snowsubl
"""


# The snow's diagnostics in PlaSim's output (nsnowout, set by rainini when nevapprec 2 runs with nsnowsub 1): pumamod
# holds their sums since the last output and the steps they count, since PlaSim reads its restart before rainmod's
# namelist; outgp writes their means as codes 501-523 and outreset clears them.
OUTPUT_MARKER = "Terluna: the snow's diagnostics"
OUTPUT_MOD_EDITS = [
    ("      real :: aprs(NHOR)  = 0. ! acculumated snow fall\n",
     "      real :: aprs(NHOR)  = 0. ! acculumated snow fall\n"
     f"      integer, parameter :: NSNOWD = {NDIAG} ! {OUTPUT_MARKER} (rainmod's snowsubl, codes 501-{500 + NDIAG})\n"
     "      integer :: nsnowout = 0  ! 1 writes them\n"
     "      integer :: nsnacc = 0    ! their steps since the last output\n"
     "      real :: asnow(NHOR,NSNOWD) = 0. ! their sums since then\n"),
]
OUTPUT_OUT_EDITS = [
    ("      call writegp(40,aprs,144,0)\n",
     "      call writegp(40,aprs,144,0)\n"
     f"      if (nsnowout == 1) then   ! {OUTPUT_MARKER}, means since the last output\n"
     "       do jsnow=1,NSNOWD\n"
     "        asnow(:,jsnow)=asnow(:,jsnow)/real(max(nsnacc,1))\n"
     "        isnowcode=500+jsnow\n"
     "        call writegp(40,asnow(1,jsnow),isnowcode,0)\n"
     "       enddo\n"
     "      endif\n"),
    ("      aprs(:)=0.\n", f"      aprs(:)=0.\n      asnow(:,:)=0.   ! {OUTPUT_MARKER}\n      nsnacc=0\n"),
]
OUTPUT_PATCHES = [('plasimmod.f90', OUTPUT_MARKER, OUTPUT_MOD_EDITS), ('outmod.f90', OUTPUT_MARKER, OUTPUT_OUT_EDITS)]

# The runner's choices: the adopted snow with its shadows (the design run), its highest-forming-layer top as a
# sensitivity, and the adopted snow alone, for its cost.
CHOICES = {'flakes': dict(nsnowtop=1, nsnowshad=1), 'flakes_highest_top': dict(nsnowtop=2, nsnowshad=1),
           'flakes_no_shadows': dict(nsnowtop=1, nsnowshad=0)}


def setting(choice: str) -> dict:
    """The snow as a run's configuration records it, for a choice of CHOICES."""
    if choice not in CHOICES:
        raise ValueError(f'the snow from its flakes takes a choice of {sorted(CHOICES)}, not {choice!r}')
    import hashlib
    prm = dict(PARAMETERS, **CHOICES[choice])
    return dict(choice=choice, namelist={k: prm[k] for k in NAMELIST_NAMES},
                product=str(PRODUCT.relative_to(ROOT)),
                product_sha256=hashlib.sha256(PRODUCT.read_bytes()).hexdigest()[:16],
                variants=[v[0] for v in VARIANTS], air_factor=AIR_FACTOR,
                diagnostics={str(code): name for code, (name, _, _) in DIAGNOSTICS.items()},
                code='climate/gcm/snow_sublimation.py (snowprep, snowsubl and their Fortran, in rain_drops.py\'s '
                     'patch; OUTPUT_PATCHES)')
