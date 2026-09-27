"""Canopy photosynthesis under the Open Moon's clear sky and Earth's, from the surface-light product.

    python -m biosphere.canopy.model        # writes results/canopy.json (a few minutes, one thread)

The light is the illumination domain's surface-light product: direct and diffuse spectra
at the ground for the design Moon (1.2 atm, 294 K, behind the titania film dimmed 5%)
and the middle atmosphere's Earth control, at Sun heights of 1-90 degrees. The canopy
reflects some of it back to the sky, which returns part (the product's reflectance from
below). Inside the canopy the light is followed by wavelength from 400 to 760 nm
(radiation.py) with PROSPECT-D leaves (leaf_optics.py), and each layer's sunlit and
shaded leaves photosynthesise by the C3 leaf model (leaf.py). A sunlit leaf's share of
the direct beam is spread over its orientations, uniform in the cosine for spherical
leaves. Rubisco capacity falls with depth as exp(-kn L), with kn from Lloyd et al.
(2010). Photons from 400 to 700 nm drive photosynthesis; the ePAR variant also counts
700-750 nm, which Zhen and Bugbee (2020) found as effective when mixed with shorter
wavelengths.

Both worlds hold O2 at Earth's partial pressure and 400 ppm of CO2: 48.6 Pa on the Moon
and 40.5 Pa on Earth. Leaves are at 22 C by day and respire at 19 C by night, the
equatorial land air of the chosen climate at its warmest and coldest (climate/gcm run
A28_dim5); Earth gets the same temperatures so that light and CO2 are compared alone.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np

from shared.constants import AVOGADRO, PLANCK, SPEED_OF_LIGHT, SYNODIC_MONTH_DAYS
from biosphere.long_night import periodic_storage_requirement
from biosphere.canopy import fetch_inputs, leaf, leaf_optics as lo, radiation as rad

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.biosphere.canopy-photosynthesis/1'
LIGHT = ROOT / 'illumination' / 'surface_light' / 'results' / 'surface_light.json'
LIGHT_SCHEMA = 'terluna.illumination.surface-light/1'

BAND_NM = (400.0, 760.0)
PAR_NM = (400.0, 700.0)
EPAR_NM = (400.0, 750.0)
O2_PA = 21227.0
WORLDS = dict(moon=dict(case='moon_1.2atm/design', co2_pa=400e-6 * 121590.0, solar_day_s=SYNODIC_MONTH_DAYS * 86400.0),
              earth=dict(case='earth_control/unfiltered', co2_pa=400e-6 * 101325.0, solar_day_s=86400.0))
DAY_C, NIGHT_C = 22.0, 19.0
LATITUDES = (0.0, 30.0, 60.0)
LAI_SWEEP = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
GRAMS_C_PER_UMOL = 12.011e-6
PALE = lo.Leaf(chlorophyll=20.0, carotenoids=4.0)
U_NODES, U_WEIGHTS = (lambda x, w: (0.5 * (x + 1), 0.5 * w))(*np.polynomial.legendre.leggauss(8))


@dataclass(frozen=True)
class Stand:
    """A canopy: its structure, the leaves in it and what counts as photosynthetic light."""
    lai: float = 5.0
    clumping: float = 1.0
    soil: float = 0.1                 # Lambertian soil reflectance, spectrally neutral
    vcmax_top: float = 60.0
    pale_top_lai: float = 0.0         # leaf area at the top of the canopy carrying the pale leaf
    band_nm: tuple = PAR_NM

    def canopy(self):
        return rad.Canopy(lai=self.lai, clumping=self.clumping)

    def optics_key(self):
        return self.lai, self.clumping, self.soil, self.pale_top_lai


def kn(vcmax_top):
    """Decline of Rubisco capacity per unit leaf area index (Lloyd et al. 2010)."""
    return math.exp(0.00963 * vcmax_top - 2.43)


def surface_light():
    """Photon spectra (umol m-2 s-1 per nm) at the ground for each world, and the product's hash."""
    raw = LIGHT.read_bytes()
    data = json.loads(raw)
    if data['schema'] != LIGHT_SCHEMA:
        raise ValueError('Unexpected surface-light product schema')
    lam = np.asarray(data['centres_nm'])
    sel = (lam > BAND_NM[0]) & (lam < BAND_NM[1])
    to_umol = lam[sel] * 1e-9 / (PLANCK * SPEED_OF_LIGHT) / AVOGADRO * 1e6
    worlds = {}
    for name, world in WORLDS.items():
        case = data['cases'][world['case']]
        keys = list(case['spectra'])
        worlds[name] = dict(
            heights=np.array([float(k) for k in keys]),
            direct=np.array([case['spectra'][k]['direct_w_m2_nm'] for k in keys])[:, sel] * to_umol,
            diffuse=np.array([case['spectra'][k]['diffuse_black_w_m2_nm'] for k in keys])[:, sel] * to_umol,
            reflectance_from_below=np.asarray(data['columns'][case['column']]['reflectance_from_below'])[sel])
    return lam[sel], worlds, hashlib.sha256(raw).hexdigest()[:16]


class Optics:
    """Canopy responses to unit direct light at each Sun height and to unit diffuse light, cached per stand."""

    def __init__(self, lam, heights):
        self.lam, self.heights = lam, heights
        coeffs = lo.coefficients()
        self.leaves = {leaf_: tuple(np.interp(lam, lo.WAVELENGTH_NM, x) for x in lo.optics(leaf_, coeffs))
                       for leaf_ in (lo.Leaf(), PALE)}
        self.cache = {}

    def leaf_spectra(self, stand: Stand):
        canopy = stand.canopy()
        depth = canopy.depth()
        top = depth < stand.pale_top_lai
        (r0, t0), (r1, t1) = self.leaves[lo.Leaf()], self.leaves[PALE]
        rho = np.where(top[:, None], r1[None], r0[None])
        tau = np.where(top[:, None], t1[None], t0[None])
        return rho, tau

    def __call__(self, stand: Stand):
        key = stand.optics_key()
        if key not in self.cache:
            rho, tau = self.leaf_spectra(stand)
            soil = np.full(self.lam.size, stand.soil)
            canopy = stand.canopy()
            self.cache[key] = dict(
                diffuse=rad.solve(canopy, rho, tau, soil),
                direct=[rad.solve(canopy, rho, tau, soil, mu0=math.sin(math.radians(h))) for h in self.heights])
        return self.cache[key]


def rates(stand: Stand, resp, lam, direct, diffuse, refl, ca_pa, t_c=DAY_C):
    """Canopy gross photosynthesis and its light at one Sun height (per unit ground area).

    direct and diffuse are the incident photon spectra over a black ground; the sky returns
    part of what the canopy reflects: F = (diffuse + R A_dir direct) / (1 - R A_dif).
    """
    d_resp, f_resp = resp
    band = (lam > stand.band_nm[0]) & (lam < stand.band_nm[1])
    par = (lam > PAR_NM[0]) & (lam < PAR_NM[1])
    sky = (diffuse + refl * d_resp['reflected'] * direct) / (1 - refl * f_resp['reflected'])
    shade = (direct * d_resp['absorbed_diffuse'] + sky * f_resp['absorbed_diffuse'])[:, band].sum(axis=1)
    facing = (direct * d_resp['direct_on_facing_leaf'])[:, band].sum(axis=1)
    sunlit = d_resp['sunlit']
    canopy = stand.canopy()
    n, dl = canopy.layers()
    vcmax25 = stand.vcmax_top * np.exp(-kn(stand.vcmax_top) * canopy.depth())
    a_shade = leaf.gross(shade, vcmax25, t_c, ca_pa, O2_PA)
    a_sun = sum(w * leaf.gross(shade + u * facing, vcmax25, t_c, ca_pa, O2_PA) for u, w in zip(U_NODES, U_WEIGHTS))
    gpp = float(np.sum(dl * ((1 - sunlit) * a_shade + sunlit * a_sun)))
    absorbed = float(np.sum(dl * (shade + 0.5 * sunlit * facing)))
    incident = direct + sky
    reflected = d_resp['reflected'] * direct + f_resp['reflected'] * sky
    return dict(gpp_umol_m2_s=gpp, absorbed_umol_m2_s=absorbed,
                incident_par_umol_m2_s=float(incident[par].sum()),
                incident_band_umol_m2_s=float(incident[band].sum()),
                sky_return_par=float(sky[par].sum() / diffuse[par].sum() - 1) if diffuse[par].sum() > 0 else 0.0,
                canopy_albedo_par=float(reflected[par].sum() / incident[par].sum()),
                sunlit_lai=float(np.sum(sunlit * dl)),
                diffuse_share_par=float(sky[par].sum() / incident[par].sum()))


def respiration(stand: Stand, t_c):
    """Leaf respiration of the canopy (umol CO2 m-2 s-1 of ground)."""
    canopy = stand.canopy()
    n, dl = canopy.layers()
    vcmax25 = stand.vcmax_top * np.exp(-kn(stand.vcmax_top) * canopy.depth())
    return float(np.sum(dl * leaf.respiration(vcmax25, t_c)))


def cycle_mean(heights_deg, values, latitude_deg):
    """Mean over a whole solar cycle, day and night, with the Sun on the celestial equator.

    Values are interpolated linearly in the sine of the Sun's height and taken to fall
    linearly to zero between the lowest height given and the horizon.
    """
    x = np.sin(np.radians(np.asarray(heights_deg, dtype=float)))
    order = np.argsort(x)
    xp = np.concatenate([[0.0], x[order]])
    fp = np.concatenate([[0.0], np.asarray(values, dtype=float)[order]])
    hour = np.linspace(-math.pi / 2, math.pi / 2, 20001)
    sin_h = math.cos(math.radians(latitude_deg)) * np.cos(hour)
    return float(np.trapezoid(np.interp(sin_h, xp, fp), hour)) / (2 * math.pi)


def carbon_cycle(stand: Stand, heights, gpp, latitude_deg, solar_day_s, steps=2880):
    """Leaf carbon over one solar cycle (g C per m2 of ground): photosynthesis, day and night leaf
    respiration, and the smallest store that carries the leaves through the cycle."""
    t = (np.arange(steps) + 0.5) / steps
    sin_h = math.cos(math.radians(latitude_deg)) * np.cos(2 * math.pi * t - math.pi)
    x = np.sin(np.radians(heights))
    order = np.argsort(x)
    g = np.interp(np.maximum(sin_h, 0.0), np.concatenate([[0.0], x[order]]), np.concatenate([[0.0], np.asarray(gpp)[order]]))
    day = sin_h > 0
    r_day, r_night = respiration(stand, DAY_C), respiration(stand, NIGHT_C)
    net = np.where(day, g - r_day, -r_night)
    dt = solar_day_s / steps
    store = periodic_storage_requirement(net, dt)
    to_g = GRAMS_C_PER_UMOL
    return dict(gpp_g_c_m2=float(np.sum(g) * dt * to_g),
                day_respiration_g_c_m2=float(np.sum(day) * dt * r_day * to_g),
                night_respiration_g_c_m2=float(np.sum(~day) * dt * r_night * to_g),
                net_g_c_m2=float(store['cycle_net'] * to_g),
                balance_closes=bool(store['cycle_balance_feasible']),
                store_for_the_night_g_c_m2=float(store['worst_single_cycle_deficit'] * to_g),
                longest_deficit_hours=float(store['deficit_duration_days'] / 3600.0))


def table(stand: Stand, optics: Optics, lam, light, co2_pa):
    """Rates at every Sun height of a light source (heights, direct, diffuse, reflectance_from_below)."""
    resp = optics(stand)
    return [rates(stand, (resp['direct'][i], resp['diffuse']), lam, light['direct'][i], light['diffuse'][i],
                  light['reflectance_from_below'], co2_pa) for i in range(len(light['heights']))]


def summary(stand: Stand, heights, rows, world):
    """Means over the whole solar cycle at each latitude, and the leaf carbon over one cycle."""
    gpp = [row['gpp_umol_m2_s'] for row in rows]
    day_s = WORLDS[world]['solar_day_s']
    return {str(lat): dict(
        gpp_mean_umol_m2_s=cycle_mean(heights, gpp, lat),
        gpp_g_c_m2_per_24h=cycle_mean(heights, gpp, lat) * 86400.0 * GRAMS_C_PER_UMOL,
        par_mean_umol_m2_s=cycle_mean(heights, [row['incident_par_umol_m2_s'] for row in rows], lat),
        absorbed_mean_umol_m2_s=cycle_mean(heights, [row['absorbed_umol_m2_s'] for row in rows], lat),
        cycle=carbon_cycle(stand, heights, gpp, lat, day_s)) for lat in LATITUDES}


EVIDENCE = ('A one-dimensional, horizontally uniform canopy under clear skies, driven by the surface-light '
            'product. Leaves are PROSPECT-D green broadleaf leaves with a spherical angle distribution; leaf '
            'photosynthesis is the steady-state C3 model with typical broadleaf capacities (Vcmax 60 umol m-2 s-1 '
            'at the top, falling with depth); internal CO2 is a fixed 0.7 of ambient, so stomata, water supply '
            'and the slower diffusion of CO2 in denser air are not modelled, and leaves are at the air '
            'temperature. There is no acclimation to continuous light, no feedback from full sugar stores, no '
            'photoinhibition, no respiration of stems and roots and no growth. The surface light understates '
            'light at Sun heights below about 20 degrees (its atlas comparison). This compares how the two skies '
            'drive the same canopy; it is not a yield forecast.')
READING_RULE = ('by_sun_height[world] holds, per Sun height, canopy gross photosynthesis (umol CO2 m-2 s-1 of ground), '
                'the photosynthetic photons it absorbs and those incident (with the sky\'s return of what the canopy '
                'reflects), the canopy\'s albedo and diffuse share in 400-700 nm, and its sunlit leaf area. '
                'summaries[world][latitude] are means over the whole solar cycle, day and night, with the Sun on '
                'the celestial equator; gpp_g_c_m2_per_24h converts the mean to grams of carbon per 24 hours. '
                'cycle is the leaves\' carbon over one solar cycle (the Moon\'s is the synodic month, Earth\'s an '
                'equinox day): photosynthesis, leaf respiration by day and night, the net, and the smallest store '
                'that carries the leaves through the cycle (biosphere/long_night.py). decomposition holds the '
                'cycle-mean photosynthesis of the base stand under light that changes from Earth\'s to the Moon\'s '
                'one factor at a time. lai and variants hold equatorial summaries.')


def run(verbose=True):
    say = print if verbose else (lambda *a, **k: None)
    lam, light, light_sha = surface_light()
    heights = light['moon']['heights']
    if not np.array_equal(heights, light['earth']['heights']):
        raise ValueError('The two worlds must share Sun heights')
    optics = Optics(lam, heights)
    base = Stand()
    out = dict(by_sun_height={}, summaries={}, decomposition={}, lai={}, variants={})
    t0 = time.time()
    rows = {w: table(base, optics, lam, light[w], WORLDS[w]['co2_pa']) for w in WORLDS}
    for w in WORLDS:
        out['by_sun_height'][w] = [dict(sun_deg=float(h), **r) for h, r in zip(heights, rows[w])]
        out['summaries'][w] = summary(base, heights, rows[w], w)
    say(f'base stand: {time.time() - t0:.0f} s', flush=True)

    # From Earth's light to the Moon's one factor at a time: CO2, the number of photons (Earth's spectrum
    # and diffuse share scaled at each Sun height), the Moon's spectrum with Earth's diffuse share at each
    # wavelength, and the Moon's own diffuse share.
    moon, earth = light['moon'], light['earth']
    par = (lam > PAR_NM[0]) & (lam < PAR_NM[1])
    moon_total, earth_total = moon['direct'] + moon['diffuse'], earth['direct'] + earth['diffuse']
    scale = (moon_total[:, par].sum(axis=1) / earth_total[:, par].sum(axis=1))[:, None]
    earth_direct_share = earth['direct'] / np.maximum(earth_total, 1e-300)
    moon_co2, earth_co2 = WORLDS['moon']['co2_pa'], WORLDS['earth']['co2_pa']
    sources = dict(
        earth=(earth, earth_co2),
        earth_light_moon_co2=(earth, moon_co2),
        earth_light_scaled_to_moon_photons=(dict(earth, direct=earth['direct'] * scale, diffuse=earth['diffuse'] * scale), moon_co2),
        moon_light_earth_diffuse_share=(dict(moon, direct=moon_total * earth_direct_share,
                                             diffuse=moon_total * (1 - earth_direct_share)), moon_co2),
        moon=(moon, moon_co2))
    for name, (src, co2) in sources.items():
        r = table(base, optics, lam, src, co2)
        out['decomposition'][name] = {str(lat): cycle_mean(heights, [x['gpp_umol_m2_s'] for x in r], lat)
                                      for lat in LATITUDES}

    for lai in LAI_SWEEP:
        t0 = time.time()
        stand = replace(base, lai=lai)
        for w in WORLDS:
            out['lai'].setdefault(w, {})[str(lai)] = summary(stand, heights, table(stand, optics, lam, light[w],
                                                                                    WORLDS[w]['co2_pa']), w)['0.0']
        say(f'leaf area index {lai:g}: {time.time() - t0:.0f} s', flush=True)

    variants = dict(extended_par=replace(base, band_nm=EPAR_NM), pale_top_third=replace(base, pale_top_lai=base.lai / 3),
                    clumped=replace(base, clumping=0.7), vcmax_40=replace(base, vcmax_top=40.0),
                    vcmax_90=replace(base, vcmax_top=90.0))
    for name, stand in variants.items():
        t0 = time.time()
        for w in WORLDS:
            out['variants'].setdefault(name, {})[w] = summary(stand, heights, table(stand, optics, lam, light[w],
                                                                                    WORLDS[w]['co2_pa']), w)['0.0']
        say(f'{name}: {time.time() - t0:.0f} s', flush=True)

    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    prospect = next(f['sha256'][:16] for f in fetch_inputs.manifest()['files'] if f['name'] == 'prospect_d_spectra.txt')
    return dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', files={p: digest(p) for p in (
            'biosphere/canopy/model.py', 'biosphere/canopy/radiation.py', 'biosphere/canopy/leaf.py',
            'biosphere/canopy/leaf_optics.py', 'biosphere/long_night.py')},
            inputs=dict(surface_light_sha256=light_sha, surface_light_schema=LIGHT_SCHEMA,
                        prospect_d_sha256=prospect)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        units=dict(rates='umol CO2 m-2 s-1 of ground', photons='umol m-2 s-1', carbon='g C m-2 of ground'),
        settings=dict(stand=base.__dict__, leaf_traits=leaf.Traits().__dict__, leaf=lo.Leaf().__dict__,
                      pale_leaf=PALE.__dict__, kn=kn(base.vcmax_top), day_leaf_c=DAY_C, night_leaf_c=NIGHT_C,
                      o2_pa=O2_PA, co2_pa={w: WORLDS[w]['co2_pa'] for w in WORLDS},
                      surface_light_cases={w: WORLDS[w]['case'] for w in WORLDS}, latitudes_deg=LATITUDES,
                      lai_sweep=LAI_SWEEP, layer_lai=base.canopy().layer_lai, streams=base.canopy().streams),
        **out)


def _round(obj, digits=6):
    if isinstance(obj, float):
        return float(f'{obj:.{digits}g}')
    if isinstance(obj, dict):
        return {k: _round(v, digits) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_round(v, digits) for v in obj]
    return obj


def main():
    product = _round(run())
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / 'canopy.json').write_text(json.dumps(product, indent=1) + '\n')
    print('wrote', RESULTS / 'canopy.json')


if __name__ == '__main__':
    main()
