"""The velocity change freight from the depot needs to reach each near-Earth source, and the water within each reach.

The survey's pools are windows of orbits. This product gives every body the survey quantifies a velocity change of
its own, so that the water a programme can reach follows from the velocity change its freight can afford.

For a catalogued body the velocity change is the approximate minimum of Shoemaker and Helin (1978, NASA CP-2053,
pp. 245-256): a transfer orbit tangent to the Earth's orbit and to the body's orbit at its aphelion, half the plane
change made at each end, speeds normalized to the Earth's orbital speed of 30 km/s and 0.5 km/s added to match computed
missions (their equations 2 to 11). Its form for a mission from low Earth orbit reproduces their Table 1 with low Earth
orbit at 200 km. The tugs' leg runs from a depot near the Earth's escape energy, its departure burned on a low pass
over the Earth: the figure from low Earth orbit without the escape from low orbit, with the same 0.5 km/s. The loaded
return is taken as its mirror, captured on a low pass. The figures assume ideal phasing and give no trajectory, launch
window or trip time.

The water is the survey's own: its Monte Carlo of the catalogued bodies is run again unchanged, and each body's water,
draw by draw, is added to the bin of its velocity change; the survey's undiscovered bodies are divided among the bins
by its orbital model, cell by cell. Quantities are contained water before recovery, as in the survey.
"""
from __future__ import annotations
import gzip, hashlib, io, json, math
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from shared.constants import EARTH_GM, EARTH_RADIUS
from . import composition, neo, population, protection

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
POPULATION = HERE / 'results/population.json.gz'
# Shoemaker and Helin's normalizing speed and correction (equation 11), and the low Earth orbit that reproduces their
# Table 1; the Earth's aphelion distance divides the Amor asteroids from the Apollos.
EARTH_SPEED_KM_S = 30.0
CORRECTION_KM_S = 0.5
LEO_ALTITUDE_KM = 200.0
EARTH_APHELION_AU = 1.017
# Upper edges of the bins of velocity change on the depot leg; the last bin holds bodies beyond the last edge.
EDGES_M_S = tuple(float(v) for v in range(100, 20001, 100))
PERCENTILES = (('p05', 0.05), ('p50', 0.5), ('p95', 0.95))
# Shoemaker and Helin's Table 1 with the elements printed in the companion paper of the same volume (pp. 161-176):
# body, a (AU), e, i (degrees), their figure of merit, their velocity change from low Earth orbit (km/s), and the
# rendezvous they assume: at perihelion for the short missions to Amors and Apollos, at aphelion for 1976AA-type bodies.
TABLE_1 = (('1943 Anteros', 1.430, 0.256, 8.70, 0.27, 8.7, 'perihelion'), ('1976AA', 0.966, 0.183, 18.91, 0.27, 8.7, 'aphelion'),
           ('433 Eros', 1.458, 0.223, 10.83, 0.29, 9.2, 'perihelion'), ('1620 Geographos', 1.244, 0.335, 13.33, 0.29, 9.3, 'perihelion'),
           ('1976UA', 0.844, 0.450, 5.85, 0.31, 9.7, 'aphelion'), ('1685 Toro', 1.367, 0.436, 9.37, 0.31, 9.7, 'perihelion'),
           ('1977VA', 1.865, 0.394, 2.97, 0.31, 9.8, 'perihelion'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def root(x):
    return np.sqrt(np.maximum(x, 0.0))


def speeds(a, e, i_deg, at='aphelion'):
    """Shoemaker and Helin's speeds in units of the Earth's orbital speed: departure from the Earth's neighbourhood U_t
    and the rendezvous impulse U_R. At aphelion, the minimum (their equations 3 to 10); at perihelion, the short
    missions to Amors and Apollos of their Table 1 (equations 12 to 17), 1976AA-type bodies keeping aphelion."""
    a, e, i = np.broadcast_arrays(np.asarray(a, float), np.asarray(e, float), np.radians(np.asarray(i_deg, float)))
    c = np.cos(i / 2)
    Q, q = a * (1 + e), a * (1 - e)
    aten = a <= 1.0
    amor = ~aten & (q > EARTH_APHELION_AU)
    if at == 'aphelion':
        ut2 = np.where(aten, 2 - 2 * root(2 * Q - Q * Q) * c, 3 - 2 / (Q + 1) - 2 * root(2 * Q / (Q + 1)) * c)
        uc2 = np.where(aten, 3 / Q - 1 - 2 / Q * root(2 - Q), 3 / Q - 2 / (Q + 1) - 2 / Q * root(2 / (Q + 1)) * np.where(amor, c, 1.0))
        ur2 = 3 / Q - 1 / a - 2 / Q * root(a / Q * (1 - e * e)) * np.where(amor, 1.0, c)
        sign = -1.0
    elif at == 'perihelion':
        if np.any(aten):
            raise ValueError('Perihelion rendezvous is defined here for Amors and Apollos only')
        ut2 = np.where(amor, 3 - 2 / (q + 1) - 2 * root(2 * q / (q + 1)) * c, 2 - 2 * root(2 * q - q * q) * c)
        uc2 = np.where(amor, 3 / q - 2 / (q + 1) - 2 / q * root(2 / (q + 1)), 3 / q - 1 - 2 / q * root(2 - q))
        ur2 = 3 / q - 1 / a - 2 / q * root(a / q * (1 - e * e)) * c
        sign = np.where(amor, 1.0, -1.0)
    else:
        raise ValueError('Rendezvous at aphelion or perihelion')
    uc, ur = root(uc2), root(ur2)
    return root(ut2), root(uc2 + sign * 2 * ur * uc * c + ur2)


def figure_of_merit(a, e, i_deg, at='aphelion'):
    """Shoemaker and Helin's figure of merit F = U_L + U_R from low Earth orbit (equations 1 and 2)."""
    ut, uR = speeds(a, e, i_deg, at)
    u0 = math.sqrt(EARTH_GM / (EARTH_RADIUS + 1000 * LEO_ALTITUDE_KM)) / 1000 / EARTH_SPEED_KM_S
    return np.sqrt(ut ** 2 + 2 * u0 ** 2) - u0 + uR


def leo_rendezvous_km_s(a, e, i_deg, at='aphelion'):
    """The approximate minimum velocity change from low Earth orbit to rendezvous (equation 11)."""
    return EARTH_SPEED_KM_S * figure_of_merit(a, e, i_deg, at) + CORRECTION_KM_S


def depot_leg_m_s(a, e, i_deg):
    """One leg between a depot near the Earth's escape energy and the body, at aphelion, with Shoemaker and Helin's
    correction. The departure burn falls on a low pass over the Earth, where a craft already at escape energy moves at
    the escape speed of low orbit, so it costs sqrt(U_t^2 + S^2) - S: the figure from low Earth orbit without the escape
    from low orbit. The change of the depot's orbit to that pass, and back, is not counted beyond the 0.5 km/s."""
    ut, uR = speeds(a, e, i_deg)
    escape = math.sqrt(2 * EARTH_GM / (EARTH_RADIUS + 1000 * LEO_ALTITUDE_KM)) / 1000 / EARTH_SPEED_KM_S
    return 1000 * (EARTH_SPEED_KM_S * (np.sqrt(ut ** 2 + escape ** 2) - escape + uR) + CORRECTION_KM_S)


def bin_of(dv_m_s):
    """The bin of each velocity change: the first upper edge not below it, or the last bin beyond the last edge."""
    return np.searchsorted(np.asarray(EDGES_M_S), np.asarray(dv_m_s, float), side='left')


def elements(record):
    values = [neo.number(record, k) for k in ('a', 'e', 'i')]
    return None if any(v is None for v in values) else values


class Recorder(population.Characterization):
    """The survey's characterization unchanged, which also adds each catalogued body's water, draw by draw, to the bin
    of its velocity change as the survey's Monte Carlo draws it."""

    def __init__(self, cfg, header, block_bins):
        super().__init__(cfg, header)
        self.block_bins = iter(block_bins)
        self.water = np.zeros((cfg['draws'], len(EDGES_M_S) + 1))

    def materials(self, diameter_km, albedo, group, rng):
        values = super().materials(diameter_km, albedo, group, rng)
        index = next(self.block_bins)
        for k in np.unique(index):
            self.water[:, k] += values[:, index == k, 1].sum(axis=1)
        return values


def known_records(records):
    """The bodies the survey quantifies, in its order: neither under a scientific hold nor without H and diameter."""
    holds = protection.load_inputs()
    return [r for r in records if not protection.exclusion_reasons(r, holds) and not (neo.number(r, 'H') is None and neo.number(r, 'diameter') is None)]


def grid_shares(raw, geometry, nodes=4):
    """For each H bin of the survey's orbital model, the share of its near-Earth orbits in each of the survey's three
    orbital categories (narrow, the rest of the wider screen, beyond) and each bin of velocity change, from Gauss
    points within each supplied cell, q < 1.3 AU as in the survey."""
    lines = raw.decode().splitlines()
    first = lambda i: lines[i].split('!')[0].split()
    nh, hmin, hmax = map(float, first(1))
    na, amin, amax = map(float, first(6))
    ne, emin, emax = map(float, first(7))
    ni, imin, imax = map(float, first(8))
    values = np.loadtxt(io.StringIO('\n'.join(lines[13:])), usecols=range(5))
    ih, ia, ie, ii = (values[:, k].astype(int) - 1 for k in range(4))
    weight = values[:, 4]
    da, de, di = (amax - amin) / na, (emax - emin) / ne, (imax - imin) / ni
    x, w = np.polynomial.legendre.leggauss(nodes)
    u, wu = (x + 1) / 2, w / 2
    out = np.zeros((int(nh), 3, len(EDGES_M_S) + 1))
    narrow, wider = geometry['narrow'], geometry['wider']
    for start in range(0, len(weight), 20000):
        s = slice(start, start + 20000)
        a = (amin + (ia[s, None, None, None] + u[None, :, None, None]) * da)
        e = (emin + (ie[s, None, None, None] + u[None, None, :, None]) * de)
        i = (imin + (ii[s, None, None, None] + u[None, None, None, :]) * di)
        a, e, i = np.broadcast_arrays(a, e, i)
        wt = weight[s, None, None, None] * wu[None, :, None, None] * wu[None, None, :, None] * wu[None, None, None, :]
        valid = (a * (1 - e) < 1.3) & (a > 0) & (e < 1)
        inner = (a >= narrow['a_min']) & (a <= narrow['a_max']) & (e <= narrow['e_max']) & (i <= narrow['i_max'])
        middle = (a >= wider['a_min']) & (a <= wider['a_max']) & (e <= wider['e_max']) & (i <= wider['i_max'])
        category = np.where(inner, 0, np.where(middle, 1, 2))
        k = bin_of(depot_leg_m_s(np.where(valid, a, 1.0), np.where(valid, e, 0.0), i))
        h = np.broadcast_to(ih[s, None, None, None], a.shape)
        np.add.at(out, (h[valid], category[valid], k[valid]), wt[valid])
    return [dict(H_low=hmin + j * (hmax - hmin) / nh, H_high=hmin + (j + 1) * (hmax - hmin) / nh, shares=out[j]) for j in range(int(nh))]


def undiscovered(survey, grid):
    """The survey's undiscovered water in each bin of velocity change, draw by draw, for its central population.

    Each orbital category's water, draw by draw, is divided among the bins in proportion to the expected undiscovered
    water of each H bin of the survey's counts and the orbital model's shares within that category, averaged over the
    model's H bins each count bin covers as the survey averages its category probabilities."""
    pools = survey['scenarios']['central']['pools']
    water = {p: np.asarray(pools[p]['undiscovered']['water_equivalent_kg']) for p in ('narrow', 'wider', 'all_nea')}
    by_category = [water['narrow'], water['wider'] - water['narrow'], water['all_nea'] - water['wider']]
    weights = np.zeros((3, len(EDGES_M_S) + 1))
    for b, moments in zip(survey['population_bins'], survey['moment_diagnostics']):
        if abs(b['H_low'] - moments['H_low']) > 1e-9:
            raise ValueError('The survey\'s count bins and moment diagnostics disagree')
        covered = [g['shares'] for g in grid if g['H_low'] >= b['H_low'] - 1e-9 and g['H_high'] <= b['H_high'] + 1e-9]
        if not covered:
            raise ValueError('Missing orbital bins')
        joint = np.mean([s / max(s.sum(), 1e-300) for s in covered], axis=0)
        conditional = joint / np.maximum(joint.sum(axis=1, keepdims=True), 1e-300)
        expected = np.asarray(b['missing_mean_by_branch'][1]) * moments['single_body_mean_kg'][1]
        weights += expected[:, None] * conditional
    fractions = weights / np.maximum(weights.sum(axis=1, keepdims=True), 1e-300)
    return sum(w[:, None] * f[None, :] for w, f in zip(by_category, fractions)), fractions


def quantiles(values):
    return {key: np.quantile(values, level, axis=0).tolist() for key, level in PERCENTILES}


def product(cfg=None):
    cfg = json.loads((HERE / 'population_inputs.json').read_text()) if cfg is None else cfg
    source = cfg['source_files']
    archive = ROOT / cfg['archive_directory']
    raw = population.checked(ROOT / source['catalogue']['path'], source['catalogue']['sha256'])
    records = composition.jpl_rows(raw)
    survey_table = population.checked(ROOT / source['taxonomy']['path'], source['taxonomy']['sha256'])
    matched, _ = composition.match_observations(records, composition.parse_table(survey_table))
    grid_raw = population.checked(archive / source['neomod_grid']['filename'], source['neomod_grid']['sha256'])
    header, _ = population.orbit_table(grid_raw, cfg['geometry'])
    survey = json.loads(gzip.decompress(POPULATION.read_bytes()))
    if survey['producer']['sha256'] != digest(HERE / 'population.py') or survey['random_seed'] != cfg['random_seed'] or survey['draws'] != cfg['draws']:
        raise ValueError('The survey product is not the one its producer and inputs give')
    # Each catalogued body's velocity change, in the survey's order and blocks.
    known = known_records(records)
    legs = np.array([depot_leg_m_s(*elements(r)) if elements(r) else np.inf for r in known])
    bins = bin_of(legs)
    size = cfg['catalogue_block_size']
    char = Recorder(cfg, header, [bins[s:s + size] for s in range(0, len(known), size)])
    rng = np.random.default_rng(cfg['random_seed'])
    pools, _ = population.catalogued(records, matched, cfg, char, rng)
    stated = survey['scenarios']['central']['pools']
    same = {p: bool(np.array_equal(pools[p][:, 1], np.asarray(stated[p]['catalogued']['water_equivalent_kg']))) for p in population.POOLS}
    if not all(same.values()):
        raise ValueError('The catalogued water differs from the survey\'s draws: ' + str(same))
    if not np.allclose(char.water.sum(axis=1), pools['all_nea'][:, 1], rtol=1e-12, atol=0):
        raise ValueError('The binned catalogued water does not add up to the survey\'s')
    grid = grid_shares(grid_raw, cfg['geometry'])
    missing, fractions = undiscovered(survey, grid)
    total = char.water + missing
    within = np.cumsum(total, axis=1)
    surveyed = np.asarray(stated['all_nea']['total']['water_equivalent_kg'])
    if not np.allclose(within[:, -1], surveyed, rtol=1e-12, atol=0):
        raise ValueError('The binned water does not add up to the survey\'s total')
    # Shoemaker and Helin's Table 1, from low Earth orbit, as the check of the formulas.
    table = []
    for name, a, e, i, F, dv, at in TABLE_1:
        f = float(figure_of_merit(a, e, i, at))
        table.append(dict(body=name, a_AU=a, e=e, i_deg=i, rendezvous_at=at, figure_of_merit=f, published_figure_of_merit=F,
                          leo_rendezvous_km_s=EARTH_SPEED_KM_S * f + CORRECTION_KM_S, published_km_s=dv))
    # NASA's table of trajectories for human round trips, joined by name, against the estimate from low Earth orbit.
    nhats_path = ROOT / 'research/runs/source_access/inputs/neo/nhats_summary.json'
    nhats = json.loads(nhats_path.read_text())
    named = {neo.name_key(d['fullname']): d for d in nhats['data']}
    pairs = [(elements(r), float(named[neo.name_key(r['full_name'])]['min_dv']['dv'])) for r in records
             if neo.name_key(r['full_name']) in named and elements(r)]
    ours = np.array([float(leo_rendezvous_km_s(*el)) for el, _ in pairs])
    theirs = np.array([dv for _, dv in pairs])
    # The largest bodies, by measured diameter or else by H at an albedo of 0.14, with their velocity changes.
    def size(r):
        d = neo.number(r, 'diameter')
        return d if d is not None else 1329 / math.sqrt(0.14) * 10 ** (-(neo.number(r, 'H') or 99.0) / 5)
    holders = []
    for j in sorted(range(len(known)), key=lambda j: -size(known[j]))[:40]:
        r = known[j]
        el = elements(r)
        holders.append(dict(name=r['full_name'].strip(), a_AU=el[0] if el else None, e=el[1] if el else None, i_deg=el[2] if el else None,
                            diameter_km=neo.number(r, 'diameter'), depot_leg_m_s=float(legs[j]) if np.isfinite(legs[j]) else None,
                            leo_rendezvous_km_s=float(leo_rendezvous_km_s(*el)) if el else None))
    deps = {str(p.relative_to(ROOT)): digest(p) for p in (Path(__file__), HERE / 'population.py', HERE / 'population_inputs.json', HERE / 'composition.py',
                                                          HERE / 'protection.py', HERE / 'protection_inputs.json', HERE / 'neo.py', POPULATION, nhats_path)}
    return {'schema': 'terluna.resources.source-reach/1', 'producer': {'path': str(Path(__file__).relative_to(ROOT)), 'sha256': digest(__file__)},
            'dependencies': deps, 'source_files': source,
            'method': {'velocity_change': 'Shoemaker and Helin (1978), NASA CP-2053 pp. 245-256, equations 2 to 11: rendezvous at aphelion, half the plane change at each end, speeds in units of 30 km/s, 0.5 km/s added. The depot leg is 30 km/s times (sqrt(U_t^2 + S^2) - S + U_R) plus 0.5 km/s, S the escape speed of low orbit: a departure burned on a low pass over the Earth from a depot near escape energy, and rendezvous; the loaded return its mirror.',
                       'low_Earth_orbit_altitude_km': LEO_ALTITUDE_KM, 'earth_speed_km_s': EARTH_SPEED_KM_S, 'correction_km_s': CORRECTION_KM_S,
                       'water': 'The survey\'s central population. Catalogued bodies: the survey\'s Monte Carlo run again unchanged, each body\'s water added to its bin draw by draw; the survey\'s own catalogued pools are reproduced draw by draw. Undiscovered bodies: each orbital category\'s water, draw by draw, divided among the bins by the orbital model\'s cells (4 Gauss points in each of a, e and i), weighted by the expected undiscovered water of each count bin. Bodies without orbital elements sit beyond the last edge.',
                       'limits': 'Ideal phasing at aphelion; no launch window, trip time, synodic period, gravity assist or aerocapture. Contained water before recovery; the survey\'s ranges are conditional, not reserves.'},
            'depot_leg_upper_edges_m_s': list(EDGES_M_S),
            'contained_water_within_kg': quantiles(within),
            'catalogued_water_within_kg': quantiles(np.cumsum(char.water, axis=1)),
            'undiscovered_water_within_kg': quantiles(np.cumsum(missing, axis=1)),
            'undiscovered_category_fractions': fractions.tolist(),
            'catalogued_bodies': len(known), 'catalogued_bodies_without_orbit': int(np.sum(~np.isfinite(legs))),
            'catalogued_pools_reproduced': same,
            'table_1_check': table,
            'nasa_human_trajectory_check': {'joined_bodies': len(pairs), 'spearman_rank_correlation': float(spearmanr(ours, theirs).correlation),
                                            'median_ratio_nasa_to_estimate': float(np.median(theirs / ours)),
                                            'share_nasa_above_estimate': float(np.mean(theirs > ours)),
                                            'reading': 'NASA\'s minimum for a human round trip of at most 450 days from low Earth orbit, with Earth return by entry, against the estimate of a one-way rendezvous from low Earth orbit at aphelion; a different mission, so the check is of rank and scale, not of value.'},
            'largest_bodies': holders,
            'evidence': 'Approximate velocity changes from a published figure of merit for ideal two-impulse transfers, and the survey\'s conditional water. No trajectory, source, reserve or route is qualified.',
            'qualified_reserve_kg': None}


def main():
    x = product()
    path = HERE / 'results/reach.json.gz'
    path.write_bytes(gzip.compress((json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode(), mtime=0))
    print(path, flush=True)


if __name__ == '__main__':
    main()
