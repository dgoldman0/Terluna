"""The solar shield's ring fleet as a light in the Moon's night: what its night-side tiles send onto the night
hemisphere, by diffuse scatter and by mirror reflection.

    python -m illumination.fleet_light.model
    # -> illumination/fleet_light/results/night_light.json

The lead ring fleet (research/studies/solar_shield_array) keeps every tile facing the Moon round its whole orbit.
Behind the Moon that face is the one the Sun lights, so the night half of the fleet shows the Moon its sunlit side.
The model builds the fleet from the shield's own products: the 1,859 rings' crossing heights over the 8,625 km
stack, each ring's radius by height and its forced eccentricity (apolune toward the Sun), the shared node line at
right angles to the Sun in the Moon's orbit plane, and the window and annulus films' visible optics. Each ring is a
continuous strip of tiles 10 km wide. The Sun reaches a night-side tile through the day-side layers that project onto
it and the night-side layers in front of it (their counts from the projected strip area), and the Moon's umbra and
penumbra darken the tiles behind it.

Two kinds of light reach the night hemisphere:
- diffuse scatter from the Moon-facing face, Lambertian with reflectance rho_d (requirement O5 sets about 0.1% for
  Sun-facing surfaces); the illuminance on level ground is integrated over every lit strip above the horizon;
- mirror reflection, which a tile facing exactly the Moon's centre returns at the impact parameter it came in on, so
  it misses the Moon; tilting the tile turns the reflected ray through twice the tilt. The model tilts each tile by a
  fixed angle in a random direction (or out of the ring's plane only) and counts the reflected light that lands.

The Moon's orbit plane stands for its equator (they differ by up to 6.7 degrees), and the observer's angle from the
anti-solar point stands for local time from midnight.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import MOON_RADIUS, SUN_RADIUS, AU

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.illumination.fleet-night-light/1'
LAYOUT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'ring_layout.json'
FILMS = ROOT / 'protection' / 'spectra' / 'annulus_film.json'
OUT = HERE / 'results' / 'night_light.json'
SUN_LUX = 128_000.0                   # sunlight above the air at 1 AU (about 1,361 W/m2 at 94 lm/W)
TILE_WIDTH_KM = 10.0
WINDOW_RADIUS_KM = 1830.0             # rings crossing below this height carry the window stack (integrated ledger)
RHO_D = 1e-3                          # O5's diffuse reflectance for Sun-facing surfaces; results scale with it
TILTS_DEG = (0.0, 1.5, 3.0, 6.0)
ROLLS_DEG = (0.0, 30.0, 60.0, 80.0)       # the shield's clearance-safe rolls outside service
RINGS_MODELLED = 465                  # every fourth ring, each weighted for four
THETA_STEPS = 1440                    # samples round each ring
ALPHA_DEG = np.arange(0.0, 90.01, 5.0)     # observer's angle from the anti-solar point
BETA_DEG = np.arange(-80.0, 80.01, 10.0)   # observer's latitude from the orbit plane
BIN_KM = 50.0
EVIDENCE = ('Geometric optics of a model fleet built from the shield branch\'s layout and film products: continuous '
            'strips of tiles, Lambertian diffuse scatter at a stated reflectance, specular reflection with fixed tilt '
            'angles, shading by projected layer counts and the Moon\'s shadow. No attitude law of the fleet has been '
            'flown on the night side; the tilts stand in for its keeping and pitch envelope.')
READING_RULE = ('Illuminance in lux on level ground. diffuse_lux[alpha][beta] is at rho_d; it scales in proportion '
                'to rho_d. alpha is the observer\'s angle from the anti-solar point (0 at local midnight, 90 at the '
                'terminator), beta the latitude from the Moon\'s orbit plane. specular[tilt] gives the share of the '
                'night half\'s mirror-reflected light that lands on the Moon and the mean illuminance it makes over '
                'the night hemisphere. Luminous fluxes in lumens; 1 PW of sunlight is about 0.94e17 lm.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def fleet():
    """Ring crossing heights (km), semi-major axes (km), eccentricities, and the visible optics of each ring's tiles."""
    layout = json.loads(LAYOUT.read_text())
    stack = layout['stack']
    top = stack['height_top_km']
    heights = np.linspace(-top, top, stack['rings'])[::stack['rings'] // RINGS_MODELLED][:RINGS_MODELLED]
    weight = stack['rings'] / heights.size
    by = stack['ideal_radius_km']['by_height']
    a = np.interp(heights, [r['height_km'] for r in by], [r['radius_km'] for r in by])
    flown = layout['flown']['matched']['rings']
    e = np.interp(heights, [r['height_km'] for r in flown], [r['forced_eccentricity'] for r in flown])
    films = json.loads(FILMS.read_text())
    window = films['window_stack']
    annulus = next(f for f in films['films'] if f['coated'] and f['titania_um'] == 0.1 and f['silica_um'] == 4.0)
    optics = {name: dict(R=f['visible_reflectance'], T=f['band_transmission']['visible_400_700'])
              for name, f in (('window', window), ('annulus', annulus))}
    is_window = np.abs(heights) < WINDOW_RADIUS_KM
    return heights, a, e, weight, is_window, optics


def strips(heights, a, e):
    """Positions (km) and unit outward normals of every ring sample; the Sun lies along -x, the node line along y."""
    sin_t = np.clip(heights / a, -1.0, 1.0)
    cos_t = np.sqrt(1.0 - sin_t ** 2)
    theta = (np.arange(THETA_STEPS) + 0.5) * 2.0 * np.pi / THETA_STEPS
    r = (a * (1 - e ** 2))[:, None] / (1.0 + e[:, None] * np.cos(theta)[None, :])      # perilune at theta = 0, night side
    u = np.stack([cos_t, np.zeros_like(cos_t), -sin_t], axis=-1)                         # night crossing direction
    y = np.array([0.0, 1.0, 0.0])
    P = r[..., None] * (np.cos(theta)[None, :, None] * u[:, None, :] + np.sin(theta)[None, :, None] * y)
    dtheta = 2.0 * np.pi / THETA_STEPS
    dP = np.gradient(P, dtheta, axis=1)
    length = np.linalg.norm(dP, axis=-1) * dtheta                                        # km of ring per sample
    n = P / np.linalg.norm(P, axis=-1, keepdims=True)
    tangent = dP / np.linalg.norm(dP, axis=-1, keepdims=True)
    return P, n, length, tangent


def coverage(P, n, length, weight, night):
    """Layers of strip projected onto the plane facing the Sun (count per area), for the day and night halves."""
    edges_y = np.arange(-23000.0, 23000.01, BIN_KM)
    edges_z = np.arange(-10000.0, 10000.01, BIN_KM)
    area = TILE_WIDTH_KM * length * np.abs(n[..., 0]) * weight
    out = {}
    for name, mask in (('day', ~night), ('night', night)):
        h, _, _ = np.histogram2d(P[..., 1][mask], P[..., 2][mask], bins=(edges_y, edges_z), weights=area[mask])
        out[name] = h / BIN_KM ** 2
    def lookup(name, yy, zz):
        iy = np.clip(np.searchsorted(edges_y, yy) - 1, 0, edges_y.size - 2)
        iz = np.clip(np.searchsorted(edges_z, zz) - 1, 0, edges_z.size - 2)
        return out[name][iy, iz]
    return lookup


def moon_shadow(P):
    """Fraction of the Sun's disc a point sees past the Moon (umbra 0, outside the penumbra 1, linear between)."""
    b = np.hypot(P[..., 1], P[..., 2])
    behind = np.maximum(P[..., 0], 0.0)
    sun_angle = SUN_RADIUS / AU
    moon_r = MOON_RADIUS / 1e3
    inner, outer = moon_r - behind * sun_angle, moon_r + behind * sun_angle
    lit = np.clip((b - inner) / np.maximum(outer - inner, 1e-9), 0.0, 1.0)
    return np.where(P[..., 0] > 0, lit, 1.0)


def model() -> dict:
    heights, a, e, weight, is_window, optics = fleet()
    P, n, length, tangent = strips(heights, a, e)
    night = P[..., 0] > 0.0
    lookup = coverage(P, n, length, weight, night)
    R = np.where(is_window, optics['window']['R'], optics['annulus']['R'])[:, None] * np.ones_like(length)
    T = np.where(is_window, optics['window']['T'], optics['annulus']['T'])[:, None] * np.ones_like(length)
    yy, zz = P[..., 1], P[..., 2]
    c_day, c_night = lookup('day', yy, zz), lookup('night', yy, zz)
    in_front = np.maximum(c_night - 1.0, 0.0) / 2.0                     # night layers between a tile and the Sun, on average
    E_in = SUN_LUX * moon_shadow(P) * T ** (c_day + in_front)              # lux on the tile's sunlit face
    mu_i = np.clip(n[..., 0], 0.0, None)                                   # the Moon-facing face's cosine to the Sun, night side
    dA_m2 = TILE_WIDTH_KM * length * weight * 1e6
    lit = night & (mu_i > 0)
    flux_in = float((E_in * mu_i * dA_m2)[lit].sum())                      # lumens onto the night half
    L = (RHO_D * E_in * mu_i / np.pi)[lit]                                 # cd/m2 of the Moon-facing face
    Pl, nl, dAl, tl = P[lit], n[lit], dA_m2[lit], tangent[lit]
    see_through = (T ** in_front)[lit]                                     # nearer night layers the scatter passes back through
    moon_km = MOON_RADIUS / 1e3
    diffuse = np.zeros((ALPHA_DEG.size, BETA_DEG.size))
    for i, al in enumerate(np.radians(ALPHA_DEG)):
        for j, be in enumerate(np.radians(BETA_DEG)):
            up = np.array([np.cos(be) * np.cos(al), np.cos(be) * np.sin(al), np.sin(be)])
            D = Pl - moon_km * up
            d = np.linalg.norm(D, axis=-1)
            mu_o = D @ up / d
            mu_e = np.einsum('ij,ij->i', D, nl) / d
            ok = (mu_o > 0) & (mu_e > 0)
            diffuse[i, j] = float((L * see_through * mu_e * mu_o * dAl / (d * 1e3) ** 2)[ok].sum())
    zenith_luminance = float(np.max(L))
    # Mirror reflection: tilt each lit night tile by a fixed angle and see where the ray goes.
    rng = np.random.default_rng(20261008)
    xhat = np.array([1.0, 0.0, 0.0])
    specular = {}
    flux_spec = (R * E_in * mu_i * dA_m2)[lit]
    hemisphere_m2 = 2.0 * np.pi * MOON_RADIUS ** 2
    for mode in ('random', 'out_of_plane'):
        for tilt in TILTS_DEG:
            tilt_r = np.radians(tilt)
            if mode == 'random':                                             # any direction about the normal
                k = rng.normal(size=nl.shape); k -= np.einsum('ij,ij->i', k, nl)[:, None] * nl
            else:                                                            # out of the ring's plane, either way
                k = np.cross(nl, tl) * rng.choice([-1.0, 1.0], size=(nl.shape[0], 1))
            k /= np.maximum(np.linalg.norm(k, axis=-1, keepdims=True), 1e-12)
            m = np.cos(tilt_r) * nl + np.sin(tilt_r) * k                     # the tilted outward normal
            face = -m                                                        # the sunlit, Moon-facing face
            ray = xhat - 2.0 * np.einsum('j,ij->i', xhat, face)[:, None] * face
            # distance of closest approach of the reflected ray to the Moon's centre, if it heads toward it
            t_close = -np.einsum('ij,ij->i', Pl, ray)
            closest = np.linalg.norm(Pl + np.maximum(t_close, 0.0)[:, None] * ray, axis=-1)
            lands = (t_close > 0) & (closest < moon_km)
            share = float(flux_spec[lands].sum() / flux_spec.sum())
            # where it lands: the first crossing of the Moon's sphere, as angle from the anti-solar point and latitude
            Pk, rk = Pl[lands], ray[lands]
            proj = -np.einsum('ij,ij->i', Pk, rk)
            hit = Pk + (proj - np.sqrt(np.maximum(moon_km ** 2 - (np.linalg.norm(Pk + proj[:, None] * rk, axis=-1)) ** 2, 0.0)))[:, None] * rk
            hit_alpha = np.degrees(np.arccos(np.clip(hit[:, 0] / moon_km, -1, 1)))
            cos_inc = np.clip(-np.einsum('ij,ij->i', rk, hit / moon_km), 0.0, 1.0)
            edges = np.arange(0.0, 90.01, 10.0)
            ring_area = 2 * np.pi * MOON_RADIUS ** 2 * (np.cos(np.radians(edges[:-1])) - np.cos(np.radians(edges[1:])))
            landed, _ = np.histogram(hit_alpha, bins=edges, weights=flux_spec[lands])
            specular[f'{mode}_{tilt:g}'] = dict(tilt_deg=tilt, mode=mode, landing_share=round(share, 6),
                                                mean_lux_over_night_hemisphere=round(float(flux_spec[lands].sum()) / hemisphere_m2, 4),
                                                lux_by_alpha_band=dict(edges_deg=edges.tolist(),
                                                                       lux=np.round(landed / ring_area, 3).tolist()),
                                                landed_mean_cos_incidence=round(float(np.average(cos_inc, weights=flux_spec[lands])) if lands.any() else 0.0, 3),
                                                tiles_that_must_tilt_away_share=round(float(
                                                    (np.hypot(Pl[:, 1], Pl[:, 2]) < moon_km + 2 * tilt_r * np.linalg.norm(Pl, axis=-1)).sum() / Pl.shape[0]), 4))
    # Rolling the night half's tiles about the line through their ring's tile centres (the shield's clearance-safe
    # turn; integrated_comparison.md) tips each Moon-facing face away from both the Sun and the Moon. Half the rings
    # roll each way.
    rolls = {}
    E_lit = E_in[lit]
    for roll in ROLLS_DEG:
        lux_eq, landed = np.zeros(ALPHA_DEG.size), 0.0
        for sign in (1.0, -1.0):
            rr = np.radians(roll) * sign
            m = nl * np.cos(rr) + np.cross(tl, nl) * np.sin(rr)
            mu_r = np.clip(m[:, 0], 0.0, None)
            Lr = RHO_D * E_lit * mu_r / np.pi
            for i, al in enumerate(np.radians(ALPHA_DEG)):
                up = np.array([np.cos(al), np.sin(al), 0.0])
                D = Pl - moon_km * up
                d = np.linalg.norm(D, axis=-1)
                mu_o, mu_e = D @ up / d, np.einsum('ij,ij->i', D, m) / d
                ok = (mu_o > 0) & (mu_e > 0) & (mu_r > 0)
                lux_eq[i] += 0.5 * float((Lr * see_through * mu_e * mu_o * dAl / (d * 1e3) ** 2)[ok].sum())
            face = -m
            ray = xhat - 2.0 * np.einsum('j,ij->i', xhat, face)[:, None] * face
            t_close = -np.einsum('ij,ij->i', Pl, ray)
            closest = np.linalg.norm(Pl + np.maximum(t_close, 0.0)[:, None] * ray, axis=-1)
            lands = (t_close > 0) & (closest < moon_km) & (mu_r > 0)
            landed += 0.5 * float((flux_spec * (mu_r / np.maximum(mu_i[lit], 1e-12)))[lands].sum())
        rolls[f'{roll:g}'] = dict(roll_deg=roll, diffuse_lux_at_equator=np.round(lux_eq, 4).tolist(),
                                  specular_mean_lux_over_night_hemisphere=round(landed / hemisphere_m2, 5))
    sun_solid_angle = np.pi * (SUN_RADIUS / AU) ** 2
    distance_km = float(np.median(np.linalg.norm(Pl, axis=-1)))
    tile_solid_angle = (TILE_WIDTH_KM / distance_km) ** 2
    glint = {name: round(o['R'] * SUN_LUX * tile_solid_angle / sun_solid_angle, 2) for name, o in optics.items()}
    mid = int(np.argmin(np.abs(BETA_DEG)))
    return dict(rings_modelled=int(heights.size), ring_weight=round(float(weight), 3), optics=optics, rho_d=RHO_D,
                night_half_luminous_flux_lm=flux_in, night_half_tile_area_km2=float(dAl.sum() / 1e6),
                diffuse_reflected_lm=float(RHO_D * flux_in), specular_reflected_lm=float(flux_spec.sum()),
                alpha_deg=ALPHA_DEG.tolist(), beta_deg=BETA_DEG.tolist(), diffuse_lux=np.round(diffuse, 3).tolist(),
                diffuse_lux_at_equator=np.round(diffuse[:, mid], 3).tolist(),
                zenith_band_luminance_cd_m2=round(zenith_luminance, 2), specular=specular, night_roll=rolls,
                glint=dict(spot_lux=glint, spot_diameter_km=round(distance_km * 2 * SUN_RADIUS / AU, 0),
                           note='one tile\'s mirror image of the Sun, where it lands'))


def main(argv=None) -> int:
    result = model()
    product = dict(schema=SCHEMA, producer=dict(domain='illumination', files={'fleet_light/model.py': digest(__file__)},
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (LAYOUT, FILMS)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, **result)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    eq = result['diffuse_lux_at_equator']
    print(f"night half: {result['night_half_tile_area_km2']:.3g} km2 lit, {result['night_half_luminous_flux_lm']:.3g} lm in")
    print(f"diffuse at rho_d {RHO_D:g}, equator: midnight {eq[0]} lux, 45 deg {eq[9]} lux, 80 deg {eq[16]} lux; "
          f"band luminance up to {result['zenith_band_luminance_cd_m2']} cd/m2")
    for k, s in result['specular'].items():
        print(f"specular {k:16s}: lands {s['landing_share']:.2e} of {result['specular_reflected_lm']:.3g} lm, "
              f"mean {s['mean_lux_over_night_hemisphere']} lux over the night hemisphere")
    print('glint spot lux:', result['glint'])
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
