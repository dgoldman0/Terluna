"""The ring fleet's own infrared glow on the Open Moon, and the ways to answer it: a screen.

    OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.solar_shield_array.glow   # results/glow.json

The ring fleet's films absorb a few percent of the sunlight they pass and radiate it from both faces, and the faces
that look at the Moon send it part of that heat. The investigated climate is the planned one, and it leaves the glow
out, so the glow has to be answered (research/decisions.md). This screen builds the fleet's stack from the shield
study's products and follows its heat onto the Moon:

- optics: each film's reflectance and transmittance from both faces at every angle, from the protection domain's
  optical constants and its layer sequence, with the stable thin-film recursion of protection/model.py taken to
  oblique incidence (s and p averaged); the sunlight they absorb by band and at every incidence, the solar infrared
  beyond 2.5 um included;
- the stack: an N-layer model of parallel, specular, partly transparent plates in the thermal infrared, solved for
  every plate's temperature, with the radiance each stack sends along its normal toward the Moon and its normal
  reflectance; the provisioning screen's two grey opaque layers are its check;
- the fleet: every ring of the stack's layout on its own radius and forced eccentricity, tiles spaced evenly in time
  and held as the dynamics hold them (perpendicular to the velocity, with the 1.5 degree shingle), the sunlight each
  tile receives through the layers in front of it and the Moon's shadow, the layers over it in the Moon's view, and
  where its glow lands, by lunar day and night;
- what the stack returns: the films' specular reflection of the Moon's own outgoing infrared and of its reflected
  sunlight, which a tile facing the Moon images back onto the Moon, traced over each tile's view of the Moon;
- the answers: films that absorb less, heat mirrors and strong emitters, attitude, their combinations, and the added
  dimming that answers what remains, with the power that dimming gives back.

It runs no dynamics and no climate model. The Moon's outgoing infrared takes its spectrum from the atmosphere domain's
correlated-k model on the middle atmosphere's titania-stack columns and its size from the GCM design case.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection import model as thin
from protection.spectra import annulus_film as af
from protection.spectra import short_wave as sw

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results'/'glow.json'
SCHEMA = 'terluna.research.fleet-glow/1'
HEAT = HERE/'results'/'array_heat.json'
LAYOUT = HERE/'results'/'ring_layout.json'
LAYOUT_ROLLED = HERE/'results'/'ring_layout_rolled.json'
LEDGER = HERE/'results'/'integrated_ledger.json'
ATTITUDE = HERE/'results'/'attitude_schemes.json'
FILMS = ROOT/'protection'/'spectra'/'annulus_film.json'
FIRST = ROOT/'provisioning'/'results'/'first_screen.json'
NIGHT = ROOT/'illumination'/'fleet_light'/'results'/'night_light.json'
SOURCES = HERE/'glow_sources.json'
SKY = ROOT/'illumination'/'sky'/'results'/'solved_sky.json'
MIDDLE = ROOT/'atmosphere'/'middle_atmosphere'/'results'
COLUMNS = ('moon_1.2atm_titania_stack', 'moon_1.2atm_titania_stack_ts298')
PRODUCTS = (HEAT, LAYOUT, LEDGER, ATTITUDE, FILMS, FIRST, NIGHT, SOURCES, LAYOUT_ROLLED)
OPTICAL = tuple(sw.SOURCES/n for n in af.INPUTS)

DESIGN = dict(
    # The films as the comparison carries them (array_heat): titania and silica in micrometres, coated.
    films=dict(annulus=[.1, 4.], window=[1., 10.]),
    # Optics grids. The thermal tables reach 2-125 um (Franta silica, Siefke titania); beyond, the last values hold.
    thermal_um=[2., 125.], thermal_points=4000, thermal_bins=60, far_um=1000., mu_nodes=8,
    sun_cosines=[1., .95, .9, .8, .7, .6, .5, .4, .3, .2, .15, .1, .06, .03],
    solar_tail_um=[2.5, 125.], solar_tail_points=1500, sun_temperature_K=5772.,
    reference_K=200.,
    # The fleet's geometry and attitude (sourced): the shingle tilt of ring_bundle_run, the service arc, roll and ramp
    # of attitude_schemes, and its pitch envelope.
    shingle_deg=1.5, service_half_angle_deg=25., roll_ramp_deg=30., roll_deg=80., pitch_envelope_deg=[-3., 6.],
    # The tilt through the nodes, as attitude_schemes sets it: in full within 10 degrees of each node, out by 20.
    node_window_deg=[10., 20.],
    # Sampling (numerical choices).
    samples_per_ring=360, sun_samples=1440, sun_bin_km=[100., 50.], return_ring_stride=8, return_sample_stride=2,
    cap_nodes=[14, 28], landing_table=181,
    # The Moon's light. The GCM design case's planetary albedo and surface temperature (climate/gcm, corrected design
    # case, years 20-29); its outgoing infrared is uniform over the Moon and Lambertian, its reflected sunlight
    # Lambertian from the day side with a spectrally flat albedo (design guesses).
    planetary_albedo=.280, surface_K=294.9, lumens_per_watt=128000./1361.,
    # Answers: films that absorb less.
    thin_titania_um=.05,                 # the lightest coated titania that meets 1e-5 (annulus_film.json)
    wide_gap_edge_nm=[220., 250.],       # design guess: an absorber with HfO2's or Si3N4's band edge in the titania's place
    wide_gap_g_m2=[1., 2.],              # design guess: 0.1-0.2 um of a dense wide-gap oxide or nitride
    uv_mirror_nm=dict(annulus=[245., 390.], window=[245., 340.]),   # bands a dielectric mirror reflects; the window
                                         # keeps the titania's own edge above 340 nm, and with it its spectrum
    uv_mirror_reflectance=[.95, .98],    # measured anchors: torchio_2002 (99% near 250 nm), optosigma (up to 98%, 245-400)
    uv_mirror_index=2.1,                 # HfO2 near 300 nm, the high-index layers (design guess)
    uv_mirror_density=[9.68, 2.2],       # HfO2 and silica, g/cm3
    # Answers: heat mirrors on Moon-facing faces and strong emitters on the outer faces (design guesses anchored to
    # measured coatings in glow_sources.json). solar: share of sunlight the coating absorbs; blocked: share it takes
    # out of the beam, absorbed share included.
    coatings=dict(
        mesh=dict(emissivity=[.1, .3], solar=[.005, .02], blocked=[.05, .2], g_m2=[.02, .2]),
        metasurface=dict(emissivity=[.15, .4], solar=[.005, .02], blocked=[.05, .2], g_m2=[.02, .2]),
        tco=dict(emissivity=[.11, .21], solar=[.03, .1], blocked=[.05, .15], g_m2=[.5, 1.5]),
        silver=dict(emissivity=[.02, .1], solar=[.04, .1], blocked=[.15, .3], g_m2=[.1, .3]),
        emitter=dict(emissivity=[.85, .85], solar=[.002, .01], blocked=[0., 0.], g_m2=[2., 10.])),
    # The climate's response and the power the added dimming gives back. The efficacy is the GCM's: in the held glow
    # branch of the gases, glow and dimming study (13 October, results/gases_glow_dimming.json, glow.held.efficacy)
    # the dimming that holds the design temperature is 0.71 of the glow's infrared watts matched as sunlight the Moon
    # absorbs, since the infrared heats the upper air. The register carries it (13 October). On the corrected
    # radiation of 17 October the pair gave about 0.59; its rerun on the final radiation waits on the target climate.
    infrared_efficacy=.7106, conversion=[.2, .3],
)

EVIDENCE = ('A screen. The films\' optics are computed from the protection domain\'s optical constants and layer '
            'sequence with the thin-film recursion taken to oblique incidence; the stack\'s thermal state from an '
            'N-layer model of parallel specular plates; the fleet from the shield study\'s layout, with tiles held '
            'as the dynamics hold them and sunlight attenuated by projected layer counts; the Moon\'s outgoing '
            'infrared from the atmosphere domain\'s correlated-k model on its titania-stack columns scaled to the GCM '
            'design case. No dynamics, attitude law, coating or climate model is run; coatings and the Moon\'s '
            'angular emission are design guesses labelled in assumptions.')
READING_RULE = ('Heat on the Moon in W/m2 averaged over the whole Moon (global), its day hemisphere (day) or its '
                'night hemisphere (night); PW of sunlight absorbed in the films; temperatures K at the sunward '
                'crossing for the ring in the stack\'s middle. The glow is the films\' own emission; the returns '
                'are the Moon\'s own outgoing infrared and reflected sunlight that the films reflect back onto it. '
                'Added dimming in percent of the design sunlight (1,173.25 W/m2 at the top of the air, 293.3 W/m2 '
                'averaged over the Moon) at the stated infrared efficacy; power from that dimming in TW of '
                'electricity at the ledger\'s 20-30% conversion. Ranges pair low and high design guesses.')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------- optics

def stack_rt(w_um, layers, thick_um, cos0):
    """Reflectance and transmittance of a layer stack between vacuum half-spaces at incidence cos0, s and p averaged.

    protection/model.py's stable recursion (reflection coefficients built from the last interface back), with each
    layer's tilted admittance: n cos(theta) for s, n / cos(theta) for p, the cosine taken on the decaying branch."""
    w = np.asarray(w_um, float)
    s0 = math.sqrt(max(0., 1.-cos0**2))
    ns = [np.ones_like(w, complex)]+[np.asarray(n, complex)*np.ones_like(w) for n in layers]+[np.ones_like(w, complex)]
    cs = []
    for n in ns:
        c = np.sqrt(1.-(s0/n)**2+0j)
        cs.append(np.where((n*c).imag < 0, -c, c))
    R = np.zeros_like(w)
    T = np.zeros_like(w)
    for pol in 's', 'p':
        def fresnel(i, j):
            a, b = (ns[i]*cs[i], ns[j]*cs[j]) if pol == 's' else (ns[j]*cs[i], ns[i]*cs[j])
            return (a-b)/(a+b), 2*ns[i]*cs[i]/(a+b)
        r, t = fresnel(len(ns)-2, len(ns)-1)
        for j in range(len(layers)-1, -1, -1):
            phase = np.exp(2j*np.pi*ns[j+1]*cs[j+1]*thick_um[j]/w)
            rj, tj = fresnel(j, j+1)
            den = 1+rj*r*phase**2
            t, r = tj*t*phase/den, (rj+r*phase**2)/den
        R += .5*np.abs(r)**2
        T += .5*np.abs(t)**2
    return R, T


def film_stack(silica, titania, titania_um, silica_um, absorber=None):
    """The coated film's layers and thicknesses from its sunward face, as annulus_film builds them; absorber replaces
    the titania's index (for a wide-gap absorber) when given."""
    ar = sw.stored_ar_layers()
    core = titania if absorber is None else absorber
    porous, mixed = thin.mix_n(np.ones_like(silica), silica, .5), thin.mix_n(silica, core, .5)
    return [porous, silica, mixed, core, mixed, silica, porous], [ar[0], ar[1], ar[2], titania_um, ar[3], silica_um,
                                                                  ar[4]]


def planck_um(w_um, T):
    """Spectral radiance (W m^-2 sr^-1 um^-1)."""
    w = np.asarray(w_um, float)*1e-6
    return 2*K.PLANCK*K.SPEED_OF_LIGHT**2/w**5/np.expm1(K.PLANCK*K.SPEED_OF_LIGHT/(w*K.BOLTZMANN*T))*1e-6


def solar_spectrum():
    """Sunlight (W m^-2 um^-1) on two grids: annulus_film's 0.1-2,500 nm (WHI 2008 below 202 nm, TSIS-1 above) and
    the solar infrared beyond 2.5 um, a 5,772 K blackbody's shape scaled so that the whole spectrum gives the solar
    constant."""
    w_nm = af.grid()
    f = af.spectrum(w_nm)*1e3
    near = float(np.trapezoid(f, w_nm/1e3))
    tail_um = np.geomspace(*DESIGN['solar_tail_um'], DESIGN['solar_tail_points'])
    shape = planck_um(tail_um, DESIGN['sun_temperature_K'])
    beyond = float(np.trapezoid(shape, tail_um))
    # The blackbody beyond 125 um is added to the last point's weight through the scaling, a few parts in 10^5.
    scale = (K.SOLAR_CONSTANT-near)/beyond
    return dict(near_um=w_nm/1e3, near=f, near_W_m2=near, tail_um=tail_um, tail=shape*scale,
                tail_W_m2=K.SOLAR_CONSTANT-near)


def film_solar(spectrum, titania_um, silica_um, absorber_edge_nm=None, uv_mirror=None):
    """Absorptance, reflectance and transmittance of a film over the whole solar spectrum, lit on its outer
    (sunward) face or its inner (Moon-facing) face, at each incidence cosine of DESIGN; and the absorbed sunlight by
    band at normal incidence on the outer face.

    absorber_edge_nm: a wide-gap absorber in the titania's place, taken as a step (design guess): everything below
    the edge that reaches it is absorbed, nothing above; the film's other layers are kept. uv_mirror: (reflectance,
    lo, hi, faces) of an all-dielectric mirror on the named faces ahead of the titania, reflecting that share of the
    band before the film acts on the rest, lossless (design guess)."""
    w1, f1, w2, f2 = spectrum['near_um'], spectrum['near'], spectrum['tail_um'], spectrum['tail']
    s1, t1 = af.indices(w1*1e3)
    s2, t2 = sw.indices(w2*1e3)
    out = dict(cosines=list(DESIGN['sun_cosines']))
    rows = {}
    for face, order in (('outer', 1), ('inner', -1)):
        A, R, T = [], [], []
        for c in DESIGN['sun_cosines']:
            parts = []
            for w, f, si, ti in ((w1, f1, s1, t1), (w2, f2, s2, t2)):
                if absorber_edge_nm is None:
                    layers, thick = film_stack(si, ti, titania_um, silica_um)
                    r, t = stack_rt(w, layers[::order], thick[::order], c)
                    a = 1-r-t
                else:
                    # The silica film with the stack's layers and no titania, then the step absorber on what passes.
                    layers, thick = film_stack(si, si, 1e-6, silica_um)
                    r, t = stack_rt(w, layers[::order], thick[::order], c)
                    below = w*1e3 < absorber_edge_nm
                    a = np.where(below, 1-r, 1-r-t)
                    t = np.where(below, 0., t)
                if uv_mirror is not None and face in uv_mirror[3]:
                    rm, lo, hi = uv_mirror[:3]
                    band = (w*1e3 >= lo) & (w*1e3 < hi)
                    r = np.where(band, rm+(1-rm)*r, r)
                    t = np.where(band, (1-rm)*t, t)
                    a = np.where(band, (1-rm)*a, a)
                parts.append((float(np.trapezoid(a*f, w)), float(np.trapezoid(r*f, w)), float(np.trapezoid(t*f, w))))
            A.append(sum(p[0] for p in parts)/K.SOLAR_CONSTANT)
            R.append(sum(p[1] for p in parts)/K.SOLAR_CONSTANT)
            T.append(sum(p[2] for p in parts)/K.SOLAR_CONSTANT)
            if c == 1. and face == 'outer':
                rows['near_A'] = parts[0][0]/spectrum['near_W_m2']
                rows['near_R'] = parts[0][1]/spectrum['near_W_m2']
                rows['near_T'] = parts[0][2]/spectrum['near_W_m2']
                rows['tail_A'] = parts[1][0]/spectrum['tail_W_m2']
        out[face] = dict(A=A, R=R, T=T)
    # Absorbed sunlight by band, normal incidence on the outer face, in points of the solar constant.
    layers, thick = (film_stack(s1, t1, titania_um, silica_um) if absorber_edge_nm is None
                     else film_stack(s1, s1, 1e-6, silica_um))
    r, t = stack_rt(w1, layers, thick, 1.)
    a = 1-r-t
    if absorber_edge_nm is not None:
        a = np.where(w1*1e3 < absorber_edge_nm, 1-r, a)
    if uv_mirror is not None and 'outer' in uv_mirror[3]:
        rm, lo, hi = uv_mirror[:3]
        a = np.where((w1*1e3 >= lo) & (w1*1e3 < hi), (1-rm)*a, a)
    bands = {}
    for lo, hi in ((0, 200), (200, 242), (242, 320), (320, 340), (340, 360), (360, 380), (380, 400), (400, 700),
                   (700, 2500)):
        sel = (w1*1e3 >= lo) & (w1*1e3 < hi)
        bands[f'{lo}_{hi}_nm'] = float(np.trapezoid((a*f1)[sel], w1[sel])/K.SOLAR_CONSTANT*100)
    bands['beyond_2500_nm'] = rows['tail_A']*spectrum['tail_W_m2']/K.SOLAR_CONSTANT*100
    out.update(rows, bands_points=bands)
    return out


def film_visible(spectrum, titania_um, silica_um):
    """Each face's normal reflectance and the transmittance for the Moon's reflected sunlight: weighted by the
    visible band (for illuminance) and by the sunlight the window stack passes to the Moon (for its energy)."""
    w1, f1 = spectrum['near_um'], spectrum['near']
    s1, t1 = af.indices(w1*1e3)
    layers, thick = film_stack(s1, t1, titania_um, silica_um)
    ri, t = stack_rt(w1, layers[::-1], thick[::-1], 1.)
    ro, _ = stack_rt(w1, layers, thick, 1.)
    rw, tw = stack_rt(w1, *film_stack(s1, t1, *DESIGN['films']['window']), 1.)
    vis = (w1 >= .4) & (w1 < .7)
    passed = f1*tw
    def mean(x, weight, sel=slice(None)):
        return float(np.trapezoid((x*weight)[sel], w1[sel])/np.trapezoid(weight[sel], w1[sel]))
    return dict(R_visible=mean(ri, f1, vis), R_visible_outer=mean(ro, f1, vis), T_visible=mean(t, f1, vis),
                R_moonlight=mean(ri, passed), R_moonlight_outer=mean(ro, passed), T_moonlight=mean(t, passed))


def mu_nodes():
    """Gauss-Legendre direction cosines and weights on (0, 1) for hemispheric fluxes, and the normal (weight 0) for
    the Moon's direction."""
    x, w = np.polynomial.legendre.leggauss(DESIGN['mu_nodes'])
    return np.r_[.5*(x+1), 1.], np.r_[.5*w, 0.]


def thermal_bins():
    """Edges (um) of the infrared bands: 0, then log-spaced from 2 to 125 um, the far tail to 1,000 um, and infinity."""
    inner = np.geomspace(*DESIGN['thermal_um'], DESIGN['thermal_bins']+1)
    return np.r_[0., inner, DESIGN['far_um'], np.inf]


def film_thermal(titania_um, silica_um):
    """Each face's reflectance and the transmittance in the thermal infrared, averaged over each band (weighted by a
    200 K blackbody within it) at each direction cosine; the first and last bands take the end values."""
    w = np.geomspace(*DESIGN['thermal_um'], DESIGN['thermal_points'])
    silica, titania = sw.indices(w*1e3)
    layers, thick = film_stack(silica, titania, titania_um, silica_um)
    mu, _ = mu_nodes()
    edges = thermal_bins()
    nb = edges.size-1
    weight = planck_um(w, DESIGN['reference_K'])
    which = np.clip(np.searchsorted(edges, w, side='right')-1, 1, nb-3)
    out = {}
    for face, order in (('out', 1), ('in', -1)):
        R = np.zeros((nb, mu.size))
        T = np.zeros((nb, mu.size))
        for j, c in enumerate(mu):
            r, t = stack_rt(w, layers[::order], thick[::order], float(c))
            for b in range(1, nb-2):
                sel = which == b
                R[b, j] = np.sum(r[sel]*weight[sel])/np.sum(weight[sel])
                T[b, j] = np.sum(t[sel]*weight[sel])/np.sum(weight[sel])
            R[0, j], T[0, j] = r[0], t[0]
            R[-2:, j], T[-2:, j] = r[-1], t[-1]
        out[face] = (R, T)
    # Reciprocity: the transmittance is the same from both faces; average the two computations.
    tau = .5*(out['out'][1]+out['in'][1])
    return dict(rho_out=out['out'][0], rho_in=out['in'][0], tau=tau, edges_um=edges,
                tau_asymmetry=float(np.max(np.abs(out['out'][1]-out['in'][1]))))


def normal_emissivity(titania_um, silica_um, T):
    """Both faces' normal emissivity, Planck-weighted over 2-125 um at T: array_heat's quantity, for the check."""
    w = np.geomspace(*DESIGN['thermal_um'], DESIGN['thermal_points'])
    silica, titania = sw.indices(w*1e3)
    layers, thick = film_stack(silica, titania, titania_um, silica_um)
    weight = planck_um(w, T)
    faces = []
    for order in (1, -1):
        r, t = stack_rt(w, layers[::order], thick[::order], 1.)
        faces.append(float(np.trapezoid((1-r-t)*weight, w)/np.trapezoid(weight, w)))
    return faces


# ------------------------------------------------------------------------------------------------------- the stack

def blackbody_fraction(u):
    """Share of a blackbody's emission at wavelengths above lambda, u = hc/(lambda k T), from two series."""
    u = np.asarray(u, float)
    out = np.zeros_like(u)
    small = u < .5
    us = u[small]
    # Below u = 0.5 the exponential series converges slowly; the integral's power series (Bernoulli numbers) is used.
    s = us**3/3-us**4/8+us**5/60-us**7/5040+us**9/272160-us**11/13305600+us**13/622702080
    out[small] = 15/np.pi**4*s
    ul = u[~small]
    total = np.zeros_like(ul)
    for n in range(1, 90):
        total += np.exp(-n*ul)*(ul**3/n+3*ul**2/n**2+6*ul/n**3+6/n**4)
    out[~small] = 1-15/np.pi**4*total
    return out


def band_radiance(T, edges_um):
    """Blackbody radiance (W m^-2 sr^-1) in each band between consecutive edges, for temperatures T (any shape)."""
    T = np.asarray(T, float)[..., None]
    lam = np.asarray(edges_um, float)*1e-6
    c2 = K.PLANCK*K.SPEED_OF_LIGHT/K.BOLTZMANN
    with np.errstate(divide='ignore', over='ignore'):
        u = np.where(lam > 0, c2/(np.where(lam > 0, lam, 1.)*T), np.inf)
    # fraction of emission at wavelengths below each edge: 1 - F(above)
    below = np.where(np.isinf(u), 0., 1.-blackbody_fraction(np.where(np.isinf(u), 1., u)))
    below = np.where(lam == 0, 0., below)
    below = np.where(np.isinf(lam), 1., below)
    return K.STEFAN_BOLTZMANN*T**4/np.pi*np.diff(below, axis=-1)


class Plates:
    """N parallel plates, inner (nearest the Moon) first: reflectance of the outer face, of the inner face, and the
    transmittance, each (N, bands, directions)."""

    def __init__(self, rho_out, rho_in, tau, edges_um):
        self.rho_out, self.rho_in, self.tau = (np.asarray(x, float) for x in (rho_out, rho_in, tau))
        self.edges_um = np.asarray(edges_um, float)
        self.n = self.rho_out.shape[0]

    @classmethod
    def repeat(cls, plate, n):
        return cls(*(np.repeat(np.asarray(plate[k], float)[None], n, axis=0) for k in ('rho_out', 'rho_in', 'tau')),
                   plate['edges_um'])

    @classmethod
    def stack(cls, plates):
        return cls(*(np.stack([np.asarray(p[k], float) for p in plates]) for k in ('rho_out', 'rho_in', 'tau')),
                   plates[0]['edges_um'])


def response(plates):
    """Linear response of a stack: per band, the net power each plate absorbs per unit radiance source in each plate
    (G), the radiance it sends along the normal toward the Moon per unit source (D), its normal reflectance for
    radiance arriving from the Moon's side (R0), and the outward and inward fluxes per unit source (Fout, Fin).

    Unknowns per band and direction: the radiance leaving each plate's outer face upward (J_up) and its inner face
    downward (J_dn). A plate passes what reaches it through tau, reflects it from either face and emits from each
    face its absorptance times the source: J_up,k = eo_k B_k + rho_out,k J_dn,k+1 + tau_k J_up,k-1 and J_dn,k =
    ei_k B_k + rho_in,k J_up,k-1 + tau_k J_dn,k+1, with nothing arriving from space above and radiance L below."""
    n = plates.n
    mu, wmu = mu_nodes()
    ro, ri, ta = plates.rho_out, plates.rho_in, plates.tau          # (n, b, m)
    eo, ei = 1-ro-ta, 1-ri-ta
    nb, nm = ro.shape[1], ro.shape[2]
    A = np.zeros((nb, nm, 2*n, 2*n))
    rhs = np.zeros((nb, nm, 2*n, n+1))
    idx = np.arange(n)
    A[:, :, 2*idx, 2*idx] = 1.
    A[:, :, 2*idx+1, 2*idx+1] = 1.
    for k in range(n):
        if k+1 < n:
            A[:, :, 2*k, 2*(k+1)+1] = -ro[k]
            A[:, :, 2*k+1, 2*(k+1)+1] = -ta[k]
        if k >= 1:
            A[:, :, 2*k, 2*(k-1)] = -ta[k]
            A[:, :, 2*k+1, 2*(k-1)] = -ri[k]
        rhs[:, :, 2*k, k] = eo[k]
        rhs[:, :, 2*k+1, k] = ei[k]
    rhs[:, :, 0, n] = ta[0]
    rhs[:, :, 1, n] = ri[0]
    X = np.linalg.solve(A, rhs)                                      # (b, m, 2n, n+1)
    up, dn = X[:, :, 0::2, :], X[:, :, 1::2, :]                      # (b, m, n, n+1)
    below = np.zeros((nb, nm, 1, n+1))
    below[..., 0, n] = 1.
    above = np.zeros((nb, nm, 1, n+1))
    incoming = np.concatenate([below, up[:, :, :-1]], axis=2)+np.concatenate([dn[:, :, 1:], above], axis=2)
    net = incoming-(up+dn)                                           # (b, m, n, n+1) radiance per unit source
    flux = 2*np.pi*(wmu*mu)[None, :, None, None]
    G = np.sum(net*flux, axis=1)                                     # (b, n, n+1)
    return dict(G=G, D=dn[:, -1, 0, :], Fout=np.sum(up[:, :, -1, :]*flux[:, :, 0], axis=1),
                Fin=np.sum(dn[:, :, 0, :]*flux[:, :, 0], axis=1), Lin=dn[:, :, 0, :], mu=mu)


def solve(plates, q, resp=None, tol=1e-11):
    """Temperatures of the plates for absorbed sunlight q (W/m^2 per plate, inner first), and what leaves: the
    radiance along the normal toward the Moon as a flux-equivalent (pi L), the inward and outward fluxes, the
    Moon-directed radiance at every direction cosine, and the energy residual."""
    resp = resp or response(plates)
    q = np.asarray(q, float)
    edges = plates.edges_um
    G = resp['G'][:, :, :plates.n]                                   # (b, k, m)
    T = (np.maximum(q, 1e-3)/(2*.5*K.STEFAN_BOLTZMANN))**.25
    for _ in range(60):
        B = band_radiance(T, edges)                                  # (m, b)
        F = q+np.einsum('bkm,mb->k', G, B)
        dB = (band_radiance(T+1e-3, edges)-band_radiance(T-1e-3, edges))/2e-3
        J = np.einsum('bkm,mb->km', G, dB)
        step = np.linalg.solve(J, -F)
        T = np.maximum(T+step, 1.)
        if np.max(np.abs(F)) < tol*max(1., float(np.sum(q))):
            break
    B = band_radiance(T, edges)
    toward = float(np.pi*np.einsum('bm,mb->', resp['D'][:, :plates.n], B))
    inward = float(np.einsum('bm,mb->', resp['Fin'][:, :plates.n], B))
    outward = float(np.einsum('bm,mb->', resp['Fout'][:, :plates.n], B))
    by_mu = np.pi*np.einsum('bjm,mb->j', resp['Lin'][:, :, :plates.n], B)
    residual = float(np.sum(q)-inward-outward)
    return dict(T=T, toward=toward, inward=inward, outward=outward, by_mu=by_mu, residual=residual,
                F=float(np.max(np.abs(q+np.einsum('bkm,mb->k', G, B)))))


def grey_plate(e_out, e_in):
    """An opaque grey plate, angle-independent, in one band spanning every wavelength."""
    mu, _ = mu_nodes()
    one = np.ones((1, mu.size))
    return dict(rho_out=(1-e_out)*one, rho_in=(1-e_in)*one, tau=0*one, edges_um=np.array([0., np.inf]))


def coated(plate, inner=None, outer=None):
    """A film plate with a grey opaque coating on its inner face (a heat mirror of emissivity inner) and/or its outer
    face (an emitter of emissivity outer). A coated plate is opaque; the uncoated face sees the coating through the
    film, its reflectance summed incoherently: R + T^2 rho_c / (1 - R' rho_c)."""
    ro, ri, ta = plate['rho_out'], plate['rho_in'], plate['tau']
    if inner is not None and outer is not None:
        new_ro, new_ri = np.full_like(ro, 1-outer), np.full_like(ri, 1-inner)
    elif inner is not None:
        rc = 1-inner
        new_ro, new_ri = ro+ta**2*rc/(1-ri*rc), np.full_like(ri, rc)
    elif outer is not None:
        rc = 1-outer
        new_ro, new_ri = np.full_like(ro, rc), ri+ta**2*rc/(1-ro*rc)
    else:
        return plate
    return dict(rho_out=new_ro, rho_in=new_ri, tau=np.zeros_like(ta), edges_um=plate['edges_um'])


def two_layer_check(heat):
    """The provisioning screen's two-layer estimate rebuilt with the N-layer solver in grey opaque mode: the share of
    the two layers' heat that leaves toward the Moon, as built and with heat mirrors, and the hottest layer."""
    films, window = heat['tiles']['films'], heat['fleet']['window_ring_share']
    weight = {n: (window if n == 'window' else 1-window) for n in films}
    def share(e_out, e_outer_moon, e_inner_moon, q):
        plates = Plates.stack([grey_plate(e_out, e_inner_moon), grey_plate(e_out, e_outer_moon)])
        res = solve(plates, [q, q])
        return res['inward']/(2*q), res
    out = {}
    as_is = sum(weight[n]*share(f['emissivity_sunward'], f['emissivity_shaded'], f['emissivity_shaded'],
                                f['absorbed_W_m2'])[0] for n, f in films.items())
    out['two_layers_as_is'] = round(as_is, 3)
    for e in (.03, .1):
        for label, strong in (('outer_faces_as_now', None), ('outer_faces_strong', .85)):
            row = {}
            for case, inner_of_outer in (('every_tile', e), ('nearest_layer_only', None)):
                total, hottest = 0., 0
                for n, f in films.items():
                    eo = strong or f['emissivity_sunward']
                    s, res = share(eo, inner_of_outer if inner_of_outer is not None else f['emissivity_shaded'], e,
                                   f['absorbed_W_m2'])
                    total += weight[n]*s
                    hottest = max(hottest, round(float(np.max(res['T']))))
                row[case] = dict(share=round(total, 3), cut=round(as_is/total, 1), hottest_layer_K=hottest)
            out[f'mirror_{e:g}_{label}'] = row
    return out


def first_screen(heat, layout):
    """The provisioning screen's glow on the Moon, rebuilt from the same products: the films' absorbed light, a share
    toward the Moon between a deep stack's half and one layer's emissivity split, and the Moon's view from the
    outermost and innermost ring radii."""
    films, window = heat['tiles']['films'], heat['fleet']['window_ring_share']
    one = sum((window if n == 'window' else 1-window)*f['emissivity_shaded']/(f['emissivity_sunward']
              + f['emissivity_shaded']) for n, f in films.items())
    rings = layout['flown']['matched']['rings']
    radii = [min(r['radius_km'] for r in rings), max(r['radius_km'] for r in rings)]
    view = [(K.MOON_RADIUS/1e3/r)**2 for r in radii[::-1]]
    area = 4*np.pi*K.MOON_RADIUS**2
    absorbed = heat['tiles']['absorbed_PW']
    return [a*1e15*s*v/area for a, s, v in zip(absorbed, (.5, one), view)]


# ------------------------------------------------------------------------------------------------------- the fleet

SUN = np.array([-1., 0., 0.])      # toward the Sun; the night side is +x, the rings' shared node line is y


def fleet_rings(layout, heat, ledger, drive=1.):
    """Every ring of the stack: crossing height, apolune radius (the layout's radius, at the sunward crossing),
    forced eccentricity and tilt, from the flown matched rings by height; window rings are the middle share of the
    stack that array_heat carries; the tiles are shared evenly among the rings. With an attitude scheme's eccentricity
    drive over the as-built one, the forced eccentricity scales with it (attitude_schemes' reading rule) and the
    sunward crossing keeps its radius."""
    stack = layout['stack']
    flown = sorted(layout['flown']['matched']['rings'], key=lambda r: r['height_km'])
    hs = np.array([r['height_km'] for r in flown])
    heights = np.linspace(-stack['height_top_km'], stack['height_top_km'], stack['rings'])
    r_apo = np.interp(heights, hs, [r['radius_km'] for r in flown])
    e = np.interp(heights, hs, [r['forced_eccentricity'] for r in flown])*drive
    tilt = np.arcsin(heights/r_apo)
    share = heat['fleet']['window_ring_share']
    return dict(heights=heights, r_apo=r_apo, e=e, a=r_apo/(1+e), tilt=tilt,
                spacing=np.abs(np.gradient(tilt)), window=np.abs(heights) <= share*stack['height_top_km'],
                tiles_per_ring=float(np.mean(heat['fleet']['tiles']))/stack['rings'],
                tile_km=ledger['requirements']['tile_side_m']/1e3)


def fleet_samples(rings, n, ring_index=None, sample_index=None):
    """Tiles spaced evenly in time round each ring: position (km), distance, the velocity's direction, the tile's
    outward normal as the dynamics set it (perpendicular to the velocity in the orbit plane, tilted by the shingle
    against the motion), the phase from the sunward crossing, the true anomaly and the local layer count in the
    Moon's view (along-track tiles over a strip's angular length times rings over its angular width)."""
    ri = np.arange(rings['heights'].size) if ring_index is None else np.asarray(ring_index)
    M = 2*np.pi*(np.arange(n)+.5)/n
    if sample_index is not None:
        M = M[sample_index]
    a, e, t = rings['a'][ri, None], rings['e'][ri, None], rings['tilt'][ri, None]
    E = M[None, :]+0*a
    for _ in range(30):
        E = E-(E-e*np.sin(E)-M[None, :])/(1-e*np.cos(E))
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    r = a*(1-e*np.cos(E))
    un = np.stack([np.cos(t), 0*t, -np.sin(t)], axis=-1)                       # perilune, on the night side
    y = np.array([0., 1., 0.])
    P = r[..., None]*(np.cos(nu)[..., None]*un+np.sin(nu)[..., None]*y)
    dr = r*e*np.sin(nu)/(1+e*np.cos(nu))
    dP = dr[..., None]*(np.cos(nu)[..., None]*un+np.sin(nu)[..., None]*y) \
        + r[..., None]*(-np.sin(nu)[..., None]*un+np.cos(nu)[..., None]*y)
    v = dP/np.linalg.norm(dP, axis=-1, keepdims=True)
    h = np.cross(un, y)
    h = h/np.linalg.norm(h, axis=-1, keepdims=True)
    perp = np.cross(v, np.broadcast_to(h, v.shape))
    tau = math.radians(DESIGN['shingle_deg'])
    normal = math.cos(tau)*perp-math.sin(tau)*v
    phase = np.angle(np.exp(1j*(nu-np.pi)))
    side = rings['tile_km']
    along = rings['tiles_per_ring']*side*r/(2*np.pi*a**2*np.sqrt(1-e**2))
    across = (side/r)/(rings['spacing'][ri, None]*np.maximum(np.abs(np.cos(nu)), 1e-6))
    area = rings['tiles_per_ring']/n*(side*1e3)**2+0*r                        # m2 of tile per sample
    return dict(P=P, r=r, v=v, normal=normal, phase=phase, nu=nu, layers=along*across, area=area,
                window=np.broadcast_to(rings['window'][ri, None], r.shape), e=np.broadcast_to(e, r.shape), ring=ri)


def ramp(x, a, b):
    """Zero below a, one above b, half a cosine between (attitude_schemes.ramp)."""
    z = np.clip((np.asarray(x, float)-a)/(b-a), 0., 1.)
    return .5-.5*np.cos(np.pi*z)


def roll_law(phase, law, chi_deg=None):
    """Roll about the line of tile centres (rad) at each phase from the sunward crossing. 'night': the night half
    only, ramped from the nodes; 'outside_service': attitude_schemes' rolled law, zero in the service arc and at the
    nodes, one sign on the day half and the other on the night half, with cosine ramps; 'outside_service_node_<deg>':
    the same roll, with the tilt through the nodes that node_tilt gives."""
    chi = math.radians(DESIGN['roll_deg'] if chi_deg is None else chi_deg)
    if law.startswith('outside_service_node_'):
        law = 'outside_service'
    a = np.abs(phase)
    s, w = math.radians(DESIGN['service_half_angle_deg']), math.radians(DESIGN['roll_ramp_deg'])
    night = chi*ramp(a-np.pi/2, 0., w)
    if law == 'as_built':
        return 0*a
    if law == 'night':
        return np.where(a <= np.pi/2, 0., night)
    if law == 'outside_service':
        day = -chi*ramp(a, s, s+w)*ramp(np.pi/2-a, 0., w)
        return np.where(a <= np.pi/2, day, night)
    raise ValueError(law)


def node_tilt(phase, law):
    """attitude_schemes' pitch through the nodes for an 'outside_service_node_<deg>' law (rad): the tilt toward the Sun
    about the orbit normal where the roll passes through radial-facing, in full within the window's first angle of
    each node and ramped out by its second; zero for every other law."""
    if not law.startswith('outside_service_node_'):
        return 0*phase
    theta = math.radians(float(law.rsplit('_', 1)[1]))
    lo, hi = (math.radians(x) for x in DESIGN['node_window_deg'])
    off = np.abs(np.abs(phase)-np.pi/2)
    return np.sign(np.sin(phase))*theta*(1-ramp(off, lo, hi))


def rotate(x, axis, angle):
    c, s = np.cos(angle)[..., None], np.sin(angle)[..., None]
    return x*c+np.cross(axis, x)*s+axis*np.sum(axis*x, axis=-1, keepdims=True)*(1-c)


def moon_shadow(P):
    """Share of the Sun's disc a point sees past the Moon: umbra 0, outside the penumbra 1, linear between
    (illumination/fleet_light)."""
    b = np.hypot(P[..., 1], P[..., 2])
    behind = np.maximum(P[..., 0], 0.)
    alpha = K.SUN_RADIUS/K.AU
    moon = K.MOON_RADIUS/1e3
    inner, outer = moon-behind*alpha, moon+behind*alpha
    lit = np.clip((b-inner)/np.maximum(outer-inner, 1e-9), 0., 1.)
    return np.where(P[..., 0] > 0, lit, 1.)


def projected_counts(rings, S):
    """Layers of film between each tile and the Sun, from the tiles' area projected on the plane facing the Sun, for
    the day and night halves (fleet_light's coverage): on average half the other tiles of its own half lie in front,
    and a night tile also has the whole day half in front. The projected area is binned from a denser sampling of
    every ring (sun_samples per ring, closer along each ring than the bins are wide)."""
    by, bz = DESIGN['sun_bin_km']
    ey, ez = np.arange(-24000., 24000.01, by), np.arange(-11000., 11000.01, bz)
    grids = dict(day=np.zeros((ey.size-1, ez.size-1)), night=np.zeros((ey.size-1, ez.size-1)))
    n = DESIGN['sun_samples']
    for start in range(0, rings['heights'].size, 200):
        D = fleet_samples(rings, n, ring_index=np.arange(start, min(start+200, rings['heights'].size)))
        proj = np.abs(D['normal'][..., 0])*D['area']/1e6
        night = D['P'][..., 0] > 0
        for name, mask in (('day', ~night), ('night', night)):
            hist, _, _ = np.histogram2d(D['P'][..., 1][mask], D['P'][..., 2][mask], bins=(ey, ez),
                                        weights=proj[mask])
            grids[name] += hist/(by*bz)
    iy = np.clip(np.searchsorted(ey, S['P'][..., 1])-1, 0, ey.size-2)
    iz = np.clip(np.searchsorted(ez, S['P'][..., 2])-1, 0, ez.size-2)
    cd, cn = grids['day'][iy, iz], grids['night'][iy, iz]
    night = S['P'][..., 0] > 0
    return np.where(night, cd+np.maximum(cn-1, 0)/2, np.maximum(cd-1, 0)/2)


def cap_nodes(r_km):
    """Quadrature over the Moon's surface that a point r_km from its centre sees: angle from the sub-point (Gauss in
    its cosine) and azimuth, with the surface element (km2)."""
    nz, nphi = DESIGN['cap_nodes']
    R = K.MOON_RADIUS/1e3
    lo = R/np.asarray(r_km, float)
    x, w = np.polynomial.legendre.leggauss(nz)
    cz = .5*(x+1)[..., None]*(1-lo)+lo                                         # cos of the angle from the sub-point
    wz = .5*w[..., None]*(1-lo)
    phi = (np.arange(nphi)+.5)*2*np.pi/nphi
    return cz, wz, phi, R


def landing_share(r_km=19000.):
    """Share of a Moon-facing tile's emission toward the Moon that lands on the day hemisphere, against the cosine of
    the tile's angle from the subsolar direction (a tile facing the Moon's centre, emitting as a Lambertian face)."""
    cpsi = np.linspace(-1, 1, DESIGN['landing_table'])
    cz, wz, phi, R = cap_nodes(r_km)
    cz, wz = cz[:, 0], wz[:, 0]
    sz = np.sqrt(1-cz**2)
    out = []
    for c in cpsi:
        p = np.array([-c, math.sqrt(max(0., 1-c*c)), 0.])                   # at angle psi from the Sun
        e1 = np.cross(p, [0., 0., 1.])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(p, e1)
        d = cz[:, None, None]*p+sz[:, None, None]*(np.cos(phi)[None, :, None]*e1+np.sin(phi)[None, :, None]*e2)
        D = R*d
        P = r_km*p
        to_tile = P-D
        dist = np.linalg.norm(to_tile, axis=-1)
        cos_d = np.sum(d*to_tile, axis=-1)/dist
        cos_t = np.sum(p*to_tile, axis=-1)/dist
        wgt = np.clip(cos_d, 0, None)*np.clip(cos_t, 0, None)/dist**2*(wz[:, None]*R**2*(2*np.pi/phi.size))
        day = (d @ SUN) > 0
        out.append(float(np.sum(wgt*day)/np.sum(wgt)))
    return cpsi, np.array(out)


def specular_geometry(S, normal, chunk=1500):
    """For each tile, the Moon's light it receives and where its specular reflection lands, per unit radiance: for
    light the Moon sends out evenly (its outgoing infrared) and for light leaving only the day side as a Lambertian
    reflector (its reflected sunlight). Each tile's view of the Moon is integrated over the cap it sees; a ray from
    a point of the Moon is reflected about the tile's normal and followed until it meets the Moon or misses it.

    Returns per tile: incident (W/m2 per W/m2/sr), landed on the day side, landed on the night side, for both kinds of
    light, and the reflected sunlight landed on the night side by angle from the antisolar point (10-degree bands)."""
    P = S['P'].reshape(-1, 3)
    nrm = normal.reshape(-1, 3)
    r = np.linalg.norm(P, axis=-1)
    R = K.MOON_RADIUS/1e3
    nz, nphi = DESIGN['cap_nodes']
    x, wx = np.polynomial.legendre.leggauss(nz)
    phi = (np.arange(nphi)+.5)*2*np.pi/nphi
    out = {k: np.zeros(r.size) for k in ('ir_in', 'ir_day', 'ir_night', 'sw_in', 'sw_day', 'sw_night')}
    bands = np.zeros((r.size, 9))
    for s0 in range(0, r.size, chunk):
        sl = slice(s0, min(r.size, s0+chunk))
        Pc, nc, rc = P[sl], nrm[sl], r[sl]
        p = Pc/rc[:, None]
        e1 = np.cross(p, np.array([0., 0., 1.]))
        bad = np.linalg.norm(e1, axis=-1) < 1e-6
        e1[bad] = np.cross(p[bad], np.array([1., 0., 0.]))
        e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
        e2 = np.cross(p, e1)
        lo = R/rc
        cz = .5*(x[None, :]+1)*(1-lo[:, None])+lo[:, None]                    # (t, z)
        wz = .5*wx[None, :]*(1-lo[:, None])
        sz = np.sqrt(1-cz**2)
        d = (cz[..., None, None]*p[:, None, None, :]
             + sz[..., None, None]*(np.cos(phi)[None, None, :, None]*e1[:, None, None, :]
                                    + np.sin(phi)[None, None, :, None]*e2[:, None, None, :]))    # (t, z, f, 3)
        D = R*d
        to_tile = Pc[:, None, None, :]-D
        dist = np.linalg.norm(to_tile, axis=-1)
        u = to_tile/dist[..., None]
        cos_d = np.clip(np.sum(d*u, axis=-1), 0, None)
        cos_t = np.clip(np.sum(u*nc[:, None, None, :], axis=-1), 0, None)
        domega = (R**2*wz[..., None]*(2*np.pi/nphi))*cos_d/dist**2
        flux = cos_t*domega                                                    # per unit radiance
        ur = u-2*np.sum(u*nc[:, None, None, :], axis=-1)[..., None]*nc[:, None, None, :]
        b = np.sum(Pc[:, None, None, :]*ur, axis=-1)
        disc = b**2-(rc**2-R**2)[:, None, None]
        lands = (b < 0) & (disc >= 0)
        sdist = -b-np.sqrt(np.clip(disc, 0, None))
        X = Pc[:, None, None, :]+sdist[..., None]*ur
        day = (X @ SUN) > 0
        bright = np.clip(d @ SUN, 0, None)
        for key, w in (('ir', flux), ('sw', flux*bright)):
            out[key+'_in'][sl] = np.sum(w, axis=(1, 2))
            out[key+'_day'][sl] = np.sum(w*(lands & day), axis=(1, 2))
            out[key+'_night'][sl] = np.sum(w*(lands & ~day), axis=(1, 2))
        cos_anti = np.clip(-(X @ SUN)/R, -1, 1)
        angle = np.degrees(np.arccos(cos_anti))
        which = np.clip((angle//10).astype(int), 0, 8)
        wn = flux*bright*(lands & ~day)
        for k in range(9):
            bands[sl, k] = np.sum(wn*(which == k), axis=(1, 2))
    shape = S['P'].shape[:-1]
    res = {k: v.reshape(shape) for k, v in out.items()}
    res['sw_night_bands'] = bands.reshape(shape+(9,))
    return res


# ---------------------------------------------------------------------------------- the Moon's outgoing infrared

def outgoing_infrared(top_of_air):
    """The Moon's outgoing infrared by band: the atmosphere domain's correlated-k model on the middle atmosphere's
    clear titania-stack columns at 288 and 298 K, interpolated to the design case's surface temperature, scaled to
    the GCM design case's outgoing infrared (its absorbed sunlight, (1 - albedo) S/4)."""
    from atmosphere.radiative_convective import ck
    from atmosphere.middle_atmosphere import lbl_check, run as rn
    plan = rn.cases()
    model = ck.CKLongwave()
    rows = []
    for name in COLUMNS:
        case = plan[name]
        prof, t, trop, o3, col, _ = lbl_check.rebuild(case, MIDDLE)
        f = model.fluxes(col, extra={'O3': o3*col['layer_column_cm2']})
        rows.append((float(case.surface_temperature_k), float(f['olr']), np.asarray(f['olr_bands'], float)))
    (t0, o0, b0), (t1, o1, b1) = rows
    x = (DESIGN['surface_K']-t0)/(t1-t0)
    bands = b0+(b1-b0)*x
    clear = float(np.sum(bands))
    total = (1-DESIGN['planetary_albedo'])*top_of_air/4
    return dict(bands_cm1=[list(b) for b in ck.BANDS], clear_sky_W_m2=clear, design_W_m2=total,
                columns={name: dict(surface_K=t, olr_W_m2=o) for name, (t, o, _) in zip(COLUMNS, rows)},
                band_W_m2=(bands*total/clear).tolist(), table_key=ck.table_key())



def infrared_weights(olr, edges_um):
    """The Moon's outgoing infrared on the thermal bands, each correlated-k band spread evenly in wavenumber."""
    edges = np.asarray(edges_um, float)
    with np.errstate(divide='ignore'):
        nu_edges = np.where(edges > 0, 1e4/np.where(edges > 0, edges, 1.), np.inf)        # cm-1, falling
    w = np.zeros(edges.size-1)
    for (lo, hi), p in zip(olr['bands_cm1'], olr['band_W_m2']):
        for b in range(edges.size-1):
            top, bottom = min(hi, nu_edges[b]), max(lo, nu_edges[b+1])
            if top > bottom:
                w[b] += p*(top-bottom)/(hi-lo)
    return w


# ------------------------------------------------------------------------------------------ films and their stacks

_THERMAL = {}


def thermal_plate(titania_um, silica_um):
    key = (titania_um, silica_um)
    if key not in _THERMAL:
        _THERMAL[key] = film_thermal(titania_um, silica_um)
    return _THERMAL[key]


LAYERS = list(range(1, 13))+[14, 16, 20, 24, 28, 32]


def stack_tables(plate, q_column, olr_w, moonlight, visible, first=None):
    """For stacks of N identical plates sharing q_column of absorbed sunlight evenly: the share of it sent toward the
    Moon as a flux-equivalent at each direction cosine (pi L / q), the hemispheric inward share, the plates'
    temperatures, the normal reflectance for the Moon's infrared (weighted by its spectrum) and, by incoherent adding,
    for its reflected sunlight and visible light."""
    mu, _ = mu_nodes()
    rows = []
    for n in LAYERS:
        plates = Plates.repeat(plate, n) if first is None else Plates.stack([first]+[plate]*(n-1))
        resp = response(plates)
        res = solve(plates, [q_column/n]*n, resp)
        r_ir = float(np.sum(resp['D'][:, n]*olr_w)/np.sum(olr_w))
        rows.append(dict(n=n, by_mu=(res['by_mu']/q_column).tolist(), inward=res['inward']/q_column,
                         T=res['T'].tolist(), residual=res['residual']/q_column, R_ir=r_ir))
    for key, (ri, ro, t) in (('R_sw', moonlight), ('R_vis', visible)):
        rn = ri if first is None else first.get(key, ri)
        values = [rn]
        for _ in range(1, LAYERS[-1]):
            rn = ri+t**2*rn/(1-ro*rn)
            values.append(rn)
        for row in rows:
            row[key] = values[row['n']-1]
    return dict(mu=mu.tolist(), rows=rows)


def table_lookup(table, n_eff, mu_view, key='by_mu'):
    """Interpolate a stack table in the layer count (linear) and, for by_mu, the direction cosine (linear)."""
    ns = np.array([row['n'] for row in table['rows']], float)
    n = np.clip(n_eff, 1, ns[-1])
    i = np.clip(np.searchsorted(ns, n)-1, 0, ns.size-2)
    f = (n-ns[i])/(ns[i+1]-ns[i])
    if key != 'by_mu':
        v = np.array([row[key] for row in table['rows']])
        return v[i]*(1-f)+v[i+1]*f
    grid = np.array([row['by_mu'] for row in table['rows']])                   # (n, mu)
    mu = np.array(table['mu'])
    order = np.argsort(mu)
    mu, grid = mu[order], grid[:, order]
    m = np.clip(mu_view, mu[0], mu[-1])
    j = np.clip(np.searchsorted(mu, m)-1, 0, mu.size-2)
    g = (m-mu[j])/(mu[j+1]-mu[j])
    lo = grid[i, j]*(1-g)+grid[i, j+1]*g
    hi = grid[i+1, j]*(1-g)+grid[i+1, j+1]*g
    return lo*(1-f)+hi*f


def uv_mirror_mass(lo_nm, hi_nm, n_high, reflectance):
    """Areal mass (g/m2) of an all-dielectric quarter-wave mirror over a band (design estimate): pairs of a
    high-index oxide and silica, enough for the stated reflectance at each stack's centre, and enough stacks, their
    stop bands tiled geometrically, to cover the band."""
    n_low = float(np.real(sw.indices(np.array([300.]))[0][0]))
    n_sub = n_low
    y = (1+math.sqrt(reflectance))/(1-math.sqrt(reflectance))
    pairs = math.ceil(math.log(y*n_sub/n_high**2)/(2*math.log(n_high/n_low)))
    dg = 2/math.pi*math.asin((n_high-n_low)/(n_high+n_low))
    ratio = (1+dg)/(1-dg)
    stacks = math.ceil(math.log(hi_nm/lo_nm)/math.log(ratio))
    centres = lo_nm*(hi_nm/lo_nm)**((np.arange(stacks)+.5)/stacks)
    rho_h, rho_l = DESIGN['uv_mirror_density']
    high_nm = (pairs+1)*np.sum(centres)/(4*n_high)
    low_nm = pairs*np.sum(centres)/(4*n_low)
    return dict(pairs=pairs, stacks=stacks, thickness_um=float((high_nm+low_nm)/1e3),
                g_m2=float((high_nm*rho_h+low_nm*rho_l)/1e3), silica_index_300nm=n_low)


# ------------------------------------------------------------------------------------------------------ the answers

def uv_mirror(ring, end, faces):
    lo, hi = DESIGN['uv_mirror_nm'][ring]
    return (DESIGN['uv_mirror_reflectance'][end], lo, hi, faces)


FILM_VARIANTS = dict(
    as_built=dict(annulus={}, window={}),
    thin_titania=dict(annulus=dict(titania_um=DESIGN['thin_titania_um']), window={}),
    wide_gap_220=dict(annulus=dict(absorber_edge_nm=DESIGN['wide_gap_edge_nm'][0]), window={}),
    wide_gap_250=dict(annulus=dict(absorber_edge_nm=DESIGN['wide_gap_edge_nm'][1]), window={}),
    uv_mirror_outer_95=dict(annulus=dict(uv=uv_mirror('annulus', 0, ('outer',))),
                            window=dict(uv=uv_mirror('window', 0, ('outer',)))),
    uv_mirror_outer_98=dict(annulus=dict(uv=uv_mirror('annulus', 1, ('outer',))),
                            window=dict(uv=uv_mirror('window', 1, ('outer',)))),
    uv_mirror_both_98=dict(annulus=dict(uv=uv_mirror('annulus', 1, ('outer', 'inner'))),
                           window=dict(uv=uv_mirror('window', 1, ('outer', 'inner')))),
    wide_gap_250_window_mirror=dict(annulus=dict(absorber_edge_nm=DESIGN['wide_gap_edge_nm'][1]),
                                    window=dict(uv=uv_mirror('window', 1, ('outer',)))),
)


def coating_spec(key):
    """A coating choice: heat mirror (inner face), emitter (outer face), the sunlight it absorbs and takes out of the
    beam, its mass; 'low' takes each design guess's favourable end, 'high' the other."""
    if key is None:
        return dict(inner=None, outer=None, solar=0., blocked=0., g_m2=0., nearest=False)
    parts = key.split('+')
    spec = dict(inner=None, outer=None, solar=0., blocked=0., g_m2=0., nearest=False)
    for part in parts:
        name, _, end = part.partition('_')
        nearest = end.endswith('nearest')
        end = end.replace('_nearest', '')
        c = DESIGN['coatings'][name]
        i = 0 if end in ('low', '') else 1
        if name == 'emitter':
            spec['outer'] = c['emissivity'][i]
        else:
            spec['inner'] = c['emissivity'][i]
            spec['nearest'] = nearest
        spec['solar'] += c['solar'][i]
        spec['blocked'] += c['blocked'][i]
        spec['g_m2'] += c['g_m2'][i]
    return spec


def film_bundle(spectrum, olr_w, q_col, ring, film_key, coat_key):
    """One ring type's films for an answer: sunlight absorbed by face and incidence, the thermal plate, the stack
    tables at that film's column absorption, the sunward crossing's temperatures, and the reflectances the returns
    use."""
    base = DESIGN['films'][ring]
    v = FILM_VARIANTS[film_key][ring]
    ti, si = v.get('titania_um', base[0]), base[1]
    coat = coating_spec(coat_key)
    sol = film_solar(spectrum, ti, si, v.get('absorber_edge_nm'), v.get('uv'))
    vis = film_visible(spectrum, ti if v.get('absorber_edge_nm') is None else 1e-6, si)
    th = thermal_plate(ti, si)
    plain = dict(rho_out=th['rho_out'], rho_in=th['rho_in'], tau=th['tau'], edges_um=th['edges_um'])
    plate = coated(plain, coat['inner'], coat['outer'])
    alpha = sol['outer']['A'][0]+coat['solar']
    # The coating's reflection of the Moon's light: the beam share it takes out, less what it absorbs, on the face
    # it covers (a heat mirror faces the Moon); the film's transmission falls by the share it takes out.
    keep = 1-coat['blocked']
    add = coat['blocked']-coat['solar']
    def triple(r_in, r_out, t, covered):
        return ((r_in+add, r_out, t*keep) if covered == 'inner' else (r_in, r_out+add, t*keep) if covered == 'outer'
                else (r_in, r_out, t))
    covered = 'inner' if coat['inner'] is not None else 'outer' if coat['outer'] is not None else None
    moon = triple(vis['R_moonlight'], vis['R_moonlight_outer'], vis['T_moonlight'], covered)
    seen = triple(vis['R_visible'], vis['R_visible_outer'], vis['T_visible'], covered)
    q = q_col[ring]*alpha/q_col['alpha_'+ring]
    if coat['nearest']:
        plain_moon = (vis['R_moonlight'], vis['R_moonlight_outer'], vis['T_moonlight'])
        plain_seen = (vis['R_visible'], vis['R_visible_outer'], vis['T_visible'])
        first = dict(plate, R_sw=moon[0], R_vis=seen[0])
        table = stack_tables(plain, q, olr_w, plain_moon, plain_seen, first=first)
        stack1, stack2 = Plates.stack([plate]), Plates.stack([plate, plain])
    else:
        table = stack_tables(plate, q, olr_w, moon, seen)
        stack1, stack2 = Plates.repeat(plate, 1), Plates.repeat(plate, 2)
    q_sun = K.SOLAR_CONSTANT*alpha
    T1 = solve(stack1, [q_sun])['T']
    T2 = solve(stack2, [q_sun, q_sun])['T']
    # Light pressure at normal incidence on the outer face, (2R + A) S/c: what the film's sail loading rests on. A
    # coating's beam share counts as reflected, less what it absorbs.
    pressure = 2*(sol['outer']['R'][0]+coat['blocked']-coat['solar'])+alpha
    return dict(solar=sol, alpha=alpha, blocked=coat['blocked'], solar_extra=coat['solar'], transmission=keep,
                pressure_coefficient=pressure,
                table=table, q_column=q, T_sunward_one=T1.tolist(), T_sunward_two=T2.tolist(), visible=vis,
                coating=coat, titania_um=ti, silica_um=si)


def attitude(S, law='as_built', chi_deg=None, pitch_deg=0.):
    """Each tile's normal under a law: pitched about the orbit normal (a constant pitch and the law's tilt through the
    nodes) and then rolled about the line of tile centres, as attitude_schemes turns its frames."""
    chi = roll_law(S['phase'], law, chi_deg)
    pitch = node_tilt(S['phase'], law)+math.radians(pitch_deg)
    normal = S['normal']
    if np.any(pitch):
        h = np.cross(S['P'], S['v'])
        h = h/np.linalg.norm(h, axis=-1, keepdims=True)
        normal = rotate(normal, h, np.broadcast_to(pitch, chi.shape))
    normal = rotate(normal, S['v'], chi)
    return normal, chi


def absorbed(S, fixed, bundles, normal):
    """Sunlight each tile absorbs (W): the incidence on its lit face, the film's absorptance there, the Moon's
    shadow and the layers of film in front of it."""
    c = normal @ SUN
    mu = np.abs(c)
    out = np.zeros_like(mu)
    for ring, mask in (('annulus', ~S['window']), ('window', S['window'])):
        b = bundles[ring]
        cos = np.array(b['solar']['cosines'])[::-1]
        a_out = np.interp(mu, cos, np.array(b['solar']['outer']['A'])[::-1])
        a_in = np.interp(mu, cos, np.array(b['solar']['inner']['A'])[::-1])
        t = np.interp(mu, cos, np.array(b['solar']['outer']['T'])[::-1])*b['transmission']
        alpha = np.where(c > 0, a_out, a_in)+b['solar_extra']
        q = K.SOLAR_CONSTANT*mu*alpha*fixed['shadow']*t**fixed['front']
        out = np.where(mask, q*S['area'], out)
    return out


def evaluate(S, fixed, bundles, geo, law='as_built', chi_deg=None, pitch_deg=0.):
    """The glow and the returns on the Moon for one answer: films (bundles) and attitude."""
    normal, chi = attitude(S, law, chi_deg, pitch_deg)
    p_abs = absorbed(S, fixed, bundles, normal)
    phat = S['P']/S['r'][..., None]
    mu_view = np.clip(np.sum(normal*phat, axis=-1), 0., 1.)
    n_eff = 1+(S['layers']-1)*np.abs(np.cos(chi))
    R = K.MOON_RADIUS/1e3
    glow = np.zeros_like(p_abs)
    for ring, mask in (('annulus', ~S['window']), ('window', S['window'])):
        s = table_lookup(bundles[ring]['table'], n_eff, mu_view)
        glow = np.where(mask, p_abs*s*mu_view*(R/S['r'])**2, glow)
    cpsi, share = fixed['landing']
    day = np.interp(phat @ SUN, cpsi, share)
    area = 4*np.pi*K.MOON_RADIUS**2
    half = area/2
    # Returns, on the subsample the specular geometry was traced for.
    sub = fixed['sub']
    w = fixed['sub_weight']
    olr = fixed['olr_W_m2']
    sun_top = fixed['top_of_air_W_m2']
    lux_top = fixed['top_of_air_lux']
    A = DESIGN['planetary_albedo']
    ret = dict(ir_day=0., ir_night=0., sw_day=0., sw_night=0., lux_bands=np.zeros(9))
    for ring, mask in (('annulus', ~S['window'][sub]), ('window', S['window'][sub])):
        tab = bundles[ring]['table']
        ne = n_eff[sub]
        r_ir = table_lookup(tab, ne, None, 'R_ir')/ne
        r_sw = table_lookup(tab, ne, None, 'R_sw')/ne
        r_vis = table_lookup(tab, ne, None, 'R_vis')/ne
        a_m2 = S['area'][sub]*w*mask
        ret['ir_day'] += float(np.sum(a_m2*r_ir*olr/np.pi*geo['ir_day']))
        ret['ir_night'] += float(np.sum(a_m2*r_ir*olr/np.pi*geo['ir_night']))
        ret['sw_day'] += float(np.sum(a_m2*r_sw*A*sun_top/np.pi*geo['sw_day']))
        ret['sw_night'] += float(np.sum(a_m2*r_sw*A*sun_top/np.pi*geo['sw_night']))
        ret['lux_bands'] += np.sum((a_m2*r_vis*A*lux_top/np.pi)[..., None]*geo['sw_night_bands'], axis=(0, 1))
    edges = np.radians(np.arange(0., 90.01, 10.))
    band_area = 2*np.pi*K.MOON_RADIUS**2*(np.cos(edges[:-1])-np.cos(edges[1:]))
    night_half = S['P'][..., 0] > 0
    row = dict(
        absorbed_PW=float(np.sum(p_abs))/1e15,
        absorbed_by_half_PW=dict(day=float(np.sum(p_abs*~night_half))/1e15, night=float(np.sum(p_abs*night_half))/1e15),
        absorbed_by_ring_PW=dict(annulus=float(np.sum(p_abs*~S['window']))/1e15,
                                 window=float(np.sum(p_abs*S['window']))/1e15),
        glow_W_m2={'global': float(np.sum(glow))/area, 'day': float(np.sum(glow*day))/half,
                   'night': float(np.sum(glow*(1-day)))/half},
        glow_by_half_W_m2=dict(day_half=float(np.sum(glow*~night_half))/area,
                               night_half=float(np.sum(glow*night_half))/area),
        glow_by_ring_W_m2=dict(annulus=float(np.sum(glow*~S['window']))/area,
                               window=float(np.sum(glow*S['window']))/area),
        infrared_return_W_m2={'global': (ret['ir_day']+ret['ir_night'])/area, 'day': ret['ir_day']/half,
                              'night': ret['ir_night']/half},
        sunlight_return_W_m2={'global': (ret['sw_day']+ret['sw_night'])/area, 'day': ret['sw_day']/half,
                              'night': ret['sw_night']/half},
        night_lux_from_returns=dict(mean=float(np.sum(ret['lux_bands']))/half,
                                    by_angle_from_midnight_deg=dict(edges=list(range(0, 91, 10)),
                                                                    lux=(ret['lux_bands']/band_area).tolist())))
    infrared = {k: row['glow_W_m2'][k]+row['infrared_return_W_m2'][k] for k in ('global', 'day', 'night')}
    row['infrared_total_W_m2'] = infrared
    row['heat_total_W_m2'] = {k: infrared[k]+row['sunlight_return_W_m2'][k] for k in ('global', 'day', 'night')}
    return row


# The answers screened: (key, film variant, coating, attitude law, roll degrees, pitch degrees).
ANSWERS = [
    ('as_built', 'as_built', None, 'as_built', None, 0.),
    ('thin_titania', 'thin_titania', None, 'as_built', None, 0.),
    ('wide_gap_220', 'wide_gap_220', None, 'as_built', None, 0.),
    ('wide_gap_250', 'wide_gap_250', None, 'as_built', None, 0.),
    ('uv_mirror_outer_95', 'uv_mirror_outer_95', None, 'as_built', None, 0.),
    ('uv_mirror_outer_98', 'uv_mirror_outer_98', None, 'as_built', None, 0.),
    ('uv_mirror_both_98', 'uv_mirror_both_98', None, 'as_built', None, 0.),
    ('mesh_low', 'as_built', 'mesh_low', 'as_built', None, 0.),
    ('mesh_high', 'as_built', 'mesh_high', 'as_built', None, 0.),
    ('metasurface_low', 'as_built', 'metasurface_low', 'as_built', None, 0.),
    ('metasurface_high', 'as_built', 'metasurface_high', 'as_built', None, 0.),
    ('tco_low', 'as_built', 'tco_low', 'as_built', None, 0.),
    ('tco_high', 'as_built', 'tco_high', 'as_built', None, 0.),
    ('silver_low', 'as_built', 'silver_low', 'as_built', None, 0.),
    ('silver_high', 'as_built', 'silver_high', 'as_built', None, 0.),
    ('emitter', 'as_built', 'emitter_low', 'as_built', None, 0.),
    ('mesh_low_emitter', 'as_built', 'mesh_low+emitter_low', 'as_built', None, 0.),
    ('mesh_low_nearest_emitter', 'as_built', 'mesh_low_nearest+emitter_low', 'as_built', None, 0.),
    ('pitch_plus_6', 'as_built', None, 'as_built', None, DESIGN['pitch_envelope_deg'][1]),
    ('pitch_minus_3', 'as_built', None, 'as_built', None, DESIGN['pitch_envelope_deg'][0]),
    ('night_roll_30', 'as_built', None, 'night', 30., 0.),
    ('night_roll_60', 'as_built', None, 'night', 60., 0.),
    ('night_roll_80', 'as_built', None, 'night', 80., 0.),
    ('rolled_outside_service_30', 'as_built', None, 'outside_service', 30., 0.),
    ('rolled_outside_service_60', 'as_built', None, 'outside_service', 60., 0.),
    ('rolled_outside_service_80', 'as_built', None, 'outside_service', 80., 0.),
    ('night_roll_80_thin_titania', 'thin_titania', None, 'night', 80., 0.),
    ('night_roll_80_wide_gap_220', 'wide_gap_220', None, 'night', 80., 0.),
    ('night_roll_80_wide_gap_250', 'wide_gap_250', None, 'night', 80., 0.),
    ('night_roll_80_wide_gap_250_window_mirror', 'wide_gap_250_window_mirror', None, 'night', 80., 0.),
    ('night_roll_80_uv_mirror_outer_98', 'uv_mirror_outer_98', None, 'night', 80., 0.),
    ('night_roll_80_emitter', 'as_built', 'emitter_low', 'night', 80., 0.),
    ('night_roll_80_mesh_low_emitter', 'as_built', 'mesh_low+emitter_low', 'night', 80., 0.),
    ('rolled_outside_service_80_wide_gap_220', 'wide_gap_220', None, 'outside_service', 80., 0.),
    ('rolled_outside_service_80_wide_gap_250', 'wide_gap_250', None, 'outside_service', 80., 0.),
    ('rolled_outside_service_80_mesh_low_emitter', 'as_built', 'mesh_low+emitter_low', 'outside_service', 80., 0.),
    ('rolled_outside_service_80_uv_mirror_outer_98', 'uv_mirror_outer_98', None, 'outside_service', 80., 0.),
]
# The 80-degree rolls with the forced eccentricity their own push leaves: attitude_schemes' scheme for each law, its
# drive scaling every ring's eccentricity about the sunward crossing's radius.
OWN_ECCENTRICITY = [
    ('night_roll_80_own_eccentricity', 'as_built', None, 'night', 80., 0., 'night_rolled_80_deg'),
    ('rolled_outside_service_80_own_eccentricity', 'as_built', None, 'outside_service', 80., 0., 'rolled_80_deg'),
    ('rolled_outside_service_80_node_3_own_eccentricity', 'as_built', None, 'outside_service_node_3', 80., 0.,
     'rolled_80_deg_node_3_deg'),
    ('rolled_outside_service_80_node_6_own_eccentricity', 'as_built', None, 'outside_service_node_6', 80., 0.,
     'rolled_80_deg_node_6_deg'),
]
# The rolled law on the layout re-flown for its push (ring_layout_rolled): the stack's radii and forced eccentricities
# as that run lays them out, with the films the next answers would add.
ROLLED_LAYOUT = [
    ('rolled_outside_service_80_rolled_layout', 'as_built', None, 'outside_service', 80., 0., 'layout'),
    ('rolled_outside_service_80_wide_gap_250_rolled_layout', 'wide_gap_250', None, 'outside_service', 80., 0., 'layout'),
    ('rolled_outside_service_80_uv_mirror_outer_98_rolled_layout', 'uv_mirror_outer_98', None, 'outside_service', 80.,
     0., 'layout'),
]
SCHEME_OF_LAW = {'night': 'night_rolled_{:g}_deg', 'outside_service': 'rolled_{:g}_deg'}


def scheme_of(law, chi):
    """attitude_schemes' scheme for a rolled law at chi degrees, or None."""
    if law.startswith('outside_service_node_'):
        return f"rolled_{chi:g}_deg_node_{law.rsplit('_', 1)[1]}_deg"
    return SCHEME_OF_LAW[law].format(chi) if law in SCHEME_OF_LAW else None


def dimming(row, design_sun):
    """Added dimming (% of the design sunlight) that answers the heat a row leaves: infrared the Moon absorbs whole,
    at the stated efficacy, against sunlight the Moon absorbs (1 - albedo) of; returned sunlight as sunlight. With it
    the records' convention (the heat over the design sunlight) and the electricity that dimming gives back."""
    a = DESIGN['planetary_albedo']
    ir, sw_ = row['infrared_total_W_m2']['global'], row['sunlight_return_W_m2']['global']
    d = 100*(DESIGN['infrared_efficacy']*ir/((1-a)*design_sun)+sw_/design_sun)
    return dict(percent=d, records_convention_percent=100*(ir+sw_)/design_sun,
                percent_at_efficacy_1=100*(ir/((1-a)*design_sun)+sw_/design_sun),
                glow_only_percent=100*DESIGN['infrared_efficacy']*row['glow_W_m2']['global']/((1-a)*design_sun))


SAIL = dict(window_tile_g_m2=26., window_dynamic_coefficient=.276, as_built=None)


def sail_loading(bundles, added):
    """One sail loading for the nested bundle: the annulus tile mass that carries the window tiles' loading (the
    comparison's 26 g/m2 x 0.137 / 0.276, here with the films' full-spectrum coefficients), with the window's dynamic
    coefficient (its redirected 13.8% carried as reflected) moved by its film's change and its tile by the mass an
    answer adds; against the annulus tile's own mass with what the answer adds to it, and the reflection the annulus
    film would need to carry its own mass at the window's loading."""
    built = SAIL['as_built']
    window = SAIL['window_dynamic_coefficient']+bundles['window']['pressure_coefficient']-built['window']
    window_tile = SAIL['window_tile_g_m2']+added['window']
    loading = window/window_tile
    base = SAIL['window_tile_g_m2']*built['annulus']/SAIL['window_dynamic_coefficient']
    annulus_tile = base+added['annulus']
    pa = bundles['annulus']['pressure_coefficient']
    return dict(annulus_coefficient=pa, window_dynamic_coefficient=window, window_tile_g_m2=window_tile,
                annulus_tile_g_m2=annulus_tile, matched_annulus_tile_g_m2=pa/loading,
                mismatch_g_m2=annulus_tile-pa/loading, as_built_annulus_tile_g_m2=base,
                reflection_to_match=max(0., (loading*annulus_tile-pa)/2))


def costs(key, film, coat_key, law, chi, pitch, bundles, attitude_p):
    """What an answer adds to a tile and asks of its attitude: mass (g/m2 and t a tile), sunlight absorbed and taken
    out of the beam, the films' temperatures at the sunward crossing, and the attitude scheme's push and store."""
    coat = coating_spec(coat_key)
    out = dict(films=film, coating=coat_key, attitude=law, roll_deg=chi, pitch_deg=pitch)
    mass = {}
    for ring, b in bundles.items():
        base_ti = DESIGN['films'][ring][0]
        m = coat['g_m2']
        v = FILM_VARIANTS[film][ring]
        if 'titania_um' in v:
            m += (v['titania_um']-base_ti)*af.DENSITY_G_UM['titania']
        if 'absorber_edge_nm' in v:
            m += float(np.mean(DESIGN['wide_gap_g_m2']))-base_ti*af.DENSITY_G_UM['titania']
        if 'uv' in v:
            rm, lo, hi, faces = v['uv']
            m += len(faces)*uv_mirror_mass(lo, hi, DESIGN['uv_mirror_index'], rm)['g_m2']
        mass[ring] = dict(g_m2=m, t_per_tile=m*1e-6*1e8, solar_absorptance=b['alpha'],
                          taken_out_of_beam=b['blocked'], pressure_coefficient=b['pressure_coefficient'],
                          T_sunward_one_K=b['T_sunward_one'], T_sunward_two_K=b['T_sunward_two'])
    out['by_ring'] = mass
    out['sail_loading'] = sail_loading(bundles, {ring: v['g_m2'] for ring, v in mass.items()})
    name = scheme_of(law, chi) if chi else None
    if name:
        scheme = attitude_p['schemes'].get(name)
        if scheme:
            centred = attitude_p['schemes'].get('balanced_'+name)
            out['attitude_study'] = dict(
                scheme=name, law=('night half only' if law == 'night' else 'outside service, both halves'
                                  + (', tilted through the nodes' if law.startswith('outside_service_node_') else '')),
                eccentricity_drive_over_as_built=scheme['eccentricity_drive']['over_as_built'],
                store_N_m_s_median_trim_25=scheme['store_N_m_s']['trim_0.25']['median'],
                store_N_m_s_median_trim_25_centred_filter=(centred['store_N_m_s']['trim_0.25']['median']
                                                           if centred else None),
                clearance_min_gap_m=attitude_p['clearance']['laws'][name]['min_gap_m'])
    return out


def compute(verbose=True):
    """Every answer through the screen; returns the parts the product is built from."""
    heat, layout, ledger, attitude_p, films_p, first_p, night_p, sources, layout_rolled = (
        json.loads(p.read_text()) for p in PRODUCTS)
    top_of_air = ledger['requirements']['climate_top_of_air_W_m2']
    design_sun = top_of_air/4
    spectrum = solar_spectrum()
    olr = outgoing_infrared(top_of_air)
    edges = thermal_bins()
    olr_w = infrared_weights(olr, edges)
    rings = fleet_rings(layout, heat, ledger)
    S = fleet_samples(rings, DESIGN['samples_per_ring'])
    sub = (slice(None, None, DESIGN['return_ring_stride']), slice(None, None, DESIGN['return_sample_stride']))
    window_vis = film_visible(spectrum, *DESIGN['films']['window'])
    fixed = dict(front=projected_counts(rings, S), shadow=moon_shadow(S['P']), landing=landing_share(), sub=sub,
                 sub_weight=S['r'].size/S['r'][sub].size, olr_W_m2=olr['design_W_m2'],
                 top_of_air_W_m2=top_of_air,
                 top_of_air_lux=DESIGN['lumens_per_watt']*K.SOLAR_CONSTANT*window_vis['T_visible']
                 * (1-ledger['requirements']['dimmer_share']))
    # The column absorption the stack tables use: the as-built fleet's sunlight per area of its footprint on the sky.
    solar0 = {ring: dict(solar=film_solar(spectrum, *DESIGN['films'][ring]), transmission=1., solar_extra=0.)
              for ring in ('annulus', 'window')}
    p0 = absorbed(S, fixed, solar0, attitude(S)[0])
    foot = S['area']/S['layers']
    q_col = {}
    for ring, mask in (('annulus', ~S['window']), ('window', S['window'])):
        q_col[ring] = float(np.sum(p0*mask)/np.sum(foot*mask))
        q_col['alpha_'+ring] = solar0[ring]['solar']['outer']['A'][0]
    bundles = {}
    def bundle(film, coat):
        key = (film, coat)
        if key not in bundles:
            bundles[key] = {ring: film_bundle(spectrum, olr_w, q_col, ring, film, coat) for ring in ('annulus', 'window')}
        return bundles[key]
    SAIL['as_built'] = {ring: b['pressure_coefficient'] for ring, b in bundle('as_built', None).items()}
    geos, answers = {}, {}
    cases = ([(key, film, coat, law, chi, pitch, None) for key, film, coat, law, chi, pitch in ANSWERS]
             + OWN_ECCENTRICITY+ROLLED_LAYOUT)
    for key, film, coat, law, chi, pitch, scheme in cases:
        if scheme is None:
            S_case, fixed_case = S, fixed
            gkey = (law, chi, pitch)
        else:
            if scheme == 'layout':                      # the re-flown layout carries the law's push in its flights
                drive = layout_rolled['variant']['eccentricity_drive_over_as_built']
                rings_case = fleet_rings(layout_rolled, heat, ledger)
            else:
                drive = attitude_p['schemes'][scheme]['eccentricity_drive']['over_as_built']
                rings_case = fleet_rings(layout, heat, ledger, drive)
            S_case = fleet_samples(rings_case, DESIGN['samples_per_ring'])
            fixed_case = dict(fixed, front=projected_counts(rings_case, S_case), shadow=moon_shadow(S_case['P']),
                              sub_weight=S_case['r'].size/S_case['r'][sub].size)
            gkey = (law, chi, pitch, scheme)
        if gkey not in geos:
            normal, _ = attitude(S_case, law, chi, pitch)
            geos[gkey] = specular_geometry({'P': S_case['P'][sub]}, normal[sub])
        b = bundle(film, coat)
        row = evaluate(S_case, fixed_case, b, geos[gkey], law, chi, pitch)
        if scheme is not None:
            row['eccentricity'] = dict(scheme=scheme if scheme != 'layout' else 'ring_layout_rolled',
                                       drive_over_as_built=drive,
                                       range=[float(rings_case['e'].min()), float(rings_case['e'].max())],
                                       perilune_km=[float((rings_case['a']*(1-rings_case['e'])).min()),
                                                    float((rings_case['a']*(1-rings_case['e'])).max())])
        row['dimming'] = dimming(row, design_sun)
        # A coating that takes light out of the window rings' beam dims the Moon itself, unextracted; the rest is
        # added dimming, which extraction can turn into electricity.
        row['dimming']['coating_provides_percent'] = 100*b['window']['blocked']
        row['dimming']['to_add_percent'] = row['dimming']['percent']-row['dimming']['coating_provides_percent']
        # A dimmer that absorbs in the window rings' films absorbs round their whole orbits: each point of absorptance
        # it adds raises their glow in proportion, against the sunlight it takes from the Moon while in the window.
        per_point = row['glow_by_ring_W_m2']['window']/(100*b['window']['alpha'])
        a_moon = (1-DESIGN['planetary_albedo'])*design_sun/100
        row['dimming']['absorbing_dimmer_glow_W_m2_per_percent'] = per_point
        row['dimming']['absorbing_dimmer_gives_back'] = DESIGN['infrared_efficacy']*per_point/a_moon
        row['dimming']['electricity_TW'] = [max(0., row['dimming']['to_add_percent'])
                                            * dimmer_light_TW_per_percent(ledger)*c for c in DESIGN['conversion']]
        row['costs'] = costs(key, film, coat, law, chi, pitch, b, attitude_p)
        answers[key] = row
        if verbose:
            print(f"{key:44s} glow {row['glow_W_m2']['global']:5.2f} (day {row['glow_W_m2']['day']:5.2f}, night "
                  f"{row['glow_W_m2']['night']:5.2f})  IR back {row['infrared_return_W_m2']['global']:5.2f}  light "
                  f"back {row['sunlight_return_W_m2']['global']:4.2f}  dimming {row['dimming']['percent']:5.2f}%",
                  flush=True)
    return dict(products=dict(heat=heat, layout=layout, ledger=ledger, attitude=attitude_p, films=films_p,
                              first=first_p, night=night_p, layout_rolled=layout_rolled), answers=answers, bundles=bundles, q_col=q_col, olr=olr,
                olr_w=olr_w, rings=rings, S=S, fixed=fixed, spectrum=spectrum, geos=geos, design_sun=design_sun)


def dimmer_light_TW_per_percent(ledger):
    """Light (TW) that 1% of the design sunlight carries over the Moon's disc: the ledger's dimmer light, 5% of the
    titania stack's, scaled to the design sunlight."""
    stack_light = ledger['electricity']['dimmer_optical_TW']/(100*ledger['requirements']['dimmer_share'])
    passed = 1-ledger['requirements']['dimmer_share']
    return stack_light*passed


CODE = ['research/studies/solar_shield_array/glow.py', 'protection/model.py', 'protection/spectra/annulus_film.py',
        'protection/spectra/short_wave.py', 'atmosphere/middle_atmosphere/escape.py',
        'atmosphere/middle_atmosphere/lbl_check.py', 'atmosphere/middle_atmosphere/equilibrium.py',
        'atmosphere/middle_atmosphere/run.py', 'atmosphere/radiative_convective/ck.py',
        'atmosphere/radiative_convective/longwave.py', 'atmosphere/radiative_convective/thermodynamics.py']


def input_files():
    from atmosphere.middle_atmosphere import fetch_inputs
    files = list(PRODUCTS)+[SKY]+list(OPTICAL)+[fetch_inputs.path('whi2008_ref_solar_irradiance_ver2.dat')]
    for name in COLUMNS:
        files += [MIDDLE/'profiles'/f'{name}.csv', MIDDLE/'cases'/f'{name}.json']
    return files


def rnd(x, n=6):
    """Round floats in nested structures to n significant figures, for a stable, readable product."""
    if isinstance(x, dict):
        return {k: rnd(v, n) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, n) for v in x]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return 0. if v == 0 or not math.isfinite(v) else float(f'{v:.{n}g}')
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.ndarray):
        return rnd(x.tolist(), n)
    return x


ASSUMPTIONS = [
    ('films', 'sourced', 'the comparison\'s films: 0.1 um of titania on 4 um of silica on the annulus rings, the window '
     'stack (1 um on 10 um) on the window rings, both with the protection domain\'s anti-reflection and matching '
     'layers (array_heat, annulus_film)'),
    ('optical constants', 'sourced', 'Franta silica and Siefke titania, CXRO below 24.8 nm, Malitson silica without '
     'absorption at 210-2,500 nm, as the protection domain evaluates its films'),
    ('oblique incidence', 'estimate', 'smooth specular films, s and p averaged; the protection domain\'s recursion '
     'reproduced to 1e-15 at normal incidence'),
    ('solar infrared beyond 2.5 um', 'estimate', 'a 5,772 K blackbody\'s shape scaled so the spectrum gives the solar '
     'constant: 45.4 W/m2'),
    ('the stack', 'estimate', 'at each point of the sky, the local layers as parallel infinite plates sharing the '
     'local absorbed sunlight evenly, no conduction between them; thermal bands 2-125 um, held beyond'),
    ('tiles and rings', 'sourced', 'every ring of the layout on the flown matched rings\' apolune radius and forced '
     'eccentricity by height, the inventory\'s 28.2 million tiles shared evenly and spaced evenly in time, held '
     'perpendicular to the velocity with the 1.5 degree shingle (ring_bundle); keeping tilts left out'),
    ('window rings', 'sourced', 'the middle 26% of the stack by height, array_heat\'s share'),
    ('sunlight on each tile', 'estimate', 'attenuated by the films the projected counts put in front (fleet_light\'s '
     'rule), with the Moon\'s penumbra linear'),
    ('layers in the Moon\'s view', 'estimate', 'along-track tiles over a tile\'s angular length times rings over a '
     'strip\'s angular width; a rolled stack opens to 1 + (n - 1)|cos roll| layers (design guess)'),
    ('the Moon\'s outgoing infrared', 'estimate', 'the atmosphere domain\'s correlated-k model on the clear titania-'
     'stack columns at 288 and 298 K, interpolated to 294.9 K and scaled to the GCM design case\'s (1 - 0.280) S/4 = '
     '211.2 W/m2; uniform over the Moon and Lambertian (design guess)'),
    ('the Moon\'s reflected sunlight', 'design guess', 'a Lambertian day side with the planetary albedo 0.28, flat in '
     'wavelength, its spectrum the sunlight the window stack passes'),
    ('infrared efficacy', 'GCM run', 'the dimming that holds the design temperature against the glow is 0.71 of its '
     'infrared watts matched as absorbed sunlight (the gases, glow and dimming study\'s held branch, 13 October); '
     'about 0.59 on the corrected radiation of 17 October; efficacy 1, the energy balance alone, is kept beside it'),
    ('wide-gap absorber', 'design guess', 'an absorber with a step edge at 220 or 250 nm in the titania\'s place, '
     '1-2 g/m2, assumed to hold the film\'s present opacity below 175 nm'),
    ('ultraviolet mirror', 'design guess', 'lossless, reflecting 95-98% of 245-390 nm (annulus) or 245-340 nm (window, '
     'keeping its spectrum) ahead of the titania; mass from quarter-wave HfO2 (n 2.1) and silica stacks'),
    ('heat mirrors and emitters', 'design guess', 'grey opaque coatings of the stated emissivity, sunlight absorbed '
     'and taken out of the beam, and mass, anchored to measured coatings (glow_sources.json)'),
    ('rolls', 'sourced', 'attitude_schemes\' law: zero in the 25 degree service arc and at the nodes, 30 degree cosine '
     'ramps; the night roll keeps the night half only'),
    ('dimming\'s electricity', 'sourced', '1% of the design sunlight over the Moon\'s disc is the ledger\'s dimmer light '
     'scaled, at 20-30% conversion'),
]


IR_BANDS_UM = [(5., 7.5), (7.5, 8.5), (8.5, 9.6), (9.6, 12.), (12., 15.), (15., 20.), (20., 25.), (25., 40.), (40., 125.)]
RECOMMENDED = 'rolled_outside_service_80_rolled_layout'
NEXT = ('rolled_outside_service_80_wide_gap_250_rolled_layout', 'rolled_outside_service_80_uv_mirror_outer_98_rolled_layout')


def film_summary(b, ring, heat, olr_w):
    sol = b['solar']
    th = thermal_plate(b['titania_um'], b['silica_um'])
    centres = np.sqrt(th['edges_um'][1:]*th['edges_um'][:-1])
    bands = {}
    for lo, hi in IR_BANDS_UM:
        sel = (centres >= lo) & (centres < hi)
        bands[f'{lo:g}_{hi:g}_um'] = dict(R_moon_face=float(np.mean(th['rho_in'][sel, -1])),
                                         T=float(np.mean(th['tau'][sel, -1])))
    f = heat['tiles']['films'][ring]
    e = normal_emissivity(b['titania_um'], b['silica_um'], f['temperature_K'])
    one = b['table']['rows'][0]
    return dict(titania_um=b['titania_um'], silica_um=b['silica_um'],
                solar_absorptance=dict(normal_outer_face=sol['outer']['A'][0], normal_inner_face=sol['inner']['A'][0],
                                       within_0p1_2500_nm=sol['near_A'], beyond_2500_nm=sol['tail_A'],
                                       array_heat=f['solar_absorptance']),
                absorbed_points_of_sunlight_by_band=sol['bands_points'],
                by_incidence=dict(cosines=sol['cosines'], outer_face=sol['outer']['A'], inner_face=sol['inner']['A']),
                reflectance=dict(solar_outer=sol['outer']['R'][0], visible_moon_face=b['visible']['R_visible'],
                                 moonlight_moon_face=b['visible']['R_moonlight']),
                pressure_coefficient=b['pressure_coefficient'],
                thermal=dict(normal_emissivity_at_film_K=dict(outer=e[0], moon_face=e[1], K=f['temperature_K']),
                             moon_infrared_normal_reflectance_one_layer=one['R_ir'], by_band=bands),
                temperature_sunward_crossing_K=dict(one_layer=b['T_sunward_one'], two_layers=b['T_sunward_two']))


def against_twilight(night_lux):
    """The night light the films return, band by band from midnight, beside the Moon's own twilight there: the solved
    sky's total horizontal illuminance for the design Moon's 1.2 atm column with the Sun at each band's edges below the
    horizon (illumination/sky). The returns' share of the twilight is taken against the band's darker edge."""
    rows = sorted(json.loads(SKY.read_text())['worlds']['moon_1.2atm']['samples'], key=lambda r: r['sun_deg'])
    suns = np.array([r['sun_deg'] for r in rows])
    lux = np.log(np.array([r['total_horizontal_lux'] for r in rows]))
    edges = night_lux['by_angle_from_midnight_deg']['edges']
    out = []
    for lo, hi, returned in zip(edges[:-1], edges[1:], night_lux['by_angle_from_midnight_deg']['lux']):
        deep, shallow = -(90.-lo), -(90.-hi)                 # the Sun's height at the band's two edges
        sky = [float(np.exp(np.interp(x, suns, lux))) if x >= suns[0] else None for x in (deep, shallow)]
        out.append(dict(from_midnight_deg=[lo, hi], sun_deg=[deep, shallow], returned_lux=returned,
                        twilight_lux=sky, share_of_twilight=(returned/sky[0] if sky[0] else None)))
    return out


def build_product(c):
    P = c['products']
    heat, layout, first = P['heat'], P['layout'], P['first']
    S, fixed, rings = c['S'], c['fixed'], c['rings']
    A = c['answers']
    built = c['bundles'][('as_built', None)]
    as_built = A['as_built']
    # The fleet as the screen builds it.
    normal = S['normal']
    phat = S['P']/S['r'][..., None]
    deviation = np.degrees(np.arccos(np.clip(np.sum(normal*phat, axis=-1), -1, 1)))
    cover = float(np.sum(S['area']/(S['r']*1e3)**2/S['layers']))
    tile_sky = float(np.sum(S['area']/(S['r']*1e3)**2))
    lune = 2*float(np.max(np.abs(rings['tilt'])))/np.pi
    geo = c['geos'][('as_built', None, 0.)]
    sub = fixed['sub']
    fw = S['area'][sub]/S['layers'][sub]
    landed = float(np.sum(fw*(geo['ir_day']+geo['ir_night']))/np.sum(fw*geo['ir_in']))
    p_abs = absorbed(S, fixed, built, normal)
    alpha = np.where(S['window'], built['window']['alpha'], built['annulus']['alpha'])
    lit = float(np.sum(p_abs)/np.sum(S['area']*K.SOLAR_CONSTANT*alpha))
    # What the Moon sees of the glow: its glow-weighted view factor and share, against the first screen's.
    n_eff = S['layers']
    mu_view = np.clip(np.sum(normal*phat, axis=-1), 0, 1)
    share = np.where(S['window'], table_lookup(built['window']['table'], n_eff, mu_view),
                     table_lookup(built['annulus']['table'], n_eff, mu_view))
    view = (K.MOON_RADIUS/1e3/S['r'])**2
    w = p_abs*share
    fleet = dict(
        rings=int(rings['heights'].size), window_rings=int(np.sum(rings['window'])), tiles=float(np.mean(heat['fleet']['tiles'])),
        apolune_km=[float(rings['r_apo'].min()), float(rings['r_apo'].max())],
        perilune_km=[float((rings['r_apo']*(1-rings['e'])/(1+rings['e'])).min()),
                     float((rings['r_apo']*(1-rings['e'])/(1+rings['e'])).max())],
        forced_eccentricity=[float(rings['e'].min()), float(rings['e'].max())],
        ring_length_over_circle_at_apolune=float(np.mean(2*np.pi*rings['a']*(1-rings['e']**2/4)/(2*np.pi*rings['r_apo']))),
        along_track_layers=dict(apolune=float(np.mean(rings['tiles_per_ring']*rings['tile_km']*rings['r_apo']
                                                      / (2*np.pi*rings['a']**2*np.sqrt(1-rings['e']**2)))),
                                designed=10./8.5),
        normal_off_moon_centre_deg=dict(max=float(deviation.max()), at_sunward_crossing=float(np.median(
            deviation[np.abs(S['phase']) < np.radians(5)])), within_service_arc_max=float(np.max(
            deviation[np.abs(S['phase']) < np.radians(DESIGN['service_half_angle_deg'])]))),
        absorbed_PW=as_built['absorbed_PW'], lit_factor=lit, array_heat_absorbed_PW=heat['tiles']['absorbed_PW'],
        array_heat_lit_factor=heat['design']['lit_factor'],
        layers_in_moon_view=dict(mean=tile_sky/cover, first_screen=first['heat']['layers_in_the_moons_view']),
        sky_covered=dict(share=cover/(4*np.pi), lune=lune),
        specular_image_landing_share=landed,
        glow_weighted=dict(view_factor=float(np.sum(w*view)/np.sum(w)), moon_share=float(np.sum(p_abs*share*view)/np.sum(p_abs*view)),
                           first_screen_view_factor=first['heat']['view_factor'],
                           first_screen_moon_share=[first['heat']['moon_facing_share']['deep_stack'],
                                                    first['heat']['moon_facing_share']['one_layer']]))
    # The stacks as built.
    stacks = {}
    for ring, b in built.items():
        stacks[ring] = dict(column_absorbed_W_m2=b['q_column'], rows=[dict(
            layers=row['n'], toward_moon_share=row['by_mu'][-1], inward_hemispheric_share=row['inward'],
            temperature_K=dict(inner=row['T'][0], outer=row['T'][-1], max=max(row['T'])),
            moon_infrared_reflectance=row['R_ir'], moonlight_reflectance=row['R_sw'], visible_reflectance=row['R_vis'])
            for row in b['table']['rows']])
    two = two_layer_check(heat)
    fs = first_screen(heat, layout)
    product_two = {k: ({kk: vv for kk, vv in v.items() if kk != 'one_layer_share'} if isinstance(v, dict) else v)
                   for k, v in first['heat']['heat_mirror'].items()}
    residual = max(abs(row['residual']) for bs in c['bundles'].values() for b in bs.values() for row in b['table']['rows'])
    films = {ring: film_summary(b, ring, heat, c['olr_w']) for ring, b in built.items()}
    fp = P['films']
    ref = dict(annulus=next(f for f in fp['films'] if f['coated'] and f['titania_um'] == .1 and f['silica_um'] == 4.)['bolometric'],
               window=fp['window_stack']['bolometric'])
    checks = dict(
        optics_reproduce_annulus_film=all(abs(built[r]['solar']['near_A']-ref[r]['A']) < 2e-6
                                          and abs(built[r]['solar']['near_T']-ref[r]['T']) < 2e-6 for r in ref),
        emissivities_reproduce_array_heat=all(
            abs(films[r]['thermal']['normal_emissivity_at_film_K']['outer']-heat['tiles']['films'][r]['emissivity_sunward']) < 1e-9
            and abs(films[r]['thermal']['normal_emissivity_at_film_K']['moon_face']-heat['tiles']['films'][r]['emissivity_shaded']) < 1e-9
            for r in ref),
        solar_spectrum_gives_solar_constant=abs(c['spectrum']['near_W_m2']+c['spectrum']['tail_W_m2']-K.SOLAR_CONSTANT) < 1e-9,
        stack_energy_balance=residual < 1e-9,
        two_layer_estimate_reproduced=two == product_two,
        first_screen_reproduced=[round(x, 2) for x in fs] == first['heat']['fleet_infrared_w_m2'],
        sky_covered_matches_lune=abs(cover/(4*np.pi)/lune-1) < .01,
        absorbed_within_array_heat_range=heat['tiles']['absorbed_PW'][0]*.99 <= as_built['absorbed_PW'] <= heat['tiles']['absorbed_PW'][1],
        outgoing_infrared_closes=abs(sum(c['olr']['band_W_m2'])-c['olr']['design_W_m2']) < 1e-9)
    rec = A[RECOMMENDED]
    mirrors = [k for k in A if k.split('_')[0] in ('mesh', 'metasurface', 'tco', 'silver')]
    mirror_back = [min(A[k]['infrared_return_W_m2']['global'] for k in mirrors),
                   max(A[k]['infrared_return_W_m2']['global'] for k in mirrors)]
    rise = [A[k]['heat_total_W_m2']['global']-as_built['heat_total_W_m2']['global'] for k in mirrors]
    mirror_rise = [min(rise), max(rise)]
    records = [
        dict(record='provisioning README finding 1 and the register row of 9 October',
             states='the Moon receives 5.7-9.9 W/m2 of glow',
             screen=f"{as_built['glow_W_m2']['global']:.1f} W/m2 of glow, and {as_built['infrared_return_W_m2']['global']:.1f} "
                    f"W/m2 of the Moon's own infrared and {as_built['sunlight_return_W_m2']['global']:.2f} W/m2 of its "
                    f"reflected sunlight returned by the films"),
        dict(record='provisioning README finding 1 and the register row of 9 October',
             states='a heat mirror on Moon-facing faces cuts the Moon\'s share 1.5-2 times, or 3-9 times on the nearest layer',
             screen=f"the cuts hold for the glow, but a mirror facing the Moon images the Moon onto itself: the candidates "
                    f"return {mirror_back[0]:.0f}-{mirror_back[1]:.0f} W/m2 of its own infrared and raise the heat the "
                    f"Moon receives by {mirror_rise[0]:.1f}-{mirror_rise[1]:.1f} W/m2"),
        dict(record='provisioning screen.heat()', states='a deep stack sends half its heat toward the Moon',
             screen=f"along the Moon's direction the semi-transparent stack sends {stacks['annulus']['rows'][-1]['toward_moon_share']:.2f} "
                    f"(annulus) and {stacks['window']['rows'][-1]['toward_moon_share']:.2f} (window) of its heat at 32 layers"),
        dict(record='provisioning screen.heat()', states='tiles 19,355-21,392 km from the Moon',
             screen=f"the layout's radii are apolune radii; the night half rides at perilune, "
                    f"{fleet['perilune_km'][0]:,.0f}-{fleet['perilune_km'][1]:,.0f} km"),
        dict(record='domain synthesis heat ledger', states='the dimming that answers heat is its W/m2 over 293 W/m2',
             screen=f"infrared the Moon absorbs whole is matched by sunlight it absorbs (1 - {DESIGN['planetary_albedo']}) of, "
                    f"so the dimming is {1/(1-DESIGN['planetary_albedo']):.2f} times that at efficacy 1 and "
                    f"{DESIGN['infrared_efficacy']/(1-DESIGN['planetary_albedo']):.2f} times at the GCM's "
                    f"{DESIGN['infrared_efficacy']:.2f}"),
        dict(record='joint synthesis, the night with the ring fleet',
             states='the fleet lights the night by diffuse scatter and glints of sunlight',
             screen=f"the films' specular reflection of the Moon's own day side also lands on the night: "
                    f"{as_built['night_lux_from_returns']['mean']:.1f} lux averaged over the night hemisphere as built"),
    ]
    return dict(
        schema=SCHEMA,
        producer=dict(files={f: digest(ROOT/f) for f in CODE}, constants=constants_used(CODE),
                      correlated_k_table=c['olr']['table_key']),
        inputs={str(Path(p).relative_to(ROOT)): digest(p) for p in input_files()},
        evidence=EVIDENCE, reading_rule=READING_RULE,
        assumptions=[dict(name=n, evidence=k, statement=t) for n, k, t in ASSUMPTIONS],
        design=json.loads(json.dumps(DESIGN)),
        results=rnd(dict(
            films=films, outgoing_infrared=c['olr'], stacks=stacks, fleet=fleet,
            first_screen=dict(rebuilt_W_m2=fs, product_W_m2=first['heat']['fleet_infrared_w_m2']),
            two_layer_check=dict(rebuilt=two, product=product_two),
            design_sunlight_W_m2=c['design_sun'],
            dimmer_light_TW_per_percent=dimmer_light_TW_per_percent(P['ledger']),
            answers=A,
            recommended=dict(answer=RECOMMENDED, then=list(NEXT), heat_W_m2=rec['heat_total_W_m2'],
                             glow_W_m2=rec['glow_W_m2'], infrared_return_W_m2=rec['infrared_return_W_m2'],
                             sunlight_return_W_m2=rec['sunlight_return_W_m2'], dimming=rec['dimming'],
                             costs=rec['costs'], night_lux=rec['night_lux_from_returns'],
                             night_lux_against_twilight=against_twilight(rec['night_lux_from_returns'])),
            records=records)),
        checks=checks)


def main():
    t0 = time.process_time()
    product = build_product(compute())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1)+'\n')
    print(json.dumps(product['checks'], indent=1))
    print(f'cpu {time.process_time()-t0:.0f} s -> {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
