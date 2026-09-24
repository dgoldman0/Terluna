"""Middle atmosphere: rate coefficients, chemistry Jacobian and transport, photolysis and correlated-k checks.

Tests that need external inputs (ultraviolet cross-sections, HITRAN lines, the
built correlated-k tables) skip with that reason when they are absent.
"""
from __future__ import annotations
import math
import unittest
import numpy as np

from atmosphere.radiative_convective import thermodynamics as th, ck
from atmosphere.radiative_convective import fetch_inputs as rc_inputs
from atmosphere.middle_atmosphere import chemistry as ch, fetch_inputs


def _present(module):
    try:
        return module.restore()['status'] == 'PASS'
    except (OSError, ValueError):
        return False


NEEDS_UV = unittest.skipUnless(_present(fetch_inputs) and _present(rc_inputs),
                               'ultraviolet or solar inputs not restored (fetch_inputs --download)')


def _tables_built():
    try:
        ck.load()
        return True
    except (OSError, ValueError, FileNotFoundError, KeyError):
        return False


def _column(planet=th.EARTH, levels=61, top=1.0):
    air = th.earthlike_air(101325.0, 400.0)
    p = np.geomspace(101325.0, top, levels)
    t = np.clip(288.0 - 30.0 * np.log(101325.0 / p), 210.0, 288.0)
    x = np.clip(0.01 * (p / 101325.0) ** 3, 4e-6, None)
    return th.column_from_levels(planet, air, p, t, x)


class RateCoefficientTests(unittest.TestCase):
    def rate(self, reactants, products, t, m=2.46e19, h2o=0.0):
        for r, p_, k, _ in ch.REACTIONS:
            if r == reactants and p_ == products:
                return float(np.asarray(k(np.array([t]), np.array([m]), np.array([h2o])))[0])
        raise KeyError(reactants)

    def test_bimolecular_values_at_298k(self):
        # JPL 19-5 recommended k(298 K), cm^3 s^-1.
        self.assertAlmostEqual(self.rate(('O', 'O3'), ('O2', 'O2'), 298.0) / 8.0e-15, 1.0, delta=0.01)
        self.assertAlmostEqual(self.rate(("NO", "O3"), ("NO2", "O2"), 298.0) / 1.9e-14, 1.0, delta=0.03)
        self.assertAlmostEqual(self.rate(('OH', 'O3'), ('HO2', 'O2'), 298.0) / 7.3e-14, 1.0, delta=0.01)
        self.assertAlmostEqual(self.rate(('O1D', 'N2'), ('O', 'N2'), 298.0) / 3.1e-11, 1.0, delta=0.01)

    def test_termolecular_limits(self):
        k = ch.termolecular(1.8e-30, 3.0, 2.8e-11, 0.0)
        t = np.array([250.0])
        # The broadening factor 0.6^(1/(1+log10(k0 M/kinf)^2)) tends to 1 only logarithmically.
        low = float(k(t, np.array([1.0]), None)[0]) / (1.8e-30 * (250 / 300) ** -3)
        high = float(k(t, np.array([1e40]), None)[0]) / 2.8e-11
        self.assertAlmostEqual(low, 1.0, delta=2e-3)
        self.assertAlmostEqual(high, 1.0, delta=2e-3)
        mid = 2.8e-11 / (1.8e-30 * (250 / 300) ** -3)          # k0 M = kinf: half the low limit times Fc
        self.assertAlmostEqual(float(k(t, np.array([mid]), None)[0]) / (0.5 * 2.8e-11 * 0.6), 1.0, delta=1e-12)
        # Chapman recombination: JPL 19-5 k0(298) = 6.1e-34 (298/300)^-2.4 cm^6 s^-1.
        self.assertAlmostEqual(self.rate(('O', 'O2'), ('O3',), 298.0, m=1.0) / 6.19e-34, 1.0, delta=0.01)

    def test_n2o5_equilibrium(self):
        t, m = np.array([250.0]), np.array([5e18])
        forward = float(ch._NO2_NO3(t, m, None)[0])
        backward = float(ch._n2o5_decomposition(t, m, None)[0])
        self.assertAlmostEqual(forward / backward / (2.7e-27 * math.exp(11000 / 250.0)), 1.0, delta=1e-12)


class ChemistryNumericsTests(unittest.TestCase):
    def test_jacobian_matches_finite_differences(self):
        rng = np.random.default_rng(1)
        nl = 3
        fixed = dict(M=np.full(nl, 1e17), O2=np.full(nl, 2e16), N2=np.full(nl, 8e16), H2O=np.full(nl, 1e12))
        # Densities of similar size keep rounding in the differences small; rates are at most
        # quadratic in any one species, so central differences are exact apart from rounding.
        n = 10 ** rng.uniform(7, 9, size=(len(ch.SPECIES), nl))
        t = np.array([200.0, 230.0, 260.0])
        j = {name: np.full(nl, 1e-4) for name in ch.PHOTOLYSIS}
        _, jac = ch._rates_and_jacobian(n, fixed, t, j)
        for i in range(len(ch.SPECIES)):
            dn = np.zeros_like(n); dn[i] = n[i] * 1e-3
            plus, _ = ch._rates_and_jacobian(n + dn, fixed, t, j, include_jac=False)
            minus, _ = ch._rates_and_jacobian(n - dn, fixed, t, j, include_jac=False)
            fd = (plus - minus) / (2 * dn[i])
            gross = (np.abs(plus) + np.abs(minus)) / (2 * dn[i])
            self.assertTrue(np.all(np.abs(jac[:, :, i].T - fd) <= 1e-7 * np.abs(fd) + 1e-11 * gross))

    def test_transport_conserves_a_tracer(self):
        col = _column(th.MOON)
        geo = ch.Column(col, 2e4, ch.Mixing(scale=36.0), ch.Boundary())
        rows, cols, vals = geo.transport()
        n = np.linspace(1.0, 3.0, geo.nl) * geo.fixed['M'] * 1e-9
        dndt = np.bincount(rows, weights=vals * n[cols], minlength=geo.nl)
        self.assertLess(abs(np.sum(dndt * geo.volume)) / np.sum(np.abs(dndt) * geo.volume), 1e-12)
        uniform = geo.fixed['M'] * 1e-9
        still = np.bincount(rows, weights=vals * uniform[cols], minlength=geo.nl)
        self.assertLess(np.max(np.abs(still)), 1e-12 * np.max(np.abs(dndt)))


class SteadyStateTests(unittest.TestCase):
    def test_solver_reaches_a_balanced_steady_state(self):
        # Fixed photolysis rates that fall off with depth stand in for the sunlight; the solver's
        # answer must make chemistry, transport, deposition and washout cancel for every species.
        col = _column(th.EARTH, levels=31, top=10.0)
        nl = col['layer_p_pa'].size
        depth = np.exp(-col['layer_p_pa'] / 3e3)

        class Sun:
            def evaluate(self, col, densities, albedo):
                rates = {name: 1e-4 * depth for name in ch.PHOTOLYSIS}
                rates['O2 -> O + O'] = 1e-10 * depth
                rates['N2O -> N2 + O1D'] = 1e-7 * depth
                return dict(rates=rates, heating=np.zeros(nl), heating_nodes=np.zeros((nl, 6)))

        mixing, boundary = ch.Mixing(), ch.Boundary()
        state = ch.solve(col, Sun(), 0.1, 2e4, mixing, boundary)
        self.assertTrue(state.info['converged'])
        geo = ch.Column(col, 2e4, mixing, boundary)
        rows, cols, vals = geo.transport()
        rates, jac = ch._rates_and_jacobian(state.n, geo.fixed, geo.t, Sun().evaluate(col, {}, 0.1)['rates'])
        for name in ('O3', 'OH', 'NO2', 'HNO3', 'N2O'):
            i = ch.INDEX[name]
            total = rates[i] + np.bincount(rows, weights=vals * state.n[i][cols], minlength=nl)
            if name in dict(boundary.deposition):
                total[0] -= dict(boundary.deposition)[name] * state.n[i, 0] * geo.area[0] / geo.volume[0]
            if name in boundary.soluble:
                total -= np.where(geo.tropo, geo.washout, 0.0) * state.n[i]
            if name in dict(boundary.fixed):
                total[0] = 0.0                                   # fixed mixing ratio at the ground
            # Against the gross chemical loss rate; the solver's Newton tolerance is 1e-3 in density.
            gross = np.abs(jac[:, i, i]) * state.n[i] + np.abs(vals).max() * state.n[i]
            self.assertLess(np.max(np.abs(total) / np.maximum(gross, 1e-30)), 1e-2, name)


@NEEDS_UV
class PhotolysisTests(unittest.TestCase):
    def test_thin_absorbing_column_sees_the_top_of_atmosphere_flux(self):
        from atmosphere.middle_atmosphere import photolysis as ph
        col = _column(top=1.0)
        sun = ph.Sunlight(None, n_mu0=4)
        tau_abs, tau_sca, _ = sun.optics(col, {})
        self.assertTrue(np.all(tau_abs >= 0) and np.all(tau_sca > 0))
        out = sun.evaluate(col, {}, 0.0)
        # In the topmost layer the O2 photolysis rate approaches the unattenuated value times the
        # global-mean geometry factor 0.5 * mean(1) = 0.5 (the actinic direct beam does not depend on mu0).
        j_top = out['rates']['O2 -> O + O'][-1]
        j_free = 0.5 * float(np.sum(sun.xs.sigma('O2', np.array([250.0]))[0] * sun.photons))
        # Light scattered up from below adds a few per cent at the top.
        self.assertGreater(j_top / j_free, 1.0)
        self.assertLess(j_top / j_free, 1.08)

    def test_energy_conservation(self):
        from atmosphere.middle_atmosphere import photolysis as ph
        col = _column(top=1.0)
        sun = ph.Sunlight(None, n_mu0=4)
        dens = {'O3': 5e12 * np.exp(-((np.log(col['layer_p_pa']) - np.log(3e3)) / 1.2) ** 2)}
        out = sun.evaluate(col, dens, 0.1)
        absorbed = out['heating_all'].sum() + out['surface_net'].sum()
        self.assertAlmostEqual(absorbed / out['toa_net'].sum(), 1.0, delta=1e-9)


@unittest.skipUnless(_tables_built(), 'correlated-k tables not built (ck.build())')
class CorrelatedKTests(unittest.TestCase):
    def test_quadrature_and_overlap(self):
        self.assertAlmostEqual(ck.W.sum(), 1.0, delta=1e-12)
        a = np.random.default_rng(2).uniform(0, 3, size=(2, ck.NB, ck.NG))
        a.sort(axis=-1)
        same = ck.overlap(a, np.zeros_like(a))
        self.assertTrue(np.allclose(same, a, rtol=1e-6, atol=1e-12))
        # Rebinning conserves the band-mean optical depth.
        b = np.sort(np.random.default_rng(3).uniform(0, 1, size=(2, ck.NB, ck.NG)), axis=-1)
        both = ck.overlap(a, b)
        self.assertTrue(np.allclose((both * ck.W).sum(-1), (a * ck.W).sum(-1) + (b * ck.W).sum(-1), rtol=1e-9))

    def test_olr_matches_line_by_line_table(self):
        # The committed line-by-line sweep holds OLR for these columns; the correlated-k scheme must
        # reproduce it (no ozone in either).
        import csv
        from pathlib import Path
        from atmosphere.radiative_convective import climate as cl
        table = Path(ck.__file__).resolve().parent / 'results' / 'inverse_climate.csv'
        rows = {r['key']: float(r['olr']) for r in csv.DictReader(open(table, newline=''))}
        model = ck.CKLongwave()
        for name, planet, dry in (('earth_1.0atm', th.EARTH, 101325.0), ('moon_1.2atm', th.MOON, 121590.0)):
            sc = cl.Scenario(name, planet, dry, 400.0)
            for ts in (270.0, 300.0):
                col = sc.column(ts, sc.longwave_levels)
                ref = rows[f'{name}|manabe_wetherald|co2=400|strat=200|albedo=0.13|ts={ts:g}|shield=none']
                self.assertLess(abs(model.fluxes(col)['olr'] - ref), 1.0, f'{name} {ts}')

    def test_band_planck_sums_to_sigma_t4(self):
        from shared.constants import STEFAN_BOLTZMANN
        total = ck.band_planck(np.array([288.0]))[0].sum() * math.pi
        # The bands stop at 3000 cm^-1, which omits 0.02% of 288-K emission.
        self.assertAlmostEqual(total / (STEFAN_BOLTZMANN * 288.0 ** 4), 1.0, delta=2e-3)
