"""The flight models keep gravity out of buoyant lift, reproduce Lamb's coefficients and beam theory, and obey the
rule of similarity: under a sixth of the gravity, the same flyer six times larger carries the same shares."""
from __future__ import annotations
import unittest
import numpy as np

from engineering.flight import buoyant as bu
from engineering.flight import winged as wi

G_EARTH, G_MOON = 9.80665, 1.6242
SIMILAR = G_EARTH / G_MOON


class Lift(unittest.TestCase):
    def test_gravity_cancels_in_buoyant_lift(self):
        earth, moon = bu.Air(G_EARTH, 1.225), bu.Air(G_MOON, 1.225)
        self.assertEqual(bu.net_lift_kg_m3(earth), bu.net_lift_kg_m3(moon))

    def test_pure_gases_in_standard_air(self):
        air = bu.Air(G_EARTH, 1.225, 288.15)
        # Pure hydrogen and helium at the air's temperature and pressure: 1.225 (1 - M_gas / M_air).
        self.assertAlmostEqual(bu.net_lift_kg_m3(air, bu.HYDROGEN, 1.0), 1.225 * (1 - 2.01588 / 28.9644), places=6)
        self.assertAlmostEqual(bu.net_lift_kg_m3(air, bu.HELIUM, 1.0), 1.225 * (1 - 4.002602 / 28.9644), places=6)
        self.assertAlmostEqual(bu.net_lift_kg_m3(air, bu.HYDROGEN, 1.0), 1.1397, places=3)

    def test_impurity_and_superheat(self):
        air = bu.Air(G_EARTH, 1.2, 280.0)
        self.assertLess(bu.net_lift_kg_m3(air, purity=0.95), bu.net_lift_kg_m3(air, purity=1.0))
        self.assertGreater(bu.net_lift_kg_m3(air, superheat_k=10.0), bu.net_lift_kg_m3(air))


class Hulls(unittest.TestCase):
    def test_lamb_coefficients(self):
        # Lamb, Hydrodynamics, section 373: prolate spheroids of fineness 4, 6 and 10.
        for f, k1, k2 in ((4.0, 0.082, 0.860), (6.0, 0.045, 0.918), (10.0, 0.021, 0.960)):
            a, b = bu.lamb_coefficients(f)
            self.assertAlmostEqual(float(a), k1, delta=0.002)
            self.assertAlmostEqual(float(b), k2, delta=0.002)
        a, b = bu.lamb_coefficients(1.0001)
        self.assertAlmostEqual(float(a), 0.5, places=3)
        self.assertAlmostEqual(float(b), 0.5, places=3)

    def test_munk_moment_for_the_los_angeles(self):
        """NACA Report 333 takes k2 - k1 = 0.9 for the Los Angeles (fineness 7.2); Munk's moment is q Vol (k2 - k1)
        sin 2 alpha."""
        k1, k2 = bu.lamb_coefficients(200.0 / 27.64)
        self.assertAlmostEqual(float(k2 - k1), 0.9, delta=0.01)
        hull, air = bu.Hull(200.0, 200.0 / 27.64), bu.Air(G_EARTH, 1.2)
        m = bu.munk_moment(hull, air, 30.0, 0.1)
        self.assertAlmostEqual(float(m), 0.5 * 1.2 * 900.0 * float(hull.volume_m3) * float(k2 - k1) * np.sin(0.2),
                               places=3)

    def test_hull_of_volume(self):
        h = bu.Hull.of_volume(200_000.0, 245.0 / 41.2)
        self.assertAlmostEqual(float(h.volume_m3), 200_000.0, places=3)
        self.assertAlmostEqual(float(h.length_m), 246.1, delta=1.0)     # LZ 129: 245 m for about 200,000 m3

    def test_similarity_of_rigid_hulls(self):
        """In air of the same density, at the same speed and gust, the Moon's hull 6.04 times longer needs the same
        share of its lift for weight-driven and aerodynamic structure, and a smaller share for its cover and cells."""
        design = bu.Design(frame=bu.DURALUMIN, girder_factor=1.0, frame_factor=0.56, frame_gas_factor=7.5)
        earth, moon = bu.Air(G_EARTH, 1.225), bu.Air(G_MOON, 1.225)
        lift = float(bu.net_lift_kg_m3(earth))
        for L in (100.0, 245.0, 800.0):
            e = bu.rigid(bu.Hull(L, 6.0), earth, lift, design)
            m = bu.rigid(bu.Hull(L * SIMILAR, 6.0), moon, lift, design)
            self.assertAlmostEqual(float(m['weight_driven_share'] / e['weight_driven_share']), 1.0, places=9)
            self.assertAlmostEqual(float(m['aerodynamic_share'] / e['aerodynamic_share']), 1.0, places=9)
            self.assertAlmostEqual(float(e['areal_share'] / m['areal_share']), SIMILAR, places=9)

    def test_woodward_gust_moment_for_the_shenandoah(self):
        """Woodward (1976): the Shenandoah at her top speed of 91 ft/s in a 35 ft/s gust, in air of 0.0021 slug/ft3,
        gets a peak gust moment of 3,950,000 lb ft (2,290,000 ft3, 680 ft)."""
        ft, slug_ft3, lbft = 0.3048, 515.379, 1.3558179
        length, volume, fineness = 680 * ft, 2.29e6 * ft ** 3, 8.6
        c = volume / (length * (length / fineness) ** 2)
        hull = bu.Hull(length, fineness, volume_coefficient=c)
        m = bu.midship_moments(hull, bu.Air(G_EARTH, 0.0021 * slug_ft3), 0.0,
                               bu.Design(speed_m_s=91 * ft, gust_m_s=35 * ft))
        self.assertAlmostEqual(float(m['gust']) / lbft / 3.95e6, 1.0, delta=0.005)

    def test_weight_driven_share_grows_with_length_and_gravity(self):
        design = bu.Design(frame=bu.DURALUMIN)
        air = bu.Air(G_EARTH, 1.225)
        lift = float(bu.net_lift_kg_m3(air))
        r = bu.rigid(bu.Hull(np.array([200.0, 400.0]), 6.0), air, lift, design)
        self.assertAlmostEqual(float(r['weight_driven_share'][1] / r['weight_driven_share'][0]), 2.0, places=9)
        self.assertAlmostEqual(float(r['aerodynamic_share'][1] / r['aerodynamic_share'][0]), 1.0, places=9)

    def test_pressure_hull_tension(self):
        """The crown's tension is the pressure plus the gas head over the diameter, times the radius."""
        air = bu.Air(G_EARTH, 1.225)
        lift = float(bu.net_lift_kg_m3(air))
        h = bu.Hull(60.0, 4.0)
        r = bu.pressure_hull(h, air, lift, bu.Fabric('test', 0.2, 1e5), bu.Design(speed_m_s=20.0, gust_m_s=0.0))
        q = 0.5 * 1.225 * 400.0
        self.assertAlmostEqual(float(r['pressure_pa']), 1.25 * q, places=6)        # no gust: the margin decides
        self.assertAlmostEqual(float(r['hoop_n_m']), (1.25 * q + lift * G_EARTH * 15.0) * 7.5, places=6)


class Wings(unittest.TestCase):
    def test_elliptic_root_moment_and_rectangular_planform(self):
        # A rectangular wing: the integral of the elliptic loading's moment over the semi-span is pi/32 of
        # (lift x span), so the planform factor is 1/32.
        self.assertAlmostEqual(wi.planform_factor(1.0), 1.0 / 32.0, places=5)
        self.assertLess(wi.planform_factor(0.3), wi.planform_factor(0.1))

    def test_pratt_gust_factor(self):
        # At large mass ratio the alleviation factor tends to 0.88; gravity cancels in the extra lift.
        big = wi.gust_load_factor(1e6, 1.0, 100.0, 10.0, 1.0, G_EARTH, 5.0)
        self.assertAlmostEqual(float(big * 2 * G_EARTH * 1e6 / (1.0 * 10.0 * 100.0 * 5.0)), 0.88, places=4)
        dn_e = wi.gust_load_factor(600.0, 8.0, 200.0, 15.0, 0.4, G_EARTH)
        dn_m = wi.gust_load_factor(600.0, 8.0, 200.0, 15.0, 0.4, G_MOON)
        self.assertAlmostEqual(float(dn_m * G_MOON), float(dn_e * G_EARTH), places=9)

    def test_similarity_of_wings(self):
        """At one wing loading in pascals and one density, the wing 6.04 times the span carries the same share."""
        for b in (30.0, 88.0):
            e, m = wi.Wing(b), wi.Wing(b * SIMILAR)
            me, mm = 6000.0 * e.area_m2 / G_EARTH, 6000.0 * m.area_m2 / G_MOON
            se = wi.wing_mass_kg(e, me, G_EARTH, 3.75, areal_kg_m2=0.0) / me
            sm = wi.wing_mass_kg(m, mm, G_MOON, 3.75, areal_kg_m2=0.0) / mm
            self.assertAlmostEqual(float(sm / se), 1.0, places=9)
            self.assertAlmostEqual(float(mm / me), SIMILAR ** 3, places=9)

    def test_hover_power_scales_as_gravity_to_one_and_a_half(self):
        e = wi.hover_power_w(100.0, G_EARTH, 1.225, 10.0)
        m = wi.hover_power_w(100.0, G_MOON, 1.225, 10.0)
        self.assertAlmostEqual(float(m / e), (G_MOON / G_EARTH) ** 1.5, places=9)


if __name__ == '__main__':
    unittest.main()
