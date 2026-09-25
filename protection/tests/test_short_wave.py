"""Short-wave transmission of the stored shield: construction, X-ray absorption and opacity checks."""
import csv
import unittest
from pathlib import Path
import numpy as np

from protection.spectra import short_wave as sw

HERE = Path(__file__).resolve().parents[1]
INPUTS_PRESENT = all((sw.SOURCES / name).is_file() for name in sw.INPUTS)
NEEDS_INPUTS = unittest.skipUnless(INPUTS_PRESENT, 'protection inputs not restored (fetch_inputs.py --download)')


@NEEDS_INPUTS
class ShortWaveTests(unittest.TestCase):
    def test_layer_sequence_matches_the_design(self):
        # With the design's own long-wave constants the rebuilt stack must equal model.optical_stack.
        from protection import model
        w = np.linspace(215.0, 1000.0, 60)
        ar = sw.stored_ar_layers()
        silica = model.si_n(w / 1000.0).astype(complex)
        titania = model.ti_nk(w / 1000.0)
        np.testing.assert_allclose(sw.design_transmission(w, silica, titania, ar),
                                   model.optical_stack(w / 1000.0, ar)[1], rtol=1e-12, atol=1e-300)

    def test_absorption_matches_the_design_xray_table(self):
        # The design's table: Beer-Lambert through 10 um of silica and 1 um of titania from CXRO f2.
        with open(HERE / 'results' / 'xray_absorption.csv', newline='') as handle:
            rows = list(csv.DictReader(handle))
        w = np.array([float(r['wavelength_nm']) for r in rows])
        mu = lambda n: 4 * np.pi * n.imag / (w * 1e-9)                   # 1/m
        mu_si, mu_ti = mu(sw.cxro_index(w, **sw.SILICA)), mu(sw.cxro_index(w, **sw.TITANIA))
        np.testing.assert_allclose(1e6 / mu_si, [float(r['silica_attenuation_length_um']) for r in rows], rtol=1e-7)
        expected = np.array([float(r['T_10um_silica_1um_titania']) for r in rows])
        np.testing.assert_allclose(np.exp(-mu_si * 10e-6 - mu_ti * 1e-6), expected, rtol=1e-6, atol=1e-300)

    def test_film_passes_only_hard_xrays(self):
        wavelength, t, silica_gap_max = sw.compute()
        self.assertTrue(np.all((t >= 0.0) & (t <= 1.0)))
        self.assertGreater(float(t[0]), 0.9)                               # 0.1 nm
        self.assertLess(float(t[wavelength >= 2.5].max()), 1e-6)
        self.assertLess(float(t[wavelength >= 5.0].max()), 1e-20)
        # Where titania rests on the CXRO estimate, the measured silica alone is already opaque.
        self.assertLess(silica_gap_max, 1e-100)


if __name__ == '__main__':
    unittest.main()
