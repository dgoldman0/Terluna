"""Band-energy ledger and optical-input audit tests."""
from __future__ import annotations
import unittest
from atmosphere.spectral_interface import Band, deposited_heat, historical_solar_bands, gap_energy_bounds, audit_optical_files


class SpectralTests(unittest.TestCase):
    def test_missing_input_rejected(self):
        with self.assertRaises(ValueError):deposited_heat([], (1,120))
        b=Band(1,120,1e-3,None,.9,.1,2,'incomplete test')
        with self.assertRaises(ValueError):deposited_heat([b],(1,120))
    def test_gap_and_overlap_rejected(self):
        for start in [49,51]:
            bs=[Band(1,50,.001,.01,1,.1,2,'test'),Band(start,120,.001,.01,1,.1,2,'test')]
            with self.assertRaises(ValueError):deposited_heat(bs,(1,120))
    def test_energy_ledger(self):
        q,rows=deposited_heat([Band(1,120,.004,.001,.8,.1,3,'synthetic unit test')],(1,120))
        self.assertAlmostEqual(q,7.2e-7)
        self.assertAlmostEqual(rows[0]['heat_W_m2']+rows[0]['other_energy_W_m2'],rows[0]['absorbed_W_m2'])
    def test_historical_band_total_and_units(self):
        self.assertAlmostEqual(sum(b.irradiance_w_m2 for b in historical_solar_bands()),.00464)
    def test_historical_bands_do_not_invent_response(self):
        with self.assertRaises(ValueError):deposited_heat(historical_solar_bands(),(.1,118.))
    def test_gap_energy_without_intraband_interpolation(self):
        r=gap_energy_bounds(historical_solar_bands(),(1239.841984/50,120.181141))
        self.assertAlmostEqual(r['lower_fraction'],.375)
        self.assertAlmostEqual(r['upper_fraction'],3.79/4.64)
    def test_changed_optical_archive_rejected(self):
        with self.assertRaises(ValueError):audit_optical_files('wavelength,irradiance\n202,1\n300,1\n','0.12 1.2 0.1\n0.3 1.3 0.2\n')
    def test_invalid_transmission(self):
        with self.assertRaises(ValueError):
            deposited_heat([Band(1,120,.01,1.1,1,.1,3,'test')],(1,120))


if __name__ == "__main__":
    unittest.main()
