"""The network vocabulary validates, enforces the passive-failure rule, and the ledger reproduces the report."""
from __future__ import annotations
import unittest

from engineering.network.model import Link, Network, Node
from engineering.network import ledger
from protection.modules import catalogue


class NetworkValidation(unittest.TestCase):
    def network(self, outcome='misses_planets'):
        nodes = [Node('mine', 'source_extraction', 'icy body'), Node('catcher', 'cislunar_capture', 'cislunar space')]
        links = [Link('mine', 'catcher', 'nitrogen_feedstock',
                      dict(packet_kg=1e7, passive_failure_outcome=outcome))]
        return Network(nodes, links)

    def test_consistent_network_validates(self):
        self.network().validate()

    def test_passive_failure_must_miss_planets(self):
        with self.assertRaises(ValueError):
            self.network('earth_impact').validate()

    def test_unknown_names_are_rejected(self):
        with self.assertRaises(ValueError):
            Network([Node('a', 'spaceport', 'x')], []).validate()
        with self.assertRaises(ValueError):
            Network([Node('a', 'factory', 'x')], [Link('a', 'b', 'propellant')]).validate()
        with self.assertRaises(ValueError):
            Network([Node('a', 'factory', 'x'), Node('b', 'factory', 'y')],
                    [Link('a', 'b', 'propellant', dict(speed=1.0))]).validate()


class Ledger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = ledger.ledger(catalogue.build())
        cls.century = [s for s in cls.result['scenarios'] if s['life_years'] == 100.0 and s['recovery'] == 0.999][0]

    def test_every_catalogue_material_has_a_commodity(self):
        for m in catalogue.build()['modules']:
            for material in m['materials_kg']:
                self.assertIn(ledger.MATERIAL_COMMODITY[material], ledger.COMMODITIES)

    def test_reproduces_the_report_replacement_flows(self):
        # protection/report.md, section 15: base filter 1,547 kg/s and all dry held screen hardware
        # 1,883 kg/s at 100-year replacement; magnets 95,000-951,000 kg/s at their planning allowance.
        by_module = self.century['gross_by_module_kg_s']
        self.assertAlmostEqual(by_module['optical_cell'] / 1547.0, 1.0, places=3)
        self.assertAlmostEqual((by_module['optical_cell'] + by_module['holding_pack']) / 1883.0, 1.0, places=3)
        low, high = self.century['magnet_planning_allowance_gross_kg_s']
        self.assertAlmostEqual(low / 95_000, 1.0, places=2)
        self.assertAlmostEqual(high / 951_000, 1.0, places=2)

    def test_reproduces_the_report_propellant_over_a_billion_years(self):
        # Section 15: 278,000 kg/s for 1 Gyr is 8.79e21 kg.
        self.assertAlmostEqual(self.century['over_service_kg']['propellant'] / 8.79e21, 1.0, places=2)


if __name__ == '__main__':
    unittest.main()
