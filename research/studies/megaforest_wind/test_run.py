"""Scope guard for the megaforest wind study runner."""
import unittest
from research.studies.megaforest_wind.run import run, ROOT


class ScopeTests(unittest.TestCase):
    def test_immersion_output_rejected(self):
        with self.assertRaises(ValueError): run(ROOT/'immersion'/'forbidden-wind-output')
        self.assertFalse((ROOT/'immersion'/'forbidden-wind-output').exists())


if __name__ == "__main__":
    unittest.main()
