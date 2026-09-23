"""Scope guard for the forest-patch study runner."""
import unittest
from research.studies.forest_patch.run import safe_output, ROOT


class ScopeTests(unittest.TestCase):
    def test_immersion_output_guard(self):
        for p in (ROOT/'immersion',ROOT/'immersion'/'x',ROOT/'unknown'/'..'/'immersion'/'x'):
            with self.assertRaises(ValueError):safe_output(p)


if __name__ == "__main__":
    unittest.main()
