"""The lane rules allow the intended dependencies and reject the rest."""
import unittest

from shared.check_layers import allowed, lane


class LaneRuleTests(unittest.TestCase):
    def test_lanes(self):
        self.assertEqual(lane("atmosphere/lower_air.py"), "domain")
        self.assertEqual(lane("immersion/engine/core.js"), "immersion")
        self.assertEqual(lane("immersion/bake/columns.mjs"), "immersion-bake")
        self.assertEqual(lane("research/studies/forest_patch/run.py"), "research")

    def test_intended_dependencies(self):
        self.assertTrue(allowed("biosphere/forest_patch.py", "atmosphere/lower_air.py"))
        self.assertTrue(allowed("atmosphere/lower_air.py", "shared/constants.py"))
        self.assertTrue(allowed("immersion/engine/core.js", "shared/constants.js"))
        self.assertTrue(allowed("immersion/bake/columns.mjs", "atmosphere/column/export_columns.cjs"))
        self.assertTrue(allowed("visualization/labs/build_accuracy.py", "illumination/references/x.py"))
        self.assertTrue(allowed("visualization/labs/x.py", "immersion/engine/cloud-renderer.js"))

    def test_rejected_dependencies(self):
        self.assertFalse(allowed("immersion/engine/cloud-renderer.js", "atmosphere/column/weather-column.js"))
        self.assertFalse(allowed("atmosphere/lower_air.py", "immersion/engine/core.js"))
        self.assertFalse(allowed("illumination/sky/solver.py", "visualization/month-of-light/cycle.js"))
        self.assertFalse(allowed("shared/constants.py", "atmosphere/lower_air.py"))
        self.assertFalse(allowed("research/check.py", "archive/training_data/x.py"))


if __name__ == "__main__":
    unittest.main()
