"""The lane rules allow the intended dependencies and reject the rest."""
import unittest
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

from shared.check_layers import allowed, lane, tracked_sources


class LaneRuleTests(unittest.TestCase):
    def test_source_inventory_handles_new_deleted_and_space_named_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "shared").mkdir()
            (root / "shared/new file.py").write_text("pass\n")
            output = "shared/new file.py\0shared/deleted.py\0archive/old.py\0"
            with patch("shared.check_layers.ROOT", root), patch("shared.check_layers.subprocess.run") as run:
                run.return_value = SimpleNamespace(stdout=output)
                self.assertEqual(tracked_sources(), ["shared/new file.py"])
                self.assertIn("--others", run.call_args.args[0])

    def test_lanes(self):
        self.assertEqual(lane("atmosphere/lower_air.py"), "domain")
        self.assertEqual(lane("immersion/engine/core.js"), "immersion-engine")
        self.assertEqual(lane("immersion/world/cove.js"), "immersion-world")
        self.assertEqual(lane("immersion/experiences/shoreline/app.js"), "immersion")
        self.assertEqual(lane("immersion/bake/columns.mjs"), "immersion-bake")
        self.assertEqual(lane("research/studies/forest_patch/run.py"), "research")
        self.assertEqual(lane("resources/ledger.py"), "domain")
        self.assertEqual(lane("provisioning/needs.py"), "domain")

    def test_intended_dependencies(self):
        self.assertTrue(allowed("biosphere/forest_patch.py", "atmosphere/lower_air.py"))
        self.assertTrue(allowed("atmosphere/lower_air.py", "shared/constants.py"))
        self.assertTrue(allowed("resources/ledger.py", "geography/water_inventory.py"))
        self.assertTrue(allowed("provisioning/needs.py", "biosphere/canopy/model.py"))
        self.assertTrue(allowed("immersion/engine/core.js", "shared/constants.js"))
        self.assertTrue(allowed("immersion/bake/columns.mjs", "atmosphere/column/export_columns.cjs"))
        self.assertTrue(allowed("visualization/labs/build_accuracy.py", "illumination/references/x.py"))
        self.assertTrue(allowed("visualization/labs/x.py", "immersion/engine/cloud-renderer.js"))
        self.assertTrue(allowed("immersion/world/cove-weather.js", "immersion/engine/core.js"))
        self.assertTrue(allowed("immersion/experiences/shoreline/main.js", "immersion/world/cove.js"))

    def test_rejected_dependencies(self):
        self.assertFalse(allowed("immersion/engine/cloud-renderer.js", "atmosphere/column/weather-column.js"))
        self.assertFalse(allowed("atmosphere/lower_air.py", "immersion/engine/core.js"))
        self.assertFalse(allowed("illumination/sky/solver.py", "visualization/month-of-light/cycle.js"))
        self.assertFalse(allowed("shared/constants.py", "atmosphere/lower_air.py"))
        self.assertFalse(allowed("research/check.py", "archive/training_data/x.py"))
        self.assertFalse(allowed("immersion/engine/scene.js", "immersion/world/cove.js"))
        self.assertFalse(allowed("archive/renderer-lab/main.js", "immersion/engine/om.js"))
        self.assertFalse(allowed("immersion/engine/om.js", "archive/renderer-lab/main.js"))


if __name__ == "__main__":
    unittest.main()
