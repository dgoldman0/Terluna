"""Constants recorded by name and value: found from imports, compared with the current table."""
import pytest

from shared import constants
from shared.provenance import constants_changed, constants_read, constants_used


def test_names_come_from_imports_and_attribute_reads(tmp_path):
    first = tmp_path / "first.py"
    first.write_text("from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS as month\n")
    second = tmp_path / "second.py"
    second.write_text("from shared import constants as k\nimport shared.constants\n"
                      "x = k.AU + shared.constants.SUN_RADIUS\n")
    (tmp_path / "table.json").write_text("{}")
    names = constants_read([first, second, tmp_path / "table.json"])
    assert names == ["AU", "MOON_RADIUS", "SUN_RADIUS", "SYNODIC_MONTH_DAYS"]
    assert constants_used([first]) == {"MOON_RADIUS": constants.MOON_RADIUS,
                                       "SYNODIC_MONTH_DAYS": constants.SYNODIC_MONTH_DAYS}


def test_the_whole_table_cannot_be_read_by_name(tmp_path):
    whole = tmp_path / "whole.py"
    whole.write_text("from shared.constants import DATA\n")
    with pytest.raises(ValueError):
        constants_read([whole])


def test_only_a_changed_value_shows():
    recorded = {"MOON_RADIUS": constants.MOON_RADIUS, "AU": constants.AU}
    assert constants_changed(recorded) == {}
    assert constants_changed(None) == {}
    changed = constants_changed({**recorded, "MOON_RADIUS": 1.0, "RETIRED": 3})
    assert changed == {"MOON_RADIUS": (1.0, constants.MOON_RADIUS), "RETIRED": (3, None)}
