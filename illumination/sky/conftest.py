"""pytest wiring for the sky solver, whose modules import each other by file name.

The solver needs numba. Without it these tests are not collected, and the pytest
header says so, instead of failing the repository-wide run.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
HAVE_NUMBA = importlib.util.find_spec("numba") is not None

if HAVE_NUMBA:
    sys.path.insert(0, str(HERE))
else:
    collect_ignore_glob = ["test_*.py"]


def pytest_report_header(config):
    if not HAVE_NUMBA:
        return "illumination/sky: numba not installed; sky solver tests not collected"
