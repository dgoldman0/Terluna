"""The lab's modules import each other by name and the domains from the repository root."""
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
for path in (str(LAB), str(LAB.parents[1])):
    if path not in sys.path:
        sys.path.insert(0, path)
