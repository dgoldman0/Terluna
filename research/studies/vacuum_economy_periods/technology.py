"""The technology of the vacuum period, 2041-2241: each quantity of the stack as a path in a slow, a central and a fast
case, the static case that holds 2041, and the technologies that arrive later, read from period_technology.json."""
from __future__ import annotations
from functools import lru_cache
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = ('static', 'slow', 'central', 'fast')


@lru_cache(maxsize=None)
def roadmap():
    return json.loads((HERE / 'period_technology.json').read_text())


@lru_cache(maxsize=None)
def value(name, case, year):
    """One quantity in a calendar year: the 2041 value in the static case, otherwise interpolated between the stated
    years, geometrically where the path says so, and held after the last."""
    t = roadmap()
    q = t['quantities'][name]
    years = t['years']
    values = q['slow' if case == 'static' else case]
    if case == 'static' or year <= years[0]:
        return float(values[0])
    if year >= years[-1]:
        return float(values[-1])
    i = next(k for k in range(1, len(years)) if year <= years[k])
    f = (year - years[i - 1]) / (years[i] - years[i - 1])
    a, b = values[i - 1], values[i]
    if q.get('interpolate') == 'geometric' and a > 0 and b > 0:
        return float(a * (b / a) ** f)
    return float(a + (b - a) * f)


def arrival(name, case):
    """The first year a technology is available in a case, or None where it does not arrive (always in the static case)."""
    if case == 'static':
        return None
    return roadmap()['arrivals'][name]['year'][case]


def arrived(name, case, year):
    first = arrival(name, case)
    return first is not None and year >= first


def parameter(name, key, case):
    """A parameter an arrival brings: a single value, or a value by case."""
    v = roadmap()['arrivals'][name][key]
    return v[case] if isinstance(v, dict) else v


class Year:
    """The technology of one calendar year in one case, read by attribute: `Year('central', 2100).array_W_kg`."""

    def __init__(self, case, year):
        self.case, self.year = case, year

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        return value(name, self.case, self.year)

    def arrived(self, name):
        return arrived(name, self.case, self.year)
