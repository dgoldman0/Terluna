"""The solar cycle of the light below 175 nm, band by band, from FISM2's daily record.

    python -m atmosphere.loss_response.cycle      # writes results/solar_cycle.json in a few seconds

R1 is a long-term average, and the loss response gives each state at quiet Sun and at solar maximum. FISM2's daily
record from February 1947, integrated over each band of the limb tables (atmosphere/middle_atmosphere/
limb_inputs.json), measures each band's own cycle. Every calendar year's mean over the WHI 2008 quiet week is a measured activity: a factor on each
band; Lyman-alpha's factor for the sky's glow, which is solar Lyman-alpha scattered by interplanetary hydrogen; and the
energy-weighted factors on the X-rays below 10 nm and on the ultraviolet from 10 to 175 nm, which the escape model's
count over the disk takes, the energy-limited bound weighting each band itself. Solar cycles 23 and 24 run from the
minimum of August 1996 to that of December 2019, which the calendar years 1997-2019 span; 2020-2025 give cycle 25 as
far as the record is taken. The author accepted the measured ultraviolet for the cycle average on 2026-10-07.

Each solar maximum is the year around the month of its smoothed sunspot maximum (SILSO, version 2), from six months
before it to five after, for cycles 19 to 25; cycle 18's year begins before the record. The loss response's solar
maximum is the strongest of them, by the energy-weighted ultraviolet, taken whole: every band, the X-rays and the glow
from the same year. The author chose the measured maxima, with the record taken back to 1947, on 2026-10-07, in place
of the escape model's 2.5 on the ultraviolet, FISM2's mean rise of the X-rays over the last three maxima and 1.5 on the
glow, which the loss response keeps as its stress case. The strongest 27 days, a solar rotation, and the mean spectrum
of cycles 19 to 24 are kept for reference.

Each year is taken as a steady state. The upper air's slowest layers take years to respond and its highest days, so
the mean of the yearly states' losses leans high against a column that follows the cycle, and the state of the whole
span's mean spectrum, with no swing at all, gives the low end. FISM2 is a model of the measured irradiance (Chamberlin
et al. 2020); before 2002 it rests on proxies of solar activity, and before 1978 on the 10.7 cm radio flux alone, so
the spectra of cycles 19 to 21 are the least certain.
"""
from __future__ import annotations
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from shared.provenance import constants_used
from atmosphere.middle_atmosphere import escape, fetch_inputs, fetch_limb_inputs
from atmosphere.loss_response import traced

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'solar_cycle.json'
SCHEMA = 'terluna.atmosphere.solar-cycle/2'
QUIET_WEEK = (dt.date(2008, 4, 10), dt.date(2008, 4, 16))        # the WHI 2008 quiet-Sun campaign
SPANS = dict(cycles_23_24=(1997, 2019), cycle_25_to_2025=(2020, 2025))
# The month of each cycle's smoothed sunspot maximum: SILSO, version 2 (sidc.be/SILSO/cyclesminmax).
SUNSPOT_MAXIMA = dict(cycle_19=(1958, 3), cycle_20=(1968, 11), cycle_21=(1979, 12), cycle_22=(1989, 11),
                      cycle_23=(2001, 11), cycle_24=(2014, 4), cycle_25=(2024, 10))
# Cycles 19 to 24, from the minimum of April 1954 to that of December 2019 (SILSO).
LONG_MEAN = ('cycles_19_24_mean_spectrum', dt.date(1954, 4, 1), dt.date(2019, 11, 30))
ROTATION_DAYS = 27
PARTS = ('1947-02_1996-07', '1996-08_2025-12')       # the record's files of each band, joined in this order
LYMAN_ALPHA_NM = 121.567
WHI = 'whi2008_ref_solar_irradiance_ver2.dat'
_RECORD = {}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def year_around(year, month):
    """The year around a month: from the first day six months before it to the last day five months after."""
    start = dt.date(year + (month - 7) // 12, (month - 7) % 12 + 1, 1)
    after = month + 6
    return start, dt.date(year + (after - 1) // 12, (after - 1) % 12 + 1, 1) - dt.timedelta(days=1)


MAXIMA = {name: year_around(*month) for name, month in SUNSPOT_MAXIMA.items()}


def band_edges():
    """The limb tables' bands (limb_heat.DESIGN), read from their product."""
    return traced.limb_product()['design']['bands_nm']


def file_name(lo, hi, part=PARTS[-1]):
    return f'fism2_band_{lo:g}-{hi:g}nm_{part}.csv'


def file_names():
    edges = band_edges()
    return [file_name(lo, hi, part) for lo, hi in zip(edges[:-1], edges[1:]) for part in PARTS]


def record():
    """FISM2's daily irradiance integrated over each band (W/m^2): the dates and an array (bands, days), each band's
    files joined and checked to run day by day without a gap."""
    if 'data' not in _RECORD:
        edges = band_edges()
        dates, values = None, []
        for lo, hi in zip(edges[:-1], edges[1:]):
            rows = [r for part in PARTS
                    for r in fetch_limb_inputs.path(file_name(lo, hi, part)).read_text().strip().splitlines()[1:]]
            days = [dt.datetime.strptime(r.split(',')[0], '%Y%j').date() for r in rows]
            if dates is None:
                dates = np.array(days)
                if np.any(np.diff(dates) != dt.timedelta(days=1)):
                    raise ValueError('the record skips or repeats a day')
            elif list(dates) != days:
                raise ValueError(f'band {lo:g}-{hi:g} nm covers other days')
            values.append([float(r.split(',')[1]) for r in rows])
        _RECORD['data'] = dates, np.array(values)
    return _RECORD['data']


def whi_energy():
    """The WHI 2008 quiet Sun's energy in each band (W/m^2), which weights the bands' factors."""
    w, f = escape.whi_quiet_sun()
    edges = band_edges()
    return np.array([float(f[(w >= lo) & (w < hi)].sum() * 0.1) for lo, hi in zip(edges[:-1], edges[1:])])


def factors(start, end):
    """Each band's mean irradiance from start to end (inclusive) over its mean in the quiet week."""
    dates, values = record()
    span = (dates >= start) & (dates <= end)
    quiet = (dates >= QUIET_WEEK[0]) & (dates <= QUIET_WEEK[1])
    if not span.any():
        raise ValueError(f'no record from {start} to {end}')
    return values[:, span].mean(axis=1) / values[:, quiet].mean(axis=1), int(span.sum())


def activity(label, start, end):
    """A measured activity: each band's factor, the glow's (Lyman-alpha's), and the energy-weighted factors on the
    X-rays and on the ultraviolet from 10 to 175 nm."""
    edges = band_edges()
    bands, days = factors(start, end)
    energy = whi_energy()
    xray = np.array([hi <= 10.0 for hi in edges[1:]])
    glow = bands[int(np.searchsorted(edges, LYMAN_ALPHA_NM, side='right') - 1)]
    return dict(label=label, start=start.isoformat(), end=end.isoformat(), days=days, bands=[float(x) for x in bands],
                glow=float(glow), xray=float(bands[xray] @ energy[xray] / energy[xray].sum()),
                uv=float(bands[~xray] @ energy[~xray] / energy[~xray].sum()))


def years(first, last):
    """The measured activity of each calendar year."""
    return [activity(str(y), dt.date(y, 1, 1), dt.date(y, 12, 31)) for y in range(first, last + 1)]


def span_mean(name):
    """The measured activity of a whole span's mean spectrum."""
    first, last = SPANS[name]
    return activity(f'{name}_mean_spectrum', dt.date(first, 1, 1), dt.date(last, 12, 31))


def periods():
    """Every year of each span with its span, and each span's mean spectrum."""
    return {name: dict(years=years(*bounds), mean_spectrum=span_mean(name)) for name, bounds in SPANS.items()}


def maxima():
    """The measured activity of the year around each solar maximum."""
    return {name: activity(name, *bounds) for name, bounds in MAXIMA.items()}


def strongest(measured):
    """The loss response's solar maximum: the year around the maximum with the most ultraviolet, whole."""
    name = max(measured, key=lambda k: measured[k]['uv'])
    return dict(measured[name], label='solar_maximum', cycle=name)


def strongest_rotation():
    """The ROTATION_DAYS in the record with the most ultraviolet, energy-weighted, as a measured activity."""
    dates, values = record()
    edges = band_edges()
    energy = whi_energy()
    xray = np.array([hi <= 10.0 for hi in edges[1:]])
    quiet = (dates >= QUIET_WEEK[0]) & (dates <= QUIET_WEEK[1])
    uv = (values[~xray] / values[~xray][:, quiet].mean(axis=1, keepdims=True)).T @ energy[~xray] / energy[~xray].sum()
    i = int(np.argmax(np.convolve(uv, np.ones(ROTATION_DAYS) / ROTATION_DAYS, 'valid')))
    return activity(f'strongest_{ROTATION_DAYS}_days', dates[i], dates[i + ROTATION_DAYS - 1])


def product():
    """The written solar-cycle product, checked against this schema."""
    data = json.loads(OUT.read_text())
    if data['schema'] != SCHEMA:
        raise ValueError(f"unexpected solar-cycle schema {data['schema']}")
    return data


def main(argv=None) -> int:
    edges = band_edges()
    energy = whi_energy()
    spans = periods()
    measured = maxima()
    files = [HERE / 'cycle.py', HERE / 'traced.py', HERE.parents[1] / 'atmosphere/middle_atmosphere/escape.py']
    code = {str(p.relative_to(HERE.parents[1])): digest(p) for p in files}
    inputs = {name: digest(fetch_limb_inputs.path(name)) for name in file_names()}
    inputs[WHI] = digest(fetch_inputs.path(WHI))
    product = dict(
        schema=SCHEMA,
        # The bands are limb_heat.DESIGN's (bands_nm below); the limb product, which reads this one, is not pinned.
        producer=dict(domain='atmosphere', files=code, constants=constants_used(files), inputs=inputs),
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=('Factors are irradiance over the mean of the WHI 2008 quiet week (10-16 April 2008), per band of '
                      'bands_nm. An activity\'s "bands" weight the limb tables\' band heats and the photon fluxes of the '
                      'exosphere step; "glow" scales the sky\'s Lyman-alpha glow; "xray" and "uv" are energy-weighted '
                      'over the bands below and above 10 nm. A span\'s years are steady states whose mean loss leans '
                      'high; its mean spectrum gives the low end. maxima are the years around each cycle\'s sunspot '
                      'maximum; solar_maximum, the strongest of them whole, is the loss response\'s solar maximum. '
                      'convention holds the escape model\'s factors, the loss response\'s stress case.'),
        bands_nm=list(edges), whi_quiet_energy_W_m2=[float(x) for x in energy], quiet_week=[d.isoformat() for d in QUIET_WEEK],
        convention=dict(solar_maximum_ultraviolet=escape.SOLAR_MAXIMUM, solar_maximum_xray=escape.XRAY_SOLAR_MAXIMUM),
        sunspot_maxima={name: f'{y}-{m:02d}' for name, (y, m) in SUNSPOT_MAXIMA.items()},
        maxima=measured, solar_maximum=strongest(measured), strongest_rotation=strongest_rotation(),
        long_mean=activity(*LONG_MEAN), spans=spans)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    for name, span in spans.items():
        m = span['mean_spectrum']
        print(f"{name}: mean spectrum ultraviolet {m['uv']:.3f}, X-rays {m['xray']:.2f}, glow {m['glow']:.3f}")
    for name, m in dict(measured, long_mean=product['long_mean'], rotation=product['strongest_rotation']).items():
        print(f"{name}: {m['start']} to {m['end']}, ultraviolet {m['uv']:.3f}, X-rays {m['xray']:.2f}, "
              f"glow {m['glow']:.3f}")
    print(f"solar maximum: {product['solar_maximum']['cycle']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
