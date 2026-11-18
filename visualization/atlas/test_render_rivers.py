"""Tests of the page's named rivers (render.river_marks): each mark's name stands on its course, and its card's lines
and figures say what the waterways record holds."""

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appearance  # noqa: E402
import render  # noqa: E402


def record(path, **over):
    w = dict(anchor='river:1.0,2.0', name='Testmayu', tongue='Runa (Quechua-rooted)', meaning='Test River',
             ends_in='Hanaqucha', ends_in_anchor='sea:Mare Moscoviense', fed_by='headwaters', works=False,
             length_km=1200, lake_km=800, lakes_crossed=5, flow_mouth_m3s=10000, lowest_year_m3s=9000,
             regions=[dict(name='Paraurqu', anchor='region:x', tongue='runa', share=1.0)],
             course=[[0.0, 10.0], [0.0, 20.0], [0.0, 50.0]])
    w.update(over)
    small = dict(w, anchor='river:5.0,6.0', name='Huchuymayu', flow_mouth_m3s=1500, lowest_year_m3s=900,
                 lakes_crossed=0, lake_km=0)
    unnamed = dict(w, anchor='river:7.0,8.0', name=None)
    Path(path).write_text(json.dumps(dict(waterways=[small, w, unnamed])))


def test_marks_are_named_ordered_and_placed_half_way():
    with tempfile.TemporaryDirectory() as d:
        record(Path(d) / 'w.json')
        marks = render.river_marks(Path(d) / 'w.json')
    assert [m['name'] for m in marks] == ['Testmayu', 'Huchuymayu']          # the largest first, the unnamed left out
    m = marks[0]
    assert abs(m['lat']) < 1e-6 and abs(m['lon'] - 30.0) < 1e-3              # half way along 40 degrees of equator
    assert m['links'] == ['sea:Mare Moscoviense', 'region:x']
    assert m['course'][0] == [0.0, 10.0]


def test_card_lines_say_what_the_record_holds():
    with tempfile.TemporaryDirectory() as d:
        record(Path(d) / 'w.json', fed_by='a branch the works cut from Hanaqucha', works=True)
        big, small = render.river_marks(Path(d) / 'w.json')
    rows = dict(big['rows'])
    assert rows['Flow where it ends'] == '10,000\u00a0m³/s' and rows['Lowest year'] == '9,000\u00a0m³/s'
    assert rows['Course'] == 'cut by the works' and rows['Water from'] == 'a branch the works cut from Hanaqucha'
    width = appearance.RIVER_WIDTH * 10000 ** 0.5
    assert rows['Width where it ends'] == f'about {width / 1000:.1f}\u00a0km'
    assert 'Its water is a branch the works cut from Hanaqucha.' in big['text']
    assert '800\u00a0km of it through 5 lakes' in big['text']
    assert 'through' not in small['text'].split('.')[0]                      # no lakes crossed, none named
    assert dict(small['rows'])['Width where it ends'].endswith('\u00a0m')
