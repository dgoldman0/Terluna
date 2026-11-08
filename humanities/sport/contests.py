"""The sport study's four contests, run together: two over a ball (a pass against an interceptor, a lobbed ball
contested) and two between bodies (a chase, shoving for position), each on the Moon beside Earth and across the
floor's grips (state.GRIPS).

    python -m humanities.sport.contests     # writes humanities/sport/results/contests.json (about 17 minutes)

The pass contest takes nearly all of that time, on three worker processes of under 200 MB each (the reach engine's
jump search, per place, pass and moment); the lob takes about 20 s, the chase about 7 s and the shove a moment.

Each contest module states its setting, its measures and its own assumptions; the shared ones are in state.py, the
players' reach in reach.py and the ball's flight and holds in windows.py. The product keeps each contest's summaries;
the per-place maps behind them are recomputed by the modules on demand.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from humanities.dance import envelope as dance
from humanities.sport import contest_chase, contest_lob, contest_pass, contest_shove
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results' / 'contests.json'
SCHEMA = 'terluna.humanities.sport-contests/1'
FILES = ('humanities/sport/contests.py', 'humanities/sport/contest_pass.py', 'humanities/sport/contest_lob.py',
         'humanities/sport/contest_chase.py', 'humanities/sport/contest_shove.py', 'humanities/sport/reach.py',
         'humanities/sport/windows.py', 'humanities/sport/state.py', 'humanities/sport/vocabulary.py',
         'humanities/sport/envelope.py',
         'humanities/dance/envelope.py', 'humanities/dance/body/model.py', 'humanities/dance/body/pose.py',
         'habitation/low_gravity/movement.py')
INPUTS = (dance.WINDS, dance.POSES)


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def results():
    return dict(
        schema=SCHEMA,
        producer=dict(domain='humanities', files={f: digest(f) for f in FILES},
                      constants=constants_used([ROOT / f for f in FILES]),
                      inputs={str(p.relative_to(ROOT)): digest(p.relative_to(ROOT)) for p in INPUTS}),
        evidence='Screens of four contests between modelled players: each player a centre of mass with the reach the '
                 'envelope\'s closed forms allow (reach.py) and balls flying in the design air (windows.py), for two '
                 'representative bodies at the shared gravity and Earth\'s, swept over three floor grips. They compare '
                 'how contests are decided; no contest has been played, solved on the whole body or measured.',
        pass_against_interceptor=contest_pass.screen(),
        contested_lob=contest_lob.screen(),
        chase=contest_chase.screen(),
        shoving_for_position=contest_shove.screen())


def main() -> int:
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
