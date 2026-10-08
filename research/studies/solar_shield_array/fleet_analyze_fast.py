"""Run the established analyzer with an equivalent vectorized encounter backend."""
import argparse
import json
from pathlib import Path

from shared.provenance import constants_used
from protection.dynamics.fleet_exclusions import swept_conflicts
from . import fleet_analyze
from .fleet_run import identities, digest
from .cycling_search import ROOT, HERE
from .cycling_analysis import compact_series_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--template', default='fleet_dense.json')
    a = p.parse_args()
    # Only this process's numerical backend is replaced. Original producers
    # and previously exported products retain their recorded source identities.
    fleet_analyze.swept_conflicts = swept_conflicts
    cases = {name: fleet_analyze.analyze(name, 21600., 16, 0.) for name in a.names}
    template = HERE/'results'/a.template
    product = json.loads(template.read_text())
    sources = {**identities(), str(Path(__file__).relative_to(ROOT)): digest(__file__),
        'protection/dynamics/fleet_exclusions.py': digest(ROOT/'protection/dynamics/fleet_exclusions.py'),
        'research/studies/solar_shield_array/fleet_analyze.py': digest(HERE/'fleet_analyze.py'),
        'research/studies/solar_shield_array/cycling_analysis.py': digest(HERE/'cycling_analysis.py')}
    product['producer'] = dict(source_hashes=sources, constants=constants_used(sources),
        template_product=a.template, template_product_sha256=digest(template),
        encounter_backend='Vectorized pair minima; tested against the reference swept collision and solar-shadow exclusions')
    product['cases'] = cases
    product['parameters'] = dict(names=a.names, cadence_hours=6., suns=16, rotation=0., output=a.output)
    (HERE/'results'/a.output).write_text(compact_series_json(product))
    print('Wrote '+a.output, flush=True)


if __name__ == '__main__':
    main()
