"""Bind bounded experiment products, preserved failures, requirements and checks."""
import hashlib
import json
from pathlib import Path
import subprocess
from .coupled_common import HERE,ROOT,RUN,digest


def main():
    names=['benchmark','optimization_attempt','refinement_attempt','optimization',
           'validation_attempt','validation','cycle_probe','restart','audit']
    rows=[];data={}
    for name in names:
        path=HERE/'results'/('coupled_'+name+'.json');p=json.loads(path.read_text());data[name]=p
        rows.append(dict(path=str(path.relative_to(ROOT)),sha256=digest(path),completed=p['completed'],resources=p['resources']))
    observed=sum(p['resources']['cpu_s'] for p in data.values())
    allowances=dict(interrupted_optimization_s=90.,failed_validation_export_s=60.,
        failed_finalize_cpu_s=data['validation']['recovery']['additional_finalize_attempt_cpu_s'],ancillary_compilation_and_inspection_s=20.)
    charged=observed+sum(allowances.values())
    if charged>1800:raise ValueError('Numerical allocation exceeded')
    requirements=[]
    for entry in json.loads((HERE/'natural_checks.json').read_text())['unchanged_pinned_requirements']:
        path=entry['path'];original=subprocess.check_output(['git','show','d8de3e8438832bd51d916d32a9adef5c6d9d1dc4:'+path],cwd=ROOT)
        sha=digest(ROOT/path)
        if sha!=hashlib.sha256(original).hexdigest():raise ValueError('Requirement changed: '+path)
        requirements.append(dict(path=path,sha256=sha,unchanged_since_d8de3e8=True))
    out=dict(schema='terluna.research.coupled-checks/1',products=rows,
        early_checkpoint='d55082e05d1622d82f1cb4121cf73a5f8f60854d',
        budget=dict(cpu_limit_s=1800,observed_producer_cpu_s=observed,allowances=allowances,charged_cpu_s=charged,
            peak_rss_MiB=max(p['resources']['max_rss_MiB'] for p in data.values()),address_space_limit_MiB=2048,
            new_raw_MiB=sum(p.stat().st_size for p in RUN.rglob('*') if p.is_file())/1024**2,new_raw_limit_MiB=512,
            numerical_threads=1,checks_separate_allowance_s_per_commit=300),
        unchanged_pinned_requirements=requirements,accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        checks=dict(early_make_check='745 passed, 55 skipped, same 13 baseline failures; later Makefile targets not reached',
            final_make_check=None,targeted=None),
        remaining='No accepted local improvement, return, next service, recurrence, hardware supply, handover or fleet power. Larger optimization/learning campaign paused.')
    (HERE/'coupled_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['budget'],indent=2))


if __name__=='__main__':main()
