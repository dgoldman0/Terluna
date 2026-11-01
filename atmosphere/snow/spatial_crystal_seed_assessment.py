"""Two coarse warm-restart sensitivities to the unmeasured initial tip rim.

Both retain the fitted mass and dimensions exactly in the prescribed cavity
family. This is an identifiability screen, not calibration of the seed shape
or a morphology-convergence test. The 2 micrometre baseline is preserved in
fd5fa37904382a39d7ef82806421027719afd608.
"""
import json
from pathlib import Path
from atmosphere.snow.spatial_crystal_assessment import run, file_id

HERE=Path(__file__).resolve().parent
RESULT=HERE/'results/spatial_crystal_seed_assessment.json'


def assess():
    rows=[]
    for rim in (1e-6,4e-6):
        row=run('earth',300.,'coarse_local',rim)
        rows.append(row)
        r=row['warm_restarts']['cases'][0]
        print(json.dumps(dict(rim_m=rim,complete=r['complete'],error=r['error'],
            ratios=r['final_to_observed_ratios_a_c_mass'])),flush=True)
    return dict(schema='terluna.atmosphere.spatial-crystal-seed-assessment/1',
        evidence='coarse conditional initial-profile sensitivity, not empirical or morphology validation',
        producer=file_id(Path(__file__).resolve()),runs=rows,limits=[
            'Only initial tip rim changes: temperature, pressure, observed initial mass and maximum dimensions, local kinetics, curvature and ventilation are held fixed.',
            'One micrometre cells underresolve the kinetic length and the thinnest tip; the purpose is to screen sensitivity, not rank profiles by physical accuracy.',
            'Both seed profiles are prescribed and unmeasured. Separate observed crystal sheaths and stochastic step sources remain absent.'])


if __name__=='__main__':
    result=assess()
    RESULT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
