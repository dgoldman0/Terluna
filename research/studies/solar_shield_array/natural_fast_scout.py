"""Replay the scout workflow with an explicitly labelled cheaper force model."""
from . import natural_scout as runner
from protection.dynamics.service_reduced import reduced_path

original_setup=runner.setup
original_write=runner.write


def setup(extra):
    return original_setup([*extra,'protection/dynamics/service_reduced.py',
        'research/studies/solar_shield_array/natural_fast_scout.py'])


def write(name,data,before,cpu):
    data['scope']='All 361 members in a centre-gravity/centre-Sun model, with body-visible fraction and no mutual shadows. Contact and service tests are proposal diagnostics. Post-contact mathematical states and unexecuted Sun-facing turns earn no physical service credit.'
    original_write('natural_fast_scout.json',data,before,cpu)


if __name__=='__main__':
    runner.setup=setup;runner.write=write;runner.independent_path=reduced_path
    # Keep the interrupted finite-area scout raw states independently.
    original_save=runner.save
    runner.save=lambda name,**arrays:original_save('fast_'+name,**arrays)
    runner.main()
