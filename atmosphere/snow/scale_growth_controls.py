"""Isolate anisotropic ventilation in the pinned SCALE ice growth routine.

Mode 0 retains every native physics expression. Mode 1 sets only the ratio of
polar to equatorial ventilation factors to one; mean mass-transfer ventilation
remains native. Both modes record the first live ice particle at every growth
call. This is a mechanism diagnostic, not an empirically calibrated repair.
"""
from pathlib import Path
import argparse
import difflib
import hashlib
import json
import shutil

from atmosphere.snow import scale_lunar as lunar, scale_build

CACHE = lunar.CACHE
BUILD_NAME = 'growth_controls'
SOURCE = 'contrib/SDM/sdm_subldep.f90'
SOURCE_SHA = '37afcfdb871cdc3da53488810ada6e8541128eccd3a0e774faf6e359408cb444'
MODES = {0: 'native anisotropic ventilation',
         1: 'unit polar/equatorial ventilation ratio; native mean ventilation'}
TRACE_FIELDS = (
    'call', 'mode', 'particle_index', 'elapsed_before_s', 'elapsed_after_s',
    'dt_s', 'temperature_K', 'pressure_Pa', 'air_density_kg_m3',
    'ice_supersaturation', 'intrinsic_growth_ratio', 'capacitance_m',
    'mean_ventilation', 'native_axis_ventilation_ratio', 'axis_ventilation_ratio',
    'effective_growth_ratio', 'deposition_density_kg_m3', 'equatorial_radius_m',
    'polar_radius_m', 'mass_kg', 'density_kg_m3', 'terminal_speed_m_s',
    'reynolds', 'schmidt', 'ventilation_x', 'diffusive_resistance',
    'thermal_resistance', 'new_equatorial_radius_m', 'new_polar_radius_m',
    'new_mass_kg', 'new_density_kg_m3', 'mass_increment_kg', 'volume_increment_m3',
)
STATE = '''
  ! Local mechanism diagnostic; the native path remains mode zero.
  logical, save :: gc_initialized = .false.
  integer, save :: gc_mode = 0, gc_trace_unit
  integer(kind=8), save :: gc_calls = 0
  real(RP), save :: gc_elapsed = 0.0_RP
'''


def setup_routine():
    header = '# ' + ' '.join(TRACE_FIELDS)
    return f'''
  subroutine gc_setup()
    logical :: exists
    integer :: config_unit, ios
    if(gc_initialized)return
    inquire(file='growth_control.txt',exist=exists)
    if(exists)then
      open(newunit=config_unit,file='growth_control.txt',status='old',action='read',iostat=ios)
      if(ios/=0)error stop 'Cannot open growth_control.txt'
      read(config_unit,*,iostat=ios)gc_mode
      close(config_unit)
      if(ios/=0)error stop 'Cannot read growth control mode'
    endif
    if(gc_mode/=0 .and. gc_mode/=1)error stop 'Unknown growth control mode'
    open(newunit=gc_trace_unit,file='growth_trace.txt',status='new',action='write',iostat=ios)
    if(ios/=0)error stop 'Growth trace already exists or cannot be written'
    write(gc_trace_unit,'(A)') '{header}'
    gc_initialized=.true.
  end subroutine gc_setup
'''


def patch(text):
    """Reject unidentified source; modify only the explicitly selected ratio."""
    if hashlib.sha256(text.encode()).hexdigest() != SOURCE_SHA:
        raise ValueError('Unidentified native ice growth source')
    replacements = [
        ('\ncontains\n', '\n'+STATE+'\ncontains\n'+setup_routine()+'\n'),
        ('    real(RP) :: fratio    ! ventilation effect ratio fc/fa\n',
         '    real(RP) :: fratio    ! ventilation effect ratio fc/fa\n'
         '    real(RP) :: gc_native_fratio\n'),
        ('    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)\n',
         '    call gc_setup()\n'
         '    gc_calls=gc_calls+1\n'
         '    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)\n'),
        ('       growth_ratio_star = growth_ratio * fratio\n',
         '       gc_native_fratio=fratio\n'
         '       if(gc_mode==1)fratio=1.0_RP\n'
         '       growth_ratio_star = growth_ratio * fratio\n'),
        ('       !### update the super-droplet ###!\n',
         '''       !### update the super-droplet ###!
       ! Record the first live ice particle. These are the inputs actually used
       ! after prescribed air and native terminal speed have been evaluated.
       ! The new state includes native geometric limiters, before the driver's
       ! separate shared-water ceiling. That ceiling is checked independently.
       if(m==1)then
         write(gc_trace_unit,'(3(I20,1X),30(ES25.16,1X))') &
              gc_calls,gc_mode,n,gc_elapsed,gc_elapsed+real(sdm_dtsbl,RP), &
              real(sdm_dtsbl,RP),t_sd,p_sd,rhom_sd,ssi,growth_ratio,capaci,barf, &
              gc_native_fratio,fratio,growth_ratio_star,rho_dep,sd_ra,sd_rc, &
              sd_mass,sd_rho,sd_tvel,n_re,n_sc,n_ventX, &
              Fac_ice_dd_invrhoi,Fac_ice_kk_invrhoi,new_sd_ra,new_sd_rc, &
              new_sd_mass,new_sd_rho,delta_mass,delta_vol
         flush(gc_trace_unit)
       endif
'''),
        ('    return\n  end subroutine sdm_subldep\n',
         '    gc_elapsed=gc_elapsed+real(sdm_dtsbl,RP)\n'
         '    return\n  end subroutine sdm_subldep\n'),
    ]
    for old, new in replacements:
        if text.count(old) != 1:
            raise ValueError('Ambiguous ice growth insertion: '+old.strip())
        text = text.replace(old, new, 1)
    return text


def check_build(cache=CACHE):
    """Require the existing source, executables and frozen producer to match."""
    target=Path(cache)/BUILD_NAME
    record=json.loads((target/'manifest.json').read_text())
    producer=lunar.base.digest(Path(__file__))
    if record['producer_sha256'] != producer:
        raise ValueError('Growth control build has a different producer')
    if lunar.base.digest(target/'producer.py') != producer:
        raise ValueError('Frozen growth control producer changed')
    if record['original_sha256'] != SOURCE_SHA:
        raise ValueError('Growth control build used an unidentified native source')
    if lunar.base.digest(target/'work'/SOURCE) != record['patched_sha256']:
        raise ValueError('Patched growth source changed after build')
    for name, digest in record['build']['executable_sha256'].items():
        if lunar.base.digest(target/'work'/scale_build.CASE_REL/name) != digest:
            raise ValueError('Growth control executable changed: '+name)
    return record


def build(cache=CACHE):
    cache=Path(cache).resolve(); target=cache/BUILD_NAME
    if target.exists():
        return check_build(cache)
    # Validate before copying or creating a new build directory.
    original=(cache/'lunar/work'/SOURCE).read_text()
    patched=patch(original)
    target.mkdir()
    work=target/'work'
    shutil.copytree(cache/'lunar/work',work,symlinks=True)
    source=work/SOURCE; source.write_text(patched)
    diff=''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),
                                   fromfile='a/'+SOURCE,tofile='b/'+SOURCE))
    (target/'growth.patch').write_text(diff)
    (target/'producer.py').write_bytes(Path(__file__).read_bytes())
    deps=cache/'deps/usr'
    result=scale_build.build(work,'-I'+str(deps/'include'),
        str(deps/'lib/x86_64-linux-gnu/libnetcdff.so.7')+
        ' /usr/lib/x86_64-linux-gnu/libnetcdf.so.19 -Wl,-rpath,'+
        str(deps/'lib/x86_64-linux-gnu')+' -lm -lz',workers=1)
    record=dict(build=result,
        original_manifest=json.loads((cache/'lunar/manifest.json').read_text()),
        hook_sha256=lunar.base.digest(lunar.HOOK),
        producer_sha256=lunar.base.digest(Path(__file__)),
        source=SOURCE,original_sha256=SOURCE_SHA,patched_sha256=lunar.base.digest(source),
        patch_sha256=lunar.base.digest(target/'growth.patch'),
        modes=MODES,trace_fields=TRACE_FIELDS,
        evidence='Controlled anisotropic-ventilation sensitivity. Mode 0 is checked against the preserved native executable; neither mode is an empirical repair.')
    lunar.boxes.save(target/'manifest.json',record)
    return check_build(cache)


def attach(directory, mode, cache=CACHE):
    """Select the control build for an unexecuted isolated-crystal experiment."""
    if type(mode) is not int or mode not in MODES:
        raise ValueError('Growth mode must be zero or one')
    directory=Path(directory)
    if (directory/'execution.json').exists() or (directory/'growth_trace.txt').exists():
        raise ValueError('Preserve the existing execution; configure a fresh run')
    record=check_build(cache)
    settings=json.loads((directory/'settings.json').read_text())
    if settings.get('collisions') or settings.get('movement',True):
        raise ValueError('This control is restricted to isolated stationary crystals')
    (directory/'growth_control.txt').write_text(str(mode)+'\n')
    settings.update(build=record,build_directory=BUILD_NAME,growth_control_mode=mode,
        growth_control_description=MODES[mode],
        growth_control_producer_sha256=lunar.base.digest(Path(__file__)),
        growth_trace_scope='First live ice particle per native growth call, after air and fall-speed evaluation and before the driver water ceiling; compare that ceiling independently.')
    lunar.boxes.save(directory/'settings.json',settings)
    inputs=json.loads((directory/'inputs.json').read_text())
    # Preserve the original explicit input set and add only the selector.
    names=sorted(set(inputs)|{'growth_control.txt'})
    lunar.boxes.save(directory/'inputs.json',
                     {name:lunar.base.digest(directory/name) for name in names})
    return directory


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['build'])
    args=parser.parse_args()
    print(json.dumps(build(),indent=2))
