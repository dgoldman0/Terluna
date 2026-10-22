"""Build a passive collision-count observer of the pinned SCALE source.

The patch inserts diagnostics around the existing multiplicity cap. It does not
change candidate pairing, rates, random draws, particle updates or time steps.
Numerical parity with the original executable must be checked after building.
Upstream BSD terms are retained in references/scale_bsd_2_clause.txt.
"""
from pathlib import Path
import argparse
import difflib
import json
import shutil

from atmosphere.snow import scale_build

SOURCE_SHA256='d6b9eed09771d1139e1533221d3c586dc10d9ac70f8c418a40429bf6c1d94ad1'
SOURCE_REL='contrib/SDM/sdm_coalescence_cold.f90'


def patch(text):
    import hashlib
    if hashlib.sha256(text.encode()).hexdigest()!=SOURCE_SHA256:
        raise ValueError('Unidentified collision source')
    changes=[
        ('    integer :: skip_flag ! flag to save whether the aggregation is canceled or not',
         '''    integer :: observer_unit
    integer, save :: observer_calls=0
    real(RP) :: observer_attempted, observer_deficit, observer_accepted
    real(RP) :: observer_pmax, observer_ratio_max
    integer :: observer_clipped
''',False),
        ('    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)',
         '''    observer_attempted=0.0_RP
    observer_deficit=0.0_RP
    observer_pmax=0.0_RP
    observer_ratio_max=0.0_RP
    observer_clipped=0
''',False),
        ('          !### set coalescence count ###!',
         '''          observer_pmax=max(observer_pmax,c_rate(tc))
          observer_ratio_max=max(observer_ratio_max,c_rate(tc)/ &
               real(max(sd_n(icptc),sd_n(icptp))/min(sd_n(icptc),sd_n(icptp)),kind=RP))
''',False),
        ('          sd_ncol = min( sd_ncol, int(sd_n1/sd_n2,kind=RP) )',
         '''          observer_attempted=observer_attempted+real(sd_ncol,kind=RP)*real(sd_n2,kind=RP)
          observer_deficit=observer_deficit+real(max(0,sd_ncol-int(sd_n1/sd_n2)),kind=RP)*real(sd_n2,kind=RP)
          if(sd_ncol>int(sd_n1/sd_n2)) observer_clipped=observer_clipped+1
''',False),
        ('    deallocate( fsort_freq )',
         '''
    observer_calls=observer_calls+1
    observer_accepted=observer_attempted-observer_deficit
    if(observer_calls==1)then
       open(newunit=observer_unit,file='collision_observer.txt',status='new',action='write')
       write(observer_unit,'(A)') '# call dt attempted accepted deficit capped_pairs max_probability max_probability_over_capacity'
    else
       open(newunit=observer_unit,file='collision_observer.txt',status='old',position='append',action='write')
    end if
    write(observer_unit,'(I10,1X,4(ES25.16,1X),I10,1X,2(ES25.16,1X))') &
         observer_calls,sdm_dtcol,observer_attempted,observer_accepted,observer_deficit, &
         observer_clipped,observer_pmax,observer_ratio_max
    close(observer_unit)
''',True),
    ]
    for anchor,addition,after in changes:
        if text.count(anchor)!=1:
            raise ValueError(f'Ambiguous observer insertion: {anchor}')
        text=text.replace(anchor,anchor+addition if after else addition+anchor)
    return text


def build(cache):
    cache=Path(cache).resolve()
    root=cache/'observer'
    root.mkdir(exist_ok=False)
    work=root/'work'
    shutil.copytree(cache/'work',work,symlinks=True)
    source=work/SOURCE_REL
    original=source.read_text();modified=patch(original)
    source.write_text(modified)
    diff=''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True),
                                     fromfile='a/'+SOURCE_REL,tofile='b/'+SOURCE_REL))
    (root/'observer.patch').write_text(diff)
    record=dict(original_sha256=SOURCE_SHA256,patched_sha256=scale_build.sha256(source),
                patch_sha256=scale_build.sha256(root/'observer.patch'),source=SOURCE_REL,
                meaning='Attempted and omitted physical collision counts under the original random draws. These do not predict a counterfactual corrected trajectory.',
                instrumentation='Read-only scalar counters, one text record per native collision call; one MPI rank only.')
    deps=cache/'deps/usr'
    result=scale_build.build(work,'-I'+str(deps/'include'),
        str(deps/'lib/x86_64-linux-gnu/libnetcdff.so.7')+' /usr/lib/x86_64-linux-gnu/libnetcdf.so.19 '+
        '-Wl,-rpath,'+str(deps/'lib/x86_64-linux-gnu')+' -lm -lz',workers=1)
    record['build']=result
    (root/'observer.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',required=True,type=Path)
    args=p.parse_args()
    print(json.dumps(build(args.cache),indent=2))
