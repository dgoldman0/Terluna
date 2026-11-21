"""Identity, missingness and non-duplicating spectral-distribution controls."""
from collections import Counter
import gzip
import json
import math
from pathlib import Path

import pytest
from resources.source_access import composition as c


def observation(number=381677, designation='2009 BJ81', tax='C;X', filename='a381677.sp250.txt'):
    fields=[(' '*380)]
    line=list(fields[0])
    for start,length,value in ((0,32,filename),(33,6,'' if number is None else str(number)),
                               (40,10,designation),(123,6,tax)):
        line[start:start+length]=f'{value:<{length}}'
    return ''.join(line)


def source_bytes(*lines):
    return ('Byte-by-byte Description of file: ajac532ft7_mrt.txt\n124-129 A6 Taxon\n'+'\n'.join(lines)+'\n').encode()


def row(identity='1',number=1,designation='2001 AA',diameter='1',tax=None):
    return {'spkid':identity,'number':number,'designation':designation,
            'full_name':f'{number} ({designation})','a':'1.2','e':'.2','i':'10',
            'diameter':diameter,'albedo':'.05','spec_B':tax,'spec_T':None}


def entry(tax):
    return {'source':'mithneos','taxonomy':tax}


def test_fixed_width_taxonomy_and_alias_identity_without_filename_doublecount():
    raw=source_bytes(observation(),observation(filename='au2009bj81.sp252n1.txt',tax='Xk:'))
    observations=c.parse_table(raw)
    assert observations[0]=={'spectrum':'a381677.sp250.txt','number':381677,'designation':'2009 BJ81','taxonomy':'C;X'}
    matched,unmatched=c.match_observations([row('20381677',381677,'2009 BJ81')],observations)
    assert not unmatched and list(matched)==['20381677']
    assert len(matched['20381677'])==2
    assert all(o['matched_by']=='both' for o in matched['20381677'])
    assert c.object_state([entry(o['taxonomy']) for o in observations])['family']=='conflicting_families'


def test_provisional_to_numbered_match_and_unmatched_remain_explicit():
    rows=[row('123',42,'2001 AB12')]
    observations=c.parse_table(source_bytes(observation(None,'2001 ab12',filename='au2001ab12.sp1.txt'),
                                            observation(999,'2014 LW21',filename='a999.sp1.txt')))
    matched,unknown=c.match_observations(rows,observations)
    assert len(matched)==1 and matched['123'][0]['matched_by']=='designation'
    assert unknown[0]['number']==999
    with pytest.raises(ValueError,match='different JPL'):
        c.match_observations([row('1',1,'2001 AA'),row('2',2,'2001 AB')],
                            [{'number':1,'designation':'2001 AB','taxonomy':'S'}])


def test_malformed_duplicate_table_and_jpl_identities_are_rejected():
    with pytest.raises(ValueError,match='duplicate'):
        c.parse_table(source_bytes(observation(),observation()))
    with pytest.raises(ValueError,match='header'):
        c.parse_table(b'<html>error</html>')
    with pytest.raises(ValueError,match='Ambiguous JPL'):
        c.match_observations([row('1',1),row('2',1)],[])


def test_taxonomy_preserves_ambiguity_uncertainty_and_low_quality_exclusion():
    assert c.taxonomy('C;X')['family']=='C_or_X'
    assert c.taxonomy('Xk:')=={'raw':'Xk:','flag':'uncertain','usable':True,'family':'X'}
    for tax in ('Sq::','::','U',''):
        assert c.taxonomy(tax)['usable'] is False
    assert c.taxonomy('S;Sr')['family']=='S_Q'
    assert c.taxonomy('Sq;Q')['family']=='S_Q'
    assert c.taxonomy('SMU','jpl_spec_T')['family'] is None
    assert c.taxonomy('new_class')['family'] is None
    assert c.object_state([entry('S'),entry('Q')])['family']=='S_Q'
    assert c.object_state([entry('C;X'),entry('B')])['family']=='conflicting_families'


def test_repeats_and_cross_catalogue_overlap_do_not_amplify_object_or_mass_weights():
    rows=[row('1',tax='S'),row('2',2,diameter='2'),row('3',3,diameter=None)]
    matched={'1':[{'taxonomy':'S'}]*4,'2':[{'taxonomy':'B'}]}
    result=c.summary(rows,matched,1500)
    assert result['coverage']=={'pool':3,'jpl':1,'mithneos':2,'both':1,'either':2,
                                'newly_observed_in_mithneos':1,'neither':1}
    d=result['integrated_classifications']
    assert d['counts']=={'C_B':1,'S_Q':1,'no_observation':1}
    assert d['known_size_mass_kg']['by_observed_family']['C_B']==pytest.approx(8*d['known_size_mass_kg']['by_observed_family']['S_Q'])
    assert d['known_size_mass_kg']['missing_measure_objects']==1
    assert d['count_inference']['conditional_observed_mixture_completion']['family_measure']=={'C_B':1.5,'S_Q':1.5}


def test_bounds_and_completion_conserve_measure_without_becoming_confidence_limits():
    states=[{'family':'C_B'},{'family':'X'},{'family':'C_or_X'},
            {'family':'no_observation'},{'family':'conflicting_families'}]
    d=c.measure_distribution(states,[2,1,1,4,2])
    assert d['total_known_measure']==10 and d['unresolved_measure']==6
    assert d['identification_bounds']=={'C_B':[.2,.8],'C_or_X':[.1,.7],'X':[.1,.7]}
    scenario=d['conditional_observed_mixture_completion']
    assert scenario['family_measure']=={'C_B':5,'C_or_X':2.5,'X':2.5}
    assert sum(scenario['family_measure'].values())==10
    assert d['unseen_family_upper_fraction']==.6
    assert c.measure_distribution([{'family':'unclassified'}],[1])['conditional_observed_mixture_completion'] is None
    with pytest.raises(ValueError):c.measure_distribution(states,[1])
    with pytest.raises(ValueError):c.measure_distribution([states[0]],[-1])


def test_geometry_limits_dark_selection_and_unknown_mass_are_independent():
    limits={'a_min_au':.8,'a_max_au':1.3,'e_max':.3,'i_max_deg':15,'albedo_max':.1,'excluded_numbers':[99]}
    r=row();assert c.selected(r,limits)
    for field,value in [('a','1.31'),('a','.79'),('e','.31'),('i','15.1'),('diameter',None),('number',99)]:
        assert not c.selected({**r,field:value},limits)
    assert not c.selected({**r,'albedo':None},limits)
    assert c.selected({**r,'albedo':None},limits,dark=False)
    assert c.selected({**r,'albedo':'.9'},limits,dark=False)
    assert c.assumed_mass(r,1500)==pytest.approx(math.pi/6*1e9*1500)
    assert c.assumed_mass({**r,'diameter':None},1500) is None
    with pytest.raises(ValueError):c.assumed_mass(r,0)


def test_hash_failure_prevents_substitution(tmp_path):
    f=tmp_path/'external';f.write_bytes(b'changed')
    with pytest.raises(ValueError,match='hash or length mismatch'):
        c.checked_bytes(f,{'sha256':c.digest(b'original'),'bytes':8})


def test_frozen_sources_product_and_existing_inventory_ids_reproduce():
    config=json.loads((c.HERE/'composition_inputs.json').read_text())
    archive=c.ROOT/config['archive_directory']
    if not (archive/config['survey']['filename']).exists() or not Path(config['jpl']['path']).exists():
        pytest.skip('Exact external catalogue inputs absent; restore the pinned originals')
    result=c.product()
    saved=json.loads(gzip.decompress((c.HERE/'results/composition.json.gz').read_bytes()))
    assert result==saved
    assert result['survey_account']['observations']==491
    assert result['survey_account']['distinct_number_designation_identities']==419
    assert result['survey_account']['distinct_filename_prefixes']==420
    assert result['survey_account']['matched_jpl_objects']==417
    assert result['survey_account']['matched_observations']==489
    assert {o['designation'] for o in result['survey_account']['unmatched_observations']}=={'2011 JT9','2014 LW21'}
    wide=result['pools']['wider_dark_screen']
    assert wide['coverage']['either']==17 and wide['coverage']['neither']==114
    assert wide['integrated_classifications']['counts']=={'C_B':6,'C_or_X':6,'S_Q':3,'X':2,'no_observation':114}
    reference=json.loads((c.ROOT/'resources/enclosed_mining/inputs/geometry_screen.json').read_text())['groups']['wider_geometry_excluding_named_science_targets']
    assert {r['spkid'] for r in result['wider_objects']}=={str(x) for x in reference['selected_object_ids']}
    mass=wide['integrated_classifications']['known_size_mass_kg']
    assert mass['total_known_measure']==pytest.approx(reference['mass_at_assumed_density_kg'])
    assert result['pools']['narrow_dark_screen']['coverage']['pool']==17
    assert result['pools']['wider_geometry_all_albedos']['coverage']['pool']==382
    for pool in result['pools'].values():
        assert pool['population_estimate'] is None and pool['statistical_confidence_interval'] is None
        for key in ('jpl_only_classifications','mithneos_only_classifications','integrated_classifications'):
            d=pool[key]
            assert sum(d['counts'].values())==d['object_count']
            assert sum(d['fractions_of_pool'].values())==pytest.approx(1)
            for measure in ('count_inference','known_size_mass_kg'):
                m=d[measure];assert sum(m['by_observed_family'].values())==pytest.approx(m['total_known_measure'])
                completion=m['conditional_observed_mixture_completion']
                if completion:assert sum(completion['family_measure'].values())==pytest.approx(m['total_known_measure'])
    assert result['grade_inference'] is None and result['qualified_reserve'] is None
