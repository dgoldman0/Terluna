import gzip
import hashlib
import json

import pytest

from resources.enclosed_mining import run as m


def test_two_tiers_do_not_double_count_the_nested_catalogue():
    inventory = {"inventory": {"groups": {
        "dark_geometry_excluding_named_science_targets": {"mass_at_assumed_density_kg": 100.},
        "dark_with_diameter": {"mass_at_assumed_density_kg": 1000.}}}}
    tiers = m.catalogue_tiers(inventory, {"narrow_extractable_fraction": .1,
                                        "broader_remainder_extractable_fraction": .1})
    assert [t["feed_allowance_kg"] for t in tiers] == [10., 90.]
    assert sum(t["feed_allowance_kg"] for t in tiers) == 100.


def test_finite_allowances_roll_over_then_leave_unallocated_demand():
    tiers = [{"id": "first", "feed_allowance_kg": 10.},
             {"id": "second", "feed_allowance_kg": 5.}]
    result = m.allocate_feed([{"year": 0, "source_feed_kg": 8.},
                              {"year": 1, "source_feed_kg": 10.}], tiers)
    assert result["annual"][1]["tier_draw_kg"] == {"first": 2., "second": 5.}
    assert result["unallocated_feed_kg"] == 3.
    assert result["requested_feed_kg"] == 18.
    assert result["supplied_under_scenario_kg"] == 15.
    assert result["first_unallocated_year"] == 1
    assert not result["within_hypothetical_mass_allowance"]
    assert result["qualified_reserve_sufficiency"] is None
    assert tiers[0]["feed_allowance_kg"] == 10.


def test_annual_access_limit_and_release_date_preserve_unmet_feed():
    tiers = [{"id": "later", "feed_allowance_kg": 100., "available_year": 1,
              "annual_feed_limit_kg": 5.}]
    result = m.allocate_feed([{"year": y, "source_feed_kg": 8.} for y in range(3)], tiers)
    assert [r["supplied_under_scenario_kg"] for r in result["annual"]] == [0., 5., 5.]
    assert result["unallocated_feed_kg"] == 14.
    assert result["annual"][-1]["remaining_by_tier_kg"]["later"] == 90.


def test_footprint_distinguishes_disturbed_and_prospective_area():
    f = m.footprint(1e12, {"bulk_density_kg_m3": 1500., "worked_depth_m": 5.,
                          "accessible_area_fraction": .25})
    assert f["excavated_area_km2"] == pytest.approx(133.33333333333334)
    assert f["prospective_area_km2"] == pytest.approx(533.3333333333334)


def test_zero_demand_neither_consumes_allowance_nor_proves_reserves():
    result = m.allocate_feed([{"year": 0, "source_feed_kg": 0.}], [])
    assert result["within_hypothetical_mass_allowance"]
    assert result["supplied_under_scenario_kg"] == 0.
    assert result["qualified_reserve_sufficiency"] is None


@pytest.mark.parametrize("value", [-1., float("inf"), float("nan"), True])
def test_invalid_feed_is_rejected(value):
    with pytest.raises(ValueError):
        m.allocate_feed([{"year": 0, "source_feed_kg": value}], [])


def test_repeated_year_or_tier_identity_is_rejected():
    with pytest.raises(ValueError):
        m.allocate_feed([{"year": 0, "source_feed_kg": 0.}]*2, [])
    with pytest.raises(ValueError):
        m.allocate_feed([], [{"id": "same", "feed_allowance_kg": 1.}]*2)


def test_changed_source_bytes_and_schema_are_rejected(tmp_path):
    p = tmp_path/"product.json"
    p.write_text(json.dumps({"schema": "example/1"}))
    spec = {"path": p.name, "sha256": m.digest(p), "schema": "example/1"}
    assert m.load_checked(spec, tmp_path)["schema"] == "example/1"
    p.write_text(json.dumps({"schema": "example/2"}))
    with pytest.raises(ValueError, match="Unreviewed source"):
        m.load_checked(spec, tmp_path)
    spec["sha256"] = m.digest(p)
    with pytest.raises(ValueError, match="Unexpected source schema"):
        m.load_checked(spec, tmp_path)


def test_pinned_catalogue_product_and_fraction_scenarios():
    q = json.loads((m.HERE/"inputs.json").read_text())
    inventory = m.load_checked(q["inventory_product"])
    narrow = m.catalogue_tiers(inventory, q["access_scenarios"]["narrow_only"])
    assert narrow[0]["feed_allowance_kg"] == pytest.approx(3.1188112554468115e11)
    assert narrow[1]["feed_allowance_kg"] == 0.
    broad = m.catalogue_tiers(inventory, q["access_scenarios"]["broader_ten_percent"])
    assert sum(t["feed_allowance_kg"] for t in broad) == pytest.approx(8.511203050095491e14)


def test_seed_is_included_once_before_actual_operation_with_nested_fields():
    case = {"seed": {"bootstrap_years": 2, "lunar_raw_kg": 100., "source_feed_kg": 60.},
            "annual": [{"year": 2, "physical": {"lunar_raw_kg_year": 7., "source_feed_kg_year": 11.}},
                       {"year": 3, "physical": {"lunar_raw_kg_year": 13., "source_feed_kg_year": 17.}}]}
    fields = {"year": "year", "lunar_virgin_feed_kg": ["physical", "lunar_raw_kg_year"],
              "source_feed_kg": ["physical", "source_feed_kg_year"]}
    rows = m.normalize_case(case, fields)
    assert [r["year"] for r in rows] == [0, 1, 2, 3]
    assert [r["period_kind"] for r in rows] == ["bootstrap_average"]*2+["operation"]*2
    assert sum(r["lunar_virgin_feed_kg"] for r in rows) == 120.
    assert sum(r["source_feed_kg"] for r in rows) == 88.
    allocation = m.allocate_feed(rows, [{"id": "finite", "feed_allowance_kg": 65.}])
    assert allocation["first_unallocated_year"] == 2
    assert allocation["unallocated_feed_kg"] == 23.


@pytest.mark.parametrize("duration", [0, -1, 2.5, True])
def test_seed_duration_must_be_positive_integer(duration):
    case = {"seed": {"bootstrap_years": duration, "lunar_raw_kg": 100., "source_feed_kg": 60.},
            "annual": [{"year": 8}]}
    with pytest.raises(ValueError, match="Bootstrap duration"):
        m.normalize_case(case, {"year": "year"})


def test_case_and_grade_scenario_identities_prevent_overwriting_paths():
    cfg = json.loads((m.HERE/"inputs.json").read_text())
    cfg["annual_fields"] = {"year": "year", "lunar_virgin_feed_kg": "lunar", "source_feed_kg": "source"}
    inventory = m.inventory_inputs(cfg)
    base = {"name": "same-market", "physical_scenario": "base",
            "seed": {"bootstrap_years": 1, "lunar_raw_kg": 2., "source_feed_kg": 3.,
                     "evidence": "synthetic seed"},
            "annual": [{"year": 1, "lunar": 4., "source": 5.}],
            "convergence": {"converged": False}, "first_capacity_shortfall": 1}
    poor = dict(base, physical_scenario="poor")
    results = m.assess_physical_cases({"cases": [base, poor]}, inventory, cfg)
    assert set(results) == {"base:same-market", "poor:same-market"}
    for result in results.values():
        assert result["cumulative_lunar_virgin_feed_kg"] == 6.
        assert result["cumulative_source_feed_kg"] == 8.
        assert not result["physical_convergence"]["converged"]
        assert result["first_physical_capacity_shortfall"] == 1
        assert result["source_access"]["narrow_only"]["within_hypothetical_mass_allowance"]
        assert result["source_access"]["narrow_only"]["qualified_reserve_sufficiency"] is None
    with pytest.raises(ValueError, match="Repeated physical"):
        m.assess_physical_cases({"cases": [base, base]}, inventory, cfg)


def test_compressed_snapshot_checks_both_archive_and_exact_source_bytes(tmp_path):
    raw = b'{"schema":"compressed-example/1","value":1}\n'
    p = tmp_path/"product.json.gz"
    p.write_bytes(gzip.compress(raw,mtime=0))
    spec = {"path": p.name, "sha256": m.digest(p), "schema": "compressed-example/1",
            "compression": "gzip", "content_sha256": hashlib.sha256(raw).hexdigest()}
    assert m.load_checked(spec,tmp_path)["value"] == 1
    spec["content_sha256"] = "0"*64
    with pytest.raises(ValueError,match="decompressed source"):
        m.load_checked(spec,tmp_path)
    spec["content_sha256"] = hashlib.sha256(raw).hexdigest()
    p.write_bytes(gzip.compress(raw.replace(b'1}',b'2}'),mtime=0))
    with pytest.raises(ValueError,match="Unreviewed source product"):
        m.load_checked(spec,tmp_path)


def test_excluding_largest_and_wider_geometry_have_distinct_finite_ceilings():
    cfg = json.loads((m.HERE/'inputs.json').read_text())
    inv = m.inventory_inputs(cfg)
    scenarios = cfg['access_scenarios']
    without = m.catalogue_tiers(inv, scenarios['broader_ten_percent_without_largest'])
    assert sum(t['feed_allowance_kg'] for t in without) == pytest.approx(3.1241570473523944e14)
    assert without[1]['excluded_objects'][0]['spkid'] == 20003552
    broad = m.catalogue_tiers(inv, scenarios['broader_ten_percent'])
    assert without[0] == broad[0]
    wide = m.catalogue_tiers(inv, scenarios['wider_geometry_all_mass_ceiling'])
    assert sum(t['feed_allowance_kg'] for t in wide) == pytest.approx(2.143151330135278e14)
    assert wide[1]['catalogue_group'] == 'wider_geometry_excluding_named_science_targets'
    exhausted = m.allocate_feed([{'year':200,'source_feed_kg':2.3e14}],wide)
    assert not exhausted['within_hypothetical_mass_allowance']
    assert exhausted['qualified_reserve_sufficiency'] is None


def test_wider_geometry_applies_every_boundary_and_excludes_unknown_diameters():
    from resources.enclosed_mining import catalogue as c
    fields=['spkid','full_name','diameter','albedo','a','e','i']
    rows=[[1,'boundary',1.,.1,1.8,.5,30.], [2,'too distant',1.,.1,1.8001,.5,30.],
          [3,'too eccentric',1.,.1,1.8,.5001,30.], [4,'too inclined',1.,.1,1.8,.5,30.0001],
          [5,'missing',None,.1,1.8,.5,30.], [6,'Bennu scientific comparison',1.,.1,1.8,.5,30.]]
    cfg=json.loads((m.HERE/'inputs.json').read_text())
    out=c.summarize({'count':len(rows),'fields':fields,'data':rows},1500.,cfg['wider_geometry_screen']['limits'])
    assert out['selected_object_ids'] == [1]
    assert out['count'] == 1
    assert out['mass_at_assumed_density_kg'] == pytest.approx(1500.*__import__('math').pi/6*1e9)


def test_additional_catalogue_product_reproduces_exact_archived_source():
    from pathlib import Path
    from resources.enclosed_mining import catalogue as c
    cfg=json.loads((m.HERE/'inputs.json').read_text())
    archive=Path(cfg['catalogue_source']['archive'])
    if not archive.exists():
        pytest.skip('Exact external JPL catalogue archive is unavailable; frozen assessment inputs remain independently pinned')
    saved=m.load_checked(cfg['geometry_product'])
    assert c.product(cfg) == saved


def test_checkpoint_condition_preserves_prior_shortage_and_unknown_qualification():
    inventory={'inventory':{'groups':{
        'dark_geometry_excluding_named_science_targets':{'mass_at_assumed_density_kg':10.},
        'dark_with_diameter':{'mass_at_assumed_density_kg':10.}}}}
    cfg={'lunar_footprint':{'bulk_density_kg_m3':1500.,'worked_depth_m':5.,'accessible_area_fraction':.25},
         'access_scenarios':{'finite':{'narrow_extractable_fraction':1.,'broader_remainder_extractable_fraction':0.}},
         'checkpoints':[1,2]}
    out=m.assess_case([{'year':y,'source_feed_kg':v,'lunar_virgin_feed_kg':0.} for y,v in [(1,12.),(2,0.)]],inventory,cfg)
    for row in out['source_access_checkpoints']:
        q=row['comparisons']['finite']
        assert q['cumulative_unallocated_feed_kg']==2.
        assert not q['within_hypothetical_mass_allowance_to_checkpoint']
        assert q['qualified_resource_and_delivery_schedule'] is None


def test_saved_assessment_reproduces_and_includes_every_seed_and_annual_feed():
    import math
    saved=json.loads(gzip.decompress((m.HERE/'results/assessment.json.gz').read_bytes()))
    assert m.product()==saved
    assert saved['producer']['sha256']==m.digest(m.HERE/'run.py')
    assert len(saved['cases'])==87
    config=json.loads((m.HERE/'inputs.json').read_text())
    for spec in config['physical_products']:
        physical=m.load_checked(spec)
        for case in physical['cases']:
            result=saved['cases'][case['physical_scenario']+':'+case['name']]
            assert result['cumulative_lunar_virgin_feed_kg']==pytest.approx(case['seed']['lunar_raw_kg']+math.fsum(r['physical']['lunar_raw_kg_year'] for r in case['annual']),rel=2e-14)
            assert result['cumulative_source_feed_kg']==pytest.approx(case['seed']['source_feed_kg']+math.fsum(r['physical']['source_feed_kg_year'] for r in case['annual']),rel=2e-14)
            assert len(result['annual'])==case['seed']['bootstrap_years']+len(case['annual'])
            for comparison in result['source_access'].values():
                assert comparison['qualified_reserve_sufficiency'] is None
                for row in comparison['annual']:
                    assert row['supplied_under_scenario_kg']+row['unallocated_feed_kg']==pytest.approx(row['requested_feed_kg'],rel=2e-14,abs=1e-6)
                    assert all(v>=0 for v in row['remaining_by_tier_kg'].values())
            for driver in result['supply_driver_checkpoints']:
                source=physical['scenarios']['source']
                expected=max(driver['source_water_kg_year']/(source['water_grade']*source['water_recovery']),
                             driver['source_nitrogen_kg_year']/(source['nitrogen_grade']*source['nitrogen_recovery']),
                             driver['source_virgin_metal_kg_year']/(source['metal_grade']*source['metal_recovery']))
                assert driver['source_feed_kg_year']==pytest.approx(expected,rel=2e-14)
                assert driver['resident_guest_and_initial_water_kg_year']>=0.
    base=[c for c in saved['cases'].values() if c['physical_scenario']=='base']
    assert len(base)==81
    for c in base:
        cp={r['year']:r['comparisons'] for r in c['source_access_checkpoints']}
        assert not cp[100]['narrow_all_mass_ceiling']['within_hypothetical_mass_allowance_to_checkpoint']
        assert cp[200]['broader_ten_percent_without_largest']['within_hypothetical_mass_allowance_to_checkpoint']
        assert not cp[220]['wider_geometry_all_mass_ceiling']['within_hypothetical_mass_allowance_to_checkpoint']
    for poor in (c for c in saved['cases'].values() if c['physical_scenario']=='poor_source'):
        reference=saved['cases']['base:'+poor['physical_case_name'].removesuffix(':poor_source')]
        assert poor['cumulative_source_feed_kg']>reference['cumulative_source_feed_kg']
        cp=next(r for r in poor['source_access_checkpoints'] if r['year']==200)['comparisons']
        assert not cp['broader_ten_percent_without_largest']['within_hypothetical_mass_allowance_to_checkpoint']
