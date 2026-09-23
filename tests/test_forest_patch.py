import unittest
from dataclasses import replace
from pathlib import Path
import math
import numpy as np
from atmosphere.lower_air import LowerAir
from biosphere.megaforest_mechanics import TreeConfig,RootConfig
from climate.forest_patch_flow import FlowConfig,PatchGrid
from biosphere.forest_patch import PatchTree,PatchDrag,make_patch,tree_mass
from biosphere.forest_patch_mechanics import InteractionConfig,PatchStructure,partition_soil,point_shape
from biosphere.forest_damage import damage_sequence
from research.run_forest_patch import safe_output,ROOT


class FlowTests(unittest.TestCase):
    def setUp(self):self.grid=PatchGrid(FlowConfig(nx=24,ny=16,nz=12))

    def test_empty_wind_tunnel_is_uniform(self):
        r=self.grid.solve(np.zeros(self.grid.shape))
        np.testing.assert_allclose(r.centered[...,0],1.,atol=1e-13)
        np.testing.assert_allclose(r.centered[...,1:],0.,atol=1e-13)
        self.assertEqual(r.diagnostics['steps'],1)

    def test_projection_removes_gradient(self):
        g=self.grid;x,y,z=np.meshgrid(*g.centers,indexing='ij')
        p=np.cos(np.pi*(x+6)/12)*np.cos(np.pi*(z)/3)
        v=g.faces()
        for a in range(3):
            s=[slice(None)]*3;s[a]=slice(1,-1)
            v[a][tuple(s)]+=np.diff(p,axis=a)/g.spacing[a]
        v,_=g.project(v)
        for actual,expected in zip(v,g.faces()):np.testing.assert_allclose(actual,expected,atol=1e-12)

    def test_projection_random_field(self):
        g=self.grid;rng=np.random.default_rng(42)
        v=[a+rng.normal(0,.1,a.shape) for a in g.faces()];g._boundary(v)
        v,_=g.project(v)
        self.assertLess(np.max(abs(g.divergence(v))),1e-12)

    def test_incompatible_flux_rejected(self):
        v=self.grid.faces();v[0][0]=2
        with self.assertRaises(ValueError):self.grid.project(v)

    def test_negative_or_nan_drag_rejected(self):
        for bad in (-1.,np.nan):
            with self.assertRaises(ValueError):self.grid.solve(np.full(self.grid.shape,bad))

    def test_drag_slab_pressure_balance(self):
        g=self.grid;b=2.
        r=g.solve(np.full(g.shape,b))
        np.testing.assert_allclose(r.centered[...,0],1.,atol=1e-12)
        np.testing.assert_allclose(np.diff(r.pressure,axis=0)/g.spacing[0],-b/2,atol=1e-11)

    def test_porosity_deflects_flow_and_closes_momentum(self):
        g=self.grid;x,y,z=np.meshgrid(*g.centers,indexing='ij');b=np.zeros(g.shape)
        b[(abs(x)<1)&(abs(y)<1)&(z>.3)&(z<1.2)]=3
        r=g.solve(b)
        self.assertLess(r.centered[b>0,0].mean(),.95)
        self.assertGreater(np.max(abs(r.centered[...,2])),.01)
        self.assertLess(r.diagnostics['momentum_relative_residual'],1e-4)
        self.assertLess(r.diagnostics['max_divergence'],1e-11)

    def test_nonconvergence_raises(self):
        g=PatchGrid(replace(self.grid.config,max_steps=4,tolerance=1e-15))
        b=np.zeros(g.shape);b[8:12,6:9,1:5]=4
        with self.assertRaises(RuntimeError):g.solve(b)


class GeometryTests(unittest.TestCase):
    def setUp(self):self.air=LowerAir();self.grid=PatchGrid(FlowConfig(nx=24,ny=16,nz=12))

    def test_all_shapes_match_initial_mass(self):
        for shape in ('flat','dome','ramp','irregular','gap'):
            trees,m=make_patch(shape)
            self.assertLess(abs(m['initial_mass_kg']/m['reference_flat_mass_kg']-1),1e-11)
        self.assertGreater(m['max_height_m'],300.)

    def test_rotation_preserves_pair_distances(self):
        a,_=make_patch(direction_deg=0);b,_=make_patch(direction_deg=67)
        da=math.hypot(a[0].x_m-a[-1].x_m,a[0].y_m-a[-1].y_m)
        db=math.hypot(b[0].x_m-b[-1].x_m,b[0].y_m-b[-1].y_m)
        self.assertAlmostEqual(da,db)

    def test_projected_area_analytic_check(self):
        p=PatchTree('one',0,0,TreeConfig(),RootConfig());d=PatchDrag(self.grid,[p]);t=p.tree
        crown=np.pi*t.crown_radius_m*(1-t.crown_base_fraction)*t.height_m/2*t.crown_frontal_fraction*t.crown_drag_coefficient
        stem=t.base_diameter_m*(1+t.tip_diameter_fraction)*t.height_m/2*t.stem_drag_coefficient
        self.assertLess(abs(d.entries[0]['target_cd_area_m2']/(crown+stem)-1),.001)
        self.assertLess(d.area_relative_error,1e-12)

    def test_tree_drag_exactly_partitions_grid(self):
        trees,_=make_patch(nrows=3,ncols=3);d=PatchDrag(self.grid,trees);r=self.grid.solve(d.coefficient)
        loads,l=d.loads(r,13.,1.45)
        self.assertEqual(len(loads),9);self.assertLess(l['force_partition_relative_error'],1e-12)

    def test_shedding_reduces_drag_and_standing_mass(self):
        p=PatchTree('one',0,0,TreeConfig(),RootConfig());q=replace(p,remaining_crown=.35)
        self.assertLess(tree_mass(q.effective_tree()),tree_mass(p.tree))
        self.assertLess(PatchDrag(self.grid,[q]).coefficient.sum(),PatchDrag(self.grid,[p]).coefficient.sum())

    def test_failure_removes_aerodynamic_obstacle(self):
        p=PatchTree('one',0,0,TreeConfig(),RootConfig(),alive=False)
        self.assertEqual(PatchDrag(self.grid,[p]).coefficient.sum(),0)

    def test_clipped_tree_rejected(self):
        p=PatchTree('one',10000,0,TreeConfig(),RootConfig())
        with self.assertRaises(ValueError):PatchDrag(self.grid,[p])

    def test_duplicate_tree_ids_rejected(self):
        p=PatchTree('one',0,0,TreeConfig(),RootConfig())
        with self.assertRaises(ValueError):PatchDrag(self.grid,[p,p])

    def test_unimplemented_streamlining_rejected(self):
        p=PatchTree('one',0,0,replace(TreeConfig(),vogel_exponent=-.8),RootConfig())
        with self.assertRaises(ValueError):PatchDrag(self.grid,[p])

    def test_stale_flow_for_new_geometry_is_rejected(self):
        p=PatchTree('one',0,0,TreeConfig(),RootConfig());d=PatchDrag(self.grid,[p])
        r=self.grid.solve(d.coefficient)
        altered=PatchDrag(self.grid,[replace(p,remaining_crown=.35)])
        with self.assertRaises(ValueError):altered.loads(r,15.,1.45)

    def test_nonfinite_reference_density_rejected(self):
        p=PatchTree('one',0,0,TreeConfig(),RootConfig());d=PatchDrag(self.grid,[p])
        r=self.grid.solve(d.coefficient)
        with self.assertRaises(ValueError):d.loads(r,15.,float('nan'))

    def test_force_zero_at_zero_wind(self):
        p=PatchTree('one',0,0,TreeConfig(),RootConfig());d=PatchDrag(self.grid,[p])
        r=self.grid.solve(d.coefficient);l,_=d.loads(r,0.,1.45)
        np.testing.assert_equal(l['one']['total_force_n'],0.)


class StructureTests(unittest.TestCase):
    def setUp(self):
        self.air=LowerAir();t=TreeConfig();r=RootConfig()
        self.trees=[PatchTree('a',-52.5,0,t,r),PatchTree('b',52.5,0,t,r)]

    def loads(self,fa=1e5,fb=0):
        return {p.identifier:dict(z_m=np.array([250.]),force_n=np.array([[f,0.,0.]]),
            crown_force_n=np.array([[f,0.,0.]]),total_force_n=np.array([f,0.,0.]),
            drag_weighted_speed_ratio=1.) for p,f in zip(self.trees,[fa,fb])}

    def test_hermite_rigid_translation_and_rotation(self):
        z=np.linspace(0,100,9);d=np.zeros(18);d[::2]=2+z*.1;d[1::2]=.1
        for h in (0,3,27,93,100):self.assertAlmostEqual(point_shape(z,h)@d,2+.1*h)

    def test_analytical_flexible_base_cantilever(self):
        t=replace(TreeConfig(),height_m=100,base_diameter_m=5,tip_diameter_fraction=1.,
            wood_density_kg_m3=1e-9,crown_wood_to_stem_mass=0,leaf_area_index=0,
            epiphyte_mass_kg_m2_footprint=0,retained_water_mm=0,elements=16)
        p=PatchTree('a',0,0,t,RootConfig());model=PatchStructure(self.air,[p])
        f=1e4;load=dict(z_m=np.array([100.]),force_n=np.array([[f,0.,0.]]),
            crown_force_n=np.array([[0.,0.,0.]]),total_force_n=np.array([f,0.,0.]),drag_weighted_speed_ratio=1.)
        r=model.solve({'a':load})['trees'][0]
        expected=f*100**3/(3*t.young_modulus_pa*np.pi*5**4/64)+f*100**2/model.info['a']['rootk']
        self.assertLess(abs(r['tip_displacement_m']/expected-1),1e-7)
        self.assertLess(abs(r['base_moment_nm']/(f*100)-1),1e-7)

    def test_contact_transfers_load_to_neighbor(self):
        a=PatchStructure(self.air,self.trees).solve(self.loads())
        b=PatchStructure(self.air,self.trees,InteractionConfig(crown_contacts=True)).solve(self.loads())
        self.assertLess(b['trees'][0]['base_moment_nm'],a['trees'][0]['base_moment_nm'])
        self.assertGreater(b['trees'][1]['base_moment_nm'],0)
        np.testing.assert_allclose(b['net_internal_contact_force_n'],0,atol=1e-12)
        self.assertLess(b['linear_equilibrium_relative_residual'],1e-8)

    def test_contacts_do_not_pull(self):
        r=PatchStructure(self.air,self.trees,InteractionConfig(crown_contacts=True)).solve(self.loads(0,1e5))
        self.assertFalse(any(l['active'] for l in r['links']))

    def test_broken_contacts_are_removed(self):
        c=InteractionConfig(crown_contacts=True)
        m=PatchStructure(self.air,self.trees,c)
        r=PatchStructure(self.air,self.trees,c,[l['id'] for l in m.links]).solve(self.loads())
        self.assertEqual(r['links'],[])

    def test_grafts_conserve_internal_moments(self):
        m=PatchStructure(self.air,self.trees,InteractionConfig(root_grafts=True));r=m.solve(self.loads())
        np.testing.assert_allclose(r['net_internal_graft_moment_nm'],0,atol=1e-12)
        self.assertGreater(r['trees'][1]['root_soil_moment_nm'],0)

    def test_soil_overlap_is_counted_once(self):
        ts=[replace(p,roots=replace(p.roots,plate_radius_m=100.)) for p in self.trees]
        a,d=partition_soil(ts,5.)
        self.assertLess(abs(d['partition_error_m2']),1e-8)
        self.assertLess(a['a']['capacity_fraction'],1.)
        self.assertLess(a['b']['capacity_fraction'],1.)

    def test_failed_roots_do_not_free_instant_anchorage(self):
        ts=[replace(p,roots=replace(p.roots,plate_radius_m=100.)) for p in self.trees]
        a,_=partition_soil(ts,5.);ts[0]=replace(ts[0],alive=False);b,_=partition_soil(ts,5.)
        self.assertEqual(a['b']['capacity_fraction'],b['b']['capacity_fraction'])

    def test_structure_retains_failed_root_reservations(self):
        ts=[replace(p,roots=replace(p.roots,plate_radius_m=100.)) for p in self.trees]
        old=PatchStructure(self.air,ts)
        ts[0]=replace(ts[0],alive=False)
        new=PatchStructure(self.air,ts)
        self.assertEqual(old.soil['b']['capacity_fraction'],new.soil['b']['capacity_fraction'])
        self.assertLess(new.soil['b']['capacity_fraction'],1.)

    def test_invalid_soil_cell_size_rejected(self):
        with self.assertRaises(ValueError):partition_soil(self.trees,0.)

    def test_contact_overload_is_reported(self):
        c=InteractionConfig(crown_contacts=True,contact_break_force_n=1.)
        r=PatchStructure(self.air,self.trees,c).solve(self.loads())
        self.assertGreater(max(l['utilization'] for l in r['links']),1.)

    def test_large_deflection_is_flagged(self):
        r=PatchStructure(self.air,self.trees).solve(self.loads(1e8,0))
        self.assertTrue(r['domain_limited'])


class DamageAndScopeTests(unittest.TestCase):
    def test_domain_limit_stops_damage_claim(self):
        air=LowerAir();g=PatchGrid(FlowConfig(nx=24,ny=16,nz=12));ts,_=make_patch(nrows=2,ncols=2)
        r=damage_sequence(air,g,ts,[100])
        self.assertEqual(r['status'],'MODEL_DOMAIN_LIMIT');self.assertEqual(r['events'],[])

    def test_zero_wind_preserves_stand_and_ledger(self):
        air=LowerAir();g=PatchGrid(FlowConfig(nx=24,ny=16,nz=12));ts,_=make_patch(nrows=2,ncols=2)
        r=damage_sequence(air,g,ts,[0])
        self.assertEqual(r['standing_tree_count'],4);self.assertEqual(r['mass_ledger_residual_kg'],0.)

    def test_decreasing_wind_rejected(self):
        with self.assertRaises(ValueError):damage_sequence(LowerAir(),None,[],[2,1])

    def test_immersion_output_guard(self):
        for p in (ROOT/'immersion',ROOT/'immersion'/'x',ROOT/'unknown'/'..'/'immersion'/'x'):
            with self.assertRaises(ValueError):safe_output(p)

    def test_hypothetical_connections_default_off(self):
        c=InteractionConfig();self.assertFalse(c.crown_contacts);self.assertFalse(c.root_grafts)

if __name__=='__main__':unittest.main()
