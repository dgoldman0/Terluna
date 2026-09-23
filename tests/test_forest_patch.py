"""Analytic, limiting, conservation, interaction and scope tests for forest patch."""
import unittest
from dataclasses import replace
from pathlib import Path
import numpy as np
from atmosphere.lower_air import LowerAir
from climate.forest_flow import Grid,SpatialFlow
from biosphere.forest_patch import Individual,make_patch,TreeDrag,biomass_proxy
from biosphere.forest_mechanics import Support,StandMechanics,share_soil,shape_row
from research.run_forest_patch import output_guard,ROOT,envelope,damage_sequence


class PatchTests(unittest.TestCase):
    def setUp(self):
        self.air=LowerAir();self.g=Grid(nx=24,ny=20,nz=16)
        self.f,_=make_patch(nx=3,ny=3,elements=12)

    def test_matched_budgets(self):
        for shape in ('flat','dome','tapered','ramp','irregular','gap'):
            f,m=make_patch(shape,nx=5,ny=5)
            self.assertAlmostEqual(m['aboveground_mass_proxy_kg']/m['target_mass_proxy_kg'],1,places=10)

    def test_equal_height_is_not_equal_mass(self):
        _,a=make_patch('flat',matched_mass=False);_,b=make_patch('dome',matched_mass=False)
        self.assertLess(b['aboveground_mass_proxy_kg'],a['aboveground_mass_proxy_kg'])
        self.assertEqual(b['maximum_height_m'],300.)

    def test_geometry_seed_reproducible(self):
        a,_=make_patch('irregular');b,_=make_patch('irregular');self.assertEqual(a,b)

    def test_invalid_shape(self):
        with self.assertRaises(ValueError):make_patch('invented')

    def test_grid_validation(self):
        for kwargs in ({'nx':2},{'top':-1},{'mixing_length_m':float('nan')},{'reservoir_fraction':.6}):
            with self.assertRaises(ValueError):SpatialFlow(replace(self.g,**kwargs))

    def test_stem_integral(self):
        d=TreeDrag(self.f,self.g,self.air)
        actual=d.area[d.kind==0].sum()
        expected=sum(f.tree.stem_drag_coefficient*f.tree.base_diameter_m*f.tree.height_m*(1+f.tree.tip_diameter_fraction)/2 for f in self.f)
        self.assertAlmostEqual(actual/expected,1,places=12)

    def test_drag_rows_normalized(self):
        d=TreeDrag(self.f,self.g,self.air)
        np.testing.assert_allclose(np.asarray(d.W.sum(axis=1)).ravel(),1,atol=2e-15)

    def test_exchange_conservation(self):
        d=TreeDrag(self.f,self.g,self.air);rng=np.random.default_rng(4)
        v=rng.normal(size=(3,*self.g.shape));_,ledger=d.forces(v,17.)
        self.assertLess(ledger['relative_residual'],1e-12)

    def test_rigid_force_scaling(self):
        d=TreeDrag(self.f,self.g,self.air);v=np.ones((3,*self.g.shape))
        a,_=d.forces(v,10);b,_=d.forces(v,20);np.testing.assert_allclose(b,4*a)

    def test_reconfiguration_and_envelope_guard(self):
        f=[replace(x,tree=replace(x.tree,vogel_exponent=-.8)) for x in self.f]
        d=TreeDrag(f,self.g,self.air);v=np.ones((3,*self.g.shape))
        a,_=d.forces(v,10);b,_=d.forces(v,20)
        self.assertLess(np.linalg.norm(b[d.kind==1]),4*np.linalg.norm(a[d.kind==1]))
        with self.assertRaises(ValueError):envelope(None,d,a)

    def test_drag_domain_guard(self):
        f=[replace(self.f[0],x=1100)]
        with self.assertRaises(ValueError):TreeDrag(f,self.g,self.air)

    def test_duplicate_ids(self):
        with self.assertRaises(ValueError):TreeDrag([self.f[0],self.f[0]],self.g,self.air)

    def test_exact_vacuum(self):
        d=TreeDrag([],self.g,self.air);v,F,D=SpatialFlow(self.g,30).solve(d)
        np.testing.assert_allclose(v[0],np.cos(np.pi/6));np.testing.assert_allclose(v[1],.5)
        self.assertEqual(F.shape,(0,3));self.assertEqual(D['outer_iterations'],0)

    def test_projection_is_idempotent(self):
        s=SpatialFlow(self.g);v=np.random.default_rng(9).normal(size=s.shape)
        h=s.project(s.fft(v));np.testing.assert_allclose(s.project(h),h,atol=1e-10)
        self.assertLess(np.max(abs(s.inv(1j*np.sum(s.k*h,axis=0)))),1e-13)

    def test_uniform_drag_analytic(self):
        class Uniform:
            def coefficient(obj,v,speed):return np.full(v.shape[1:],.002)
            def forces(obj,v,speed):return np.empty((0,3)),{'relative_residual':0.}
        s=SpatialFlow(self.g);s.reservoir=np.full(s.shape[1:],.01)
        v,_,D=s.solve(Uniform());expected=2*.01/(.01+np.sqrt(.01**2+2*.002*.01))
        self.assertAlmostEqual(float(v[0].mean()),expected,places=5)
        self.assertLess(D['momentum_relative_residual'],self.g.tolerance)

    def test_spatial_flow_conservation(self):
        d=TreeDrag(self.f,self.g,self.air);v,F,D=SpatialFlow(self.g).solve(d)
        self.assertLess(D['max_divergence_per_m'],1e-11)
        self.assertLess(D['mean_momentum_relative_residual'],1e-4)
        self.assertLess(D['force_ledger']['relative_residual'],1e-12)
        self.assertGreater(D['max_speed_over_reference'],1.) # bypass permitted

    def test_reversed_wind_symmetry(self):
        d=TreeDrag(self.f,self.g,self.air)
        a,_,_=SpatialFlow(self.g,0).solve(d);b,_,_=SpatialFlow(self.g,180).solve(d)
        np.testing.assert_allclose(a[0],-b[0,::-1],atol=2e-5)

    def test_nonconvergence_rejected(self):
        gg=replace(self.g,max_outer=4,tolerance=1e-10)
        with self.assertRaises(RuntimeError):SpatialFlow(gg).solve(TreeDrag(self.f,gg,self.air))

    def test_soil_disjoint(self):
        f,L=share_soil(self.f);np.testing.assert_equal(f,1)
        self.assertAlmostEqual(L['allocated_area_m2'],L['union_area_m2'])

    def test_soil_overlap_not_double_counted(self):
        a=self.f[0];b=replace(a,key='second')
        f,L=share_soil([a,b]);np.testing.assert_equal(f,.5)
        self.assertAlmostEqual(L['allocated_area_m2'],L['union_area_m2'])

    def test_point_load_force_and_moment(self):
        z=np.linspace(0,300,17);trans=np.tile([1.,0.],len(z));rotation=np.empty(2*len(z));rotation[::2]=z;rotation[1::2]=1
        for h in np.linspace(0,300,61):
            s=shape_row(z,h);self.assertAlmostEqual(s@trans,1);self.assertAlmostEqual(s@rotation,h)

    def test_zero_wind_structures(self):
        d=TreeDrag(self.f,self.g,self.air);F=np.zeros((len(d.area),3));R=StandMechanics(self.f,self.air).evaluate(d,F)
        self.assertTrue(all(r['tip_displacement_m']==0 for r in R['trees']))
        self.assertTrue(all(r['root_utilization']==0 for r in R['trees']))

    def test_linear_structure_scaling(self):
        d=TreeDrag(self.f,self.g,self.air);v=np.zeros((3,*self.g.shape));v[0]=1.;F,_=d.forces(v,5.)
        m=StandMechanics(self.f,self.air);a=m.evaluate(d,F);b=m.evaluate(d,4*F)
        np.testing.assert_allclose([r['tip_displacement_m'] for r in b['trees']],4*np.array([r['tip_displacement_m'] for r in a['trees']]),rtol=1e-9)
        self.assertLess(b['horizontal_force_balance_relative_residual'],1e-8)

    def _pair(self):
        a=self.f[0];return [replace(a,key='left',x=-40,y=0),replace(a,key='right',x=40,y=0)]

    def test_contact_transfers_load(self):
        f=self._pair();d=TreeDrag(f,self.g,self.air);F=np.zeros((len(d.area),3))
        row=np.where((d.owner==0)&(d.kind==1))[0][-1];F[row,0]=1e5
        a=StandMechanics(f,self.air).evaluate(d,F)
        b=StandMechanics(f,self.air,Support(crown_mode='contact')).evaluate(d,F)
        self.assertEqual(a['trees'][1]['root_utilization'],0)
        self.assertGreater(b['trees'][1]['root_utilization'],0)
        self.assertLess(b['horizontal_force_balance_relative_residual'],1e-8)
        self.assertGreater(b['active_crown_contacts'],0)

    def test_contact_does_not_pull(self):
        f=self._pair();d=TreeDrag(f,self.g,self.air);F=np.zeros((len(d.area),3));F[np.where((d.owner==0)&(d.kind==1))[0][-1],0]=-1e5
        r=StandMechanics(f,self.air,Support(crown_mode='contact')).evaluate(d,F)
        self.assertEqual(r['active_crown_contacts'],0)

    def test_clonal_root_transfers_moment(self):
        f=self._pair();d=TreeDrag(f,self.g,self.air);F=np.zeros((len(d.area),3));F[np.where((d.owner==0)&(d.kind==1))[0][-1],0]=1e5
        r=StandMechanics(f,self.air,Support(root_link_nm_rad=2e9)).evaluate(d,F)
        self.assertGreater(r['trees'][1]['root_utilization'],0)
        self.assertLess(r['equilibrium_relative_residual'],2e-8)

    def test_broken_link_removes_support(self):
        f=self._pair();s=Support(crown_mode='bonded');m=StandMechanics(f,self.air,s)
        self.assertGreater(len(m.links),0)
        n=StandMechanics(f,self.air,s,broken_links=[x['key'] for x in m.links]);self.assertEqual(len(n.links),0)

    def test_domain_flag(self):
        d=TreeDrag(self.f,self.g,self.air);v=np.ones((3,*self.g.shape));F,_=d.forces(v,80.)
        R=StandMechanics(self.f,self.air).evaluate(d,F)
        self.assertTrue(any(not r['domain_valid'] for r in R['trees']))

    def test_lost_tree_removes_drag(self):
        d=TreeDrag(self.f,self.g,self.air);f=[replace(self.f[0],alive=False),*self.f[1:]];e=TreeDrag(f,self.g,self.air)
        self.assertLess(e.area.sum(),d.area.sum());self.assertFalse(np.any(e.owner==0))

    def test_crown_loss_reduces_drag_and_mass(self):
        a=self.f[0];b=replace(a,crown_remaining=.5)
        d=TreeDrag([a],self.g,self.air);e=TreeDrag([b],self.g,self.air)
        self.assertLess(e.area[e.kind==1].sum(),d.area[d.kind==1].sum())
        self.assertLess(biomass_proxy(b.effective()),biomass_proxy(a.effective()))

    def test_damage_low_wind_control(self):
        R,_=damage_sequence(self.f,self.air,self.g,Support(),[0.,1.])
        self.assertEqual(R['status'],'RAMP_COMPLETED');self.assertEqual(R['events'],[])
        self.assertEqual(R['final_standing_count'],len(self.f))

    def test_immersion_output_guard(self):
        with self.assertRaises(ValueError):output_guard(ROOT/'immersion'/'results')
        with self.assertRaises(ValueError):output_guard(ROOT/'immersion')
        self.assertEqual(output_guard(ROOT/'research'/'runs'),(ROOT/'research'/'runs').resolve())


if __name__=='__main__':unittest.main()
