"""The simple boundary layer's density patch: it applies once to CM1's own turb.F where the driver's block ends and
divides each tendency the block fills by the density at that tendency's own points; only the new build carries it,
the GCM-fed boxes take that build, and they hold their air to the GCM's from above the GCM's convection."""
import pytest

from climate.crm import cm1_run as c

SOURCE = c.CM1_HOME / c.SOURCE['version'] / 'src'
needs_source = pytest.mark.skipif(not (SOURCE / 'turb.F').exists(), reason=f'CM1 is not unpacked in {SOURCE}')


def test_only_the_new_build_carries_the_patch_and_the_site_boxes_take_it():
    assert c.BUILD_PATCHES['moon_omp_rain_law_pbl'] == c.BUILD_PATCHES['moon_omp_rain_law'] + [c.PBL_DENSITY_PATCH]
    assert c.BUILDS['moon_omp_rain_law_pbl'] == c.BUILDS['moon_omp']
    assert c.PBL_DENSITY_PATCH not in c.PATCHES
    # the builds of 28 October stand on it: the wind's own time scale, the gust over water, the pull from aloft and
    # the inflow of the reference's air
    assert all(c.PBL_DENSITY_PATCH not in patches for label, patches in c.BUILD_PATCHES.items()
               if label not in ('moon_omp_rain_law_pbl', 'moon_omp_rain_law_pbl_wind', 'moon_omp_rain_law_pbl_wind_gust',
                                'moon_omp_rain_law_pbl_wind_gust_aloft',
                                'moon_omp_rain_law_pbl_wind_gust_aloft_inflow'))
    for name in ('box_near_side_sea', 'box_0e_design'):
        assert c.CASES[name]['build'] == 'moon_omp_rain_law_pbl'
        assert c.CASES[name]['nudge']['held_from'] == 'convection_top'


@needs_source
def test_the_patch_divides_each_tendency_by_the_density_at_its_own_points_once():
    text = (SOURCE / 'turb.F').read_text(encoding='latin-1')
    patched = c.apply_patch(text, *c.PBL_DENSITY_PATCH)
    assert c.apply_patch(patched, *c.PBL_DENSITY_PATCH) == patched                 # idempotent
    assert patched.count(c.PBL_DENSITY_MARKER) == 1
    block = patched[patched.index('dobr09:'):patched.index('endif  dobr09')]
    for field in ('thpten', 'qvpten', 'qcpten', 'qipten'):
        assert f'{field}(i,j,k) = {field}(i,j,k)*rr(i,j,k)' in block, field
    assert 'upten(i,j,k) = upten(i,j,k)/(0.5*(rho(i-1,j,k)+rho(i,j,k)))' in block
    assert 'vpten(i,j,k) = vpten(i,j,k)/(0.5*(rho(i,j-1,k)+rho(i,j,k)))' in block
    # after every call that fills them, in both of the block's branches
    assert block.count('call     turbvz') == 2 and block.rindex('call     turbvz') < block.index(c.PBL_DENSITY_MARKER)
    # nothing else in the file moves
    assert patched.replace(block, '') == text.replace(text[text.index('dobr09:'):text.index('endif  dobr09')], '')


def test_the_air_is_held_from_above_the_gcms_convection():
    nudge = dict(tau_s=259200.0, z_start_m=8000.0, z_full_m=16000.0, held_from='convection_top')
    out = c.held_heights(nudge, dict(held_from_m=93000.0, held_fully_from_m=129000.0), 150000.0, 1.0)
    assert out == dict(nudge, z_start_m=93000.0, z_full_m=129000.0)
    # held fully from the model's top when the GCM has no level above the first, the heights over the case's scale
    out = c.held_heights(nudge, dict(held_from_m=129000.0, held_fully_from_m=None), 150000.0, 2.0)
    assert out['z_start_m'] == 64500.0 and out['z_full_m'] == 75000.0
    with pytest.raises(ValueError):
        c.held_heights(nudge, dict(held_from_m=None, held_fully_from_m=None), 150000.0, 1.0)
