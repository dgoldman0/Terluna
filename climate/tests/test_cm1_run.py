"""Checks of the CM1 runner's pure pieces: patching, namelist editing, grids, surfaces and restarts."""
import hashlib
import json

import numpy as np
import pytest

from climate.crm import cm1_run as c

GCM_CONFIGURATION = c.HERE.parent / 'gcm' / 'products' / 'moon_gcm_configuration.json'
needs_planet = pytest.mark.skipif(not GCM_CONFIGURATION.exists(), reason='climate/gcm/products/moon_gcm_configuration.json '
                                  'is missing; python -m climate.gcm.boundary writes it')


def test_patch_applies_once_and_is_idempotent():
    text = 'a = 1\nb = 2\nb = 2\n'
    edits = [('', '! head\n', 1), ('a = 1\n', 'a = 1  ! Terluna: marked\n', 1), ('b = 2', 'b = 3', 2)]
    once = c.apply_patch(text, 'demo', 'Terluna: marked', edits)
    assert once == '! head\na = 1  ! Terluna: marked\nb = 3\nb = 3\n'
    assert c.apply_patch(once, 'demo', 'Terluna: marked', edits) == once


def test_patch_refuses_a_different_source():
    with pytest.raises(RuntimeError):
        c.apply_patch('b = 2\n', 'demo', 'Terluna: marked', [('b = 2', 'b = 3  ! Terluna: marked', 2)])
    with pytest.raises(RuntimeError):
        c.apply_patch('b = 2\n', 'demo', 'Terluna: never added', [('b = 2', 'b = 3', 1)])


def test_every_patch_adds_its_marker():
    for name, marker, edits in [*c.MAKEFILE.values(), *c.PATCHES]:
        assert any(marker in new for _, new, _ in edits), name


def test_openmp_builds_compile_at_o3_for_this_processor_unless_given_their_own_flags():
    makefile = '#FC   = gfortran\n#OPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions ' \
               '--param=max-vartrack-size=0 -fopenmp\n#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token ' \
               '-ffreestanding\n#OMP  = -DOPENMP\n.F.o:\n\t$(CPP) $(DM) $(OMP) $(DP) $(ADV) $(OUTPUTOPT) $*.F > $*.f90\n'
    name, marker, edits = c.MAKEFILE['omp']
    base = c.apply_patch(makefile, name, marker, edits)
    for label, flags in c.OPTIMIZATION.items():
        assert c.BUILDS[label] == c.BUILDS['moon_omp']
        name, marker, edits = c.optimization_patch(label, flags)
        fast = c.apply_patch(base, name, marker, edits)
        opts = [line for line in fast.splitlines() if line.startswith('OPTS')]
        assert len(opts) == 1 and f' {flags} -finline-functions ' in opts[0] and (' -O2 ' in opts[0]) == (flags == '-O2')
        assert fast.replace(opts[0], '') == base.replace(c.OMP_OPTS.strip('\n'), '')
    assert c.optimization_for('moon_omp', 'omp') == c.optimization_for('moon_omp_earth_fall', 'omp') == '-O3 -march=native'
    assert c.optimization_for('moon_omp_o2', 'omp') == '-O2' and c.optimization_for('moon', 'mpi') is None


def test_the_rain_test_build_takes_earth_fall_speeds_and_keeps_the_host_gravity_elsewhere():
    source = ("     REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP\n"
              'SUBROUTINE GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst,cm1db)\n         G = 9.81\n'
              '            ACN(K) = G*RHOW/(18.*MU(K))\n           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)\n')
    for old in ('UMS=MIN(UMS,1.2*dum)', 'UNS=MIN(UNS,1.2*dum)', 'UMG=MIN(UMG,20.*dum)', 'UNG=MIN(UNG,20.*dum)'):
        source += (old + '\n') * 3
    for old in ('UMR=MIN(UMR,9.1*dum)', 'UNR=MIN(UNR,9.1*dum)'):
        source += (old + '\n') * 5
    source += 'UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)\nUNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)\n'
    lunar = c.apply_patch(source, *next(p for p in c.PATCHES if p[0] == 'morrison.F'))
    earth = c.apply_patch(lunar, *c.EARTH_FALL_PATCH)
    assert '(G/9.81)' not in earth and earth.count('(GFALL/9.81)') == 8 and 'GFALL = 9.81' in earth
    assert 'ACN(K) = GFALL*RHOW' in earth and 'ALPHA = G*MW' in earth and 'G = terluna_g_host' in earth
    assert c.BUILD_PATCHES['moon_omp_earth_fall'] == [c.EARTH_FALL_PATCH]



def test_the_wind_build_relaxes_the_wind_alone_on_var10_and_keeps_lsnudge_tau_without_it():
    source = ('      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad,terluna_w\n'
              '      if( donudge .and. do_lsnudge_u )then\n'
              '          tem1 = -lsnudgefac*( (umove+uavg(k))-lsnudge_u(k,lsnudge_count) )/(lsnudge_tau)\n'
              '          tem2 = -lsnudgefac*( (vmove+vavg(k))-lsnudge_v(k,lsnudge_count) )/(lsnudge_tau)\n'
              '          tem1 = -lsnudgefac*terluna_w*( thavg(k)-lsnudge_th(k,lsnudge_count) )/(lsnudge_tau)\n')
    out = c.apply_patch(source, *c.WIND_TAU_PATCH)
    assert out.count('/(terluna_tau_uv)') == 2 and out.count('/(lsnudge_tau)') == 1          # the held air keeps its own
    assert 'terluna_tau_uv = lsnudge_tau' in out and 'if( var10.gt.0.0 ) terluna_tau_uv = var10' in out
    assert out.index('terluna_tau_uv = lsnudge_tau') < out.index('if( donudge .and. do_lsnudge_u )then')
    assert c.BUILD_PATCHES['moon_omp_rain_law_pbl_wind'] == c.BUILD_PATCHES['moon_omp_rain_law_pbl'] + [c.WIND_TAU_PATCH]
    assert c.BUILDS['moon_omp_rain_law_pbl_wind'] == c.BUILDS['moon_omp']
    for day in ('day3', 'day4'):
        case, base = c.CASES[f'box_near_side_sea_l15_wind_{day}'], c.CASES[f'box_near_side_sea_l15_{day}']
        assert case['build'] == 'moon_omp_rain_law_pbl_wind' and case['restart_namelist'] == dict(param8=dict(var10=7200.0))
        same = lambda d: {k: v for k, v in d.items() if k not in ('build', 'restart_namelist', 'purpose', 'restart_from')}
        assert same(case) == same(base) and case['restart_from']['day'] == base['restart_from']['day']


def test_the_gust_build_takes_beljaars_gust_over_water_and_changes_nothing_else():
    source = ('!      Use Beljaars over land, old MM5 (Wyngaard) formula over water\n'
              '        if (xland(i).lt.1.5) then\n'
              '        fluxc = max(hfx(i)/rhox(i)/cp                    &\n'
              '              + ep1*tskv*qfx(i)/rhox(i),0.)\n'
              '        VCONV = vconvc*(g/tgdsa(i)*pblh(i)*fluxc)**.33\n'
              '        else\n        VCONV = SQRT(DTHVM)\n        endif\n')
    out = c.apply_patch(source, *c.GUST_PATCH)
    assert 'if (xland(i).lt.1.5 .or. .true.) then' in out and c.GUST_MARKER in out
    assert out.replace('if (xland(i).lt.1.5 .or. .true.) then   ! ' + c.GUST_MARKER + ", at the build's gravity",
                       'if (xland(i).lt.1.5) then') == source
    assert c.BUILD_PATCHES['moon_omp_rain_law_pbl_wind_gust'] == \
        c.BUILD_PATCHES['moon_omp_rain_law_pbl_wind'] + [c.GUST_PATCH]
    assert all(c.GUST_PATCH not in patches for label, patches in c.BUILD_PATCHES.items()
               if label != 'moon_omp_rain_law_pbl_wind_gust')
    for day in ('day3', 'day4'):
        case, wind = c.CASES[f'box_near_side_sea_l15_wind_gust_{day}'], c.CASES[f'box_near_side_sea_l15_wind_{day}']
        assert case['build'] == 'moon_omp_rain_law_pbl_wind_gust'
        same = lambda d: {k: v for k, v in d.items() if k not in ('build', 'purpose', 'restart_from')}
        assert same(case) == same(wind) and case['restart_from']['day'] == wind['restart_from']['day']
    assert c.CASES['box_near_side_sea_l15_wind_gust_day4']['restart_from']['case'] == 'box_near_side_sea_l15_wind_gust_day3'

def test_the_twilight_build_adds_its_hook_after_the_sun_and_only_there():
    name, marker, edits = c.TWILIGHT_PATCH
    assert name == 'radiation_driver.F' and any(marker in new for _, new, _ in edits)
    assert c.BUILD_PATCHES['moon_omp_twilight'] == [c.TWILIGHT_PATCH] and c.BUILDS['moon_omp_twilight'] == c.BUILDS['moon_omp']
    # the rings', Jura's box, and the boxes on the design climate with the rain law (and the boundary layer's fix, the
    # wind's own nudging time scale and the gust over water at lunar gravity)
    with_hook = {'moon_omp_twilight', 'moon_omp_terrain_sea_twilight', 'moon_omp_rain_law', 'moon_omp_rain_law_pbl',
                 'moon_omp_rain_law_pbl_wind', 'moon_omp_rain_law_pbl_wind_gust'}
    assert all(c.TWILIGHT_PATCH not in patches for label, patches in c.BUILD_PATCHES.items() if label not in with_hook)
    sun = next(p for p in c.PATCHES if p[0] == 'radiation_driver.F')
    added = ''.join(new for _, new, _ in sun[2])
    # It extends the lines the Sun's patch writes, which the build applies first.
    for old in ('use input, only : var12,var16,var17,var18,var19', 'real :: terluna_lat,terluna_s'):
        assert old in added and any(old in o for o, _, _ in edits)
    for tag in ('45n', '45s'):
        case = c.CASES[f'ring_{tag}_lsw_twilight']
        assert case['twilight'] in c.TWILIGHT_MODES and case['build'] == 'moon_omp_twilight'
        assert {k: v for k, v in case.items() if k not in ('twilight', 'build', 'purpose')} == \
               {k: v for k, v in c.CASES[f'ring_{tag}_lsw'].items() if k not in ('build', 'purpose')}


def test_the_corrected_rings_take_their_own_gcm_run_with_and_without_the_twilight():
    run, years = c.CORRECTED_GCM
    assert run != c.GCM_RUN and years[1] - years[0] == 9
    for tag in ('45n', '45s'):
        plain, lit = c.CASES[f'ring_{tag}_corrected'], c.CASES[f'ring_{tag}_corrected_twilight']
        assert plain['gcm'] == lit['gcm'] == c.CORRECTED_GCM and 'gcm' not in c.CASES[f'ring_{tag}_lsw']
        assert plain['build'] == 'moon_omp' and 'twilight' not in plain
        assert lit['build'] == 'moon_omp_twilight' and lit['twilight'] == 'full'
        # Apart from the GCM run, the build and the twilight they are the earlier 45-degree rings.
        same = lambda case: {k: v for k, v in case.items() if k not in ('gcm', 'twilight', 'build', 'purpose')}
        assert same(plain) == same(lit) == same(c.CASES[f'ring_{tag}_lsw'])
        assert '3% dimmer shield' in plain['purpose'] and '5% dimmer' not in lit['purpose']


def test_the_twilight_table_reads_as_the_hook_expects(tmp_path):
    product = c.HERE.parent / 'results' / 'gcm' / 'twilight.json'
    if not product.exists():
        pytest.skip('twilight.json not built')
    tables = {}
    for mode in c.TWILIGHT_MODES:
        record = c.write_twilight(tmp_path / c.TWILIGHT_TABLE, mode)
        rows = (tmp_path / c.TWILIGHT_TABLE).read_text().splitlines()
        n, shells = (int(x) for x in rows[0].split())
        assert (n, shells) == (record['elevations'], record['shells']) and len(rows) == n + 4
        edges = np.array(rows[1].split(), float)
        assert edges.size == shells + 1 and edges[0] == 0.0 and np.all(np.diff(edges) > 0)
        table = np.array([r.split() for r in rows[2:n + 2]], float)
        assert table.shape == (n, 2 + shells) and np.all(np.diff(table[:, 0]) > 0)
        escape, temperature = (np.array(r.split(), float) for r in rows[n + 2:])
        assert escape.size == temperature.size == shells + 1 and np.all(temperature > 100.0)
        tables[mode] = table, escape
    # The full table adds the shells' widening, which takes more sunlight and lets out more infrared.
    assert tables['full'][1].sum() > tables['curved'][1].sum() > 0.0
    with pytest.raises(ValueError):
        c.write_twilight(tmp_path / c.TWILIGHT_TABLE, 'both')


def test_a_case_set_up_from_another_takes_its_inputs_and_only_a_new_executable(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'RUNS', tmp_path / 'runs')
    monkeypatch.setattr(c, 'CM1_HOME', tmp_path / 'cm1')
    source = tmp_path / 'runs' / 'box_0e_small'
    source.mkdir(parents=True)
    for name in ('input_sounding', 'namelist.template', 'terluna_surface.txt', 'progress.json', 'namelist.input',
                 'cm1out_t000001_s.dat', 'cm1rst_t000002_s.dat', 'cm1_segment_001.log'):
        (source / name).write_text(name)
    (source / 'RRTMG_LW_DATA').symlink_to(source / 'input_sounding')
    (source / 'case.json').write_text(json.dumps(dict(case='box_0e_small', build=dict(label='moon_omp', executable_sha256='old'))))
    build = tmp_path / 'cm1' / 'build' / 'moon_omp_earth_fall'
    build.mkdir(parents=True)
    (build / 'cm1.exe').write_text('exe')
    (build / 'build.json').write_text(json.dumps(dict(label='moon_omp_earth_fall', executable_sha256='new')))
    case = c.setup('box_0e_small_earth_fall')
    assert sorted(p.name for p in case.iterdir()) == ['RRTMG_LW_DATA', 'case.json', 'cm1.exe', 'input_sounding',
                                                      'namelist.template', 'terluna_surface.txt']
    assert (case / 'cm1.exe').resolve() == (build / 'cm1.exe').resolve() and (case / 'RRTMG_LW_DATA').is_symlink()
    record = json.loads((case / 'case.json').read_text())
    assert record['case'] == 'box_0e_small_earth_fall' and record['build']['label'] == 'moon_omp_earth_fall'
    assert record['inputs_from'] == dict(case='box_0e_small', build='moon_omp', executable_sha256='old')


def test_morrison_fall_speeds_follow_gravity():
    # V ~ g^((B+1)/3): Stokes droplets (B = 2) fall g times slower, hail (B = 0.5) as sqrt(g)
    assert (1 / 6) ** ((2 + 1) / 3) == pytest.approx(1 / 6)
    assert (1 / 6) ** ((0.5 + 1) / 3) == pytest.approx((1 / 6) ** 0.5)
    assert 'AR = AR*(G/9.81)**((BR+1.)/3.)' in c.MORRISON_BLOCK


def test_set_namelist_changes_only_the_named_entry():
    text = ' &param1\n dx     =  3000.0,\n dx_inner = 5.0,\n /\n &param10\n dx     =  1.0,\n /\n'
    out = c.set_namelist(text, 'param1', 'dx', 6011.2)
    assert ' dx     =  6011.2,' in out and ' dx_inner = 5.0,' in out and out.endswith(' dx     =  1.0,\n /\n')
    assert '.true.' in c.set_namelist(' &p\n flag = .false.,\n /\n', 'p', 'flag', True)
    with pytest.raises(KeyError):
        c.set_namelist(text, 'param1', 'dy', 1.0)


def test_vertical_grid_is_fine_below_and_even_aloft():
    z = np.array(c.vertical_grid(dz0=100.0, z_fine=1000.0, stretch=1.08, dz_max=2000.0, ztop=150000.0))
    dz = np.diff(z)
    assert z[0] == 0.0 and z[-1] == pytest.approx(150000.0)
    assert np.all(dz > 0) and dz[:10] == pytest.approx(100.0)
    assert np.all(dz[1:] / dz[:-1] < 1.081) and dz.max() <= 2000.0 * 1.5


def test_ring_surface_segments_cover_the_ring(monkeypatch):
    lon = (np.arange(1440) + 0.5) * 0.25
    lat = 89.875 - 0.25 * np.arange(720)
    water = np.zeros((720, 1440))
    water[:, (lon > 90) & (lon < 180)] = 1.0                           # a quarter of the equator is sea
    from climate.gcm import boundary
    monkeypatch.setattr(boundary, 'lakes_product', lambda *a: dict(lat_deg=lat, lon_deg=lon, water_fraction=water,
                                                                   sha256='test'))
    segments, share, _ = c.ring_surface(nx=360, dx=1000.0, sea_k=299.0, land_k=295.0)
    assert share == pytest.approx(0.25, abs=0.01)
    assert segments[0][0] < 0 and segments[-1][1] > 360 * 1000.0
    assert [s[2] for s in segments] == [1, 2, 1]
    assert all(a[1] == b[0] for a, b in zip(segments, segments[1:]))


def test_latest_restart_needs_a_complete_set(tmp_path):
    for part in 'isuvwx':
        (tmp_path / f'cm1rst_t000001_{part}.dat').write_bytes(b'')
    for part in 'isuv':
        (tmp_path / f'cm1rst_t000002_{part}.dat').write_bytes(b'')
    assert c.latest_restart(tmp_path) == 1
    empty = tmp_path / 'empty'
    empty.mkdir()
    assert c.latest_restart(empty) == 0


def test_a_run_on_from_a_restart_carries_the_wind_windows_on_by_whole_lunar_days(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'solar_day_s', lambda: 3600.0)               # 36 windows of 100 s a day
    source, case = tmp_path / 'source', tmp_path / 'case'
    source.mkdir()
    case.mkdir()
    for i in range(72):                                                   # two days, each phase's profile repeating
        t2 = 1.0e30 if i == 71 else (i + 1) * 100.0
        (source / f'lsnudge_{i + 1:04d}.dat').write_text(
            f' *  Header\n  {i * 100.0:22.1f} {t2:22.6e}\n *  Profile\n *  z\n  0.0 300.0 10.0 {i % 36:.1f} 0.0\n')
    record = c.restart_windows(source, case, 7150.0, 10650.0)            # from inside the 72nd window
    files = sorted(case.glob('lsnudge_*.dat'))
    assert record['windows'] == len(files) == 36 and record['first_window_of_source'] == 72
    first, second, last = (f.read_text().splitlines() for f in (files[0], files[1], files[-1]))
    assert first[1].split() == ['7100.0', '7.200000e+03'] and first[4].split()[3] == '35.0'
    assert second[1].split()[0] == '7200.0' and second[4].split()[3] == '0.0'      # the next day's first phase
    assert last[1].split()[1] == '1.000000e+30' and last[4].split()[3] == '34.0'
    (source / 'lsnudge_0040.dat').write_text(' *  Header\n  3900.0 4.0e+03\n *  Profile\n *  z\n  0.0 300.0 10.0 9.0 0.0\n')
    with pytest.raises(RuntimeError):                                     # a day's windows that disagree
        c.restart_windows(source, case, 7150.0, 10650.0)
    single = tmp_path / 'single'
    single.mkdir()
    (single / 'lsnudge_0001.dat').write_text(' *  Header\n  0.0 1.0e+30\n')
    assert c.restart_windows(single, case, 0.0, 100.0) is None


def test_a_run_on_from_a_run_on_case_takes_its_windows_from_the_first_case(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'solar_day_s', lambda: 3600.0)               # 36 windows of 100 s a day
    root, third, fourth = (tmp_path / n for n in ('root', 'third', 'fourth'))
    for folder in (root, third, fourth):
        folder.mkdir()
    for i in range(72):
        t2 = 1.0e30 if i == 71 else (i + 1) * 100.0
        (root / f'lsnudge_{i + 1:04d}.dat').write_text(
            f' *  Header\n  {i * 100.0:22.1f} {t2:22.6e}\n *  Profile\n *  z\n  0.0 300.0 10.0 {i % 36:.1f} 0.0\n')
    (root / 'case.json').write_text(json.dumps(dict(case='root')))
    (third / 'case.json').write_text(json.dumps(dict(case='third', restart_from=dict(case='root', day=0.083))))
    c.restart_windows(root, third, 7150.0, 10650.0)                       # the third day's windows, from phase 35
    assert c.window_root(third) == root and c.window_root(root) == root
    c.restart_windows(c.window_root(third), fourth, 10800.0, 14400.0)    # the fourth day from the third's restart
    first = sorted(fourth.glob('lsnudge_*.dat'))[0].read_text().splitlines()
    assert first[1].split()[0] == '10800.0' and first[4].split()[3] == '0.0'      # the day's first phase
    c.restart_windows(third, fourth, 10800.0, 14400.0)                    # the run-on case's own windows
    shifted = sorted(fourth.glob('lsnudge_*.dat'))[0].read_text().splitlines()
    assert shifted[4].split()[3] == '35.0'                                # a window out of phase


def test_storms_wrap_around_the_ring():
    from climate.crm import ring_analysis as ra
    rain = np.zeros(20)
    rain[[0, 1, 18, 19]] = 5.0                                         # one storm across the seam
    rain[[8, 9, 10]] = 2.0
    found = sorted(ra.storms(rain, threshold=1.0))
    assert found == [(8, 3), (18, 4)]
    assert ra.storms(np.full(5, 3.0)) == [(0, 5)] and ra.storms(np.zeros(5)) == []


def test_hour_angle_grows_eastward_and_with_time():
    from climate.crm import ring_analysis as ra
    geo = dict(x=np.array([0.0, 250.0, 500.0]), sun=dict(h0=-90.0, length=1000.0))
    assert ra.hour_angle(geo, 0.0, 100.0) == pytest.approx([-90.0, 0.0, 90.0])
    assert ra.hour_angle(geo, 25.0, 100.0) == pytest.approx([0.0, 90.0, -180.0])


def test_read_snapshot_follows_the_grads_layout(tmp_path):
    from climate.crm import ring_analysis as ra
    ctl = ('dset ^cm1out_t00%y4_s.dat\nxdef 3 linear 0 1\nydef 1 linear 0 1\nzdef 2 levels\n1\n2\nvars 2\n'
           'rain 0 99 rain\nth 2 99 theta\nendvars\n')
    (tmp_path / 'cm1out_s.ctl').write_text(ctl)
    np.arange(9, dtype='<f4').tofile(tmp_path / 'cm1out_t000001_s.dat')
    d = ra.read_snapshot(tmp_path, 1)
    assert d['rain'].tolist() == [0, 1, 2] and d['th'].tolist() == [[3, 4, 5], [6, 7, 8]]


def test_wet_bulb_matches_the_psychrometric_chart():
    from climate.crm import ring_analysis as ra
    t = 303.15
    assert ra.wet_bulb(t, 0.5 * ra.saturation_pa(t), 101325.0) - 273.15 == pytest.approx(22.0, abs=0.2)
    assert ra.wet_bulb(298.15, ra.saturation_pa(298.15), 121600.0) == pytest.approx(298.15, abs=1e-6)
    assert ra.vapour_pa(0.622, 2000.0) == pytest.approx(1000.0)                  # r = 0.622 is half the air


def test_hour_angle_table_keeps_land_and_water_apart():
    from climate.crm import ring_analysis as ra
    bins = ra.HOUR_BINS
    comp = {k: np.zeros((2, bins)) for k in ('t2', 'tsk', 'rh2', 'tw2', 'prate', 'evap', 's10', 'hpbl', 'lcl', 'cape',
                                              'cloud', 'cloud_high', 'cloud_flight', 'fog', 'cloud_top_m', 'pwat')}
    comp['t2'][0], comp['t2'][1] = 300.0, 290.0
    comp['wind_levels'] = np.ones((2, bins, len(ra.WIND_HEIGHTS_M)))
    table = ra.hour_angle_table(np.linspace(-175.0, 175.0, bins), comp)
    assert table['land']['air_2m_c'][0] == pytest.approx(26.9) and table['water']['air_2m_c'][0] == pytest.approx(16.9)
    assert len(table['land']['wind_by_height_m_s']) == len(ra.WIND_HEIGHTS_M) == len(table['wind_heights_km'])
    assert table['land']['fog'] == [0.0] * bins


def test_storm_tracks_follow_overlap_and_measure_drift():
    from climate.crm import ring_analysis as ra
    found = [[(10, 4)], [(12, 4), (40, 3)], [(14, 4)], []]            # one storm drifting east, one brief
    tracks = ra.track_storms(found, nx=50, reach=1)
    assert sorted(len(t) for t in tracks) == [1, 3]
    summary = ra.track_summary(tracks, output_s=3600.0, dx=1000.0, nx=50)
    assert summary['drift_m_s']['median'] == pytest.approx(4 * 1000.0 / (2 * 3600.0))   # 4 columns in 2 hours
    assert summary['lifetime_h']['max'] == pytest.approx(3.0) and summary['seen_once_share'] == pytest.approx(0.5)
    across = ra.track_storms([[(48, 3)], [(0, 3)]], nx=50, reach=1)   # across the seam of the ring
    assert len(across) == 1 and ra.track_summary(across, 3600.0, 1000.0, 50)['drift_m_s']['median'] > 0


@needs_planet
def test_rain_spells_count_runs_at_each_place():
    from climate.crm import ring_analysis as ra
    prate = np.zeros((6, 2))
    prate[[1, 2], 0] = 5.0                                              # one 2-snapshot spell on land
    prate[[0, 3, 5], 1] = 5.0                                           # three single spells on water
    out = ra.rain_spells(prate, np.array([True, False]), output_s=3600.0, threshold=1.0)
    assert out['land']['duration_h']['median'] == pytest.approx(2.0)
    assert out['water']['spells_per_lunar_day'] is None                # six hours is not a lunar day
    month = ra.rain_spells(prate, np.array([True, False]), output_s=29.530589 * 86400.0 / 6, threshold=1.0)
    assert month['water']['spells_per_lunar_day'] == pytest.approx(3.0, rel=1e-3)


def test_a_case_ends_at_its_last_whole_restart_interval():
    two_lunar_days = dict(days=59.06, restart_s=43200.0)
    assert c.case_length_s(two_lunar_days, 1.0) == pytest.approx(59.0 * 86400.0)
    scaled = dict(days=20.0, restart_s=86400.0)                         # the Earth-gravity pair case
    s = 0.16556767345166287
    assert c.case_length_s(scaled, s) == pytest.approx(20.0 * 86400.0 * s)


def test_noon_lies_where_the_hour_angle_is_zero():
    from climate.crm import ring_analysis as ra
    geo = dict(x=np.array([0.0]), sun=dict(h0=-90.0, length=1000.0))
    assert ra.noon_x_km(geo, 0.0, 100.0) == pytest.approx(0.25)
    assert ra.noon_x_km(geo, 25.0, 100.0) == pytest.approx(0.0)                 # a quarter day later noon reached x = 0
    x = ra.noon_x_km(geo, 60.0, 100.0) * 1000.0
    assert ra.hour_angle(dict(x=np.array([x]), sun=geo['sun']), 60.0, 100.0)[0] == pytest.approx(0.0, abs=1e-9)


def test_field_product_pools_rain_and_keeps_sections_to_90_km(tmp_path):
    import json
    from climate.crm import ring_analysis as ra
    nz, nx, nt, w = 6, 8, 3, 2 * ra.STORM_HALF_WIDTH
    heights = np.array([1.0, 20.0, 60.0, 89.0, 95.0, 120.0]) * 1000.0
    rain = np.zeros((nt, nx)); rain[1, 5] = 7.0
    arrays = dict(heights_m=heights, hour_angle_deg=np.array([-5.0, 5.0]), x_km=np.arange(nx) * 6.0,
                  section_u=np.ones((2, nz)), section_w=np.zeros((2, nz)), section_cloud=np.zeros((2, nz)),
                  hov_prate_mm_h=rain, hov_time_day=np.arange(nt) * 0.125, hov_noon_x_km=np.zeros(nt),
                  land=np.array([1, 1, 0, 0, 1, 0, 1, 1], bool), storm_day=30.0, storm_hour_angle_deg=45.0,
                  storm_rain_max_mm_h=7.0, storm_condensate_g_kg=np.ones((nz, w)), storm_rain_water_g_kg=np.zeros((nz, w)),
                  storm_updraft_m_s=np.zeros((nz, w)), storm_rain_mm_h=np.zeros(w), storm_wind_10m_m_s=np.zeros(w),
                  storm_air_2m_c=np.zeros(w), storm_land=np.ones(w, bool))
    path = tmp_path / 'fields.npz'
    ra.write_fields(path, dict(case='t', span_days=[0, 1], evidence='e'), arrays)
    f = np.load(path)
    assert json.loads(str(f['metadata']))['schema'] == ra.FIELDS_SCHEMA
    assert f['section_eastward_m_s'].shape == (2, 4) and f['storm_condensate_g_kg'].shape == (5, w)
    assert f['rain_map_max_mm_h'].shape == (nt, 2) and float(f['rain_map_max_mm_h'][1, 1]) == 7.0
    assert f['rain_map_land_share'].tolist() == [0.5, 0.75]


def test_water_mean_profile_skips_land_columns():
    land = np.array([True, False, False])
    snap = dict(th=np.array([[300.0, 290.0, 292.0]]), qv=np.array([[0.02, 0.01, 0.012]]),
                psfc=np.array([1.0e5, 1.2e5, 1.22e5]), t2=np.array([305.0, 298.0, 300.0]), q2=np.array([0.02, 0.014, 0.016]))
    prof = c.water_mean_profile([snap, snap], land)
    assert prof['theta_k'][0] == pytest.approx(291.0) and prof['qv_kg_kg'][0] == pytest.approx(0.011)
    assert prof['surface_pa'] == pytest.approx(1.21e5) and prof['air_2m_k'] == pytest.approx(299.0)


@needs_planet
def test_rings_off_the_equator_are_shorter_and_turn_the_wind():
    nx, dx, length = c.ring_grid(80.0, 6000.0, 8)
    assert length == pytest.approx(2 * np.pi * c.planet()['radius_m'] * np.cos(np.radians(80.0)))
    assert nx % 8 == 0 and nx * dx == pytest.approx(length) and abs(dx - 6000.0) < 400.0
    omega = c.planet()['rotation_rate_rad_s']
    assert c.coriolis(80.0) == pytest.approx(2 * omega * np.sin(np.radians(80.0)))
    assert c.coriolis(-80.0) == pytest.approx(-c.coriolis(80.0)) and c.coriolis(0.0) == 0.0


@needs_planet
def test_only_rings_off_the_equator_turn_the_wind_and_lower_the_sun():
    air = dict(sunlight_w_m2=1300.0)
    south = c.case_settings(c.CASES['ring_80s'], 312, 6075.7, 111, 150000.0, air, 1.6242)
    assert south['param2']['icor'] == 1 and south['param2']['lspgrad'] == 1 and south['param3']['fcor'] < 0
    assert south['param11']['ctrlat'] == -80.0 and south['param19']['do_lsnudge_v'] and south['param9']['output_vinterp'] == 1
    equator = c.case_settings(c.CASES['ring'], 1816, 6011.0, 111, 150000.0, air, 1.6242)
    assert equator['param2']['icor'] == 0 and equator['param2']['lspgrad'] == 0 and equator['param3']['fcor'] == 0.0
    assert equator['param11']['ctrlat'] == 0.0 and not equator['param19']['do_lsnudge_v']
    assert equator['param9']['output_vinterp'] == 0
    assert equator['param8']['var13'] == 0.0 and south['param8']['var13'] == 0.0
    sinking = c.case_settings(c.CASES['ring_70s_lsw'], 624, 5983.4, 111, 150000.0, air, 1.6242)
    assert sinking['param8']['var13'] == 1.0 and sinking['param11']['ctrlat'] == -70.0


def test_ring_surface_off_the_equator_averages_wide_columns(monkeypatch):
    lon = (np.arange(1440) + 0.5) * 0.25
    lat = 89.875 - 0.25 * np.arange(720)
    water = np.zeros((720, 1440))
    water[np.ix_(np.abs(lat - 80.0) < 1.0, lon < 90.0)] = 1.0               # a quarter of 80 N is sea
    water[np.ix_(np.abs(lat) < 1.0, lon < 180.0)] = 1.0                     # half of the equator, which 80 N ignores
    water[np.ix_(np.abs(lat - 80.0) < 1.0, (lon > 180.0) & (lon < 180.5))] = 1.0   # a sliver narrower than a column
    from climate.gcm import boundary
    monkeypatch.setattr(boundary, 'lakes_product', lambda *a: dict(lat_deg=lat, lon_deg=lon, water_fraction=water,
                                                                   sha256='test'))
    segments, share, _ = c.ring_surface(nx=312, dx=6075.7, sea_k=293.0, land_k=291.0, latitude=80.0)
    assert share == pytest.approx(0.25, abs=0.01)
    assert [s[2] for s in segments] == [2, 1]


def test_land_evaporates_as_the_gcm_bucket_allows(tmp_path):
    assert c.gcm_soil(tmp_path) == c.PLASIM_SOIL
    (tmp_path / 'landmod_namelist').write_text(' &landmod_nl\n WSMAX = 1.0\n /END\n')
    assert c.gcm_soil(tmp_path) == dict(wsmax=1.0, drhsfull=0.4)
    assert c.land_moisture(dict(land_moisture='gcm'), dict(land_wetness=0.03)) == pytest.approx(0.03)
    assert c.land_moisture(dict(land_moisture=0.5), dict(land_wetness=0.03)) == pytest.approx(0.5)


def test_circulation_reports_the_turned_wind_only_where_there_is_one():
    from climate.crm import ring_analysis as ra
    zh = np.array(ra.WIND_HEIGHTS_M)
    section = dict(u=np.ones((2, zh.size)), w=np.zeros((2, zh.size)), cloud=np.zeros((2, zh.size)))
    assert 'northward_m_s' not in ra.circulation_table(np.array([-5.0, 5.0]), zh, section)
    section['v'] = np.full((2, zh.size), -2.0)
    table = ra.circulation_table(np.array([-5.0, 5.0]), zh, section)
    assert table['northward_m_s'][0] == [-2.0, -2.0]
    assert 'S' in ra.evidence(-80.0) and ra.evidence(0.0) == ra.EVIDENCE


def test_plasim_half_levels_put_each_full_level_midway():
    half = np.array([0.0, 0.2, 0.5, 1.0])
    assert c.plasim_half_levels(0.5 * (half[1:] + half[:-1])) == pytest.approx(half)


def test_band_vertical_wind_sinks_where_air_converges_aloft():
    lat = np.array([85.0, 80.0, 75.0, 70.0])                             # rows north to south
    edges = np.array([0.0, 0.5, 1.0])
    flux = np.array([[3.0] * 4, [-1.0] * 4])                             # poleward aloft, equatorward below
    radius, gravity = 1.0e6, 1.6
    w, (south, north) = c.band_vertical_wind(lat, edges, flux, lat == 80.0, radius, gravity)
    assert (south, north) == (77.5, 82.5)
    above = 2.0 * 0.5 / gravity                                          # the column mean (1.0) removed: 3 - 1 = 2
    into = 2 * np.pi * radius * (np.cos(np.radians(77.5)) - np.cos(np.radians(82.5))) * above
    area = 2 * np.pi * radius ** 2 * (np.sin(np.radians(82.5)) - np.sin(np.radians(77.5)))
    assert w[1] == pytest.approx(-into / area) and w[1] < 0             # converging aloft, so it sinks
    assert w[0] == 0.0 and w[2] == pytest.approx(0.0, abs=1e-15)


def test_vertical_wind_file_ends_at_zero_at_the_model_top(tmp_path):
    path = tmp_path / 'terluna_wls.txt'
    c.write_vertical_wind(path, dict(z_m=[0.0, 2000.0, 60000.0], w_m_s=[0.0, -0.003, -0.002]), 150000.0)
    lines = path.read_text().splitlines()
    assert lines[0] == '4' and lines[-1].split() == ['150000.000', '0.000000e+00']
    assert [float(l.split()[1]) for l in lines[1:3]] == [0.0, -0.003]


def test_tilted_path_crosses_the_equator_at_its_nodes_and_mirrors():
    nx = 1816
    a, a_prime = c.tilted_path(nx, 45.0, 0.0), c.tilted_path(nx, 45.0, 180.0)
    assert abs(a['lat'][0]) < 0.2 and a['lon'][0] < 0.2                     # heading north across 0 E
    assert a['lat'].max() == pytest.approx(45.0, abs=0.01) and a['lat'].min() == pytest.approx(-45.0, abs=0.01)
    assert a['lon'][np.argmax(a['lat'])] == pytest.approx(90.0, abs=0.2)
    assert np.all(np.diff(a['dlon']) > 0) and a['dlon'][-1] == pytest.approx(360.0, abs=0.3)
    # A-prime is A mirrored across the equator: at the same longitude, the opposite latitude
    lat_mirror = np.interp(a['lon'], np.sort(a_prime['lon']), a_prime['lat'][np.argsort(a_prime['lon'])])
    assert np.allclose(lat_mirror, -a['lat'], atol=0.1)
    # the Coriolis parameter CM1 builds along the ring is that of each column's latitude
    s = 2 * np.pi * (np.arange(nx) + 0.5) / nx
    omega2 = 2 * 2.6617e-6
    assert np.allclose(omega2 * np.sin(np.radians(45.0)) * np.sin(s), omega2 * np.sin(np.radians(a['lat'])))


def test_wetness_classes_and_their_land_use_rows():
    assert c.wetness_class(0.04) == 20 and c.wetness_class(0.47) == 27 and c.wetness_class(1.0) == 30
    assert [c.wetness_class(w) for w in (0.001, 0.005, 0.013, 0.018)] == [31, 31, 32, 33]   # ground drier than 0.05
    rows = c.landuse_rows(((20, 0.05), (30, 0.90), (31, 0.003)))
    assert rows[20].startswith('20,') and '   .05,' in rows[20] and 'wetness 0.05' in rows[20]
    assert rows[30].startswith('30,') and '   .90,' in rows[30]
    assert rows[31].startswith('31,') and '   .003,' in rows[31] and 'wetness 0.003' in rows[31]
    indices = [i for i, _ in c.WETNESS_CLASSES]
    assert len(set(indices)) == len(indices) and not {16, 24} & set(indices)   # water and snow keep their rows


def test_path_surface_samples_the_product_along_the_path(monkeypatch):
    lon = (np.arange(1440) + 0.5) * 0.25
    lat = 89.875 - 0.25 * np.arange(720)
    water = np.zeros((720, 1440))
    water[:, lon < 90.0] = 1.0                                               # water west of 90 E
    from climate.gcm import boundary
    monkeypatch.setattr(boundary, 'lakes_product', lambda *a: dict(lat_deg=lat, lon_deg=lon, water_fraction=water, sha256='t'))
    path = c.tilted_path(64, 45.0, 0.0)
    n = len(path['lat'])
    segs, share, _ = c.path_surface(path['lat'], path['lon'], 1000.0, np.full(n, 297.0), np.full(n, 293.0), np.full(n, 0.4))
    wet = np.array([s[2] == 2 for s in segs])
    assert np.array_equal(wet, path['lon'] < 90.0)
    assert all(s[3] == 16 for s in segs if s[2] == 2) and all(s[3] == 26 for s in segs if s[2] == 1)
    assert segs[0][0] < 0 and segs[-1][1] > 64 * 1000.0 and share == pytest.approx(wet.mean())


def test_path_heights_keep_water_at_its_level_and_smooth_the_land(monkeypatch):
    lon = (np.arange(1440) + 0.5) * 0.25
    lat = 89.875 - 0.25 * np.arange(720)
    level = -1000.0
    water = np.zeros((720, 1440))
    water[:, lon < 90.0] = 1.0                                               # sea west of 90 E
    water[:, (lon > 180.0) & (lon < 200.0)] = 1.0                           # a raised lake
    height = np.where(lon < 90.0, level - 500.0, level + 3000.0)[None, :] * np.ones((720, 1))
    surface = np.where((lon > 180.0) & (lon < 200.0), level + 2000.0, level)[None, :] * np.ones((720, 1))
    from climate.gcm import boundary
    monkeypatch.setattr(boundary, 'lakes_product', lambda *a: dict(lat_deg=lat, lon_deg=lon, water_fraction=water,
                                                                    mean_height_m=height, water_surface_m=surface,
                                                                    meta=dict(level_m=level), sha256='t'))
    path = c.tilted_path(360, 0.0, 0.0)                                      # the equator, one column a degree
    z = c.path_heights(path['lat'], path['lon'], passes=2)
    sea, lake = path['lon'] < 89.5, (path['lon'] > 180.5) & (path['lon'] < 199.5)
    assert np.all(z[sea] == 0.0) and np.allclose(z[lake], 2000.0)            # water at its own level
    inland = (path['lon'] > 95.0) & (path['lon'] < 175.0)
    assert np.allclose(z[inland], 3000.0) and np.all(z >= 0.0)
    coast = np.argmin(np.abs(path['lon'] - 90.5))
    assert 0.0 < z[coast] < 3000.0                                           # the land's edge is smoothed


@needs_planet
def test_the_terrain_ring_is_the_flat_ring_with_ground_heights():
    air = dict(sunlight_w_m2=1300.0)
    flat, rough = (c.case_settings(c.CASES[n], 1816, 6011.0, 111, 150000.0, air, 1.6242)
                   for n in ('ring_70_45e', 'ring_70_45e_terrain'))
    assert rough['param0']['terrain_flag'] and rough['param2']['itern'] == 4 and rough['param2']['irdamp'] == 1
    assert not flat['param0']['terrain_flag'] and flat['param2']['itern'] == 0 and flat['param2']['irdamp'] == 2
    differ = {(s, k) for s in flat for k in flat[s] if flat[s][k] != rough[s].get(k)}
    assert differ == {('param0', 'terrain_flag'), ('param2', 'itern'), ('param2', 'irdamp')}
    assert c.CASES['ring_70_45e_terrain']['build'] == 'moon_omp_terrain' and 'moon_omp_terrain' in c.BUILD_PATCHES
    patched = {name for name, _, _ in c.BUILD_PATCHES['moon_omp_terrain']}
    assert patched == {'param.F', 'lsnudge.F', 'solve1.F', 'init_surface.F',   # every column's ground, for boxes
                       'input.F', 'base.F', 'adv_routines.F', 'adv.F'}      # and its own wind


def test_vertical_wind_table_has_one_row_per_column(tmp_path):
    path = tmp_path / 'terluna_wls2d.txt'
    c.write_vertical_wind_2d(path, [0.0, 1000.0, 2000.0], [[0.0, 0.002, 0.0], [0.0, -0.003, 0.0]])
    lines = path.read_text().splitlines()
    assert lines[0] == '2 3' and lines[1].split() == ['0.000', '1000.000', '2000.000']
    assert float(lines[3].split()[1]) == -0.003


@needs_planet
def test_tilted_rings_use_the_beta_plane_and_their_own_sun():
    air = dict(sunlight_w_m2=1300.0)
    a_prime = c.case_settings(c.CASES['ring_a_prime'], 1816, 6011.0, 111, 150000.0, air, 1.6242)
    assert a_prime['param2']['betaplane'] == 1 and a_prime['param2']['icor'] == 1 and a_prime['param2']['lspgrad'] == 0
    assert a_prime['param8']['var12'] == 45.0 and a_prime['param8']['var13'] == 2.0
    assert a_prime['param8']['var18'] == -90.0 + 180.0                     # local time follows longitude on every ring
    assert a_prime['param3']['fcor'] == pytest.approx(c.coriolis(90.0), rel=1e-5)
    equator = c.case_settings(c.CASES['ring'], 1816, 6011.0, 111, 150000.0, air, 1.6242)
    assert equator['param2']['betaplane'] == 0 and equator['param8']['var12'] == 0.0


@needs_planet
def test_a_great_circle_of_no_tilt_is_the_equator_with_no_coriolis_force():
    air = dict(sunlight_w_m2=1300.0)
    cfg = c.CASES['ring_equator']
    equator = c.case_settings(cfg, 1816, 6011.0, 111, 150000.0, air, 1.6242)
    assert equator['param2']['icor'] == 0 and equator['param2']['betaplane'] == 0 and equator['param3']['fcor'] == 0.0
    assert not equator['param19']['do_lsnudge_v'] and equator['param11']['ctrlat'] == 0.0
    assert equator['param8']['var12'] == 0.0 and equator['param8']['var13'] == 2.0     # the GCM's vertical wind by column
    assert c.turning_rate(cfg, 0.0) == 0.0 and c.turning_rate(c.CASES['ring_70_45e'], 0.0) == c.coriolis(90.0)
    path = c.tilted_path(360, 0.0, 0.0)
    assert np.allclose(path['lat'], 0.0) and np.allclose(path['lon'], np.arange(360) + 0.5)
    from climate.crm import ring_analysis as ra
    text = ra.evidence(0.0, dict(tilt_deg=0.0, node_deg=0.0), 'A28_dim5_moon')
    assert 'along the equator' in text and 'Coriolis' not in text and '(A28_dim5_moon, the 5% dimmer shield)' in text
    assert 'tilted 70 degrees' in ra.evidence(0.0, dict(tilt_deg=70.0, node_deg=45.0), 'A28_dim5_moon')


def test_hour_angle_along_a_tilted_ring_follows_longitude():
    from climate.crm import ring_analysis as ra
    path = c.tilted_path(360, 45.0, 180.0)
    geo = dict(x=(np.arange(360) + 0.5) * 1000.0, sun=dict(h0=-90.0 + 180.0, length=360000.0, dlon=path['dlon']))
    h = ra.hour_angle(geo, 0.0, 100.0)
    expected = (path['lon'] - 90.0 + 180.0) % 360.0 - 180.0                  # local time = longitude - 90 at t = 0
    assert np.allclose(h, expected, atol=1e-6)


@needs_planet
def test_a_box_is_three_dimensional_with_one_sun_and_its_rings_forcing():
    air = dict(sunlight_w_m2=1300.0)
    box = c.case_settings(c.CASES['box_0e'], 64, 6011.0, 111, 150000.0, air, 1.6242)
    assert box['param0']['nx'] == 64 and box['param0']['ny'] == 64 and box['param1']['dy'] == box['param1']['dx']
    assert box['param8']['var16'] == 0.0                                     # the Sun the same across the box
    assert box['param8']['var18'] == -90.0 + c.CASES['box_0e']['site']['lon_deg']   # the site's local time
    assert box['param8']['var11'] == 1.0 and box['param8']['var13'] == 1.0 and box['param8']['var12'] == 0.0
    assert box['param2']['icor'] == 0 and box['param2']['betaplane'] == 0 and box['param11']['ctrlat'] == 0.0
    assert box['param19']['do_lsnudge_v'] and box['param9']['output_vinterp'] == 1
    assert box['param9']['output_u'] == 0 and box['param9']['output_radten'] == 0
    ring = c.case_settings(c.CASES['ring_a'], 1816, 6011.0, 111, 150000.0, air, 1.6242)
    assert 'var11' not in ring['param8'] and 'output_u' not in ring['param9'] and ring['param0']['ny'] == 1


def test_day_night_table_reads_as_the_patch_expects(tmp_path):
    z = np.array([0.0, 500.0, 1000.0])
    th = np.arange(6, dtype=float).reshape(2, 3) * 1e-5
    qv = -np.arange(6, dtype=float).reshape(2, 3) * 1e-8
    c.write_day_night(tmp_path / 'terluna_lsadv.txt', dict(z_m=z, theta_k_s=th, qv_kg_kg_s=qv))
    lines = (tmp_path / 'terluna_lsadv.txt').read_text().splitlines()
    assert lines[0].split() == ['2', '3'] and np.allclose(np.array(lines[1].split(), float), z)
    for b in range(2):                                                      # per bin: a theta line, then a vapour line
        assert np.allclose(np.array(lines[2 + 2 * b].split(), float), th[b])
        assert np.allclose(np.array(lines[3 + 2 * b].split(), float), qv[b])
    assert 'terluna_nb.eq.0' in c.ADVECTION_BLOCK and 'var11.gt.0.5' in c.ADVECTION_BLOCK


def test_land_overrides_reach_the_right_columns_of_the_table():
    plain = c.landuse_rows(c.WETNESS_CLASSES)
    assert c.landuse_rows(c.WETNESS_CLASSES, None) == plain
    rough = c.landuse_rows(c.WETNESS_CLASSES, dict(roughness_cm=200.0, therin=4.97, wetness=1.0))
    fields = rough[30].split(',')                                           # index ALBD SLMO SFEM SFZ0 THERIN SCFX SFHC name
    assert [float(f) for f in fields[1:6]] == [20.0, 1.0, 0.95, 200.0, 4.97] and fields[0] == '30'
    assert float(plain[30].split(',')[4]) == 10.0 and float(plain[30].split(',')[2]) == 0.9
    assert 'roughness 200 cm' in rough[30] and 'wetness 1.00' in rough[30]


def test_the_dry_ground_land_ring_differs_from_its_control_in_the_land_alone():
    ring, control = c.CASES['ring_45n_corrected_gcm_land'], c.CASES['ring_45n_corrected']
    assert set(ring) - {'land'} == set(control) - {'land'}
    assert {k: v for k, v in ring.items() if k not in ('land', 'purpose')} == \
           {k: v for k, v in control.items() if k not in ('land', 'purpose')}
    assert ring['land'] == dict(roughness_cm=200.0, therin=4.97) and control.get('land') is None
    # a circle of latitude's one land row: the control's row (its wetness written in two places), roughness and soil
    index, base = c.PLACEHOLDER_LAND
    plain = base.replace('.50', f'{0.014:.2f}'.lstrip('0')).split(',')
    given = c.landuse_rows(((index, round(0.014, 2)),), ring['land'])[index].split(',')
    assert given[:4] == plain[:4] and given[6] == plain[6] and given[7].split("'")[0] == plain[7].split("'")[0]
    assert float(given[4]) == 200.0 and float(given[5]) == 4.97 and float(plain[4]) == 10.0 and float(plain[5]) == 4.0


def test_the_sideways_ring_differs_from_its_control_in_the_transport_alone():
    ring, control = c.CASES['ring_45n_corrected_sideways'], c.CASES['ring_45n_corrected']
    assert set(ring) - {'sideways'} == set(control) - {'sideways'} and ring['sideways'] is True
    assert {k: v for k, v in ring.items() if k not in ('sideways', 'purpose')} == \
           {k: v for k, v in control.items() if k not in ('sideways', 'purpose')}


def test_sideways_forcing_stands_ground_first_in_theta_and_vapour(monkeypatch):
    from climate.crm import ring_heat_budget as hb
    gcm = dict(sigma=[0.1, 0.5, 0.9], rows_deg=[52.6, 47.1, 41.5],
               mean_meridional_k_day=[-0.1, 0.0, -0.2], eddies_k_day=[0.0, 0.1, -0.1],
               vapour=dict(mean_meridional_g_kg_day=[0.0, 0.0, -0.1], eddies_g_kg_day=[0.0, 0.0, 0.05]))
    monkeypatch.setattr(hb, 'gcm_sideways', lambda *args: gcm)
    ref = dict(surface_pa_sea=1.2e5, profile=dict(z_m=[30000.0, 10000.0, 1000.0]))
    f = c.sideways_forcing('run', (10, 19), 45.0, 3.0, ref)
    assert f['z_m'].tolist() == [0.0, 1000.0, 10000.0, 30000.0] and f['theta_k_s'].shape == (1, 4)
    exner = (0.9 * 1.2e5 / 1.0e5) ** (287.04 / 1005.7)
    assert f['theta_k_s'][0, 0] == f['theta_k_s'][0, 1] == pytest.approx(-0.3 / 86400.0 / exner)
    assert f['qv_kg_kg_s'][0, 1] == pytest.approx(-0.05 / 8.64e7) and f['qv_kg_kg_s'][0, 3] == 0.0


@needs_planet
def test_the_sideways_ring_switches_on_the_large_scale_advection():
    air = dict(sunlight_w_m2=1198.0)
    on = c.case_settings(c.CASES['ring_45n_corrected_sideways'], 1288, 5993.1, 111, 150000.0, air, 1.6242)
    off = c.case_settings(c.CASES['ring_45n_corrected'], 1288, 5993.1, 111, 150000.0, air, 1.6242)
    assert on['param8']['var11'] == 1.0 and 'var11' not in off['param8']


@needs_planet
def test_the_land_tests_are_quarter_boxes_that_write_the_radiation_at_the_ground():
    air = dict(sunlight_w_m2=1300.0)
    small = c.case_settings(c.CASES['box_0e_small_rough'], 32, 6011.0, 111, 150000.0, air, 1.6242)
    assert small['param0']['nx'] == 32 and small['param0']['ny'] == 32 and small['param9']['output_sfcparams'] == 1
    big = c.case_settings(c.CASES['box_0e'], 64, 6011.0, 111, 150000.0, air, 1.6242)
    assert big['param9']['output_sfcparams'] == 0
    assert c.CASES['box_0e_small_gcm_land']['land'] == c.GCM_LAND and c.CASES['box_0e_small']['land'] is None


def test_box_columns_lie_along_the_heading_with_y_to_its_left(monkeypatch):
    monkeypatch.setattr(c, 'planet', lambda: dict(radius_m=1.0e6))
    lat, lon = c.box_columns(0.0, 10.0, 0.0, 3, 3, 1.0e4)               # heading north on the equator
    assert lat.shape == (3, 3) and lat[1, 1] == pytest.approx(0.0) and lon[1, 1] == pytest.approx(10.0)
    assert lat[1, 2] > 0.0 and lon[1, 2] == pytest.approx(10.0)          # x runs north
    assert lon[2, 1] < 10.0 and lat[2, 1] == pytest.approx(0.0, abs=1e-9)   # y runs west, to the left
    step = np.radians(lat[1, 2]) * 1.0e6
    assert step == pytest.approx(1.0e4, rel=1e-6)


def test_box_heights_join_at_the_edges_and_keep_the_interior():
    raw = np.full((40, 40), 1000.0)
    raw[15:25, 15:25] = 3000.0                                           # a plateau in the middle
    raw[0, :] = 500.0                                                    # a low edge
    out = c.box_heights(raw, 6000.0, passes=0, taper_m=30.0e3)
    h = out['heights']
    assert np.abs(h[0] - h[-1]).max() < 0.05 * 500.0                     # the 500 m step where the box wraps is gone
    assert np.abs(h[:, 0] - out['edge_m']).max() < 0.05 * 500.0
    assert np.allclose(h[17:23, 17:23], 3000.0)                          # beyond the taper the ground is the atlas's
    assert out['max_m'] == pytest.approx(3000.0) and out['max_slope'] > 0.0


def test_the_terrain_build_reads_every_column_s_ground():
    name, marker, edits = next(p for p in c.TERRAIN_PATCHES if p[0] == 'init_surface.F')
    text = ('      real :: sx1,sx2,sxl,stsk,stmn\n        enddo\n        close(unit=97)\n\n'
            '      ELSEIF( initsfc.ne.1 .and. initsfc.ne.2 )THEN\n')
    out = c.apply_patch(text, name, marker, edits)
    assert "file='terluna_surface2d.txt'" in out and 'logical :: terluna_two' in out
    assert out.index('terluna_surface2d.txt') < out.index('ELSEIF( initsfc.ne.1')


def test_the_highland_box_is_box_0e_over_terrain_on_ring_70_45e():
    box, base = c.CASES['box_highland'], c.CASES['box_0e']
    assert {k for k in box if box[k] != base.get(k)} == {'site', 'day_night', 'build', 'terrain', 'surface_output', 'purpose'}
    assert box['site']['ring'] == 'ring_70_45e' and box['day_night']['rings'] == ('ring_70_45e',)
    assert box['build'] == 'moon_omp_terrain' and (box['nx'], box['ny']) == (64, 64)


def test_the_terrain_build_forces_each_column_at_its_own_height():
    # every edit the terrain build makes to text the common patches wrote finds its anchor exactly once there
    assert c.ADVECTION_BLOCK.count(c.ADVECTION_OLD) == 1 and 'zh(i,j,k)' in c.ADVECTION_COLUMNS
    assert c.LSW_BLOCK.count('    IF( var13.gt.1.5 )THEN\n') == 1
    assert all(c.WSUB_BLOCK.count(old) == 1 for name, _, edits in c.TERRAIN_PATCHES if name == 'adv_routines.F'
               for old, _, _ in edits if 'terluna_w2(i' in old)
    assert 'zh(i,j,k)-var14' in c.TERRAIN_NUDGE and '{ten} = {ten}' in c.TERRAIN_NUDGE      # the hold by each column's height
    assert 'zf(i,j,k).ge.terluna_z' in c.W3_BLOCK
    air = dict(sunlight_w_m2=1300.0)
    setting = lambda n: c.case_settings(c.CASES[n], 64, 6011.0, 111, 150000.0, air, 1.6242)['param8']['var13']
    assert setting('box_highland') == 3.0 and setting('box_0e') == 1.0 and setting('ring_70_45e') == 2.0


FINE_TEMPLATE = (' &param1\n dx     =  6000.0,\n dy     =  6000.0,\n dtl    =  40.0,\n timax  = 5102784,\n'
                 ' run_time =  -999.9,\n tapfrq =   10800.0,\n rstfrq =  43200.0,\n /\n'
                 ' &param2\n irst      =  0,\n rstnum    =  1,\n ptype     =  27,\n iptra     =  1,\n npt       =  7,\n'
                 ' pdtra     =  0,\n /\n'
                 ' &param8\n var1      =   0.0,\n var2      =   0.0,\n var3      =   0.0,\n var4      =   5000.0,\n var5      =   6.8,\n var6      =   3.0,\n var7      =   12.0,\n'
                 ' var8      =   3.0,\n var9      =   0.0,\n var10     =   12000.0,\n var18     =   -90.0,\n'
                 ' var19     =   2551443.0,\n /\n'
                 ' &param14\n diagfrq        =     10800.0,\n /\n')


def coarse_run(folder):
    """A coarse box of two by two columns and three levels with one output, at day 7.5."""
    folder.mkdir(parents=True)
    np.savetxt(folder / 'input_grid_z', [0.0, 100.0, 300.0, 700.0])
    for name in ('lsnudge_0001.dat', 'terluna_wls.txt', 'terluna_lsadv.txt', 'LANDUSE.TBL'):
        (folder / name).write_text(name)
    (folder / 'terluna_surface.txt').write_text('1\n-1000000000.0 1000000000.0 1.0 30 294.163 294.163\n')
    (folder / 'namelist.template').write_text(FINE_TEMPLATE)
    (folder / 'case.json').write_text(json.dumps(dict(
        case='coarse', configuration=dict(output_s=10800.0, site=dict(lat_deg=0.0, lon_deg=0.0)),
        grid=dict(nx=2, ny=2, dx_m=6000.0, length_m=12000.0, nz=3, ztop_m=700.0),
        build=dict(label='moon_omp_elec', executable_sha256='coarse'), site=dict(ring='ring_a'))))
    names = [('th', 3), ('qv', 3), ('uinterp', 3), ('vinterp', 3), ('psfc', 0), ('t2', 0), ('q2', 0), ('tsk', 0)]
    (folder / 'cm1out_s.ctl').write_text('xdef 2 linear 0 1\nydef 2 linear 0 1\nzdef 3 levels\n1\n2\n3\n'
                                         f'vars {len(names)}\n' + ''.join(f'{n} {k} 99 {n}\n' for n, k in names) +
                                         'endvars\n')
    th = np.array([[300.0, 302.0, 304.0, 306.0], [310.0, 310.0, 310.0, 310.0], [320.0, 321.0, 322.0, 323.0]])
    fields = [th, np.full((3, 4), 0.01), np.full((3, 4), 2.0), np.full((3, 4), -1.0), np.full(4, 1.2e5),
              np.array([294.0, 295.0, 296.0, 297.0]), np.full(4, 0.012), np.array([300.0, 302.0, 304.0, 306.0])]
    n = round(7.5 * 86400.0 / 10800.0) + 1
    np.concatenate([f.ravel() for f in fields]).astype('<f4').tofile(folder / f'cm1out_t{n:06d}_s.dat')


def test_a_fine_box_starts_from_the_coarse_run_s_averaged_air_at_its_hour(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'RUNS', tmp_path / 'runs')
    monkeypatch.setattr(c, 'CM1_HOME', tmp_path / 'cm1')
    tree = tmp_path / 'tree'
    (tree / 'run').mkdir(parents=True)
    for name in ('RRTMG_LW_DATA', 'RRTMG_SW_DATA'):
        (tree / 'run' / name).write_text(name)
    monkeypatch.setattr(c, 'fetch', lambda: tree)
    build = tmp_path / 'cm1' / 'build' / 'moon_omp_elec'
    build.mkdir(parents=True)
    (build / 'cm1.exe').write_text('exe')
    (build / 'build.json').write_text(json.dumps(dict(label='moon_omp_elec', executable_sha256='fine')))
    coarse_run(tmp_path / 'runs' / 'coarse')
    monkeypatch.setitem(c.CASES, 'fine', dict(c.CASES['box_0e_elec_fine'],
                                              fine_from=dict(case='coarse', day=None, refine=3)))
    with pytest.raises(RuntimeError, match='choose the day'):
        c.setup('fine')
    monkeypatch.setitem(c.CASES, 'fine', dict(c.CASES['box_0e_elec_fine'],
                                              fine_from=dict(case='coarse', day=7.5, refine=3)))
    case = c.setup('fine')
    lines = (case / 'input_sounding').read_text().splitlines()
    assert [float(v) for v in lines[0].split()] == pytest.approx([1200.0, 295.5 * (1.0e5 / 1.2e5) ** (287.04 / 1005.7),
                                                                 12.0], abs=1e-3)
    rows = np.array([[float(v) for v in line.split()] for line in lines[1:]])
    assert rows[:, 0].tolist() == [50.0, 200.0, 500.0, 700.0]                    # scalar levels, then the top
    assert rows[:3, 1].tolist() == pytest.approx([303.0, 310.0, 321.5])         # the means over the columns
    assert rows[3, 1] == pytest.approx(321.5 + 11.5 / 300.0 * 200.0, abs=1e-3)   # carried on the top two levels' slope
    assert rows[:, 2] == pytest.approx(10.0) and rows[:, 3] == pytest.approx(2.0) and rows[:, 4] == pytest.approx(-1.0)
    assert (case / 'terluna_surface.txt').read_text().split()[5:] == ['303.000', '294.163']   # mean skin, deep ground
    text = (case / 'namelist.template').read_text()
    hour = (-90.0 + 360.0 * 7.5 * 86400.0 / 2551443.0 + 180.0) % 360.0 - 180.0
    assert c.namelist_value(text, 'param1', 'dx') == '2000.0' and c.namelist_value(text, 'param1', 'dtl') == '13.333'
    assert float(c.namelist_value(text, 'param8', 'var18')) == pytest.approx(hour, abs=1e-4)
    assert c.namelist_value(text, 'param1', 'tapfrq') == '900.0' and c.namelist_value(text, 'param8', 'var5') == '6.8'
    assert c.namelist_value(text, 'param1', 'timax') == str(round(4.0 * 86400.0))
    assert (case / 'lsnudge_0001.dat').read_text() == 'lsnudge_0001.dat' and (case / 'cm1.exe').is_symlink()
    record = json.loads((case / 'case.json').read_text())
    assert record['grid']['dx_m'] == pytest.approx(2000.0) and record['grid']['length_m'] == pytest.approx(4000.0)
    assert record['initial']['snapshot'] == 61 and record['initial']['day'] == 7.5
    assert record['configuration']['start_hour_angle_deg'] == pytest.approx(hour, abs=1e-4)
    assert record['build']['executable_sha256'] == 'fine' and record['site'] == dict(ring='ring_a')


def test_a_rerun_from_a_restart_takes_the_run_s_inputs_and_restart_and_its_own_electricity(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'RUNS', tmp_path / 'runs')
    monkeypatch.setattr(c, 'CM1_HOME', tmp_path / 'cm1')
    tree = tmp_path / 'tree'
    (tree / 'run').mkdir(parents=True)
    for name in ('RRTMG_LW_DATA', 'RRTMG_SW_DATA'):
        (tree / 'run' / name).write_text(name)
    monkeypatch.setattr(c, 'fetch', lambda: tree)
    build = tmp_path / 'cm1' / 'build' / 'moon_omp_elec'
    build.mkdir(parents=True)
    (build / 'cm1.exe').write_text('exe')
    (build / 'build.json').write_text(json.dumps(dict(label='moon_omp_elec', executable_sha256='rerun')))
    source = tmp_path / 'runs' / 'coarse'
    coarse_run(source)
    for part in 'isuvwx':
        (source / f'cm1rst_t000015_{part}.dat').write_text(part)
    (source / 'progress.json').write_text(json.dumps(dict(segments=[
        dict(segment=7, model_s=604800.0, executable='first'), dict(segment=8, model_s=691200.0, executable='second')])))
    monkeypatch.setitem(c.CASES, 'again', dict(c.CASES['box_0e_elec_uncapped'],
                                               restart_from=dict(case='coarse', day=7.0)))
    with pytest.raises(RuntimeError, match='no restart at day 7.0'):
        c.setup('again')
    monkeypatch.setitem(c.CASES, 'again', dict(c.CASES['box_0e_elec_uncapped_corona'],
                                               restart_from=dict(case='coarse', day=7.5)))
    with pytest.raises(RuntimeError, match='vertical field'):
        c.setup('again')
    (source / 'terluna_ez_648000.bin').write_text('ez')
    case = c.setup('again')
    assert all((case / f'cm1rst_t000015_{part}.dat').read_text() == part for part in 'isuvwx')
    assert (case / 'terluna_ez_648000.bin').read_text() == 'ez' and (case / 'cm1.exe').is_symlink()
    assert (case / 'lsnudge_0001.dat').read_text() == 'lsnudge_0001.dat' and (case / 'cm1out_s.ctl').exists()
    text = (case / 'namelist.template').read_text()
    assert c.namelist_value(text, 'param8', 'var8') == '4.0' and c.namelist_value(text, 'param8', 'var2') == '3000.0'
    assert c.namelist_value(text, 'param1', 'dx') == '6000.0'                  # the run's own grid
    progress = json.loads((case / 'progress.json').read_text())
    assert progress['segments'] == [] and progress['start_s'] == 648000.0
    record = json.loads((case / 'case.json').read_text())
    assert record['restart_from'] == dict(case='coarse', day=7.5, restart=15, model_s=648000.0, made_by='second')
    assert record['electricity']['corona_v_m'] == 3000.0 and record['build']['executable_sha256'] == 'rerun'
    assert c.status('again').startswith('again: 7.50 of 14.50 days')
    with pytest.raises(RuntimeError, match='has been set up'):
        c.setup('again')


def test_vapour_above_the_gcms_top_keeps_that_levels_share_of_ice_saturation():
    from atmosphere.radiative_convective.thermodynamics import saturation_pressure
    zg = np.array([1e3, 20e3, 60e3, 100e3, 130e3])
    th = np.array([300.0, 310.0, 330.0, 370.0, 420.0])
    qv = np.array([1.2e-2, 3e-3, 2e-4, 4e-6, 3e-7])
    slope = (th[-1] - th[-2]) / (zg[-1] - zg[-2])
    theta_at = lambda z: np.where(z <= zg[-1], np.interp(z, zg, th), th[-1] + slope * (z - zg[-1]))
    qv_below = lambda z: np.exp(np.interp(z, zg, np.log(qv)))
    qv_at = c.held_share_above(zg, theta_at, qv_below, 121400.0, 1.6242, 150e3)
    below = np.linspace(0.0, 130e3, 50)
    assert np.array_equal(qv_at(below), qv_below(below))
    # Above the top the vapour falls, and its share of ice saturation stays the top level's.
    z = np.linspace(130e3, 150e3, 9)
    q = qv_at(z)
    assert np.all(np.diff(q) < 0)
    cp, rd, g = 1005.7, 287.04, 1.6242
    zz = np.linspace(0.0, 150e3, 30001)
    inv = 1.0 / (theta_at(zz) * (1.0 + 0.61 * qv_at(zz)))
    exner = (121400.0 / 1e5) ** (rd / cp) - g / cp * np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(zz))])
    p, t = 1e5 * np.interp(z, zz, exner) ** (cp / rd), theta_at(z) * np.interp(z, zz, exner)
    share = q * p / (0.622 + q) / saturation_pressure(t)
    assert np.allclose(share, share[0], rtol=0.01)


def test_the_sea_build_makes_the_two_dimensional_water_columns_sea():
    name, marker, edits = c.SEA_COLUMNS_PATCH
    assert name == 'init_surface.F' and c.BUILD_PATCHES['moon_omp_terrain_sea'] == c.TERRAIN_PATCHES + [c.SEA_COLUMNS_PATCH]
    old, new, count = edits[0]
    assert count == 1 and old in c.SURFACE_2D_BLOCK and marker in new and marker not in old
    assert 'xland(i,j) = merge(2.0,1.0,' in new and '.eq.16)' in new          # land use 16, water, is sea
    assert c.BUILDS['moon_omp_terrain_sea'] == c.BUILDS['moon_omp_terrain']


def test_a_box_keeps_its_waters_at_their_level_and_smooths_its_land():
    raw = np.zeros((20, 20))
    raw[8:12, 8:12] = 2000.0                                                  # an island in a sea
    wet = raw == 0.0
    out = c.box_heights(raw, 6000.0, passes=2, taper_m=30000.0, wet=wet)
    assert np.all(out['heights'][wet] == 0.0)                                  # the sea stays at sea level
    assert out['heights'][~wet].max() < 2000.0 and out['heights'][~wet].min() > 0.0
    dry = c.box_heights(raw, 6000.0, passes=2, taper_m=30000.0)                # without the mask the island spreads
    assert (dry['heights'][wet] > 0.0).any()


def test_a_box_on_a_parallel_ring_takes_the_rows_reference_and_the_columns_round_it(monkeypatch):
    profile = dict(sea_surface_k=296.0, land_ground_k=295.0, land_wetness=0.02, lapse_k_m=0.0016, rows_deg=[45.0],
                   profile=dict(z_m=[0.0, 1000.0]))
    monkeypatch.setattr(c, 'gcm_equator_profile', lambda run, years, band, latitude: dict(profile, latitude_deg=latitude))
    monkeypatch.setattr(c, 'gcm_vertical_wind', lambda lat, band, run, years: dict(z_m=[0.0, 5.0e4], w_m_s=[0.0, -0.005]))
    cfg = c.CASES['box_jura']
    zw = [0.0, 1000.0, 2.5e4, 2.0e5]
    site = c.parallel_box_site(cfg, zw, 1288, c.CASES['ring_45n_corrected'], *c.CORRECTED_GCM)
    assert site['heading_deg'] == 90.0 and site['sea_k'] == 296.0 and site['wetness'] == 0.02
    assert site['w_m_s'][0] == 0.0 and site['w_m_s'][-1] == 0.0 and site['w_m_s'][2] == pytest.approx(-0.0025)
    # the ring runs east from 0 E: the columns within 100 km of 323.875 E
    length = 2.0 * np.pi * c.planet()['radius_m'] * np.cos(np.radians(45.0))
    centres = (np.array(site['columns']) + 0.5) * length / 1288 / length * 360.0
    assert np.all(np.abs(centres - 323.875) <= 100.0e3 / length * 360.0 + 1e-9) and len(site['columns']) >= 30
    with pytest.raises(ValueError):
        c.parallel_box_site(dict(cfg, site=dict(cfg['site'], lat_deg=20.0)), zw, 1288, c.CASES['ring_45n_corrected'])


def test_the_jura_box_is_the_island_and_its_sea_on_the_corrected_design_case():
    box = c.CASES['box_jura']
    assert box['kind'] == 'box' and (box['nx'], box['ny']) == (100, 90) and box['gcm'] == c.CORRECTED_GCM
    assert box['terrain']['sea'] and box['build'] == 'moon_omp_terrain_sea_twilight' and box['twilight'] == 'full'
    assert box['site']['ring'] == 'ring_45n_corrected' and box['day_night']['rings'] == ('ring_45n_corrected_twilight',)
    assert c.BUILD_PATCHES[box['build']] == c.TERRAIN_PATCHES + [c.SEA_COLUMNS_PATCH, c.TWILIGHT_PATCH]
    # the hook places each shell by its column's own face heights, which over terrain stand on the column's ground
    assert any('zf(i,j,k+1)' in new for _, new, _ in c.TWILIGHT_PATCH[2])


NIGHT_TEMPLATE = (' &param0\n nx           =     2,\n ny           =     2,\n nz           =      3,\n /\n'
                  ' &param1\n dx     =  6000.0,\n dy     =  6000.0,\n dtl    =  40.0,\n timax  = 5102784,\n'
                  ' run_time =  -999.9,\n tapfrq =   10800.0,\n rstfrq =  43200.0,\n statfrq =   3600.0,\n /\n'
                  ' &param2\n irst      =  0,\n rstnum    =  1,\n /\n'
                  ' &param8\n var11     =   1.0,\n var14     =   93132.6,\n var18     =   -90.0,\n'
                  ' var19     =   2551443.0,\n var20     =   1.0,\n /\n'
                  ' &param9\n output_radten    = 1,\n output_thbudget  = 1,\n output_qvbudget  = 1,\n /\n'
                  ' &param14\n diagfrq        =     10800.0,\n /\n')


def night_source(folder):
    """A sea box of two by two columns and three levels with one output at day 7.5 and three wind windows, the
    second of which holds that time."""
    coarse_run(folder)
    (folder / 'namelist.template').write_text(NIGHT_TEMPLATE)
    (folder / 'terluna_surface.txt').write_text('1\n-1000000000.0 1000000000.0 2.0 16 302.036 302.036\n')
    for name in ('terluna_twilight.txt', 'terluna_rain_law.txt'):
        (folder / name).write_text(name)
    (folder / 'lsnudge_0001.dat').unlink()
    edges = [0.0, 600000.0, 700000.0, 1.0e30]
    for k in range(3):
        (folder / f'lsnudge_{k + 1:04d}.dat').write_text(
            f' *  Header\n  {edges[k]:22.1f} {edges[k + 1]:22.6e}\n *  Profile\n *  columns\n'
            f'        0.0000     284.0000    13.000000     {k:.4f}     1.0000\n')
    record = json.loads((folder / 'case.json').read_text())
    record['configuration']['site'] = dict(lat_deg=-5.0, lon_deg=-45.0, surface='sea')
    record.update(restart_from=dict(case='earlier', day=0.0), nudging_wind=dict(windows=3))
    (folder / 'case.json').write_text(json.dumps(record))


def test_a_held_night_starts_from_the_source_s_mean_air_and_holds_that_hour_s_forcing(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'RUNS', tmp_path / 'runs')
    monkeypatch.setattr(c, 'CM1_HOME', tmp_path / 'cm1')
    tree = tmp_path / 'tree'
    (tree / 'run').mkdir(parents=True)
    for name in ('RRTMG_LW_DATA', 'RRTMG_SW_DATA'):
        (tree / 'run' / name).write_text(name)
    monkeypatch.setattr(c, 'fetch', lambda: tree)
    build = tmp_path / 'cm1' / 'build' / 'moon_omp_rain_law_pbl'
    build.mkdir(parents=True)
    (build / 'cm1.exe').write_text('exe')
    (build / 'build.json').write_text(json.dumps(dict(label='moon_omp_rain_law_pbl', executable_sha256='night')))
    source = tmp_path / 'runs' / 'coarse'
    night_source(source)
    monkeypatch.setitem(c.CASES, 'night', dict(c.CASES['box_near_side_sea_night_2km'],
                                               night_from=dict(case='coarse', day=7.0, refine=3)))
    with pytest.raises(RuntimeError, match='no output at day 7.0'):
        c.setup('night')
    monkeypatch.setitem(c.CASES, 'night', dict(c.CASES['box_near_side_sea_night_2km'],
                                               night_from=dict(case='coarse', day=7.5, refine=3)))
    case = c.setup('night')
    rows = np.array([[float(v) for v in line.split()] for line in (case / 'input_sounding').read_text().splitlines()[1:]])
    assert rows[:3, 1].tolist() == pytest.approx([303.0, 310.0, 321.5])          # the source's means over its columns
    # the wind window holding day 7.5 (648,000 s), its profile kept and held from the start to the end
    nudge = (case / 'lsnudge_0001.dat').read_text().splitlines()
    assert [float(v) for v in nudge[1].split()] == [0.0, 1.0e30] and nudge[4].split()[3] == '1.0000'
    assert nudge[:1] + nudge[2:] == (source / 'lsnudge_0002.dat').read_text().splitlines()[:1] + \
        (source / 'lsnudge_0002.dat').read_text().splitlines()[2:]
    assert not (case / 'lsnudge_0002.dat').exists()
    for name in ('terluna_lsadv.txt', 'terluna_wls.txt', 'terluna_surface.txt', 'terluna_twilight.txt',
                 'terluna_rain_law.txt', 'LANDUSE.TBL', 'input_grid_z'):
        assert (case / name).read_bytes() == (source / name).read_bytes(), name   # the source's own, as written
    text = (case / 'namelist.template').read_text()
    hour = (-90.0 + 360.0 * 7.5 * 86400.0 / 2551443.0 + 180.0) % 360.0 - 180.0
    value = lambda section, key: c.namelist_value(text, section, key)
    assert float(value('param8', 'var18')) == pytest.approx(hour, abs=1e-4)
    assert float(value('param8', 'var19')) == c.STOPPED_CLOCK_S                  # the Sun and the table held
    assert value('param8', 'var14') == '93132.6' and value('param8', 'var11') == '1.0'
    assert value('param0', 'nx') == '64' and value('param0', 'ny') == '64'
    assert value('param1', 'dx') == '2000.0' and value('param1', 'dtl') == '13.333'
    assert value('param1', 'timax') == str(round(2.0 * 86400.0)) and value('param1', 'tapfrq') == '21600.0'
    assert value('param1', 'rstfrq') == '43200.0' and value('param1', 'statfrq') == '900.0'
    assert value('param14', 'diagfrq') == '1800.0' and value('param2', 'irst') == '0'
    assert [value('param9', k) for k in ('output_thbudget', 'output_qvbudget')] == ['0', '0']   # the budgets' 3-D fields
    assert value('param9', 'output_radten') == '1'                     # as the source's: CM1 writes it with radiation
    record = json.loads((case / 'case.json').read_text())
    assert record['restart_from'] is None and record['nudging_wind']['windows'] == 1
    assert record['grid']['nx'] == 64 and record['grid']['ny'] == 64 and record['grid']['dx_m'] == pytest.approx(2000.0)
    assert record['initial']['snapshot'] == 61 and record['held']['hour_angle_deg'] == pytest.approx(hour, abs=1e-4)
    assert record['held']['wind_window']['source_file'] == 'lsnudge_0002.dat'
    assert record['held']['sun_elevation_deg'] == pytest.approx(
        np.degrees(np.arcsin(np.cos(np.radians(5.0)) * np.cos(np.radians(hour)))), abs=0.01)
    assert record['configuration']['start_hour_angle_deg'] == pytest.approx(hour + 45.0, abs=1e-4)
    assert record['build']['executable_sha256'] == 'night'
    assert record['runner'] == hashlib.sha256(open(c.__file__, 'rb').read()).hexdigest()[:16]


def test_the_held_nights_differ_in_their_grid_alone():
    names = [f'box_near_side_sea_night_{tag}' for tag in ('6km', '3km', '2km')]
    boxes = [c.CASES[n] for n in names]
    grid = ('nx', 'ny', 'night_from', 'purpose')
    for box in boxes[1:]:
        assert {k: v for k, v in box.items() if k not in grid} == {k: v for k, v in boxes[0].items() if k not in grid}
    assert [(b['nx'], b['night_from']['refine']) for b in boxes] == [(32, 1), (64, 2), (64, 3)]
    assert all(b['night_from']['case'] == 'box_near_side_sea_day3' and b['night_from']['day'] == 84.5 for b in boxes)
    sea = c.CASES['box_near_side_sea']
    shared = ('build', 'nudge', 'large_scale_w', 'day_night', 'site', 'twilight', 'twilight_product', 'rain_law',
              'rain_limits', 'z_grid', 'droplets_cm3', 'gcm', 'land', 'surface_output', 'wind_by_local_time')
    assert all(boxes[0][k] == sea[k] for k in shared)                            # the sea box's build and forcing
    assert boxes[0]['days'] == 2.0 and boxes[0]['output_s'] == 6 * 3600.0 and boxes[0]['restart_s'] == 12 * 3600.0
    assert c.case_length_s(boxes[0], 1.0) == 2.0 * 86400.0


def test_the_earth_length_cases_differ_from_their_sources_in_the_length_and_what_they_write_alone():
    day3, linf = c.CASES['box_near_side_sea_day3'], c.CASES['box_near_side_sea_day3_linf75']
    own = ('output_s', 'keep_restarts', 'restart_namelist', 'purpose')
    assert {k: v for k, v in linf.items() if k not in own} == {k: v for k, v in day3.items() if k not in own}
    assert linf['restart_from'] == dict(case='box_near_side_sea', day=59.0) and linf['days'] == 88.5
    assert linf['restart_namelist'] == dict(param3=dict(l_inf=75.0), param9=dict(output_thbudget=0, output_qvbudget=0))
    assert linf['output_s'] == 6 * 3600.0 and linf['restart_s'] == day3['restart_s'] and linf['keep_restarts'] == 2
    night, night75 = c.CASES['box_near_side_sea_night_6km'], c.CASES['box_near_side_sea_night_6km_linf75']
    assert {k: v for k, v in night75.items() if k not in ('namelist', 'purpose')} == \
        {k: v for k, v in night.items() if k not in ('namelist', 'purpose')}
    assert night75['namelist'] == dict(c.NIGHT_OUTPUT, param3=dict(l_inf=75.0))
    # no case set up before them takes either new key, so the runner sets them up and runs them as before (the
    # wind's own time scale of 28 October, and the gust beside it, take restart_namelist after them)
    new = {'box_near_side_sea_day3_linf75', 'box_near_side_sea_night_6km_linf75', 'box_near_side_sea_l15_wind_day3',
           'box_near_side_sea_l15_wind_day4', 'box_near_side_sea_l15_wind_gust_day3', 'box_near_side_sea_l15_wind_gust_day4'}
    assert not [n for n, v in c.CASES.items() if n not in new and ('restart_namelist' in v or 'keep_restarts' in v)]


def test_pruning_keeps_the_starting_restart_the_newest_and_any_incomplete_set(tmp_path):
    for n in (118, 119, 120, 121, 122):
        for part in 'isuvwx':
            (tmp_path / f'cm1rst_t{n:06d}_{part}.dat').write_bytes(b'')
    for part in 'isuv':                                                   # being written
        (tmp_path / f'cm1rst_t000123_{part}.dat').write_bytes(b'')
    assert c.prune_restarts(tmp_path, 2, start=118) == [119, 120]
    left = sorted({p.name[8:14] for p in tmp_path.glob('cm1rst_t*_s.dat')})
    assert left == ['000118', '000121', '000122', '000123'] and c.latest_restart(tmp_path) == 122
    assert len(list(tmp_path.glob('cm1rst_t000121_*.dat'))) == 6 and c.prune_restarts(tmp_path, 2, start=118) == []


LENGTH_TEMPLATE = FINE_TEMPLATE + (' &param3\n fcor    = -1.155e-06,\n l_inf   =     453.0,\n /\n'
                                   ' &param9\n output_radten    = 1,\n output_thbudget  = 1,\n output_qvbudget  = 1,\n /\n')


def test_a_run_on_with_earth_s_length_changes_the_source_s_namelist_there_alone(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'RUNS', tmp_path / 'runs')
    monkeypatch.setattr(c, 'CM1_HOME', tmp_path / 'cm1')
    tree = tmp_path / 'tree'
    (tree / 'run').mkdir(parents=True)
    for name in ('RRTMG_LW_DATA', 'RRTMG_SW_DATA'):
        (tree / 'run' / name).write_text(name)
    monkeypatch.setattr(c, 'fetch', lambda: tree)
    monkeypatch.setattr(c, 'apply_rain_limits', lambda case, cfg, g: None)   # the source's limits, written elsewhere
    build = tmp_path / 'cm1' / 'build' / 'moon_omp_rain_law_pbl'
    build.mkdir(parents=True)
    (build / 'cm1.exe').write_text('exe')
    (build / 'build.json').write_text(json.dumps(dict(label='moon_omp_rain_law_pbl', executable_sha256='linf')))
    source = tmp_path / 'runs' / 'coarse'
    coarse_run(source)
    (source / 'namelist.template').write_text(LENGTH_TEMPLATE)
    for part in 'isuvwx':
        (source / f'cm1rst_t000015_{part}.dat').write_text(part)
    (source / 'progress.json').write_text(json.dumps(dict(segments=[dict(segment=8, model_s=691200.0,
                                                                          executable='source')])))
    monkeypatch.setitem(c.CASES, 'linf', dict(c.CASES['box_near_side_sea_day3_linf75'],
                                              restart_from=dict(case='coarse', day=7.5), days=14.5))
    case = c.setup('linf')
    text = (case / 'namelist.template').read_text()
    value = lambda section, key: c.namelist_value(text, section, key)
    assert value('param3', 'l_inf') == '75.0' and value('param3', 'fcor') == '-1.155e-06'
    assert [value('param9', k) for k in ('output_thbudget', 'output_qvbudget', 'output_radten')] == ['0', '0', '1']
    assert value('param1', 'tapfrq') == '21600.0' and value('param14', 'diagfrq') == '10800.0'   # 3-hourly diagnostics
    assert value('param1', 'rstfrq') == '43200.0'                       # restarts numbered as the source's
    changed = [(a, b) for a, b in zip(LENGTH_TEMPLATE.splitlines(), text.splitlines()) if a != b]
    assert len(changed) == 4 and len(text.splitlines()) == len(LENGTH_TEMPLATE.splitlines())
    assert all((case / f'cm1rst_t000015_{part}.dat').read_text() == part for part in 'isuvwx')
    record = json.loads((case / 'case.json').read_text())
    assert record['configuration']['restart_namelist']['param3'] == dict(l_inf=75.0)
    assert record['restart_from']['restart'] == 15 and record['build']['executable_sha256'] == 'linf'


def test_the_held_night_with_earth_s_length_takes_it_over_the_source_s(tmp_path, monkeypatch):
    monkeypatch.setattr(c, 'RUNS', tmp_path / 'runs')
    monkeypatch.setattr(c, 'CM1_HOME', tmp_path / 'cm1')
    tree = tmp_path / 'tree'
    (tree / 'run').mkdir(parents=True)
    for name in ('RRTMG_LW_DATA', 'RRTMG_SW_DATA'):
        (tree / 'run' / name).write_text(name)
    monkeypatch.setattr(c, 'fetch', lambda: tree)
    build = tmp_path / 'cm1' / 'build' / 'moon_omp_rain_law_pbl'
    build.mkdir(parents=True)
    (build / 'cm1.exe').write_text('exe')
    (build / 'build.json').write_text(json.dumps(dict(label='moon_omp_rain_law_pbl', executable_sha256='night')))
    source = tmp_path / 'runs' / 'coarse'
    night_source(source)
    (source / 'namelist.template').write_text(NIGHT_TEMPLATE + ' &param3\n l_inf   =     453.0,\n /\n')
    for name, base in (('night75', 'box_near_side_sea_night_6km_linf75'), ('night', 'box_near_side_sea_night_6km')):
        monkeypatch.setitem(c.CASES, name, dict(c.CASES[base], night_from=dict(case='coarse', day=7.5, refine=1)))
    long_, short = (c.namelist_value((c.setup(n) / 'namelist.template').read_text(), 'param3', 'l_inf')
                    for n in ('night', 'night75'))
    assert (long_, short) == ('453.0', '75.0')
    a, b = ((c.RUNS / n / 'namelist.template').read_text().splitlines() for n in ('night', 'night75'))
    assert [(x, y) for x, y in zip(a, b) if x != y] == [(' l_inf   =     453.0,', ' l_inf   =     75.0,')]
