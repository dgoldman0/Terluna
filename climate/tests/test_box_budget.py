"""The box budget's reading of CM1's domain diagnostics and of the prescribed day-night heating (no CM1 run)."""
import numpy as np
import pytest

from climate.crm import box_budget as bb
from climate.crm import cm1_output as co


def test_the_day_night_heating_follows_the_patch_round_the_day():
    z = np.array([0.0, 1000.0, 2000.0])
    th = np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0], [4.0, 4.0, 4.0]])   # four bins of 90 degrees
    zh = np.array([500.0, 1500.0, 2500.0])
    # var18 = -180 puts the hour angle at -180 degrees at time 0: halfway between the last bin's centre and the first's
    assert np.allclose(bb.day_night_heating(z, th, zh, 0.0, -180.0, 360.0), [2.5, 2.5, 0.0])
    # 45 degrees later the hour angle stands at the first bin's centre
    assert np.allclose(bb.day_night_heating(z, th, zh, 45.0, -180.0, 360.0), [1.0, 1.0, 0.0])
    # a quarter of the way on from there, between the first and second bins
    assert np.allclose(bb.day_night_heating(z, th, zh, 67.5, -180.0, 360.0), [1.25, 1.25, 0.0])


def test_layer_sums_weigh_each_level_by_its_density_and_depth():
    zw = np.array([0.0, 1000.0, 3000.0, 14000.0, 30000.0])
    rho = np.array([1.4, 1.3, 1.0, 0.5])
    rate = np.array([1.0, 2.0, -1.0, 0.5])                               # K/day
    out = bb.layer_sums(rate, rho, zw)
    assert out['0-1.77 km'] == pytest.approx(1.0 * 1.4 * co.CP * 1000.0 / bb.DAY)
    assert out['1.77-6.53 km'] == pytest.approx(2.0 * 1.3 * co.CP * 2000.0 / bb.DAY)
    assert out['column'] == pytest.approx(sum(r * p * co.CP * d for r, p, d in zip(rate, rho, np.diff(zw))) / bb.DAY)


def test_layer_means_weigh_each_level_by_its_mass():
    zw = np.array([0.0, 1000.0, 3000.0, 14000.0, 30000.0])
    rho = np.array([1.4, 1.3, 1.0, 0.5])
    t = np.array([300.0, 290.0, 280.0, 250.0])
    out = bb.layer_means(t, rho, zw)
    assert out['0-1.77 km'] == pytest.approx(300.0) and out['1.77-6.53 km'] == pytest.approx(290.0)
    assert out['6.53-13.74 km'] == pytest.approx(280.0) and out['13.74-23.06 km'] == pytest.approx(250.0)
    two = bb.layer_means(t, rho, zw, layers=((0.0, 3.0),))               # the lowest two levels together, by mass
    assert two['0-3 km'] == pytest.approx((300.0 * 1400.0 + 290.0 * 2600.0) / 4000.0)


def test_a_diagnostic_output_reads_by_its_descriptor(tmp_path):
    (tmp_path / 'cm1out_diag_s.ctl').write_text('dset ^cm1out_diag_00%y4_s.dat\nxdef 1 linear 0 1\nydef 1 linear 0 1\n'
                                                'zdef 3 levels\n0.05\n0.15\n0.25\nvars 2\nmtime 0 99 time\n'
                                                'th 3 99 theta\nendvars\n')
    path = tmp_path / 'cm1out_diag_000002_s.dat'
    np.array([10800.0, 300.0, 301.0, 302.0], dtype='<f4').tofile(path)
    d = bb.read_diag(tmp_path, path)
    assert d['mtime'] == 10800.0 and np.allclose(d['th'], [300.0, 301.0, 302.0])
    np.array([10800.0, 300.0, 301.0], dtype='<f4').tofile(path)
    with pytest.raises(ValueError):
        bb.read_diag(tmp_path, path)


def test_the_water_budget_closes_and_reads_the_boundary_layer_against_the_evaporation(tmp_path, monkeypatch):
    """Two levels of 1000 kg/m2 each a day apart: the boundary layer brings 3 mm/day at the ground, the imposed
    ascent and the nudging take 0.5 each, 1 mm/day condenses and rains out, and the column keeps the other 1."""
    import json
    monkeypatch.setattr(bb.ra, 'RUNS', tmp_path)
    case = tmp_path / 'box'
    case.mkdir()
    (case / 'case.json').write_text(json.dumps(dict(grid=dict(nx=1, dx_m=1000.0, length_m=1000.0),
                                                    configuration=dict(start_hour_angle_deg=0.0))))
    (case / 'input_grid_z').write_text('0\n1000\n3000\n')
    (case / 'terluna_surface.txt').write_text('1\n0 1000 2 16 300 300\n')
    (case / 'namelist.input').write_text(' var11 = 0.0,\n var18 = 0.0,\n var19 = 1.0,\n')
    names = ('qv', 'qc', 'rho', 'qvb_pbl', 'qvb_mp', 'qvb_lsw', 'qvb_nudge', 'qvb_hadv', 'qvb_vadv')
    (case / 'cm1out_diag_s.ctl').write_text('xdef 1 linear 0 1\nydef 1 linear 0 1\nzdef 2 levels\n500\n2000\n'
                                            f'vars {3 + len(names)}\nmtime 0 99 t\nqfx 0 99 e\nprate 0 99 p\n'
                                            + ''.join(f'{n} 2 99 x\n' for n in names) + 'endvars\n')
    per = 1.0 / (bb.DAY * 1000.0)                                         # 1 mm/day in one level's 1000 kg/m2
    for n, (day, qv) in enumerate(((29.5, (0.010, 0.005)), (30.5, (0.0105, 0.0055))), 1):
        fields = dict(qv=qv, qc=(0.0, 0.0), rho=(1.0, 0.5), qvb_pbl=(3 * per, 0.0), qvb_mp=(0.0, -per),
                      qvb_lsw=(-0.5 * per, 0.0), qvb_nudge=(0.0, -0.5 * per), qvb_hadv=(0.2 * per, -0.2 * per),
                      qvb_vadv=(0.0, 0.0))
        values = [day * bb.DAY, 3.0 / bb.DAY, 1.0 / bb.DAY] + [v for name in names for v in fields[name]]
        np.array(values, dtype='<f4').tofile(case / f'cm1out_diag_{n:06d}_s.dat')
    w = bb.water('box', 29.5)
    assert w['outputs'] == 2 and w['span_days'] == pytest.approx([29.5, 30.5])
    c = w['column_mm_day']
    assert c['boundary_layer'] == pytest.approx(3.0, rel=1e-5) and c['microphysics'] == pytest.approx(-1.0, rel=1e-5)
    assert c['imposed_ascent'] == pytest.approx(-0.5, rel=1e-5) and c['nudging'] == pytest.approx(-0.5, rel=1e-5)
    assert c['own_transport'] == pytest.approx(0.0, abs=1e-6)
    assert w['bands_mm_day']['nudging']['0-16 km'] == pytest.approx(-0.5, rel=1e-5)
    assert w['evaporation_mm_day'] == pytest.approx(3.0, rel=1e-5) and w['rain_mm_day'] == pytest.approx(1.0, rel=1e-5)
    assert w['storage_mm_day'] == pytest.approx(1.0, rel=1e-4) and w['rest_mm_day'] == pytest.approx(0.0, abs=1e-4)
    assert w['boundary_layer_over_evaporation'] == pytest.approx(1.0, rel=1e-5)
    assert bb.water('box', 29.5, to_day=30.5)['storage_mm_day'] == pytest.approx(1.0, rel=1e-4)   # the span's end
    with pytest.raises(RuntimeError, match='fewer than two'):
        bb.water('box', 29.5, to_day=30.0)


def test_the_inflow_of_the_reference_s_air_is_reckoned_from_its_written_rates(tmp_path, monkeypatch):
    """The water test's two levels on the inflow build: the lower draws the reference's air in at 1e-5 /s, the upper
    lets the box's own out; the box stands 1 then 0.5 g/kg below the reference there, so the inflow brings
    1e-5 * 0.75e-3 kg/kg * 1000 kg/m2 a second, 0.648 mm/day, and the rest takes it back out of the closure."""
    import json
    monkeypatch.setattr(bb.ra, 'RUNS', tmp_path)
    case = tmp_path / 'box'
    case.mkdir()
    (case / 'case.json').write_text(json.dumps(dict(grid=dict(nx=1, dx_m=1000.0, length_m=1000.0),
                                                    configuration=dict(start_hour_angle_deg=0.0))))
    (case / 'input_grid_z').write_text('0\n1000\n3000\n')
    (case / 'terluna_surface.txt').write_text('1\n0 1000 2 16 300 300\n')
    (case / 'namelist.input').write_text(' var7 = 1.0,\n var11 = 0.0,\n var18 = 0.0,\n var19 = 1.0,\n')
    (case / bb.INFLOW_RATES).write_text(' 2\n   500.0   1.0E-05\n   2000.0   0.0\n')
    (case / 'lsnudge_0001.dat').write_text(' *  Header:   lsnudge_time1 (s)   lsnudge_time2 (s)\n   0.0   7.0e4\n'
                                           ' *  Profile.\n *    z (m)   theta (K)    qv (g/kg)  u  v\n'
                                           '  0.0  300.0  11.0  1.0  0.0\n  500.0  301.0  11.0  1.0  0.0\n'
                                           '  2000.0  305.0  5.5  1.0  0.0\n  3000.0  308.0  5.5  1.0  0.0\n')
    names = ('qv', 'qc', 'rho', 'qvb_pbl', 'qvb_mp', 'qvb_lsw', 'qvb_nudge', 'qvb_hadv', 'qvb_vadv')
    (case / 'cm1out_diag_s.ctl').write_text('xdef 1 linear 0 1\nydef 1 linear 0 1\nzdef 2 levels\n500\n2000\n'
                                            f'vars {3 + len(names)}\nmtime 0 99 t\nqfx 0 99 e\nprate 0 99 p\n'
                                            + ''.join(f'{n} 2 99 x\n' for n in names) + 'endvars\n')
    per = 1.0 / (bb.DAY * 1000.0)
    for n, (day, qv) in enumerate(((29.5, (0.010, 0.005)), (30.5, (0.0105, 0.0055))), 1):
        fields = dict(qv=qv, qc=(0.0, 0.0), rho=(1.0, 0.5), qvb_pbl=(3 * per, 0.0), qvb_mp=(0.0, -per),
                      qvb_lsw=(-0.5 * per, 0.0), qvb_nudge=(0.0, -0.5 * per), qvb_hadv=(0.0, 0.0), qvb_vadv=(0.0, 0.0))
        values = [day * bb.DAY, 3.0 / bb.DAY, 1.0 / bb.DAY] + [v for name in names for v in fields[name]]
        np.array(values, dtype='<f4').tofile(case / f'cm1out_diag_{n:06d}_s.dat')
    inflow = bb.read_inflow(case, (case / 'namelist.input').read_text())
    assert inflow['z'] == pytest.approx([500.0, 2000.0]) and inflow['rate'] == pytest.approx([1.0e-5, 0.0])
    assert inflow['theta'] == pytest.approx([301.0, 305.0]) and inflow['qv'] == pytest.approx([0.011, 0.0055])
    w = bb.water('box', 29.5)
    assert w['column_mm_day']['inflow'] == pytest.approx(1.0e-5 * 0.75e-3 * 1000.0 * bb.DAY, rel=1e-4)
    assert w['rest_mm_day'] == pytest.approx(-w['column_mm_day']['inflow'], rel=1e-3)   # 3 - 0.5 - 0.5 - 1 + inflow + rest = 1
    assert bb.read_inflow(case, ' var7 = 0.0,\n') is None and bb.read_inflow(case, ' var11 = 1.0,\n') is None


def test_a_term_only_the_case_holds_counts_as_nil_in_the_reference():
    a = dict(layers_w_m2=dict(radiation={'column': -59.0}, inflow={'column': 22.5}),
             microphysics_layers_w_m2=dict(condensation={'column': 30.0}), surface=dict(t2_k=299.6))
    b = dict(layers_w_m2=dict(radiation={'column': -57.0}), microphysics_layers_w_m2=dict(condensation={'column': 31.0}),
             surface=dict(t2_k=299.3))
    d = bb.difference(a, b)
    assert d['layers_w_m2']['inflow']['column'] == pytest.approx(22.5)
    assert d['layers_w_m2']['radiation']['column'] == pytest.approx(-2.0)
    assert d['microphysics_layers_w_m2']['condensation']['column'] == pytest.approx(-1.0)
    assert d['surface']['t2_k'] == pytest.approx(0.3)
