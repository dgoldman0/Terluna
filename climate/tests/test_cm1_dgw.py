"""Tests of the damped gravity wave's build (cm1_run.DGW_PATCHES, 5 November 2026).

The patch adds code only behind its switch, so with var9 at 0 the build is the inflow build's arithmetic; its wind
(cm1_run.dgw_wind, a copy of the Fortran's arithmetic) is nil when the reference is the box's own mean, rises over a warm
anomaly and sinks under a cool one, and relaxes the first baroclinic mode at the rate var5 sets; the inflow's
bookkeeping keeps heat and water, drawing in the reference's air where the wind converges and sending out the box's
own; and box_budget reads the build's record of its interval means."""
import math

import numpy as np
import pytest

from climate.crm import box_budget as bb
from climate.crm import cm1_run as cr

G = 1.6242188765608128            # the Moon's GM/R^2 (cm1_run.planet), as the moon builds take it
CP, RD, P00 = 1005.7, cr.RD_CM1, 1.0e5


def pristine(name):
    path = cr.CM1_HOME / cr.SOURCE['version'] / 'src' / name
    if not path.exists():
        pytest.skip('the pristine CM1 source is not on this machine')
    return path.read_text(encoding='latin-1')


def patched(name, label):
    text = pristine(name)
    for file, marker, edits in [cr.MAKEFILE[cr.BUILDS[label][1]], *cr.PATCHES, *cr.BUILD_PATCHES[label]]:
        if file == name:
            text = cr.apply_patch(text, file, marker, edits)
    return text


def test_switch_off_is_the_inflow_build():
    for name in ('input.F', 'turb.F', 'solve1.F'):
        undone = patched(name, cr.DGW_BUILD)
        for file, marker, edits in cr.DGW_PATCHES:
            if file == name:
                for old, new, count in reversed(edits):
                    assert undone.count(new) == count
                    undone = undone.replace(new, old)
        assert undone == patched(name, cr.INFLOW_BUILD)
    for block in (cr.DGW_BLOCK, cr.DGW_INFLOW_SUMS, cr.DGW_FLUSH, cr.DGW_TABLE, cr.DGW_HPBL):
        first = next(line for line in block.splitlines() if line.strip() and not line.strip().startswith('!'))
        assert 'var9.gt.0.5' in first                        # every statement it adds runs only with the switch on
    added = [line for _, _, edits in cr.DGW_PATCHES for _, text, _ in edits for line in text.splitlines()]
    assert max(len(line) for line in added) <= 132                          # free form's limit


def test_settings():
    assert cr.dgw_settings(cr.DGW_BUILD, 418.0e3, 40.0, 72.5e3) == dict(var9=1.0, var5=10450.0, var6=40.0, var4=72500.0)
    with pytest.raises(ValueError):
        cr.dgw_settings(cr.INFLOW_BUILD, 418.0e3, 40.0, 72.5e3)
    site = dict(lat_deg=0.0, half_lat_deg=5.6, half_lon_deg=8.5)
    a = cr.planet()['radius_m']
    assert cr.footprint_m(site) == pytest.approx(math.sqrt(4.0 * a * a * math.radians(5.6) * math.radians(8.5)))


def column(n_bv=0.002, theta0=300.0, ps=1.2e5):
    """A column of uniform buoyancy frequency at lunar gravity on the design-years boxes' grid: scalar and w heights,
    the base state's pressure, Exner function and potential temperature on scalar levels, the density on w levels."""
    zw = np.asarray(cr.vertical_grid(**cr.CASES['box_0e_design_years']['z_grid']), float)
    zh = 0.5 * (zw[1:] + zw[:-1])
    theta = lambda z: theta0 * np.exp(n_bv ** 2 * z / G)
    fine = np.linspace(0.0, zw[-1], 200001)
    exner = (ps / P00) ** (RD / CP) - np.concatenate([[0.0], np.cumsum(G / (CP * theta(fine[1:])) * np.diff(fine))])
    pi_at = lambda z: np.interp(z, fine, exner)
    rho_at = lambda z: P00 * pi_at(z) ** (CP / RD) / (RD * theta(z) * pi_at(z))
    return dict(zw=zw, zh=zh, prs0=P00 * pi_at(zh) ** (CP / RD), pi0=pi_at(zh), theta=theta(zh), rf0=rho_at(zw),
                rho=rho_at(zh))


def wind(col, theta_mean, tau=1.0e4, c=40.0, hpbl=0.0, lid=72.5e3, qv_mean=None, qv_ref=None):
    nil = np.zeros_like(col['zh'])
    return cr.dgw_wind(col['zh'], col['zw'], col['prs0'], col['pi0'], col['rf0'], theta_mean,
                       nil if qv_mean is None else qv_mean, col['theta'], nil if qv_ref is None else qv_ref,
                       tau, c, 300.0, hpbl, lid, G)


def test_no_wind_for_the_reference_itself():
    col = column()
    q = 0.012 * np.exp(-col['zh'] / 8000.0)
    assert np.all(wind(col, col['theta'], qv_mean=q, qv_ref=q) == 0.0)


def first_mode(col, lid=72.5e3, floor=300.0):
    """The discrete first baroclinic mode of dgw_wind's operator on the column: omega's shape on the w levels between
    the ground and the lid, and its speed (m/s), from M v = -mu D v with D = Rd S_p / p at the w levels."""
    zw, zh = col['zw'], col['zh']
    nk = zh.size
    pf = np.empty(nk + 1)
    pf[0] = col['prs0'][0] * (col['prs0'][0] / col['prs0'][1]) ** ((zh[0] - zw[0]) / (zh[1] - zh[0]))
    pf[1:nk] = np.sqrt(col['prs0'][:-1] * col['prs0'][1:])
    pf[nk] = col['prs0'][-1]
    lidk = next(k for k in range(1, nk + 1) if zw[k] >= lid)
    i = np.arange(1, lidk)
    h1, h2 = pf[i] - pf[i - 1], pf[i + 1] - pf[i]
    a, c = 2.0 / (h1 * (h1 + h2)), 2.0 / (h2 * (h1 + h2))
    m = np.diag(-(a + c)) + np.diag(c[:-1], 1) + np.diag(a[1:], -1)
    pif = 0.5 * (col['pi0'][i - 1] + col['pi0'][i])
    sp = -pif * (col['theta'][i] - col['theta'][i - 1]) / (col['prs0'][i] - col['prs0'][i - 1])   # -(T/theta) dtheta/dp
    d = RD * sp / pf[i]
    mu, vec = np.linalg.eig(np.linalg.solve(np.diag(d), -m))
    k = int(np.argmin(mu.real))
    return dict(i=i, v=vec[:, k].real, c=1.0 / math.sqrt(mu[k].real), sp=sp, pf=pf, lid=lidk)


def test_warm_anomaly_rises_at_the_derived_rate():
    col = column()
    mode = first_mode(col)
    i, sp = mode['i'], mode['sp']
    t_w = sp * mode['v']                                    # the temperature anomaly the mode's omega answers
    t_w *= 0.5 / t_w[np.argmax(np.abs(t_w))]                # a warm anomaly peaking at 0.5 K
    t_scalar = np.interp(col['zh'], col['zw'][i], t_w, left=0.0, right=0.0)
    tau = 1.0e4
    for sign in (1.0, -1.0):
        w = wind(col, col['theta'] + sign * t_scalar / col['pi0'], tau=tau, c=mode['c'])
        peak = i[np.argmax(np.abs(t_w))]
        assert np.sign(w[peak]) == sign                     # warm rises, cool sinks
        tendency = -col['rf0'][i] * G * w[i] * sp           # dT'/dt = omega S_p, omega = -rho g w
        rate = -tendency[np.argmax(np.abs(t_w))] / (sign * t_w[np.argmax(np.abs(t_w))])
        assert rate == pytest.approx(1.0 / tau, rel=0.05)
    assert 20.0 < mode['c'] < 120.0                         # a lunar column's first mode, about Earth's speed


def test_the_inflow_keeps_heat_and_water():
    col = column()
    rng = np.random.default_rng(5)
    anomaly = 0.3 * np.sin(np.pi * np.clip(col['zh'] / 72.5e3, 0.0, 1.0)) + 0.05 * rng.standard_normal(col['zh'].size)
    w = wind(col, col['theta'] + anomaly / col['pi0'], hpbl=2000.0)
    dz = np.diff(col['zw'])
    div = (col['rf0'][1:] * w[1:] - col['rf0'][:-1] * w[:-1]) / dz         # horizontal mass convergence, kg/m3/s
    assert abs((div * dz).sum()) < 1e-12 * np.abs(div * dz).sum()           # the air drawn in equals the air sent out
    for a, a_ref in ((col['theta'] + anomaly, col['theta']), (0.012 * np.exp(-col['zh'] / 8e3),
                                                              0.013 * np.exp(-col['zh'] / 9e3))):
        face = np.r_[0.0, 0.5 * (a[1:] + a[:-1]), 0.0]                      # any face values: the flux part telescopes
        flux = col['rf0'] * w * face
        advect = (-(flux[1:] - flux[:-1]) / dz + a * div) / col['rho']       # CM1's wsub: the flux form less a times div
        inflow = np.maximum(div, 0.0) / col['rho'] * (a_ref - a)            # INFLOW_BLOCK at its rate max(div, 0)/rho0
        gained = ((advect + inflow) * col['rho'] * dz).sum()
        expected = (dz * np.where(div > 0.0, div * a_ref, div * a)).sum()  # the reference's air in, the box's out
        assert gained == pytest.approx(expected, rel=1e-9, abs=1e-12 * np.abs(div * dz * a).sum())
        same = ((-(col['rf0'][1:] * w[1:] - col['rf0'][:-1] * w[:-1]) / dz + div) * dz).sum()
        assert abs(same) < 1e-12 * np.abs(div * dz).sum()                   # uniform air: nothing gained or lost


def test_box_budget_reads_the_record(tmp_path):
    (tmp_path / 'namelist.input').write_text(' &param8\n var9 = 1.0,\n var7 = 1.0,\n /\n &param14\n diagfrq = 10800.0,\n /\n')
    rows = []
    for end, w in ((10800.0, 0.01), (21600.4, -0.02)):
        rows.append(f'   {end:.10E}   {270.0:.7E}   {1500.0:.7E}           2')
        for zh, zf in ((50.0, 0.0), (150.0, 100.0)):
            rows.append(f'   {zh:.7E}   {zf:.7E}   {w * zf / 100.0:.7E}   {1e-5:.7E}   {2e-8:.7E}   {3e-6:.7E}')
    (tmp_path / bb.DGW_RECORD).write_text('\n'.join(rows) + '\n')
    record = bb.read_dgw(tmp_path, (tmp_path / 'namelist.input').read_text())
    assert [r['end_s'] for r in record] == [10800.0, 21600.4]
    assert record[1]['w'].tolist() == pytest.approx([0.0, -0.02])
    assert record[0]['inflow_qv'].tolist() == pytest.approx([2e-8, 2e-8])
    assert record[0]['table_theta'].tolist() == pytest.approx([3e-6, 3e-6])
    assert bb.dgw_at(record, 21600.0, 10800.0)['zb_m'] == 1500.0
    with pytest.raises(ValueError):
        bb.dgw_at(record, 40000.0, 10800.0)
    assert bb.read_dgw(tmp_path, ' &param8\n var9 = 0.0,\n /\n') is None
