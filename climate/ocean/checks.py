"""The ocean models against geography's linear tide solver.

geography/tide_dynamics.py solves the linear shallow-water equations on the
atlas's quarter-degree cells, with the Moon's rotation and a linear bed
friction: it gives each sea's seiches and its response to each tidal line.
The same cells, depths and settings in an ocean model should give the same
answers, which tests the model's lunar gravity, radius, rotation and tide
forcing together:

- a channel on the equator, whose slowest seiche is 2L / sqrt(g H);
- Smythii-Marginis's three slowest seiches (geography: 40, 21 and 17 hours),
  from a free oscillation started in the solver's own mode shapes, without
  rotation or friction (SCHISM also on a uniform floor);
- MITgcm: Mare Fecunditatis's response to the 9.61-day line 2D + M' with ten
  and a hundred times the base friction (geography: 5.4 and 0.6 times the
  equilibrium height over the mare), and the Near-side Sea's response, with
  the hundredfold case again at half the time step and record interval;
- SCHISM: Smythii-Marginis under the monthly line M' through the lunar
  patch's tide reader.

`python -m climate.ocean.checks mitgcm` runs MITgcm's checks (about half an hour,
the forced cases on four processes) and writes results/ocean/checks_mitgcm.json
under climate/; `schism` runs SCHISM's. Each check starts from a clear folder.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from climate.ocean import forcing, lunar, mitgcm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "climate" / "results" / "ocean"
NEARSIDE, SMYTHII = 433, 874


def quarter_degree(body):
    """The solver's sea (its largest face-joined part) on the atlas cells, as an MITgcm grid and the solver's SeaGrid."""
    from geography import tide_dynamics as dynamics, tides
    g = tides.load_seas()
    sea = dynamics.SeaGrid(g["labels"], g["depth"], g["lat"], g["lon"], body)
    lat, lon = g["lat"], np.where(g["lon"] > 180, g["lon"] - 360, g["lon"])
    rows, cols = np.nonzero(sea.mask)
    lon_c = lon[cols]
    if lon_c.max() - lon_c.min() > 180:          # a sea astride 180 E would need the other convention
        raise ValueError("sea straddles the 180 E meridian")
    j0, j1 = rows.min() - 2, rows.max() + 2      # rows run north to south in the atlas
    lon_lo, lon_hi = lon_c.min() - 0.5, lon_c.max() + 0.5
    col_order = np.argsort(lon)
    keep = col_order[(lon[col_order] > lon_lo) & (lon[col_order] < lon_hi)]
    depth = np.where(sea.mask, np.maximum(g["depth"], dynamics.MIN_DEPTH_M), 0.0)
    sub = depth[j0:j1 + 1][:, keep][::-1]       # south first
    lat_sub = lat[j0:j1 + 1][::-1]
    lon_sub = lon[keep]
    xg = np.r_[lon_sub - 0.125, lon_sub[-1] + 0.125]
    yg = np.r_[lat_sub - 0.125, lat_sub[-1] + 0.125]
    index = dict(rows=(lat_sub, lon_sub))
    grid = mitgcm.Grid(xg=xg, yg=yg, depth=sub, name=f"atlas body {body}")
    # each solver cell's (j, i) on the MITgcm grid, from the atlas indices (rows reversed, columns in `keep` order)
    position = np.full(len(lon), -1)
    position[keep] = np.arange(len(keep))
    jj, ii = j1 - sea.rows, position[sea.cols]
    if (ii < 0).any() or not (np.array_equal(lat_sub[jj], lat[sea.rows]) and np.array_equal(lon_sub[ii], lon[sea.cols])
                              and (sub[jj, ii] > 0).all()):
        raise RuntimeError("solver cells do not map onto the grid's water cells")
    return grid, sea, (jj, ii), index


def fit_lines(t_s, series, frequencies_cpd):
    """Least-squares mean and complex amplitudes Re[A exp(i w t)] of `series` (time, ...) at the frequencies."""
    w = 2 * np.pi * np.asarray(frequencies_cpd) / 86400.0
    design = np.column_stack([np.ones_like(t_s)] + [np.cos(x * t_s) for x in w] + [np.sin(x * t_s) for x in w])
    sol, *_ = np.linalg.lstsq(design, series.reshape(len(t_s), -1), rcond=None)
    k = len(w)
    amp = sol[1:k + 1] - 1j * sol[k + 1:2 * k + 1]
    return sol[0].reshape(series.shape[1:]), amp.reshape((k,) + series.shape[1:])


def dominant_frequency(t_s, x, guess_cpd):
    """A single oscillation's frequency: the least-squares fit of a mean and one sinusoid with the smallest residual,
    within 20% of a guess, by a scan finer than the record's resolution refined there. (The fitted amplitude's peak
    sits slightly off the true frequency in a finite record, by about 0.1% over a dozen cycles; the residual's
    minimum does not.)"""
    from scipy.optimize import minimize_scalar
    def cost(f):
        mean, a = fit_lines(t_s, x, [f])
        w = 2 * np.pi * f / 86400.0
        fitted = mean + np.real(a[0] * np.exp(1j * w * t_s))
        return float(np.sum((x - fitted) ** 2))
    span = (t_s[-1] - t_s[0]) / 86400.0
    scan = np.arange(guess_cpd * 0.8, guess_cpd * 1.2, 0.1 / span)
    best = scan[int(np.argmin([cost(f) for f in scan]))]
    r = minimize_scalar(cost, bounds=(best - 0.1 / span, best + 0.1 / span), method="bounded", options=dict(xatol=1e-8))
    return r.x


def seiches_mitgcm(count=3):
    from geography import tide_dynamics as dynamics
    grid, sea, (jj, ii), _ = quarter_degree(SMYTHII)
    periods_h, vectors = dynamics.seiche_periods(sea, count)
    eta = np.zeros_like(grid.depth)
    shapes = []
    for k in range(count):
        v = vectors[:, k] / np.max(np.abs(vectors[:, k]))
        shapes.append(v)
        eta[jj, ii] += 0.1 * v                       # 10 cm in each mode
    dt = 120.0
    exp = mitgcm.Experiment(name="check_seiches_smythii", grid=grid, deltat=dt, rotation=False, tide=None,
                            output_s=900.0, output_fields=("ETAN",), implicit=0.5, initial=dict(eta=eta))
    exp.write(fresh=True)
    days = 14.0
    exp.run(days * 86400.0, days * 86400.0)
    its = exp.snapshots()
    t = np.array([i * dt for i in its])
    series = np.array([exp.state(i)["ETAN"][jj, ii] for i in its])          # (time, cells)
    w = sea.area / sea.area.sum()
    out = []
    for k in range(count):
        modal = series @ (w * shapes[k]) / np.sum(w * shapes[k] ** 2)
        f = dominant_frequency(t, modal, 24.0 / periods_h[k])
        out.append(dict(mode=k + 1, solver_hours=float(periods_h[k]), mitgcm_hours=float(24.0 / f),
                        relative_difference=float(24.0 / f / periods_h[k] - 1)))
    return dict(cells=int(sea.cells), deltat_s=dt, days=days, modes=out,
                mean_level_range_m=monitor_range(exp.folder / "run.log", "dynstat_eta_mean"))


def monitor_range(log, name):
    """The spread of one of MITgcm's monitor statistics over a run (max minus min)."""
    import re
    values = [float(v) for v in re.findall(rf"%MON {name}\s*=\s*(\S+)", Path(log).read_text())]
    return float(max(values) - min(values)) if values else None


def fecunditatis_region(sea):
    from geography import nomenclature
    named = nomenclature.by_name(nomenclature.features())
    centre = named["Mare Fecunditatis"]
    pts = forcing.unit_points(np.degrees(sea.lat), np.degrees(sea.lon))
    c = forcing.unit_points(centre["lat"], centre["lon"])
    angle = np.arccos(np.clip(pts @ c, -1, 1))
    return angle * lunar.RADIUS / 1000 < centre["diameter_km"] / 2


def forced_mitgcm(argument="2D +M'", cases=((10, 110.0, 50.0, 1800.0, 21600.0), (100, 60.0, 30.0, 1800.0, 21600.0),
                                            (100, 60.0, 30.0, 900.0, 10800.0))):
    """MITgcm against the solver for one line; each case is (friction factor, days run, days fitted at the end, time
    step, forcing record interval). The third case halves the step and the record interval: the forcing's linear
    interpolation between records carries a line at 1 - (w T)^2 / 12 of its amplitude, and the step's timing lags it by
    up to half a step, so both residuals should fall with it."""
    from geography import tide_dynamics as dynamics
    entry, amp = forcing.line(argument)
    f_cpd = entry["frequency_cpd"]
    grid, sea, (jj, ii), _ = quarter_degree(NEARSIDE)
    grid = mitgcm.pad_even(grid)
    region = fecunditatis_region(sea)
    w_all = sea.area / sea.area.sum()
    values = forcing.monomials(np.degrees(sea.lat), np.degrees(sea.lon))
    f_cells = values @ amp
    eq = sea.equilibrium(f_cells)
    spec = mitgcm.tide_forcing_grid(grid, step=1.0)
    out = []
    for factor, days, fit_days, deltat, period in cases:
        seconds = np.arange(0.0, days * 86400.0 + period, period)
        coefs = forcing.line_series(seconds, amp, f_cpd) * forcing.ramp(seconds, 10.0)[:, None]
        r = factor * dynamics.FRICTION_M_S
        solver = dynamics.forced_response(sea, f_cells, f_cpd, friction=r)
        suffix = "" if deltat == 1800.0 else f"_dt{deltat:g}"
        exp = mitgcm.Experiment(name=f"check_forced_{argument.replace(' ', '').replace(chr(39), 'p')}_x{factor}{suffix}",
                                grid=grid, deltat=deltat, linear_friction=r, rotation=True,
                                tide=dict(spec, period_s=period, file="tide.bin"), output_s=21600.0,
                                output_fields=("ETAN",), output_from_s=(days - fit_days) * 86400.0, implicit=0.5,
                                tiles=(2, 2))
        exp.write(fresh=True)
        forcing.write_mitgcm(exp.folder / "tide.bin", spec, seconds, coefs, attraction=False)
        exp.run(days * 86400.0, (days - fit_days) * 86400.0)
        its = [i for i in exp.snapshots() if i * exp.deltat >= (days - fit_days) * 86400.0 - 1]
        t = np.array([i * exp.deltat for i in its])
        series = np.array([exp.state(i)["ETAN"][jj, ii] for i in its])
        _, a = fit_lines(t, series, [f_cpd])
        model = a[0]
        np.savez_compressed(exp.folder / "fit.npz", model=model, solver=solver, equilibrium=eq, area=sea.area,
                            depth=sea.depth, lat_deg=np.degrees(sea.lat), lon_deg=np.degrees(sea.lon), fecunditatis=region)
        def ratio(x, m):
            wr = sea.area[m]
            return float(np.sqrt(np.sum(wr * np.abs(x[m]) ** 2) / np.sum(wr * np.abs(eq[m]) ** 2)))
        whole = np.ones(sea.cells, bool)
        out.append(dict(friction_factor=factor, friction_m_s=r, days=days, fit_days=fit_days, deltat_s=deltat,
                        forcing_record_s=period, forcing_interpolation_amplitude=float(
                            1 - (2 * math.pi * f_cpd / 86400.0 * period) ** 2 / 12),
                        fecunditatis_ratio=dict(solver=ratio(solver, region), mitgcm=ratio(model, region)),
                        nearside_ratio=dict(solver=ratio(solver, whole), mitgcm=ratio(model, whole)),
                        fecunditatis_phase_deg=dict(
                            solver=float(np.degrees(np.angle(np.sum(sea.area[region] * solver[region]) / np.sum(sea.area[region] * eq[region])))),
                            mitgcm=float(np.degrees(np.angle(np.sum(sea.area[region] * model[region]) / np.sum(sea.area[region] * eq[region]))))),
                        rms_difference_over_rms_solver=float(np.sqrt(np.sum(w_all * np.abs(model - solver) ** 2) /
                                                                     np.sum(w_all * np.abs(solver) ** 2)))))
    return dict(line=dict(argument=entry["argument"], period_days=entry["period_days"]), cells=int(sea.cells),
                fecunditatis_cells=int(region.sum()), results=out)


def cells_and_nodes(grid, mesh):
    """Sparse maps between the grid's cells and a from_cells mesh's nodes: node = mean of its wet cells,
    cell = mean of its four corner nodes."""
    import scipy.sparse as sparse
    step = grid.yg[1] - grid.yg[0]
    lon_step = grid.xg[1] - grid.xg[0]
    jn = np.rint((mesh.lat - grid.yg[0]) / step).astype(int)
    inn = np.rint((mesh.lon - grid.xg[0]) / lon_step).astype(int)
    ny, nx = grid.depth.shape
    wet = grid.depth > 0
    rows, cols = [], []
    for n, (j, i) in enumerate(zip(jn, inn)):
        for jj in (j - 1, j):
            for ii in (i - 1, i):
                if 0 <= jj < ny and 0 <= ii < nx and wet[jj, ii]:
                    rows.append(n)
                    cols.append(jj * nx + ii)
    incidence = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(mesh.nodes, ny * nx))
    node_from_cell = sparse.diags(1.0 / np.maximum(incidence.sum(axis=1).A1, 1)) @ incidence
    cell_from_node = sparse.diags(1.0 / np.maximum(incidence.sum(axis=0).A1, 1)) @ incidence.T
    return node_from_cell, cell_from_node


def seiches_schism(count=3, days=14.0, uniform_depth=None):
    """SCHISM's free oscillation of Smythii-Marginis; with `uniform_depth` every cell takes that depth, so the
    solver's face depths and SCHISM's corner depths agree and only the lunar constants and the numerics differ."""
    from geography import tide_dynamics as dynamics
    from climate.ocean import mesh as meshes, schism
    grid, sea, (jj, ii), _ = quarter_degree(SMYTHII)
    if uniform_depth:
        grid.depth = np.where(grid.depth > 0, float(uniform_depth), 0.0)
        sea.depth[:] = uniform_depth
        sea.u_depth[:] = uniform_depth
        sea.v_depth[:] = uniform_depth
    periods_h, vectors = dynamics.seiche_periods(sea, count)
    m = meshes.from_cells(grid.depth, grid.yc, grid.xc)
    node_from_cell, cell_from_node = cells_and_nodes(grid, m)
    eta = np.zeros_like(grid.depth)
    shapes = []
    for k in range(count):
        v = vectors[:, k] / np.max(np.abs(vectors[:, k]))
        shapes.append(v)
        eta[jj, ii] += 0.1 * v
    exp = schism.Experiment(name="check_seiches_smythii" + (f"_uniform{int(uniform_depth)}" if uniform_depth else ""), mesh=m, dt=120.0, drag=0.0, advection=False, rotation=False,
                            output_s=900.0, file_s=86400.0, initial_eta=node_from_cell @ eta.ravel(), compute_ranks=2)
    exp.write(fresh=True)
    exp.run(days * 86400.0, days * 86400.0)
    t, fields = exp.series(("elevation",))
    cells = np.array([(cell_from_node @ e).reshape(grid.depth.shape)[jj, ii] for e in fields["elevation"]])
    w = sea.area / sea.area.sum()
    out = []
    for k in range(count):
        modal = cells @ (w * shapes[k]) / np.sum(w * shapes[k] ** 2)
        f = dominant_frequency(t, modal, 24.0 / periods_h[k])
        out.append(dict(mode=k + 1, solver_hours=float(periods_h[k]), schism_hours=float(24.0 / f),
                        relative_difference=float(24.0 / f / periods_h[k] - 1)))
    return dict(cells=int(sea.cells), nodes=m.nodes, elements=m.elements, deltat_s=exp.dt, days=days, modes=out,
                uniform_depth_m=uniform_depth,
                reading_rule=("Every cell at one depth: the solver's face depths and SCHISM's corner depths are the same "
                              "floor, so the periods test gravity, radius and the numerics." if uniform_depth else
                              "SCHISM's depths sit at the cell corners (the mean of the cells sharing each corner) and vary "
                              "linearly across each triangle; the solver gives each face the shallower cell's depth, which "
                              "slows waves across Smythii's craters. Its periods with the mean or the deeper cell's depth at "
                              "each face are 34.9 and 32.5 hours for the slowest mode, bracketing SCHISM's."))


def forced_schism(argument="M'", body=SMYTHII, days=110.0, fit_days=55.0):
    """SCHISM under one monthly line through the lunar patch's tide reader (ilunartide = 1), without bed friction,
    against the solver with its base friction: Smythii-Marginis follows its equilibrium within 1.4% at this line,
    so the gain and phase test the forcing's amplitude, pattern and sign."""
    from geography import tide_dynamics as dynamics
    from climate.ocean import mesh as meshes, schism
    entry, amp = forcing.line(argument)
    f_cpd = entry["frequency_cpd"]
    grid, sea, (jj, ii), _ = quarter_degree(body)
    m = meshes.from_cells(grid.depth, grid.yc, grid.xc)
    _, cell_from_node = cells_and_nodes(grid, m)
    values = forcing.monomials(np.degrees(sea.lat), np.degrees(sea.lon))
    f_cells = values @ amp
    eq = sea.equilibrium(f_cells)
    solver = dynamics.forced_response(sea, f_cells, f_cpd)
    period = 10800.0
    seconds = np.arange(0.0, days * 86400.0 + period, period)
    coefs = forcing.line_series(seconds, amp, f_cpd) * forcing.ramp(seconds, 10.0)[:, None]
    exp = schism.Experiment(name=f"check_forced_{argument.replace(' ', '').replace(chr(39), 'p')}_{body}", mesh=m,
                            dt=600.0, drag=0.0, advection=False, rotation=True, output_s=10800.0,
                            file_s=(days - fit_days) * 86400.0 / 5, compute_ranks=2)
    exp.write(seconds, coefs, attraction=False, fresh=True)
    exp.run(days * 86400.0, (days - fit_days) * 86400.0)
    t, fields = exp.series(("elevation",))
    keep = t >= (days - fit_days) * 86400.0 - 1
    cells = np.array([(cell_from_node @ e).reshape(grid.depth.shape)[jj, ii] for e in fields["elevation"][keep]])
    _, a = fit_lines(t[keep], cells, [f_cpd])
    model = a[0]
    w = sea.area / sea.area.sum()
    def regress(x, y):
        g = np.sum(w * x * np.conj(y)) / np.sum(w * np.abs(y) ** 2)
        return dict(gain=float(abs(g)), phase_deg=float(np.degrees(np.angle(g))))
    return dict(line=dict(argument=entry["argument"], period_days=entry["period_days"]), body=body, days=days,
                fit_days=fit_days, nodes=m.nodes, elements=m.elements,
                schism_on_solver=regress(model, solver), schism_on_equilibrium=regress(model, eq),
                solver_on_equilibrium=regress(solver, eq),
                reading_rule="Complex least-squares gain of one response on another over the sea's cells (area-weighted); "
                             "a sign error would show as a 180-degree phase, a wrong pattern as a gain far from one.")


def channel(model):
    """A 10-degree channel on the equator, 1,000 m deep and 2 degrees wide: its slowest seiche is 2L / sqrt(g H)
    in the continuum, which tests a model's lunar gravity and radius with nothing else in play."""
    from geography import tide_dynamics as dynamics
    from climate.ocean import mesh as meshes, schism
    lat = np.arange(1.875, -2.0, -0.25)
    lon = np.arange(-0.875, 11.0, 0.25)
    inside = (np.abs(lat) < 1.0)[:, None] & ((lon > 0) & (lon < 10))[None, :]
    depth = np.where(inside, 1000.0, 0.0)
    sea = dynamics.SeaGrid(inside.astype(int), depth, lat, lon, 1)
    periods_h, vectors = dynamics.seiche_periods(sea, 1)
    analytic_h = 2 * lunar.RADIUS * math.radians(10.0) / math.sqrt(lunar.GRAVITY * 1000.0) / 3600.0
    yc, xc, d = lat[::-1], lon, depth[::-1]
    grid = mitgcm.Grid(xg=np.r_[xc - 0.125, xc[-1] + 0.125], yg=np.r_[yc - 0.125, yc[-1] + 0.125], depth=d)
    jj, ii = len(lat) - 1 - sea.rows, sea.cols           # rows reversed (south first); columns as they are
    v = vectors[:, 0] / np.abs(vectors[:, 0]).max()
    eta = np.zeros_like(d)
    eta[jj, ii] = 0.1 * v
    if model == "mitgcm":
        exp = mitgcm.Experiment(name="check_channel", grid=grid, deltat=30.0, rotation=False, output_s=300.0,
                                output_fields=("ETAN",), initial=dict(eta=eta))
        exp.write(fresh=True)
        exp.run(2 * 86400.0, 2 * 86400.0)
        its = exp.snapshots()
        t = np.array([i * exp.deltat for i in its])
        series = np.array([exp.state(i)["ETAN"][jj, ii] for i in its])
    else:
        m = meshes.from_cells(d, yc, xc)
        node_from_cell, cell_from_node = cells_and_nodes(grid, m)
        exp = schism.Experiment(name="check_channel", mesh=m, dt=30.0, drag=0.0, advection=False, rotation=False,
                                output_s=300.0, file_s=86400.0, initial_eta=node_from_cell @ eta.ravel(), compute_ranks=2)
        exp.write(fresh=True)
        exp.run(2 * 86400.0, 2 * 86400.0)
        t, fields = exp.series(("elevation",))
        series = np.array([(cell_from_node @ e).reshape(d.shape)[jj, ii] for e in fields["elevation"]])
    w = sea.area / sea.area.sum()
    modal = series @ (w * v) / np.sum(w * v ** 2)
    f = dominant_frequency(t, modal, 24.0 / periods_h[0])
    return dict(analytic_hours=analytic_h, solver_hours=float(periods_h[0]), model_hours=float(24.0 / f),
                model_over_analytic=float(24.0 / f / analytic_h - 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("model", choices=["mitgcm", "schism"])
    parser.add_argument("--only", choices=["seiches", "forced"])
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f"checks_{args.model}.json"
    record = json.loads(path.read_text()) if path.exists() else dict(schema="terluna.climate.ocean-checks/1", model=args.model)
    if args.only in (None, "seiches"):
        record["channel"] = channel(args.model)
        path.write_text(json.dumps(record, indent=1) + "\n")
        record["seiches_smythii"] = seiches_mitgcm() if args.model == "mitgcm" else seiches_schism()
        if args.model == "schism":
            record["seiches_smythii_uniform"] = seiches_schism(uniform_depth=1000.0)
        path.write_text(json.dumps(record, indent=1) + "\n")
    if args.only in (None, "forced"):
        if args.model == "mitgcm":
            record["forced_fecunditatis"] = forced_mitgcm()
        else:
            record["forced_monthly_smythii"] = forced_schism()
        path.write_text(json.dumps(record, indent=1) + "\n")
    print(path)


if __name__ == "__main__":
    main()
