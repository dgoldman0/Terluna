"""Isolate an existing atmospheric restart and recover instantaneous wave inputs.

The copied executable and physical inputs stay fixed. Output cadence changes
in the copy; the original climate archive and installed model stay intact.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time

import numpy as np

from climate.gcm.wave_coverage import ROOT, digest
from shared.constants import DATA, MOON_SURFACE_GRAVITY

DEFAULT_SOURCE = ROOT / "climate/gcm/runs/A28_dim5_moon"
DEFAULT_OUTPUT = ROOT / "research/runs/waves/weather"


def edit_namelist(text, settings):
    """Require exactly one existing assignment for every requested setting."""
    for key, value in settings.items():
        pattern = r"(?im)^(\s*" + re.escape(key) + r"\s*=\s*)[^\n]+"
        text, count = re.subn(pattern, lambda m: m[1] + str(value) + " ", text)
        if count != 1:
            raise ValueError(f"Expected one namelist assignment for {key}; found {count}")
    return text


def scalar(text, key):
    matches = re.findall(r"(?im)^\s*" + re.escape(key) + r"\s*=\s*([\d.eEdD+-]+)", text)
    if len(matches) != 1:
        raise ValueError(f"Expected one value for {key}")
    return float(matches[0].replace("d", "e").replace("D", "e"))


def run(source=DEFAULT_SOURCE, output=DEFAULT_OUTPUT, *, steps=1446, cadence=6, snapshots=True):
    """Run the saved eight-rank binary, using copied inputs and a bounded duration."""
    if cadence < 1 or steps < cadence or steps % cadence:
        raise ValueError("Steps must be a positive multiple of the snapshot cadence")
    source, output = source.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(f"Run folder already exists; retained for inspection: {output}")
    work = source / "model"
    executable = work / "most_plasim_t21_l10_p8.x"
    restart = work / "MOST_REST.00029"
    inputs = sorted(work.glob("*_namelist")) + sorted(work.glob("*.sra"))
    inputs += [work / "sun.dat", work / "sun_hr.dat", executable, restart]
    inputs += [work / name for name in ("restart_dsnow", "restart_xsnow") if (work / name).is_file()]
    original = {str(p): digest(p) for p in inputs}
    config = json.loads((source / "progress.json").read_text())["configuration"]
    if config["model"]["resolution"] != "T21" or config["ncpus"] != 8:
        raise ValueError("This runner uses the archived T21 eight-rank configuration")
    planet = (work / "planet_namelist").read_text()
    if not np.isclose(scalar(planet, "GA"), MOON_SURFACE_GRAVITY, rtol=1e-10):
        raise ValueError("Restart gravity differs from the lunar scenario")
    output.mkdir(parents=True)
    for p in inputs:
        shutil.copy2(p, output / ("plasim_restart" if p == restart else p.name))
    shutil.copy2(Path(__file__), output / "producer.py")
    settings = dict(N_RUN_STEPS=steps, N_RUN_YEARS=0, N_RUN_MONTHS=0,
                    N_RUN_DAYS=0, NSNAPSHOT=int(snapshots), NSTPS=cadence)
    namelist = output / "plasim_namelist"
    namelist.write_text(edit_namelist(namelist.read_text(), settings))
    identity = {p.name: digest(p) for p in output.iterdir()}
    manifest = dict(schema="terluna.climate.wave-snapshot-run/1", state="prepared",
                    producer_sha256=digest(Path(__file__)), source_inputs=original,
                    copied_inputs=identity, source_configuration=config, output_settings=settings,
                    timestep_s=config["model"]["timestep_min"] * 60,
                    source_progress_sha256=digest(source / "progress.json"),
                    source_restart=str(restart), mpi_ranks=8, timeout_s=600)
    record = output / "run.json"
    record.write_text(json.dumps(manifest, indent=2) + "\n")
    start = time.monotonic()
    with (output / "console.log").open("w") as log:
        process = subprocess.Popen(["mpiexec", "-np", "8", "./" + executable.name], cwd=output,
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                   env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"))
        try:
            returncode = process.wait(timeout=600)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise
    required = ["plasim_output", "plasim_status", "plasim_diag"]
    if snapshots:
        required += ["plasim_snapshot"]
    if returncode or (output / "Abort_Message").exists() or any(
            not (output / name).is_file() for name in required):
        raise RuntimeError(f"Atmospheric run failed; inspect {output}")
    if any(digest(Path(p)) != value for p, value in original.items()):
        raise RuntimeError("Original climate inputs changed during the run")
    manifest.update(state="complete", elapsed_wall_s=time.monotonic() - start,
                    original_inputs_unchanged=True,
                    output_sha256={name: digest(output / name) for name in required})
    record.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def export(run_dir=DEFAULT_OUTPUT / "snapshots", output=DEFAULT_OUTPUT / "atmosphere.npz"):
    """Decode spectra with the installed model's postprocessor, retaining native grids."""
    from exoplasim import pyburn

    producer = run_dir / "export_producer.py"
    shutil.copy2(Path(__file__), producer)
    producer_sha256 = digest(producer)
    record = json.loads((run_dir / "run.json").read_text())
    if record["state"] != "complete" or record["output_settings"]["NSNAPSHOT"] != 1:
        raise ValueError("A completed snapshot run is required")
    for name, value in record["output_sha256"].items():
        if digest(run_dir / name) != value:
            raise ValueError(f"Changed atmospheric output: {name}")
    planet = (run_dir / "planet_namelist").read_text()
    gascon, radius = scalar(planet, "GASCON"), scalar(planet, "PLARAD")
    # pyburn's API takes radius in its own reference Earth radii. This cancels
    # the reference-radius conversion in that installed implementation.
    reference_radius = DATA["model_closures"]["exoplasim_3_4_2"]["postprocessor_reference_radius_m"]
    data = pyburn.dataset(str(run_dir / "plasim_snapshot"),
                          [str(n) for n in (130, 131, 132, 133, 134, 139, 159, 172, 173, 180, 181, 210)],
                          radius=radius / reference_radius, gravity=scalar(planet, "GA"), gascon=gascon,
                          logfile=str(run_dir / "decode.log"))
    # Retain the decoded intermediate so later audits can inspect units and
    # coordinates independently of this exporter's field selection.
    np.savez_compressed(run_dir / "decoded.npz", **{k: np.asarray(v[0]) for k, v in data.items()},
                        metadata=np.array(json.dumps({k: v[1] for k, v in data.items()})))
    get = lambda key: np.asarray(data[key][0])
    steps = get("time")
    cadence = record["output_settings"]["NSTPS"]
    expected = np.arange(cadence - 1, record["output_settings"]["N_RUN_STEPS"], cadence)
    if not np.array_equal(steps, expected):
        raise ValueError("Unexpected snapshot times; preserve and inspect the step convention")
    # Resolve model-native names through the installed code table, including
    # friction velocity cubed, roughness and the signed surface stresses.
    code = lambda n: get(pyburn.ilibrary[str(n)][0])
    arrays = dict(time_s=(steps - steps[0]) * record["timestep_s"],
                  restart_elapsed_s=steps * record["timestep_s"],
                  longitude_deg=get("lon"), latitude_deg=get("lat"),
                  eastward_lowest_m_s=get("ua")[:, -1], northward_lowest_m_s=get("va")[:, -1],
                  temperature_lowest_K=get("ta")[:, -1], specific_humidity_lowest=get("hus")[:, -1],
                  surface_pressure_Pa=code(134) * 100., surface_temperature_K=code(139),
                  friction_velocity_cubed_m3_s3=code(159), land_fraction=code(172),
                  roughness_m=code(173), stress_eastward_Pa=code(180), stress_northward_Pa=code(181),
                  sea_ice_fraction=code(210))
    # pyburn labels surface pressure hPa after dividing the raw Pa by 100.
    if data[pyburn.ilibrary["134"][0]][1][2] != "hPa":
        raise ValueError("Unexpected postprocessor pressure units")
    shape = (len(steps), len(arrays["latitude_deg"]), len(arrays["longitude_deg"]))
    for key, values in arrays.items():
        if not np.isfinite(values).all() or (values.ndim > 1 and values.shape != shape):
            raise ValueError(f"Invalid atmospheric field: {key}")
    if np.any(arrays["surface_pressure_Pa"] < 1000) or np.any(arrays["roughness_m"] <= 0):
        raise ValueError("Invalid pressure or roughness")
    rho = arrays["surface_pressure_Pa"] / (gascon * arrays["surface_temperature_K"])
    arrays["stress_closure_density_kg_m3"] = rho
    stress = np.hypot(arrays["stress_eastward_Pa"], arrays["stress_northward_Pa"])
    recovered = rho * np.maximum(arrays["friction_velocity_cubed_m3_s3"], 0)**(2/3)
    relative = np.abs(stress - recovered) / np.maximum(stress, 1e-8)
    # Fluxmod uses the previous physics surface temperature for this closure;
    # the snapshot temperature is after the surface step. Keep the residual.
    metadata = dict(schema="terluna.climate.gcm-wave-snapshots/1",
                    evidence="Instantaneous atmospheric fields from a copied settled GCM restart, with unchanged climate physics.",
                    reading_rule=f"Global native T21 grid, sampled every {cadence * record['timestep_s']/3600:g} Earth hours. Winds are at the lowest sigma level. Surface stress is the model's downward momentum transfer; its density closure uses surface temperature and the configured dry-gas constant. Wave forcing requires an explicit surface coupling rule.",
                    producer=dict(domain="climate", source_sha256=producer_sha256,
                                  postprocessor_sha256=digest(Path(pyburn.__file__))),
                    run_record_sha256=digest(run_dir / "run.json"), raw_sha256=record["output_sha256"]["plasim_snapshot"],
                    source_configuration=record["source_configuration"], sigma_lowest=float(get("lev")[-1]),
                    gas_constant_J_kg_K=gascon, gravity_m_s2=scalar(planet, "GA"), radius_m=radius,
                    sample_count=len(steps), sample_interval_s=cadence * record["timestep_s"],
                    duration_s=float(arrays["time_s"][-1]), first_snapshot_restart_step=float(steps[0]),
                    pressure_conversion="pyburn hPa multiplied by 100", stress_closure_relative_max=float(relative.max()),
                    units=dict(time_s="s since first snapshot",restart_elapsed_s="model steps since restart multiplied by timestep in seconds",
                               longitude_deg="degrees east",latitude_deg="degrees north",
                               eastward_lowest_m_s="m/s",northward_lowest_m_s="m/s",
                               temperature_lowest_K="K",specific_humidity_lowest="kg/kg",
                               surface_pressure_Pa="Pa",surface_temperature_K="K",
                               friction_velocity_cubed_m3_s3="m^3/s^3",land_fraction="fraction",
                               roughness_m="m",stress_eastward_Pa="Pa",stress_northward_Pa="Pa",
                               sea_ice_fraction="fraction",stress_closure_density_kg_m3="kg/m^3"),
                    limitations=["T21 resolves regional circulation; coastal breezes and gusts require finer atmospheric grids.",
                                 "Earth-derived boundary-layer and sea-roughness laws operate with lunar gravity.",
                                 "The climate land mask and atlas coastline have different resolutions.",
                                 "The separate climate programme retains its recorded thermal and humidity uncertainties."])
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **arrays, metadata=np.array(json.dumps(metadata)))
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "control", "export"))
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.action == "export":
        result = export(args.output / "snapshots", args.output / "atmosphere.npz")
    else:
        result = run(args.source, args.output / ("snapshots" if args.action == "run" else "control"),
                     snapshots=args.action == "run")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
