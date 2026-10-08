"""Verify and export compact numerical snapshots from ignored run files."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from .run import ROOT, HERE
from shared.provenance import constants_changed


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def publish(run_dir):
    holding_path=run_dir/"holding.json";validation_path=run_dir/"validation.json"
    holding=json.loads(holding_path.read_text());validation=json.loads(validation_path.read_text())
    if validation["holding_product_sha256"]!=digest(holding_path):
        raise ValueError("Validation refers to a different holding run")
    for relative,expected in holding["producer"]["source_hashes"].items():
        if digest(ROOT/relative)!=expected:
            raise ValueError(f"Study source changed after run: {relative}")
    if constants_changed(holding["producer"]["constants"]):
        raise ValueError("A shared constant used in the holding run has changed")
    if digest(Path(__file__).with_name("validate.py"))!=validation["producer"]["runner_sha256"]:
        raise ValueError("Validation runner changed after run")
    q=holding["aperture_quadrature"]
    error=abs(q["variable_distance"]["mean_power_TW"]/q["variable_distance_finer"]["mean_power_TW"]-1)
    if error>0.001:
        raise ValueError("Aperture quadrature differs by over 0.1%")
    for label in ("baseline","variable_distance"):
        fine=validation["cadence"][label]["step_600_s_closure"]["mean_power_TW"]
        coarse=(holding["baseline"]["ideal"] if label=="baseline" else holding["trajectory_search"][label]["fine_closure"])["mean_power_TW"]
        if abs(fine/coarse-1)>0.001:
            raise ValueError(f"Time sampling differs by over 0.1%: {label}")
    if validation["propagation"]["feedback"]["max_position_error_m"]>50:
        raise ValueError("Reference controller exceeds the formation margin")
    output=HERE/"results";output.mkdir(exist_ok=True)
    holding["units"]={"times":"seconds in propagation; TDB calendar days in eclipse records",
                      "distance":"km unless _m suffix", "acceleration":"mm/s^2 unless _m_s2 suffix",
                      "power":"TW unless _W suffix", "mass":"kg", "areal_mass":"kg/m^2",
                      "propellant_flow":"kg/s", "coefficients":"km in Sun-following basis, as defined in model.py"}
    holding["reading_rule"]="Use aperture_quadrature for final array comparisons; baseline and fixed_distance_scan use the centre-force surrogate. ideal is a lower bound from relaxed optical control. variable_distance_10g uses the 50 g/m2 trajectory with an unbuilt 10 g/m2 replacement, without re-optimising it."
    holding["producer"].update({"snapshot_exporter":"research/studies/solar_shield_array/publish.py",
                               "exporter_sha256":digest(Path(__file__)),"raw_run_sha256":digest(holding_path)})
    final=output/"holding.json";final.write_text(json.dumps(holding,indent=2,allow_nan=False)+"\n")
    validation["raw_holding_product_sha256"]=validation["holding_product_sha256"]
    validation["holding_product_sha256"]=digest(final)
    validation["producer"].update({"snapshot_exporter_sha256":digest(Path(__file__)),"raw_validation_sha256":digest(validation_path)})
    (output/"validation.json").write_text(json.dumps(validation,indent=2,allow_nan=False)+"\n")
    return {"aperture_quadrature_relative_difference":error,"files":[str(final),str(output/"validation.json")]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir",type=Path,default=ROOT/"research/runs/solar_shield_array")
    print(json.dumps(publish(parser.parse_args().run_dir),indent=2))
