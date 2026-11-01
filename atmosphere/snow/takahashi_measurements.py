"""Read the public Takahashi (1991) data copy with only Python's stdlib.

The dataset is Fuchs et al. (2025b), doi:10.5281/zenodo.15804055,
1991_takahashi.xlsx, supplied with Takahashi's permission per Readme.txt.
Rows retain their original Excel row identity. Missing cax values stay None.
"""
from pathlib import Path
import hashlib
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

EXPECTED_SHA256 = "bfb05429b98aa4c5857f07a42f9594ddc8d1c7b9969baa025a1c34e448d72a63"
EXPECTED_MD5 = "ecb41c07eab7982cadeff694c4f45418"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def read_rows(path):
    data = Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != EXPECTED_SHA256:
        raise ValueError("Unidentified laboratory workbook SHA256")
    if hashlib.md5(data).hexdigest() != EXPECTED_MD5:
        raise ValueError("Workbook does not match the published Zenodo MD5")
    with zipfile.ZipFile(path) as z:
        strings = ["".join(n.itertext()) for n in
                   ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)]
        xmlrows = ET.fromstring(z.read("xl/worksheets/sheet1.xml")).findall(".//m:row", NS)
        rows = []
        for row in xmlrows:
            values = {}
            for cell in row.findall("m:c", NS):
                v = cell.find("m:v", NS)
                if v is None:
                    continue
                column = re.match(r"[A-Z]+", cell.attrib["r"]).group()
                values[column] = strings[int(v.text)] if cell.attrib.get("t") == "s" else float(v.text)
            if not values:
                continue
            if row.attrib["r"] == "1":
                if values != {"A": "temp", "B": "time", "C": "aax", "D": "cax", "E": "mass_old", "F": "mass"}:
                    raise ValueError("Unexpected workbook header")
                continue
            for key in ("A", "B", "C", "F"):
                if not math.isfinite(values[key]):raise ValueError("Nonfinite measurement")
            if any(values[key]<=0 for key in ("B", "C", "F")):
                raise ValueError("Nonpositive measured age, width or mass")
            if "D" in values and (not math.isfinite(values["D"]) or values["D"]<=0):
                raise ValueError("Invalid measured polar dimension")
            old_mass = float(values["E"].replace("x10^", "e"))
            if not math.isclose(old_mass, values["F"], rel_tol=1e-12):
                raise ValueError("Disagreement between text and numeric mass")
            result = {"excel_row": int(row.attrib["r"]), "temperature_C": values["A"],
                      "age_min": values["B"], "a_full_mm": values["C"],
                      "c_full_mm": values.get("D"), "mass_g": values["F"]}
            result["mass_kg"] = result["mass_g"] * 1e-3
            result["maximum_dimension_mm"] = max(x for x in (result["a_full_mm"], result["c_full_mm"]) if x is not None)
            result["spheroid_density_kg_m3"] = None
            if result["c_full_mm"] is not None:
                v = math.pi / 6 * (result["a_full_mm"] * 1e-3)**2 * (result["c_full_mm"] * 1e-3)
                result["spheroid_density_kg_m3"] = result["mass_kg"] / v
            flags=[]
            if result['c_full_mm'] is None:flags.append('c_axis_missing_in_source')
            # Material value from the pinned SCALE source, not a fitted lunar constant.
            if result['spheroid_density_kg_m3'] is not None and result['spheroid_density_kg_m3']>916.8:
                flags.append('mass_and_axes_exceed_native_ice_density_under_spheroid_mapping')
            excluded=[]
            if result['excel_row']==72:
                flags.append('source_c_axis_decimal_anomaly_unresolved')
                excluded=['c_full_mm','maximum_dimension_mm','spheroid_density_kg_m3','c_to_a_axis_ratio']
            result['quality_flags']=flags
            result['calibration_excluded_observables']=excluded
            result['maximum_dimension_basis']='Maximum of reported axes; missing polar dimensions are not reconstructed.'
            rows.append(result)
    if len(rows)!=169:raise ValueError("Unexpected measured-row count")
    return rows


def diagnostic_windows(rows):
    """Explicit ±0.5 C windows, NOT author-assigned groups or uncertainty intervals."""
    out = {}
    for target in (-3.7, -5.3, -8.6, -10.6, -14.4, -20.1, -22.0):
        selected = [r for r in rows if abs(r["temperature_C"] - target) <= 0.50000000001]
        groups = {}
        for age in sorted({r["age_min"] for r in selected}):
            subset = [r for r in selected if r["age_min"] == age]
            entry = {"rows": [r["excel_row"] for r in subset], "n": len(subset), "calibration_exclusions": {str(r["excel_row"]):r["calibration_excluded_observables"] for r in subset if r["calibration_excluded_observables"]}}
            for key in ("temperature_C", "a_full_mm", "c_full_mm", "mass_g", "maximum_dimension_mm"):
                v = [r[key] for r in subset if r[key] is not None]
                entry[key + "_range"] = [min(v), max(v)] if v else None
            groups[str(age)] = entry
        out[str(target)] = {"selection": "abs(measured_temperature_C-target_C)<=0.5", "by_age": groups,
                            "three_minute_rows": [r for r in selected if r["age_min"] == 3]}
    return out


if __name__ == "__main__":
    path = Path(sys.argv[1])
    rows = read_rows(path)
    result = {"schema":"terluna.atmosphere.takahashi-measurements/1",
              "evidence":"Measured crystal data transcribed from the identified public workbook; geometric conversions and temperature windows are our diagnostics, not confidence intervals.",
              "attribution":"Takahashi et al. (1991), doi:10.2151/jmsj1965.69.1_15; workbook supplied with Tsuneya Takahashi's permission in Fuchs et al. (2025b), doi:10.5281/zenodo.15804055.",
              "source_license":"CC-BY-4.0, https://creativecommons.org/licenses/by/4.0/; source dataset metadata checked. This product converts units, adds row identifiers, derived spheroid densities and diagnostic groupings.",
              "parser_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "source": "https://zenodo.org/records/15804055", "source_file": path.name,
              "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "source_md5": EXPECTED_MD5,
              "units": {"age": "min", "full_axes": "mm", "mass": "g"},
              "unit_evidence": "Fuchs2025 f08.py lines304-314 converts axes*1e3/(minutes*60) to micrometers/s; mass strings agree with numeric grams and original1991 paper.",
              "geometry_caution": "Derived spheroid density is a representation diagnostic, not measured bulk density. Missing planar thickness remains absent.",
              "data_quality": {
                  "flag_counts":{flag:sum(flag in r['quality_flags'] for r in rows) for flag in sorted({f for r in rows for f in r['quality_flags']})},
                  "row72":{
                      "source_cell":"xl/worksheets/sheet1.xml:D72",
                      "source_value_mm":0.8923,
                      "source_comparison":"The original paper's Figure 7 (printed p.22) places the corresponding polar-axis series near 0.09 mm at 12 minutes. Table 2 predicts 0.08584 mm at nominal -10.6 C; row 72 is at -10.3 C. A decimal anomaly is suspected, not corrected.",
                      "author_plot":"Fuchs f08.py includes this age/temperature selection and has no explicit row exclusion.",
                      "handling":"Preserve the source entry and exclude its polar dimension and dependent size/shape quantities from calibration pending source verification. Raw window ranges retain the flag."},
                  "geometry_caution":"Exceeding bulk ice density after spheroid mapping indicates representation incompatibility, not necessarily bad measurements."},
              "row_count": len(rows), "rows": rows, "diagnostic_windows": diagnostic_windows(rows)}
    print(json.dumps(result, indent=2))
