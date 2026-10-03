"""Display the computed sea response to a recovered atmospheric weather sequence."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from visualization.waves.coastal import digest

ROOT = Path(__file__).resolve().parents[2]


def render(source=ROOT / "climate/waves/results/weather.json", case_name=None,
           output=ROOT / "visualization/waves/results/weather.png"):
    product = json.loads(source.read_text())
    if product["schema"] != "terluna.climate.weather-wave-study/1":
        raise ValueError("Unsupported weather-wave product")
    case_name = case_name or product["displayed_case"]
    case = next(c for c in product["cases"] if c["name"] == case_name)
    grid, forcing = product["basin_grid"], product["forcing"]
    wet = np.asarray(grid["wet"])
    values = np.zeros(wet.shape)
    values[wet] = case["node_max_hs_m"]
    offshore = np.asarray(case["offshore_records"])
    time = np.asarray(case["hourly_time_hours"])/24
    fig, axes = plt.subplots(2,2,figsize=(13,10))
    fig.subplots_adjust(top=.85,bottom=.16,left=.08,right=.93,hspace=.43,wspace=.35)
    fig.suptitle("Waves through a changing lunar weather sequence",x=.08,y=.975,ha="left",fontsize=20,weight="bold")
    fig.text(.08,.926,"Smythii–Marginis · 7 Earth days · recovered GCM surface stress → SWAN",fontsize=12)
    fig.text(.08,.893,f"1° basin grid · {case['step_s']} s timestep · {case['displayed_from_run_hour']:g} h preceding weather · lunar gravity",fontsize=11)
    ax = axes[0,0]
    vmax = np.ceil(values.max()*2)/2
    mesh = ax.pcolormesh(grid["longitude_deg"],grid["latitude_deg"],np.ma.array(values,mask=~wet),
                         shading="nearest",cmap="viridis",vmin=0,vmax=vmax)
    np.testing.assert_array_equal(np.ma.compressed(mesh.get_array()),case["node_max_hs_m"])
    ax.set_facecolor("#eeeeeb")
    ax.set_aspect(1/np.cos(np.radians(np.mean(grid["latitude_deg"]))))
    ax.plot(offshore[0,1],offshore[0,2],"*",ms=11,color="white",mec="#263443",label="Reference point")
    ax.set(xlabel="Longitude (°E)",ylabel="Latitude (°N)",title="Largest Hs at each node, days 1–7")
    ax.legend(loc="upper left",fontsize=9)
    fig.colorbar(mesh,ax=ax,pad=.025,label="Significant wave height (m)")
    ax = axes[0,1]
    line, = ax.plot(forcing["monthly_time_days"],forcing["monthly_basin_mean_stress_Pa"],color="#446f83",lw=1.6)
    np.testing.assert_array_equal(line.get_ydata(),forcing["monthly_basin_mean_stress_Pa"])
    ax.axvspan(forcing["start_day_from_first_snapshot"],forcing["end_day_from_first_snapshot"],color="#d59b45",alpha=.2,label="Simulated week")
    ax.set(xlabel="Earth days since first atmospheric snapshot",ylabel="Basin-mean surface stress (Pa)",title="The 30-day atmospheric continuation")
    ax.legend(fontsize=9)
    ax = axes[1,0]
    for field,label,color in ((case["hourly_max_hs_m"],"Largest Hs across the sea","#9d522c"),
                               (offshore[:,3],"Hs at the reference point","#236b8e")):
        line, = ax.plot(time,field,label=label,color=color,lw=1.8)
        np.testing.assert_array_equal(line.get_ydata(),field)
    ax.axvspan(0,1,color="grey",alpha=.12,label="First 24 h of selected week")
    ax.set(xlabel="Earth days into simulated week",ylabel="Significant wave height (m)",title="The sea builds and eases through the week",xlim=(0,7))
    ax.legend(fontsize=8.5)
    ax = axes[1,1]
    period = np.where(offshore[:,3] > 0,offshore[:,5],np.nan)
    line, = ax.plot(time,period,color="#236b8e",lw=1.8,label="Mean wave period")
    np.testing.assert_array_equal(line.get_ydata(),period)
    ax.set(xlabel="Earth days into simulated week",ylabel="Mean period Tm01 (s)",title="Period and forcing at the reference point",xlim=(0,7))
    right = ax.twinx()
    equivalent = np.linalg.norm(np.asarray(case["offshore_equivalent_wind_m_s"]),axis=1)
    line2, = right.plot(time,equivalent,color="#9d522c",lw=1.4,ls="--",label="Stress-equivalent input")
    np.testing.assert_array_equal(line2.get_ydata(),equivalent)
    right.set_ylabel("Stress-equivalent input (m/s)")
    ax.legend([line,line2],[line.get_label(),line2.get_label()],loc="lower right",fontsize=8.5)
    for ax in axes.flat:
        ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.text(.08,.082,"Selected strong-forcing episode. Wave occurrence across months and seasons requires a longer sample.",fontsize=10)
    fig.text(.08,.051,"Regional atmospheric forcing and coarse basin geometry; lunar wind input, dissipation and coastal breaking remain research questions.",fontsize=9.5)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=160,facecolor="white")
    plt.close(fig)
    metadata = dict(schema="terluna.visualization.weather-waves/1",evidence=product["evidence"],
                    reading_rule=product["reading_rule"],producer_sha256=digest(Path(__file__)),
                    input_sha256=digest(source),output_sha256=digest(output),
                    matplotlib_version=matplotlib.__version__,displayed_case=case["name"],
                    plotted_arrays_checked_exact=True,map_limits_m=[0,vmax])
    output.with_suffix(".json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(output)


if __name__ == "__main__":
    render()
