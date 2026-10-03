"""Render the second lunar solar cycle from the climate domain's wave product."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from visualization.waves.coastal import digest

ROOT = Path(__file__).resolve().parents[2]


def render(source=ROOT/"climate/waves/results/cycle.json",
           output=ROOT/"visualization/waves/results/cycle.png"):
    product = json.loads(source.read_text())
    if product["schema"] != "terluna.climate.wave-cycle/1":
        raise ValueError("Unsupported wave-cycle product")
    grid, stats, series, clock = [product[k] for k in ("grid","statistics","series","clock")]
    wet = np.asarray(grid["wet"])
    time = (np.asarray(series["time_hours"])-clock["report_interval_hours"][0])/24
    np.testing.assert_allclose(time[[0,-1]],[0,clock["cycle_days"]],rtol=0,atol=1e-12)
    metre = next(row for row in stats["thresholds"] if row["hs_m"] == 1)
    offshore = np.asarray(series["offshore"])
    fig, axes = plt.subplots(3,2,figsize=(13.5,13.5))
    fig.subplots_adjust(top=.88,bottom=.11,left=.08,right=.93,hspace=.46,wspace=.33)
    fig.suptitle("A full cycle of waves in Smythii–Marginis",x=.08,y=.975,ha="left",fontsize=20,weight="bold")
    fig.text(.08,.939,f"{clock['cycle_days']:.5f} Earth days of wave spin-up, followed by {clock['cycle_days']:.5f} days of analysis",fontsize=12)
    fig.text(.08,.915,f"Global GCM surface stress → SWAN · 1° basin grid · {product['case']['step_s']} s timestep · one reporting cycle",fontsize=11)
    map_specs = ((stats["node_max_hs_m"],"Largest Hs during the reporting cycle","Significant wave height (m)",
                  0,np.ceil(stats["max_hs_m"]*2)/2,"viridis"),
                 (np.asarray(metre["node_fraction"])*100,"Time with Hs ≥ 1 m","Share of reporting cycle (%)",0,100,"cividis"))
    for ax,(field,title,label,low,high,palette) in zip(axes[0],map_specs):
        z=np.zeros(wet.shape)
        z[wet]=field
        mesh=ax.pcolormesh(grid["longitude_deg"],grid["latitude_deg"],np.ma.array(z,mask=~wet),
                           shading="nearest",cmap=palette,vmin=low,vmax=high)
        np.testing.assert_array_equal(np.ma.compressed(mesh.get_array()),field)
        ax.set_facecolor("#eeeee9")
        ax.set_aspect(1/np.cos(np.radians(np.mean(grid["latitude_deg"]))))
        ax.plot(offshore[0,1],offshore[0,2],"*",ms=11,mec="#223444",color="white",label="Reference point")
        ax.set(xlabel="Longitude (°E)",ylabel="Latitude (°N)",title=title)
        fig.colorbar(mesh,ax=ax,pad=.025,label=label)
    axes[0,0].legend(loc="upper left",fontsize=8.5)
    ax=axes[1,0]
    for field,label,color in ((series["max_hs_m"],"Largest Hs across the sea","#ab603d"),
                              (series["area_mean_hs_m"],"Basin mean Hs","#668050"),
                              (offshore[:,3],"Hs at the reference point","#226f91")):
        line,=ax.plot(time,field,label=label,color=color,lw=1.5)
        np.testing.assert_array_equal(line.get_ydata(),field)
    ax.set(xlabel="Earth days into the second cycle",ylabel="Significant wave height (m)",title="Heights through the full cycle",xlim=(0,clock["cycle_days"]))
    ax.legend(fontsize=8.5)
    ax=axes[2,0]
    period=np.where(offshore[:,3]>=.1,offshore[:,5],np.nan)
    line,=ax.plot(time,period,color="#276f91",lw=1.5)
    np.testing.assert_array_equal(line.get_ydata(),period)
    ax.set(xlabel="Earth days into the second cycle",ylabel="Mean period Tm01 (s)",
           title="Wave period at the reference point · Hs ≥ 0.1 m",xlim=(0,clock["cycle_days"]))
    ax=axes[2,1]
    forcing=product["forcing"]
    forcing_days=np.asarray(forcing["monthly_time_days"])-clock["cycle_days"]
    line,=ax.plot(forcing_days,forcing["monthly_basin_mean_stress_Pa"],color="#ab603d",lw=1.5)
    np.testing.assert_array_equal(line.get_ydata(),forcing["monthly_basin_mean_stress_Pa"])
    ax.set(xlabel="Earth days into the second cycle",ylabel="Basin-mean surface stress (Pa)",
           title="Atmospheric momentum supplied to the sea",xlim=(0,clock["cycle_days"]))
    ax=axes[1,1]
    for threshold,color in zip(stats["thresholds"],["#778c48","#276f91","#ac603b","#76558e"]):
        field=np.asarray(threshold["sampled_area_fraction"])*100
        line,=ax.plot(time,field,label=f"Hs ≥ {threshold['hs_m']:g} m",color=color,lw=1.5)
        np.testing.assert_array_equal(line.get_ydata(),field)
    ax.set(xlabel="Earth days into the second cycle",ylabel="Resolved basin area (%)",title="How widely the sea is exposed",xlim=(0,clock["cycle_days"]),ylim=(0,100))
    ax.legend(fontsize=8.5)
    for ax in axes.flat:
        ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.text(.08,.059,"Occurrence applies to this simulated cycle. Lunar wave physics, weak-wave accuracy and finer coastal geography remain open.",fontsize=10)
    fig.text(.08,.035,"The first cycle is excluded from every map and statistic. Threshold durations interpolate hourly heights and exact cycle boundaries.",fontsize=9.5)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=160,facecolor="white")
    plt.close(fig)
    metadata=dict(schema="terluna.visualization.wave-cycle/1",producer_sha256=digest(Path(__file__)),
                  input_sha256=digest(source),output_sha256=digest(output),
                  evidence=product["evidence"],reading_rule=product["reading_rule"],
                  matplotlib_version=matplotlib.__version__,plotted_arrays_checked_exact=True,
                  clock=clock,displayed_case=product["case"]["name"])
    output.with_suffix(".json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(output)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,default=ROOT/"climate/waves/results/cycle.json")
    parser.add_argument("--output",type=Path,default=ROOT/"visualization/waves/results/cycle.png")
    args=parser.parse_args()
    render(args.source,args.output)
