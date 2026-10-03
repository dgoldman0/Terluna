"""Display the computed shore exposure and its fixed-station refinement."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
LABELS=dict(west_face="W",east_face="E",southern_shore="S",offshore="O")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_case(summary):
    path=ROOT/summary["path"]
    if digest(path)!=summary["sha256"]:
        raise ValueError("Coastal case differs from the analysed product")
    case=json.loads(path.read_text())
    if case["schema"]!="terluna.climate.shore-case/1":
        raise ValueError("Unsupported coastal case")
    return case,path


def maps(report,output):
    report_path=ROOT/"climate/waves/results/shore.json"
    candidates=[c for c in report["cases"] if c["name"]==f"h852_s{c['stride']}_d36"]
    summary=min(candidates,key=lambda c:c["spacing_m"])
    case,path=load_case(summary)
    terrain_path=ROOT/case["terrain_path"]
    if digest(terrain_path)!=case["terrain_sha256"]:
        raise ValueError("Terrain differs from the wave calculation")
    with np.load(terrain_path,allow_pickle=False) as source:
        terrain={key:source[key].copy() for key in ("longitude_deg","latitude_deg","depth_m","wet")}
    data,wet=np.asarray(case["values"]),np.asarray(case["wet"])
    hs=data[...,2]
    power=np.linalg.norm(data[...,8:10],axis=-1)*case["energy_flux_conversion"]
    fig,axes=plt.subplots(1,3,figsize=(15,6.5))
    fig.subplots_adjust(left=.065,right=.975,bottom=.27,top=.8,wspace=.24)
    fig.suptitle("Where eastern Smythii's waves reach the coast",x=.065,y=.965,
                 ha="left",fontsize=21,weight="bold")
    fig.text(.065,.908,f"Second solar cycle, day {case['day_in_second_cycle']:.2f} · "
             f"{case['spacing_m']:.0f} m wave grid · 118 m flooded-rock terrain",fontsize=12)
    fig.text(.065,.865,"Steady coastal response to the sampled directional sea and winds",fontsize=11)
    fields=[(terrain["longitude_deg"],terrain["latitude_deg"],terrain["depth_m"],terrain["wet"],"Water depth (m)","Blues",1500),
            (case["longitude_deg"],case["latitude_deg"],hs,wet,"Significant wave height (m)","viridis",float(np.ceil(hs[wet].max()*5)/5)),
            (case["longitude_deg"],case["latitude_deg"],power,wet,"Net energy transport (W/m)","magma",float(np.ceil(power[wet].max()/100)*100))]
    for ax,(x,y,values,mask,title,cmap,upper) in zip(axes,fields):
        mesh=ax.pcolormesh(x,y,np.ma.array(values,mask=~mask),shading="nearest",cmap=cmap,vmin=0,vmax=upper)
        np.testing.assert_array_equal(np.ma.compressed(mesh.get_array()),values[mask])
        ax.contour(terrain["longitude_deg"],terrain["latitude_deg"],terrain["depth_m"],levels=[0],colors="#57554f",linewidths=.6)
        ax.set_facecolor("#dfdbd2")
        ax.set(xlabel="Longitude (°E)",ylabel="Latitude (°N)",title=title)
        ax.set_aspect(1/np.cos(np.deg2rad(np.mean(case["latitude_deg"]))))
        ax.set_xlim(case["longitude_deg"][0],case["longitude_deg"][-1])
        ax.set_ylim(case["latitude_deg"][0],case["latitude_deg"][-1])
        fig.colorbar(mesh,ax=ax,pad=.025,shrink=.8)
        for name,site in summary["sites"].items():
            ax.plot(site["longitude_deg"],site["latitude_deg"],"o",ms=4,mec="black",mfc="white",mew=.6)
            ax.annotate(LABELS[name],(site["longitude_deg"],site["latitude_deg"]),xytext=(5,5),
                        textcoords="offset points",fontsize=9,weight="bold",color="black",
                        bbox=dict(facecolor="white",alpha=.8,edgecolor="none",pad=1))
    step=max(1,len(wet)//13)
    xx,yy=np.meshgrid(case["longitude_deg"],case["latitude_deg"])
    take=np.zeros_like(wet);take[::step,::step]=True;take &= wet & (hs>=.1)
    vector=data[...,8:10]*case["energy_flux_conversion"]
    axes[2].quiver(xx[take],yy[take],vector[take,0]/np.cos(np.deg2rad(yy[take])),vector[take,1],color="white",alpha=.8,
                   angles="xy",scale_units="xy",scale=max(power[wet].max(),1)/.055,width=.004)
    fig.text(.065,.16,"W: western face    E: eastern face    S: southern shore    O: offshore reference",fontsize=11)
    fig.text(.065,.112,"Arrows show the computed direction and relative size of net energy transport.\n"
             "The shoreward station values integrate incoming directions separately. Grid and boundary checks accompany the study.",fontsize=10)
    fig.text(.065,.052,"Static flooded rock at the atlas water level; beach formation, reflection, diffraction and run-up remain separate calculations.",fontsize=9,color="#444444")
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=170,facecolor="white")
    plt.close(fig)
    return [report_path,path,terrain_path],dict(case=case["name"],wave_nodes=int(wet.sum()),
        native_terrain_shape=list(terrain["depth_m"].shape),displayed_hs_max_m=float(hs[wet].max()),
        displayed_transport_max_w_m=float(power[wet].max()),exact_array_checks=3)


def phase_maps(report,output):
    names=[f"h{hour}_s2_d36" for hour in (852,1020,1068)]
    summaries={c["name"]:c for c in report["cases"]}
    if any(name not in summaries for name in names):
        return None
    loaded=[load_case(summaries[name]) for name in names]
    terrain_path=ROOT/loaded[0][0]["terrain_path"]
    if any(c["terrain_sha256"]!=digest(terrain_path) for c,p in loaded):
        raise ValueError("Weather phases must share the same terrain")
    with np.load(terrain_path,allow_pickle=False) as terrain:
        tx,ty,depth=terrain["longitude_deg"],terrain["latitude_deg"],terrain["depth_m"]
    powers=[]
    for c,path in loaded:
        data=np.asarray(c["values"])
        powers.append(np.linalg.norm(data[...,8:10],axis=-1)*c["energy_flux_conversion"])
    maximum=max(float(power[np.asarray(c["wet"])].max()) for power,(c,path) in zip(powers,loaded))
    upper=np.ceil(maximum/100)*100
    fig,axes=plt.subplots(1,3,figsize=(14,6.5))
    fig.subplots_adjust(left=.065,right=.88,bottom=.24,top=.77,wspace=.22)
    fig.suptitle("The exposed coast changes as the sea turns",x=.065,y=.962,ha="left",fontsize=21,weight="bold")
    fig.text(.065,.9,"Three fixed weather phases from the second solar cycle · common 237 m wave grid",fontsize=12)
    headings=("East-northeast sea","North-northeast sea","West-southwest sea")
    for ax,(c,path),power,heading in zip(axes,loaded,powers,headings):
        wet=np.asarray(c["wet"])
        data=np.asarray(c["values"])
        mesh=ax.pcolormesh(c["longitude_deg"],c["latitude_deg"],np.ma.array(power,mask=~wet),
                           shading="nearest",cmap="magma",vmin=0,vmax=upper)
        np.testing.assert_array_equal(np.ma.compressed(mesh.get_array()),power[wet])
        ax.contour(tx,ty,depth,levels=[0],colors="#57554f",linewidths=.6)
        ax.set_facecolor("#dfdbd2")
        ax.set(xlabel="Longitude (°E)",ylabel="Latitude (°N)",title=f"Day {c['day_in_second_cycle']:.2f}\n{heading}")
        ax.set_aspect(1/np.cos(np.deg2rad(np.mean(c["latitude_deg"]))))
        ax.set_xlim(c["longitude_deg"][0],c["longitude_deg"][-1])
        ax.set_ylim(c["latitude_deg"][0],c["latitude_deg"][-1])
        xx,yy=np.meshgrid(c["longitude_deg"],c["latitude_deg"])
        step=max(1,len(wet)//12)
        take=np.zeros_like(wet);take[::step,::step]=True;take &= wet & (data[...,2]>=.1) & (power>=1)
        vector=data[...,8:10]*c["energy_flux_conversion"]
        heading_vectors=vector[take]/np.linalg.norm(vector[take],axis=-1)[:,None]
        ax.quiver(xx[take],yy[take],heading_vectors[:,0]/np.cos(np.deg2rad(yy[take])),heading_vectors[:,1],
                  color="white",alpha=.8,angles="xy",scale_units="xy",scale=1/.035,width=.004)
        for name,site in summaries[c["name"]]["sites"].items():
            ax.plot(site["longitude_deg"],site["latitude_deg"],"o",ms=4,mec="black",mfc="white",mew=.6)
            ax.annotate(LABELS[name],(site["longitude_deg"],site["latitude_deg"]),xytext=(4,5),
                        textcoords="offset points",fontsize=9,weight="bold",
                        bbox=dict(facecolor="white",alpha=.8,edgecolor="none",pad=1))
    fig.colorbar(mesh,cax=fig.add_axes([.915,.30,.014,.43]),label="Net wave-energy transport (W/m)")
    fig.text(.065,.14,"Names describe basin wave travel. Equal-length arrows show local energy-flow direction; colours compare power.\n"
             "W and E mark the two headland faces; S marks the southern shore; O is the offshore reference.",fontsize=10)
    fig.text(.065,.06,"Selected stationary responses on flooded rock. Timing, breaking-zone structure and run-up require further calculations.",fontsize=10,color="#444444")
    fig.savefig(output,dpi=170,facecolor="white")
    plt.close(fig)
    inputs=[ROOT/"climate/waves/results/shore.json",terrain_path]+[p for c,p in loaded]
    return inputs,dict(cases=names,common_spacing_m=loaded[0][0]["spacing_m"],
                        colour_scale_w_m=[0,float(upper)],arrows="Equal length; local net energy-flow direction where Hs>=0.1 m and power>=1 W/m",exact_array_checks=3)


def save_record(output,inputs,display):
    record=dict(schema="terluna.visualization.shore-exposure/1",
        evidence="Exact array rendering of computed wave height, native flooded terrain and resolved net energy transport.",
        reading_rule="Maps show selected fixed weather phases from the second solar cycle. The domain report supplies numerical comparisons and the boundary interpolation assumption.",
        producer=dict(path="visualization/waves/shore.py",sha256=digest(Path(__file__))),
        inputs=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p)) for p in inputs],
        output=dict(path=output.name,sha256=digest(output)),display=display,
        matplotlib_version=matplotlib.__version__)
    output.with_suffix(".json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(display))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=HERE/"results/shore.png")
    args=parser.parse_args()
    report=json.loads((ROOT/"climate/waves/results/shore.json").read_text())
    if report["schema"]!="terluna.climate.shore-exposure/1":
        raise ValueError("Unsupported coastal exposure report")
    save_record(args.output,*maps(report,args.output))
    phase_output=args.output.with_name(args.output.stem+"_phases"+args.output.suffix)
    result=phase_maps(report,phase_output)
    if result:
        save_record(phase_output,*result)


if __name__=="__main__":
    main()
