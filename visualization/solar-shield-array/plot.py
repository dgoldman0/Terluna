"""Display the committed holding product; physics stays in protection/research."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"research/studies/solar_shield_array/results/holding.json"
OUT=Path(__file__).with_name("results")


def main():
    product=json.loads(SOURCE.read_text())
    if product["schema"]!="terluna.protection.photogravitational-holding/1":
        raise ValueError("Unsupported holding product")
    OUT.mkdir(exist_ok=True)
    q=product["aperture_quadrature"]
    names=["Fixed\nno sail thrust","Fixed\nideal control","Moving\nideal control","Moving\n10 g/m² film"]
    keys=["baseline_without_solar","baseline","variable_distance_finer","variable_distance_10g"]
    values=[q[k]["mean_power_TW"] for k in keys]
    fig,axes=plt.subplots(1,3,figsize=(16,6.4))
    fig.subplots_adjust(left=.055,right=.985,bottom=.27,top=.82,wspace=.36)
    fig.patch.set_facecolor("#f4f5f7")
    colors=["#6c7582","#4379a6","#183e6c","#bf802a"]
    ax=axes[0];bars=ax.bar(np.arange(4),values,color=colors,width=.68)
    for bar,value in zip(bars,values):ax.text(bar.get_x()+bar.get_width()/2,value+7,f"{value:.1f}",ha="center",fontsize=10)
    ax.set_xticks(np.arange(4),names,fontsize=9);ax.set_ylim(0,350)
    ax.set_ylabel("Mean electrical holding power, TW")
    ax.set_title("Whole aperture · 4 lunar radii",loc="left",fontsize=12)
    ax=axes[1]
    coeff=np.array(product["trajectory_search"]["variable_distance"]["coefficients_km"])[0]
    phase=np.linspace(0,2*np.pi,721);basis=[np.ones_like(phase)]
    for n in range(1,4):basis.extend([np.cos(n*phase),np.sin(n*phase)])
    distance=np.stack(basis,axis=1)@coeff
    ax.plot(np.degrees(phase),distance/1000,color=colors[2],lw=2)
    ax.axhline(78,color=colors[0],ls="--",lw=1)
    ax.set(xlabel="Synodic longitude difference, degrees",ylabel="Sunward distance, thousand km",xlim=(0,360))
    ax.set_xticks([0,90,180,270,360]);ax.set_title("Bounded monthly trajectory",loc="left",fontsize=12)
    ax=axes[2]
    rows=[r for r in product["mass_sensitivity"] if r["closure"]["closed"] and r["sigma_kg_m2"]<=.2]
    x=np.array([r["sigma_kg_m2"]*1000 for r in rows]);y=np.array([r["closure"]["mean_power_TW"] for r in rows])
    ax.loglog(x,y,"o-",color=colors[1],lw=2)
    ax.axvline(25.9,color=colors[3],ls=":",lw=1.5,label="Stored oxide layers: 25.9 g/m²")
    ax.legend(loc="lower right",frameon=False,fontsize=8)
    ax.set(xlabel="Base film + support mass, g/m²",ylabel="Mean electrical holding power, TW")
    ax.set_title("Fixed 78,000 km · centre-force scan",loc="left",fontsize=12)
    for ax in axes:
        ax.set_facecolor("white");ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",alpha=.18);ax.set_axisbelow(True)
    fig.suptitle("Terluna shield holding — DE440 ephemeris study",fontsize=16,ha="left",x=.02,y=.96)
    fig.text(.02,.07,"Ideal optical control is a lower bound. First three bars use 50 g/m² base mass. A 10 g/m² filter is a new-material target.\nPropulsion: 30 km/s exhaust, 70% efficiency, 45° cant, 300 W/kg hardware, seven-day propellant buffer; storage excluded.",fontsize=9,color="#454c56")
    image=OUT/"holding.png";fig.savefig(image,dpi=180,bbox_inches="tight");plt.close(fig)
    (OUT/"provenance.json").write_text(json.dumps({"schema":"terluna.visualization.holding/1","product_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "renderer_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"bar_values_TW":values,
        "figure":"holding.png","evidence":"Displays the committed numerical product; no additional physical model."},indent=2)+"\n")
    print(image)


if __name__=="__main__":main()
