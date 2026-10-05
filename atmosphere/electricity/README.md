# Electricity

The physics behind the [atmospheric-electricity study](../../research/studies/atmospheric_electricity/README.md):
how storms separate charge, how much the air conducts, and how their thunder
carries.

| File | What it does |
|---|---|
| [charging.py](charging.py) | Graupel's bouncing collisions with ice crystals or snow integrated over both exponential size distributions, and the charge each bounce separates under four laboratory laws: Saunders and Peck (1998) and Takahashi (1978) as WRF-ELEC implements them (MicroTed/wrf4-elec at commit e43041b, public domain under the WRF notice), and the charges measured near lunar impact speeds by Pradeep Kumar et al. (2024) and Ávila et al. (2013) |
| [muons.py](muons.py) | Cosmic-ray muons through a deep air column with MCEq, and the ionization they and their decay electrons make; for Earth's standard atmosphere and the Open Moon's column |
| [conductivity.py](conductivity.py) | Ion production from the CRAC:CRII tables and the MCEq muons; small ions' recombination, mobility and attachment to aerosol and droplets; the conductivity and charge relaxation time of the Open Moon's column in clear air and in cloud |
| [thunder.py](thunder.py) | Thunder heard at the ground: rays through a layered atmosphere with wind (the effective sound speed), as circular arcs within each layer, turning, reflecting from a hard ground and gathered by range and arrival time; absorption by octave band after ISO 9613-1, scaled to the air's oxygen and nitrogen; Few's (1969) spectral peak; levels against the threshold of hearing (ISO 226) and A-weighted. Checked against ISO 9613's table, the inverse-square law, exact arcs in a linear sound-speed profile and the silent zone's edge |
| [inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py) | Takahashi's charging table, the CRAC:CRII v2 tables and the PDG muon energy-loss table for dry air, fetched and hash-checked, kept out of Git |

**Checks.** The ion chain reproduces Earth's measured fair-weather positive
conductivity (Gringel 1978, as fitted in Nicoll 2012) within 2–8 % from 5 to 30
km at Germany's 3 GV cutoff. MCEq (SIBYLL 2.3e, Hillas–Gaisser H3a) gives Earth's
sea-level vertical muon intensity within 6–13 % of the Bogdanova et al. (2006)
fit above 10 GeV/c; above 1 GeV/c it gives 51 m⁻² s⁻¹ sr⁻¹ against the fit's 76
(measurements 60–70), short where 3-D effects matter. The muons that reach the
lunar troposphere start above 7–20 GeV. At Earth's 1,030 g/cm² its muons and their
decay electrons make 1,050 ion pairs per gram per second, against CRAC:CRII's
1,500 there with the hadronic and electromagnetic cascade included. The charging
laws reproduce WRF-ELEC's critical rime accretion rates and Takahashi's tabulated
charges.

**Results.** The Open Moon's column (`results/conductivity_moon.json`): muons and
their decay electrons make 0.04 ion pairs per cm³ per second at the ground and
0.05–0.06 at the storms' charging zone (30–50 km), where clear air conducts
1.2–1.9 × 10⁻¹⁴ S/m and cloud 1 × 10⁻¹⁶ S/m. The muon products
(`results/muon_ionization_{earth,moon}.json`) give the ionization by depth.

MCEq runs in its own environment on the data drive, outside the repository:

```sh
python3 -m venv /media/projectspace/terluna-research/venvs/mceq
/media/projectspace/terluna-research/venvs/mceq/bin/pip install MCEq==1.4.2   # downloads its 277 MB database on first use
python -m atmosphere.electricity.fetch_inputs --download
/media/projectspace/terluna-research/venvs/mceq/bin/python -m atmosphere.electricity.muons earth   # 4 minutes
/media/projectspace/terluna-research/venvs/mceq/bin/python -m atmosphere.electricity.muons moon    # 16 minutes
climate/gcm/.venv/bin/python -m atmosphere.electricity.conductivity
```
