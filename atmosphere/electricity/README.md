# Electricity

The physics behind the [atmospheric-electricity study](../../research/studies/atmospheric_electricity/README.md):
how storms separate charge, how much the air conducts, and how their thunder
carries.

| File | What it does |
|---|---|
| [charging.py](charging.py) | Graupel's bouncing collisions with ice crystals or snow integrated over both exponential size distributions, and the charge each bounce separates under four laboratory laws: Saunders and Peck (1998) and Takahashi (1978) as WRF-ELEC implements them (MicroTed/wrf4-elec at commit e43041b, public domain under the WRF notice), and the charges measured near lunar impact speeds by Pradeep Kumar et al. (2024) and Ávila et al. (2013) |
| [muons.py](muons.py) | Cosmic-ray muons through a deep air column with MCEq, and the ionization they and their decay electrons make; for Earth's standard atmosphere and the Open Moon's column |
| [conductivity.py](conductivity.py) | Ion production from the CRAC:CRII tables and the MCEq muons; small ions' recombination, mobility and attachment to aerosol and droplets, the aerosol swollen by the air's humidity; the conductivity and charge relaxation time of the Open Moon's column in clear air and in cloud |
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

**Humidity** (2026-10-05). The lunar air is humid: in the equatorial box's clear
air the relative humidity is 84 % at the ground, 64 % at 10 km and 52–60 %
from 20 to 50 km (`climate/results/crm/mixed_phase_box_0e.json`). Aerosol takes
up water: by κ-Köhler theory with its Kelvin term (Petters and Kreidenweis 2007)
and κ 0.3, a continental mean, particles of 0.05 µm dry radius grow 1.35-fold
at the ground and 1.10-fold at 34 km, with κ 0.1 and 1.0 (organic-rich to
marine) given as a bracket, since nothing yet sets the Open Moon's aerosol.
Ions attach to the grown particles by Hõrrak et al.'s (2008) coefficient,
scaled with the ions' mobility, so that it grows as the pressure falls, where
the coefficient used before held its sea-level value at every height. With 100
particles per cm³ the clear air then conducts 3.7 × 10⁻¹⁵ S/m at the ground and
7.4 × 10⁻¹⁵ S/m at 34 km, against 4.2 and 9.4 × 10⁻¹⁵ before: the humidity
takes 12 % at the ground (5–24 % across the κ bracket) and 5 % aloft, the
pressure scaling 18–30 % from 30 to 50 km. Charge in that air relaxes in 40
minutes at the ground and 20 at 34 km. The scaling with mobility is the
continuum regime's: it holds while the particles are much larger than the
ions' mean free path, and as the air thins the attachment grows more slowly
toward its kinetic limit. With Fuchs and Sutugin's transition factor relative
to the coefficient's sea-level calibration, for an ion mean free path of 15 nm
(Tammet et al.'s transition length) to 50 nm (derived from the ions' mobility
and mass) at sea level, the clear air conducts 5–12 % more at 34 km, 12–25 %
more at 50 km, 23–45 % more at 70 km and 64–71 % more at the column's top
(the `_transition_` keys), and 1–3 % less at the ground. The column's values
aloft are therefore the low end, and they stand until stage 3 rebuilds the
column with the regional aerosol (the author's decision, 2026-10-06). The small ions' mobility and
recombination take no humidity term: the reference mobilities were measured in
humid boreal air, at its 80–85 % relative humidity on average (Hõrrak 2001),
field and laboratory results disagree on any further effect at 1–3 mol % of
water, and recombination at 53–85 % relative humidity lies within its
measurement uncertainty of the value used (Franchin et al. 2015). Without
aerosol, and in cloud, the column is as it was. The literature behind these
choices is in `humid_conductivity.md` on the data drive. The column's aerosol is
one assumed state; the [aerosol study](../../research/studies/open_moon_aerosol/README.md)
builds the particles each region would hold from the project's designs and puts
the ground's conductivity at 2.2–2.7 × 10⁻¹⁵ S/m over the Moon (0.6–4.7 across its
cases, with Earth checks that lean high by up to two to three times), and at 2–4 % of
that in fog.

MCEq runs in its own environment on the data drive, outside the repository:

```sh
python3 -m venv /media/projectspace/terluna-research/venvs/mceq
/media/projectspace/terluna-research/venvs/mceq/bin/pip install MCEq==1.4.2   # downloads its 277 MB database on first use
python -m atmosphere.electricity.fetch_inputs --download
/media/projectspace/terluna-research/venvs/mceq/bin/python -m atmosphere.electricity.muons earth   # 4 minutes
/media/projectspace/terluna-research/venvs/mceq/bin/python -m atmosphere.electricity.muons moon    # 16 minutes
climate/gcm/.venv/bin/python -m atmosphere.electricity.conductivity
```
