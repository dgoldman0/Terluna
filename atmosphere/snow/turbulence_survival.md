# How large a lunar snowflake can survive

A lunar snowflake carries no more load at its junctions than an Earth flake of a given span until it is 2.2–2.6 times
larger, when both fly through the same turbulence. Earth's ordinary large flakes of about 2 cm therefore correspond to
lunar flakes of about 4–5 cm, and the 4–6 cm aggregates observed in a heavy, sheared Earth storm to lunar flakes of
about 8–15 cm. The flake's own weight sets almost all of this: at terminal speed its drag equals its weight, which is a
sixth of Earth's for the same ice. Turbulent velocity differences across the flake add between 1% and a third of the
load over the range of turbulence tested, and they shrink the lunar span by up to a fifth only where lunar turbulence is
ten times Earth's. This is a scaling estimate by analogy with observed flakes; no measured strength law transfers to a
decimetre aggregate.

## The loads

A flake falling at terminal speed V carries its weight m g as drag, and that weight passes through the junctions that
hold the flake together. Turbulence adds a fluctuating drag where the air's velocity differs across the flake. For a
velocity increment du, small beside V, the fluctuation is 2 m g du / V, with du = (εD)^(1/3) across a span D in the
inertial range at dissipation rate ε, and the viscous form below the Kolmogorov length. The fall speed follows
Heymsfield and Westbrook's (2010) law in each world's air at its freezing level: 0.83 kg/m³ on the Moon's 301 K ladder
and 1.08 kg/m³ in the Earth winter column. The mass law is McSnow's aggregate law, Mitchell's (1996) aggregates of side
planes, columns and bullets, which the aggregation of plan step 4 uses; Leinonen's dendrite aggregates and Brown and
Francis (1995) give the sensitivity.

Earth's largest observed flakes are taken as surviving at their limit: aggregates of 4–6 cm in the 9 mm/h CASP II
storm, imaged from aircraft and collected on the ground, with breakup inferred to be inhibited near freezing while
aggregation continued (Lawson et al. 1993, 1998), and flakes of a centimetre or two in ordinary snowfall. For each
Earth flake the calculation finds the lunar span with the same peak load, weight plus turbulent fluctuation, for
dissipation rates of 10⁻⁵ to 10⁻² m²/s³ and lunar turbulence from a tenth to ten times Earth's.

| Earth flake | Lunar span with the same load | Settling load alone |
|---|---|---|
| 2 cm, ordinary large flakes | 3.9–5.2 cm | 4.7–5.2 cm |
| 4 cm | 7.7–10.3 cm | 9.4–10.3 cm |
| 6 cm, the largest CASP II aggregates | 11.4–15.5 cm | 14.1–15.5 cm |

The settling-only span is the gravity ratio raised to one over the mass exponent: 6.04^(1/2.1) = 2.36 for Mitchell's
aggregates, 2.46 for Leinonen's and 2.58 for Brown and Francis's. A lunar flake falls about half as fast as an Earth
flake of the same span, so the same turbulence is a larger fraction of its speed: du/V reaches 0.13 on the Moon at
ε = 10⁻² m²/s³, against 0.05 on Earth. At those rates turbulence carries up to 36% of the lunar flake's load if the
Moon's air is ten times as turbulent as Earth's, which brings the 6 cm equivalent down to about 12 cm.

## What the estimate assumes

The weakest junction carries a load that grows with the flake, and its strength does not depend on the flake's size,
so a lunar flake built from similar crystals and contacts fails at the same junction load as an Earth flake. Earth's
largest observed flakes may instead be limited by the time or water available to grow them; the lunar spans would then
be conservative. Lunar flakes spend one to four days aloft, which gives their contacts time to sinter and strengthen
(Demmenie et al. 2025) and also exposes them to repeated sublimation, which loosened contacts in laboratory dendrites
(Oraltay and Hallett 1989); neither effect is quantified here. The lunar dissipation rate is unmeasured: the archived
cloud simulations save no turbulence field ([turbulence](turbulence.md)), so the estimate spans two decades around
Earth's.

## Beside the aggregation

McSnow's aggregation, which reproduces the median size of Earth's snowflakes, grows the lunar column's flakes past
these spans within the first one and a half kilometres below the cloud top at 3 mm/h
([aggregation](mcsnow_aggregation.md)). Where aggregation outruns survival in this way, breakup decides the size that
falls, and these load-equivalent spans indicate where breakup begins. That inference rests on the analogy above and on
McSnow's lack of a whole-aggregate breakup law.

## Reproduction

```sh
python -m atmosphere.snow.turbulence_survival
```

The [runner](turbulence_survival.py) writes the compact [product](results/turbulence_survival.json) with every span,
load and air state.
