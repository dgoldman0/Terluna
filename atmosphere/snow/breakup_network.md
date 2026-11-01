# Junction loads in a falling snowflake: a first contact-network calculation

The force on the most loaded junction of a snowflake aggregate falling at terminal speed grows more slowly than the
flake's weight, and the bending moment there grows faster than the weight but more slowly than weight times span. In
model aggregates whose mass grows as D^1.90, the peak force grows as D^1.5–1.6 when the upper crystals sit in the shadow
of the lower ones and as D^1.0 when every crystal takes drag on its own projected area; the peak moment grows as
D^2.5–2.6 and D^2.0 under the same two rules, while m·g·D grows as D^2.9. With the junction strength fixed by Earth
flakes of 4–6 cm at their limit, a junction that fails in bending gives a Moon/Earth breaking-span ratio of 1.95 with
shadowing and 2.25 without for Mitchell's aggregates (b = 2.1), where the load analogy of
[the survival note](turbulence_survival.md) gives 2.35. Earth's 4–6 cm flakes then correspond to lunar flakes of
7.8–11.7 cm with shadowing and 9.0–13.5 cm without, against the analogy's 9.4–14.1 cm. A junction failing under direct
force would give ratios of 2.8–4.3 and lunar spans of 11–26 cm. All of this is a calculation on idealised aggregates
with geometric drag rules; no part of it is a measured breakup law.

For a narrow ice neck, bending sets the stress. At the calibration flakes the largest moment is 6–10 mm times the
largest force, so in a circular neck of radius r the bending stress 4M/(πr³) exceeds the direct stress F/(πr²) by
4M/(Fr) ≈ 25–41 mm / r, at least 25 for any neck narrower than a millimetre. The moment criterion is therefore the
result carried forward, with the force criterion as the other bound; the unknown neck radius cancels from the ratio.

The loads stay small because each part of a falling flake carries most of its own weight as drag; a junction carries
only what drag and weight on one side of it leave unbalanced, with the inertial loads of the flake's acceleration.
Without shadowing the imbalance comes from the crystals' random orientations and grows as the square root of their
number, so the normalised loads fall as D^−0.9. Shadowing adds a systematic imbalance, the lowest crystals taking the
drag and the upper ones hanging their weight on them; it grows as the flake fills in (the unshadowed share of the
crystals' area falls from 0.87 at 64 crystals to 0.72 at 8192), and the normalised loads fall only as D^−0.3 to D^−0.4.
Natural flakes are denser in projection than these (area ratio 0.35 against 0.08–0.14); if their moments grew as m·g·D,
the ratio would fall to 6.04^(1/(b+1)) = 1.79–1.86, or 7.1–10.7 cm for Mitchell's aggregates.

| Load growth, assumed or computed | Ratio for b = 2.1 / 2.0 / 1.9 | Lunar span from 4–6 cm, b = 2.1 |
|---|---|---|
| Load analogy: load ∝ m g | 2.35 / 2.46 / 2.58 | 9.4–14.1 cm |
| Dense limit: moment ∝ m g D | 1.79 / 1.82 / 1.86 | 7.1–10.7 cm |
| Moment, opaque shadow | 1.95 / 2.00 / 2.05 | 7.8–11.7 cm |
| Moment, no shadow | 2.25 / 2.34 / 2.44 | 9.0–13.5 cm |
| Force, opaque shadow | 2.84 / 3.03 / 3.26 | 11.4–17.1 cm |
| Force, no shadow | 4.28 / 4.87 / 5.67 | 17.1–25.7 cm |

Settling loads, broadside orientation, fits over 64–8192 crystals; b is the exponent of Mitchell's, Leinonen's or Brown
and Francis's mass law (1.9 is also the model's own). The partial shadow gives moment ratios of 2.04 / 2.10 / 2.17, fits
from 256 crystals shift moment ratios by 0.06 or less, and random orientations give 1.90–2.25.

Turbulence changes these spans by a few per cent. A vertical air velocity varying linearly across the flake by
du = (εD)^(1/3) changes each crystal's drag by twice its local velocity over the fall speed V. With equal dissipation on
both worlds and Mitchell's mass law, du/V at ε = 10⁻² m²/s³ is 0.044 for a 4 cm Earth flake and 0.12 for its lunar
counterpart. From ε = 10⁻⁴ to 10⁻² the lunar moment-criterion span shrinks from 7.7–11.6 to 7.6–11.3 cm with shadowing
and from 8.8–13.1 to 8.1–11.8 cm without; the survival note's whole-flake form (settling load times 1 + 2du/V) gives
7.4–11.1 cm with shadowing. The force-criterion spans without shadowing fall from 16.4–24.3 to 13.9–19.9 cm.

When the junction with the largest moment fails, the flake sheds about a third of its mass. Over flakes of 256 crystals
and more the larger piece keeps a median 0.67 with shadowing and 0.64 without (half the flakes 0.57–0.76, nine in ten at
least 0.52, at most 3% above 0.9), little changed by random orientations (0.63–0.64) or turbulence at du/V = 0.1
(0.65–0.66). The largest-force junction lies farther out (median 0.80–0.84 broadside, 0.72–0.75 random). The cap in
McSnow's sub-cloud references ([sub-cloud test](mcsnow_subcloud.md)) holds aggregates at 14 cm, the analogy's equivalent
of Earth's 6 cm flakes, and splits a flake into equal pieces; under the moment criterion that equivalent is 11.7 cm with
shadowing and 13.5 cm without, and the equal split, 0.5, sits at the bottom edge of the network's distribution, where a
split near two thirds and one third matches it. Only the first failure is followed, and the pieces' own spans are not
computed.

## How the loads were computed

Hierarchical ballistic cluster–cluster aggregation in 3D builds the aggregates from thin circular plates of one diameter
d0 and one mass, each oriented at random: pairs of equal clusters, both turned at random, move straight together with a
random offset until two plates' circumscribing spheres touch, so each merge makes one contact and the contacts form a
spanning tree. Twenty realisations of 4096 plates give 20 aggregates of that size and 40 at each power of two from 8 to
2048; pairing them gives 10 of 8192. The maximum dimension D (largest centre-to-centre distance plus d0) runs from 14 d0
at 64 plates to 185 d0 at 8192, and the mass grows as D^(1.90 ± 0.02) (bootstrap standard deviation, 64–8192 plates).

Each flake falls with its thinnest principal axis vertical (vertical extent 0.41–0.47 of D), with a random orientation
as the sensitivity. Its drag equals its weight and is shared by projected area seen from below on 0.1 d0 pixels: each
pixel goes to the lowest plate over it (opaque), to the plates stacked over it in shares of 0.5^rank (partial), or to
every plate over it (no shadow: each plate's own projected area). The flake's linear and angular acceleration from rest
(d'Alembert) carries the torque that drag and weight leave and the net force turbulence adds, so each junction's loads
balance from either side; the force and moment at a contact, about the midpoint of the two plates' centres, are sums
over one side. The turbulent shear runs along a random horizontal direction and takes the package's Heymsfield–Westbrook
fall speed (turbulence_survival.load) in freezing-level air of 1.08 kg/m³ (Earth) and 0.83 kg/m³ (Moon), with d0 = 1 mm.

| Plates | D / d0 | F/(m g), opaque | F/(m g), no shadow | M/(m g D), opaque | M/(m g D), no shadow | Unshadowed |
|---|---|---|---|---|---|---|
| 64 | 14 | 0.042 | 0.036 | 0.0067 | 0.0067 | 0.87 |
| 256 | 30 | 0.026 | 0.019 | 0.0039 | 0.0032 | 0.81 |
| 1024 | 61 | 0.020 | 0.011 | 0.0032 | 0.0019 | 0.76 |
| 4096 | 125 | 0.015 | 0.0058 | 0.0027 | 0.0012 | 0.72 |
| 8192 | 185 | 0.014 | 0.0041 | 0.0024 | 0.0007 | 0.72 |

Medians, broadside, settling only. Over 64–8192 plates F/(m g) and M/(m g D) fall as D^(−0.38 ± 0.03) and
D^(−0.40 ± 0.04) (opaque), D^(−0.53 ± 0.03) and D^(−0.58 ± 0.04) (partial) and D^(−0.86 ± 0.03) and D^(−0.88 ± 0.04) (no
shadow); from 256 plates the opaque exponents are −0.31 ± 0.05 and −0.31 ± 0.07.

## Assumptions and limits

The aggregates are idealised: equal plates, contacts decided by circumscribing spheres, and one contact per merge, which
makes a tree, where natural flakes have several contacts and interlocked branches that share the moment at a neck. They
grow by joining equal clusters, while natural flakes grow largely by sweeping much smaller crystals and flakes; that
history shapes the contact tree and may move the two-thirds split. They are sparser in projection than natural flakes,
which weakens their shadowing; the dense-limit row brackets that side. The fits cover 14–185 d0: with d0 = 1 mm the
moment-criterion lunar spans (78–146 d0) lie inside, and the force-criterion spans without shadowing (171–257 d0) reach
beyond it. The drag rules are geometric. The flakes fall at Reynolds numbers of 3000–6700 on their span and 40–110 on a
millimetre crystal, where wakes, flow through and around the porous flake and extra drag at its edges shape the loads;
none is resolved, and edge effects lie outside both bounding rules. The orientation is fixed and the loads are
instantaneous, leaving out tumbling, flutter and collisions.

Every junction is equally strong at every flake size; sintering aloft and loosening by sublimation remain unquantified,
as in the survival note. The calibration takes Earth's 4–6 cm CASP II flakes to be breakup-limited; for Mitchell's mass
law it implies junction strengths of 13–25 µN and 0.8–2.4 × 10⁻⁷ N·m with shadowing (8–13 µN and 0.5–1.3 × 10⁻⁷ N·m
without), beside flake weights of 0.50–1.18 mN. Turbulence is a linearised shear at the flake's scale, equally strong on
both worlds.

## What would test it

A resolved flow calculation over scanned natural aggregates, or pressure and velocity measurements on printed replicas
at Reynolds numbers of 10³–10⁴, would fix how drag is shared between lower and upper crystals and so the exponents.
Natural aggregates held at terminal speed in a vertical wind tunnel while they grow would show whether they fail at the
junction of largest moment and keep about two thirds of their mass. Micromechanical tests of necks in natural aggregates
would give the bending strength to set beside the implied moments, testing whether Earth's 4–6 cm flakes are
breakup-limited at all.

## Files and checks

[breakup_network.py](breakup_network.py), run as `OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.breakup_network`,
computes every number here in about 30 CPU-seconds and 190 MB and writes [the product](results/breakup_network.json)
with the inputs, seeds and the sha256 of the module and of the imported turbulence_survival.py and flakes.py; no outside
source was consulted. A self-test on a 128-plate aggregate finds contacts at d0 within 2 × 10⁻¹⁵, no overlaps, loads
balanced within 8 × 10⁻¹⁴, tree sums equal to brute-force sums within 2 × 10⁻¹⁴ and zero loads under drag in proportion
to weight; the largest imbalance over the run is 7 × 10⁻¹⁰ in plate units.
[test_snow_breakup_network.py](../tests/test_snow_breakup_network.py) repeats the self-test and checks that the product
was written by this module and that its ratios follow from its fits. These checks test the arithmetic; the physics rests
on the assumptions above.
