# Sailing and shape

**The Open Moon's mean wind blows east at every height above about 2 km, so no layer carries a floater west with
the Sun; sailing instead moves floaters north and south, 30–70° of latitude in a lunar day, and the high eastward
winds shorten their nights.** The design GCM's layer-mean eastward wind is −0.17 m/s at its lowest level, about
0.7 km above the ground, 0.45 m/s at 4 km, 1.75 m/s at 10 km, 3.5 m/s at 18 km, 6.8 m/s at 29 km and 12 m/s at 42 km.
Holding the Sun needs a westward ground speed of 4.3 cos(latitude) m/s ([navigation](navigation.md)). A wing on a long
tether turns the shear between two layers into a cross-wind force, as a sail does between air and water, and a body
that rides the eastward wind aloft sees the Sun return sooner. Numbers are *screen* values (`sailing` in the
product). A wind product resolved by latitude and solar hour, being built on another branch, will replace the single
global-mean profile and the five CM1 rings used here.

The GCM's sigma levels are placed by the wind product's pressure profile from 2.5 km up, where it agrees with
ρgΔz to 1%; its 0 km pressure is not hydrostatic with the levels above, so the lowest level is placed by extending
the 2.5–5 km log-pressure profile downward, to 0.74 km.

## Riding the eastward wind shortens the night

**A floater riding the mean wind at 10 km sees a 21-day solar cycle at the equator and a 16-day one at 60°; at 20 km,
15 and 10 days.** Moving east at u against the Sun's westward pace v speeds the passage of local time by (v + u)/v.
Its longest geometric night falls from the fixed ground's 14.8 days to 10.5 days at 10 km and 7.6 days at 20 km at
the equator, and to 8.2 and 5.2 days at 60°. At 40 km it would be 4.1 days, in air colder than the 25.7 km freezing
height. Riding the CM1 rings' wind round the full cycle, counting the hours where it runs slow, takes 21.5 days at
10 km and 16.8 at 20 km at the equator, 21.3 days at 45° N and 18.7 days at 80° N for 10 km.

| Latitude | 5 km | 10 km | 20 km | 30 km |
|---|---:|---:|---:|---:|
| 0° | 25.6 days | 21.0 | 15.3 | 11.1 |
| 30° | 25.1 | 20.1 | 14.2 | 10.1 |
| 45° | 24.2 | 18.8 | 12.7 | 8.8 |
| 60° | 22.6 | 16.3 | 10.3 | 6.8 |

The latitude rows use the global-mean wind at every latitude, the simplification the coming product removes.

## The evening return flow nearly holds the Sun

**At 1 km the low evening flow runs west at almost the Sun's pace, so floaters linger near sunset for three to six
weeks in the rings, and for one to two weeks if the flow is a little weaker.** In the CM1 rings the air near the
ground flows toward the afternoon from both sides. Its westward branch, back from the evening and night, reaches 4.16
m/s at 1 km at the equator against the Sun's 4.28 m/s, at 85° past noon with the Sun 5° up: a body there drifts only
0.12 m/s behind the Sun and spends 28 days within 10° of that hour, and a full cycle at that height takes 75 days. At
45° N the branch reaches 2.95 m/s against 3.03 m/s, also near 85° past noon, and holds a body 39 days within 10° of
it, the Sun 3.5° up. At 45° S the slowest drift is 0.22 m/s at 65° past noon, the Sun 17° up (21 days). At 80° S the
flow exceeds the Sun's pace and holds bodies at 51° past noon, the Sun 6° up, drawing them in with an e-folding time
of 8.7 days.

**The margins are small, so the lingering is a possibility to test.** The equatorial body drifts at 3% of the wind's
speed relative to the Sun. With the westward branch 5% weaker the equatorial dwell falls from 28 to 15 days, at 10%
to 10 days and at 20% to 6.5; at 45° N from 39 to 18, 12 and 7 days; at 45° S from 21 to 13, 10 and 6 days. The 80° S
hold survives a 5% weaker branch and becomes a 34-day dwell at 10%. The rings are two-dimensional: air cannot
converge on their storms from the sides, which likely strengthens the low westward branch. The coming wind product
will settle it.

**The lingering hour is beautiful and stormy.** The equatorial sunset band at 1 km sits at the evening end of the
afternoon storm sector, where land rain still falls at 0.4–0.55 mm/h: a low floater there lives in long sunset
light, close to the forests, under the last storms of the day. At 80° S the held hour is in the dry, calm polar air.

## Sailing between two layers

**A wing hung below a body sails it across the shear at 0.4–0.8 m/s: 30–70° of latitude in a lunar day.** The body
sits in one wind, the wing in another; the wing's lift across its relative wind and the body's drag balance. In the
global-mean profile the shear is 1.9 m/s between 10 km and the lowest level and 1.1 m/s over 5–10 km, about 0.21 m/s
per kilometre below 20 km and 0.37 m/s per kilometre between 20 and 40 km. With a wing of a tenth of the body's
projected area (lift coefficient 1, drag 0.1, design choices) the steady balance gives:

| Wing hangs at | Shear | Across the shear | Along the shear | Latitude per lunar day | Tether and ballast, 28 MPa | 92 MPa |
|---|---:|---:|---:|---:|---:|---:|
| 0.74 km, from a body at 10 km | 1.9 m/s | 0.58–0.80 m/s | 0.42–1.30 m/s | 49–67° | 4 mm–39 cm, 0.20–0.65 kg/m² | 2 mm–15 cm, 0.10–0.29 kg/m² |
| 5 km, from a body at 10 km | 1.1 m/s | 0.37–0.46 m/s | 0.23–0.73 m/s | 31–39° | 2 mm–16 cm, 0.04–0.10 kg/m² | 1 mm–8 cm, 0.03–0.08 kg/m² |

The ranges span a sphere's drag (C_d 0.47) and a streamlined body's (0.05), and bodies of 40 m to 2 km, which sail
alike because body drag, wing area and tether force all scale with area. A wing of 30% of the projected area sails
0.34–0.86 m/s across. **A kilometres-long tether must carry its own weight.** Fibre at 28 MPa and 1,500 kg/m³ can hang
at most 11.6 km under lunar gravity, and a weightless wing on it would trail only a few kilometres below the body.
Each tether is therefore sized with a weighted wing so that the wing hangs 90% of the tether's length below the body,
with the top stress at the allowable: the 10.3 km tether to the lowest level needs a wing weighted at about half the
body's pull, and the 5.6 km tether one weighted at 1.4 times it, giving the hanging masses in the table, which the
body carries. 300 MPa fibre at a century's 92 MPa halves the 10-km hanging mass. The tether's drag is included. The
along-shear speed slows the body's eastward drift by up to the shear itself, which lengthens its solar cycle; it
cannot carry the body west of the westernmost layer's wind.

**The same physics drives the StratoSail concept.** A wing hung on a tether up to 15 km below a stratospheric balloon
uses the wind difference between the two heights; the balloon acts as the keel and the wing as the sail, giving the
drift a bias of a few metres per second (Aaron, Heun & Nock 2002). Earth's mean zonal wind differences between 35 and
20 km are 11–26 m/s in their examples; the Open Moon's lower air is gentler, and between 20 and 40 km its 7.4 m/s
shear would give cross speeds about four times those at 10 km (*derived*, cross speed in proportion to the shear).

## Sailors between air and water

**The Portuguese man-of-war sails up to 50–55° from downwind with its float as the sail and its tentacles and polyps
as the sea anchor.** Its sailing speed is a few tens of centimetres a second, at most about 0.4 m/s, and the colonies
come in mirror forms that sail on opposite tacks (Iosilevskii & Weihs 2009, with Totton & Mackie's observed 40–45°).
The by-the-wind sailor *Velella* carries a low, flat sail that works mostly by drag at attack angles of 28–87°, its
hull broadside to the water's flow; selection favoured stability over windward performance (Francis 1991).

**A floater near the sea can sail the same way.** A low body dragging a drogue in the water of a thousandth of its
projected area sails 48° off downwind at 0.6 of the wind's speed; a drogue of 1–5% of its area holds it 70–78° off
downwind at 0.11–0.24 of the wind (the drogue's lift coefficient 1 and drag 0.1 are design choices). This suits
coastal and sea floaters; a high floater's tether would have to reach 10–20 km to touch water.

## Shape

**Streamlining cuts drag 9.4 times (C_d 0.47 to 0.05).** Powered motion costs fall in proportion
([navigation](navigation.md)); in a sailing balance it lets the wing pull the body further along the shear (1.3 m/s
against 0.42 m/s at 10 km with a tenth-area wing) while the cross-shear speed stays near 0.6–0.8 m/s. Flat colonies
and quilts are broadside to vertical gusts and edge-on to horizontal ones; tilted, a quilt is itself a wing.

## What the coming wind product will settle

- North–south winds: the rings and the global-mean profile carry none, so the mean meridional drift is open.
- The westward branch's strength and hour at every latitude and height, which decides where sunset lingering and
  Sun-held hours exist.
- Shear by latitude and hour, which sets the sailing speed along real paths.

## Sources and checks

[sailing.py](sailing.py) holds the functions; [sailing_sources.json](sailing_sources.json) records Aaron, Heun & Nock
(2002), Iosilevskii & Weihs (2009), Francis (1991) and Munro et al. (2019) with what was read. Inputs are the
[design climatology](../../../climate/gcm/products/climatology_A28_dim5_moon.npz) (`ua_layer_mean`), the
[GCM wind product](../../../climate/results/gcm/global_winds_A28_dim5_moon.json) for the sigma-to-height mapping, and
the CM1 rings at the equator, 45° N and S and 80° N and S (the GCM's vertical wind imposed).
[test_sailing.py](test_sailing.py) checks the relative period's limits, the circuit period against the steady period,
a stable and an unstable hour attractor on a constructed wind, the dwell time, the sailing balance's residual, its
drag-only and equal-area limits and its tack symmetry, the hanging tether's depth and stress, and the extrapolated
lowest level. An independent read-only review on 9 October 2026 found the tethers sized without their own weight,
the ring periods taken at the hour-mean wind, the lowest sigma level placed by a non-hydrostatic surface pressure and
the sunset dwell given without its sensitivity; this version corrects all four, and the sailing speeds hold.
