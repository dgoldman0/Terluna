# Aimed projectiles on the Moon

**8 November 2026.** A stage of the sport study ([../README.md](../README.md)) of its own: how aimed projectiles fly in the Open Moon's air, what shooting asks of a shooter's body there, and whether a paintball can be dodged. A paintball game needs these figures first. The study flies five projectiles (an arrow, an air-rifle pellet, a small-bore bullet, a paintball and a clay target) in the design run's air at the sea-level datum (1.40 kg/m³) and at the shared gravity, beside Earth's air and gravity. It sets the shooter's recoil and sway against the two worlds, and lets a target dodge a seen paintball under the sport study's reach engine. Everything is in [aimed.py](aimed.py), with its product [results/projectiles.json](results/projectiles.json) and its sources in [sources.json](sources.json). The author's directions are in the decisions register, section "Sport from the body, 7 November 2026".

## The flights

The projectiles are points with quadratic drag (the clay is a spinning disc with lift), launched 1.4 m above the floor and aimed at a chest 1.2 m up. Two things change on the Moon. The drop over a given flight is a sixth of Earth's, so the sight barely lifts. The air is 15% denser, so light projectiles lose speed a little faster. Over sporting distances the flight times hardly change, because they are set by the muzzle speed and the drag.

| Projectile | Distance | Flight (Moon, Earth) | Arrives at (Moon, Earth) | Aim above level (Moon, Earth) |
|---|---|---|---|---|
| Paintball, 3.2 g at 91 m/s | 20 m | 0.28 s, 0.27 s | 56, 60 m/s | -0.4°, 0.3° |
| | 30 m | 0.48, 0.46 s | 44, 49 m/s | -0.1°, 1.2° |
| | 50 m | 1.07, 0.97 s | 27, 32 m/s | 0.5°, 3.7° |
| Arrow, 19 g at 60 m/s | 70 m | 1.28, 1.27 s | 50, 51 m/s | 0.9°, 6.0° |
| Air-rifle pellet, 0.53 g at 175 m/s | 10 m | 0.061, 0.061 s | 153, 155 m/s | (drop 3 mm, 17 mm) |
| .22 bullet, 2.6 g at 326 m/s | 50 m | 0.162, 0.160 s | 294, 298 m/s | (drop 2 cm, 12 cm) |

A paintball in the lunar air loses half its speed every 29 m (Earth 33 m) and falls at most 8.2 m/s (Earth 21.5). Shot level, it drops only 0.15 m by 30 m, where Earth's has dropped 0.85 m and reached the floor before 40 m. A crosswind of 1 m/s moves it 0.15 m by 30 m (0.13 on Earth). It arrives at 30 m with 3.1 J (3.8 on Earth) and at 60 m with 0.7 J (1.1), so it stings much as on Earth at matching range.

**A stray shot carries farther, by less than gravity alone suggests.** The best elevation from 1.4 m carries a paintball 125 m (95 on Earth), an arrow 531 m (232; 815 m and 283 m for an aligned arrow with laminar flow), an air-rifle pellet 249 m (200) and a .22 bullet 1.5 km (1.1 km), with up to 40 s in the air. The denser air and drag, which hold every long flight near its terminal speed, keep the lunar ranges to 1.2-2.3 times Earth's (2.9 for the laminar arrow) rather than six. A range's safety zone still has to grow by those factors, and upward by more: a .22's apex reaches 380 m.

**Clay targets glide.** A clay thrown along its own plane carries lift as a spinning disc, at a sixth of the weight in denser air. Started from Andert, Freudenthal and Levedag's coefficients, the model gives a skeet high house's long throw 62.6 m on Earth, inside the measured 61-63 m; the low house's comes out 61 m against a measured 68-70, an 11% shortfall. The same machine on the Moon throws a clay 135 m, rising 20-27 m and staying up 18-20 s. To land at the field's 62-69 m a lunar machine would throw at 12.6-14.1 m/s, and the clay would then hang 10-12 s, three times Earth's 4 s. Clay shooting on the Moon is a different sport: slower, longer targets in a longer window.

## The shooter

**Recoil hardly matters, except for a shotgun.** A shot's impulse gives the shooter the same speed change on both worlds, but the feet's grip and the floor's moment at the support's edge are a sixth as strong on the Moon. The impulse is 0.29 N s for a paintball marker, 1.2 for a bow, 0.09 for an air rifle and 1.0 for a .22. A 24 g trap load at 400 m/s, 11.2 N s with its powder gases, rocks a standing lunar shooter back 1.5° and takes 0.56 s to settle (Earth 0.26° and 0.09 s; the larger body 0.55° and 0.35 s). A shotgun's second shot comes half a second after the first, so on the Moon it is fired from a body still settling.

**A lunar shooter stands still without trying.** A standing body is an inverted pendulum on its ankles. On Earth the ankles' intrinsic stiffness is 91 ± 23% of what gravity's toppling needs, so standing still takes continual active control (Loram and Lakie 2002). Carried to the Moon, the same stiffness is 5.5 times what the toppling needs (4.1-6.9 across that spread), so the body stands passively stable and sways, if disturbed, with a period of about 2.4 s. This is a screen of one mechanism. How a shooter's nervous system retunes its stiffness at a sixth of the weight is unmeasured, and holding a sight steady also depends on the arms and the breath.

## Dodging a paintball

A target that sees the shot starts to move 0.2 s after it leaves the muzzle, under the reach engine's floor law, and escapes once its centre is the torso's half-width plus the ball's radius (0.21 m) off the ball's line. On the Moon at a grip of 0.8, that takes a shot fired from at least 41 m at a target already moving with a foot down. A target in the ready crouch needs 45.5 m, one mid-bound 48.7 m and one standing upright 49.7 m. On Earth the same targets need 29.7, 31.5, 40.5 and 40.1 m. The lunar ball is slower, but a lunar player starts to move more slowly still: by 30 m the ball arrives in 0.48 s, and a player in the ready crouch has moved 3 cm.

**A player in the air is a sure hit.** An arc cannot change, so the shooter leads it exactly. A lunar jumper rising 1 m is in the air 2.2 s, and from 30 m any shot fired in the first 1.7 s of its flight hits. A 2 m rise gives 2.7 s. Earth's highest jump stays up 0.49 s, so from 30 m the window is 0.03 s.

**Ducking is slow.** Sinking 0.5 m behind cover takes 0.83 s on the Moon (0.41 s on Earth), since a body sinks no faster than g. Pulling down on the cover with the arms halves it, 0.36 s (the sport README's "Moving by obstacles", with the walls, ledges and bars the reach engine now has).

For a paintball game this means: at the usual ranges of 10-40 m a lunar shot cannot be dodged once fired, and a jump makes a target. The play goes into not being seen and not being where a shooter expects. That means cover, the moment of breaking from it, a path an opponent cannot predict, and obstacles to move by.

## For a game

The product's `for_games` block gathers what a game engine needs to fly a paintball and move a player within these limits. For the paintball it gives the mass, diameter, drag coefficient and muzzle speed, and the law a = -k |v - w| (v - w) - g z with k = ρ C_D (π d²/4) / 2m, w the wind (k = 0.0242 per metre on the Moon, 0.0211 on Earth). For each world it gives gravity and air density. For the player it gives the reaction time, the floor's push limit at each grip (1.30 m/s² at 0.8 on the Moon), the sprint profile, the top speed, take-off speeds and the hands' reach, with pointers to the reach engine's floor law and jump timing for the rest.

## Evidence and limits

Each projectile has its source or assumption beside it in [sources.json](sources.json).
- **Paintball:** the velocity limit and calibre come from a marker manual's statement of the ASTM standard; the mass and the smooth sphere's drag are assumptions.
- **Arrow:** the drag coefficients are Park's.
- **.22 bullet:** the drag is fitted to the maker's velocity table, to within 2% at 50 and 100 yards.
- **Clay target:** the coefficients are the starting values of a published fit, and the throws are checked against its measurements.
- **Pellet, bow, shotgun load and powder gases:** assumptions.

Everything else is mechanics on the sport study's bodies and reach engine. The flights leave out spin's sideways force (a paintball's curve), a fletched arrow's oscillation and the drag crisis. The dodge leaves out the shooter's aim error and a paintball's own scatter, which spread hits and make long shots chancier on both worlds. Nothing was fired, thrown or dodged; these are screens.

## Reproduction and checks

```sh
python -m humanities.sport.projectiles.aimed      # results/projectiles.json, in about 15 s
OPENBLAS_NUM_THREADS=1 python -m pytest -q humanities/tests/test_sport_projectiles.py
```

The tests check the following:
- that the product regenerates and binds its code, constants and inputs;
- flights without air against the vacuum's closed form, and a fall from rest against the tanh law to the terminal speed;
- the .22 against its maker's table, and the clay against its measured skeet throw;
- that a lunar drop is a sixth of Earth's over the same short flight;
- recoil's momentum, and its rock-back scaling with 1/g;
- that Earth's ankle stiffness stands a shooter still on the Moon but not on Earth;
- that a longer flight never shortens a target's escape, a player in the air never escapes, lunar dodges need more distance than Earth's, and ducking keeps to gravity's clock.

On 8 November: 7 passed.
