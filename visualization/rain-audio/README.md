# Rain audio

An audio lab for the Open Moon's rain: it renders how rain sounds on concrete, puddles and grass with leaves in the
three target climates, the hot world, the resort and the lock, beside the same rain on Earth, in mono and in stereo,
near and far, and how thunder sounds below its flash. The sound starts from Earth rain textures fitted to CC0 field
recordings, which the author accepted as the listening baseline, and carries them to the Moon with the drops the
atmosphere domain computes arriving at the ground.

It is the rain work of 11–13 October 2026, done in ChatGPT and recovered on 17 October (the assessment of that day on
the research drive, `integration-history-review-2026-10-17/assessment.md`, finding F02), brought in cleaned up. Its
renderer and models reproduce the archive's renders bit for bit; its scenes now read the domain's products in place of
the archive's prescribed drop sizes and its air of 20 °C and 80% humidity; and the storm study's thunder layer of
11 October comes with it.

```sh
python3 visualization/rain-audio/render.py --output <folder>   # about three minutes on one thread, 1.2 GB
python3 -m pytest visualization/rain-audio/tests                # about a minute
python3 visualization/rain-audio/fetch_inputs.py --source <folder>   # the recordings, to refit the models
```

Its tests run on their own, like the lane's other Python tools'. Renders go to the folder given, here ignored by Git;
those of 17 October are on the research drive in `atmosphere-rain-audio-2026-10-17/renders/`.

## What it reads

| Input | What the lab takes from it |
|---|---|
| [atmosphere/rain/results/rain_arrivals.json](../../atmosphere/rain/results/rain_arrivals.json) (`terluna.atmosphere.rain-arrivals/1`) | For each target, on the Moon and on Earth: the drops of each bounded cloud-base spectrum that reach the ground, class by class, with their size and impact speed and their number per mm/h of ground rain; drips from leaves; 8–12 mm drops; the air at 2 m |
| [atmosphere/electricity/thunder.py](../../atmosphere/electricity/thunder.py) | The air's absorption by ISO 9613-1, its oxygen and nitrogen relaxation scaled to the air's own fractions (`absorption_db_per_m`), and the sound speed in moist air (`sound_speed`) |
| [research/studies/atmospheric_electricity/results/thunder.json](../../research/studies/atmospheric_electricity/results/thunder.json) (`terluna.research.atmospheric-electricity-thunder/1`) | Each flash's first sound, the span of its sound and the loudest second of each octave band below it |
| [models/](models/) and [models/models.json](models/models.json) (`terluna.visualization.rain-texture-models/1`) | The three texture models, 2.3 MB, with their sources' hashes, credits and fit intervals, the archive's recorded samples, and the stereo spacing's basis |
| [inputs.json](inputs.json) | The three CC0 recordings behind the models, restored by hash from the research drive or the rain work's archive; needed only to refit |

The lab keeps the models in Git as a small product, in keeping with the lane's committed assets of a few megabytes,
because every render needs them and they regenerate only from recordings the repository does not carry. The
recordings (27 MB) stay on the drive and come back by `fetch_inputs.py`, whose hash check stops on any changed byte, in
the manner of `protection/fetch_inputs.py`. Rendered audio stays on the drive.

The recordings, each CC0 1.0 on its source page as checked on 17 October 2026: *Rain on concrete* (BigSoundBank 1289)
and *Rain on Puddle* (1290) by Joseph Sardin, BigSoundBank / LaSonotheque, recorded with a SoundDevices MixPre-3 and a
Sennheiser ME66; and *Rain on grass and leaves* by SamsterBirdies (Freesound 584268), recorded with a Tascam
DR-100mkIII, of which the public high-quality MP3 preview is used. CC0 asks no credit; this gives it.

| File | Holds |
|---|---|
| [texture.py](texture.py) | The Earth textures: band envelopes on the ERB-number scale, their Gaussian-copula statistics with smoothed complex modulation cross-spectra, version 3's intensity correction, fresh carriers; fitting from a recording |
| [scenes.py](scenes.py) | A scene's rain from the arrivals product: direct rain, drips and large drops sharing one water budget, and the air |
| [transfer.py](transfer.py) | The impact-force spectra and their ratio, the minimum-phase filter that applies it, and the large drops' events |
| [propagation.py](propagation.py) | Direct paths to a spaced pair of microphones with the thunder model's absorption; the rain sheet's 2 x 2 covariance |
| [thunder_layer.py](thunder_layer.py) | Thunder below its flash, from the thunder product's band levels, delay and span |
| [levels.py](levels.py) | Loudness by ITU-R BS.1770-4, peaks, fades and WAV files |
| [render.py](render.py) | The auditions for the three targets and their manifest |

## How the sound is made

**The Earth texture.** Each model holds how one recording's sound fluctuates in 64 bands: 257 quantiles of each band's
amplitude envelope, measured every 0.25 ms, and the smoothed cross-spectra of the bands' modulations, which keep how
they rise and fall together and drop the order in which they did. A render draws fresh innovations with those
statistics, maps them through the quantiles, corrects the summed intensity with a monotone map calibrated on a pilot
render, and lays the envelopes on fresh noise. It reads no recording.

**The rain at the ground.** A scene takes the arriving drops of one bounded cloud-base spectrum, the rain product's
central one (mass-weighted mean diameter 1.5 mm, shape 3) by default, at a rain rate chosen at the ground, 5 mm/h by
default, as in the rain work. On grass and leaves a quarter of the water arrives as 3 mm drips released at rest 1 m up.
Each target pairs the Moon's rain with Earth's rain falling through air of the same temperature and humidity at 1 atm,
as the rain product pairs them; the Earth texture stands for that Earth rain.

**The transfer.** Each impact delivers its drop's momentum J = mU as a force pulse lasting a fraction of D/U, a gamma
pulse with a square-root rise and an exponential tail (after Soto et al. 2014 and Mitchell et al. 2019), wider by 1.6
on puddles. The Moon's rain is the Earth texture filtered by the square root of the ratio of the two scenes' summed
force spectra: a causal minimum-phase filter that keeps the texture's timing and its material's ringing. Slow drops make
long, gentle pulses, so the Moon's rain loses most at high frequencies.

**Stereo and distance.** Rain falls on a sheet from 0.2 to 40 m around two omnidirectional microphones 8 cm apart and
1.6 m up. Each of 4,096 cells reaches them with 1/r spreading, a delay of r/c and the thunder model's absorption in the
target's air, and the sum is a 2 x 2 cross-spectral covariance by frequency; the stereo rain is two texture renders
mixed by its square root. The source is calibrated once, to Earth's field, and kept for the Moon. A single 1 m² patch
can also be placed at a point, with its absolute delay.

**Thunder.** The thunder study traced the electrified box's flashes to the ground through the box's air. The layer
renders a flash heard below it: noise in each octave band that reaches the threshold of hearing, and in the 16 Hz band,
under one envelope stretched to the product's span, each band scaled so its loudest second has the product's level.

## What it renders

The renders of 17 October, at 5 mm/h on the ground in each target, are on the drive with a manifest of every file's
hash, level, gain and source. Per target there are Earth-then-Moon comparisons in mono and stereo, each 12 s passage at
the rain work's listening level of −35.4 LUFS (`*-matched.wav`), and the whole 30 s renders with one gain, keeping the
model's relative level (`*-shared-gain.wav`); the Moon's rain on concrete with 1% of its water in 8–12 mm drops; and a
grass patch near, middle and far. Across the targets each surface has the Earth texture followed by the three Moons at
one gain, and four thunder events play at one gain and each at −23 LUFS. A gain that would lift a true peak above
−1 dBFS stops there.

The rain arriving at 5 mm/h on the ground, central spectrum, from the arrivals product:

| | Hot world | Resort | Lock | Earth, in the resort's air |
|---|---|---|---|---|
| Drops arriving per m² per second | 757 | 1,066 | 1,016 | 1,542 |
| Their median size and speed | 1.23 mm, 1.40 m/s | 1.03 mm, 1.18 m/s | 1.09 mm, 1.25 m/s | 0.91 mm, 3.52 m/s |
| Rain leaving the cloud base that it takes, independent drops | 289 mm/h | 72 mm/h | 27 mm/h | 6.5 mm/h |
| A 3 mm drip released from 1 m strikes at | 1.60 m/s | 1.59 m/s | 1.59 m/s | 4.01 m/s |

The transfer and the levels it gives, the Moon against Earth in the same target's air:

| Surface | Target | Transfer at 100 Hz; 1 kHz; 10 kHz | Moon minus Earth, mono: level; loudness |
|---|---|---|---|
| Concrete | Hot world | −6.2; −13.5; −20.1 dB | −17.7 dB; −18.7 LU |
| | Resort | −6.9; −14.2; −20.7 dB | −18.4 dB; −19.4 LU |
| | Lock | −7.5; −14.8; −21.2 dB | −18.9 dB; −20.0 LU |
| Puddles | Hot world | −6.6; −16.1; −20.2 dB | −18.0 dB; −19.0 LU |
| | Resort | −7.3; −16.7; −20.8 dB | −18.7 dB; −19.7 LU |
| | Lock | −7.9; −17.3; −21.3 dB | −19.3 dB; −20.2 LU |
| Grass and leaves | Hot world | −6.9; −13.9; −20.0 dB | −16.7 dB; −17.9 LU |
| | Resort | −7.4; −14.5; −20.6 dB | −17.3 dB; −18.5 LU |
| | Lock | −7.9; −15.1; −21.2 dB | −17.9 dB; −19.1 LU |

At the same rain on the ground the Moon's rain is 17–19 dB quieter than Earth's and darker, its highs falling 13–14 dB
more than its lows, because its drops strike at about a third of Earth's speed: each carries less momentum and presses
three to four times as long. The hot world is the loudest of the three targets, since its drier, deeper fall leaves only
larger, faster drops. The stereo renders give the same differences within 0.1 dB. The Earth texture stands for a
slightly different Earth rain in each target's air, and the three Earth scenes' force spectra differ by at most 0.45 dB.

The Moon's air, with less oxygen and more nitrogen than Earth's, absorbs about 3.5% more than Earth's mixture would
at 125 Hz and 5–14% less at 4–8 kHz; at 1.2 atm the same humidity holds less vapour. Over 40 m it takes 0.2–0.3 dB at
1 kHz and 3.2–4.2 dB at 10 kHz, against 3.7–4.8 dB in Earth's air of the same temperature and humidity. In 30 s on a
4 m² patch the large-drop audition carries two drops of 7.8–8.9 mm striking at 4.8–5.1 m/s in each target, where
3.8–4.1 are expected; each carries 1,300–2,600 times the impulse of the median arriving drop and stands far above the
quiet lunar bed.

The thunder events, below their flashes, with central acoustic efficiency:

| Event | First sound after the flash | Span of 90% of its energy | Loudest second, A-weighted |
|---|---|---|---|
| The Moon, a typical flash in cloud at 30–42 km | 87 s | 34 s | 64.0 dBA, its energy at 16–63 Hz |
| The Moon, a typical ground strike | 2 s | 27 s | 95.4 dBA |
| Earth's benchmark flash in cloud at 4.8–10.8 km | 15 s | 17 s | 76.5 dBA |
| Earth's benchmark ground strike | 2 s | 14 s | 88.3 dBA |

In the shared-gain thunder file full scale stands for 23.2 Pa, 121 dB re 20 µPa, so its digital levels are pressures.

## Checks

| Check | Result | Evidence state |
|---|---|---|
| The texture port against the archive's raw 30 s masters, one per surface, and the lunar renderer's Earth render | Bit for bit; a 3 s render matches the recovered code in this environment | Reproduction of the original |
| Refitting the concrete and puddle models from their recordings | The same arrays, the puddles' modulation factors within 1.1e-16 (the grass recording is an MP3, which needs ffmpeg, absent here) | Reproduction from the sources |
| Envelopes generated over 60 s against each model's fitted quantiles | Every band's 10th, 50th and 90th percentiles within 10%, their medians across bands within 3% | Statistical reproduction |
| The transfer port, given the archive's own scenes | The archive's saved transfer curves bit for bit | Reproduction of the original |
| The transfer's scalings | Zero frequency at 2 Σ j (mU)²; speed doubled gives 4 at low and 32 at high frequency; fewer, larger drops at fixed rain force more at low frequency (s³); flux and width as they must; identical scenes leave the texture untouched | Analytic limits of the model |
| The propagation's controls | Mirror and centreline symmetry, distances, the arrival order, 1/r, √area, the annuli's power and the field's covariance (Hermitian, positive, rank one for coincident microphones, converged in azimuth to 1e-6) | Geometric and numerical checks |
| Absorption | The thunder model's own, whose ISO 9613-1 table check lies in atmosphere/tests | Identity with the domain |
| The thunder layer | Each band's loudest second within 0.01 dB of the product, within 0.5 dB when filtered back out of the whole event; the envelope's span within 1% | Faithfulness to the product |
| The loudness meter | A full-scale 997 Hz sine at −3.01 LKFS and EBU Tech 3341's first case at −23.0 LUFS | Standard reference readings |
| The arrivals product (atmosphere/tests/test_rain_arrivals.py) | Survival shares equal the rain product's to 3e-16; each spectrum's fluxes deliver 1 mm/h; drips below the vacuum and terminal speeds | Numerical checks of a runnable calculation |

The lunar sound itself is a rendering of an assumed transfer: no lunar rain has been heard, and listening judges it.

## What the listening of 11–13 October settled

The rain work's first generators of 11 October, one built from short impact kernels and one from recorded fragments,
were rejected: they sounded static, with tones like insects, and the author asked for an Earth-accurate generator before
any lunar drop work. A reconstruction that kept the recordings' envelopes on fresh carriers was accepted as a step. The
independent generator of the same day, fitted from statistics alone, was accepted for concrete and puddles, tentatively;
its grass and leaves first sounded slightly sharper than the recording, and the measurements found the generated passage
6.5 LU louder than the original; once the original, the reconstruction and the generator were matched at −35.4 LUFS the
author judged it improved and opened the lunar work. These three textures are the Earth baseline (the register's row of
17 October), and comparisons are made at matched loudness, with shared-gain versions kept for relative level.

On 13 October the author set bounded, assumed drop-size distributions in place of a CM1 or cloud-growth campaign, and
asked for stereo and distance next, Earth before the Moon (both in the register). The lunar auditions and the stereo
and distance auditions of 13 October carry no recorded verdict, and these renders none yet. The storm study of
11 October was heard in pieces: its drop demonstration was inaudible on the author's playback and its thunder could not
be picked out, so separate auditions followed, with no recorded verdict. Its thunder layer comes in here because the
product it rendered is unchanged (the same sha256) and the layer renders that product's levels faithfully; its rain
layer gave way to the rain generator. The drive catalogues these records
(`atmosphere-rain-audio-2026-10-17/LISTENING-RECORDS.md`).

## Assumptions and limits

- **The drops.** The bounded cloud-base spectra are scenario boundary conditions, and the drops fall independently in
  the rain product's parcel column. At 5 mm/h on the ground that column asks 27–289 mm/h of the cloud base, where
  coalescence and the rain's own moistening of the air would let more through; the rain product's encounter diagnostic
  says collisions would reshape the arriving spectrum.
- **The humidity.** Each target's air at 2 m comes through the arrivals product from the climate targets' wet-bulb
  temperature at 2 m, which pairs the GCM's 2 m air temperature over each 3-day window with the lowest level's specific
  humidity. In the GCM's output files that humidity (`hus`) is a single-step snapshot at each window's boundary, and
  `domain/climate` is correcting it. The humidity sets the cloud base, how much of each drop evaporates on the way down
  and the air's absorption; the products keep it as it stands until the integration.
- **The rain rate.** 5 mm/h at the ground is the scene's choice, the same on Earth and the Moon; the rain at any other
  rate scales from the same spectrum, and the texture's own timing keeps its Earth statistics.
- **The force pulse.** Its gamma shape, its widths of 0.24, 0.32 and 0.44 D/U and the puddles' factor of 1.6 are
  declared approximations. Incident momentum stands for the impulse delivered, and the surface's response to force and
  the path to the microphone are taken as the same on Earth and the Moon, so they cancel in the ratio; this holds best
  for hard ground, less well for water and leaves, and puddle bubbles are left as the Earth texture has them.
- **The Earth texture.** It stands for Earth rain of the product's central spectrum in each target's air; the
  recordings' own rain, drops and microphone distances are unmeasured. Its level is digital, so the rain has no
  absolute sound pressure.
- **Drips and large drops.** The quarter of the water dripping as 3 mm drops from 1 m awaits biosphere/canopy's
  interception and storage. The 1% of water in 8–12 mm drops, its 4 m² patch and its events' decays are a sensitivity.
- **Stereo.** The 8 cm spacing and the envelope correlation of 0.85 are effective values fitted to the grass
  recording's coherence (models/models.json); ground reflection, refraction by wind and temperature, scattering, an
  impact's radiation pattern and a listener's head are left out.
- **Absorption.** ISO 9613-1's relaxation frequencies are Earth's, with only their strengths scaled to the air's oxygen
  and nitrogen, as the thunder model takes them.
- **Thunder.** The product traces the electrified box's flashes through the box's own air, the old design climate's; no
  target climate's thunder is computed. The envelope's fine shape and the noise are synthesized; the bands, the delay
  and the span are the product's.
