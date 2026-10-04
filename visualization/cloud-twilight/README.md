# Cloud-twilight diagnostics

`python visualization/cloud-twilight/plot.py` reads the committed
`research/studies/cloud_twilight/results/cloud_twilight.json` product and writes
an ignored four-panel PNG with a product/plot hash manifest. Panels show incident
cloud illumination, separate ground-light context, the spectrally combined
Sun/cloud/eye transmission factor, and corrected cloud occurrence by height/time.

This first plot draws no cloud radiance or visual contrast. It keeps the
current design-column beam separate from the old unfiltered sky background.
The source fraction in panel C has no cloud scattering law; panel D is a
per-cell occurrence fraction, not a probability of seeing a cloud.

`MPLCONFIGDIR=/tmp/terluna-mpl python visualization/cloud-twilight/evening.py`
draws the regional follow-up from the study's computed cloud radiance: the
evening cloud regions (`regions.png`), the selected views in colour and
brightness (`cloud_views.png`), a cloud system through 72 hours with its
contrast (`evening_history.png`) and a distant cloud in the deep evening
(`deep_evening.png`), with a local gallery (`index.html`) whose exposure can be
adjusted for the faint late views. It reads `regional_evenings.json`,
`evening_appearance.json`, `deep_evenings.json` and the `evening_scenes_*.json`
products from the study's results and writes to the ignored
`research/runs/optical_comfort/evening/figures/` with a hash manifest. Display
smoothing and exposure act on the saved XYZ values.
