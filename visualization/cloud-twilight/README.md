# Cloud-twilight diagnostics

`python visualization/cloud-twilight/plot.py` reads the committed
`research/studies/cloud_twilight/results/cloud_twilight.json` product and writes
an ignored four-panel PNG with a product/plot hash manifest. Panels show incident
cloud illumination, separate ground-light context, the spectrally combined
Sun/cloud/eye transmission factor, and corrected cloud occurrence by height/time.

Cloud radiance and visual contrast are not calculated. The plot keeps the
current design-column beam separate from the old unfiltered sky background.
The source fraction in panel C has no cloud scattering law; panel D is a
per-cell occurrence fraction, not a probability of seeing a cloud.
