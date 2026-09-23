# Accuracy labs

Offline, self-contained pages for inspecting the A1–A3 light-transport references
in [illumination/references](../../illumination/references/) beside the
experience's real-time approximations. The labs display computed results; they do
not compute physics of their own.

| Lab | Shows | Built from |
|---|---|---|
| A1 (`index.accuracy-lab.html`, `src/accuracy-lab.js`) | The shared atmospheric profiles, first-order spectral skies and the frozen-cloud benchmark | `illumination/references/data/a1/generated/` |
| A2 (`index.a2-accuracy-lab.html`, `src/a2-accuracy-lab.js`) | Molecular multiple-scattering sky and cloud convergence | `illumination/references/data/a2/generated/`, `A2_VALIDATION.json` |
| A3 (`index.a3-accuracy-lab.html`, `src/a3-accuracy-lab.js`) | Coupled gas-plus-cloud benchmark cases | `illumination/references/data/a3/generated/`, `A3_VALIDATION.json` |

Generate the reference data first (see the references README), then:

```sh
python build_accuracy.py          # Open_Moon_A1_Accuracy_Lab.html
python build_a2_accuracy.py       # Open_Moon_A2_Accuracy_Lab.html
python build_a3_accuracy.py       # Open_Moon_A3_Accuracy_Lab.html
python tools/a1/check_lab.py      # headless browser checks (Playwright + Chromium)
```

The built pages are ignored by Git; `validation/a*/lab-browser.json` records the
browser checks that were run on each published build.
