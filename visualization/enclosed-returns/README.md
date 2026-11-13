# Staged construction and export returns

This renderer reads the admitted [return study](../../research/studies/enclosed_returns/README.md) and displays its operating cash, pooled funder cash returns and additional funding requirements. The figure uses the breakthrough customer scenario and identifies the real-demand-growth assumption. It carries the source-access and deployment qualifications with the comparison.

Run `python visualization/enclosed-returns/render.py`. The PNG, SVG and numerical plotting manifest go to ignored `research/runs/enclosed_returns/figures/`. The manifest records source and renderer hashes and every plotted value. The renderer does not calculate production, financing or socioeconomic benefits.

The [bottleneck renderer](bottlenecks.py) reads the separate stratified product and writes `research/runs/enclosed_returns/bottlenecks/bottlenecks.png`, an SVG and a manifest. It displays the discounted receipt-to-payout reconciliation, the timing of earnings, intervention NPVs and source-tug replacement sensitivity. Workforce failures appear as unresolved cells. Reproduce with `python visualization/enclosed-returns/bottlenecks.py`; the renderer performs no economic calculation beyond grouping already reported monetary strata for display.

The [revised renderer](revised.py) displays the matched recovery/refurbishment comparison, operating trajectories and cumulative payouts from the revised financial products. Run `python visualization/enclosed-returns/revised.py`; PNG, SVG and the source/value manifest go to ignored `research/runs/enclosed_returns/revised_figures/`. The figure keeps 0.5% and 1% real demand growth distinct and excludes the separate accelerated-rollout interventions. It calculates no additional economic return.
