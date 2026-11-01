# Independent collision references

The fixed-pivot collision operator spreads the mass distribution strongly when a large collector captures much smaller particles, even while its mean growth agrees closely with an independent continuous-mass calculation. This is a demonstrated numerical effect in an aggregation-only control. It prevents attributing the full column's large tail entirely to the inherited geometry. The present controls do not determine how much of that coupled tail is numerical, and they do not establish a physical maximum snowflake size.

The producer is [`collision_reference.py`](collision_reference.py); [`results/collision_reference.json`](results/collision_reference.json) binds its sources, shared constants and atmospheric inputs. It leaves the earlier producers and products intact. Both sides use the same Brown–Francis mass law, fixed projected-area ratio, drag law and temperature-dependent sticking at −10 °C in the lunar column. The prescribed air has pressure 52.485 kPa and density 0.69020 kg/m³. Deposition, settling and chipping are disabled on both sides to isolate coalescence. This is a numerical comparison of an inherited physical closure, not empirical validation of that closure.

Each independent finite-volume realization contains actual particles. An unordered pair of masses \(m_i,m_j\) coalesces at rate \(K(m_i,m_j)/V\), with an exponential waiting time and an exactly sampled pair. Its product mass is \(m_i+m_j\); its monomer count is the sum of the parents' counts. There are no superparticle weights, collision batches, mass pivots or time steps in this reference. The kernel includes the same sticking factor as the fixed-pivot control,

\[
K(m_i,m_j)=e^{0.09T_C}\,q\,\frac{\pi}{4}(D_i+D_j)^2|v_i-v_j|.
\]

Here \(q=0.26019797\), \(T_C=-10\), and diameter and velocity follow the frozen geometry and drag code. Equal-mass particles have zero differential-settling collision rate. Actual particle positions, turbulent encounters and hydrodynamic interception beyond the inherited area factor are outside this homogeneous model.

The finite-volume process is a Marcus–Lushnikov coalescent. Its generator acting on the empirical number measure gives the discrete counterpart of the Smoluchowski weak form,

\[
\frac{d}{dt}\int f(m)n(dm)=\frac12\iint K(m,m')
[f(m+m')-f(m)-f(m')]n(dm)n(dm').
\]

This identifies the intended kinetic comparison; it is not a proof that every late-time infinite-volume limit exists for the extrapolated snow law. [Norris (1998)](https://arxiv.org/abs/math/9801145) establishes stochastic-to-kinetic convergence under sufficient conditions. [Fournier and Laurençot (2007)](https://arxiv.org/abs/0706.2057) show, for a particular family of strongly gelling kernels, why the treatment of giant particles can distinguish Smoluchowski and Flory limits. Their kernels are not this snow kernel. Their abstracts were read here to establish the scope of that caution; no theorem from those papers is claimed to apply to the inherited law. Finite-volume convergence must therefore be assessed for each reported observable and time.

A separate prescribed-bath calculation makes the unequal-size error particularly clear. One collector starts at diameter 15.8608 cm and captures particles of diameter 0.20934 mm from a maintained bath of 100,000 particles/m³. The exact reference is a pure-birth chain: after \(k\) captures its mass is \(M_0+ks\), and its next-event rate is \(n_bK(M_0+ks,s)\). Each capture supplies exactly \(s\) of external bath mass. This replenished-bath budget is distinct from the closed-box mass budget.

The final observation at 272.627 s corresponds to 10% initial-rate mass growth; the actual changing rates give mean growth of 10.5480%. The 8,192 independent collector histories give the following comparison. Variances are normalized by the initial collector mass squared. The tail is a **number fraction** above \(1.2M_0\), equivalent to diameter above 17.4582 cm.

| Method | Mean mass / \(M_0\) | Variance / \(M_0^2\) | Fraction above \(1.2M_0\) |
|---|---:|---:|---:|
| Continuous mass, 8,192 histories | 1.10548019 | 0.000000395456 | 0 observed |
| Fixed pivots, 96 masses | 1.10551106 | 0.03761784 | 28.3383% |
| Fixed pivots, 191 masses | 1.10549200 | 0.01741612 | 16.1279% |
| Fixed pivots, 381 masses | 1.10548272 | 0.00838767 | 18.6133% |

The 96-pivot mean differs by only 0.00279%, while its variance is about 95,125 times the sampled physical variance. At 381 pivots the variance remains about 21,210 times larger. The variance decreases with mass-grid refinement; a tail fraction at a particular threshold need not decrease monotonically because the distributions remain discrete. The independent 4,096-history reference gives normalized mean 1.10547544 and variance 0.000000385680. The two mean estimates differ by 0.40 combined standard errors. Their number-CDF 95% Dvoretzky–Kiefer–Wolfowitz bands have halfwidths 2.122% and 1.501%, respectively. Thus zero observed tail is not a claim of zero true probability, but sampling error at this confidence level is far too small to explain the fixed-grid tail. These bands apply to number CDFs, not to mass-weighted CDFs.

The pure-birth calculation admits at most 118,324 events per history; the largest observed count was 31,933 across the two ensembles. The maximum rate over all pre-cap states gives a dominating Poisson mean of 42,165.64 over the observation interval. The Chernoff upper bound on hitting the computational event cap has logarithm −45,930.81. Its numeric exponential underflows to zero in the JSON; the retained logarithm states the bound and must not be read as a proof of exactly zero probability. A history reaching the cap would reject the calculation. Both fixed-pivot bath controls integrate a linear generator by a matrix exponential, including the separate bath source and upper-domain sinks. No time-step error explains their spreading: two half exponentials agree with one full exponential to at most 2.23×10⁻¹⁶ in the exported state, and all exported upper sinks are zero to floating-point precision. The maximum relative bath-water residual is 2.34×10⁻¹⁶.

The ordinary closed boxes start from ten equally populated masses between 0.1 and 0.4 mm, all exactly present in the original 96-pivot grid. The 191- and 381-pivot nested refinements retain those original masses bitwise. This declared discrete distribution ensures identical initial masses and concentrations; it is not an exact representation of the continuous log-uniform column injection. Total initial number is 100,000/m³ and ice mass is 0.142016 g/m³. Each reference size has 64 independent realizations. Initial particle counts 320, 640, 1,280 and 2,560 correspond to volumes 0.0032, 0.0064, 0.0128 and 0.0256 m³. Enlarging the volume at fixed initial concentrations tests finite-particle effects separately from ensemble sampling error.

At six hours, the second mass moment and the mass fraction above 1 mm are:

| Method | Second mass moment, kg²/m³ | Ice mass above 1 mm |
|---|---:|---:|
| Continuous mass, 1,280 initial particles | (2.676 ± 0.234)×10⁻¹² | 19.29% ± 0.72 percentage points |
| Continuous mass, 2,560 initial particles | (2.939 ± 0.163)×10⁻¹² | 19.78% ± 0.56 percentage points |
| Fixed pivots, 96 masses | 8.952×10⁻¹² | 23.19% |
| Fixed pivots, 191 masses | 3.581×10⁻¹² | 21.13% |
| Fixed pivots, 381 masses | 2.999×10⁻¹² | 19.87% |

The displayed uncertainties are one standard error across whole independent realizations. Particles within one realization are correlated and are not treated as independent samples. The finer grid approaches the largest-volume reference for these observables. The far tail remains less well determined: the largest-volume reference has 1.614% ± 0.418 percentage points of mass above 3 mm, whereas no mass above 5 mm was observed in its 64 realizations. The latter observation supplies no tight rare-tail bound; the product records a conservative pointwise 95% Hoeffding halfwidth of 16.98 percentage points for realization-level mass CDFs. Agreement in the bulk moments therefore does not validate the far tail.

At twelve hours the independent reference has **not converged in volume**. Its mean second mass moment rises from 2.763×10⁻¹¹ to 1.529×10⁻¹⁰ kg²/m³ across the four increasing volumes. The largest particle allowed by the total mass in those finite boxes ranges from 4.687 to 14.002 mm. These bounds follow from finite simulated volume, not from snow physics. Meanwhile the fixed-pivot second moment falls from 5.944×10⁻⁸ to 1.007×10⁻⁹ kg²/m³ across the three refinements. Both approximations remain sensitive to their numerical controls; the twelve-hour stress point cannot select the physical large-tail outcome or be called a converged kinetic reference. A low upper-domain loss alone does not resolve that failure: the 381-pivot mass exit fraction is only 4.49×10⁻¹⁸ at twelve hours.

Halving the ordinary 191-pivot time step from 30 to 15 s changes the six-hour second moment by 0.0137% and its mass CDF by at most 4.38×10⁻⁶ in absolute fraction. At twelve hours the corresponding changes are 0.1274% and 6.83×10⁻⁵. These are small compared with the mass-grid differences. Maximum relative closed-box mass residuals are 1.87×10⁻¹⁵ for the exact-event ensembles and 7.64×10⁻¹⁶ for the fixed-pivot controls. The largest discrepancy between incrementally maintained and recomputed event-rate sums is 7.00×10⁻¹⁶. Those bookkeeping checks support implementation correctness; they do not establish physical closure accuracy.

The focused tests independently check two-particle exponential waiting, a three-particle analytic constant-kernel death process, unequal pair selection, constant-rate Poisson capture including a tail probability, conservation, identical seed representation, bath-source/domain ledgers, exponential composition and rejection of an inadequate event cap. Product checks bind the saved sources and inputs and verify the exported population/ledger structure. The producer completed in 78.6 s with one numerical worker. Its native fixed-grid populations, independent-collector mass/count arrays, ensemble CDFs, whole-realization standard errors and finite-volume maxima are retained for source-bound visualization.

This establishes a practical repair requirement: the collision algorithm must represent many small captures without converting them into rare jumps across widely spaced mass pivots. An accurately resolved capture or collision treatment can be tested against the present prescribed-bath distribution before returning to the full column. The continuous-mass reference supplies that independent target, but its direct event cost and finite-volume effects make it a diagnostic rather than a ready replacement for the full atmospheric solver. Chipping remains a separate extension: the inherited smaller-parent channel includes fractional expected chip counts, so an event-based implementation must specify their stochastic meaning and preserve the exact remnant rather than silently rounding them. Geometry, interception, contact survival and climate-supported growth histories remain physical questions after this numerical repair.

**Provenance, 25 October 2026.** [results/collision_reference.json](results/collision_reference.json) was built with the
rain solver as it stood before 25 October, whose Chapman–Enskog vapour diffusivity stood 17.6% below Massman's (1998)
measured value. It keeps its values; its record names `atmosphere/rain/drops.py` at that solver's hash, kept byte for
byte as [atmosphere/rain/baseline/drops_2026-10-24.py](../rain/baseline/drops_2026-10-24.py) and listed in
`research/provenance.json`, and its identity checks read that copy ([shared/baselines.py](../../shared/baselines.py)).
