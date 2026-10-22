# Collision events in the full snow column

This calculation extends the independent collision reference through the complete prescribed lunar cloud. Closed aggregation removes donors at each event; vapour growth, fragment production, seed injection and settling share one water and number ledger. Particle masses remain continuous. The [control matrix](event_column_results.md) determines which computed distributions survive numerical refinement. The inherited aggregate geometry, moist environment and chipping law retain their physical limitations.

The completed finite-time bath method remains a separate tested component. Its replenished donor population does not supply a closed nonlinear collision update. The new [Python interface](event_column.py) and [C++ event kernel](event_column.cpp) instead extend the earlier finite-volume reference, using a faster sampler for the same collision generator. This choice makes the full coupled calculation tractable while introducing a volume-convergence requirement. It does not install a validated deterministic replacement into PlaSim or CM1.

## Collision generator and positive rate selection

Each particle represents one particle in a finite computational volume V. For two different particles, the inherited encounter kernel is

\[
K_{ij}=C(D_i+D_j)^2|v_i-v_j|,\qquad C=q\pi/4.
\]

An unordered pair meets at rate K/V. It sticks with probability E = exp(0.09 T_C), with that decision made separately for every encounter. A sticking event removes both parents and inserts their actual summed mass. With chipping disabled, E multiplies the event rate directly. Equal-mass particles have zero differential-settling encounter rate. No collision batch, statistical particle weight, population cutoff or maximum flake dimension enters the calculation.

The local thinning construction follows from the algebraic bound

\[
K_{ij}\le \overline K_{ij}
=2C(D_i^2+D_j^2)(v_i+v_j).
\]

Write S₂ = ΣD², Sᵥ = Σv, and S₂ᵥ = ΣD²v. Including self-pair proposals, the total candidate rate is `2C(N S₂ᵥ + S₂ Sᵥ)/V`. The first term chooses one uniform particle and one weighted by D²v; the second chooses particles weighted by D² and v. Adding both orderings produces the stated unordered-pair bound. Self-pairs are rejected. A distinct pair is accepted with probability K/\(\overline K\). Exponential candidate waits followed by these null events therefore retain the desired pair rate. The moment trees are rebuilt after growth and changed locally after each collision. Their parent sums are computed from positive children, avoiding subtractive rate drift.

[Eibeck and Wagner (2000)](https://doi.org/10.1137/S1064827599353488) provide the published lineage of fictitious-jump and acceptance–rejection coalescence algorithms. Only their abstract and bibliographic record were read for this extension. The separable bound, tree implementation and column coupling above are local derivations, checked against the independently coded direct sampler and analytic finite-population controls. No theorem from that paper is asserted for the full lunar model.

At fixed masses, air properties and layer membership, the collision substep samples the finite-volume generator without a collision timestep approximation, apart from floating-point arithmetic and random-number generation. The full operator split still has timestep error. Moreover, agreement with a finite-volume reference does not establish an infinite-volume kinetic limit for the extrapolated snow kernel. The earlier twelve-hour closed boxes already exposed unresolved volume dependence.

## Fractional chipping and its allocation remainder

The frozen smaller-parent law prescribes a real-valued loss ℓ per nonsticking encounter. The smaller parent's mass decreases by ℓ; its fragment-number production is ℓ/m_f, where m_f is the inherited 0.1 mm fragment mass. This expectation generally is not an integer. Rounding it independently at each encounter would change the specified mean loss or water budget.

The event calculation accumulates ℓ/m_f in the layer's fragment allocation remainder. Every whole particle is emitted at mass m_f; the remaining fraction stays in that layer's number and mass ledgers. Thus each remainder lies in [0, 1), and the total unallocated mass per square metre is bounded by `layers × m_f / area`. The corresponding particles join growth, settling and collisions after allocation. Their delay is an approximation, rather than exact sampling of a fully specified stochastic breakup law. They are never discarded.

The product also accumulates the unrestricted growth demand of one m_f particle in every layer over the integration interval, divided by area. This bounds the instantaneous demand omitted from fractional remainders accumulated in time. It does not bound the nonlinear effect of delayed fragment growth, transport, nucleation or collisions. Increasing volume reduces the remainder's relative representation; the complete volume controls include this approximation. The code retains the prescribed chip law's millimetre size restriction and maximum chip count. It supplies no whole-aggregate contact failure or mechanical survival calculation.

## Growth, source and vertical transport

The prescribed air is still and saturated over liquid water. It has no resolved updraughts, dry gaps, evolving vapour reservoir or latent-heating feedback. The layer air properties, seed-size interval, cloud-top reference injection speed, vapour supply and integration duration match the earlier prescribed scenarios. Injection is a Poisson stream with expected number flux `N_seed × u_in`; individual seed diameters are sampled from the original continuous log-uniform interval. Consequently the source agrees with the deterministic boundary flux in expectation. Each realization records the actual injected water, including its sampling variation.

For each growth substep, positive Heun stages advance actual masses. At each stage the common layer multiplier is the smaller of one and available vapour production divided by particle demand. Each stage therefore respects the prescribed layer supply. The ledger records the mass actually added after floating-point storage. No mass-space advection grid is used.

The primary matrix retains the earlier upwind height-transport equation. A particle in layer k moves to k−1 with rate v(m,k)/Δz; its waiting time is exponential. The mean occupation thus follows the same finite-volume settling generator as the sectional baseline. The transport substep may cross several layers and refreshes velocity after each crossing. A separate ballistic mode follows a particle's position through piecewise-constant layer air. That mode changes vertical transport and is excluded from the primary matched controls. Initial ballistic pilots were exploratory and are superseded by the source-bound upwind matrix.

A symmetric composition applies half injection, half growth, half settling, a full collision step, half settling, half growth and half injection. Collision masses and rates update after every accepted event. The macrostep ends exactly at the filling boundary and four gathering-block boundaries; separate 60 s and 30 s controls test the split. Layer occupation is accumulated at macrostep endpoints, so these time-mean diagnostic profiles also need timestep refinement. They are not exact residence-time integrals.

The complete ledger is

\[
M_{\rm injected}+M_{\rm deposited}
=M_{\rm freezing\ exit}+M_{\rm remaining}+M_{\rm fractional\ fragments}.
\]

The number ledger adds chipping production and subtracts one particle per sticking event. There is no upper mass-domain exit. A largest particle in one realization is bounded by that realization's available water, which increases with simulated volume. This numerical bound cannot establish a physical maximum snowflake size. The chosen horizontal area specifies the statistical volume of each well-mixed layer; it is not a literal narrow tube with walls confining a flake. Increasing it adds represented particles at fixed concentration. This is a finite-population refinement, distinct from the earlier audit that rescaled weights while keeping the same particle packets.

## Controls, uncertainties and product interpretation

The closed-box matrix uses the existing independent reference's exact initial species, four volumes and five observation times. Each new volume has 128 independent realizations, compared with the reference's 64. Number, second mass moment and every exported mass-CDF query must lie within five combined realization standard errors, with a numerical tolerance. This is a compatibility screen across the stored queries, not a joint confidence guarantee or rare-tail certificate.

The full-column matrix contains four independent realizations per setting. Lunar area increases from 0.0001 to 0.0004 to 0.0016 m² at 30 layers and 60 s maximum step. A separate 30 s setting holds area at 0.0004 m². The 60-layer setting doubles area to 0.0008 m², preserving each layer's statistical volume while halving vertical spacing. Earth controls use the same geometry and 100 seeds/L, with 1 mm/h vapour supply and two areas. These generic Earth columns are transfer controls. The earlier published-simulation comparisons remain separate evidence; this matrix is not another reproduction of those archived Earth profiles.

Declared numerical screens require all cases to complete, relative water and number residuals below 10⁻¹⁰, and each listed refinement to change the mean realization mass CDF by at most 0.05 and mean flux and mass quantiles by at most 5%. These practical tolerances do not establish physical correctness or bound finite-sample uncertainty. Failed controls remain in the report. In particular, agreement in the water ledger cannot compensate for a distribution that shifts with volume.

Standard errors are computed across complete independent realizations. Aggregating particles interact and compete for a common water supply, so individual particles do not provide independent uncertainty samples. Reported mean quantiles are means of realization quantiles; they differ from quantiles of a pooled distribution. The four arrival blocks show changes during gathering but do not themselves establish a stationary cloud or four independent samples.

The versioned products retain each case's settings, completed duration, ledgers, cumulative arrival number/mass/second-moment/projected-area distributions, height profiles, source and climate-input hashes, compiler flags and compiled-binary identity. The final null diameter in the query list denotes +infinity and includes every arrival. Large dimensions are never deleted from the water or optical-area moments. Numerical budgets reject an unfinished calculation with an explicit failure; they are not physical cutoffs.

Focused tests check physical-formula parity, analytic finite-volume death chains, two-particle exponential waiting, direct-sampler agreement, chipping-loss parity, noninteracting vertical transport against an independent matrix exponential, source statistics, source-only retention, complete shallow-column ledgers, supply limits and reproducibility. The [results report](event_column_results.md) records which checks and matrix runs were actually completed. Neither the source hashes, successful accounting nor these numerical controls validate giant aggregate geometry, fracture, real cloud exposure or snowstorm optics.
