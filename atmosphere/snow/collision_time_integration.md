# Finite-time capture updates for the prescribed bath

The maintained-bath controls now have two positive finite-time capture implementations: a direct [Poisson mixture](collision_jumps.py) for constant rates and a [uniformized birth process](collision_birth.py) for size-dependent rates. Both add actual donor masses before reconstructing the three cell moments. This note gives their mathematical basis, approximation bounds and source lineage; the [results report](collision_jump_results.md) evaluates the completed controls. Their scope is a prescribed, replenished bath. Closed aggregation additionally consumes donors and changes rates with the evolving population.

Let the donor mass be \(d\), the constant capture rate \(\lambda\), the step length \(h\), and \(\mu=\lambda h\). Starting from a point mass, forward Euler gives weights \(1-\mu,\mu\) for zero and one captures. Its formal capture-count variance is \(\mu(1-\mu)\), so realizability fails for \(\mu>1\). Multiplying the population by a small positive number changes neither this conditional distribution nor its restriction: at \(100\,\mathrm{s}^{-1}\), the bound is \(h\le0.01\,\mathrm{s}\). An SSP method requiring realizable Euler stages inherits that limitation even when the affected cell contains extremely little population. This derivation explains a restriction; it does not identify every rejection in a particular run.

Positive implicit stepping can alter the variance instead. Backward Euler applied to the unbounded constant-rate birth ladder, again starting from a point mass, gives a geometric capture-count distribution with mean \(\mu\) and variance \(\mu(1+\mu)\). The physical finite-time distribution is Poisson, whose mean and variance both equal \(\mu\). Thus nonnegative probabilities alone do not establish acceptable time integration for this narrow-distribution problem.

The constant-rate implementation uses the transition directly. For initial number density \(n(m)\), extend it by zero outside its represented support and write

\[
n_{\mathrm{new}}(m)=\sum_{k=0}^{\infty}p_k(\mu)n(m-kd),
\qquad p_k(\mu)=e^{-\mu}\frac{\mu^k}{k!}.
\]

Every shifted component is nonnegative. With \(M_r=\int m^r n(m)\,dm\) and \(N=M_0\), before stopping particles at a finite-domain exit the identities are

\[
N'=N,\qquad M_1'=M_1+Nd\mu,
\]
\[
M_2'=M_2+2d\mu M_1+Nd^2(\mu^2+\mu).
\]

The variance of mass per collector increases by \(d^2\mu\). These are local derivations from the Poisson moments. The added mass belongs to the prescribed donor input; the bath is replenished, so collector mass is not conserved by itself.

For each shifted component, positive Gauss quadrature integrates number, mass and squared mass over its overlaps with destination cells, splitting at those boundaries and at intervals resolving the entropy density. Atoms enter at their actual summed masses. Reconstruction occurs after the macrostep and can still change higher moments and subsequent CDFs. The upper boundary is absorbing: the exit ledger records mass at the first crossing and stops subsequent captures. For the constant-rate process, the probability of exiting on capture number K is evaluated from the full Poisson survival probability for each initial mass; the unrestricted identities above therefore apply only before stopping.

For the size-dependent bath, each starting mass defines a ladder \(m_k=m_0+kd\), with rates \(\lambda_k=N_bK(m_k,d)\). A column-oriented generator has \(Q_{k+1,k}=\lambda_k\) and \(Q_{k,k}=-\lambda_k\). On a finite state set, including appropriate absorbing exits, choose \(\nu>0\) with \(\nu\ge\max_k\lambda_k\); an identically zero generator needs only the identity update. Then

\[
e^{Qh}=e^{-\nu h}\sum_{r=0}^{\infty}\frac{(\nu h)^r}{r!}
\left(I+Q/\nu\right)^r.
\]

The matrix inside the power is stochastic, so this uniformization represents the finite-time transition as a positive mixture. [Gross and Miller (1984)](https://doi.org/10.1287/opre.32.2.343) describe the underlying randomization method for transient Markov probabilities and event statistics. Applying it to this snow bath and its sectional reconstruction is a local adaptation. Holding \(\lambda_k\) at its initial value throughout a leap would instead introduce a rate-freezing approximation requiring separate refinement.

## Probability and moment truncation

[Fox and Glynn (1988), §§1–4](https://web.stanford.edu/~glynn/papers/1988/FG88.pdf), provide stable Poisson weights and explicit tail bounds. The local constant-rate code uses a mode-centered recurrence with a directly evaluated peak probability and Chernoff bounds to choose its count window. It does not implement their complete algorithm. Retained probabilities are never renormalized. For a Poisson count of mean μ, the local Chernoff bound at a threshold k on the appropriate side of μ is `exp(−μ + k[1 + log(μ/k)])`, with the k = 0 limit evaluated separately. Long-double arithmetic and upward conversion preserve a positive bound through ordinary double underflow; this remains floating-point evaluation of an analytic bound, rather than interval arithmetic.

For retained counts ℓ through u, the omitted count moments follow by shifting the Poisson sum:

- `q0 = P(K ≤ ℓ−1) + P(K ≥ u+1)`;
- `q1 = μ[P(K ≤ ℓ−2) + P(K ≥ u)]`;
- `q2 = μ²[P(K ≤ ℓ−3) + P(K ≥ u−1)] + q1`.

Replacing the tail probabilities with their Chernoff bounds gives separate upper bounds on omitted number, first mass moment and second mass moment: `N q0`, `M1 q0 + N d q1`, and `M2 q0 + 2 d M1 q1 + N d² q2`. Subsequent constant-rate propagation uses the same positive moment identities as the untruncated transition. Physical exit probabilities restored from the full survival function can reduce the actual omission; the bounds still cover the omitted unrestricted mixture.

[Munsky and Khammash (2006), §II](https://doi.org/10.1063/1.2145882), establish finite-state projection bounds on retained probabilities. The local birth implementation instead uses a direct coupling bound: reaching artificial ladder count C requires at least C events of the dominating Poisson process, so its probability is at most `P(Poisson(νh) ≥ C)`. Rates are evaluated at every actual ladder mass up to C−1. At a physical first exit they become zero. Uniformization retains event counts through R, with R < C, and records its omitted Poisson tail separately from the artificial-ladder bound. Thus an artificial last state is never silently treated as a physical collector.

The physical upper mass boundary U and donor mass d bound every stopped mass by U+d. Accumulated omitted and ladder-hitting probabilities therefore bound each later raw mass moment by multiplication with `(U+d)^r`, for r = 0, 1, 2. This bound is conservative when the represented collectors lie far below U. The size-dependent controls use a per-source Poisson tail tolerance of 10⁻²⁴; the default ladder has 256 events, with a separate 384-event control. A work budget of νh ≤ 64 keeps the forward Poisson recurrence away from underflow. Exceeding it halves the trial macrostep, without accepting an incomplete transition. This is a numerical work budget, not the point-population Euler restriction μ ≤ 1.

These bounds account for transition truncation conditional on the reconstructed source. They exclude entropy reconstruction, source quadrature, floating-point arithmetic and nonlinear amplification of those errors. The code separately records quadrature moment discrepancies, floating-point storage losses and number and water ledger residuals. The constant-rate implementation propagates its storage bounds through the exact transition; the size-dependent implementation records the largest per-step storage bounds. Neither quantity certifies a total solver error. Atomic CDFs are evaluated against physical atom masses so that an independent coordinate subtraction cannot move the comparison to the wrong side of a jump.

## Numerical controls and scope

The [benchmark runner](collision_jump_benchmarks.py) keeps the analytic Poisson and physical-kernel controls separate. The analytic cases sample both sides of discrete jumps within twelve count standard deviations, with omitted reference-tail probability reported. Each mass grid samples its own first crossing, while the original 191-cell case retains the earlier comparison time. Separate runs halve the macrostep, double quadrature order, vary grid spacing and shift grid faces. A further comparison uses one macrostep and two half steps to expose the effect of intermediate reconstruction.

The physical cases retain the original imposed collector, donor, bath density, snow kernel and observation times. Their independent reference consists of 8,192 exact-event histories. Timestep, quadrature, mass-grid and ladder-capacity controls compare against that reference's sampling uncertainty. Number and mass accounting, mean, variance, CDF discrepancy and moment truncation have independent admission criteria. In particular, the truncation check evaluates `(B2 + 2 mean B1 + mean² B0)/variance`; a small raw-second-moment error alone would not resolve this narrow distribution.

The focused tests compare size-dependent transitions and physical first exits with an independent dense matrix exponential, and linear birth rates with negative-binomial moments. Constant-rate tests exercise a macrostep with 20,000 expected captures, continuous source measures crossing faces, and an analytic continuous-source overshoot distribution. Both implementations are tested with population scaling by 1, 10⁻¹⁰⁰ and 10⁻³⁰⁰. They introduce no population cutoff.

The [results report](collision_jump_results.md) supplies measured accuracy and runtime. Those controls establish a maintained-bath component at their tested conditions. Closed aggregation requires paired event accounting, donor depletion and rates updated with the evolving population, followed by independent number-loss, variance and distribution comparisons. The subsequent [full-column event study](event_column_results.md) carries out that population comparison through a separate solver and finds unresolved lunar volume dependence. Physical aggregate survival and predictive snowstorm optics remain unestablished. The [source ledger](collision_time_integration_sources.json) preserves the targeted reading and retrieval limits.
