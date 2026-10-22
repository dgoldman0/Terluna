# Moment reconstruction for snow collisions

A sectional beta reconstruction is a candidate repair for the excessive spreading identified by the [continuous-mass collision reference](collision_reference.md). Each mass cell would carry particle number, mass and second mass moment, allowing a narrow population to remain narrow while its mean moves between the old pivots. This note records the mathematical proposal and its source lineage. Implementation results and adoption decisions belong in the separate collision-repair report; this note supplies neither a convergence result nor physical validation of large lunar flakes.

[Yuan, Laurent and Fox (2012)](https://doi.org/10.1016/j.jaerosci.2012.04.003) introduced beta and gamma extended quadrature methods of moments (EQMOM) for population balances, with positive reconstructed distributions and additional quadrature nodes for source integration. [Pigou et al. (2018), Appendix B.6](https://doi.org/10.1016/j.jcp.2018.03.027), describe beta reconstruction and moment-preserving Gauss–Jacobi quadrature. The proposed use of one beta distribution in each fixed mass cell is a sectional adaptation of that method family. The sources examined here do not establish this particular collision discretization's accuracy.

For a cell with bounds \(a,b\), define \(M_r=\int_a^b m^r n(m)\,dm\), \(N=M_0\), mean \(\mu=M_1/N\), and variance \(s^2=M_2/N-\mu^2\). On the normalized coordinate \(y=(m-a)/(b-a)\), let \(u=(\mu-a)/(b-a)\) and \(v=s^2/(b-a)^2\). For an interior realizable moment set,

\[
0<u<1,\qquad 0<v<u(1-u),\qquad
q=\frac{u(1-u)}{v}-1,\quad
\alpha=uq,\quad \beta=(1-u)q.
\]

A beta distribution with these parameters, multiplied by \(N\), reproduces the three moments. These formulas are the direct single-component specialization, rather than a fitted physical width. Empty cells, zero variance and endpoint means require explicit treatment. Zero variance gives a Dirac mass at the mean; maximal variance gives atoms at the endpoints. No minimum width is implied. Local scaled or centered moment arithmetic should protect the small variance from cancellation. Published beta-EQMOM develops atomic limits and secondary quadrature in detail, although [Vikas et al. (2013), §3.2](https://cfd.ku.edu/papers/2013-jcp-radiation.pdf), demonstrate them for radiation transport, not snow aggregation.

For a symmetric collision kernel, the weak collision equation for each destination cell \(c\) is

\[
\frac{dM_r^c}{dt}=
\frac12\iint K(m,m')n(m)n(m')
\left[
(m+m')^r{\bf1}_c(m+m')
-m^r{\bf1}_c(m)
-{m'}^r{\bf1}_c(m')
\right]\,dm\,dm'.
\]

Thus one binary event changes global \((M_0,M_1,M_2)\) by \((-1,0,2mm')\). The repair must evolve those budgets: particle number decreases, mass is conserved and the second moment increases. Quadrature products should enter the cell containing their actual sum mass, with consistent birth and loss weights. Self-pair factors, cell-edge ownership and mass leaving the represented domain need explicit conventions.

Moment preservation alone does not establish distribution accuracy. Gauss–Jacobi quadrature integrates the reconstruction's low polynomial moments exactly in exact arithmetic; nonlinear kernels and discontinuous destination-cell indicators need separate quadrature controls. Reconstructing a beta after each update also closes unresolved higher moments, which can affect subsequent selective losses and tail growth. Time integration must keep every cell realizable without replacing failed moments by an artificial spread. The decisive controls therefore include the maintained-bath mean and variance, narrow populations crossing cell boundaries, unequal-size mergers, and comparisons of distributions as cell width, quadrature order and timestep change independently.

The stochastic alternatives provide distinct checks. [Shima et al. (2009), §§4.1.4 and 5.1.3](https://arxiv.org/html/physics/0701103), preserve expected collision counts but derive packet collision-count variance approximately \(\min(\xi_i,\xi_j)\) times the physical variance. Equal multiplicities therefore do not by themselves certify finite-volume fluctuations. Their multiple-collision count, floor of \(p\) plus a Bernoulli draw of its fractional part, also differs from a Poisson count; the authors identify small \(p\) as the rigorous timestep condition. A stochastic candidate for the present narrow-width test should retain physical donor-mass increments and demonstrate timestep and multiplicity convergence.

An individual collector's mass distribution is a different observable from fluctuations of whole-population counts. [Li et al. (2022), §§2–4](https://www.pvv.ntnu.no/~nilshau/Own_Papers/2022_Li_etal.pdf), obtain close lucky-collector growth-time distributions under several tested multiplicities, while examining depletion and interactions between collectors. Their implementation evaluates all pairs with small collision probabilities. Those results support a continuous-mass stochastic comparison but do not certify arbitrary random-pair settings. [Dziekan and Pawlowska (2017), §§3–5](https://acp.copernicus.org/articles/17/13509/2017/acp-17-13509-2017.pdf), separately show that multiplicity and assumed mixing volume affect fluctuation statistics; their numerical thresholds are specific to their warm-cloud cases.

Constant-number Monte Carlo maintains a sample through cloning or removal while changing the represented volume. [Lin, Lee and Matsoukas (2002)](https://doi.org/10.1016/S0009-2509(02)00114-8) explicitly distinguish this procedure from physical finite-volume number fluctuations. It remains a possible mean-field population solver, with resampling error to control, rather than an automatic variance reference.

The [source ledger](collision_moments_sources.json) records inspected sections and retrieval failures. This targeted reading supports the method choices and their restrictions. It does not resolve aggregate strength, evolving geometry, or the atmospheric conditions required to produce a large flake.
