# Bounded maximum-entropy reconstruction for snow collisions

The alternative to the [beta candidate](collision_moments_methods.md) reconstructs each mass cell with a bounded exponential quadratic while retaining its number, mass and second mass moment. It addresses the choice of unresolved distribution shape. The [collision-repair report](collision_repair.md) owns the comparisons with exact capture distributions, including the beta candidate's cell-crossing failure; neither moment agreement nor the maximum-entropy prescription establishes the new candidate's accuracy.

This family has established mathematical and numerical precedents. [Mead and Papanicolaou (1984), §II](https://bayes.wustl.edu/Manual/MeadPapanicolaou.pdf), derive exponential-polynomial reconstruction from constrained entropy and its convex dual. [Massot et al. (2010), §2.2.2 of the inspected preprint](https://www.researchgate.net/publication/40534668_A_Robust_Moment_Method_for_Evaluation_of_the_Disappearance_Rate_of_Evaporating_Sprays), use bounded maximum-entropy reconstruction within a sectional population method for evaporating sprays. The proposed application to binary snow collisions, using the existing event accounting and a three-moment distribution in each mass cell, is a local adaptation. These papers do not validate that particular solver.

For normalized mass coordinate \(x=(m-a)/(b-a)\in[0,1]\), write

\[
p(x;A,B)=\frac{\exp(Ax+Bx^2)}{Z(A,B)},\qquad
Z(A,B)=\int_0^1\exp(Ax+Bx^2)\,dx.
\]

Particle density in mass coordinates is \(N p(x)/(b-a)\). Given normalized mean \(u\) and variance \(v\), choose \(A,B\) to satisfy \(E[x]=u\) and \(E[x^2]=u^2+v\). Equivalently, minimize

\[
\Phi(A,B)=\log Z-Au-B(u^2+v).
\]

Its gradient is the pair of moment residuals and its Hessian is the covariance matrix of \((x,x^2)\). This is the three-moment specialization of the published entropy construction. Entropy here uses a uniform reference measure in mass. Affine centering and scaling preserve that choice; maximizing differential entropy in logarithmic mass would define a different closure. No thermodynamic selection law for falling snow is asserted.

The quadratic coefficient must remain unrestricted. Every finite \(A,B\) gives an integrable positive density on this bounded interval. Positive \(B\) makes the log density convex; when \(-2B<A<0\), the density has an interior minimum and rises toward both endpoints. For the symmetric case \(u=1/2\), \(A=-B\); \(B=0\) gives the uniform variance \(1/12\), and positive \(B\) permits larger variance. A restriction to negative \(B\), as in a truncated-normal parameterization, would exclude valid moment states. These shape properties follow directly from the proposed formula.

The realizable interior is \(0<u<1\), \(0<v<u(1-u)\). Its boundaries require atomic measures: \(v=0\) gives a point mass at \(u\), while \(v=u(1-u)\) gives weights \(1-u,u\) at \(0,1\). Empty cells and endpoint means also need explicit treatment. Finite coefficients cannot reproduce these atoms. [Alldredge et al. (2014), §§2–4](https://arxiv.org/pdf/1306.2881), distinguish interior entropy solutions from atomic boundary measures and examine the resulting ill-conditioned optimization. The implementation retains explicit atomic limits without imposing a minimum variance.

Numerical inversion must be assessed separately from closure accuracy. The Alldredge study shows that a finite quadrature has its own realizability set and that moving moments inward through regularization can impede convergence. Its fixed-quadrature preference concerns that algorithm. The [local implementation](collision_entropy.py) instead partitions the support to resolve interior and endpoint concentrations and shifts exponent evaluation by the peak.

The implementation fits in standardized coordinates \(z=(x-u)/\sqrt v\), targeting \(E[z]=0\) and \(E[z^2]=1\). Its optimizer target is \(2\times10^{-11}\). Final parameters, including a stalled optimizer's result, are admitted only when a separately refined rule gives standardized moment residuals at most \(2\times10^{-10}\), and changing from 24- to 48-point Gauss rules on the partition changes log normalization by at most the same tolerance. Termination status remains recorded. Each cell's preceding finite-density fit is transformed into the new centered/scaled coordinates as an initial guess; failed warm starts fall back to the cold guess. Exact target-moment pairs share a fit cache. These initial guesses and cache entries do not change the target moments or admission criteria.

For collision quadrature, the code compresses a positive sampled measure by twice-reorthogonalized polynomial iteration and a symmetric eigensystem. It checks positive finite weights, support and sampled moments through degree \(2q-1\), where \(q\) is the requested quadrature order. These checks concern the sampled rule; quadrature refinement must still assess the underlying density integrals and nonlinear kernel.

At subnormal storage magnitudes, a separate reconstruction check propagates half the smallest floating-point spacing through normalization of the stored number and two local moments. It permits only defects within the resulting roundoff bounds, leaves the stored state unchanged, and records affected population, physical moment bounds and reconstruction changes. A feasible positive variance is not removed merely because this additional storage uncertainty includes zero. No population cutoff is imposed. Diagnostic counts include cache reuse, repeated validations and rejected trial stages; they are not counts of independent fits or physical collisions.

Time integration retains SSPRK2 and validates both Euler stages and their final average. Writing that average as \(Y+\tfrac12(Y^{**}-Y)\), with \(Y^{**}\) the second Euler stage, preserves identical inputs even when separately rounded subnormal halves would not. It does not remove the Euler-stage realizability restriction. The [finite-time jump update](collision_time_integration.md) remains a separate bath-only proposal. These numerical safeguards do not establish distribution accuracy or physical snow behaviour.

The collision operator still has the binary-event budget \(\Delta(N,M_1,M_2)=(-1,0,2mm')\). A successful entropy inversion does not guarantee accurate nonlinear kernel integrals or mass transferred across cell boundaries. Quadrature must resolve those boundaries, and time integration must keep all stages realizable. Three moments also cannot recover arbitrary multiple peaks or the discrete capture ladder. Repeated reconstruction changes higher moments and can alter later collision rates even when the low moments remain accurate.

The decisive comparisons therefore retain independently known distributions: the constant-rate Poisson bath, the original size-dependent capture reference, crossings at several positions relative to mass-cell edges, and closed-box aggregation. They should vary timestep, mass-cell width and quadrature independently, reporting mean, variance and CDF errors through the crossing rather than only after recovery. Reconstruction-only tests should include uniform, narrow interior, skewed, positive-\(B\) and atomic cases. An inversion check on moments generated from the same exponential family establishes numerical consistency; comparison with the independent capture distributions tests the closure.

The [source ledger](collision_entropy_sources.json) records the exact reading extent, accessed versions and limitations. No external paper bytes or numerical production results accompany this note.
