# Method: fixed-component ERA-GEE with model-based recursive partitioning

This document states the specification and the reasoning behind each design
choice. `VALIDATION.md` reports what was checked.

---

## 1. The base model

Person $i$ is observed at waves $t = 0, 1, 2$; each wave carries $Q$ continuous
outcomes. The base model is a population-averaged multivariate growth model:

$$g(\mu_{it,q}) = \beta_{1q} + \beta_{2q}\,t, \qquad q = 1,\dots,Q$$

Written in ERA form, the two components are

$$f_{\text{const}} = 1, \qquad f_{\text{time}} = t$$

with both component weights **constrained** to 1. The estimated quantity is the
$2 \times Q$ coefficient matrix

$$B = \begin{pmatrix} \beta_{11} & \cdots & \beta_{1Q} \\ \beta_{21} & \cdots & \beta_{2Q}\end{pmatrix}$$

whose first row holds outcome-specific intercepts and second row
outcome-specific slopes. The parameter vector is

$$\theta = \operatorname{vec}(B^\top) = (\beta_{11},\dots,\beta_{1Q},\ \beta_{21},\dots,\beta_{2Q})$$

— all intercepts, then all slopes. Score columns follow the same order, which
is what lets the diagnostics in §4 index into blocks of it.

### Why the components are fixed, not estimated

Under this specification $F = [\,\mathbf{1}\ \ t\,]$ is a **given design
matrix**. Two things a general ERA-GEE implementation would do must therefore
be switched off.

**The weight update.** $w_{\text{const}} = w_{\text{time}} = 1$ are
constraints. "Converges to approximately 1" and "is fixed at 1" are different
models, and only the second is the one specified here.

**The normalisation.** ERA normally standardises the predictor block, the
outcomes, and the components. Each of those breaks the specification:

* standardising $X$ centres the constant column to identically zero, which
  **destroys the intercept component**. In ordinary regression one standardises
  the predictors but never the intercept column of the design matrix; the same
  applies here.
* standardising $F$ turns $t$ into $(t - \bar t)/s_t$, after which
  $\beta_{2q}$ is no longer "change per wave".
* re-standardising the outcomes **inside each node** makes $\beta$ incomparable
  across terminal nodes, which is fatal when the whole point is to compare
  subgroups.

`R/era_gee_fixed.R` therefore asserts $\max|f_{\text{const}} - 1| = 0$ and
$\max|f_{\text{time}} - t| = 0$ on every call, and exposes outcome scaling as a
switch (`raw`, or `global` = scaled once on the full sample with the constants
held fixed at every node). Node-specific re-standardisation is not offered.

### Why raw time

Because $t = 0$ is baseline, $\beta_{1q}$ is the **baseline expected outcome**
and $\beta_{2q}$ the **expected change per wave**. Both rows of $B$ read
directly off a terminal node, which is what makes a fitted tree interpretable
without further transformation.

### Why the intercepts stay in the model

The substantive question concerns rates of change, so it is tempting to fit
slopes only and report levels descriptively. That is wrong twice over.

MOB tests whether **model parameters** are unstable across a partitioning
variable. A quantity reported descriptively beside the model is not in
$\theta$, so it is never tested. Worse, dropping the intercept forces
$g(\mu) = \beta_2 t$, so every trajectory passes through a common origin at
$t = 0$; where subgroups differ in baseline, the slope absorbs part of that
difference. Including the intercept is a **precondition for detecting slope
heterogeneity correctly**, not an optional extra.

---

## 2. Clustering, and what the working covariance does

Rows are person-waves. The working covariance $V_i$ is $Q \times Q$: it models
association among the $Q$ outcomes **within** a person-wave. Wave-to-wave
correlation is not parameterised.

That is a deliberate specification, not an oversight, and the wording matters:

> Wave-to-wave correlation is not explicitly parameterized in the working
> covariance; the fact that repeated observations belong to the same person is
> instead reflected in the cluster-level estimating functions and in robust
> instability inference.

Not "we assume the waves are independent". The supporting standard result: in
GEE, when the working correlation is misspecified, provided the mean model is
correctly specified and regularity conditions hold, the regression parameter
estimator remains consistent and the robust sandwich variance protects
inference.

### Two equivalent implementations

In GEE the independent sampling unit is the subject, and the natural empirical
estimating-function contribution is

$$U_i(\theta) = D_i^\top V_i^{-1}\{Y_i - \mu_i(\theta)\}, \qquad \sum_{i=1}^N U_i(\theta) = 0$$

The fact that the long data has $NT$ rows and the fact that the independent
score units number $NT$ are different statements. Two implementations follow:

| | **A** — used for all results | **B(a)** — theoretical reference |
|---|---|---|
| `mob()` partitions | the person-wave frame | a person-level frame |
| `estfun` | $NT \times 2Q$ | $N \times 2Q$, $U_i = \sum_t u_{it}$ |
| clustering | `mob(..., cluster = id)` | intrinsic; no `cluster` needed |

These are two implementations of one specification, so the question is not
which is better but whether A reproduces B(a). It does — see `VALIDATION.md`
§2. The reason is structural: `partykit:::mob_partynode` computes the OPG
"meat" of the fluctuation test as

```r
meat <- if (is.null(cluster)) crossprod(process)
        else crossprod(as.matrix(apply(process, 2L, tapply, as.numeric(cluster), sum)))
```

so with `cluster` supplied the covariance is built from **cluster-summed
scores** — exactly $U_i = \sum_t u_{it}$.

The precise claim is narrower than "clustered MOB is B(a)". The covariance
calculation matches B(a) exactly. The fluctuation *process* is still built over
the $NT$ rows. For a person-level **categorical** partitioning variable the two
collapse algebraically — within-category score sums are the same whether taken
over rows or persons, and the category weights are unchanged. For a
**numeric** partitioning variable the statistic is a supLM over an ordered
process in which method A walks a grid three times finer, so the two are not
algebraically identical; empirically the difference was in the third decimal.

The decomposition $U_i = \sum_t u_{it}$ is exact and unique here because $V_i$
is block-diagonal across waves. Both `era_gee_mob_fit()` (A) and
`era_gee_mob_fit_subject()` (B(a)) ship in `R/era_gee_fixed.R` so the
equivalence stays checkable.

---

## 3. Node membership

`era_gee_mob_fit()` reads node membership from `x[, "rowid"]`, which requires
`rowid` on the regressor side of the formula:

```r
cbind(y1, ..., yQ) ~ rowid + X_time | z1 + z2 + ...
```

The reason is a partykit interface detail with real consequences. At every node
below the root, partykit passes `weights` **already subset to the node**, all
strictly positive. A fit function that computes `which(weights > 0)` therefore
gets `1:n_node` and selects the first `n_node` rows of the *whole* dataset
rather than the node's members. Only the root would be fitted correctly.
Carrying an explicit row index sidesteps this: `x` **is** correctly subset at
every node.

---

## 4. The three diagnostics

All three fit the same $2 \times Q$ model on raw time. Only the parameters
entering the instability test differ.

| tree | question | test |
|---|---|---|
| joint | *which subgroups have different trajectories?* | all $2Q$ parameters |
| slope | *which subgroups change at different rates?* | $Q$ nuisance-adjusted slope scores |
| intercept | *which subgroups start at different levels?* | the $Q$ intercepts, `parm = 1:Q` |

The joint tree is the **primary** analysis. It is invariant to the centring of
time, as it must be — centring is a non-singular linear reparameterisation of
the same model space.

### The nuisance-adjusted slope score

A naive slope-focused test — simply restricting the instability test to the
slope columns — does **not** test slope heterogeneity under a raw time origin.
The slope score for person $i$ is $\sum_t t\,V^{-1} r_{it}$. A subgroup with a
pure baseline shift has residuals displaced by a constant $\delta$ at every
wave, and

$$\sum_t t\,\delta = \delta(0+1+2) = 3\delta \neq 0$$

so a pure intercept difference mechanically produces a non-zero slope score.
Measured on the reference data, the two score blocks correlate at about
**0.948** under raw time. No column selection can separate blocks that
collinear, the more so because partykit whitens the process *before* subsetting
its columns, making each selected coordinate a combination of all $2Q$ raw
scores.

Centring time removes the leakage (block correlation drops to ~0.02), but that
is a workaround rather than the answer. The hypothesis of interest,

$$H_0:\ \beta_{21}(Z),\dots,\beta_{2Q}(Z)\ \text{are stable in } Z,$$

concerns a quantity that **does not depend on the time origin**:
$\beta_1 + \beta_2 t = (\beta_1 + c\beta_2) + \beta_2(t - c)$. A properly
defined test of slope instability conditional on nuisance intercepts should
give the same answer either way. That the naive result depends on the time
coding says the naive test is not that test.

The adopted diagnostic partitions the subject-level score as
$U_i = (U_{I,i}, U_{S,i})$ and uses

$$U_{S\cdot I} = U_S - I_{SI} I_{II}^{-1} U_I,$$

projecting the nuisance intercept direction out of the slope score. Under
$t \to t - c$ the design transforms as $F_{\text{new}} = FA$ with
$A = \begin{pmatrix}1 & -c\\ 0 & 1\end{pmatrix}$, so the scores transform as
$(A^\top \otimes I_Q)$:

$$U_{I,\text{new}} = U_I, \qquad U_{S,\text{new}} = U_S - c\,U_I$$

The slope score changes only by a multiple of the intercept score — exactly the
direction the projection removes. **$U_{S\cdot I}$ is therefore invariant to
the time origin by construction**, and that is a falsifiable prediction, tested
in `run/03_time_origin_invariance.R`.

The term used here is *nuisance-adjusted* or *nuisance-orthogonalized*,
deliberately not *efficient score*: GEE rests on estimating equations rather
than a full likelihood, so "efficient" would require further justification
about which empirical covariance or sensitivity matrix is in play. Time-origin
invariance is a purely algebraic property and is testable as it stands.

Implementation notes. The projection is estimated from the **subject-level**
scores within the current node and applied row-wise to the wave-level scores;
because it is linear, summing the adjusted wave-level scores over a person
recovers the adjusted subject-level score, so method A with `cluster = id`
stays coherent with B(a). `parm` is **not** used for this diagnostic — the
$Q$-dimensional adjusted process is handed to `mob()` directly and whitened on
its own covariance.

### The intercept diagnostic

Run on raw time, unadjusted. There is no invariance to demand: $\beta_{1q}$
genuinely depends on where $t = 0$ sits, and "do baseline levels differ" is
well posed once $t = 0$ is fixed at baseline. A symmetric $U_{I\cdot S}$
adjustment is available in principle but nothing in the validation calls for
it. This diagnostic is a supplementary check, not a methodological
contribution.

### Reporting depth

The root split is the primary finding; splits below it are interpreted more
cautiously. In simulation the root variable and root cut were recovered in 100%
of replications, while ~7% of trees carried an additional, by-construction
spurious split below a correct root. This is not a rule that deeper splits are
always spurious — an empirical second split can be both strong and
theoretically meaningful — but a deeper split needs its own justification
rather than inheriting the root split's credibility.

---

## 5. Scope

The implementation supports Gaussian outcomes with the identity link and
`corstr = "independence"`, and raises an error outside that. That is the scope
in which the method was validated; see `VALIDATION.md` §6.
