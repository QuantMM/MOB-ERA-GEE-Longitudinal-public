# MOB-ERA-GEE for multivariate longitudinal outcomes

Code and validation for the method used in:

> *Modeling Intersectional Heterogeneity in Multidimensional Physical Activity
> Trajectories Across Immigrant Groups: Evidence from the Canadian Longitudinal
> Study on Aging.* Sunmee Kim, University of Manitoba.

The paper is an applied study and does not develop the methodology in its text.
This repository is where the method is specified, validated, and made
reproducible.

---

## Everything here runs without any data

The empirical analysis uses the Canadian Longitudinal Study on Aging (CLSA),
which cannot be redistributed: access requires an approved application, and the
CLSA Access Agreement does not permit sharing data beyond the research team.

**No CLSA data is needed to run anything in this repository.** Every
methodological claim is checked against a seeded synthetic dataset built to
resemble the empirical design, so any reader can reproduce all of it, exactly,
from a clean R session.

```r
install.packages(c("partykit", "MASS", "matrixcalc"))
# from the repository root:
source("run_all.R")          # everything except the 84-minute simulation
source("run_all.R"); FULL <- TRUE   # or: Rscript run_all.R --full
```

---

## The model

Each person contributes three waves; each wave carries `Q = 4` continuous
outcomes. The base model is a population-averaged multivariate growth model:

```
g(mu_itq) = beta_1q + beta_2q * t          q = 1, ..., Q
```

with **fixed components** `f_const = 1` and `f_time = t`, both component
weights **constrained to 1**, giving a `2 x Q` coefficient matrix:

| | walk | lsport | msport | ssport |
|---|---|---|---|---|
| **beta_1q** | baseline expected outcome (t = 0) | | | |
| **beta_2q** | expected change per wave | | | |

Time is the **raw** wave index 0, 1, 2, so both rows are directly
interpretable. `R/era_gee_fixed.R` asserts `max|f_const - 1| = 0` and
`max|f_time - t| = 0` on every call — these are constraints, not estimates.

Rows are person-waves. The GEE working covariance is `Q x Q` and models
association among the four outcomes *within* a person-wave. Wave-to-wave
correlation is not parameterised; the fact that repeated observations belong to
the same person is reflected in the cluster-level estimating functions and in
robust instability inference, via `mob(..., cluster = id)`.

## The three diagnostics

All three fit the **same** `2 x Q` model on raw time. They differ only in which
parameters enter the instability test.

| tree | question | test |
|---|---|---|
| **joint** | *Which subgroups have different trajectories?* | all `2Q` parameters |
| **slope** | *Which subgroups change at different rates?* | `Q` nuisance-adjusted slope scores `U_S.I = U_S − I_SI I_II^-1 U_I` |
| **intercept** | *Which subgroups start at different levels?* | the `Q` intercepts |

The slope diagnostic orthogonalises the slope score against the nuisance
intercept score. This matters: with a raw time origin the two score blocks are
correlated at about **0.95**, so simply selecting the slope columns does not
test slope heterogeneity — a pure baseline difference leaks into it. The
adjusted score is **exactly invariant to where time is coded as zero**
(verified to 1e-13), which is the property a test of slope heterogeneity should
have, since `beta_1 + beta_2*t = (beta_1 + c*beta_2) + beta_2*(t − c)` leaves
the slope unchanged.

---

## What each script does, and what it produces

| script | claim it supports | runtime |
|---|---|---|
| `run/00_regression_test.R` | the fitted model is the specified one, and the reference case recovers the planted answer | 15 s |
| `run/01_worked_example.R` | a transparent end-to-end example, with a trajectory figure | 30 s |
| `run/02_cluster_equivalence.R` | `cluster = id` implements subject-level estimating functions | 30 s |
| `run/03_time_origin_invariance.R` | the slope diagnostic does not depend on the time origin | 60 s |
| `run/04_diagnostics_table.R` | all three diagnostics behave as intended across four truth conditions | 3 min |
| `run/05_repeated_simulation.R` | calibration, recovery and parameter accuracy over 1,000 replications | 84 min (32 cores) |
| `run/05b_summarise_simulation.R` | summarises the above — runs in seconds against the shipped result file | 5 s |

`results/simulation_R1000.rds` is included, so the simulation tables can be
reproduced without re-running the 84-minute job:

```r
o <- readRDS("results/simulation_R1000.rds"); res <- o$res; N_REP <- o$n_rep
source("run/05b_summarise_simulation.R")
```

### `run/01_worked_example.R`

3,000 persons × 3 waves. Truth: baselines identical across immigrant groups,
recent immigrants declining much faster on all four outcomes, no other variable
having any effect.

```
Parameter instability tests (root):
              IMM       SEX      ETHN      EDU   INCNEED  HOMEOWN    URBAN
statistic  715.66    9.7120   38.2600   9.2153   25.1568   8.0897   9.0385
p.value   4.2e-141   0.9053    0.9962   1.0000    0.9711   1.0000   0.9449

Best splitting variable: IMM
Selected split: recent | non, established
```

| node | n | group | outcome | est. intercept | true | est. slope | true |
|---|---|---|---|---|---|---|---|
| 2 | 371 | recent | walk | 2.838 | 2.95 | **−0.431** | −0.45 |
| 2 | | | lsport | 1.397 | 1.40 | **−0.284** | −0.30 |
| 2 | | | msport | 1.125 | 1.20 | **−0.224** | −0.25 |
| 2 | | | ssport | 1.464 | 1.50 | **−0.344** | −0.35 |
| 3 | 2629 | non + established | walk | 2.946 | 2.95 | **−0.060** | −0.06 |
| 3 | | | lsport | 1.399 | 1.40 | **−0.033** | −0.03 |
| 3 | | | msport | 1.186 | 1.20 | **−0.019** | −0.02 |
| 3 | | | ssport | 1.517 | 1.50 | **−0.044** | −0.04 |

![fitted vs true trajectories](results/fig_synthetic_trajectories.png)

The intercepts land near 2.95 / 1.40 / 1.20 / 1.50 in **both** nodes, which is
the right answer — the truth had no baseline difference. The slopes separate by
roughly a factor of seven. This is only visible because the intercepts are in
the parameter vector: a slope-only model could not have shown that the
baselines agree, and any baseline difference that did exist would have been
absorbed into the slopes.

### `run/04_diagnostics_table.R`

Split variable, observed = expected in all twelve cells. `IMM` statistic in
parentheses.

| truth | joint | slope | intercept |
|---|---|---|---|
| slope-only | IMM (715.66) | IMM (484.43) | (none) (9.91) |
| intercept-only | IMM (260.87) | (none) (7.04) | IMM (186.94) |
| intercept + slope | IMM (1028.58) | IMM (337.60) | IMM (77.63) |
| null | (none) (20.62) | (none) (5.91) | (none) (11.70) |

Each diagnostic fires when, and only when, the kind of heterogeneity it asks
about is present.

### `run/05_repeated_simulation.R` — headline results

1,000 replications × 4 conditions, 9,000 trees, 0 errors. Percentages carry
**Monte Carlo standard errors**, not confidence intervals.

**Calibration under the global null** (nominal α = 0.05, Bonferroni over 7
partitioning variables) — approximately nominal for all three:

| joint | slope | intercept |
|---|---|---|
| 5.6% (MCSE 0.7) | 4.3% (MCSE 0.6) | 6.0% (MCSE 0.8) |

**Root recovery**, non-null conditions — the root variable *and* the root
category cut were correct in **100%** of replications in every cell (5,000
trees). Note this is ceiling power at a deliberately large planted signal, not
a general claim about sensitivity.

**Slope diagnostic specificity** — the key result. Under intercept-only truth
the slopes are homogeneous by construction, so any split is a false positive:

| | P(any false split) |
|---|---|
| pure null | 4.3% (MCSE 0.6) |
| intercept-only, strong nuisance intercept heterogeneity present | **4.1%** (MCSE 0.6) |

Intercept heterogeneity does not inflate the slope diagnostic's false-positive
rate — the two are essentially indistinguishable at this Monte Carlo
resolution. For contrast, a naive slope test on raw time produced a statistic
of 73.9 (p = 5.7e-12) in exactly that condition.

**Parameter recovery**, conditional on a correct root cut — essentially
unbiased, with absolute RMSE of 0.016–0.018 for slopes in the smaller
371-person node and 0.006 in the larger node.

**Tree depth.** 7.1% of trees carried an extra split below a correct root,
spread across all seven null partitioning variables. Root heterogeneity
detection was highly stable; deeper recursive splits carry the usual
multiple-testing and tree-growth variability. These results support treating
the root split as the primary finding and interpreting deeper splits more
cautiously — not a universal rule that deeper splits are spurious, but a reason
that a deeper split needs its own justification rather than inheriting the
root's credibility.

---

## What has *not* been established

The validation above is for one reference design:

- **balanced** repeated measurements — every person contributes all three waves
- one 3-level true partitioning variable, with the affected group at 12% of the sample
- a single, deliberately large effect size
- Gaussian outcomes on the raw scale
- working independence across the four outcomes
- n = 3,000 persons

Unbalanced or missing-wave designs, smaller effects, other cohort compositions,
and non-independence working correlation structures are **not** covered.
`R/era_gee_fixed.R` rejects families and working correlation structures outside
the validated scope with an error rather than approximating them silently.

---

## Layout

```
R/                            the method
  era_gee_fixed.R             fixed-component ERA-GEE fitter + two mob() fit factories
  nuisance_adjusted_slope.R   the slope diagnostic
  synthetic_data.R            the seeded reference-case generator
run/                          one script per claim, in order
results/                      figure, and the 1,000-replication result object
docs/
  METHOD.md                   specification and design rationale
  VALIDATION.md               full validation report
run_all.R
```

## Method background

The MOB-ERA-GEE framework was developed in Samantha Sichkaruk's Master's
thesis (University of Manitoba, MSpace) and presented at ISDSA, where it
received the best-paper award in that track. A methodological manuscript is in
preparation. This repository implements the fixed-component longitudinal
specialisation used in the applied paper above, together with its validation;
it is self-contained and does not depend on that codebase.

## Environment

R 4.6.1, `partykit` 1.3.0, `MASS`, `matrixcalc`. See `sessionInfo.txt`.

`partykit`'s version matters. Two behaviours the results rely on are internal
to it: how `cluster` enters the covariance of the instability test, and that
`parm` selects columns of the process *after* whitening. Both were verified
against 1.3.0 and are documented in `docs/METHOD.md`.

## License

MIT. See `LICENSE`.
