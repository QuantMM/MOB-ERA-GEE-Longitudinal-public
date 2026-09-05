# Validation

Everything below reproduces from a clean R session with no data access. Script
names refer to `run/`.

| § | claim | script |
|---|---|---|
| 1 | the fitted model is the specified one, and recovers a planted answer | `00_regression_test.R`, `01_worked_example.R` |
| 2 | `cluster = id` implements subject-level estimating functions | `02_cluster_equivalence.R` |
| 3 | the slope diagnostic does not depend on the time origin | `03_time_origin_invariance.R` |
| 4 | the three diagnostics behave as intended across four truth conditions | `04_diagnostics_table.R` |
| 5 | calibration, recovery and parameter accuracy over 1,000 replications | `05_repeated_simulation.R` |
| 6 | what has **not** been established | — |

The reference dataset (`R/synthetic_data.R`, seed 20260903) is 3,000 persons ×
3 balanced waves, `Q = 4` continuous outcomes, and seven person-level
partitioning variables whose marginals loosely mirror the CLSA comprehensive
cohort. Only `IMM` (non 59.4%, recent 12.4%, established 28.2%) ever carries an
effect; the other six are null by construction.

---

## 1. Specification and recovery

`00_regression_test.R` asserts, and passes:

```
B is 2 x Q                              PASS
rownames of B = (intercept, slope)      PASS
estfun is n_row x 2Q                    PASS
scores vanish at optimum                PASS   max |colSum| = 1.14e-10
f_const == 1 exactly                    PASS
f_time  == t exactly                    PASS
first split variable is IMM             PASS
one node is exactly {recent}            PASS
other node is exactly {established, non} PASS
node parameters within tolerance         PASS
```

`01_worked_example.R` runs one transparent end-to-end case. Truth: baselines
identical across immigrant groups, recent immigrants declining much faster on
all four outcomes, no other variable having any effect.

```
Parameter instability tests (root):
              IMM       SEX      ETHN      EDU   INCNEED  HOMEOWN    URBAN
statistic  715.66    9.7120   38.2600   9.2153   25.1568   8.0897   9.0385
p.value   4.2e-141   0.9053    0.9962   1.0000    0.9711   1.0000   0.9449

Best splitting variable: IMM
Selected split: recent | non, established
```

| node | n | outcome | est. intercept | true | est. slope | true |
|---|---|---|---|---|---|---|
| 2 (recent) | 371 | walk | 2.838 | 2.95 | −0.431 | −0.45 |
| | | lsport | 1.397 | 1.40 | −0.284 | −0.30 |
| | | msport | 1.125 | 1.20 | −0.224 | −0.25 |
| | | ssport | 1.464 | 1.50 | −0.344 | −0.35 |
| 3 (non + established) | 2629 | walk | 2.946 | 2.95 | −0.060 | −0.06 |
| | | lsport | 1.399 | 1.40 | −0.033 | −0.03 |
| | | msport | 1.186 | 1.20 | −0.019 | −0.02 |
| | | ssport | 1.517 | 1.50 | −0.044 | −0.04 |

The intercepts land near the common truth in **both** nodes, correctly, while
the slopes separate by roughly a factor of seven. Only a model carrying the
intercepts can show that.

---

## 2. Method A reproduces method B(a)

| check | result |
|---|---|
| root coefficients | `max |coef_A − coef_B(a)| = 0.000e+00` |
| root objective | 18558.531963 both, difference `0.000e+00` |
| **score aggregation** | `max |rowsum(estfun_A, by person) − estfun_B(a)| = 0.000e+00` |
| root instability tests | identical to all displayed digits |
| split variable / index | `IMM`, `2,1,2` — both |
| terminal-node membership | agreement 1.0000 (371 / 2629, zero off-diagonal) |
| terminal-node parameters | `max |Δintercept| = 0`, `max |Δslope| = 0` |

**A is exactly equivalent to B(a) for the score aggregation and the clustered
covariance, and reproduces B(a) exactly for the categorical partitioning
variables in this reference case. For numeric person-level partitioning
variables the fluctuation processes are not algebraically identical**, although
the observed discrepancy was negligible: adding a numeric person-level
covariate with no effect gave `max |statistic_A − statistic_B(a)| = 0.0998` and
`max |p_A − p_B(a)| = 6.1e-04`, with the same split variable and the same tree.

That distinction matters because `SDC_DRES_COM` (years since immigration) is
expected to enter the empirical analysis as an integer covariate, so the
exactness claim must stay confined to the categorical case.

**Why `cluster = id` is not merely a convenience.** Dropping it distorts the
root test statistics in both directions, inflating the null partitioning
variables by up to 1.53×:

```
ratio A(no cluster) / A(cluster = id):
    IMM     SEX    ETHN     EDU INCNEED HOMEOWN   URBAN
  0.966   0.898   1.090   1.531   1.350   1.414   0.898
```

`INCNEED`'s p-value moves from 0.9711 to 0.4643 — for a variable with no effect.

---

## 3. Time-origin invariance of the slope diagnostic

Pass criterion, fixed in advance: the adjusted slope instability statistic on
raw $t = (0,1,2)$ and on centred $t = (-1,0,1)$ must agree to numerical
tolerance.

| condition | \|stat(raw) − stat(centred)\| | all-Z max difference |
|---|---|---|
| slope-only | **0.000e+00** | 8.9e-14 |
| intercept-only | **1.73e-13** | 1.7e-13 |

The invariance holds at the level of the score itself, not merely the
statistic: `max |U_S·I(raw) − U_S·I(centred)|` was 2–3e-14 at wave level and
5–7e-14 at subject level. The slopes themselves are origin-free to ~1e-14, as
they must be.

### Why this was needed

Block correlation between the intercept and slope scores, matched outcomes:

| | correlation |
|---|---|
| raw t = 0,1,2 | 0.948 0.944 0.943 0.948 |
| centred t = −1,0,1 | 0.017 −0.023 0.023 0.000 |

Under intercept-only truth, where the slopes are identical by construction:

```
RAW time,     naive parm = slope      split on IMM     stat = 73.92   p = 5.7e-12
CENTRED time, naive parm = slope      split on (none)  stat =  6.70   p = 0.997
```

The joint statistic, by contrast, is exactly invariant (715.66 and 260.87 under
both codings), as a joint test over the whole model space must be.

### Adopted (B) versus the centred-time workaround (A)

| condition | B: adjusted | A: naive on centred t | difference | same conclusion |
|---|---|---|---|---|
| slope-only | 484.43 → IMM | 508.88 → IMM | 5.0% | yes |
| intercept-only | 7.04 → (none) | 6.70 → (none) | 4.8% | yes |
| **intercept + slope** | 337.60 → IMM | 392.21 → IMM | **16.2%** | yes |
| null | 5.91 → (none) | 5.86 → (none) | 0.8% | yes |

The centred-time workaround gets the same substantive answer everywhere, so it
is not a bare hack. But the gap is largest where **both** kinds of
heterogeneity are present — the realistic case — because centring orthogonalises
the design columns, not the empirical scores. The residual block correlation
under centred time is itself condition-dependent (0.017 … 0.000 under
intercept-only truth, but 0.112 … 0.058 under slope-only truth), whereas the
adjusted score's invariance is structural.

---

## 4. Four truth conditions × three diagnostics

Seed, covariates and noise are identical across conditions; only the truth for
the `recent` group changes. Observed = expected in all twelve cells, `IMM`
statistic in parentheses:

| truth | joint | slope | intercept |
|---|---|---|---|
| slope-only | IMM (715.66) | IMM (484.43) | (none) (9.91) |
| intercept-only | IMM (260.87) | (none) (7.04) | IMM (186.94) |
| intercept + slope | IMM (1028.58) | IMM (337.60) | IMM (77.63) |
| null | (none) (20.62) | (none) (5.91) | (none) (11.70) |

The "(none)" cells sit at 5.9–20.6 while every "IMM" cell is 78–1029, so no
verdict is near a boundary. Every split was `recent | established + non` at two
terminal nodes; no split on any of the six null partitioning variables.

---

## 5. Repeated simulation

1,000 replications × 4 conditions, 9,000 trees, **0 errors**, 84 minutes on 32
workers. All percentages carry **Monte Carlo standard errors**, not confidence
intervals.

### 5.1 Calibration under the global null

Nominal α = 0.05, Bonferroni over 7 partitioning variables. Approximately
nominal for all three:

| diagnostic | P(any split) |
|---|---|
| joint | 5.6% (MCSE 0.7) |
| slope | 4.3% (MCSE 0.6) |
| intercept | 6.0% (MCSE 0.8) |

By variable, as % of all replications — no variable dominates:

| | IMM | SEX | ETHN | EDU | INCNEED | HOMEOWN | URBAN |
|---|---|---|---|---|---|---|---|
| joint | 0.7 | 0.6 | 1.3 | 1.1 | 0.6 | 0.4 | 0.9 |
| slope | 0.4 | 0.4 | 0.9 | 0.8 | 0.5 | 0.7 | 0.6 |
| intercept | 0.7 | 0.5 | 0.5 | 1.5 | 0.9 | 0.8 | 1.1 |

### 5.2 Root recovery, non-null conditions

| condition | diagnostic | root variable = IMM | wrong variable | no split | **root cut exact** |
|---|---|---|---|---|---|
| slope-only | joint | 100.0% | 0.0 | 0.0 | **100.0%** |
| slope-only | slope | 100.0% | 0.0 | 0.0 | **100.0%** |
| intercept-only | joint | 100.0% | 0.0 | 0.0 | **100.0%** |
| intercept+slope | joint | 100.0% | 0.0 | 0.0 | **100.0%** |
| intercept+slope | slope | 100.0% | 0.0 | 0.0 | **100.0%** |

In every one of the 5,000 trees fitted in these cells, both the root variable
and the root category partition were correct, the partition being recorded as
`recent | established|non` in 1,000 of 1,000 replications per cell.

This is **ceiling power at a deliberately large planted signal**, not a general
claim about sensitivity. The `intercept-only / slope` cell is absent by design:
the slopes are homogeneous there, so it is a specificity result and appears in
§5.4.

### 5.3 Tree complexity

Conditional on a correct root cut:

| condition | diagnostic | n | exact two-node tree | over-split | mean nodes |
|---|---|---|---|---|---|
| slope-only | joint | 1000 | 92.9% | 7.1% | 2.072 |
| slope-only | slope | 1000 | 92.9% | 7.1% | 2.071 |
| intercept+slope | joint | 1000 | 92.9% | 7.1% | 2.072 |
| intercept+slope | slope | 1000 | 92.9% | 7.1% | 2.071 |

After a correct root split each child is effectively a null node and is then
subjected to a further round of seven instability tests, which is what an
over-split rate of this size predicts.

Depth-2 split variables among the 284 over-splitting replications (286 splits):

| EDU | INCNEED | SEX | ETHN | IMM | URBAN | HOMEOWN |
|---|---|---|---|---|---|---|
| 72 | 50 | 40 | 38 | 36 | 28 | 22 |

All seven appear; the most frequent accounts for 25.2%. There is no reason to
expect an even spread — the seven variables differ in number of levels, in
marginal distribution, and hence in how many candidate binary splits the test
sees. `EDU` has exactly zero effect in the data-generating process, so its
frequency cannot be a real secondary signal. The reportable point is
qualitative: secondary splits spread across all seven null variables with no
variable close to dominating.

### 5.4 Slope diagnostic specificity — the key result

Under intercept-only truth the slopes are homogeneous by construction, so for
the slope diagnostic that condition is itself a null and any split is a false
positive.

| | P(any false split) |
|---|---|
| pure null | **4.3%** (MCSE 0.6) |
| intercept-only, strong nuisance intercept heterogeneity present | **4.1%** (MCSE 0.6) |

of which 0.2% landed on IMM and 3.9% elsewhere.

> **Intercept heterogeneity does not inflate the slope diagnostic's
> false-positive rate: 4.1% versus 4.3% — essentially indistinguishable at the
> Monte Carlo resolution of 1,000 replications.**

For contrast, the joint diagnostic splits on IMM in 100.0% of that same
condition, correctly — the trajectories genuinely differ, in baseline. And a
naive raw-time slope test produced a statistic of 73.9 (p = 5.7e-12) there,
which is what motivated the adjusted score.

### 5.5 Parameter recovery

Coefficients are refitted on the **root's two branches**, so they are defined
whenever the root cut is correct, over-splitting or not. **Absolute RMSE is the
primary quantity**: the true slopes span two orders of magnitude across
conditions (−0.45 down to −0.02), so a ratio explodes wherever the truth is
near zero.

| parameter | conditioning | bias (range) | **absolute RMSE** |
|---|---|---|---|
| β₁ recent | correct root cut (n = 1000) | −0.0012 … 0.0000 | 0.0364 … 0.0372 |
| β₁ recent | exact two-node tree (n = 929) | −0.0008 … +0.0002 | 0.0365 … 0.0375 |
| β₂ recent | correct root cut | 0.0000 … +0.0009 | 0.0162 … 0.0177 |
| β₂ recent | exact two-node tree | −0.0001 … +0.0008 | 0.0164 … 0.0178 |
| β₁ others | correct root cut | −0.0009 … +0.0003 | 0.0130 … 0.0140 |
| β₂ others | correct root cut | −0.0002 … 0.0000 | 0.0062 … 0.0064 |

Essentially unbiased. Absolute slope RMSE was ~0.016–0.018 in the smaller
371-person node and ~0.006 in the larger node. For the deliberately steep
recent-immigrant slopes in the slope-heterogeneity conditions
(−0.45, −0.30, −0.25, −0.35) that is roughly 4–7% of the true magnitude; under
intercept-only truth the recent slopes are −0.06 to −0.02 and the same absolute
RMSE is a large fraction of them, which is why the ratio is not used as a
summary.

The two conditionings agree to the third decimal: over-splitting below the root
does not damage root-level parameter recovery.

The four conditions share replication seeds (common random numbers) and differ
only by a deterministic shift of the `recent` group's truth, so bias and RMSE
come out identical across conditions. That is expected — it removes Monte Carlo
variation between conditions — not a bug.

---

## 6. What has *not* been established

One reference design:

- **balanced** repeated measurements — every person contributes all three waves
- one 3-level true partitioning variable, affected group at 12% of the sample
- a single, deliberately large effect size
- Gaussian outcomes on the raw scale
- working independence across the four outcomes
- n = 3,000 persons

Unbalanced or missing-wave designs, smaller effects, other cohort compositions,
and non-independence working correlation structures are **not** covered.
`R/era_gee_fixed.R` raises an error outside the validated family and working
correlation rather than approximating silently.

---

## Summary

> These simulations provide strong evidence that the implementation is
> correctly calibrated and recovers the intended heterogeneity under the
> reference design used for the empirical application.

The clearest structural finding is the separation between the two levels of the
tree: **root heterogeneity detection was highly stable, while deeper recursive
splits carry the usual multiple-testing and tree-growth variability.** These
results support treating the root split as the primary finding and interpreting
deeper splits more cautiously. They do not establish that every deeper split is
spurious — but a deeper split needs its own justification rather than
inheriting the root's credibility.
