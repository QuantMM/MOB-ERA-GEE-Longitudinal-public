# MOB-ERA-GEE for multivariate longitudinal outcomes

Code accompanying the article

> Kim, S. *Migration, Social Position, and Multidimensional Physical Activity
> Over Time: Findings from the Canadian Longitudinal Study on Aging.*
> Frontiers in Sociology (to be submitted).

The method combines extended redundancy analysis via generalized estimating
equations (ERA-GEE) with model-based recursive partitioning (MOB) to detect
subgroups that differ in the fitted Baseline levels and the average change per
assessment of several longitudinal outcomes analysed jointly. For the
development and evaluation of the MOB-ERA-GEE framework, see
Sichkaruk (2026), *Advancing extended redundancy analysis via generalized
estimating equations (ERA-GEE): Integrating model-based recursive partitioning
(MOB) for subgroup detection* [Master's thesis, University of Manitoba],
http://hdl.handle.net/1993/39920.

> **Status.** Research code accompanying an applied article, not a
> general-purpose R package. The release cited in the article is **v1.0.1**.

---

## What is in this repository

| Folder | Contents | Needs CLSA data? |
|---|---|---|
| `R/`, `run/`, `run_all.R`, `docs/`, `results/` | Method implementation, implementation checks, and simulation-based validation on seeded synthetic data | **No** |
| `empirical/` | The complete analysis pipeline of the article | **Yes** (not included) |

`docs/METHOD.md` gives the specification and the reasoning behind each design
choice; `docs/VALIDATION.md` is the full validation report, including what the
validation does not establish.

---

## Part 1 — Method and validation (runs without any data)

Every methodological check can be reproduced from seeded synthetic data built
to resemble the empirical design.

```r
install.packages(c("partykit", "MASS", "matrixcalc", "sandwich"))

# from the repository root:
source("run_all.R")                  # everything except the long simulation
FULL <- TRUE; source("run_all.R")    # including the 1,000-replicate simulation
```

| Script | What it checks |
|---|---|
| `run/00_regression_test.R` | The synthetic reference case and fitted results are unchanged |
| `run/01_worked_example.R` | One complete fit on the synthetic reference case |
| `run/02_cluster_equivalence.R` | Person-wave scores aggregated by participant reproduce the participant-level formulation |
| `run/03_time_origin_invariance.R` | The change diagnostic does not depend on where time is coded as zero |
| `run/04_diagnostics_table.R` | Joint, level, and change diagnostics under four planted conditions |
| `run/05_repeated_simulation.R`, `run/05b_summarise_simulation.R` | Calibration, recovery, and specificity over 1,000 replicates per condition |

Main simulation results (3,000 participants, three assessments, four outcomes,
seven candidate partitioning variables, α = .05 with Bonferroni adjustment):

- False-split rate under the null: joint 5.6%, level 6.0%, change 4.3%
  (Monte Carlo SE 0.6–0.8 points).
- Correct root variable and category grouping: 100% in every non-null condition.
- Change diagnostic under heterogeneity in levels only: 4.1% false splits,
  versus 4.3% under the null.
- Further splits below a correct root: 7.1%, spread across the null variables.

The simulation contains no age terms; the empirical model adds them as global
coefficients (see below).

---

## Part 2 — Empirical analysis (CLSA data required, not included)

### Data

The article uses the Comprehensive cohort of the Canadian Longitudinal Study
on Aging (CLSA): Baseline Comprehensive dataset version 7.2, Follow-up 1
Comprehensive dataset version 5.1, and Follow-up 2 Comprehensive dataset
version 4.0. CLSA data cannot be redistributed; access requires an approved
application to the CLSA.

Place the exported CSV files in `empirical/CLSA Datafiles for analysis/`. The
file names in `DataPrep/01_build_integrated_dataset.R` are those of the
author's approved extract and will need to be changed for another extract.

### Running the pipeline

Use `empirical/` as the working directory and run the scripts in order:

```sh
cd empirical
for f in DataPrep/0*.R; do Rscript "$f"; done
for f in Analysis/1[5-9]_*.R Analysis/2[0-3]_*.R; do Rscript "$f"; done
```

| Script | Produces |
|---|---|
| `DataPrep/01`–`04` | Integrated long data set; recoding of the physical activity items and covariates |
| `DataPrep/05` | Descriptive report |
| `DataPrep/06` | Balanced and all-available analytic data sets |
| `DataPrep/07` | Attrition and retention |
| `DataPrep/08` | Immigration Flag and years in Canada |
| `Analysis/15` | Primary multivariate longitudinal model; numerical implementation checks |
| `Analysis/16` | Partitioning trees in the full sample and among immigrants |
| `Analysis/17` | Immigrant–non-immigrant comparison (Table 2) |
| `Analysis/18` | Finer settlement-duration bands |
| `Analysis/19` | Three-group comparison, within-group trees, sex and racialization contrasts |
| `Analysis/20` | Table 1 and Figures 1–3 |
| `Analysis/21` | Flag-determinate refit; between-group comparison of the racialization contrast |
| `Analysis/22` | Trees refitted with α = .01 (Supplementary Table S11) |
| `Analysis/23` | Supplementary Tables S1–S10 with numerical cross-checks and a small-cell guard |

`Analysis/23` writes the tables into the Supplementary Results LaTeX file; set
`PKG` at the top of the script to the folder that holds that file.

**All outputs are written to `empirical/DataPrep/out/`, contain CLSA-derived
results, and must stay local** (the folder is listed in `.gitignore`). Every
output that leaves the pipeline suppresses counts below 6, including counts
that could be recovered by subtraction, as the CLSA Publication Policy
requires.

### Model in brief

For participant *i*, assessment *t* = 0, 1, 2 (Baseline, Follow-up 1,
Follow-up 2), and outcome *q*,

```
Y_itq = β1q + β2q·t + γq·A_i + δq·(t × A_i) + ε_itq,
```

where `A_i` is age at the baseline visit centered at the sample mean. The four
design terms are ERA-GEE components with weights fixed at 1. Estimation uses a
Gaussian family, identity link, and working independence, with a
participant-clustered sandwich covariance for multivariate Wald tests.

In the partitioning trees, `β1q` and `β2q` are node-specific, while `γq` and
`δq` are global and estimated by alternation with the tree (`palm_era_gee.R`).
Three diagnostics test different sets of scores:

- **Joint:** all eight node-specific parameters.
- **Level:** the four levels.
- **Change:** the four change parameters, after projecting out the level scores.

Tree settings: minimum node size 100 participants, maximum depth 4,
test-based pre-pruning at α = .05, no post-pruning.

---

## Software

R 4.6.1 with partykit 1.3.0, sandwich 3.1.3, MASS 7.3.65, and matrixcalc.
See `sessionInfo.txt`.

## Data availability

Data are available from the Canadian Longitudinal Study on Aging
(www.clsa-elcv.ca) for researchers who meet the criteria for access to
de-identified CLSA data.

## Citation

Please cite the article and this repository (see `CITATION.cff`). For the
methodology, please also cite Sichkaruk (2026).

## License

MIT (see `LICENSE`).
