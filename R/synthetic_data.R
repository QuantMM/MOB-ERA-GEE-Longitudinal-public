# ---------------------------------------------------------------------------
# synthetic_data.R
#
# THE REFERENCE CASE.
#
# One seeded synthetic dataset resembling the intended CLSA application, used
# by every script in this folder so that all comparisons are made on exactly
# the same data. Do not change the generator without re-running
# run/00_regression_test.R.
#
# TRUTH (deliberately simple and obvious):
#   * baseline levels IDENTICAL across immigrant groups
#   * recent immigrants decline substantially faster on ALL FOUR outcomes
#   * SEX, ETHN, EDU, INCNEED, HOMEOWN, URBAN generate NO heterogeneity
#   -> MOB should recover IMM, splitting "recent" from the rest, and it should
#      do so because of SLOPE instability, not intercept instability.
#
# This is, informally, the "slope-only" condition of the four planned
# simulation conditions.
# ---------------------------------------------------------------------------

REF_YN <- c("walk_freq", "lsport_freq", "msport_freq", "ssport_freq")
REF_ZN <- c("IMM", "SEX", "ETHN", "EDU", "INCNEED", "HOMEOWN", "URBAN")

REF_TRUTH <- list(
  b1        = c(walk_freq = 2.95, lsport_freq = 1.40, msport_freq = 1.20, ssport_freq = 1.50),
  b2_others = c(walk_freq = -0.06, lsport_freq = -0.03, msport_freq = -0.02, ssport_freq = -0.04),
  b2_recent = c(walk_freq = -0.45, lsport_freq = -0.30, msport_freq = -0.25, ssport_freq = -0.35)
)

make_reference_data <- function(n_person = 3000, n_wave = 3, seed = 20260903,
                                b1        = REF_TRUTH$b1,
                                b2_others = REF_TRUTH$b2_others,
                                b2_recent = REF_TRUTH$b2_recent,
                                b1_recent = NULL) {

  # b1_recent lets later conditions (intercept-only, intercept+slope) move the
  # baseline of the "recent" group. NULL means "same baseline as everyone else".
  if (is.null(b1_recent)) b1_recent <- b1

  set.seed(seed)
  Q  <- length(REF_YN)
  YN <- REF_YN

  pick <- function(lv, p) factor(sample(lv, n_person, TRUE, p), levels = lv)
  person <- data.frame(
    id  = seq_len(n_person),
    IMM = pick(c("non", "recent", "established"),           c(.60, .12, .28)),
    SEX = pick(c("Female", "Male"),                         c(.53, .47)),
    ETHN = pick(c("White","Black","EastAsian","SouthSEAsian","MiddleEastern","Multiple"),
                c(.78, .03, .06, .07, .03, .03)),
    EDU = pick(c("LessThanSecondary","SecondaryOrSomePost","PostSecDegree"), c(.14, .30, .56)),
    INCNEED = pick(c("VeryWell","Adequately","SomeDifficulty","NotVeryWell"), c(.38, .42, .15, .05)),
    HOMEOWN = pick(c("Own","Rent","Other"),                 c(.76, .19, .05)),
    URBAN   = pick(c("Urban","Rural"),                      c(.82, .18))
  )

  # Person-specific deviations, constant across waves and correlated across
  # outcomes, induce within-person correlation. This is a property of the
  # DATA-GENERATING PROCESS, not a random effect in the model: the ERA-GEE
  # base model has no random intercept.
  cs  <- function(rho, sd) sd^2 * (diag(Q) * (1 - rho) + rho)
  U   <- MASS::mvrnorm(n_person, rep(0, Q), cs(0.35, 0.55))
  idx <- rep(seq_len(n_person), each = n_wave)
  tt  <- rep(0:(n_wave - 1), times = n_person)
  E   <- MASS::mvrnorm(n_person * n_wave, rep(0, Q), cs(0.20, 0.45))

  is_rec <- person$IMM[idx] == "recent"
  Y <- sapply(seq_len(Q), function(q)
    ifelse(is_rec, b1_recent[q], b1[q]) +
    ifelse(is_rec, b2_recent[q], b2_others[q]) * tt +
    U[idx, q] + E[, q])
  colnames(Y) <- YN

  dat <- cbind(data.frame(id = idx, wave = tt, X_time = tt), Y,
               person[idx, setdiff(names(person), "id")])
  rownames(dat) <- NULL
  dat$rowid <- seq_len(nrow(dat))          # person-WAVE row index (method A)

  # Person-level frame for the subject-level interface (method B(a)).
  # The Q outcome columns are placeholders that satisfy mob()'s model frame;
  # the fit function never reads them, it looks the person up in `dat`.
  pdat <- person
  pdat$prow <- seq_len(n_person)           # PERSON row index (method B(a))
  ph <- dat[dat$wave == 0, YN, drop = FALSE]   # baseline values as placeholders
  rownames(ph) <- NULL
  pdat <- cbind(pdat, ph)

  list(dat = dat, pdat = pdat, person = person,
       YN = YN, ZN = REF_ZN, Q = Q, n_wave = n_wave, n_person = n_person,
       truth = list(b1 = b1, b1_recent = b1_recent,
                    b2_others = b2_others, b2_recent = b2_recent))
}
