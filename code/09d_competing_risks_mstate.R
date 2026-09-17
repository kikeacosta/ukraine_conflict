# ==============================================================================
# STEP 09D - COMPETING RISKS MULTI-STATE MODEL
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Replaces the separate per-outcome Cox models in Phase 2 with a single
# competing-risks multi-state model. The separate-model approach (Phase 2)
# fits P(dead), P(alive), P(prisoner) independently, which does NOT
# guarantee they sum to 1 and does not properly share information across
# outcomes (someone who transitions to "alive" is removed from risk of
# "dead" at that instant - competing risks, not independent events).
#
# STATE STRUCTURE
# ----------------
#   State 1: missing        (starting state, at risk)
#   State 2: dead            (absorbing)
#   State 3: alive           (absorbing)
#   State 4: prisoner        (absorbing)
#
# Transitions: 1->2, 1->3, 1->4 (no transitions out of 2/3/4 are modeled;
# see LIMITATIONS at the bottom for why "resolved" states are treated as
# absorbing here even though real people occasionally get reclassified).
#
# METHOD
# ------
# Uses mstate::msprep() to build the long-format transition dataset, fits
# cause-specific Cox models stratified by transition type (this is the
# standard "cause-specific hazards" competing risks approach), then derives
# cumulative incidence functions (CIFs) with mstate::probtrans(). CIFs are
# the competing-risks-correct version of "P(dead by time t)" - unlike a
# naive Kaplan-Meier complement, they properly account for the fact that
# people who die can no longer become "alive" or "prisoner", and vice versa.
#
# REQUIRES
# --------
# - Phase 1 output: data_inter/ukr_ualosses_pooled_transitions.rds
# - Phase 2 output: data_inter/ukr_ualosses_rates_by_duration.rds (for
#   comparison against the separate-model estimates)
# - R package: mstate
#
# OUTPUT
# ------
# - data_inter/ukr_ualosses_mstate_cifs.rds       (cumulative incidence by
#                                                   cohort and time)
# - data_inter/ukr_ualosses_mstate_vs_phase2.rds  (side-by-side comparison)
# - figures/exploratory/mstate_cif_by_cohort.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

if (!require("mstate", quietly = TRUE)) {
  install.packages("mstate")
  library(mstate)
}

if (!file.exists("data_inter/ukr_ualosses_pooled_transitions.rds")) {
  stop(
    "Phase 1 output not found. Run 09_ualosses_missing_analysis_improved.R first.",
    call. = FALSE
  )
}

pooled_transitions <- read_rds("data_inter/ukr_ualosses_pooled_transitions.rds")

# ==============================================================================
# 1. BUILD PER-PERSON TIME-TO-EVENT DATA
# ==============================================================================
# Same duration construction as Phase 2 (09b): time = register date of the
# window's endpoint minus the disappearance date. Kept here rather than
# reading Phase 2's output directly because msprep() needs person-level rows
# with a single time-to-event per transition type, not the pre-aggregated
# duration bands 09b produces.

surv_data <- pooled_transitions %>%
  filter(status == "missing") %>%
  mutate(
    # Actual register export dates, from the filenames themselves
    # (250916, 251204, 260530, 260919), not estimates. Three windows now
    # that v16 (Dec 2025) sits between v14 and v18.
    register_date = case_when(
      window == "v14→v16" ~ as.Date("2025-12-04"),
      window == "v16→v18" ~ as.Date("2026-05-30"),
      window == "v18→v19" ~ as.Date("2026-09-19")
    ),
    time = as.numeric(register_date - date_evnt),
    cohort = factor(year)
  ) %>%
  filter(time > 0) %>%
  # one row per person per window; status2 is the outcome at that window's
  # endpoint. "prisoner" is folded into "alive" here (not dead, which is all
  # that matters for the mortality imputation this feeds - steps 10/11 only
  # need P(dead) vs P(not dead)). This is also a numerical necessity, not
  # just a simplification: prisoner-status RECORDING jumped sharply between
  # windows in the raw data (e.g. 2025 cohort: 12 events in v14->v18 vs 840
  # in v18->v19 - a recording-practice change, not a hazard pattern), which
  # produced severe (quasi-)complete separation when "prisoner" was fit as
  # its own transition with cohort-specific coefficients (coefficients blew
  # up to ~1e13 on the hazard scale). A 3-way split isn't needed for the
  # question this model answers; 2-way removes the unstable rare-cell terms
  # entirely rather than papering over them with penalization.
  transmute(
    id = row_number(),
    cohort,
    time,
    to_dead  = as.integer(status2 == "dead"),
    to_alive = as.integer(status2 %in% c("alive", "prisoner"))
    # censored (status2 == "missing") implicitly has both = 0
  )

message(sprintf("Competing risks dataset: %d person-windows", nrow(surv_data)))

# ==============================================================================
# 2. PREPARE MSTATE LONG-FORMAT DATA
# ==============================================================================

# trans.comprisk() is purpose-built for exactly this structure (one starting
# state, K competing absorbing outcomes, no intermediate states) rather than
# a hand-built transMat(). It also fixes state naming to "eventfree" for the
# start state, which downstream code must match.
tmat <- trans.comprisk(2, names = c("dead", "alive"))
print(tmat)

# msprep()'s wide format requires time/status to have ncol == nrow(trans) ==
# number of STATES (3: eventfree/missing + 2 outcomes), not ncol == number of
# transitions (2). The first column is the starting state itself; its status
# is structurally NA (msprep requires this - see its dimension check against
# dim(trans)[1]). Each person then gets exactly 2 output rows in the long
# format (one per possible transition out of "missing"), with status = 1 on
# whichever row matches their actual outcome, or both 0 if censored.
ms_long <- msprep(
  time = matrix(rep(surv_data$time, 3), ncol = 3),
  status = cbind(NA, as.matrix(surv_data[, c("to_dead", "to_alive")])),
  trans = tmat,
  # msprep()'s internal keep[ord1] uses base-R single-bracket column
  # semantics that tibble's stricter [.tbl_df rejects for long recycled
  # index vectors (confirmed with a minimal repro) - force a base data.frame.
  keep = as.data.frame(surv_data["cohort"]),
  id = surv_data$id
)

ms_long <- expand.covs(ms_long, "cohort", longnames = FALSE)

message(sprintf("mstate long-format rows: %d (one row per person per possible transition)",
                nrow(ms_long)))

# ==============================================================================
# 3. FIT CAUSE-SPECIFIC HAZARDS MODEL, STRATIFIED BY TRANSITION
# ==============================================================================
# Standard competing-risks approach: one Cox model, stratified so each
# transition (missing->dead, missing->alive, missing->prisoner) gets its own
# baseline hazard, but cohort effects are estimated within each stratum.
# This is what makes it "competing risks" rather than three unrelated Cox
# fits: the risk set and event indicator now correctly reflect that a person
# who becomes "alive" at time t leaves the risk set for "dead" at time t.

cohort_dummy_cols <- grep("^cohort", names(ms_long), value = TRUE)
cohort_dummy_cols <- setdiff(cohort_dummy_cols, "cohort")

fit_formula <- as.formula(
  paste0(
    "Surv(Tstart, Tstop, status) ~ ",
    paste(cohort_dummy_cols, collapse = " + "),
    " + strata(trans)"
  )
)

cr_model <- coxph(fit_formula, data = ms_long, method = "breslow")

message("\n=== COMPETING RISKS MODEL SUMMARY ===\n")
print(summary(cr_model))

write_rds(cr_model, "data_inter/ukr_ualosses_mstate_model.rds")

# ==============================================================================
# 4. CUMULATIVE INCIDENCE FUNCTIONS BY COHORT
# ==============================================================================
# This is the payoff: CIFs give P(dead by time t) and P(not dead by time t)
# that are internally consistent (they + P(still missing) sum to 1 at every
# t), unlike the separate Phase 2 Cox models.

# expand.covs(..., longnames = FALSE) names transition-specific dummies
# "cohort<levelIndex>.<transIndex>", where levelIndex is the cohort factor's
# level position EXCLUDING the reference level (level 1 = 2022 gets no
# dummy at all), NOT the literal year - confirmed against a toy dataset
# before touching this, since "cohort2023.1" (my first, wrong guess) is not
# what expand.covs() actually produces.
cohort_levels <- levels(surv_data$cohort)

cif_by_cohort <- function(cohort_year) {
  # One row per transition (1=dead, 2=alive); msfit/probtrans need both to
  # compute the joint transition-probability trajectory, not one row
  # repeated - a single shared row (an earlier version of this function used
  # data.frame(..., nrow = 1)) does not vary `trans` per row.
  newdata <- as.data.frame(matrix(0, nrow = 2, ncol = length(cohort_dummy_cols)))
  names(newdata) <- cohort_dummy_cols
  newdata$trans <- 1:2
  newdata$strata <- newdata$trans

  lvl_idx <- match(as.character(cohort_year), cohort_levels)
  if (!is.na(lvl_idx) && lvl_idx > 1) {
    for (k in 1:2) {
      col <- sprintf("cohort%d.%d", lvl_idx - 1, k)
      if (col %in% names(newdata)) newdata[newdata$trans == k, col] <- 1
    }
  }

  msf_c <- msfit(cr_model, newdata = newdata, trans = tmat)
  pt <- probtrans(msf_c, predt = 0)[[1]]

  # probtrans() evaluates every cohort's CIF over the UNION of observed
  # event/censoring times pooled across ALL cohorts (up to ~1500+ days, set
  # by 2022's long follow-up), not just that cohort's own range. For 2025,
  # whose actual maximum possible follow-up (disappearance to the Sep 2026
  # register) is under ~9 months, this means the model is extrapolated 4+
  # years past any data it could ever have, using a Cox PH assumption that
  # has no business holding that far out - and its (numerically extreme:
  # HR ~150,000x baseline for the alive transition) coefficient overflows
  # cleanly to Inf/NaN out there. The fix is not numerical - it's to never
  # report a cohort's CIF beyond its own observable horizon in the first
  # place, confirmed here rather than assumed.
  max_follow_up <- max(surv_data$time[surv_data$cohort == cohort_year])
  pt <- pt[pt$time <= max_follow_up, ]

  if (any(!is.finite(as.matrix(pt[, -1])))) {
    warning(sprintf(
      "Non-finite CIF values remain for cohort %s even within its observed follow-up range (max %.0f days) - inspect the Cox fit for this cohort before trusting its imputation.",
      cohort_year, max_follow_up
    ))
  }

  # 3 states now (eventfree/missing, dead, alive) -> pstate1/pstate2/pstate3.
  # No pstate4: "prisoner" was folded into "alive" upstream (section 1).
  as_tibble(pt) %>%
    transmute(
      cohort = cohort_year,
      time_days = time,
      time_months = time / 30.44,
      p_missing = pstate1,
      p_dead = pstate2,
      p_alive = pstate3,
      # 95% CI on the dead CIF via the delta-method SE mstate provides
      p_dead_se = se2,
      p_dead_lo = pmax(0, p_dead - 1.96 * p_dead_se),
      p_dead_hi = pmin(1, p_dead + 1.96 * p_dead_se)
    )
}

cifs <- map_dfr(levels(surv_data$cohort), cif_by_cohort)

message("\n=== CIF CHECK: probabilities sum to 1 at each timepoint ===")
check_sum <- cifs %>%
  mutate(total = p_missing + p_dead + p_alive) %>%
  summarise(max_deviation = max(abs(total - 1)))
print(check_sum)
stopifnot(check_sum$max_deviation < 1e-6)
message("✓ PASS: CIFs are internally consistent (sum to 1 within floating-point tolerance)")

write_rds(cifs, "data_inter/ukr_ualosses_mstate_cifs.rds")

# ==============================================================================
# 5. COMPARE TO PHASE 2 (SEPARATE COX MODELS)
# ==============================================================================

if (file.exists("data_inter/ukr_ualosses_rates_by_duration.rds")) {
  rates_by_duration <- read_rds("data_inter/ukr_ualosses_rates_by_duration.rds")

  duration_marks <- tibble(
    duration_band = c("0–5 mo", "6–11 mo", "12–23 mo", "24+ mo"),
    time_months = c(2.5, 8.5, 17.5, 30)
  )

  cif_at_marks <- cifs %>%
    inner_join(
      duration_marks %>%
        crossing(cohort = unique(cifs$cohort)) %>%
        mutate(cohort = as.character(cohort)),
      by = "cohort",
      relationship = "many-to-many"
    ) %>%
    mutate(time_diff = abs(time_months.x - time_months.y)) %>%
    slice_min(time_diff, by = c(cohort, duration_band), n = 1) %>%
    select(cohort, duration_band, p_dead_mstate = p_dead,
           p_dead_lo_mstate = p_dead_lo, p_dead_hi_mstate = p_dead_hi)

  comparison <- rates_by_duration %>%
    mutate(cohort = as.character(cohort)) %>%
    select(cohort, duration_band, p_dead_phase2 = p_dead,
           ci_lo_phase2 = ci_dead_lo, ci_hi_phase2 = ci_dead_hi) %>%
    left_join(cif_at_marks, by = c("cohort", "duration_band")) %>%
    mutate(
      ci_width_phase2 = ci_hi_phase2 - ci_lo_phase2,
      ci_width_mstate = p_dead_hi_mstate - p_dead_lo_mstate,
      pct_narrower = (ci_width_phase2 - ci_width_mstate) / ci_width_phase2 * 100
    )

  message("\n=== PHASE 2 (separate Cox) vs COMPETING RISKS MODEL ===\n")
  print(comparison %>%
        select(cohort, duration_band, p_dead_phase2, p_dead_mstate,
               ci_width_phase2, ci_width_mstate, pct_narrower))

  message(sprintf("\nMean CI narrowing from competing risks model: %.1f%%",
                  mean(comparison$pct_narrower, na.rm = TRUE)))

  write_rds(comparison, "data_inter/ukr_ualosses_mstate_vs_phase2.rds")
} else {
  message("Phase 2 output not found; skipping comparison (run 09b first for the full picture)")
}

# ==============================================================================
# 6. PLOT: CIF BY COHORT
# ==============================================================================

p <- cifs %>%
  ggplot(aes(x = time_months, y = p_dead, color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = p_dead_lo, ymax = p_dead_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "Time Since Disappearance (months)",
    y = "Cumulative Incidence: P(Confirmed Dead)",
    title = "Competing Risks Model: Cumulative Incidence of Death by Cohort",
    subtitle = "Shaded band: 95% CI. Unlike separate Cox models, this CIF properly\naccounts for the competing 'not dead' (alive/prisoner) outcome removing people from risk.",
    color = "Cohort", fill = "Cohort"
  ) +
  theme_bw() + theme(legend.position = "top")

ggsave("figures/exploratory/mstate_cif_by_cohort.png", p, w = 8, h = 5.5)
message("\n✓ Saved: figures/exploratory/mstate_cif_by_cohort.png")

# ==============================================================================
# LIMITATIONS OF THIS IMPLEMENTATION
# ==============================================================================
# 1. "Prisoner" is folded into "alive" (section 1), not modeled as its own
#    outcome. A first version fit 3 competing outcomes (dead/alive/prisoner)
#    and hit severe (quasi-)complete separation on the prisoner transition:
#    cohort-specific coefficients diverged to ~1e13 on the hazard scale.
#    Diagnosis: raw prisoner-status EVENT COUNTS jump sharply between windows
#    (e.g. 2025 cohort: 12 in v14->v18 vs 840 in v18->v19), consistent with a
#    change in recording practice around v19 rather than a real hazard shift;
#    that timing clustering, not a literal zero-event cell, is what broke the
#    fit. Since this model only needs to feed a dead/not-dead imputation into
#    steps 10-11, collapsing prisoner into alive removes the unstable term
#    rather than papering over it with penalization, and is substantively
#    fine for the question asked here - but it does mean this file no longer
#    distinguishes "confirmed POW" from "otherwise resolved alive", which the
#    Phase 1 pooled-rate tables (not model-based) still do descriptively.
# 2. Absorbing-state assumption: real registers occasionally reclassify
#    someone from "dead" back to "missing" (data correction). This model
#    treats dead/alive as absorbing once reached.
# 3. Cause-specific hazards (this approach) answer "instantaneous risk of
#    each outcome among those still missing" and their CIFs correctly sum
#    to 1. A Fine-Gray subdistribution hazard model answers a different
#    question (effect of covariates on the CIF directly) and would be a
#    further refinement if the cohort covariate's effect on the CIF *shape*
#    (not just level) is of interest - not implemented here.
# 4. Two-window pooling: like Phase 1/2, this treats v14->v18 and v18->v19
#    observations as exchangeable id-level records rather than a proper
#    left-truncated recurring-window design. Person-windows from the same
#    individual (if observed missing at both v14 and v18) are treated as
#    independent rows, which understates correlation slightly.
# ==============================================================================

message("\n✓ Competing risks model complete.")
message("  Compare figures/exploratory/mstate_cif_by_cohort.png against")
message("  figures/exploratory/phase2_duration_resolution_patterns.png")
message("  to see the CI-narrowing effect directly.")
