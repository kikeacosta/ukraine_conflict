# ==============================================================================
# STEP 09G - COMPETING RISKS MODEL ON TRUE PER-PERSON TRAJECTORIES
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Refits the dead-vs-alive competing-risks model from 09d, but on the
# counting-process intervals from 09f (Tstart, Tstop] per person, chained
# across all 4 register versions) instead of 09d's per-window "person-window"
# rows, which treated the same individual's re-observation in consecutive
# windows as independent Surv(0, time) draws.
#
# mstate::msprep() is built for single-episode-per-subject wide input (every
# subject starts at Tstart=0); it isn't a fit for pre-existing left-truncated
# multi-interval subject histories. Since 09f's intervals already have
# exactly the fields mstate's long format needs (Tstart, Tstop, status,
# trans, subject id, covariates), this script builds that long format
# directly rather than forcing msprep() to do something outside its design -
# msfit()/probtrans() only care that the columns are right, not how they got
# built. Multiple non-overlapping (Tstart,Tstop] intervals sharing one
# subject id is a standard Andersen-Gill counting-process Cox model.
#
# INPUT   data_inter/ukr_ualosses_intervals_countingprocess.rds   (from 09f)
# OUTPUTS data_inter/ukr_ualosses_traj_mstate_model.rds
#         data_inter/ukr_ualosses_traj_mstate_cifs.rds
#         data_inter/ukr_ualosses_traj_vs_window_comparison.rds
#         figures/exploratory/traj_cif_by_cohort.png
#         figures/exploratory/traj_vs_window_hazard_comparison.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

if (!require("mstate", quietly = TRUE)) {
  install.packages("mstate")
  library(mstate)
}

if (!file.exists("data_inter/ukr_ualosses_intervals_countingprocess.rds")) {
  stop("Run 09f_person_trajectories.R first.", call. = FALSE)
}

intervals <- read_rds("data_inter/ukr_ualosses_intervals_countingprocess.rds")
message(sprintf("Loaded %d counting-process intervals from %d people",
                nrow(intervals), n_distinct(intervals$person_id)))

# ==============================================================================
# 1. BUILD MSTATE-STYLE LONG FORMAT DIRECTLY (bypassing msprep())
# ==============================================================================
# Two rows per interval (one per competing transition), Tstart/Tstop carried
# through unchanged (real elapsed time, left-truncated for a person's later
# intervals - NOT reset to 0), status = 1 only on the row matching the
# interval's actual outcome.

tmat <- trans.comprisk(2, names = c("dead", "alive"))
print(tmat)

ms_long <- bind_rows(
  intervals %>% transmute(id = person_id, Tstart, Tstop, trans = 1,
                          status = to_dead, cohort),
  intervals %>% transmute(id = person_id, Tstart, Tstop, trans = 2,
                          status = to_alive, cohort)
) %>%
  arrange(id, Tstart, trans) %>%
  mutate(strata = trans)

# expand.covs() expects an mstate "msdata"-classed object with a `trans`
# attribute pointing at tmat - attach it manually since we built ms_long by
# hand rather than via msprep().
attr(ms_long, "trans") <- tmat
class(ms_long) <- c("msdata", class(ms_long))

ms_long <- expand.covs(ms_long, "cohort", longnames = FALSE)

message(sprintf("Long-format rows: %d", nrow(ms_long)))

# ==============================================================================
# 2. FIT THE COX MODEL
# ==============================================================================

cohort_dummy_cols <- setdiff(grep("^cohort", names(ms_long), value = TRUE), "cohort")

fit_formula <- as.formula(
  paste0("Surv(Tstart, Tstop, status) ~ ", paste(cohort_dummy_cols, collapse = " + "),
        " + strata(trans)")
)

cr_model <- coxph(fit_formula, data = ms_long, method = "breslow")

message("\n=== TRAJECTORY-BASED COMPETING RISKS MODEL ===\n")
print(summary(cr_model))

write_rds(cr_model, "data_inter/ukr_ualosses_traj_mstate_model.rds")

# ==============================================================================
# 3. CUMULATIVE INCIDENCE FUNCTIONS BY COHORT (same approach as 09d)
# ==============================================================================

cohort_levels <- levels(intervals$cohort)

cif_by_cohort <- function(cohort_year) {
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

  max_follow_up <- max(intervals$Tstop[intervals$cohort == cohort_year])
  pt <- pt[pt$time <= max_follow_up, ]

  if (any(!is.finite(as.matrix(pt[, -1])))) {
    warning(sprintf("Non-finite CIF values for cohort %s within its own follow-up range.", cohort_year))
  }

  as_tibble(pt) %>%
    transmute(
      cohort = cohort_year, time_days = time, time_months = time / 30.44,
      p_missing = pstate1, p_dead = pstate2, p_alive = pstate3,
      p_dead_se = se2,
      p_dead_lo = pmax(0, p_dead - 1.96 * p_dead_se),
      p_dead_hi = pmin(1, p_dead + 1.96 * p_dead_se)
    )
}

cifs <- map_dfr(cohort_levels, cif_by_cohort)

check_sum <- cifs %>% mutate(total = p_missing + p_dead + p_alive) %>%
  summarise(max_deviation = max(abs(total - 1)))
message(sprintf("\nCIF sum-to-1 check: max deviation = %.2e", check_sum$max_deviation))
stopifnot(check_sum$max_deviation < 1e-6)
message("✓ PASS")

write_rds(cifs, "data_inter/ukr_ualosses_traj_mstate_cifs.rds")

# ==============================================================================
# 4. COMPARE TO THE PER-WINDOW (INDEPENDENCE-VIOLATING) MODEL
# ==============================================================================

if (file.exists("data_inter/ukr_ualosses_mstate_model.rds")) {
  window_model <- read_rds("data_inter/ukr_ualosses_mstate_model.rds")

  window_hr <- broom::tidy(window_model, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(term, approach = "Per-window (09d)", hr = estimate, lo = conf.low, hi = conf.high)
  traj_hr <- broom::tidy(cr_model, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(term, approach = "Trajectory (09g)", hr = estimate, lo = conf.low, hi = conf.high)

  comparison <- bind_rows(window_hr, traj_hr) %>%
    mutate(ci_width_log = log(hi) - log(lo))  # width on log scale = symmetric, comparable across huge HR range

  message("\n=== HAZARD RATIO COMPARISON: per-window vs trajectory-based ===\n")
  print(comparison %>% arrange(term, approach) %>%
        mutate(across(c(hr, lo, hi), ~round(., 1))))

  width_comparison <- comparison %>%
    select(term, approach, ci_width_log) %>%
    pivot_wider(names_from = approach, values_from = ci_width_log) %>%
    mutate(pct_narrower = (`Per-window (09d)` - `Trajectory (09g)`) / `Per-window (09d)` * 100)

  message("\n=== CI WIDTH (log scale) COMPARISON ===\n")
  print(width_comparison)
  message(sprintf("\nMean CI narrowing from trajectory-based model: %.1f%%",
                  mean(width_comparison$pct_narrower, na.rm = TRUE)))

  write_rds(comparison, "data_inter/ukr_ualosses_traj_vs_window_comparison.rds")

  p_compare <- comparison %>%
    mutate(term = factor(term, levels = unique(term))) %>%
    ggplot(aes(x = term, y = hr, ymin = lo, ymax = hi, color = approach)) +
    geom_pointrange(position = position_dodge(width = 0.4)) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = c("Per-window (09d)" = "#95a5a6", "Trajectory (09g)" = "#c0392b")) +
    labs(x = "Cohort x Transition", y = "Hazard Ratio (log scale)",
        title = "Does Respecting True Per-Person Time Change the Estimates?",
        subtitle = "Per-window (independence-violating) vs. trajectory-based (proper counting-process) model",
        color = NULL) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

  ggsave("figures/exploratory/traj_vs_window_hazard_comparison.png", p_compare, w = 9, h = 6)
  message("✓ Saved: figures/exploratory/traj_vs_window_hazard_comparison.png")
} else {
  message("09d model not found; skipping comparison")
}

# ==============================================================================
# 5. PLOT: CIF BY COHORT (trajectory-based)
# ==============================================================================

p <- cifs %>%
  ggplot(aes(x = time_months, y = p_dead, color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = p_dead_lo, ymax = p_dead_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "Time Since Disappearance (months)", y = "Cumulative Incidence: P(Confirmed Dead)",
    title = "Trajectory-Based Competing Risks: Cumulative Incidence of Death by Cohort",
    subtitle = "Built from true per-person (Tstart,Tstop] intervals chained across v14/v16/v18/v19,\nnot independent per-window observations. Shaded band: 95% CI.",
    color = "Cohort", fill = "Cohort"
  ) +
  theme_bw() + theme(legend.position = "top")

ggsave("figures/exploratory/traj_cif_by_cohort.png", p, w = 8, h = 5.5)
message("✓ Saved: figures/exploratory/traj_cif_by_cohort.png")

message("\n✓ Trajectory-based competing risks model complete.")
