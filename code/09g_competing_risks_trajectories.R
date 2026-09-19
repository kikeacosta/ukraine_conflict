# ==============================================================================
# STEP 09G - COMPETING RISKS ON PER-PERSON TRAJECTORIES, AND THE ADOPTION TEST
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Fits a competing-risks Cox model to 09f's counting-process intervals: from
# "missing", a person can move to dead, prisoner or no longer listed (alive,
# by the rule production uses), with time measured in days since the event
# and event-year cohort as covariate.
# The multi-state coxph() of the survival package fits one baseline hazard and
# one set of cohort effects per transition, handles left truncation (a
# person's later interval starts where the earlier one ended), and its
# survfit() gives Aalen-Johansen state probabilities by cohort.
#
# THE TEST, agreed before the model was re-run (2026-09-18). The model may be
# reconsidered for production - at least for cohort-specific bounds - only if
# BOTH hold, with prisoners and people no longer listed kept apart (released
# prisoners, held rather than disappeared, are outside the population at risk):
#   1. the longest-missing cohort (2022) has the highest probability of
#      death, read from its cumulative incidence at its own maximum follow-up;
#   2. the cohort hazard ratios are plausible. Fixed here, before looking:
#      every cohort hazard ratio lies between 1/50 and 50.
# Otherwise the model stays a documented investigation and the share of the
# missing who are alive stays an explicit, evidence-based range (09).
#
# INPUT   data_inter/ukr_ualosses_intervals_countingprocess.rds (09f, gitignored)
# OUTPUTS data_inter/ukr_ualosses_traj_mstate_model.rds          (gitignored)
#         data_inter/ukr_ualosses_traj_mstate_cifs.rds            (aggregate)
#         data_inter/ukr_ualosses_multistate_test.rds             (aggregate)
#         figures/exploratory/traj_cif_by_cohort.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")
library(survival)

if (!file.exists("data_inter/ukr_ualosses_intervals_countingprocess.rds")) {
  stop("Run 09f_person_trajectories.R first.", call. = FALSE)
}
intervals <-
  read_rds("data_inter/ukr_ualosses_intervals_countingprocess.rds") |>
  mutate(event = factor(
    case_when(to_dead == 1 ~ "dead", to_prisoner == 1 ~ "prisoner",
              to_unlisted == 1 ~ "unlisted",
              .default = "censor"),
    levels = c("censor", "dead", "prisoner", "unlisted")
  ))
message(sprintf("Loaded %d intervals from %d people", nrow(intervals), n_distinct(intervals$person_id)))

# ==============================================================================
# 1. MULTI-STATE COX MODEL
# ==============================================================================
fit <- coxph(Surv(Tstart, Tstop, event) ~ cohort, data = intervals, id = person_id)
print(fit)
write_rds(fit, "data_inter/ukr_ualosses_traj_mstate_model.rds")

# ==============================================================================
# 2. STATE PROBABILITIES BY COHORT, within each cohort's own follow-up
# ==============================================================================
cohort_levels <- levels(intervals$cohort)
cif_by_cohort <- function(cy) {
  sf <- survfit(fit, newdata = data.frame(cohort = factor(cy, levels = cohort_levels)))
  p <- sf$pstate
  if (length(dim(p)) == 3) p <- p[, 1, ]
  colnames(p) <- sf$states
  follow_max <- max(intervals$Tstop[intervals$cohort == cy])
  as_tibble(p) |>
    mutate(time_days = sf$time) |>
    filter(time_days <= follow_max) |>
    transmute(cohort = cy, time_days, p_missing = .data[[sf$states[1]]],
              p_dead = dead, p_prisoner = prisoner,
              p_unlisted = unlisted)
}
cifs <- map_dfr(cohort_levels, cif_by_cohort)
stopifnot(all(abs(cifs$p_missing + cifs$p_dead + cifs$p_prisoner +
                    cifs$p_unlisted - 1) < 1e-6))
write_rds(cifs, "data_inter/ukr_ualosses_traj_mstate_cifs.rds")

# ==============================================================================
# 3. THE TEST
# ==============================================================================
at_max_follow_up <-
  cifs |>
  slice_max(time_days, n = 1, by = cohort, with_ties = FALSE) |>
  mutate(days_at_risk_min = map_dbl(cohort, \(cy) min(intervals$Tstart[intervals$cohort == cy])),
         days_at_risk_max = time_days)

co <- summary(fit)$coefficients
hr <- tibble(term = rownames(co), hr = exp(co[, "coef"]),
             lo = exp(co[, "coef"] - 1.96 * co[, "se(coef)"]),
             hi = exp(co[, "coef"] + 1.96 * co[, "se(coef)"]))

criterion_1 <- at_max_follow_up$cohort[which.max(at_max_follow_up$p_dead)] == "2022"
criterion_2 <- all(hr$hr >= 1 / 50 & hr$hr <= 50)

test <- list(
  at_max_follow_up = at_max_follow_up,
  hazard_ratios = hr,
  longest_missing_highest_p_dead = criterion_1,
  hazard_ratios_plausible = criterion_2,
  pass = criterion_1 && criterion_2
)
write_rds(test, "data_inter/ukr_ualosses_multistate_test.rds")

cat("\n=== MULTISTATE ADOPTION TEST ===\n")
print(as.data.frame(at_max_follow_up |> mutate(across(starts_with("p_"), \(x) round(x, 3)))))
print(as.data.frame(hr |> mutate(across(c(hr, lo, hi), \(x) signif(x, 3)))))
cat(sprintf("\n1. 2022 cohort has the highest P(dead): %s\n", criterion_1))
cat(sprintf("2. every cohort hazard ratio within [1/50, 50]: %s\n", criterion_2))
cat(sprintf("RESULT: %s\n", if (test$pass) "PASS - reconsider the model for cohort-specific bounds"
            else "FAIL - the model stays a documented investigation"))

# ==============================================================================
# 4. PLOT
# ==============================================================================
cifs |>
  pivot_longer(c(p_dead, p_prisoner, p_unlisted), names_to = "outcome", values_to = "p") |>
  mutate(outcome = factor(outcome, levels = c("p_dead", "p_prisoner", "p_unlisted"),
                          labels = c("Dead", "Prisoner", "No longer listed"))) |>
  ggplot(aes(time_days / 30.44, p, colour = cohort)) +
  geom_line(linewidth = 1) +
  facet_wrap(~outcome, nrow = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Months since disappearance", y = "Cumulative incidence", colour = "Cohort",
       title = "Competing risks from missing, by event-year cohort (within each cohort's follow-up)") +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/exploratory/traj_cif_by_cohort.png", w = 12, h = 4.5)

message("Done.")
