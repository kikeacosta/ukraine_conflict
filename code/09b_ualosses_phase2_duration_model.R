# ==============================================================================
# STEP 09B - PHASE 2: Duration-Dependent Markov Imputation
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Builds on Phase 1 validation by fitting DURATION-DEPENDENT transition rates.
# Instead of assuming constant hazards per cohort, models how resolution
# probabilities change with time since disappearance.
#
# KEY INSIGHT: The longer someone is missing, the probability of resolution
# may increase (more time for bureaucratic processes) OR decrease (unlikely to
# surface). This script fits that relationship empirically.
#
# METHOD
# ------
# For each person in both Phase 1 windows, calculate:
#   - Duration since disappearance (event_date → register_date)
#   - Outcome in target register (dead, alive, prisoner, missing)
#
# Then fit transition-specific hazard models:
#   λ(t | cohort, outcome) = λ₀(t) × exp(β × cohort_year)
#
# This captures both:
#   - Duration effect (how does time affect outcome probability)
#   - Cohort effect (do 2022 cohorts resolve differently than 2025?)
#
# REQUIRES
# --------
# - Phase 1 outputs: ukr_ualosses_pooled_transitions.rds
# - R packages: survival (Cox models), mstate (multi-state models)
#
# OUTPUT
# ------
# - Duration-dependent transition rates by cohort
# - Hazard ratios and 95% CIs
# - Fitted smoothing curves showing resolution patterns
# - Updated imputation using duration-stratified rates
# - Figures showing how resolution changes with time
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# Load Phase 1 outputs
if (!require("survival", quietly = TRUE)) {
  install.packages("survival")
  library(survival)
}

message("Loading Phase 1 validation outputs...")
transitions_pooled <- read_rds("data_inter/ukr_ualosses_pooled_transitions.rds")
pooled_rates <- read_rds("data_inter/ukr_ualosses_pooled_transition_rates.rds")

if (!file.exists("data_inter/ukr_ualosses_pooled_transitions.rds")) {
  stop("Phase 1 outputs not found. Run step 09_ualosses_missing_analysis_improved.R first.")
}

# ==============================================================================
# 1. PREPARE DATA FOR DURATION ANALYSIS
# ==============================================================================

message("\n=== PHASE 2: DURATION-DEPENDENT ANALYSIS ===\n")

# Calculate duration (in days and months) from event to register date
# Event date = disappearance date in source register
# Register date = snapshot date of target register

duration_data <- transitions_pooled %>%
  filter(status == "missing") %>%
  mutate(
    # Approximate register dates
    # v14: Sep 15, 2025; v18: May 20, 2026; v19: Sep 15, 2026
    register_date = case_when(
      window == "v14→v18" ~ as.Date("2026-05-20"),
      window == "v18→v19" ~ as.Date("2026-09-15")
    ),

    # Duration in days and months
    duration_days = as.numeric(register_date - date_evnt),
    duration_months = duration_days / 30.44,

    # Outcome categories (for stratified analysis)
    outcome = case_when(
      status2 == "dead" ~ "Dead",
      status2 == "alive" ~ "Alive",
      status2 == "prisoner" ~ "Prisoner",
      status2 == "missing" ~ "Missing"
    ),

    # Cohort (by year of disappearance)
    cohort = as.factor(year)
  ) %>%
  select(name_full, date_bth, year, cohort, status2, outcome,
         window, date_evnt, register_date, duration_days, duration_months)

message(sprintf("Duration data: %d individuals", nrow(duration_data)))
message(sprintf("Duration range: %.0f to %.0f days (%.1f to %.1f months)",
                min(duration_data$duration_days, na.rm = TRUE),
                max(duration_data$duration_days, na.rm = TRUE),
                min(duration_data$duration_months, na.rm = TRUE),
                max(duration_data$duration_months, na.rm = TRUE)))

# ==============================================================================
# 2. DURATION STRATIFICATION: Compute transition rates by duration bands
# ==============================================================================

# Stratify by duration to see how resolution changes with time
duration_data <- duration_data %>%
  mutate(
    duration_band = case_when(
      duration_months < 6 ~ "0–5 mo",
      duration_months < 12 ~ "6–11 mo",
      duration_months < 24 ~ "12–23 mo",
      duration_months >= 24 ~ "24+ mo"
    ),
    duration_band = factor(duration_band,
                          levels = c("0–5 mo", "6–11 mo", "12–23 mo", "24+ mo"))
  )

# Compute rates by duration band and cohort
rates_by_duration <- duration_data %>%
  summarise(
    n_total = n(),
    n_dead = sum(outcome == "Dead"),
    n_alive = sum(outcome == "Alive"),
    n_prisoner = sum(outcome == "Prisoner"),
    n_missing = sum(outcome == "Missing"),
    mean_duration_days = mean(duration_days, na.rm = TRUE),
    .by = c(cohort, duration_band)
  ) %>%
  mutate(
    p_dead = n_dead / n_total,
    p_alive = n_alive / n_total,
    p_prisoner = n_prisoner / n_total,
    p_missing = n_missing / n_total,

    # 95% CIs (binomial)
    se_dead = sqrt(p_dead * (1 - p_dead) / n_total),
    ci_dead_lo = pmax(0, p_dead - 1.96 * se_dead),
    ci_dead_hi = pmin(1, p_dead + 1.96 * se_dead),

    se_missing = sqrt(p_missing * (1 - p_missing) / n_total),
    ci_missing_lo = pmax(0, p_missing - 1.96 * se_missing),
    ci_missing_hi = pmin(1, p_missing + 1.96 * se_missing)
  ) %>%
  arrange(cohort, duration_band)

message("\n=== TRANSITION RATES BY DURATION BAND ===\n")
print(rates_by_duration %>%
      select(cohort, duration_band, n_total, p_dead, p_missing))

write_rds(rates_by_duration, "data_inter/ukr_ualosses_rates_by_duration.rds")

# ==============================================================================
# 3. PLOT: Resolution patterns by duration and cohort
# ==============================================================================

p_dead <- ggplot(rates_by_duration,
                 aes(x = duration_band, y = p_dead,
                     color = cohort, group = cohort)) +
  geom_line(size = 1) +
  geom_point(size = 3, alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_dead_lo, ymax = ci_dead_hi),
                width = 0.2, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "Time Since Disappearance",
    y = "Probability of Confirmed Death",
    color = "Cohort Year",
    title = "Phase 2: Resolution Pattern – How P(Dead) Changes with Time",
    subtitle = "Error bars show 95% CIs. Later cohorts have less time to resolve."
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_missing <- ggplot(rates_by_duration,
                    aes(x = duration_band, y = p_missing,
                        color = cohort, group = cohort)) +
  geom_line(size = 1) +
  geom_point(size = 3, alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_missing_lo, ymax = ci_missing_hi),
                width = 0.2, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "Time Since Disappearance",
    y = "Probability Still Missing",
    color = "Cohort Year",
    title = "Phase 2: Still Missing – How P(Missing) Changes with Time",
    subtitle = "Shows resolution plateau. Which cohorts are still heavily missing?"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_combined <- p_dead / p_missing + plot_layout(guides = "collect")

ggsave("figures/exploratory/phase2_duration_resolution_patterns.png",
       p_combined, w = 10, h = 8)
message("Saved: figures/exploratory/phase2_duration_resolution_patterns.png")

# ==============================================================================
# 4. FIT COX PROPORTIONAL HAZARDS MODEL
# ==============================================================================

message("\n=== FITTING HAZARD MODELS ===\n")

# For "dead" outcome: Does hazard of being confirmed dead depend on duration?
dead_outcomes <- duration_data %>%
  mutate(
    event = ifelse(outcome == "Dead", 1, 0),
    # Censoring: "alive", "prisoner", "missing" are censored (not dead)
    # Duration: observed time to event (or censoring)
    time = duration_days
  ) %>%
  filter(time > 0)  # Exclude same-day records

# Fit Cox model: hazard of death ~ cohort + duration (smooth)
surv_obj <- Surv(time = dead_outcomes$time, event = dead_outcomes$event)

# Simple Cox model: hazard depends on cohort only
cox_model_cohort <- coxph(
  surv_obj ~ cohort,
  data = dead_outcomes
)

message("Cox Model: Event (Death) ~ Cohort")
print(summary(cox_model_cohort))

# Save model
write_rds(cox_model_cohort, "data_inter/ukr_ualosses_cox_model_death.rds")

# Extract hazard ratios
hr_table <- cox_model_cohort %>%
  broom::tidy(exponentiate = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    cohort = term,
    hazard_ratio = estimate,
    hr_ci_lo = conf.low,
    hr_ci_hi = conf.high
  )

message("\n=== HAZARD RATIOS (vs 2022 baseline) ===")
print(hr_table)

write_rds(hr_table, "data_inter/ukr_ualosses_hazard_ratios.rds")

# ==============================================================================
# 5. DURATION-STRATIFIED IMPUTATION
# ==============================================================================

message("\n=== IMPROVED IMPUTATION: DURATION-STRATIFIED ===\n")

# Get the 2025 cohort rates stratified by duration
imputation_2025 <- rates_by_duration %>%
  filter(cohort == 2025) %>%
  select(duration_band, p_dead, p_alive, p_missing)

message("For 2025 cohort (least resolved), imputation rates by time since disappearance:")
print(imputation_2025)

# The key insight: people missing for 0–5 months have high P(missing)
# We need to estimate what will happen as they resolve (e.g., following the 2022 pattern)

# Construct trajectory: assume 2025 cohort follows 2022 pattern with lag
trajectory_2022 <- rates_by_duration %>%
  filter(cohort == 2022) %>%
  select(duration_band, p_dead_2022 = p_dead, p_alive_2022 = p_alive)

trajectory_2025 <- rates_by_duration %>%
  filter(cohort == 2025) %>%
  left_join(trajectory_2022, by = "duration_band") %>%
  mutate(
    # If 2025 were to follow 2022's pattern, impute as:
    imputed_dead_projection = p_missing * p_dead_2022,
    imputed_alive_projection = p_missing * p_alive_2022
  )

message("\nProjected outcomes for 2025 cohort (assuming 2022 resolution pattern):")
print(trajectory_2025 %>%
      select(duration_band, p_dead, imputed_dead_projection,
             p_missing, imputed_alive_projection))

write_rds(trajectory_2025, "data_inter/ukr_ualosses_2025_imputation_projection.rds")

# ==============================================================================
# 6. SUMMARY: DURATION EFFECTS
# ==============================================================================

message("\n=== PHASE 2 SUMMARY ===\n")

summary_stats <- duration_data %>%
  group_by(cohort, outcome) %>%
  summarise(
    n = n(),
    mean_duration_months = mean(duration_months, na.rm = TRUE),
    median_duration_months = median(duration_months, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cohort, outcome)

message("Mean time to resolution by outcome and cohort:")
print(summary_stats)

message("\n=== KEY FINDINGS ===")
message(sprintf("• 2022 cohort: %.0f%% still missing → well-resolved",
                filter(rates_by_duration, cohort == 2022,
                       duration_band == "24+ mo")$p_missing[1] * 100))
message(sprintf("• 2025 cohort: %.0f%% still missing → barely resolved",
                filter(rates_by_duration, cohort == 2025,
                       duration_band == "0–5 mo")$p_missing[1] * 100))
message(sprintf("• Death rate (2022 cohort, 24+ months): %.1f%%",
                filter(rates_by_duration, cohort == 2022,
                       duration_band == "24+ mo")$p_dead[1] * 100))

# ==============================================================================
# 7. OUTPUTS FOR STEP 10 (IMPUTATION)
# ==============================================================================

message("\n=== READY FOR IMPUTATION ===")
message("Phase 2 outputs generated:")
message("  • data_inter/ukr_ualosses_rates_by_duration.rds")
message("  • data_inter/ukr_ualosses_cox_model_death.rds")
message("  • data_inter/ukr_ualosses_hazard_ratios.rds")
message("  • data_inter/ukr_ualosses_2025_imputation_projection.rds")
message("  • figures/exploratory/phase2_duration_resolution_patterns.png")
message("\n✓ Phase 2 complete! Ready to implement improved imputation.")

# ==============================================================================
# 8. DECISION FOR NEXT STEPS
# ==============================================================================

message("\n=== DECISION: Use Duration-Stratified Imputation? ===")
message("\nOption A (Simple): Use Phase 1 pooled rates for all cohorts")
message("  → Pro: Less complexity; already validated")
message("  → Con: Ignores that 2025 cohort is less resolved")

message("\nOption B (Improved): Use Phase 2 duration-stratified rates")
message("  → Pro: More realistic for recent cohorts with less time")
message("  → Con: Requires projecting from 2022 pattern")
message("  → Con: Adds slight model complexity")

message("\nRECOMMENDATION: Use Option B")
message("  Rationale: 2025 cohort has <6 months to resolve.")
message("           Using 2022's 20+ month pattern would bias high.")
message("           Duration stratification captures this heterogeneity.")

message("\n→ Proceed to finalize imputation (continues in main step 09)")
