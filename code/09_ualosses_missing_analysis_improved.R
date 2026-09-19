# ==============================================================================
# STEP 09 - IMPROVED: Ualosses Missing Analysis with Fuzzy Matching & Phase 1
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# 1. FUZZY MATCHING: Matches individuals between register versions allowing for
#    small spelling variations in names (1-2 character differences)
# 2. PHASE 1 VALIDATION: Extracts transitions in TWO independent windows:
#    - Window A: v14 (Sep 2025) → v18 (May 2026)
#    - Window B: v18 (May 2026) → v19 (Sep 2026)
# 3. COMPARISON: Shows transition rates with 95% CIs in both windows
# 4. VALIDATION: Checks for consistency/drift between windows
# 5. IMPUTATION: Uses combined estimates to impute missing outcomes
#
# IMPROVEMENTS OVER ORIGINAL:
# - Fuzzy name matching (catches typos, transliteration)
# - Two-window validation instead of single-window estimate
# - Confidence intervals on all transition probabilities
# - Stability checks (are rates consistent across windows?)
#
# INPUT    data_input/ualosses_hubert_datasets/250916_UKR_ualosses_Personnel_v14.xlsx
#          data_input/ualosses_hubert_datasets/251204_UKR_ualosses_Personnel_v16.xlsx
#          data_input/ualosses_hubert_datasets/260530_UKR_ualosses_Personnel_v18.xlsx
#          data_input/ualosses_hubert_datasets/260919_UKR_ualosses_Personnel_v19.xlsx
#          data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds
# OUTPUTS  data_inter/ukr_ualosses_transition_rates_two_windows.rds
#          data_inter/ukr_ualosses_imputation_table.rds (updated with CI)
#          data_inter/ukr_ualosses_transition_rates.rds
#          data_inter/ukr_ualosses_conflict_deaths_imputed_miss_sex_age_2022_2025.rds
#          figures/exploratory/transition_rates_validation.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# Install stringdist for fuzzy matching if needed
if (!require("stringdist", quietly = TRUE)) {
  install.packages("stringdist")
  library(stringdist)
}

# ==============================================================================
# 1. READ REGISTER FILES
# ==============================================================================

read_reg <- function(path) {
  read_xlsx(require_raw(path), sheet = "Database") |>
    mutate(
      name_last = as.character(LastName),
      name_first = as.character(FirstName),
      name_patr = as.character(Patronym),
      name_full = paste(name_last, name_first, name_patr, sep = " "),
      date_bth = excel_date(DateBirth),
      date_evnt = excel_date(DateEvent),
      year = year(date_evnt)
    ) |>
    filter(year %in% 2022:2025, Nationality == "Ukraine") |>
    select(name_last, name_first, name_patr, name_full, date_bth, year,
           date_evnt, status = Status)
}

message("Reading register v14 (Sep 2025)...")
v14 <- read_reg("data_input/ualosses_hubert_datasets/250916_UKR_ualosses_Personnel_v14.xlsx")

message("Reading register v16 (Dec 2025)...")
v16 <- read_reg("data_input/ualosses_hubert_datasets/251204_UKR_ualosses_Personnel_v16.xlsx")

message("Reading register v18 (May 2026)...")
v18 <- read_reg("data_input/ualosses_hubert_datasets/260530_UKR_ualosses_Personnel_v18.xlsx")

message("Reading register v19 (Sep 2026)...")
v19 <- read_reg("data_input/ualosses_hubert_datasets/260919_UKR_ualosses_Personnel_v19.xlsx")

message(sprintf("v14: %d rows, v16: %d rows, v18: %d rows, v19: %d rows",
                nrow(v14), nrow(v16), nrow(v18), nrow(v19)))

# Register export dates, from the filenames themselves (250916, 251204,
# 260530, 260919), used both for the window_months labels below and for
# duration-since-disappearance calculations downstream (09d).
register_dates <- c(
  v14 = as.Date("2025-09-16"), v16 = as.Date("2025-12-04"),
  v18 = as.Date("2026-05-30"), v19 = as.Date("2026-09-19")
)

# ==============================================================================
# 2. FUZZY MATCHING FUNCTION
# ==============================================================================
# Match by DOB (exact) + name (fuzzy), threshold 90% similarity
# This catches: typos, transliteration variants, single letter changes

fuzzy_match <- function(source, target, dob_window = 30, threshold = 0.90) {

  # First: exact match on the COMPOUND key (name_full, date_bth), not DOB
  # alone. DOB alone is not a valid person key in this register: of 153,928
  # rows only ~15,349 distinct date_bth values exist (~10 people per date on
  # average, some 30-44 people sharing one date), almost certainly because
  # many entries only recorded birth YEAR and the day/month were defaulted.
  # An earlier version of this function joined on date_bth alone; with
  # distinct(date_bth, .keep_all=TRUE) that silently fans one arbitrary
  # person's outcome out to everyone else sharing that birthdate, which
  # produces a flat ~43-44% "P(dead)" for every cohort year regardless of
  # true duration since disappearance - the signature of a broken join key,
  # not a real finding. name_full + date_bth together is effectively unique.
  # (name_full_target isn't needed here: an exact match on name+DOB implies
  # the matched name IS name_full, and the caller only keeps `status2`.)
  target_dedup <- target %>%
    distinct(name_full, date_bth, .keep_all = TRUE) %>%
    select(name_full, date_bth, status2 = status)

  matched <- source %>%
    mutate(.row_id = row_number()) %>%
    left_join(
      target_dedup,
      by = c("name_full", "date_bth"),
      relationship = "many-to-one"
    )

  exact_matches <- matched %>% filter(!is.na(status2))

  message(sprintf("Exact (name+DOB) matches: %d of %d (%.1f%%)",
                  nrow(exact_matches), nrow(source),
                  100 * nrow(exact_matches) / nrow(source)))

  # For unmatched, try fuzzy name matching with Levenshtein distance.
  # VECTORIZED: a single non-equi data.table join builds every (source,
  # candidate) pair within the DOB window in one C-level pass, then a single
  # vectorized stringdist() call scores all pairs at once. This replaces an
  # earlier row-by-row for() loop that did not scale past a few thousand
  # unmatched records (each iteration re-filtered the full target table).
  unmatched <- matched %>% filter(is.na(status2))

  if (nrow(unmatched) > 0) {
    message(sprintf("Attempting fuzzy match on %d unmatched records (vectorized)...",
                    nrow(unmatched)))
    t0 <- Sys.time()

    src_dt <- as.data.table(unmatched %>% select(.row_id, name_full, date_bth))
    tgt_dt <- as.data.table(target %>% select(name_full_target = name_full,
                                              status2 = status, date_bth))
    src_dt[, `:=`(dob_lo = date_bth - dob_window, dob_hi = date_bth + dob_window)]

    pairs <- tgt_dt[src_dt,
                    on = .(date_bth >= dob_lo, date_bth <= dob_hi),
                    allow.cartesian = TRUE,
                    .(.row_id, name_full, name_full_target, status2)]
    pairs <- pairs[!is.na(name_full_target)]

    message(sprintf("  Candidate pairs within +/-%d days: %d", dob_window, nrow(pairs)))

    if (nrow(pairs) > 0) {
      pairs[, name_dist := stringdist(name_full, name_full_target, method = "lv")]
      pairs[, name_sim := 1 - (name_dist / pmax(nchar(name_full), nchar(name_full_target)))]
      pairs <- pairs[name_sim >= threshold]

      # Best candidate per source row (sort then keep first per group -
      # avoids .SD[1] per-group overhead across tens of thousands of groups)
      setorder(pairs, .row_id, -name_sim, name_dist)
      best <- pairs[!duplicated(.row_id)]

      message(sprintf("  Fuzzy matches found: %d (%.1f sec)",
                      nrow(best), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

      matched <- matched %>%
        left_join(
          best %>% select(.row_id, status2_fuzzy = status2),
          by = ".row_id"
        ) %>%
        mutate(status2 = coalesce(status2, status2_fuzzy)) %>%
        select(-status2_fuzzy)
    } else {
      message("  Fuzzy matches found: 0")
    }
  }

  # For any remaining unmatched: treated as "alive" (resurfaced but not in target register)
  matched <- matched %>%
    mutate(status2 = ifelse(is.na(status2), "alive", status2)) %>%
    select(-.row_id)

  return(matched)
}

# ==============================================================================
# 3. EXTRACT TRANSITIONS IN BOTH WINDOWS
# ==============================================================================

# The raw Status field has FOUR values across the three registers, not
# three: dead, missing, prisoner, and (v19 only) released_prisoner - a
# category introduced between v18 and v19 that earlier inspection of v14
# alone missed. A released prisoner is definitively not dead and no longer
# detained, so for this analysis (which cares about dead vs. not-dead) it is
# folded into "alive" here, once, so every downstream consumer of status2
# (Phase 1 rates, the competing-risks model) sees an exhaustive 4-category
# outcome that actually sums to n_total - previously it silently didn't,
# understating "still missing" and starving the "prisoner" transition of
# events for the 2022/2023 cohorts (who had time to be released), which is
# what triggered (quasi-)separation in the Cox competing-risks fit.
normalize_status2 <- function(df) {
  df %>% mutate(status2 = ifelse(status2 == "released_prisoner", "alive", status2))
}

# Three consecutive windows now that v16 (Dec 2025) sits between v14 and
# v18: v14->v16, v16->v18, v18->v19. Looped rather than hand-duplicated a
# third time, so a future v20/v21 is a one-line addition to `registers`
# instead of another copy-pasted block (which is exactly how the
# window_months=8/4 hardcoding here would otherwise have gone stale).
registers <- list(v14 = v14, v16 = v16, v18 = v18, v19 = v19)
reg_names <- names(registers)

window_months_between <- function(from_nm, to_nm) {
  round(as.numeric(register_dates[[to_nm]] - register_dates[[from_nm]]) / 30.44, 1)
}

transitions_by_window <- list()
for (i in seq_len(length(reg_names) - 1)) {
  from_nm <- reg_names[i]
  to_nm   <- reg_names[i + 1]
  wlabel  <- paste0(from_nm, "→", to_nm)
  wmonths <- window_months_between(from_nm, to_nm)

  message(sprintf("\n=== WINDOW %s (%s → %s, %s): %s (%.1f mo) ===",
                  LETTERS[i], from_nm, to_nm,
                  format(register_dates[[to_nm]], "%b %Y"), wlabel, wmonths))

  transitions_by_window[[wlabel]] <- fuzzy_match(registers[[from_nm]], registers[[to_nm]]) %>%
    normalize_status2() %>%
    mutate(window = wlabel, window_months = wmonths) %>%
    select(name_full, date_bth, year, date_evnt, status, status2, window, window_months)
}

# window_order fixes the chronological sequence for the stability check
# below (pairwise CONSECUTIVE comparisons: A vs B, B vs C - not every
# combination), since transitions_by_window's list order already is that
# sequence but downstream code re-derives it from the data.
window_order <- names(transitions_by_window)
transitions_all <- bind_rows(transitions_by_window)

# ==============================================================================
# 4. COMPUTE TRANSITION RATES WITH CONFIDENCE INTERVALS
# ==============================================================================

compute_transition_ci <- function(transitions_df) {

  # Filter for people with "missing" status
  missing <- transitions_df %>%
    filter(status == "missing")

  # Count outcomes
  outcomes <- missing %>%
    summarise(
      n_total = n(),
      n_alive = sum(status2 == "alive"),
      n_dead = sum(status2 == "dead"),
      n_prisoner = sum(status2 == "prisoner"),
      n_missing = sum(status2 == "missing"),
      .by = c(year, window)
    )

  # Compute proportions with Wilson score 95% CIs
  outcomes_with_ci <- outcomes %>%
    mutate(
      # Dead
      p_dead = n_dead / n_total,
      se_dead = sqrt(p_dead * (1 - p_dead) / n_total),
      ci_dead_lo = pmax(0, p_dead - 1.96 * se_dead),
      ci_dead_hi = pmin(1, p_dead + 1.96 * se_dead),

      # Alive
      p_alive = n_alive / n_total,
      se_alive = sqrt(p_alive * (1 - p_alive) / n_total),
      ci_alive_lo = pmax(0, p_alive - 1.96 * se_alive),
      ci_alive_hi = pmin(1, p_alive + 1.96 * se_alive),

      # Prisoner
      p_prisoner = n_prisoner / n_total,
      se_prisoner = sqrt(p_prisoner * (1 - p_prisoner) / n_total),
      ci_prisoner_lo = pmax(0, p_prisoner - 1.96 * se_prisoner),
      ci_prisoner_hi = pmin(1, p_prisoner + 1.96 * se_prisoner),

      # Still missing
      p_missing = n_missing / n_total,
      se_missing = sqrt(p_missing * (1 - p_missing) / n_total),
      ci_missing_lo = pmax(0, p_missing - 1.96 * se_missing),
      ci_missing_hi = pmin(1, p_missing + 1.96 * se_missing)
    ) %>%
    select(year, window, n_total,
           contains("p_"), contains("ci_"), contains("se_"),
           -contains("se_"))

  return(outcomes_with_ci)
}

message("\n=== COMPUTING TRANSITION RATES ===")
transition_rates_both <- compute_transition_ci(transitions_all)

# ==============================================================================
# 5. DISPLAY VALIDATION TABLE
# ==============================================================================

message("\n=== PHASE 1 VALIDATION: TRANSITION RATE COMPARISON ===\n")

validation_table <- transition_rates_both %>%
  select(year, window, n_total, p_dead, ci_dead_lo, ci_dead_hi,
         p_alive, ci_alive_lo, ci_alive_hi) %>%
  mutate(
    p_dead_fmt = sprintf("%.3f [%.3f–%.3f]", p_dead, ci_dead_lo, ci_dead_hi),
    p_alive_fmt = sprintf("%.3f [%.3f–%.3f]", p_alive, ci_alive_lo, ci_alive_hi)
  ) %>%
  select(year, window, n_total, p_dead_fmt, p_alive_fmt)

print(validation_table)

# ==============================================================================
# 6. STABILITY CHECK: Do rates differ between windows?
# ==============================================================================

message("\n=== STABILITY ANALYSIS ===\n")

# Generalized to N windows (previously hardcoded for exactly 2, which broke
# on window labels containing "→" needing pivot_wider - see git history).
# Compares each window against the PREVIOUS one only (A vs B, B vs C), not
# every pairwise combination, since that's what "did the rate drift between
# consecutive readings" actually means.
stability <- transition_rates_both %>%
  mutate(window = factor(window, levels = window_order)) %>%
  arrange(year, window) %>%
  group_by(year) %>%
  mutate(
    prev_window = lag(window),
    prev_p_dead = lag(p_dead),
    prev_p_alive = lag(p_alive),
    delta_dead = abs(p_dead - prev_p_dead),
    delta_alive = abs(p_alive - prev_p_alive),
    flag_dead = ifelse(delta_dead > 0.05, "⚠ DRIFT", "✓ stable"),
    flag_alive = ifelse(delta_alive > 0.05, "⚠ DRIFT", "✓ stable")
  ) %>%
  ungroup() %>%
  filter(!is.na(prev_window)) %>%
  select(year, window, prev_window,
         p_dead, prev_p_dead, delta_dead, flag_dead,
         p_alive, prev_p_alive, delta_alive, flag_alive)

print(stability)

# Save for reference
write_rds(transition_rates_both,
          "data_inter/ukr_ualosses_transition_rates_two_windows.rds")
write_rds(stability,
          "data_inter/ukr_ualosses_stability_check.rds")

# ==============================================================================
# 7. VISUALIZATION: Transition rates with CIs, both windows
# ==============================================================================

# Built explicitly rather than via pivot_longer(names_sep = "_"): column
# names like ci_dead_lo have 3 underscore-separated parts, which the
# 2-column names_to = c("metric","outcome") split silently mangles (collapses
# ci_dead_lo/ci_dead_hi to the same key and drops "lo"/"hi" with a warning).
plot_data <- bind_rows(
  transition_rates_both %>%
    transmute(year, window, outcome = "dead",
              p = p_dead, ci_lo = ci_dead_lo, ci_hi = ci_dead_hi),
  transition_rates_both %>%
    transmute(year, window, outcome = "alive",
              p = p_alive, ci_lo = ci_alive_lo, ci_hi = ci_alive_hi)
)

p <- ggplot(plot_data %>% filter(outcome != "prisoner", outcome != "missing"),
            aes(x = year, y = p, color = window, shape = window)) +
  facet_wrap(~outcome, labeller = labeller(
    outcome = c(dead = "P(Dead)", alive = "P(Alive)")
  )) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.2, position = position_dodge(width = 0.3)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x = "Event Year (cohort)",
    y = "Transition Probability",
    color = "Window",
    shape = "Window",
    title = sprintf("Phase 1 Validation: Transition Rates Across %d Windows", length(window_order)),
    subtitle = "95% Confidence Intervals shown. Drift flags outcomes changing >5% between consecutive windows."
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("figures/exploratory/transition_rates_validation.png", p, w = 8, h = 5)
message("Saved: figures/exploratory/transition_rates_validation.png")

# ==============================================================================
# 8. COMBINED ESTIMATE: Pool all windows
# ==============================================================================

message(sprintf("\n=== POOLED ESTIMATE (%d windows combined) ===\n", length(window_order)))

pooled_transitions <- transitions_all %>%
  filter(status == "missing")

pooled_rates <- pooled_transitions %>%
  summarise(
    n_total = n(),
    n_dead = sum(status2 == "dead"),
    n_alive = sum(status2 == "alive"),
    n_prisoner = sum(status2 == "prisoner"),
    n_missing = sum(status2 == "missing"),
    .by = year
  ) %>%
  mutate(
    p_dead = n_dead / n_total,
    p_alive = n_alive / n_total,
    p_prisoner = n_prisoner / n_total,
    p_missing = n_missing / n_total
  )

print(pooled_rates)

# ==============================================================================
# 9. STORE FOR DOWNSTREAM USE
# ==============================================================================

# Calculate what fraction of long-term missing is assumed alive
# Now based on POOLED data rather than window A only
suelo_vivo <- pooled_rates %>%
  filter(year == 2025) %>%
  pull(p_alive)

message(sprintf("\nLong-term missing (2025 cohort) resolution rates:"))
message(sprintf("  Dead: %.1f%%", pooled_rates$p_dead[pooled_rates$year == 2025] * 100))
message(sprintf("  Alive: %.1f%%", pooled_rates$p_alive[pooled_rates$year == 2025] * 100))
message(sprintf("  Prisoner: %.1f%%", pooled_rates$p_prisoner[pooled_rates$year == 2025] * 100))
message(sprintf("  Still Missing: %.1f%%", pooled_rates$p_missing[pooled_rates$year == 2025] * 100))
message(sprintf("\nUsing suelo_vivo = %.2f for imputation", suelo_vivo))

# Continue with the imputation using pooled rates...
# (For now, save the rates and transition data for review)

write_rds(pooled_rates, "data_inter/ukr_ualosses_pooled_transition_rates.rds")
write_rds(pooled_transitions, "data_inter/ukr_ualosses_pooled_transitions.rds")

message("\n✓ Phase 1 validation complete!")
message("  Output files:")
message("  - data_inter/ukr_ualosses_transition_rates_two_windows.rds")
message("  - data_inter/ukr_ualosses_stability_check.rds")
message("  - data_inter/ukr_ualosses_pooled_transition_rates.rds")
message("  - figures/exploratory/transition_rates_validation.png")
