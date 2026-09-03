# ==============================================================================
# STEP 11 - Probabilistic estimation of conflict and all-cause mortality
# Ukraine, 2022-2025
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Projects the Ukrainian population forward from 1 January 2022 one year at a
# time, and in each year splits total mortality into:
#
#   expected  - what mortality would have been without the war (Lee-Carter
#               forecast from the 2000-2019 trend, see 04)
#   conflict  - civilian deaths (UCDP totals, OHCHR age-sex profile) plus
#               combatant deaths (ualosses register), the latter further
#               split into confirmed and imputed-from-missing
#
# Conflict death TOTALS are uncertain, so they are drawn 5,000 times from PERT
# distributions defined by the min/mode/max in the parameter table built in 10.
# Every draw is projected through the full 2022-2025 accounting, which is what
# gives the results their credible intervals.
#
# INPUTS   data_inter/ukr_param_table.rds                  (from 10)
#          data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds (from 04)
#          data_inter/ukr_pop_sssu.rds                     (from 01)
#          data_inter/ukr_migrants_...rds                  (from 06)
#          data_inter/ukr_asfr_wpp_2022_2025.rds           (from 05)
#          data_inter/ukr_ohchr_civilian_casualties.rds    (from 07_ohchr)
#          data_inter/ukr_ualosses_...rds                  (from 08)
#
# OUTPUTS  data_inter/ukr_sim_draws_2022_2025.parquet  (all 5,000 draws)
#          data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds (summary)
#
# The life expectancy decomposition that used to live at the bottom of this
# script now has its own step, 14, which runs it inside every draw instead of
# once on the medians.
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# 1. SIMULATION CONTROL ========================================================
n_sim <- 5000

# 2. BASELINE DEMOGRAPHIC INPUTS ===============================================

# --- parameter table: min / mode / max conflict deaths per year and role -----
param_table <- read_rds("data_inter/ukr_param_table.rds")

# the combatant minimum is the count of INDIVIDUALLY CONFIRMED deaths in the
# ualosses register; anything drawn above it is attributable to the imputation
# of the missing, and that is the split reported as confirmed vs imputed
combatant_mins <- param_table %>%
  filter(role == "combatants") %>%
  select(year, min_cmb = min)

# --- counterfactual ("expected") mortality, no war --------------------------
exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") %>%
  filter(year %in% 2022:2025, source == "frcst") %>%
  select(-source)

# --- starting population: 1 January 2022, government-controlled mainland -----
# region "cnt" excludes the Donetsk and Luhansk regions
pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") %>%
  filter(reg == "cnt", year == 2022) %>%
  select(-reg)

# --- net emigration (refugee outflow), positive = people leaving ------------
migs2 <- read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds") %>%
  mutate(ems = -mix) |>
  select(-mix)

# --- fertility --------------------------------------------------------------
# WPP2024 only publishes observed ASFR through 2023, so the most recent
# observed schedule (2023) is carried forward to 2024 and 2025.
asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds")
asfr2 <- asfr %>%
  bind_rows(
    asfr %>% filter(year == 2023) %>% mutate(year = 2024),
    asfr %>% filter(year == 2023) %>% mutate(year = 2025)
  ) %>%
  arrange(year, age)

stopifnot(setequal(unique(asfr2$year), 2022:2025))

# 3. AGE-SEX PROFILES OF CONFLICT DEATHS =======================================
# The draws give yearly TOTALS; these profiles say how those totals are spread
# over age and sex. Both are normalised to sum to 1 within each year.

# civilians: OHCHR casualty records, ungrouped to single years of age
ohchr2 <- read_rds("data_inter/ukr_ohchr_civilian_casualties.rds") %>%
  select(year, sex, age, prop_cvs = cx)

# combatants: age-sex distribution of the confirmed dead in the ualosses register
ual2 <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") %>%
  filter(status == "dead", year %in% 2022:2025) |>
  select(-status) |>
  complete(year = 2022:2025, sex, age = 0:100, fill = list(dx = 0)) |>
  mutate(prop_cmb = dx / sum(dx), .by = year) |>
  select(-dx)

profile_props <- ual2 |> left_join(ohchr2, by = c("year", "sex", "age"))

# 4. STATIC INPUT GRID =========================================================
# One row per year-sex-age cell, carrying everything that does NOT vary
# between simulations.
static_inputs <- expand_grid(
  year = 2022:2025,
  sex = c("f", "m"),
  age = 0:100
) %>%
  left_join(migs2, by = c("year", "sex", "age")) %>%
  left_join(exp_mort2, by = c("year", "sex", "age")) %>%
  left_join(asfr2, by = c("year", "sex", "age")) %>%
  left_join(profile_props, by = c("year", "sex", "age")) %>%
  replace_na(list(ems = 0, mx = 0, fx = 0, prop_cvs = 0, prop_cmb = 0))

# every cell must have a real mortality rate: an mx of 0 here would mean the
# join silently failed and that cohort would be projected as immortal
stopifnot(
  nrow(static_inputs) == 4 * 2 * 101,
  all(static_inputs$mx > 0),
  !any(is.na(static_inputs))
)

# 5. THE PROJECTION FOR ONE DRAW ===============================================
# run_single_sim() is defined in 00_setup.R, so this script and 13 share a
# single implementation of the cohort-component projection.

# 6. RUN (OR REUSE) THE 5,000 DRAWS ============================================
# ~7 minutes. Cached, so re-running this script is instant unless the parquet
# is deleted. The file is large and is NOT tracked in git - see .gitignore.
sim_output_raw <- cache_parquet("data_inter/ukr_sim_draws_2022_2025.parquet", {
  set.seed(42) # inside the cache block, so the cache is reproducible

  draws_df <- param_table %>%
    group_by(role, year) %>%
    reframe(
      sim_id = 1:n_sim,
      draw = rpert(n_sim, min = min, mode = mode, max = max)
    ) %>%
    pivot_wider(names_from = role, values_from = draw) %>%
    rename(draw_cmb = combatants, draw_cvs = civilians) %>%
    left_join(combatant_mins, by = "year")

  draws_list <- draws_df %>% group_split(sim_id)

  message("  running ", n_sim, " simulations ...")
  map_df(1:n_sim, function(i) {
    run_single_sim(
      sim_id = i,
      draws_this_sim = draws_list[[i]],
      static_inputs = static_inputs,
      pop22_ini = pop22_ini
    )
  })
})

# 7. DERIVED CAUSES AND RATES ==================================================
# Stored wide (components only) and expanded here, so the cached file is six
# times smaller than the long form and nothing can drift out of sync.
sim_long <-
  sim_output_raw %>%
  mutate(
    conflict = civilian + combatant_confirmed + combatant_imputed,
    all = expected + conflict
  ) %>%
  pivot_longer(
    c(
      all, conflict, expected, civilian,
      combatant_confirmed, combatant_imputed
    ),
    names_to = "cause",
    values_to = "dx"
  ) %>%
  mutate(mx = dx / pop)

# 8. SUMMARY ACROSS DRAWS ======================================================
sce_all_probabilistic <- sim_long %>%
  summarise(
    mx_median = median(mx, na.rm = TRUE),
    mx_lower = quantile(mx, 0.025, na.rm = TRUE),
    mx_upper = quantile(mx, 0.975, na.rm = TRUE),
    dx_median = median(dx, na.rm = TRUE),
    pop_median = median(pop, na.rm = TRUE),
    .by = c(year, sex, age, cause)
  ) %>%
  arrange(year, sex, age, cause)

write_rds(
  sce_all_probabilistic,
  "data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds"
)

# quick sanity read-out
sim_output_raw %>%
  summarise(
    civilian = sum(civilian),
    combatant = sum(combatant_confirmed + combatant_imputed),
    .by = c(sim_id, year)
  ) %>%
  summarise(
    civilian = median(civilian),
    combatant = median(combatant),
    .by = year
  ) %>%
  print()

message("Done. Run 14 for the life expectancy decomposition by cause.")
