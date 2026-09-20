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
# Every input that is uncertain and can be given a distribution is drawn n_sim
# times (simulation_draws(), 00_setup.R), and every draw is projected through
# the full 2022-2025 accounting, which is what gives the results their
# uncertainty intervals: civilian deaths, the evidence on the missing alive,
# the sampling error of the resolution model and of the registration-lag
# factors, the Lee-Carter forecast of the counterfactual, and net migration.
#
# They are drawn differently, because they are uncertain in different ways.
# Each year's civilian total is its own unknown, so those draws are
# independent. The evidence on the missing alive describes the missing, not a
# year, so it is drawn once per simulation and sets every year's military
# total together. Migration is uncertain in its COVERAGE - how much of the
# outflow the sources see - and that bias runs the same way in every year, so
# each migration component is drawn once per simulation and held across all
# four years. The forecast index follows its random walk from year to year.
#
# INPUTS   data_inter/ukr_param_table.rds                  (from 10)
#          data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds (from 04)
#          data_inter/ukr_lc_forecast_error.rds            (from 04)
#          data_inter/ukr_pop_sssu.rds                     (from 01)
#          data_inter/ukr_migrants_...rds                  (from 06)
#          data_inter/ukr_asfr_wpp_2022_2025.rds           (from 05)
#          data_inter/ukr_ohchr_civilian_casualties.rds    (from 07_ohchr)
#          data_inter/ukr_ualosses_...rds                  (from 08)
#          data_inter/ukr_military_inputs.rds              (from 09)
#
# OUTPUTS  data_inter/ukr_sim_draws_2022_2025_n<n_sim>.rds  (every draw)
#          data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds (summary)
#          data_inter/ukr_sim_param_draws.rds      (the draws themselves, for 15)
#          data_inter/ukr_sim_alive_draws.rds      (the evidence on the missing
#                                                   alive in every draw, for 15)
#
# The life expectancy decomposition that used to live at the bottom of this
# script now has its own step, 14, which runs it inside every draw instead of
# once on a summary of them.
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# 1. SIMULATION CONTROL ========================================================
# n_sim comes from 00_setup.R (1,000 while the pipeline is being built, 20,000
# for the reported estimates), or from UKR_N_SIM for a one-off run. The cache
# below is keyed by it, so a run at one size can never be served another
# size's draws.
message("simulation size: n_sim = ", n_sim)

# 2. BASELINE DEMOGRAPHIC INPUTS ===============================================

# --- parameter table: min / mode / max conflict deaths per year and role -----
param_table <- read_rds("data_inter/ukr_param_table.rds")

# --- counterfactual ("expected") mortality, no war --------------------------
# the forecast means, and what each draw needs to move them along its own
# path of the forecast index: the loading bx and the forecast variance vk
exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") %>%
  filter(year %in% 2022:2025, source == "frcst") %>%
  select(-source)
lc_error <- read_rds("data_inter/ukr_lc_forecast_error.rds")

# --- starting population: 1 January 2022, continental Ukraine ----------------
# region "cnt" is the whole country without Crimea and Sevastopol; it includes
# the Donetsk and Luhansk regions
pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") %>%
  filter(reg == "cnt", year == 2022) %>%
  select(-reg)

# --- net migration (outflow), positive = people leaving ----------------------
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

# combatants: TWO age-sex distributions, not one. The register's confirmed dead
# and its missing do not share an age profile - the missing are older, with
# 6.7 points more of their mass in ages 40-49 and 5.6 points less in 20-29, and
# a mean age about a year higher. Spreading the whole combatant total over the
# confirmed-dead profile therefore placed the imputed deaths at ages where they
# did not occur, which matters because the Arriaga decomposition in 14 is
# age-weighted and younger deaths cost more life expectancy each.
#
# So registered deaths take the profile of the register's dead, and imputed
# deaths take the profile of the missing they are imputed from.
ual_raw <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")

cmb_profile <- function(keep_status, nm) {
  ual_raw |>
    filter(status == keep_status, year %in% 2022:2025) |>
    select(-status) |>
    complete(year = 2022:2025, sex, age = 0:100, fill = list(dx = 0)) |>
    summarise(dx = sum(dx), .by = c(year, sex, age)) |>
    mutate("{nm}" := dx / sum(dx), .by = year) |>
    select(-dx)
}

ual2 <-
  cmb_profile("dead", "prop_cmb_dead") |>
  left_join(cmb_profile("missing", "prop_cmb_miss"), by = c("year", "sex", "age"))

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
  left_join(lc_error$bx, by = c("sex", "age")) %>%
  left_join(lc_error$variance |> select(sex, year, vk), by = c("sex", "year")) %>%
  replace_na(list(ems = 0, w = 0, mx = 0, fx = 0, prop_cvs = 0,
                  prop_cmb_dead = 0, prop_cmb_miss = 0))

# every cell must have a real mortality rate: an mx of 0 here would mean the
# join silently failed and that cohort would be projected as immortal
stopifnot(
  nrow(static_inputs) == 4 * 2 * 101,
  all(static_inputs$mx > 0),
  !any(is.na(static_inputs)),
  # each draw is spread along this departures profile, so it must sum to 1
  all(abs(tapply(static_inputs$w, static_inputs$year, sum) - 1) < 1e-9),
  # and each combatant profile must sum to 1 within a year, or the drawn
  # totals would not be preserved once they are spread over age
  all(abs(tapply(static_inputs$prop_cmb_dead, static_inputs$year, sum) - 1) < 1e-9),
  all(abs(tapply(static_inputs$prop_cmb_miss, static_inputs$year, sum) - 1) < 1e-9)
)

# 5. THE PROJECTION FOR ONE DRAW ===============================================
# run_single_sim() is defined in 00_setup.R, so this script and 13 share a
# single implementation of the cohort-component projection.

# 6. DRAW THE PARAMETERS =======================================================
# Drawing is cheap and deterministic; the projection below is neither. Keeping
# the two apart means the draws are on disk for the figures even when the
# projection is read back from cache, and the projection stays a pure function
# of them.
# simulation_draws() (00_setup.R) holds the sampling, so 13g can repeat it
# with another PERT shape on the same random numbers. What each input's draw
# stands for:
#
# CIVILIANS are drawn independently in every year. What is uncertain is each
# year's own count, coming from event reporting specific to that year's
# events, and 2022 being at its ceiling says nothing about 2023.
#
# COMBATANTS are not like that, although they look like a count. What is
# uncertain about them is how many of the missing are alive (00_setup.R,
# alive_evidence(): the evidence, the sources and the ranges): the prisoners
# of war among them and the share of the rest alive for other reasons. Those
# are properties of the missing, not four separate facts, so they are drawn
# once per simulation and every year's total follows from them
# (military_draws()), together with a draw of the resolution model's
# estimates and of the registration-lag factors. Drawing the years
# independently would let one simulation put 2022 at "none of the missing
# alive" and 2025 at "a fifth of them alive", which is not a scenario anyone
# could defend, and the errors would partly cancel in the four-year total.
#
# THE COUNTERFACTUAL is a forecast, and its error grows with the horizon: each
# simulation follows its own path of each sex's Lee-Carter index, the two
# sexes' steps correlated as 04 measured.
#
# MIGRATION is uncertain in how much of the outflow the sources capture, and a
# register that keeps people after they return does so in every year alike.
# The WESTERN component is bracketed by source: one weight w per simulation,
# drawn from a symmetric Beta-PERT on [0, 1] with mode 0.5, blends the two
# readings' whole paths - w = 0 is the crossings reading in every year, w = 1
# the register reading - so its four-year total runs between the two sources'
# own totals, and the mode is their midpoint in every year (see 10). RUSSIA AND
# BELARUS is one 2022 entry with its own PERT bounds.
mil <- read_rds("data_inter/ukr_military_inputs.rds")
sampled <- simulation_draws(param_table, mil, lc_error, n = n_sim)
draws_df <- sampled$draws_df
param_draws <- sampled$param_draws
write_rds(sampled$alive_draws, "data_inter/ukr_sim_alive_draws.rds")

# the draws of the evidence alone, the model and the lag at their point
# values, must stay inside 10's bounds, which are the same total at the ends
# of the evidence
stopifnot(
  param_draws |>
    filter(role == "mil_alive") |>
    left_join(param_table |> filter(role == "combatants") |> select(year, min, max), by = "year") |>
    with(all(draw >= min - 1e-6 & draw <= max + 1e-6)),
  !any(is.na(draws_df))
)

# One tidy row per simulation, year and role. 15 draws its PERT figures from
# this rather than re-deriving the draws from the projection output, and 16
# uses it to attribute the spread in the results back to its inputs.
write_rds(param_draws, "data_inter/ukr_sim_param_draws.rds")

# 7. RUN (OR REUSE) THE PROJECTIONS ============================================
# About 7 minutes at 5,000 draws and 28 at 20,000. Cached, so re-running this
# script is instant unless the .rds file is deleted.
#
# THE FILE NAME CARRIES THE DRAW COUNT. Without it, a cache built during
# construction at 1,000 draws would be handed straight back to a production run
# asking for 20,000, and every interval would be reported at a size it was
# never computed at. Keying the name also means switching between the two sizes
# costs nothing, because each keeps its own cache.
#
# Large, and NOT tracked in git - see .gitignore.
sim_draws_file <- sprintf("data_inter/ukr_sim_draws_2022_2025_n%d.rds", n_sim)
sim_output_raw <- cache_rds(sim_draws_file, {
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

# 8. DERIVED CAUSES AND RATES ==================================================
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

# 9. SUMMARY ACROSS DRAWS ======================================================
# THE ESTIMATE IS THE MEAN OF THE DRAWS, here and in every step that summarises
# them, with the 2.5th and 97.5th percentiles as its interval. The mean is what
# the central projection of the deterministic steps stands for - every input at
# the centre of its distribution - so the two agree; and means add up, across
# ages, sexes, years and components, which medians do not, so every table sums
# to its own total.
sce_all_probabilistic <- sim_long %>%
  summarise(
    mx_mean = mean(mx, na.rm = TRUE),
    mx_lower = quantile(mx, 0.025, na.rm = TRUE),
    mx_upper = quantile(mx, 0.975, na.rm = TRUE),
    dx_mean = mean(dx, na.rm = TRUE),
    pop_mean = mean(pop, na.rm = TRUE),
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
    civilian = mean(civilian),
    combatant = mean(combatant),
    .by = year
  ) %>%
  print()

message("Done. Run 14 for the life expectancy decomposition by cause.")
