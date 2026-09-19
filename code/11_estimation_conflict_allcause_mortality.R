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
# Conflict death totals and net migration are both uncertain, so both are
# drawn n_sim times (00_setup.R) from PERT distributions built in 10. Every
# draw is projected through the full 2022-2025 accounting, which is what gives
# the results their uncertainty intervals.
#
# They are drawn differently, because they are uncertain in different ways.
# Each year's civilian total is its own unknown, so those draws are
# independent. The military total is uncertain through ONE parameter, alpha,
# the share of the never-resolved missing who are alive, so alpha is drawn
# once per simulation. Migration is uncertain in its COVERAGE - how much of the
# outflow the sources see - and that bias runs the same way in every year, so
# each migration component is drawn once per simulation and held across all
# four years.
#
# INPUTS   data_inter/ukr_param_table.rds                  (from 10)
#          data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds (from 04)
#          data_inter/ukr_pop_sssu.rds                     (from 01)
#          data_inter/ukr_migrants_...rds                  (from 06)
#          data_inter/ukr_asfr_wpp_2022_2025.rds           (from 05)
#          data_inter/ukr_ohchr_civilian_casualties.rds    (from 07_ohchr)
#          data_inter/ukr_ualosses_...rds                  (from 08)
#          data_inter/ukr_alpha_missing.rds, ukr_military_alpha_lines.rds (09)
#
# OUTPUTS  data_inter/ukr_sim_draws_2022_2025_n<n_sim>.rds  (every draw)
#          data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds (summary)
#          data_inter/ukr_sim_param_draws.rds      (the draws themselves, for 15)
#          data_inter/ukr_sim_alpha_draws.rds      (alpha in every draw, for 15)
#
# The life expectancy decomposition that used to live at the bottom of this
# script now has its own step, 14, which runs it inside every draw instead of
# once on the medians.
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

# the register's INDIVIDUALLY CONFIRMED deaths; anything drawn above them is
# attributable to the imputation of the missing, and that is the split
# reported as confirmed vs imputed
combatant_confirmed <- param_table %>%
  filter(role == "combatants") %>%
  select(year, conf_cmb = confirmed)

# --- counterfactual ("expected") mortality, no war --------------------------
exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") %>%
  filter(year %in% 2022:2025, source == "frcst") %>%
  select(-source)

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
set.seed(42)

# CIVILIANS are drawn independently in every year. What is uncertain is each
# year's own count, coming from event reporting specific to that year's
# events, and 2022 being at its ceiling says nothing about 2023.
draws_cvs_long <- param_table %>%
  filter(role == "civilians") %>%
  group_by(role, year) %>%
  reframe(
    sim_id = 1:n_sim,
    draw = rpert(n_sim, min = min, mode = mode, max = max)
  )

# COMBATANTS are not like that, although they look like a count. What is
# uncertain about them is alpha, the share of the never-resolved missing who
# are alive (00_setup.R, alpha_evidence(): the evidence, the sources and the
# range). That is ONE property of the missing, not four
# separate facts, so alpha is drawn once per simulation and every year's total
# follows from it: at_alpha0 - alpha x residual, the line 09's chain traces.
# Drawing the years independently would let one simulation put 2022 at "none
# of the missing alive" and 2025 at "a fifth of them alive", which is not a
# scenario anyone could defend, and the errors would partly cancel in the
# four-year total.
alpha_range <- read_rds("data_inter/ukr_alpha_missing.rds")
alpha_lines <- read_rds("data_inter/ukr_military_alpha_lines.rds")

alpha_draws <- tibble(
  sim_id = 1:n_sim,
  alpha = qpert(
    runif(n_sim),
    min = alpha_range$alpha_min,
    mode = alpha_range$alpha_mode,
    max = alpha_range$alpha_max
  )
)
write_rds(alpha_draws, "data_inter/ukr_sim_alpha_draws.rds")

draws_cmb_long <-
  expand_grid(alpha_draws, alpha_lines |> select(year, at_alpha0, residual)) |>
  transmute(role = "combatants", year, sim_id, draw = at_alpha0 - alpha * residual)

# the draws must stay inside 10's bounds, which are the same line evaluated
# at the ends of the range
stopifnot(
  draws_cmb_long |>
    left_join(param_table |> filter(role == "combatants") |> select(year, min, max), by = "year") |>
    with(all(draw >= min - 1e-6 & draw <= max + 1e-6))
)

draws_conflict_long <- bind_rows(draws_cvs_long, draws_cmb_long)

draws_conflict <- draws_conflict_long %>%
  pivot_wider(names_from = role, values_from = draw) %>%
  rename(draw_cmb = combatants, draw_cvs = civilians)

# Migration is not like that. What is uncertain is how much of the outflow
# the sources capture, and a register that keeps people after they return
# does so in every year alike. Drawing the years independently would let a
# simulation put 2022 on the register reading and 2023 on the crossing
# reading, which is not a scenario anyone could defend. So each component
# gets ONE quantile per simulation, held across all four years.
mig_component_draws <- function(rl) {
  u <- runif(n_sim)
  param_table %>%
    filter(role == rl) %>%
    group_by(role, year) %>%
    reframe(
      sim_id = 1:n_sim,
      draw = qpert(u, min = min, mode = mode, max = max)
    )
}

draws_mig_long <-
  bind_rows(
    mig_component_draws("mig_west"),
    mig_component_draws("mig_ru_by")
  )

draws_mig <- draws_mig_long %>%
  summarise(draw_mig = sum(draw), .by = c(sim_id, year))

draws_df <- draws_conflict %>%
  left_join(draws_mig, by = c("sim_id", "year")) %>%
  left_join(combatant_confirmed, by = "year")

# One tidy row per simulation, year and role. 15 draws its PERT figures from
# this rather than re-deriving the draws from the projection output, and 16
# uses it to attribute the spread in the results back to its inputs.
param_draws <- bind_rows(draws_conflict_long, draws_mig_long)
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
