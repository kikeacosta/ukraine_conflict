# ==============================================================================
# STEP 13c - Sensitivity of the life expectancy loss to net migration
# ==============================================================================
#
# 13 answers how much of the loss is a denominator effect of net migration, by
# comparing the mode with a scenario in which nobody leaves. This answers the
# companion question: how far does the loss move across the plausible range of
# the net migration estimate itself?
#
# Net migration enters the simulation as two components (06, 10), each drawn once
# per simulation and held across the four years (11):
#   mig_west   the western outflow, bracketed by source: each draw blends the
#              crossings reading's whole path with the register reading's by
#              one Beta-PERT(0, 0.5, 1) weight, so the mode is their midpoint
#   mig_ru_by  displacement into Russia and Belarus, a single 2022 entry drawn
#              from a PERT
#
# Each component is swept across its whole range, along the quantile u that
# 11 draws (so every year moves together, as in the simulation), while the
# other component and every conflict input stay at their mode: the western
# component from the crossings reading's path to the register reading's,
# Russia and Belarus across its PERT bounds. One deterministic projection per
# grid point, as in 13b.
#
# The sweep stays inside each component's range. Below it, the additive
# allocation along the departures profile that run_single_sim() uses would
# start to push age cells negative, so a curve extended towards zero would be
# an artefact of the allocation rather than a result. The loss with
# no migration at all is 13's scenario, reported alongside.
#
# Section 5 then checks three specification choices of the migration input,
# one at a time at the mode (documents/migration_methodology.md, section 10):
# where the western mode sits in the bracket (A2), the age-sex profile of the
# Russia and Belarus flow (A9), and additive against proportional allocation
# (A5). 15 reports them as table A8.
#
# INPUTS   the same as 13
#          data_inter/ukr_migration_decomposition.rds (from 13)
#          data_inter/ukr_pop_sssu.rds, Donetsk and Luhansk (from 01)
# OUTPUTS  data_inter/ukr_migration_sensitivity_e0.rds
#          data_inter/ukr_migration_specification_e0.rds (section 5)
#          figures/exploratory/migration_sensitivity_e0.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# 1. INPUTS, AS IN 13 ==========================================================
param_table <- read_rds("data_inter/ukr_param_table.rds")

exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") %>%
  filter(year %in% 2022:2025, source == "frcst") %>%
  select(-source)

pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") %>%
  filter(reg == "cnt", year == 2022) %>%
  select(-reg)

migs2 <- read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds") %>%
  mutate(ems = -mix) |>
  select(-mix)

asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds")
asfr2 <- asfr %>%
  bind_rows(
    asfr %>% filter(year == 2023) %>% mutate(year = 2024),
    asfr %>% filter(year == 2023) %>% mutate(year = 2025)
  ) %>%
  arrange(year, age)

ohchr2 <- read_rds("data_inter/ukr_ohchr_civilian_casualties.rds") %>%
  select(year, sex, age, prop_cvs = cx)

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

profile_props <-
  cmb_profile("dead", "prop_cmb_dead") |>
  left_join(cmb_profile("missing", "prop_cmb_miss"), by = c("year", "sex", "age")) |>
  left_join(ohchr2, by = c("year", "sex", "age"))

static_inputs <-
  expand_grid(year = 2022:2025, sex = c("f", "m"), age = 0:100) %>%
  left_join(migs2, by = c("year", "sex", "age")) %>%
  left_join(exp_mort2, by = c("year", "sex", "age")) %>%
  left_join(asfr2, by = c("year", "sex", "age")) %>%
  left_join(profile_props, by = c("year", "sex", "age")) %>%
  replace_na(list(ems = 0, w = 0, mx = 0, fx = 0, prop_cvs = 0,
                  prop_cmb_dead = 0, prop_cmb_miss = 0))

stopifnot(all(static_inputs$mx > 0), !any(is.na(static_inputs)))

# conflict inputs held at their mode
conflict_at_mode <-
  param_table |>
  filter(role %in% c("civilians", "combatants")) |>
  select(year, role, mode) |>
  pivot_wider(names_from = role, values_from = mode) |>
  rename(draw_cvs = civilians, draw_cmb = combatants) |>
  left_join(
    param_table |> filter(role == "combatants") |> select(year, conf_cmb = confirmed),
    by = "year"
  )

mig_bounds <- param_table |> filter(str_starts(role, "mig_")) |>
  select(year, role, min, mode, max, crossings, register)

# 2. THE GRID ==================================================================
# u is the PERT quantile of the swept component, shared by its years as in 11.
# The 2.5% and 97.5% points bound the 95% interval of the simulated draws, and
# the mode is added exactly (for the western component it is the median,
# because its blend weight is symmetric; for Russia and Belarus it is not).
# rounded so the equality tests below are exact
u_grid <- round(seq(0, 1, by = 0.025), 3)

# as 11 draws: the western blend of the two readings' paths, Russia and
# Belarus through its own PERT
component_at <- function(rl, u) {
  rows <- mig_bounds |> filter(role == rl)
  if (rl == "mig_west") {
    w <- qpert(u, min = 0, mode = 0.5, max = 1)
    rows |> transmute(year, value = crossings + w * (register - crossings))
  } else {
    rows |> transmute(year, value = qpert(u, min = min, mode = mode, max = max))
  }
}
component_at_mode <- function(rl) {
  mig_bounds |> filter(role == rl) |> transmute(year, value = mode)
}

scenarios <-
  expand_grid(component = c("mig_west", "mig_ru_by"), u = u_grid) |>
  mutate(point = case_when(
    u == 0 ~ "min", u == 0.025 ~ "p2.5", u == 0.975 ~ "p97.5", u == 1 ~ "max",
    .default = NA_character_
  )) |>
  bind_rows(tibble(component = c("mig_west", "mig_ru_by"), u = NA_real_, point = "mode"))

run_scenario <- function(comp, u) {
  other <- setdiff(c("mig_west", "mig_ru_by"), comp)
  swept <- if (is.na(u)) component_at_mode(comp) else component_at(comp, u)
  draw_mig <-
    bind_rows(swept, component_at_mode(other)) |>
    summarise(draw_mig = sum(value), .by = year)

  draws_df <-
    conflict_at_mode |>
    left_join(draw_mig, by = "year") |>
    mutate(sim_id = 1)

  sim <-
    run_single_sim(1, draws_df, static_inputs, pop22_ini) |>
    mutate(
      conflict = civilian + combatant_confirmed + combatant_imputed,
      mx_all = (expected + conflict) / pop,
      mx_bsn = expected / pop
    )

  e0_of <- function(d, col) {
    d |>
      select(year, sex, age, mx = all_of(col)) |>
      group_by(year, sex) |>
      do(lifetable(dt_in = .data)) |>
      ungroup() |>
      filter(age == 0) |>
      select(year, sex, ex)
  }

  e0_of(sim, "mx_bsn") |>
    rename(ex_bsn = ex) |>
    left_join(e0_of(sim, "mx_all") |> rename(ex_all = ex), by = c("year", "sex")) |>
    mutate(
      loss = ex_bsn - ex_all,
      cumulative_net_migration = sum(swept$value),
      min_exposure = min(sim$pop)
    )
}

# 3. RUN =======================================================================
migration_sensitivity <-
  scenarios |>
  mutate(res = map2(component, u, run_scenario)) |>
  unnest(res)

# a curve through negative exposure would be meaningless
stopifnot(all(migration_sensitivity$min_exposure >= 0))

# the mode scenario must reproduce 13's loss at the mode
check_mode <-
  migration_sensitivity |>
  filter(point == "mode") |>
  left_join(read_rds("data_inter/ukr_migration_decomposition.rds") |> select(year, sex, loss13 = loss),
            by = c("year", "sex"))
stopifnot(isTRUE(all.equal(check_mode$loss, check_mode$loss13, tolerance = 1e-8)))

write_rds(
  migration_sensitivity |> select(component, u, point, cumulative_net_migration, year, sex, loss),
  "data_inter/ukr_migration_sensitivity_e0.rds"
)

cat("\n=== 2025 LOSS ACROSS EACH COMPONENT'S RANGE ===\n")
print(as.data.frame(
  migration_sensitivity |>
    filter(!is.na(point), year == 2025) |>
    select(component, point, cumulative_net_migration, sex, loss) |>
    pivot_wider(names_from = sex, values_from = loss, names_prefix = "loss_") |>
    arrange(component, cumulative_net_migration)
))

# 4. DIAGNOSTIC PLOT (the manuscript figure is assembled in 15) ================
migration_sensitivity |>
  filter(!is.na(u)) |>
  mutate(sex = if_else(sex == "f", "Females", "Males")) |>
  ggplot(aes(cumulative_net_migration / 1e6, loss, colour = factor(year))) +
  geom_line(linewidth = 1) +
  facet_wrap(component ~ sex, scales = "free") +
  labs(x = "Cumulative net migration of the component, 2022-2025 (millions)",
       y = "Life expectancy loss (years)", colour = "Year") +
  theme_bw()
ggsave("figures/exploratory/migration_sensitivity_e0.png", w = 9, h = 6)

message("Done. data_inter/ukr_migration_sensitivity_e0.rds written for 15.")

# ==============================================================================
# 5. SPECIFICATION CHECKS
# ==============================================================================
# Three assumptions in documents/migration_methodology.md (section 10) are
# tested by changing one of them at a time, deterministically, with every
# conflict input at its mode:
#
#   A2  the western mode is the midpoint of the two readings. Alternative: a
#       third of the way from the crossings reading to the register reading.
#   A9  Russia and Belarus share the 2022 western departures profile, which the
#       exit ban on men makes young and female. Alternative: the 1 January
#       2022 age-sex structure of Donetsk and Luhansk oblasts, as if whole
#       households left - older and more male.
#   A5  a draw away from the mode is added along the departures profile.
#       Alternative: the signed age-sex profile scaled in proportion, which
#       also scales returns. At the mode the two coincide, so the check is run
#       at both ends of the western bracket.
loss_of <- function(static, draw_mig) {
  sim <-
    run_single_sim(1, conflict_at_mode |> left_join(draw_mig, by = "year") |> mutate(sim_id = 1),
                   static, pop22_ini) |>
    mutate(conflict = civilian + combatant_confirmed + combatant_imputed,
           mx_all = (expected + conflict) / pop, mx_bsn = expected / pop)
  stopifnot(min(sim$pop) >= 0)
  e0 <- function(col) {
    sim |> select(year, sex, age, mx = all_of(col)) |> group_by(year, sex) |>
      do(lifetable(dt_in = .data)) |> ungroup() |> filter(age == 0) |> select(year, sex, ex)
  }
  e0("mx_bsn") |> rename(ex_bsn = ex) |>
    left_join(e0("mx_all") |> rename(ex_all = ex), by = c("year", "sex")) |>
    transmute(year, sex, loss = ex_bsn - ex_all)
}

west_at <- function(w) {
  mig_bounds |> filter(role == "mig_west") |> transmute(year, value = crossings + w * (register - crossings))
}
totals_with_west <- function(w) {
  bind_rows(west_at(w), component_at_mode("mig_ru_by")) |> summarise(draw_mig = sum(value), .by = year)
}
mode_totals <- totals_with_west(0.5)

# A9: move the Russia and Belarus flow (2022) from the departures profile w to
# the Donetsk-Luhansk population structure. The population on 1 January 2022
# at age a is the cohort that sits in row a + 1 of 2022 (00_setup.R).
ru_by_mode <- mig_bounds |> filter(role == "mig_ru_by") |> pull(mode)
dl_profile <-
  read_rds("data_inter/ukr_pop_sssu.rds") |>
  filter(reg %in% c("dnk", "luk"), year == 2022) |>
  summarise(pop = sum(pop), .by = c(sex, age)) |>
  mutate(age = pmin(age + 1, 100)) |>
  summarise(pop = sum(pop), .by = c(sex, age)) |>
  mutate(p_dl = pop / sum(pop), year = 2022L) |>
  select(year, sex, age, p_dl)
static_a9 <-
  static_inputs |>
  left_join(dl_profile, by = c("year", "sex", "age")) |>
  mutate(ems = if_else(year == 2022, ems - ru_by_mode * w + ru_by_mode * coalesce(p_dl, 0), ems)) |>
  select(-p_dl)
stopifnot(isTRUE(all.equal(sum(static_a9$ems), sum(static_inputs$ems))))

# A5: the signed profile scaled in proportion to the draw, instead of the
# difference added along the departures profile
proportional <- function(totals) {
  static_inputs |>
    left_join(totals |> rename(target = draw_mig), by = "year") |>
    left_join(mode_totals |> rename(base = draw_mig), by = "year") |>
    mutate(ems = ems * target / base) |>
    select(-target, -base)
}

specification <- bind_rows(
  loss_of(static_inputs, mode_totals) |> mutate(scenario = "Mode"),
  loss_of(static_inputs, totals_with_west(1 / 3)) |> mutate(scenario = "A2: western mode a third of the way from the crossings reading"),
  loss_of(static_a9, mode_totals) |> mutate(scenario = "A9: Russia and Belarus on the Donetsk-Luhansk population structure"),
  loss_of(static_inputs, totals_with_west(0)) |> mutate(scenario = "A5: crossings end, additive (as used)"),
  loss_of(proportional(totals_with_west(0)), mode_totals) |> mutate(scenario = "A5: crossings end, proportional"),
  loss_of(static_inputs, totals_with_west(1)) |> mutate(scenario = "A5: register end, additive (as used)"),
  loss_of(proportional(totals_with_west(1)), mode_totals) |> mutate(scenario = "A5: register end, proportional")
)

# the first row must reproduce 13's loss at the mode
stopifnot(isTRUE(all.equal(
  specification |> filter(scenario == "Mode") |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(specification, "data_inter/ukr_migration_specification_e0.rds")

cat("\n=== SPECIFICATION CHECKS: LOSS BY YEAR, MEN, AND 2025 WOMEN ===\n")
print(as.data.frame(
  specification |>
    mutate(col = if_else(sex == "m", paste0("m_", year), paste0("f_", year))) |>
    filter(sex == "m" | year == 2025) |>
    select(scenario, col, loss) |>
    pivot_wider(names_from = col, values_from = loss) |>
    mutate(across(where(is.numeric), \(x) round(x, 3)))
))

message("Done. data_inter/ukr_migration_specification_e0.rds written for 15.")
