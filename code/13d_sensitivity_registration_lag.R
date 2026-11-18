# ==============================================================================
# STEP 13d - Sensitivity to the registration-lag correction
# ==============================================================================
#
# WHY
# ---
# 08b brings v19's dead and missing to the completeness events reach at four
# years, and 09 imputes on those completed counts. Four years is where the
# register's growth is last measured on enough event months, not where it
# stops: at 42-48 months the dead still grow by about 0.07% a month and the
# missing by about 0.2%. So the horizon makes the correction a lower bound on
# the lag. This step shows how far the results move without the correction, and
# with the growth extrapolated beyond four years at a constant rate; 13l fits a
# curve to the tail instead. It also carries the two alternatives 08b keeps to
# the way the growth is measured.
#
# SCENARIOS
# ---------
#   as registered            v19's counts, no correction
#   to 48 months (used)      08b's chain-ladder, the production estimate
#   to 60 / 72 months        the same, with the rate of the 42-48 month band
#                            carried on beyond 48 months
#   growth per person        to 48 months, with the growth of each status
#   listed                   measured per person the register lists for the
#                            event month, dead or missing, not per record of
#                            the status's own stock
#   releases dated by        to 48 months, with lags and windows counted from
#   Kaggle version           the day each release was published, not from the
#                            data it holds
#   the missing's rates      to 48 months, with the missing's band rates as 08b
#   as measured              measures them, rising again after 24 months, in
#                            place of the rates held non-increasing by lag
#   the dead's rates for     to 48 months, with the missing given the dead's
#   the missing beyond       rate of the same band from 24 months on
#   24 months
#
# The last two change the missing's long-lag rates. Pooled over the windows,
# the missing's rates as measured rise again after 24 months while the dead's
# decay: those bands are measured on the 2022-23 cohorts alone, whose missing
# lists were rebuilt when the register imported the Interior Ministry's
# database from February 2025, so they carry that import, not time since the
# event. 08b therefore holds the rates non-increasing by lag (the estimate); the
# first alternative undoes that, the second goes further and gives the missing
# the dead's decay.
#
# Each scenario goes through production's imputation at the central values
# of the evidence on the missing alive, and through the deterministic
# projection of 13b, every other input at its mode. The scenario used must
# reproduce 13b. A grid then crosses the four horizons with five points of
# that evidence - its floor, the ends of the 95% interval of the draws, its
# central values and its ceiling, as 13b's sweep sets them - to show how the
# two assumptions move the results together.
#
# INPUTS   data_inter/ukr_registration_completion.rds (08b)
#          data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds (08)
#          data_inter/ukr_ualosses_resolution_model.rds, ukr_alive_missing.rds (09)
#          data_inter/ukr_alive_sensitivity_military.rds, ukr_alive_sensitivity_e0.rds (13b)
#          the static inputs of 13b
# OUTPUTS  data_inter/ukr_registration_lag.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# 1. COMPLETION FACTORS UNDER EACH HORIZON =====================================
completion <- read_rds("data_inter/ukr_registration_completion.rds")
lambda_band <- completion$lambda_band
L_ref <- completion$L_ref

# the band table carried beyond L_ref at the rate of the last band below it
extend_bands <- function(horizon, lb = lambda_band) {
  if (horizon <= L_ref) return(lb |> filter(hi <= L_ref))
  last <- lb |> filter(hi == L_ref)
  bind_rows(lb |> filter(hi <= L_ref),
            last |> mutate(lo = L_ref, hi = horizon, band = NA))
}
chain_factor <- function(L_from, s, horizon, lb = lambda_band) {
  lb <- extend_bands(horizon, lb) |> filter(status == s)
  exp(sum(pmax(0, pmin(lb$hi, horizon) - pmax(lb$lo, L_from)) * lb$lambda))
}

# the missing's long-lag rates otherwise, the dead's rates as they are:
#   measured   08b's rates as measured, before they are held non-increasing
#   dead_tail  the dead's rate of the same band from 24 months on
lambda_constrained <- list(
  measured = completion$lambda_measured |> arrange(status, lo),
  dead_tail = lambda_band |>
    left_join(lambda_band |> filter(status == "dead") |> select(lo, lambda_dead = lambda), by = "lo") |>
    mutate(lambda = if_else(status == "missing" & lo >= 24, lambda_dead, lambda)) |>
    select(-lambda_dead)
)
# (the dead's rate is above the missing's at 24-30 months, so the second
# alternative raises that band and lowers the three after it)
stopifnot(all(map_lgl(lambda_constrained, \(lb) all(lb$lambda >= -1e-6))),
          all(map_lgl(lambda_constrained, \(lb) isTRUE(all.equal(lb$lambda[lb$status == "dead"], lambda_band$lambda[lambda_band$status == "dead"])))))
cat("\n=== MONTHLY GROWTH OF THE MISSING BY LAG BAND: USED, AS MEASURED AND AT THE DEAD'S RATES (%) ===\n")
print(as.data.frame(
  lambda_band |> filter(status == "missing") |> select(band, lo, hi, used = lambda) |>
    left_join(lambda_constrained$measured |> filter(status == "missing") |> select(lo, measured = lambda), by = "lo") |>
    left_join(lambda_constrained$dead_tail |> filter(status == "missing") |> select(lo, dead_tail = lambda), by = "lo") |>
    mutate(across(c(used, measured, dead_tail), \(x) round(100 * x, 3)))
))

# each event month's completion factor under a horizon (none: 1), under one of
# 08b's alternatives to the way the growth is measured, or under one of the
# constraints on the missing's rates
factors_under <- function(horizon, alternative = NA) {
  if (!is.na(alternative) && alternative %in% names(lambda_constrained)) {
    lb <- lambda_constrained[[alternative]]
    return(completion$month |>
             mutate(f = map2_dbl(lag, status, \(l, s) chain_factor(l, s, horizon, lb))))
  }
  if (!is.na(alternative)) {
    return(completion$month |>
             left_join(completion$alternatives |> filter(alternative == !!alternative) |> select(m, status, f),
                       by = c("m", "status")))
  }
  completion$month |>
    mutate(f = if (is.na(horizon)) 1 else map2_dbl(lag, status, \(l, s) chain_factor(l, s, horizon)))
}
ratio_under <- function(horizon, alternative = NA) {
  factors_under(horizon, alternative) |>
    summarise(n_completed = sum(n_v19 * f), n_v19 = sum(n_v19), .by = c(year, status)) |>
    transmute(year, status, ratio = n_completed / n_v19)
}
# kind: the horizons, which the grid below crosses with the missing alive; the
# alternative specifications of the growth, at the horizon used; and the
# constraints on the missing's long-lag rates, at the horizon used
scenarios <- tibble(
  scenario = c("as registered", "to 48 months (used)", "to 60 months", "to 72 months",
               "to 48 months, growth per person listed", "to 48 months, releases dated by Kaggle version",
               "to 48 months, the missing's rates as measured, rising after 24 months",
               "to 48 months, the dead's rates for the missing beyond 24 months"),
  kind = c(rep("horizon", 4), rep("specification", 2), rep("long_lag", 2)),
  horizon = c(NA, 48, 60, 72, 48, 48, 48, 48),
  alternative = c(rep(NA, 4), "exposure", "dating", "measured", "dead_tail")
)
stopifnot(setequal(scenarios$alternative[scenarios$kind == "specification"], completion$alternatives$alternative),
          setequal(scenarios$alternative[scenarios$kind == "long_lag"], names(lambda_constrained)))
ratios <- scenarios |> mutate(r = map2(horizon, alternative, ratio_under)) |> unnest(r)

# the production ratios are 08b's own
stopifnot(isTRUE(all.equal(
  ratios |> filter(scenario == "to 48 months (used)") |> arrange(status, year) |> pull(ratio),
  completion$year |> arrange(status, year) |> pull(ratio)
)))
# the two leave the dead's ratios as used, and the rates as measured, which
# the estimate only ever lowers, cannot give fewer missing
used <- ratios |> filter(scenario == "to 48 months (used)") |> select(year, status, used = ratio)
stopifnot(ratios |> filter(kind == "long_lag") |> left_join(used, by = c("year", "status")) |>
            summarise(ok = all(if_else(status == "dead", abs(ratio - used) < 1e-9,
                                       alternative != "measured" | ratio >= used - 1e-9))) |> pull(ok))

cat("\n=== COMPLETION RATIOS BY SCENARIO ===\n")
print(as.data.frame(ratios |> mutate(ratio = round(ratio, 4)) |>
                      pivot_wider(names_from = c(status, year), values_from = ratio) |>
                      select(-horizon, -kind, -alternative)))

# 2. THROUGH THE IMPUTATION CHAIN ==============================================
stocks_registered <-
  read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") |>
  summarise(n = sum(dx), .by = c(year, status)) |>
  mutate(status = as.character(status))
mil <- read_rds("data_inter/ukr_military_inputs.rds")
bound <- composition_bound(mil$composition)[1, ]
evidence <- read_rds("data_inter/ukr_alive_missing.rds")
# the scenarios of the missing (00_setup.R): the resolution-based one, the
# estimate, through every horizon; it and the two bounds in the grid below
scenario_args <- function(s) {
  switch(s,
    bound_low = list(captives = 0, other = 1),
    central = list(captives = 0, other = 1 - mil$scenarios$p_central),
    bound_high = list(captives = 0, other = 0),
    bodies = list(captives = 0, other = 1 - unname(mil$scenarios$p_bodies_year[as.character(sort(unique(mil$month$year)))])),
    rule = list(captives = evidence$captives_central, other = bound / 2)
  )
}

# each year's registered count spread over its months as v19 lists them, as
# in 09, and each month completed by its factor under the scenario
military_at <- function(horizon, scenario = "central", alternative = NA) {
  st <- factors_under(horizon, alternative) |>
    mutate(share = n_v19 / sum(n_v19), .by = c(year, status)) |>
    left_join(stocks_registered |> rename(n_year = n), by = c("year", "status")) |>
    mutate(completed = n_year * share * f)
  dead <- st |> filter(status == "dead") |> summarise(confirmed = sum(completed), .by = year)
  miss <- st |> filter(status == "missing") |> select(month = m, year, missing_stock = completed)
  sa <- scenario_args(scenario)
  impute_composition(miss, sa$captives, sa$other, mil$captive_share) |>
    left_join(dead, by = "year") |>
    transmute(year, confirmed, missing = missing_stock, imputed_dead, total = confirmed + imputed_dead)
}
military <-
  scenarios |>
  mutate(res = map2(horizon, alternative, \(h, a) military_at(h, alternative = a))) |>
  select(-alternative) |>
  unnest(res)

# 3. THE DETERMINISTIC PROJECTION, AS IN 13b ===================================
param_table <- read_rds("data_inter/ukr_param_table.rds")
exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") |>
  filter(year %in% 2022:2025, source == "frcst") |> select(-source)
pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") |> filter(reg == "cnt", year == 2022) |> select(-reg)
migs2 <- read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds") |>
  mutate(ems = -mix) |> select(-mix)
asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds")
asfr2 <- bind_rows(asfr, asfr |> filter(year == 2023) |> mutate(year = 2024),
                   asfr |> filter(year == 2023) |> mutate(year = 2025))
ohchr2 <- read_rds("data_inter/ukr_ohchr_civilian_casualties.rds") |> select(year, sex, age, prop_cvs = cx)
ual_raw <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")
cmb_profile <- function(keep_status, nm) {
  ual_raw |> filter(status == keep_status, year %in% 2022:2025) |> select(-status) |>
    complete(year = 2022:2025, sex, age = 0:100, fill = list(dx = 0)) |>
    summarise(dx = sum(dx), .by = c(year, sex, age)) |>
    mutate("{nm}" := dx / sum(dx), .by = year) |> select(-dx)
}
static_inputs <-
  expand_grid(year = 2022:2025, sex = c("f", "m"), age = 0:100) |>
  left_join(migs2, by = c("year", "sex", "age")) |>
  left_join(exp_mort2, by = c("year", "sex", "age")) |>
  left_join(asfr2, by = c("year", "sex", "age")) |>
  left_join(cmb_profile("dead", "prop_cmb_dead") |>
              left_join(cmb_profile("missing", "prop_cmb_miss"), by = c("year", "sex", "age")) |>
              left_join(ohchr2, by = c("year", "sex", "age")), by = c("year", "sex", "age")) |>
  replace_na(list(ems = 0, w = 0, mx = 0, fx = 0, prop_cvs = 0, prop_cmb_dead = 0, prop_cmb_miss = 0))
stopifnot(all(static_inputs$mx > 0), !any(is.na(static_inputs)))

base_draws <-
  central_values(param_table) |> filter(role == "civilians") |> select(year, draw_cvs = mode) |>
  left_join(central_values(param_table) |> filter(str_starts(role, "mig_")) |> summarise(draw_mig = sum(mode), .by = year),
            by = "year")

e0_loss <- function(mil) {
  sim <- run_single_sim(1, base_draws |> left_join(mil |> select(year, draw_cmb = total, conf_cmb = confirmed), by = "year") |>
                          mutate(sim_id = 1), static_inputs, pop22_ini) |>
    mutate(conflict = civilian + combatant_confirmed + combatant_imputed,
           mx_all = (expected + conflict) / pop, mx_bsn = expected / pop)
  e0 <- function(col) sim |> select(year, sex, age, mx = all_of(col)) |> group_by(year, sex) |>
    do(lifetable(dt_in = .data)) |> ungroup() |> filter(age == 0) |> select(year, sex, ex)
  e0("mx_bsn") |> rename(ex_bsn = ex) |> left_join(e0("mx_all") |> rename(ex_all = ex), by = c("year", "sex")) |>
    transmute(year, sex, loss = ex_bsn - ex_all)
}
loss <- military |> nest(.by = c(scenario, kind, horizon)) |>
  mutate(res = map(data, e0_loss)) |> select(-data) |> unnest(res)

# the horizon used must reproduce the central projection
mi <- mode_projection_inputs()
check <- loss_at_mode(mi$draws, mi$static_inputs, mi$pop22_ini) |> arrange(year, sex)
stopifnot(isTRUE(all.equal(loss |> filter(scenario == "to 48 months (used)") |> arrange(year, sex) |> pull(loss),
                           check$loss, tolerance = 1e-6)))

# 4. THE SCENARIOS AND THE HORIZON TOGETHER ===================================
# The military total and the loss under the main text's three - the bounds and
# the resolution-based scenario (00_setup.R) - and the four horizons, to show how
# the two choices move the results jointly. share_alive: the share of the
# missing counted alive under the scenario, at the horizon used.
sc_tab <- read_rds("data_inter/ukr_scenarios.rds")$by_year
scenario_points <-
  sc_tab |>
  filter(scenario %in% scenario_names) |>
  summarise(share_alive = 1 - sum(imputed_dead) / sum(missing), .by = scenario) |>
  rename(point = scenario)
grid <-
  expand_grid(scenarios |> filter(kind == "horizon") |> select(scenario, horizon), scenario_points) |>
  mutate(mil = map2(horizon, point, \(h, p) military_at(h, scenario = p)),
         # summed over rounded years, as table 3 builds its totals, so the two agree
         military = map_dbl(mil, \(m) sum(round(m$total))),
         res = map(mil, e0_loss)) |>
  select(-mil) |>
  unnest(res)
stopifnot(isTRUE(all.equal(
  grid |> filter(scenario == "to 48 months (used)", point == "central") |> arrange(year, sex) |> pull(loss),
  loss |> filter(scenario == "to 48 months (used)") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(lambda_band = lambda_band, L_ref = L_ref, ratios = ratios,
               military = military, loss = loss, grid = grid),
          "data_inter/ukr_registration_lag.rds")

cat("\n=== MILITARY TOTAL AND 2025 MALE LOSS, THE MISSING ALIVE BY HORIZON ===\n")
print(as.data.frame(
  grid |> filter(year == 2025, sex == "m") |>
    transmute(scenario, point, cell = sprintf("%s / %.2f", scales::comma(round(military)), loss)) |>
    pivot_wider(names_from = scenario, values_from = cell)
))

cat("\n=== MILITARY TOTAL BY SCENARIO, AT THE CENTRAL VALUES ===\n")
print(as.data.frame(military |> summarise(across(c(confirmed, missing, imputed_dead, total), \(x) round(sum(x))),
                                          .by = scenario)))
cat("\n=== LOSS, MEN ===\n")
print(as.data.frame(loss |> filter(sex == "m") |> mutate(loss = round(loss, 2)) |> select(-horizon, -kind) |>
                      pivot_wider(names_from = year, values_from = loss)))

message("Done. data_inter/ukr_registration_lag.rds written.")
