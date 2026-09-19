# ==============================================================================
# STEP 13d - Sensitivity to the registration-lag correction
# ==============================================================================
#
# WHY
# ---
# 08b brings v19's dead and missing to the completeness events reach at four
# years, and 09 imputes on those completed counts. Four years is where the
# register's growth is last measured on enough event months, not where it
# stops: the dead still grow by about 0.2% a month at 42-48 months. So the
# correction is a lower bound on the lag. This step shows how far the results
# move without the correction, and with the growth extrapolated beyond four
# years.
#
# SCENARIOS
# ---------
#   as registered            v19's counts, no correction
#   to 48 months (used)      08b's chain-ladder, the production estimate
#   to 60 / 72 months        the same, with the rate of the 42-48 month band
#                            carried on beyond 48 months
#
# Each scenario goes through production's imputation chain at the central
# alpha its own evidence gives, and through the deterministic projection of
# 13b, every other input at its mode. The scenario used must reproduce 13b.
#
# INPUTS   data_inter/ukr_registration_completion.rds (08b)
#          data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds (08)
#          data_inter/ukr_ualosses_resolution_model.rds, ukr_alpha_missing.rds (09)
#          data_inter/ukr_alpha_sensitivity_e0.rds (13b)
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
extend_bands <- function(horizon) {
  if (horizon <= L_ref) return(lambda_band |> filter(hi <= L_ref))
  last <- lambda_band |> filter(hi == L_ref)
  bind_rows(lambda_band |> filter(hi <= L_ref),
            last |> mutate(lo = L_ref, hi = horizon, band = NA))
}
chain_factor <- function(L_from, s, horizon) {
  lb <- extend_bands(horizon) |> filter(status == s)
  exp(sum(pmax(0, pmin(lb$hi, horizon) - pmax(lb$lo, L_from)) * lb$lambda))
}
# each event month's completion factor under a horizon (none: 1)
factors_under <- function(horizon) {
  completion$month |>
    mutate(f = if (is.na(horizon)) 1 else map2_dbl(lag, status, \(l, s) chain_factor(l, s, horizon)))
}
ratio_under <- function(horizon) {
  factors_under(horizon) |>
    summarise(n_completed = sum(n_v19 * f), n_v19 = sum(n_v19), .by = c(year, status)) |>
    transmute(year, status, ratio = n_completed / n_v19)
}
scenarios <- tibble(
  scenario = c("as registered", "to 48 months (used)", "to 60 months", "to 72 months"),
  horizon = c(NA, 48, 60, 72)
)
ratios <- scenarios |> mutate(r = map(horizon, ratio_under)) |> unnest(r)

# the production ratios are 08b's own
stopifnot(isTRUE(all.equal(
  ratios |> filter(scenario == "to 48 months (used)") |> arrange(status, year) |> pull(ratio),
  completion$year |> arrange(status, year) |> pull(ratio)
)))

cat("\n=== COMPLETION RATIOS BY SCENARIO ===\n")
print(as.data.frame(ratios |> mutate(ratio = round(ratio, 4)) |>
                      pivot_wider(names_from = c(status, year), values_from = ratio) |> select(-horizon)))

# 2. THROUGH THE IMPUTATION CHAIN ==============================================
stocks_registered <-
  read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") |>
  summarise(n = sum(dx), .by = c(year, status)) |>
  mutate(status = as.character(status))
register_alive <- sum(stocks_registered$n[stocks_registered$status %in% c("prisoner", "released_prisoner")])
model <- read_rds("data_inter/ukr_ualosses_resolution_model.rds")
alpha_used <- read_rds("data_inter/ukr_alpha_missing.rds")

# each year's registered count spread over its months as v19 lists them, as
# in 09, and each month completed by its factor under the horizon
military_under <- function(horizon) {
  st <- factors_under(horizon) |>
    mutate(share = n_v19 / sum(n_v19), .by = c(year, status)) |>
    left_join(stocks_registered |> rename(n_year = n), by = c("year", "status")) |>
    mutate(completed = n_year * share * f)
  dead <- st |> filter(status == "dead") |> summarise(confirmed = sum(completed), .by = year)
  miss <- st |> filter(status == "missing") |> select(month = m, year, missing_stock = completed)
  impute_fn <- \(a) impute_missing(a, model, miss)
  a <- alpha_evidence(impute_fn, model$resolved, register_alive,
                      net_projected_captivity = alpha_used$net_projected_captivity)
  impute_fn(a$alpha_mode) |>
    left_join(dead, by = "year") |>
    transmute(year, alpha = a$alpha_mode, confirmed, missing = missing_stock, imputed_dead,
              total = confirmed + imputed_dead)
}
military <- scenarios |> mutate(res = map(horizon, military_under)) |> unnest(res)

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
  param_table |> filter(role == "civilians") |> select(year, draw_cvs = mode) |>
  left_join(param_table |> filter(str_starts(role, "mig_")) |> summarise(draw_mig = sum(mode), .by = year),
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
loss <- military |> nest(.by = c(scenario, horizon)) |>
  mutate(res = map(data, e0_loss)) |> select(-data) |> unnest(res)

# the scenario used must reproduce 13b at the central alpha
check <- read_rds("data_inter/ukr_alpha_sensitivity_e0.rds") |> filter(range_point %in% "mode") |> arrange(year, sex)
stopifnot(isTRUE(all.equal(loss |> filter(scenario == "to 48 months (used)") |> arrange(year, sex) |> pull(loss),
                           check$loss, tolerance = 1e-8)))

write_rds(list(lambda_band = lambda_band, L_ref = L_ref, ratios = ratios,
               military = military, loss = loss),
          "data_inter/ukr_registration_lag.rds")

cat("\n=== MILITARY TOTAL BY SCENARIO, EACH AT ITS OWN CENTRAL ALPHA ===\n")
print(as.data.frame(military |> summarise(alpha = first(alpha),
                                          across(c(confirmed, missing, imputed_dead, total), \(x) round(sum(x))),
                                          .by = scenario)))
cat("\n=== LOSS, MEN ===\n")
print(as.data.frame(loss |> filter(sex == "m") |> mutate(loss = round(loss, 2)) |> select(-horizon) |>
                      pivot_wider(names_from = year, values_from = loss)))

message("Done. data_inter/ukr_registration_lag.rds written.")
