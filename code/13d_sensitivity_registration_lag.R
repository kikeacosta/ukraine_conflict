# ==============================================================================
# STEP 13d - Sensitivity of the military total to registration lag
# ==============================================================================
#
# WHY
# ---
# The register keeps adding deaths and disappearances long after the event:
# for 2025 events, v14 (16 Sep 2025) listed 9,001 dead, v19 (19 Sep 2026)
# 16,457, and even 2022 events keep growing. The imputation only reclassifies
# people already listed, so a loss not yet registered is not counted.
#
# METHOD - a chain-ladder on the four releases
# --------------------------------------------
# For each event month m, each status s (dead, missing) and each consecutive
# pair of releases r -> r+1:
#   net additions = records listed with status s in r+1 whose key is absent
#                   from r in any status
#                   - records with status s in r whose key is absent from r+1
#                   in any status
# A key present in both releases is left out whatever its status, so a missing
# person who is found dead, taken prisoner or released is a resolution, not a
# registration or a de-registration. The monthly growth rate
#   lambda = log(1 + net additions / N_r) / (months between the releases)
# is attached to the lag since the event (midpoint of the pair), and pooled
# by lag band, weighting by N_r x months.
#
# Chaining those rates from each month's lag at v19 up to 48 months, the
# longest lag measured on enough event months (section 3), gives a factor that
# brings every month to the completeness events reach at four years. The
# register still grows beyond that lag, and that growth is not corrected: the
# factor is a lower bound on the lag, not a full correction.
#
# The corrected dead and missing stocks then go through production's own
# imputation chain at the central alpha and the deterministic projection
# of 13b, every other input at its mode. This is a sensitivity; production
# does not apply it.
#
# INPUTS   the four register releases (gitignored), through a cache of
#          aggregate counts: data_inter/ualosses_registration_by_month.rds
#          data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds (08)
#          data_inter/ukr_ualosses_transition_rates.rds, ukr_alpha_missing.rds (09)
#          the static inputs of 13b
# OUTPUTS  data_inter/ukr_registration_lag.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

releases <- ual_releases

# 1. AGGREGATE COUNTS FROM THE RELEASES (cached; no names leave this block) ====
reg_counts <- cache_rds("data_inter/ualosses_registration_by_month.rds", {
  rd <- function(p) {
    ual_read_release(p) |>
      mutate(m = floor_date(date_evnt, "month")) |>
      filter(ukrainian, !is.na(m), m >= as.Date("2022-02-01"), m <= as.Date("2025-12-01")) |>
      ual_one_per_key() |>
      select(key, m, status)
  }
  d <- set_names(map(releases$path, rd), releases$release)
  losses <- \(x) filter(x, status %in% c("dead", "missing"))

  # one long table: stocks by release, and additions and drops between
  # consecutive releases. A record counts as added or dropped only if its key
  # is absent from the other release in ANY status, so a missing person who
  # becomes a prisoner or is released is a resolution, not a de-registration.
  stock <- imap_dfr(d, \(x, r) count(losses(x), m, status, name = "n") |>
                      mutate(kind = "stock", from = r, to = NA_character_))
  flows <- map_dfr(1:3, function(i) {
    a <- d[[i]]; b <- d[[i + 1]]
    bind_rows(
      losses(b) |> anti_join(a, by = "key") |> count(m, status, name = "n") |> mutate(kind = "added"),
      losses(a) |> anti_join(b, by = "key") |> count(m, status, name = "n") |> mutate(kind = "dropped")
    ) |> mutate(from = releases$release[i], to = releases$release[i + 1])
  })
  bind_rows(stock, flows)
})

stock <- reg_counts |> filter(kind == "stock") |> transmute(m, status, n, release = from)
flows <- reg_counts |> filter(kind != "stock") |> select(m, status, n, kind, from, to)

# 2. GROWTH RATES BY LAG ======================================================
lag_months <- function(release_date, m) as.numeric(release_date - (m + 14)) / 30.44
rel_date <- set_names(releases$date, releases$release)

growth <-
  flows |>
  pivot_wider(names_from = kind, values_from = n, values_fill = 0) |>
  left_join(stock |> rename(from = release, n_from = n), by = c("m", "status", "from")) |>
  mutate(
    L0 = lag_months(rel_date[from], m), L1 = lag_months(rel_date[to], m),
    dL = L1 - L0, Lmid = (L0 + L1) / 2,
    lambda = log1p((added - dropped) / n_from) / dL
  ) |>
  filter(L0 >= 1, n_from >= 50)

bands <- c(1, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48, 56)
lambda_band <-
  growth |>
  mutate(band = cut(Lmid, bands, right = TRUE)) |>
  filter(!is.na(band)) |>
  summarise(lambda = weighted.mean(lambda, n_from * dL), n_obs = n(), .by = c(status, band)) |>
  mutate(lo = bands[as.integer(band)], hi = bands[as.integer(band) + 1]) |>
  arrange(status, lo)

cat("\n=== MONTHLY GROWTH OF THE REGISTER BY LAG SINCE THE EVENT ===\n")
print(as.data.frame(lambda_band |> mutate(lambda_pct = round(100 * lambda, 3)) |>
                      select(status, band, lambda_pct, n_obs)))

# 3. CHAIN TO THE LONGEST WELL-MEASURED LAG ===================================
# The chain runs only as far as the rates rest on at least a year's worth of
# event months (12). The last band, 48-56 months, rests on the five event
# months of February-June 2022 in a single pair of releases, and there v19
# deleted a batch of February 2022 missing records in one clean-up: a one-off,
# not a registration lag.
L_ref <- lambda_band |>
  summarise(ok = all(n_obs >= 12), .by = c(lo, hi)) |>
  filter(ok) |>
  pull(hi) |>
  max()
stopifnot(L_ref == 48)
chain_factor <- function(L_from, s) {
  lb <- lambda_band |> filter(status == s)
  exp(sum(pmax(0, pmin(lb$hi, L_ref) - pmax(lb$lo, L_from)) * lb$lambda))
}

v19_month <- stock |> filter(release == "v19") |> mutate(L19 = lag_months(rel_date[["v19"]], m))
completion <-
  v19_month |>
  mutate(factor = map2_dbl(L19, status, chain_factor), n_corrected = n * factor) |>
  mutate(year = year(m)) |>
  summarise(n_v19 = sum(n), n_corrected = sum(n_corrected), .by = c(year, status)) |>
  mutate(ratio = n_corrected / n_v19) |>
  arrange(status, year)

cat("\n=== V19 COUNTS BROUGHT TO THE COMPLETENESS OF THE OLDEST EVENTS ===\n")
print(as.data.frame(completion |> mutate(across(c(n_v19, n_corrected), round), ratio = round(ratio, 4))))

# 4. THROUGH THE IMPUTATION CHAIN AND THE PROJECTION ===========================
ual_raw <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")
stocks <- ual_raw |> summarise(n = sum(dx), .by = c(year, status)) |> mutate(status = as.character(status))
tasas_long <- read_rds("data_inter/ukr_ualosses_transition_rates.rds")
alpha_mode <- read_rds("data_inter/ukr_alpha_missing.rds")$alpha_mode
ratio_of <- function(s) completion |> filter(status == s) |> select(year, ratio)

military_under <- function(scale) {
  dead <- stocks |> filter(status == "dead") |> select(year, dead = n) |>
    left_join(if (scale) ratio_of("dead") else tibble(year = 2022:2025, ratio = 1), by = "year") |>
    mutate(dead = dead * ratio) |> select(-ratio)
  miss <- stocks |> filter(status == "missing") |> select(year, missing_stock = n) |>
    left_join(if (scale) ratio_of("missing") else tibble(year = 2022:2025, ratio = 1), by = "year") |>
    mutate(missing_stock = missing_stock * ratio) |> select(-ratio)
  impute_missing(alpha_mode, tasas_long, miss) |>
    left_join(dead, by = "year") |>
    transmute(year, confirmed = dead, missing = missing_stock, imputed_dead,
              total = dead + imputed_dead)
}
military <- bind_rows(
  military_under(FALSE) |> mutate(scenario = "v19 as registered"),
  military_under(TRUE) |> mutate(scenario = "corrected for registration lag")
)

# deterministic projection, as in 13b
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
loss <- military |> nest(.by = scenario) |> mutate(res = map(data, e0_loss)) |> select(-data) |> unnest(res)

# the uncorrected scenario must reproduce 13b at the central alpha
check <- read_rds("data_inter/ukr_alpha_sensitivity_e0.rds") |> filter(range_point == "mode") |> arrange(year, sex)
stopifnot(isTRUE(all.equal(loss |> filter(scenario == "v19 as registered") |> arrange(year, sex) |> pull(loss),
                           check$loss, tolerance = 1e-8)))

write_rds(list(lambda_band = lambda_band, completion = completion, L_ref = L_ref,
               military = military, loss = loss),
          "data_inter/ukr_registration_lag.rds")

cat("\n=== MILITARY TOTAL AT THE CENTRAL ALPHA ===\n")
print(as.data.frame(military |> summarise(across(c(confirmed, missing, imputed_dead, total), \(x) round(sum(x))), .by = scenario)))
cat("\n=== LOSS, MEN ===\n")
print(as.data.frame(loss |> filter(sex == "m") |> mutate(loss = round(loss, 2)) |>
                      pivot_wider(names_from = year, values_from = loss)))

message("Done. data_inter/ukr_registration_lag.rds written.")
