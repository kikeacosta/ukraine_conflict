# ==============================================================================
# STEP 13g - Sensitivity to the counterfactual window and the PERT shape
# ==============================================================================
#
# 1. THE LEE-CARTER WINDOW
# ------------------------
# The counterfactual is a Lee-Carter forecast fitted to 2000-2019 (04). Other
# defensible windows start earlier, in the post-Soviet recovery, or later, or
# keep the COVID years in. Each is fitted as 04 fits its own, forecast to
# 2025, and put through the deterministic projection at the central value of
# every conflict and migration input. The window used must reproduce 13's loss.
# The forecast's own uncertainty, which 11 carries by drawing the forecast
# index, is shown here on its own as well: the rates one standard deviation
# of the index below and above its mean.
#
# 2. THE SHAPE OF THE PERT DISTRIBUTIONS
# --------------------------------------
# Every disputed input - civilian deaths, the three inputs on the missing
# alive (the share of the unrecorded prisoners among the missing, the
# prisoners held, the share of the unresolved alive for other reasons), the
# western blend weight, Russia and Belarus - is a Beta-PERT with shape 4, the
# standard choice. A flatter shape (2) puts more weight towards the bounds, a
# sharper one (8) less. The simulation is re-run at each shape with the same
# random numbers (simulation_draws(), 00_setup.R), so the differences in the
# intervals come from the shape alone.
#
# INPUTS   data_inter/ukr_life_tables_1989_2021.rds (03), data_input/DataDxEx.csv
#          the static inputs at the mode; the parameter table (10), the
#          military inputs (09) and the forecast error (04)
# OUTPUTS  data_inter/ukr_baseline_window_e0.rds, ukr_forecast_spread_by_window.rds
#          data_inter/ukr_pert_shape_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

# 1. LEE-CARTER WINDOWS ========================================================
mx_all <- read_rds("data_inter/ukr_life_tables_1989_2021.rds")
exposures <-
  read.csv("data_input/DataDxEx.csv", header = TRUE) |>
  rename_with(tolower) |>
  filter(data == "Ex") |>
  pivot_longer(starts_with("age"), names_to = "age", values_to = "pop") |>
  mutate(age = as.numeric(gsub("age", "", age)), sex = str_sub(sex, 1, 1)) |>
  select(year, sex, age, pop)
dt <-
  mx_all |>
  select(year, sex, age, mx) |>
  mutate(sex = str_sub(sex, 1, 1)) |>
  left_join(exposures, by = c("year", "sex", "age")) |>
  mutate(deaths = mx * pop) |>
  drop_na(mx, pop)

# as 04: adjusted on e0, jumping off the FITTED rates of the last year. The
# Lee-Miller variant jumps off the OBSERVED rates instead, which is the last row
# below: it moves the level of the counterfactual more than it moves the loss.
forecast_window <- function(years, jump = "fit") {
  dt |>
    filter(year %in% years) |>
    as_vital(index = year, key = c(sex, age), .age = "age", .sex = "sex",
             .deaths = "deaths", .population = "pop") |>
    model(lc = LC(log(mx), adjust = "e0", jump_choice = jump)) |>
    forecast(h = 2025 - max(years)) |>
    as_tibble() |>
    select(year, sex, age, mx = .mean) |>
    filter(year %in% 2022:2025)
}
windows <- tibble(
  window = c("2000-2019 (used)", "1995-2019", "2005-2019", "2010-2019",
             "2000-2021, with the COVID years",
             "2000-2019, jumping off the observed rates (Lee-Miller)"),
  first = c(2000, 1995, 2005, 2010, 2000, 2000),
  last = c(2019, 2019, 2019, 2019, 2021, 2019),
  jump = c("fit", "fit", "fit", "fit", "fit", "actual")
)
with_mx <- function(fc) mi$static_inputs |> select(-mx) |> left_join(fc, by = c("year", "sex", "age"))

baseline <-
  windows |>
  mutate(res = pmap(list(first, last, jump), function(f, l, j) {
    fc <- forecast_window(f:l, j)
    st <- with_mx(fc)
    stopifnot(all(st$mx > 0), !any(is.na(st$mx)))
    loss_at_mode(mi$draws, st, mi$pop22_ini)
  })) |>
  unnest(res)

# the window used must reproduce 04's forecast and 13's loss
stopifnot(isTRUE(all.equal(
  baseline |> filter(first == 2000, last == 2019, jump == "fit") |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-6
)))

# The forecast's own uncertainty on its own: the rates one standard deviation
# of the forecast index k below and above its mean, window 2000-2019. 11
# draws the index, so its intervals carry this. log(mx) is linear in k, so the 15.9% and 84.1%
# quantiles of each rate's forecast distribution are the rates at k -/+ 1 SD.
fc_used <-
  dt |>
  filter(year %in% 2000:2019) |>
  as_vital(index = year, key = c(sex, age), .age = "age", .sex = "sex",
           .deaths = "deaths", .population = "pop") |>
  model(lc = LC(log(mx), adjust = "e0", jump_choice = "fit")) |>
  forecast(h = 6) |>
  as_tibble() |>
  filter(year %in% 2022:2025)
k_band <-
  tibble(window = c("2000-2019, forecast mortality one SD lower",
                    "2000-2019, forecast mortality one SD higher"),
         p = c(pnorm(-1), pnorm(1))) |>
  mutate(res = map(p, function(pp) {
    fc <- fc_used |> transmute(year, sex, age, mx = quantile(mx, pp))
    st <- with_mx(fc)
    stopifnot(all(st$mx > 0), !any(is.na(st$mx)))
    loss_at_mode(mi$draws, st, mi$pop22_ini)
  })) |>
  select(-p) |>
  unnest(res)
baseline <- bind_rows(baseline, k_band)
write_rds(baseline, "data_inter/ukr_baseline_window_e0.rds")

# The forecast's spread by window: each window's forecast rates 1.96 standard
# deviations of its index either way, as the counterfactual e0 and the loss in
# each year, every other input at its mode. This is the width the drawn
# forecast gives the intervals, window by window: a window that takes in the
# 2005-2008 crisis has a more erratic index, so a wider band.
spread <-
  windows |>
  filter(jump == "fit") |>   # the jump-off is not a window; it has no band of its own
  mutate(res = map2(first, last, function(f, l) {
    fc <-
      dt |>
      filter(year %in% f:l) |>
      as_vital(index = year, key = c(sex, age), .age = "age", .sex = "sex",
               .deaths = "deaths", .population = "pop") |>
      model(lc = LC(log(mx), adjust = "e0", jump_choice = "fit")) |>
      forecast(h = 2025 - l) |>
      as_tibble() |>
      filter(year %in% 2022:2025)
    map_dfr(c(lower = 0.025, upper = 0.975), function(pp) {
      st <- with_mx(fc |> transmute(year, sex, age, mx = quantile(mx, pp)))
      stopifnot(all(st$mx > 0), !any(is.na(st$mx)))
      loss_at_mode(mi$draws, st, mi$pop22_ini)
    }, .id = "bound")
  })) |>
  unnest(res)
write_rds(spread, "data_inter/ukr_forecast_spread_by_window.rds")
cat("\n=== THE FORECAST'S 95% BAND BY WINDOW, 2025 ===\n")
print(as.data.frame(
  spread |> filter(year == 2025) |>
    mutate(across(c(e0_bsn, loss), \(x) round(x, 2))) |>
    select(window, bound, sex, e0_bsn, loss) |>
    pivot_wider(names_from = c(sex, bound), values_from = c(e0_bsn, loss))
))

cat("\n=== COUNTERFACTUAL e0 AND LOSS BY LEE-CARTER WINDOW ===\n")
print(as.data.frame(
  baseline |> filter(year == 2025) |>
    mutate(across(c(e0_bsn, loss), \(x) round(x, 2))) |>
    select(window, sex, e0_bsn, loss) |>
    pivot_wider(names_from = sex, values_from = c(e0_bsn, loss))
))

# 2. PERT SHAPES ===============================================================
param_table <- read_rds("data_inter/ukr_param_table.rds")
mil <- read_rds("data_inter/ukr_military_inputs.rds")
lc_error <- read_rds("data_inter/ukr_lc_forecast_error.rds")
# 1,000 draws per shape at the production size; fewer in a quick test run
n_shape <- min(1000, n_sim)
# as 11: the draws carry their own counterfactual, so the static inputs carry
# the forecast's loading and variance
static_shape <-
  mi$static_inputs |>
  left_join(lc_error$bx, by = c("sex", "age")) |>
  left_join(lc_error$variance |> select(sex, year, vk), by = c("sex", "year"))

loss_by_draw <- function(draws_df) {
  sims <- map_dfr(split(draws_df, draws_df$sim_id), function(d) {
    run_single_sim(d$sim_id[1], d, static_shape, mi$pop22_ini)
  })
  d <- as.data.table(sims)
  setorder(d, sim_id, year, sex, age)
  n_age <- 101
  key <- d[seq(1L, .N, by = n_age), .(sim_id, year, sex)]
  POP <- matrix(d$pop, nrow = n_age)
  DEX <- matrix(d$expected, nrow = n_age)
  DCF <- matrix(d$civilian + d$combatant_confirmed + d$combatant_imputed, nrow = n_age)
  key[, loss := lt_cols(DEX / POP, key$sex)$ex[1, ] - lt_cols((DEX + DCF) / POP, key$sex)$ex[1, ]]
  as_tibble(key)
}

shapes <- c(2, 4, 8)
pert_shape <-
  map_dfr(shapes, function(s) {
    draws <- simulation_draws(param_table, mil, lc_error, n = n_shape, shape = s)$draws_df
    loss_by_draw(draws) |>
      summarise(mean = mean(loss), lo = quantile(loss, 0.025), hi = quantile(loss, 0.975),
                .by = c(year, sex)) |>
      mutate(shape = s)
  }) |>
  mutate(width = hi - lo)
write_rds(pert_shape, "data_inter/ukr_pert_shape_e0.rds")

cat("\n=== LOSS BY PERT SHAPE: MEAN AND 95% INTERVAL ===\n")
print(as.data.frame(pert_shape |> mutate(across(c(mean, lo, hi, width), \(x) round(x, 2))) |>
                      arrange(sex, year, shape)))
message("Done. data_inter/ukr_baseline_window_e0.rds and ukr_pert_shape_e0.rds written.")
