# ==============================================================================
# STEP 13g - Sensitivity to the counterfactual window and the PERT shape
# ==============================================================================
#
# 1. THE LEE-CARTER WINDOW
# ------------------------
# The counterfactual is a Lee-Carter forecast fitted to 2000-2019 (04). Other
# defensible windows start earlier, in the post-Soviet recovery, or later, or
# keep the COVID years in. Each is fitted as 04 fits its own, forecast to
# 2025, and put through the deterministic projection at the mode of every
# conflict and migration input. The window used must reproduce 13's loss.
# The forecast's own uncertainty is shown as well: the rates one standard
# deviation of the forecast index below and above its mean.
#
# 2. THE SHAPE OF THE PERT DISTRIBUTIONS
# --------------------------------------
# Every drawn input - civilian deaths, alpha, the western blend weight, Russia
# and Belarus - is a Beta-PERT with shape 4, the standard choice. A flatter
# shape (2) puts more weight towards the bounds, a sharper one (8) less. The
# simulation is re-run at each shape with the same random numbers
# (simulation_draws(), 00_setup.R), so the differences in the intervals come
# from the shape alone.
#
# INPUTS   data_inter/ukr_life_tables_1989_2021.rds (03), data_input/DataDxEx.csv
#          the static inputs at the mode; the parameter table, the range of
#          alpha and the military lines (09, 10)
# OUTPUTS  data_inter/ukr_baseline_window_e0.rds
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

# as 04: Lee-Miller variant, adjusted on e0, jumping off the fitted rates
forecast_window <- function(years) {
  dt |>
    filter(year %in% years) |>
    as_vital(index = year, key = c(sex, age), .age = "age", .sex = "sex",
             .deaths = "deaths", .population = "pop") |>
    model(lc = LC(log(mx), adjust = "e0", jump_choice = "fit")) |>
    forecast(h = 2025 - max(years)) |>
    as_tibble() |>
    select(year, sex, age, mx = .mean) |>
    filter(year %in% 2022:2025)
}
windows <- tibble(
  window = c("2000-2019 (used)", "1995-2019", "2005-2019", "2010-2019", "2000-2021, with the COVID years"),
  first = c(2000, 1995, 2005, 2010, 2000),
  last = c(2019, 2019, 2019, 2019, 2021)
)
with_mx <- function(fc) mi$static_inputs |> select(-mx) |> left_join(fc, by = c("year", "sex", "age"))

baseline <-
  windows |>
  mutate(res = map2(first, last, function(f, l) {
    fc <- forecast_window(f:l)
    st <- with_mx(fc)
    stopifnot(all(st$mx > 0), !any(is.na(st$mx)))
    loss_at_mode(mi$draws, st, mi$pop22_ini)
  })) |>
  unnest(res)

# the window used must reproduce 04's forecast and 13's loss
stopifnot(isTRUE(all.equal(
  baseline |> filter(first == 2000, last == 2019) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-6
)))

# The forecast's own uncertainty, which the simulation does not carry: the
# rates one standard deviation of the forecast index k below and above its
# mean, window 2000-2019. log(mx) is linear in k, so the 15.9% and 84.1%
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

cat("\n=== COUNTERFACTUAL e0 AND LOSS BY LEE-CARTER WINDOW ===\n")
print(as.data.frame(
  baseline |> filter(year == 2025) |>
    mutate(across(c(e0_bsn, loss), \(x) round(x, 2))) |>
    select(window, sex, e0_bsn, loss) |>
    pivot_wider(names_from = sex, values_from = c(e0_bsn, loss))
))

# 2. PERT SHAPES ===============================================================
param_table <- read_rds("data_inter/ukr_param_table.rds")
alpha_range <- read_rds("data_inter/ukr_alpha_missing.rds")
alpha_lines <- read_rds("data_inter/ukr_military_alpha_lines.rds")
n_shape <- 1000

loss_by_draw <- function(draws_df) {
  sims <- map_dfr(split(draws_df, draws_df$sim_id), function(d) {
    run_single_sim(d$sim_id[1], d, mi$static_inputs, mi$pop22_ini)
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
    draws <- simulation_draws(param_table, alpha_range, alpha_lines, n = n_shape, shape = s)$draws_df
    loss_by_draw(draws) |>
      summarise(median = median(loss), lo = quantile(loss, 0.025), hi = quantile(loss, 0.975),
                .by = c(year, sex)) |>
      mutate(shape = s)
  }) |>
  mutate(width = hi - lo)
write_rds(pert_shape, "data_inter/ukr_pert_shape_e0.rds")

cat("\n=== LOSS BY PERT SHAPE: MEDIAN AND 95% INTERVAL ===\n")
print(as.data.frame(pert_shape |> mutate(across(c(median, lo, hi, width), \(x) round(x, 2))) |>
                      arrange(sex, year, shape)))
message("Done. data_inter/ukr_baseline_window_e0.rds and ukr_pert_shape_e0.rds written.")
