# ==============================================================================
# STEP 13n - Pandemic mortality carried into the 2022 counterfactual
# ==============================================================================
#
# WHY
# ---
# The counterfactual is fitted to 2000-2019 and forecast from there, so it holds
# no pandemic mortality at all (04). 2021 was Ukraine's worst pandemic year, and
# in the neighbouring countries excess mortality did not end with it: 2022 still
# ran above trend, by a fraction of 2021's excess. Without the war some of that
# would have been Ukraine's too, so the 2022 counterfactual is, if anything, too
# low, and the 2022 loss too high. This is one of the few biases of the estimate
# that run upwards.
#
# No source says how much would have carried over. The scenarios take a quarter
# and a half of 2021's excess - the observed rates of 2021 over the forecast for
# 2021, by age and sex, on the log scale, smoothed over eleven years of age -
# and raise the 2022 counterfactual by it; 2023-2025 are left as forecast. Every
# other input is at its mode. The excess measured this way is the pandemic's
# only if the forecast is right for 2021, which is what the counterfactual
# assumes throughout.
#
# INPUTS   data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds (04), the static
#          inputs at the mode
# OUTPUTS  data_inter/ukr_covid_carryover_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

mxs <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") |> mutate(mx = unname(mx))
excess_2021 <-
  mxs |>
  filter(year == 2021) |>
  pivot_wider(names_from = source, values_from = mx) |>
  arrange(sex, age) |>
  mutate(log_excess = log(obs / frcst),
         log_excess = as.numeric(stats::filter(log_excess, rep(1 / 11, 11), sides = 2)),
         .by = sex) |>
  # the ends of the age range, where the window runs out, take the nearest value
  group_by(sex) |>
  fill(log_excess, .direction = "downup") |>
  ungroup() |>
  select(sex, age, log_excess)

cat("\n=== 2021: OBSERVED OVER FORECAST DEATH RATES, SMOOTHED ===\n")
print(as.data.frame(
  excess_2021 |> filter(age %in% seq(0, 90, 10)) |>
    mutate(ratio = round(exp(log_excess), 3)) |> select(-log_excess) |>
    pivot_wider(names_from = sex, values_from = ratio)
))

loss_with <- function(fraction) {
  st <- mi$static_inputs |>
    left_join(excess_2021, by = c("sex", "age")) |>
    mutate(mx = if_else(year == 2022, mx * exp(fraction * log_excess), mx)) |>
    select(-log_excess)
  stopifnot(all(st$mx > 0), !any(is.na(st$mx)))
  loss_at_mode(mi$draws, st, mi$pop22_ini)
}
scenarios <- tibble(
  scenario = c("None (as estimated)", "A quarter of 2021's excess carried into 2022",
               "Half of 2021's excess carried into 2022"),
  fraction = c(0, 0.25, 0.5)
)
loss <- scenarios |> mutate(res = map(fraction, loss_with)) |> unnest(res)

stopifnot(isTRUE(all.equal(
  loss |> filter(fraction == 0) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(excess_2021 = excess_2021, loss = loss), "data_inter/ukr_covid_carryover_e0.rds")

cat("\n=== LOSS AND COUNTERFACTUAL e0 WITH PANDEMIC MORTALITY CARRIED INTO 2022 ===\n")
print(as.data.frame(
  loss |> filter(year == 2022) |>
    transmute(scenario, sex, loss = round(loss, 3), e0_bsn = round(e0_bsn, 2)) |>
    pivot_wider(names_from = sex, values_from = c(loss, e0_bsn))
))
message("Done. data_inter/ukr_covid_carryover_e0.rds written.")
