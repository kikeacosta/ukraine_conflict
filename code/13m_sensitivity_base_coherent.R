# ==============================================================================
# STEP 13m - A smaller population base, with the counterfactual kept coherent
# ==============================================================================
#
# WHY
# ---
# SSSU's 1 January 2022 population carries the 2001 census forward with
# registered migration only, so it still holds people who had left without
# deregistering, most of them of working age. 13e makes the base outside Donetsk
# and Luhansk 0.5 and 1.0 million smaller at ages 20-64 and leaves the
# counterfactual alone. But the counterfactual is fitted to death rates whose
# exposures come from the same series for the same territory - the pre-war
# series covers the regions with complete registration, which is continental
# Ukraine without Donetsk and Luhansk - so if the base is too large, so were the
# exposures, by an amount that grew from nothing at the 2001 census. Rates
# computed on smaller exposures are higher, the forecast of the counterfactual is
# higher with them, and part of what the smaller base adds to the loss is taken
# back.
#
# This step makes the two coherent. For a base X million smaller at ages 20-64
# outside Donetsk and Luhansk, the exposures of the Lee-Carter window at those
# ages are made smaller by a share that rises linearly from zero in 2001 to the
# share X is of the 2022 base, the rates are recomputed on them, the model is
# fitted again as 04 and 13g fit it, and the loss is projected on the smaller
# base with the new counterfactual. X runs to 2.0 million, beyond 13e's range:
# the 2019 electronic census and SSSU differ by under a million on a like-for-
# like territory, but other assessments put the gap higher.
#
# Donetsk and Luhansk need no such step: they are not in the pre-war series, so
# their base moves the conflict rates alone (13e, 14d).
#
# INPUTS   data_inter/ukr_life_tables_1989_2021.rds (03), data_input/DataDxEx.csv,
#          data_inter/ukr_pop_sssu.rds (01), the static inputs at the mode,
#          data_inter/ukr_denominator_donetsk_luhansk.rds (13e)
# OUTPUTS  data_inter/ukr_base_coherent_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

# 1. THE BASE, AS 13e CUTS IT ==================================================
pop_dl <-
  read_rds("data_inter/ukr_pop_sssu.rds") |>
  filter(reg %in% c("dnk", "luk"), year == 2022) |>
  summarise(pop_dl = sum(pop), .by = c(sex, age))
working_rest <-
  mi$pop22_ini |>
  left_join(pop_dl, by = c("sex", "age")) |>
  mutate(rest = pop - coalesce(pop_dl, 0), working = age >= 20 & age <= 64)
working_total <- sum(working_rest$rest[working_rest$working])
base_without <- function(million) {
  w <- working_rest |> mutate(cut = if_else(working, million * 1e6 * rest / working_total, 0))
  stopifnot(all(w$cut <= w$rest))
  w |> mutate(pop = pop - cut) |> select(all_of(names(mi$pop22_ini)))
}

# 2. THE COUNTERFACTUAL ON EXPOSURES SMALLER BY THE SAME SHARE =================
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

# the overcount as a share of the exposure: none at the 2001 census, the share
# the cut is of the 2022 base by 2022, linear between
forecast_with <- function(million) {
  share_2022 <- million * 1e6 / working_total
  dt |>
    filter(year %in% 2000:2019) |>
    mutate(over = if_else(age >= 20 & age <= 64, share_2022 * pmax(0, year - 2001) / (2022 - 2001), 0),
           pop = pop * (1 - over),
           mx = deaths / pop) |>
    select(-over) |>
    as_vital(index = year, key = c(sex, age), .age = "age", .sex = "sex",
             .deaths = "deaths", .population = "pop") |>
    model(lc = LC(log(mx), adjust = "e0", jump_choice = "fit")) |>
    forecast(h = 6) |>
    as_tibble() |>
    select(year, sex, age, mx = .mean) |>
    filter(year %in% 2022:2025)
}
with_mx <- function(fc) mi$static_inputs |> select(-mx) |> left_join(fc, by = c("year", "sex", "age"))

scenarios <- expand_grid(million = c(0, 0.5, 1.0, 2.0),
                         counterfactual = c("as fitted", "refitted on the smaller exposures")) |>
  filter(!(million == 0 & counterfactual != "as fitted"))
loss <-
  scenarios |>
  mutate(res = map2(million, counterfactual, function(x, cf) {
    st <- if (cf == "as fitted") mi$static_inputs else with_mx(forecast_with(x))
    stopifnot(all(st$mx > 0), !any(is.na(st$mx)))
    loss_at_mode(mi$draws, st, base_without(x))
  })) |>
  unnest(res)

# with nothing cut the loss is 13's at the mode; with the counterfactual as
# fitted, 13e's
stopifnot(isTRUE(all.equal(
  loss |> filter(million == 0) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))
stopifnot(isTRUE(all.equal(
  loss |> filter(million == 1, counterfactual == "as fitted") |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_denominator_donetsk_luhansk.rds")$migrants |>
    filter(str_detect(scenario, "1.0 million")) |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(loss = loss, working_total = working_total), "data_inter/ukr_base_coherent_e0.rds")

cat("\n=== LOSS AND COUNTERFACTUAL e0, MEN, BY THE CUT AT AGES 20-64 OUTSIDE DONETSK AND LUHANSK ===\n")
print(as.data.frame(
  loss |> filter(sex == "m", year %in% c(2022, 2025)) |>
    transmute(million, counterfactual, year, loss = round(loss, 2), e0_bsn = round(e0_bsn, 2)) |>
    pivot_wider(names_from = year, values_from = c(loss, e0_bsn))
))
message("Done. data_inter/ukr_base_coherent_e0.rds written.")
