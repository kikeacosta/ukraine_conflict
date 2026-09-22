# ==============================================================================
# STEP 04b - The territory of the pre-war series, and a level shift in k at 2014-15
# ==============================================================================
#
# WHY
# ---
# The counterfactual is a Lee-Carter forecast fitted to 2000-2019 (04). From 2014
# the occupation of Crimea and of parts of Donetsk and Luhansk changed what the
# State Statistics Service could register. If the series behind the fit changed
# territory then, the index k could step at 2014-15 for a reason that has
# nothing to do with mortality: the step would enter the drift, and it would
# inflate the variance of the yearly changes, which is what sets the forecast
# band that dominates the male interval from 2023.
#
# WHAT IT DOES
# ------------
#   1. the deaths and exposures behind the fit, all ages, by year: whether the
#      population they cover steps at 2014 or 2015
#   2. the fitted index of each sex (the same model as 04) and its yearly
#      changes; a level shift at year s is a change at s out of line with the
#      others, tested by a dummy for s in a regression of the changes on a
#      constant (the drift), for 2014, for 2015 and for the two together
#   3. the drift and the standard deviation of the changes with those years
#      left out, and the standard deviation of a forecast six years ahead
#      (2025) under each: what a step would do to the band if it were one
#
# Not part of the estimates: it reports, and changes nothing downstream.
#
# INPUTS   data_input/DataDxEx.csv, data_inter/ukr_life_tables_1989_2021.rds (03)
# OUTPUTS  data_inter/ukr_lc_level_shift.rds
# ==============================================================================

rm(list = ls(all = TRUE))
source("code/00_setup.R")

# 1. THE TERRITORY OF THE SERIES ==============================================
dxex <-
  read.csv("data_input/DataDxEx.csv", header = TRUE) |>
  rename_with(tolower) |>
  pivot_longer(starts_with("age"), names_to = "age", values_to = "value") |>
  summarise(value = sum(value, na.rm = TRUE), .by = c(year, region, data))
coverage <-
  dxex |>
  pivot_wider(names_from = data, values_from = value) |>
  arrange(year) |>
  mutate(exposure_change = Ex / lag(Ex) - 1)
stopifnot(n_distinct(coverage$region) == 1)

# 2. THE INDEX AND ITS YEARLY CHANGES, AS 04 FITS THEM ========================
mx_all <- read_rds("data_inter/ukr_life_tables_1989_2021.rds")
pop <-
  read.csv("data_input/DataDxEx.csv", header = TRUE) |>
  rename_with(tolower) |>
  filter(data == "Ex") |>
  pivot_longer(starts_with("age"), names_to = "age", values_to = "pop") |>
  mutate(age = as.numeric(gsub("age", "", age)), sex = str_sub(sex, 1, 1)) |>
  select(year, sex, age, pop)
ukr_vital <-
  mx_all |>
  select(year, sex, age, mx) |>
  mutate(sex = str_sub(sex, 1, 1)) |>
  left_join(pop, by = c("year", "sex", "age")) |>
  mutate(deaths = mx * pop) |>
  drop_na(mx, pop) |>
  filter(year %in% 2000:2019) |>
  as_vital(index = year, key = c(sex, age), .age = "age", .sex = "sex",
           .deaths = "deaths", .population = "pop")
kt <-
  ukr_vital |>
  model(lc = LC(log(mx), adjust = "e0", jump_choice = "fit")) |>
  time_components() |>
  as_tibble() |>
  arrange(sex, year) |>
  mutate(dk = kt - lag(kt), .by = sex) |>
  drop_na(dk)
# the fit must be 04's
walk04 <- read_rds("data_inter/ukr_lc_forecast_error.rds")$walk
stopifnot(isTRUE(all.equal(
  kt |> summarise(drift = mean(dk), sigma = sd(dk), .by = sex) |> arrange(sex) |> select(drift, sigma),
  walk04 |> arrange(sex) |> select(drift, sigma),
  check.attributes = FALSE, tolerance = 1e-8)))

# 3. THE TEST, AND THE BAND WITH AND WITHOUT THE STEP YEARS ===================
shifts <- list("2014" = 2014, "2015" = 2015, "2014 and 2015" = c(2014, 2015))
test <- map_dfr(names(shifts), function(nm) {
  kt |>
    mutate(step = as.numeric(year %in% shifts[[nm]])) |>
    group_by(sex) |>
    group_modify(function(d, key) {
      fit <- lm(dk ~ step, data = d)
      co <- summary(fit)$coefficients["step", ]
      rest <- d$dk[d$step == 0]
      n <- length(rest) + 1                    # fitted years behind the changes kept
      tibble(shift_at = nm, step = co[["Estimate"]], t = co[["t value"]], p = co[["Pr(>|t|)"]],
             drift_all = mean(d$dk), sigma_all = sd(d$dk),
             drift_without = mean(rest), sigma_without = sd(rest),
             sd_2025_all = sd(d$dk) * sqrt(6 * (1 + 6 / (nrow(d)))),
             sd_2025_without = sd(rest) * sqrt(6 * (1 + 6 / (n - 1))))
    }) |>
    ungroup()
})

cat("\n=== DEATHS AND EXPOSURES BEHIND THE FIT, ALL AGES ===\n")
print(as.data.frame(coverage |> filter(year >= 2010) |>
                      mutate(across(c(Dx, Ex), round), exposure_change = round(100 * exposure_change, 2))))
cat("\n=== A LEVEL SHIFT IN THE LEE-CARTER INDEX ===\n")
print(as.data.frame(test |> mutate(across(where(is.numeric), \(v) signif(v, 3)))))

write_rds(list(coverage = coverage, kt = kt, test = test), "data_inter/ukr_lc_level_shift.rds")
message("Done. data_inter/ukr_lc_level_shift.rds written.")
