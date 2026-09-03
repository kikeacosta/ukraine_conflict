# ==============================================================================
# STEP 04 - Counterfactual mortality: Lee-Carter forecast for 2022-2025
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Produces the "expected" mortality schedule: what death rates would have been
# had the war not happened. This is the baseline everything else is measured
# against.
#
# A Lee-Carter model (Lee-Miller variant: adjusted on e0, jumping off the
# fitted rates) is fitted to 2000-2019 and forecast six years ahead.
#
# 2020 and 2021 are deliberately EXCLUDED from the fitting window. They are
# COVID years, and including them would build pandemic excess mortality into
# a counterfactual that is supposed to represent normal conditions.
#
# INPUTS   data_inter/ukr_life_tables_1989_2021.rds   (from 03)
#          data_input/DataDxEx.csv                    (exposures)
# OUTPUT   data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds   <- used by 11, 13
#          Contains BOTH observed and forecast rows, tagged by `source`;
#          2020 and 2021 appear under both. Downstream scripts filter on
#          source == "frcst".
# ==============================================================================

rm(list = ls(all = TRUE))

## loading LC functions
source("code/00_setup.R")

# loading data
mx_all <- read_rds("data_inter/ukr_life_tables_1989_2021.rds")

# Load population data
pop_cnt <- read.csv("data_input/DataDxEx.csv", header = TRUE)

pop_cnt2 <-
  pop_cnt %>%
  rename_with(tolower) %>%
  filter(data == "Ex") %>%
  pivot_longer(starts_with("age"), names_to = "age", values_to = "pop") %>%
  mutate(
    age = as.numeric(gsub("age", "", age)),
    sex = str_sub(sex, 1, 1)
  ) %>%
  select(year, sex, age, pop)

dt <-
  mx_all %>%
  select(year, sex, age, mx) %>%
  mutate(sex = str_sub(sex, 1, 1)) %>%
  left_join(pop_cnt2, by = c("year", "sex", "age")) %>%
  mutate(deaths = mx * pop) %>%
  drop_na(mx, pop)

# Create vital object
ukr_vital <-
  dt %>%
  filter(year %in% 2000:2019) %>%
  as_vital(
    index = year,
    key = c(sex, age),
    .age = "age",
    .sex = "sex",
    .deaths = "deaths",
    .population = "pop"
  )

# Forecast using Lee-Carter model (Lee-Miller variant based on life expectancy) for 2022-2025
fc_rates <-
  ukr_vital %>%
  model(lc = LC(log(mx), adjust = "e0", jump_choice = "fit")) %>%
  forecast(h = 6)

# Extract forecasted rates
forecast_df <-
  fc_rates %>%
  as_tibble() %>%
  select(year, sex, age, mx = .mean) |>
  mutate(source = "frcst")


tst <-
  bind_rows(
    dt |>
      mutate(source = "obs"),
    forecast_df
  )

tst_out <- tst %>%
  select(-pop, -deaths)

# saving forecast estimates
write_rds(tst_out, "data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# visualizing observed and projected life expectancy ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tst |>
  mutate(
    sex = case_when(sex == "m" ~ "Male", sex == "f" ~ "Female"),
    sex = factor(sex, levels = c("Male", "Female")),
    source = case_when(
      source == "obs" ~ "Observed",
      source == "frcst" ~ "Forecast"
    ),
    source = factor(source, levels = c("Observed", "Forecast"))
  ) |>
  filter(year %in% 2019:2025) |>
  ggplot(aes(x = age, y = mx, lty = source)) +
  geom_line() +
  scale_y_log10() +
  facet_grid(sex ~ year) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave(
  "figures/exploratory/labtalk/ukr_mx_forecasted.png",
  w = 8,
  h = 4
)

# Calculate life expectancy for forecasted years
e0_forecast <-
  fc_rates %>%
  as_tibble() %>%
  select(year, sex, Age = age, mx = .mean) %>%
  as_vital(
    index = year,
    key = c(sex, Age),
    .age = "Age",
    .sex = "sex"
  ) %>%
  life_expectancy() |>
  filter(Age == 0) |>
  mutate(source = "frcst")

e0_obs <-
  dt %>%
  filter(year %in% 1995:2019) %>%
  select(year, sex, Age = age, mx) %>%
  as_vital(
    index = year,
    key = c(sex, Age),
    .age = "Age",
    .sex = "sex"
  ) %>%
  life_expectancy() |>
  filter(Age == 0) |>
  mutate(source = "obs")

e0 <- bind_rows(e0_obs, e0_forecast)

e0 |>
  mutate(
    sex = case_when(sex == "m" ~ "Male", sex == "f" ~ "Female"),
    sex = factor(sex, levels = c("Male", "Female")),
    source = case_when(
      source == "obs" ~ "Observed",
      source == "frcst" ~ "Forecast"
    ),
    source = factor(source, levels = c("Observed", "Forecast"))
  ) |>
  ggplot(aes(x = year, y = ex, lty = source)) +
  facet_wrap(~sex) +
  geom_line() +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave(
  "figures/exploratory/labtalk/ukr_le_forecasted.png",
  w = 8,
  h = 3
)
