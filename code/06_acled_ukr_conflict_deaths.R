rm(list = ls())
source("code/00_setup.R")

dts <-
  read_xlsx(
    "data_input/acled/number_of_reported_fatalities_by_country-year_as-of-22Aug2025.xlsx"
  )
cvs <-
  read_xlsx(
    "data_input/acled/number_of_reported_civilian_fatalities_by_country-year_as-of-22Aug2025.xlsx"
  )

dt <-
  dts %>%
  rename_with(tolower) %>%
  rename(dts = fatalities) %>%
  left_join(
    cvs %>%
      rename_with(tolower) %>%
      rename(civilians = fatalities)
  ) %>%
  replace_na(list(dts = 0, civilians = 0)) %>%
  mutate(combatants = dts - civilians) %>%
  select(-dts) %>%
  gather(civilians, combatants, key = role, value = dts)

dt %>%
  filter(year %in% 2014:2025) %>%
  group_by(country, year) %>%
  summarise(dts = sum(dts))

dt %>%
  filter(year %in% 2014:2025) %>%
  group_by(country, role, year) %>%
  summarise(dts = sum(dts))

ukr <-
  dt %>%
  filter(country == "Ukraine", year %in% 2022:2025) %>%
  select(-country)

ukr %>%
  filter(year >= 2022) %>%
  reframe(dts = sum(dts), .by = role)

write_rds(ukr, "data_inter/ukr_acled.rds")

# looking at conflict mortality in Russia
rus <-
  dt %>%
  filter(country == "Russia", year %in% 2022:2025) %>%
  select(-country)

rus %>%
  filter(year >= 2022) %>%
  reframe(dts = sum(dts), .by = role)
