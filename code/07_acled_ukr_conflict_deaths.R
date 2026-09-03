# ==============================================================================
# STEP 07 (ACLED) - Alternative conflict death totals, for comparison only
# ==============================================================================
#
# ACLED fatality counts are read here as an independent cross-check on the
# UCDP totals used by the pipeline. NOTHING DOWNSTREAM CONSUMES THIS OUTPUT -
# the estimates all run on UCDP (07_ucdp_ukr_conflict_deaths_invals.R).
#
# INPUTS   data_input/acled/number_of_reported_fatalities_*.xlsx
#          data_input/acled/number_of_reported_civilian_fatalities_*.xlsx
# OUTPUT   data_inter/ukr_acled.rds   (comparison only)
# ==============================================================================

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
  reframe(dts = sum(dts), .by = c(role))

# Ukraine and Russia totals for the source reconciliation table assembled
# in 16, alongside the UCDP figures.
write_rds(
  dt |>
    filter(country %in% c("Ukraine", "Russia"), year %in% 2022:2025) |>
    summarise(dts = sum(dts), .by = c(country, role)) |>
    pivot_wider(names_from = role, values_from = dts) |>
    mutate(unknown = NA_real_, source = "ACLED", type = "unadjusted"),
  "data_inter/ukr_rus_acled_source_comparison.rds"
)


dt %>%
  filter(country == "Colombia") |>
  spread(role, dts)
