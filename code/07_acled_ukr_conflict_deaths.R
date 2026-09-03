# ==============================================================================
# STEP 07 (ACLED) - Alternative conflict death totals, for comparison only
# ==============================================================================
#
# ACLED fatality counts are read here as an independent cross-check on the
# UCDP totals used by the pipeline. NOTHING DOWNSTREAM CONSUMES THIS OUTPUT -
# the estimates all run on UCDP (07_ucdp_ukr_conflict_deaths.R).
#
# The two ACLED workbooks are not tracked in git. As with the other external
# sources, the tidy country-year-role extract below is cached, so the raw
# files only ever have to be read once and a fresh clone runs without them.
# See cache_rds() in 00_setup.R.
#
# INPUTS   data_input/acled/number_of_reported_fatalities_*.xlsx
#          data_input/acled/number_of_reported_civilian_fatalities_*.xlsx
#          -> cached as data_inter/acled_fatalities_slim.rds
# OUTPUTS  data_inter/ukr_acled.rds                    (comparison only)
#          data_inter/ukr_rus_acled_source_comparison.rds  (used by 15)
# ==============================================================================

rm(list = ls())
source("code/00_setup.R")

acled_dir <- "data_input/acled"

dt <- cache_rds("data_inter/acled_fatalities_slim.rds", {
  dts <- read_xlsx(require_raw(file.path(
    acled_dir,
    "number_of_reported_fatalities_by_country-year_as-of-22Aug2025.xlsx"
  )))
  cvs <- read_xlsx(require_raw(file.path(
    acled_dir,
    "number_of_reported_civilian_fatalities_by_country-year_as-of-22Aug2025.xlsx"
  )))

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
})

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
# in 15, alongside the UCDP figures.
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
