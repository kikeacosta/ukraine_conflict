# ==============================================================================
# STEP 06 - Net migration by age and sex, 2022-2025
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Builds the net outflow that the projection subtracts from exposure, and the
# PERT bounds that 10 turns into the migration rows of the parameter table.
# Every number comes from a file in data_input/; the reasoning behind each
# step is in documents/migration_methodology.md. The functions it calls are in
# 00_setup.R (mig_*), shared with 13r.
#
#   1. REGISTER STOCK of refugees from Ukraine recorded outside Russia, at the
#      end of each year: UNHCR's published totals of refugees from Ukraine
#      recorded globally (Factsheet #3, May 2026), which count every country
#      but Russia, Belarus included. They are used as published: no netting
#      of the pre-war 2021 stock, no subtraction of Belarus, no correction of
#      the clean-ups of invalid statuses. UNHCR counts people who have fled
#      Ukraine since 24 February 2022, hold international protection,
#      temporary protection or another form of stay, and are present in the
#      host country.
#
#   2. AGE-SEX NET FLOW BY COHORT DIFFERENCING. Eurostat's temporary
#      protection stock for EU + EFTA at each year end, as published, gives
#      the age-sex structure; it is scaled to the register stock, ungrouped to
#      single ages separately by year and sex, and differenced along the
#      cohort diagonal:
#        D_t(x) = S_t(x) - S_{t-1}(x-1) + deaths abroad
#      Age 0 is set to zero: it is mostly children born abroad, who were never
#      in Ukraine's population. D is negative where a cohort shrank abroad,
#      i.e. net return. The flows are smoothed lightly across age.
#
#   3. TWO READINGS OF THE WESTERN FLOW, which bracket it:
#        register  = sum of D_t (step 2)
#        crossings = CES (2026) table 1 net crossings of citizens, plus the
#                    277k who left through Russia or Belarus and are now
#                    elsewhere (they are in the registers but crossed no
#                    western border)
#      The two readings bracket the flow as whole paths: 11 draws one weight
#      per simulation and blends the two readings' paths with it, so the
#      four-year total runs between the two sources' own totals. The mode is
#      their midpoint in every year. The per-year min and max written below
#      are the range of each year's two readings, kept for display and checks.
#
#   4. RUSSIA, one 2022 entry: mode = Data Finder stock of refugees and asylum
#      seekers at end 2022 net of 2021; min and max set below with their
#      reasons.
#
#   5. ALLOCATION AT THE MODE. The mode total differs from the register
#      reading, and the difference is spread over the DEPARTURES profile
#      (positive part of D), leaving the return cells as the register has
#      them:  ems(x) = D(x) + (mode - sum D) * w(x).  11 applies the same
#      additive rule to every draw.
#
# 13r builds the register stock bottom-up instead (Eurostat with the
# clean-ups back-cast, the Data Finder for other countries, Canada and the USA
# from programme counts), with Russia and Belarus as the separate component,
# and reports how far the loss moves.
#
# INPUTS   data_input/migration/unhcr/unhcr_factsheet3_refugees_recorded_globally.csv
#          data_input/refugees_eurostat/migr_asytpsm__custom_21093689_page_spreadsheet.xlsx
#            (stock at 31 December by country, sex and age group)
#          data_input/migration/unhcr/unhcr_population_coo_UKR_2021_2025.csv
#            (Russia only)
#          data_input/migration/ces_2026_figures.csv
#          data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds     (from 04)
# OUTPUTS  data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds
#            year, sex, age, mix (negative = leaving, at the mode), w
#          data_inter/ukr_migration_bounds.rds       <- used by 10
#          data_inter/ukr_migration_readings.rds     <- used by 10, 13r, 15
#            and the docs (register_stock is UNHCR's total, the stock the
#            register reading differences; ru is the Russian stock net of
#            2021)
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

yrs <- mig_years

# ==============================================================================
# 1. REGISTER STOCK: UNHCR'S TOTALS, 31 DECEMBER
# ==============================================================================
# The chart prints the totals rounded to 0.1 million, at the end of each June
# and December; the December values are the year-end stocks.
register_stock <-
  read_csv("data_input/migration/unhcr/unhcr_factsheet3_refugees_recorded_globally.csv",
           show_col_types = FALSE) |>
  mutate(date = as.Date(date)) |>
  filter(month(date) == 12, year(date) %in% yrs) |>
  transmute(year = as.integer(year(date)), register_stock = refugees_recorded_globally) |>
  arrange(year)

stopifnot(identical(register_stock$year, yrs),
          isTRUE(all.equal(register_stock$register_stock, c(4.7e6, 5.2e6, 5.7e6, 5.9e6))))

print(as.data.frame(register_stock))

# ==============================================================================
# 2. AGE-SEX NET FLOW BY COHORT DIFFERENCING
# ==============================================================================
# The totals carry no age-sex detail, so the register stock takes the
# structure of Eurostat's year-end extract for EU + EFTA, as published. The
# deaths abroad are computed with the forecast death rates.
eurostat_bands <- mig_eurostat_structure(mig_eurostat_extract(yrs))

flows <- mig_register_flows(eurostat_bands, register_stock, mig_mx_exp(yrs), yrs)

stopifnot(!any(is.na(flows$D)))

# The flows before smoothing are kept as D_raw (00_setup.R).
flows <- mig_smoothed(flows)

# ==============================================================================
# 3-4. READINGS AND BOUNDS
# ==============================================================================
crossings <- mig_crossings()

# Russia: the Data Finder's refugees and asylum seekers in the Russian
# Federation, net of the country's pre-war 2021 stock; a year with no record
# carries the year before forward
unhcr_ru <-
  read_csv("data_input/migration/unhcr/unhcr_population_coo_UKR_2021_2025.csv",
           show_col_types = FALSE) |>
  filter(coa_iso == "RUS") |>
  transmute(year = as.integer(year),
            ref_asy = coalesce(refugees, 0) + coalesce(asylum_seekers, 0)) |>
  complete(year = 2021:2025, fill = list(ref_asy = NA)) |>
  arrange(year) |>
  mutate(ref_asy = if_else(year == 2025 & is.na(ref_asy), lag(ref_asy), ref_asy),
         ru = base::pmax(coalesce(ref_asy, 0) - ref_asy[year == 2021], 0)) |>
  filter(year %in% yrs) |>
  select(year, ru)

readings <-
  mig_flow_totals(flows) |>
  left_join(crossings, by = "year") |>
  left_join(register_stock, by = "year") |>
  left_join(unhcr_ru, by = "year")

bounds_west <-
  readings |>
  transmute(
    year,
    component = "west",
    min = base::pmin(register, crossings),
    mode = (register + crossings) / 2,
    max = base::pmax(register, crossings)
  )

# Russia. The mode is UNHCR's statistical stock in the Russian Federation at
# end 2022 net of 2021. That all of it entered in 2022, and that it stayed,
# are assumed, not measured. The floor takes out the onward movement that is
# measured, the 277k CES finds left through Russia or Belarus and are now
# elsewhere, assuming all of them were counted in Russia before they moved on;
# overstatement by the Russian administrative data beyond that is not
# quantified and not allowed for, and neither are Ukrainian citizens who lived
# in Russia before February 2022, whom the mode nets only as the refugees and
# asylum seekers of 2021. The ceiling adds the groups the mode's basis
# (rental contracts, by CES's account) leaves out: about 303k granted Russian
# citizenship and 98k granted temporary asylum in 2022 bring it to about
# 1.64M, and the rest is a margin for informal renters and people in temporary
# accommodation, up to the 2.0M that Russian officials cited in autumn 2022.
# The sources are in documents/migration_methodology.md, section 6.
ru_onward <- crossings |> filter(year == min(yrs)) |> pull(via_ru_by)
ru_max <- 2000000

bounds_ru <-
  readings |>
  filter(year == min(yrs)) |>
  transmute(year, component = "ru", min = ru - ru_onward, mode = ru, max = ru_max)

stopifnot(ru_onward == 277000, bounds_ru$mode == 1266878)

migration_bounds <- bind_rows(bounds_west, bounds_ru)

stopifnot(
  all(migration_bounds$min <= migration_bounds$mode),
  all(migration_bounds$mode <= migration_bounds$max)
)

# ==============================================================================
# 5. ALLOCATION AT THE MODE
# ==============================================================================
mode_tot <- migration_bounds |> summarise(mode_tot = sum(mode), .by = year)

dt5 <- mig_allocate(flows, mode_tot)

write_rds(dt5, "data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds")
write_rds(migration_bounds, "data_inter/ukr_migration_bounds.rds")
write_rds(readings, "data_inter/ukr_migration_readings.rds")

print(as.data.frame(readings |> select(year, register_stock, register, crossings, departures,
                                       returns, deaths_abroad, births_abroad_dropped, ru)))
print(as.data.frame(migration_bounds))

# ==============================================================================
# DIAGNOSTIC FIGURES (exploratory)
# ==============================================================================
nice <- function(d) d |> mutate(sex = if_else(sex == "f", "Females", "Males"))

# net flow by age and sex at the mode, returns below zero
dt5 |>
  nice() |>
  ggplot(aes(age, -mix / 1e3, colour = factor(year))) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line() +
  facet_wrap(~sex) +
  scale_colour_manual(values = c("black", "grey20", "grey45", "grey70")) +
  labs(x = "Age at end of year", y = "Net outflow (thousands, negative = return)",
       colour = NULL) +
  theme_bw()
ggsave("figures/exploratory/migration_net_flow_by_age_mode.png", w = 8, h = 3.5)

# the same, from 2023, where the scale lets the returns show
dt5 |>
  filter(year > 2022) |>
  nice() |>
  ggplot(aes(age, -mix / 1e3, colour = factor(year))) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line() +
  facet_wrap(~sex) +
  labs(x = "Age at end of year", y = "Net outflow (thousands)", colour = NULL) +
  theme_bw()
ggsave("figures/exploratory/migration_net_flow_by_age_mode_2023_2025.png", w = 8, h = 3.5)
