# ==============================================================================
# STEP 06 - Net emigration by age and sex, 2022-2025
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Builds the net outflow that the projection subtracts from exposure, and the
# PERT bounds that 10 turns into the migration rows of the parameter table.
# Every number comes from a file in data_input/; the reasoning behind each
# step is in documents/migration_methodology.md.
#
#   1. REGISTER STOCK of people displaced outside Russia and Belarus, at the
#      end of each year, built bottom-up so that each part is counted in its
#      own source's definitions:
#        Eurostat temporary protection, EU + EFTA (the native source, and the
#          only one with age and sex)
#      + UNHCR Refugee Data Finder, refugees + asylum seekers, for every other
#          country of asylum except Russia, Belarus, Canada and the USA, net
#          of each country's pre-war 2021 stock
#      + Canada (CUAET) and USA (Uniting for Ukraine) from national arrival
#          figures. The Data Finder cannot be used for either: it counts only
#          ~20k in the USA (parolees are not refugees) and holds Canada flat at
#          75k before a +98k reclassification step in 2025.
#          >>> PLACEHOLDER: approximated by linear interpolation between
#          >>> published anchor points until year-end series are obtained.
#
#   2. AGE-SEX NET FLOW BY COHORT DIFFERENCING. Eurostat's age-sex structure
#      is scaled to the register stock, ungrouped to single ages separately by
#      year and sex, and differenced along the cohort diagonal:
#        D_t(x) = S_t(x) - S_{t-1}(x-1) + deaths abroad
#      Age 0 is set to zero: it is mostly children born abroad, who were never
#      in Ukraine's population. D is negative where a cohort shrank abroad,
#      i.e. net return.
#
#   3. TWO READINGS OF THE WESTERN FLOW, which bracket it:
#        register  = sum of D_t (step 2)
#        crossings = CES (2026) table 1 net crossings of citizens, plus the
#                    277k who left through Russia or Belarus and are now
#                    elsewhere (they are in the registers but crossed no
#                    western border)
#      Bounds are the two readings; the mode is their midpoint.
#
#   4. RUSSIA AND BELARUS, one 2022 entry: mode = Data Finder stock at end
#      2022 net of 2021; min and max set below with their reasons.
#
#   5. ALLOCATION AT THE MODE. The mode total differs from the register
#      reading, and the difference is spread over the DEPARTURES profile
#      (positive part of D), leaving the return cells as the register has
#      them:  ems(x) = D(x) + (mode - sum D) * w(x).  11 applies the same
#      additive rule to every draw.
#
# INPUTS   data_input/refugees_eurostat/migr_asytpsm__custom_21093689_page_spreadsheet.xlsx
#          data_input/migration/unhcr/unhcr_population_coo_UKR_2021_2025.csv
#          data_input/migration/national_programme_arrivals_PLACEHOLDER.csv
#          data_input/migration/ces_2026_figures.csv
#          data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds     (from 04)
# OUTPUTS  data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds
#            year, sex, age, mix (negative = leaving, at the mode), w
#          data_inter/ukr_migration_bounds.rds       <- used by 10
#          data_inter/ukr_migration_readings.rds     <- used by 15 and the docs
# ==============================================================================

# ==============================================================================
# !!! TODO-PLACEHOLDER: CANADA AND USA ARRIVALS (search for this tag) !!!
# ==============================================================================
# The register stock includes Ukrainians who went to Canada and the USA. The
# UNHCR Data Finder does not count them properly (see section 1c), so they
# come from national programme figures. Year-end series were not available
# when this was written, so section 1c APPROXIMATES them by linear
# interpolation between a few published cumulative totals stored in
#   data_input/migration/national_programme_arrivals_PLACEHOLDER.csv
#
# DATA NEEDED TO REPLACE IT
#
#   CANADA - Canada-Ukraine Authorization for Emergency Travel (CUAET)
#     What:  CUMULATIVE number of CUAET holders who ARRIVED in Canada since
#            17 March 2022 (entries, not applications or approvals), at
#            31 Dec 2022, 31 Dec 2023, 31 Dec 2024 and 31 Dec 2025, or at the
#            nearest dates available. Monthly would be better still.
#            Useful if it exists: number of CUAET holders still residing in
#            Canada at each year end (a stock rather than arrivals).
#     Known: 298,128 arrivals 17 Mar 2022 - 1 Apr 2024 (entry deadline was
#            31 Mar 2024, so later years change little).
#     Where: (1) IRCC page "Canada-Ukraine authorization for emergency travel:
#                Key figures"
#                https://www.canada.ca/en/immigration-refugees-citizenship/services/immigrate-canada/ukraine-measures/key-figures.html
#                It was updated about weekly while CUAET ran. Open its history
#                in the Wayback Machine, https://web.archive.org/ , and read
#                the arrivals figure on the snapshots nearest each 31 Dec.
#            (2) IRCC open data, https://open.canada.ca , search "CUAET".
#            (3) If neither works: an access-to-information request to IRCC
#                for monthly CUAET arrivals, March 2022 - December 2025.
#
#   USA - Uniting for Ukraine (U4U) parole
#     What:  CUMULATIVE number of Ukrainians PAROLED INTO THE USA under U4U
#            (arrivals, not supporter applications or travel authorisations)
#            at 31 Dec 2022, 2023, 2024 and 2025. Ideally ALSO Ukrainians
#            paroled at ports of entry between 24 Feb and 25 Apr 2022, before
#            U4U existed, which the current anchors leave out.
#     Known: 93,928 U4U arrivals by 13 Dec 2022 (USCIS letter, FOIA library);
#            about 260,000 Ukrainians in the USA in Nov 2025 (CES 2026 citing
#            Reuters; a count present, not arrivals).
#     Where: (1) DHS Office of Homeland Security Statistics, https://ohss.dhs.gov ,
#                parole and immigration statistics tables, if they break out U4U.
#            (2) USCIS Uniting for Ukraine pages and the USCIS FOIA reading
#                room (congressional correspondence often quotes totals).
#            (3) If neither works: a FOIA request to USCIS/CBP for monthly U4U
#                parole arrivals, April 2022 - December 2025.
#
# HOW TO ADD IT (no code changes needed)
#   1. Append one row per country and year end to the CSV above:
#        country (CAN or USA), date (YYYY-12-31), cumulative_arrivals,
#        qualifier ("exact"), source, url, accessed, note.
#      Delete any anchor the new rows make redundant.
#   2. Save a copy of each source page or file (PDF/HTML/XLSX) in
#        data_input/migration/web_snapshots/
#   3. Rename the CSV (drop "_PLACEHOLDER") and update the path in section 1c.
#   4. Re-run 06, 10, 11, 12, 13, 14, 15, and update
#        documents/migration_methodology.md (section "Placeholder: Canada
#        and USA") and documents/data_sources.md.
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

yrs <- 2022:2025

# ==============================================================================
# 1a. EUROSTAT TEMPORARY PROTECTION STOCK, EU + EFTA, 31 DECEMBER
# ==============================================================================
# The extract holds December of each year for EU-27 countries plus Iceland,
# Liechtenstein, Norway and Switzerland. The EU-27 aggregate row is dropped so
# nothing is counted twice.
dt <-
  read_xlsx(
    "data_input/refugees_eurostat/migr_asytpsm__custom_21093689_page_spreadsheet.xlsx",
    skip = 10
  )

eu_row <- "European Union - 27 countries (from 2020)"
cts_rmv <- c("Special value", ":", "Observation flags:", "dp", "d")

cols_names <-
  expand_grid(date = yrs, sex = c("t", "m", "f", "u"), dat = c("data", "flag")) |>
  mutate(nms = paste(date, sex, dat, sep = "_")) |>
  pull(nms)

dt2 <-
  dt |>
  set_names(c("ctr", "age", cols_names)) |>
  filter(!ctr %in% cts_rmv) |>
  drop_na(ctr) |>
  select(-ends_with("flag")) |>
  pivot_longer(-c(ctr, age), names_to = "vars", values_to = "mix") |>
  mutate(
    year = as.integer(str_sub(vars, 1, 4)),
    sex = str_sub(vars, 6, 6),
    mix = suppressWarnings(as.double(mix))
  ) |>
  replace_na(list(mix = 0)) |>
  select(ctr, year, sex, age, mix)

eurostat_tot <-
  dt2 |>
  filter(ctr != eu_row, sex == "t", age == "Total") |>
  summarise(eurostat = sum(mix), .by = year)

# unknown sex and age are allocated with the distribution of the known cases
eurostat_bands <-
  dt2 |>
  filter(ctr != eu_row, sex %in% c("m", "f"), !age %in% c("Total", "Unknown")) |>
  summarise(mix = sum(mix), .by = c(year, sex, age)) |>
  mutate(cx = mix / sum(mix), .by = year) |>
  mutate(age = case_when(
    age == "Less than 14 years" ~ 0,
    age == "From 14 to 17 years" ~ 14,
    age == "From 18 to 34 years" ~ 18,
    age == "From 35 to 64 years" ~ 35,
    age == "65 years or over" ~ 65
  )) |>
  select(year, sex, age, cx)

# ==============================================================================
# 1b. UNHCR REFUGEE DATA FINDER, OUTSIDE EU + EFTA, RUSSIA, BELARUS, CAN, USA
# ==============================================================================
eu_efta_iso <- c(
  "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU",
  "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT",
  "ROU", "SVK", "SVN", "ESP", "SWE", "ISL", "LIE", "NOR", "CHE"
)

unhcr <-
  read_csv("data_input/migration/unhcr/unhcr_population_coo_UKR_2021_2025.csv",
           show_col_types = FALSE) |>
  transmute(
    year = as.integer(year), coa_iso,
    ref_asy = coalesce(refugees, 0) + coalesce(asylum_seekers, 0)
  )

# each country's 2021 stock predates the invasion and is not war displacement;
# the few countries with no 2025 record carry their 2024 value forward
unhcr_net <-
  unhcr |>
  complete(year = 2021:2025, coa_iso, fill = list(ref_asy = NA)) |>
  arrange(coa_iso, year) |>
  group_by(coa_iso) |>
  mutate(
    ref_asy = if_else(year == 2025 & is.na(ref_asy), lag(ref_asy), ref_asy),
    base_2021 = coalesce(ref_asy[year == 2021], 0),
    net = base::pmax(coalesce(ref_asy, 0) - base_2021, 0)
  ) |>
  ungroup() |>
  filter(year %in% yrs)

unhcr_other <-
  unhcr_net |>
  filter(!coa_iso %in% c(eu_efta_iso, "RUS", "BLR", "CAN", "USA", "UKR")) |>
  summarise(unhcr_other = sum(net), .by = year)

unhcr_ru_by <-
  unhcr_net |>
  filter(coa_iso %in% c("RUS", "BLR")) |>
  summarise(ru_by = sum(net), .by = year)

# ==============================================================================
# 1c. CANADA AND USA - PLACEHOLDER
# ==============================================================================
# !!! TODO-PLACEHOLDER - see the note at the top of this script for exactly
# which data replace this and where to find them.
#
# Why not the Data Finder: it counts ~15-22k Ukrainians in the USA (asylum
# seekers and refugees only; U4U parolees are not refugees, and US arrivals
# were ~260k), and holds Canada flat at ~75k in 2022-2024 before jumping to
# 174k in 2025, a reporting change rather than a flow.
#
# Current approximation: linear interpolation between dated cumulative-arrival
# anchors, evaluated at 31 December, flat after the last anchor. Arrivals are
# used as if they were a stock, which overstates it by whoever has since left
# Canada or the USA.
anchors <-
  read_csv("data_input/migration/national_programme_arrivals_PLACEHOLDER.csv",
           show_col_types = FALSE) |>
  mutate(date = as.Date(date))

programme <-
  anchors |>
  group_by(country) |>
  reframe(
    year = yrs,
    arrivals = approx(
      x = as.numeric(date), y = cumulative_arrivals,
      xout = as.numeric(as.Date(paste0(yrs, "-12-31"))),
      rule = 2
    )$y
  ) |>
  pivot_wider(names_from = country, values_from = arrivals) |>
  rename(canada_placeholder = CAN, usa_placeholder = USA)

register_stock <-
  eurostat_tot |>
  left_join(unhcr_other, by = "year") |>
  left_join(programme, by = "year") |>
  mutate(register_stock = eurostat + unhcr_other + canada_placeholder + usa_placeholder)

print(as.data.frame(register_stock))

# ==============================================================================
# 2. AGE-SEX NET FLOW BY COHORT DIFFERENCING
# ==============================================================================
# The destinations outside EU + EFTA carry no age-sex detail, so the register
# stock takes Eurostat's structure for that year. Ungrouping is done on the
# stocks, which are positive, and only then differenced.
ung_stock <- function(d) {
  y <- base::pmax(d$mix * 1e5, 1)
  fit <- suppressWarnings(pclm(x = d$age, y = y, nlast = 36)$fitted)
  tibble(age = 0:100, S = fit / 1e5)
}

stock_sa <-
  eurostat_bands |>
  left_join(register_stock |> select(year, register_stock), by = "year") |>
  mutate(mix = cx * register_stock) |>
  group_by(year, sex) |>
  group_modify(~ ung_stock(.x)) |>
  ungroup()

mx_exp <-
  read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") |>
  filter(source == "frcst", year %in% yrs) |>
  select(year, sex, age, mx)

# Rows are indexed by age at the END of year t, the cohort born in t - x, which
# is how 11 carries its population. The open interval 100+ collects the
# previous year's 99 and 100+. Deaths abroad are added back, otherwise a
# cohort that merely died abroad would read as having returned.
cohort_diff <- function(d) {
  map_dfr(yrs, function(t) {
    now <- d |> filter(year == t) |> arrange(age) |> pull(S)
    if (t == min(yrs)) {
      return(tibble(year = t, age = 0:100, S_prev = 0, S = now))
    }
    b <- d |> filter(year == t - 1) |> arrange(age) |> pull(S)
    tibble(year = t, age = 0:100, S_prev = c(NA, b[1:99], b[100] + b[101]), S = now)
  })
}

flows <-
  stock_sa |>
  group_by(sex) |>
  group_modify(~ cohort_diff(.x)) |>
  ungroup() |>
  left_join(mx_exp, by = c("year", "sex", "age")) |>
  mutate(
    deaths_abroad = if_else(year == min(yrs), 0, coalesce(S_prev * mx, 0)),
    D = if_else(age == 0, 0, S - coalesce(S_prev, 0) + deaths_abroad)
  )

stopifnot(!any(is.na(flows$D)))

# ==============================================================================
# 3-4. READINGS AND BOUNDS
# ==============================================================================
ces <- read_csv("data_input/migration/ces_2026_figures.csv", show_col_types = FALSE)

crossings <-
  ces |>
  filter(key == "net_outflow_citizens") |>
  select(year, ces_net = value) |>
  left_join(
    ces |>
      filter(key == "left_via_russia_belarus_now_elsewhere") |>
      select(year, via_ru_by = value),
    by = "year"
  ) |>
  mutate(year = as.integer(year), via_ru_by = coalesce(via_ru_by, 0),
         crossings = ces_net + via_ru_by)

readings <-
  flows |>
  summarise(
    register = sum(D),
    departures = sum(base::pmax(D, 0)),
    returns = sum(base::pmin(D, 0)),
    deaths_abroad = sum(deaths_abroad),
    births_abroad_dropped = sum(S[age == 0]),
    .by = year
  ) |>
  left_join(crossings, by = "year") |>
  left_join(register_stock, by = "year") |>
  left_join(unhcr_ru_by, by = "year")

bounds_west <-
  readings |>
  transmute(
    year,
    component = "west",
    min = base::pmin(register, crossings),
    mode = (register + crossings) / 2,
    max = base::pmax(register, crossings)
  )

# Russia and Belarus. The mode is UNHCR's own end-2022 statistical stock, which
# stays level through 2023 (so the channel filled in the first year) and is
# consistent with the 1.3M of April 2025. The floor allows for Russian
# administrative data overstating and for onward movement; the ceiling for
# the rental-contract basis missing non-renters, informal renters and the
# naturalised, and for the 1.5-2.0M Russian officials cited in autumn 2022.
ru_by_min <- 700000
ru_by_max <- 2000000

bounds_ru_by <-
  readings |>
  filter(year == min(yrs)) |>
  transmute(year, component = "ru_by", min = ru_by_min, mode = ru_by, max = ru_by_max)

migration_bounds <- bind_rows(bounds_west, bounds_ru_by)

stopifnot(
  all(migration_bounds$min <= migration_bounds$mode),
  all(migration_bounds$mode <= migration_bounds$max)
)

# ==============================================================================
# 5. ALLOCATION AT THE MODE
# ==============================================================================
mode_tot <- migration_bounds |> summarise(mode_tot = sum(mode), .by = year)

dt5 <-
  flows |>
  mutate(w = base::pmax(D, 0) / sum(base::pmax(D, 0)), .by = year) |>
  left_join(readings |> select(year, register), by = "year") |>
  left_join(mode_tot, by = "year") |>
  mutate(
    ems = D + (mode_tot - register) * w,
    mix = -ems # negative = people leaving, the convention 11 and 13 read
  ) |>
  select(year, sex, age, mix, w) |>
  arrange(year, sex, age)

stopifnot(
  isTRUE(all.equal(
    dt5 |> summarise(v = -sum(mix), .by = year) |> arrange(year) |> pull(v),
    mode_tot |> arrange(year) |> pull(mode_tot),
    tolerance = 1e-8
  )),
  all(abs(dt5 |> summarise(s = sum(w), .by = year) |> pull(s) - 1) < 1e-9)
)

write_rds(dt5, "data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds")
write_rds(migration_bounds, "data_inter/ukr_migration_bounds.rds")
write_rds(readings, "data_inter/ukr_migration_readings.rds")

print(as.data.frame(readings |> select(year, register, crossings, departures, returns,
                                       deaths_abroad, births_abroad_dropped, ru_by)))
print(as.data.frame(migration_bounds))

# ==============================================================================
# DIAGNOSTIC FIGURES (exploratory)
# ==============================================================================
nice <- function(d) d |> mutate(sex = if_else(sex == "f", "Females", "Males"))

# register stock and its parts
register_stock |>
  select(year, eurostat, unhcr_other, canada_placeholder, usa_placeholder) |>
  pivot_longer(-year, names_to = "part", values_to = "stock") |>
  mutate(part = factor(part, levels = c("usa_placeholder", "canada_placeholder",
                                        "unhcr_other", "eurostat"))) |>
  ggplot(aes(factor(year), stock / 1e6, fill = part)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(eurostat = "#0A9396", unhcr_other = "#94D2BD",
                               canada_placeholder = "#EE9B00", usa_placeholder = "#CA6702")) +
  labs(x = NULL, y = "Stock at 31 December (millions)", fill = NULL,
       caption = "Canada and USA are placeholders interpolated between published anchors.") +
  theme_minimal()
ggsave("figures/exploratory/migration_register_stock_components.png", w = 7, h = 4)

# net flow by age and sex at the mode, returns below zero
dt5 |>
  nice() |>
  ggplot(aes(age, -mix / 1e3, colour = factor(year))) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line() +
  facet_wrap(~sex) +
  scale_colour_manual(values = c("black", "grey20", "grey45", "grey70")) +
  labs(x = "Age at end of year", y = "Net emigration (thousands, negative = return)",
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
  labs(x = "Age at end of year", y = "Net emigration (thousands)", colour = NULL) +
  theme_bw()
ggsave("figures/exploratory/migration_net_flow_by_age_mode_2023_2025.png", w = 8, h = 3.5)
