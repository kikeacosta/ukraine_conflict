# ==============================================================================
# STEP 06 - Net migration by age and sex, 2022-2025
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
#          only one with age and sex), with the clean-ups of invalid statuses
#          that Eurostat flags as breaks in the monthly series treated as
#          breaks (section 1a): the Netherlands (January 2023), Bulgaria
#          (February 2024), Germany (July 2024), Lithuania (October 2024) and
#          Italy (June 2025). Each step is measured in the monthly series by
#          sex and age group and spread back, at a constant rate, over the
#          months from the country's first month with a stock to the break,
#          so the year-end stocks before the break are lowered and those
#          after it stand as published
#      + UNHCR Refugee Data Finder, refugees + asylum seekers, for every other
#          country of asylum except Russia, Belarus, Canada and the USA, net
#          of each country's pre-war 2021 stock
#      + Canada (CUAET) and USA (Uniting for Ukraine) from national programme
#          figures, interpolated between published counts (section 1c). The
#          Data Finder cannot be used for either: it counts only ~20k in the
#          USA (parolees are not refugees) and holds Canada flat at 75k before
#          a +98k reclassification step in 2025.
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
#      The two readings bracket the flow as whole paths: 11 draws one weight
#      per simulation and blends the two readings' paths with it, so the
#      four-year total runs between the two sources' own totals. The mode is
#      their midpoint in every year. The per-year min and max written below
#      are the range of each year's two readings, kept for display and checks.
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
#   6. THE CLEAN-UPS READ AS RETURNS. The register reading is also computed
#      from the stock as published, where each clean-up is a fall in its own
#      year, and written beside the reading as used, with the measured steps
#      and the year-end stocks before and after the correction.
#
# INPUTS   data_input/refugees_eurostat/migr_asytpsm__custom_21093689_page_spreadsheet.xlsx
#            (stock at 31 December by country, sex and age group)
#          data_input/refugees_eurostat/migr_asytpsm_api_UA_monthly_by_country_sex_age.csv
#            (the same table, monthly, March 2022 to December 2025, from
#            Eurostat's SDMX 2.1 API: data/migr_asytpsm/M.PER.UA...?startPeriod=
#            2022-03&endPeriod=2025-12&format=SDMX-CSV)
#          data_input/migration/unhcr/unhcr_population_coo_UKR_2021_2025.csv
#          data_input/migration/national_programme_counts.csv
#          data_input/migration/ces_2026_figures.csv
#          data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds     (from 04)
# OUTPUTS  data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds
#            year, sex, age, mix (negative = leaving, at the mode), w
#          data_inter/ukr_migration_bounds.rds       <- used by 10
#          data_inter/ukr_migration_readings.rds     <- used by 15 and the docs
#            (eurostat and register_stock are the stocks with the breaks
#            back-cast, the stocks the register reading differences)
#          data_inter/ukr_migration_breaks.rds       <- the breaks, for the docs
#            list: flagged  country, month, change, ordinary, step, break
#                           (every month Eurostat flags, country totals)
#                  breaks   country, month, sex, age, change, ordinary, step
#                  stock    country, year, published, corrected (all sexes
#                           and ages; "EU + EFTA" is the sum over countries)
#                  register year, sex, register (as used), register_published
#                           (the clean-ups differenced as returns of the year)
# ==============================================================================

# ==============================================================================
# CANADA AND USA
# ==============================================================================
# The register stock includes Ukrainians who went to Canada and the USA. The
# UNHCR Data Finder does not count them properly (see section 1c), so they
# come from national programme figures. Neither country publishes a year-end
# count of Ukrainians present, so section 1c interpolates between the counts
# that are published, stored with their sources in
#   data_input/migration/national_programme_counts.csv
# 13c moves the two series 30% up and down to show what they carry.
#
# A better series, if it appears: CUAET arrivals and holders present at each
# year end (IRCC key-figures page history in the Wayback Machine, or IRCC open
# data), and Uniting for Ukraine parolees by month (DHS Office of Homeland
# Security Statistics). New rows in the CSV need no code change; re-run 06,
# 10, 11 onwards.
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

dt2_published <-
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

# --- breaks: clean-ups of invalid statuses ---------------------------------------
# A temporary protection status stays on the books of the country that granted
# it until that country ends it, and people who left, for Ukraine or for
# another country, were not always taken off. Some countries removed such
# invalid statuses in one clean-up. Eurostat flags each clean-up as a break in
# the monthly series (flag "b" on every cell of the country in that month),
# and its metadata (migr_asytp_esms) gives the reason: invalid statuses kept
# on the books until the clean-up (stated there for the Netherlands,
# Bulgaria, Germany and Italy; for Lithuania the source is the flag alone).
# Differenced as published, a clean-up reads as people returning to Ukraine
# in the year it happened. It is not a flow of that year: the people it
# removes had left the country of protection at unknown earlier dates and
# stayed on its books.
#
# The breaks follow from Eurostat's flags. A break is a month Eurostat flags
# in a country's series at which the country's total stock falls, once its
# ordinary monthly change (the median of its six monthly changes before that
# month) is taken off. Five are found:
#   Netherlands  January 2023
#   Bulgaria     February 2024
#   Germany      July 2024
#   Lithuania    October 2024
#   Italy        June 2025
# Denmark's flagged month, September 2025, is a rise and gets no correction.
# The flagged months and the breaks are asserted below, so a download with
# other flags stops the step until someone has looked at them.
#
# Each step is measured in the monthly series by sex and age group: the
# change from the month before the break to the break month, less the cell's
# ordinary monthly change, the median of its six monthly changes before the
# break. Only falls count; a cell that does not fall gets no correction. The
# invalid statuses are taken to have built up at a constant rate from the
# country's first month with a stock, so the stock at an earlier month t is
# lowered by
#   step x (t - first month) / (break month - first month)
# and from the break month on it stands as published. The year-end stocks are
# corrected cell by cell, by sex and age group: end 2022 for all five
# countries, end 2023 for all but the Netherlands, and end 2024 for Italy.
# The stock after the break is unchanged, so what moves is the timing of the
# outflow between years; the four-year total of the register reading moves
# only through the deaths abroad and the children aged 0, which are computed
# on the lowered stocks (section 6).
#
# Two large falls carry no break flag and stand as published: Poland in
# October 2022 and Czechia in March 2023. Poland's lies inside 2022, before
# the first year end, so it moves no end-of-year difference. Czechia's, were
# it treated as a break, would lower the 2022 outflow and raise 2023's.
#
# Spain, Greece and Cyprus still keep invalid statuses (313,090 people under
# temporary protection in the three at end 2025). They are not corrected: no
# clean-up has happened there, so no step measures them, and their stocks
# stand as published.

# months counted from year 0, so that differences are numbers of months
month_index <- function(m) as.integer(str_sub(m, 1, 4)) * 12L + as.integer(str_sub(m, 6, 7))

# the monthly series' countries, sexes and age groups, matched to the annual
# extract's; Y_LT18 overlaps the first two age groups and is not used
tp_geo_labels <- c(
  AT = "Austria", BE = "Belgium", BG = "Bulgaria", CH = "Switzerland",
  CY = "Cyprus", CZ = "Czechia", DE = "Germany", DK = "Denmark", EE = "Estonia",
  EL = "Greece", ES = "Spain", FI = "Finland", FR = "France", HR = "Croatia",
  HU = "Hungary", IE = "Ireland", IS = "Iceland", IT = "Italy",
  LI = "Liechtenstein", LT = "Lithuania", LU = "Luxembourg", LV = "Latvia",
  MT = "Malta", NL = "Netherlands", NO = "Norway", PL = "Poland",
  PT = "Portugal", RO = "Romania", SE = "Sweden", SI = "Slovenia", SK = "Slovakia"
)
tp_sex_codes <- c(F = "f", M = "m", UNK = "u")
tp_age_labels <- c(
  Y_LT14 = "Less than 14 years", `Y14-17` = "From 14 to 17 years",
  `Y18-34` = "From 18 to 34 years", `Y35-64` = "From 35 to 64 years",
  Y_GE65 = "65 years or over", UNK = "Unknown"
)

tp_monthly <-
  read_csv("data_input/refugees_eurostat/migr_asytpsm_api_UA_monthly_by_country_sex_age.csv",
           show_col_types = FALSE) |>
  filter(citizen == "UA", unit == "PER", freq == "M", geo != "EU27_2020") |>
  transmute(geo, ctr = unname(tp_geo_labels[geo]), sex, age, month = TIME_PERIOD,
            t = month_index(TIME_PERIOD), value = OBS_VALUE, flag = coalesce(OBS_FLAG, ""))

stopifnot(!anyNA(tp_monthly$ctr), all(tp_geo_labels %in% dt2_published$ctr))

# the change at a given month and the ordinary monthly change before it, by
# country and cell
step_at <- function(series, months) {
  series |>
    inner_join(months, by = "geo") |>
    filter(t >= t_break - 7, t <= t_break) |>
    arrange(geo, sex, age, t) |>
    mutate(change = value - lag(value), .by = c(geo, sex, age)) |>
    summarise(
      n_months = n(),
      change_at_break = change[t == t_break],
      ordinary = median(change[t >= t_break - 6 & t < t_break]),
      .by = c(geo, ctr, break_month, t_break, sex, age)
    ) |>
    rename(change = change_at_break) |>
    mutate(step = change - ordinary)
}

# the months Eurostat flags as breaks, one per country
tp_flagged <-
  tp_monthly |>
  filter(str_detect(flag, "b")) |>
  distinct(geo, break_month = month) |>
  mutate(t_break = month_index(break_month)) |>
  arrange(t_break, geo)

stopifnot(identical(
  paste(tp_flagged$geo, tp_flagged$break_month),
  c("NL 2023-01", "BG 2024-02", "DE 2024-07", "LT 2024-10", "IT 2025-06", "DK 2025-09")
))

tp_flag_steps <-
  step_at(tp_monthly |> filter(sex == "T", age == "TOTAL"), tp_flagged) |>
  mutate(is_break = step < 0) |>
  arrange(t_break, geo)

print(as.data.frame(tp_flag_steps |> select(ctr, break_month, change, ordinary, step, is_break)))

tp_breaks <-
  tp_flag_steps |>
  filter(is_break) |>
  select(geo, ctr, break_month, t_break)

stopifnot(identical(
  paste(tp_breaks$geo, tp_breaks$break_month),
  c("NL 2023-01", "BG 2024-02", "DE 2024-07", "LT 2024-10", "IT 2025-06")
))

# first month with a stock, from each country's total (Germany's series
# starts in August 2022, Italy's stands at zero in March 2022)
tp_first <-
  tp_monthly |>
  filter(geo %in% tp_breaks$geo, sex == "T", age == "TOTAL", value > 0) |>
  summarise(t_first = min(t), .by = geo)

tp_steps <-
  step_at(
    tp_monthly |> filter(sex %in% names(tp_sex_codes), age %in% names(tp_age_labels)),
    tp_breaks |> select(geo, break_month, t_break)
  ) |>
  mutate(
    step = if_else(change < 0, base::pmin(step, 0), 0),
    sex = unname(tp_sex_codes[sex]),
    age = unname(tp_age_labels[age])
  ) |>
  left_join(tp_first, by = "geo")

stopifnot(
  nrow(tp_steps) == nrow(tp_breaks) * length(tp_sex_codes) * length(tp_age_labels),
  all(tp_steps$n_months == 8),
  !anyNA(tp_steps)
)

# The annual extract's December cells are the monthly series' December cells
# for the countries with a break, so the steps lower the stocks they were
# measured on. The one exception is the Netherlands at end 2022, which
# Eurostat revised after the extract was taken: the extract holds 22,340
# people of unknown sex there and the monthly series 10, most of them
# assigned to men, and its total is 3,165 higher. The steps are applied to
# the extract's cells as they stand.
tp_december <-
  tp_monthly |>
  filter(geo %in% tp_breaks$geo, str_sub(month, 6, 7) == "12",
         sex %in% names(tp_sex_codes), age %in% names(tp_age_labels)) |>
  transmute(ctr, year = as.integer(str_sub(month, 1, 4)),
            sex = unname(tp_sex_codes[sex]), age = unname(tp_age_labels[age]),
            monthly = value) |>
  left_join(dt2_published, by = c("ctr", "year", "sex", "age")) |>
  filter(!(ctr == "Netherlands" & year == 2022))
stopifnot(nrow(tp_december) == (nrow(tp_breaks) * length(yrs) - 1) * 18,
          isTRUE(all.equal(tp_december$monthly, tp_december$mix)))

# share of a step that had built up by month t: none at the first month, all
# of it at the break month; from the break month on the stock is as published
built_up <- function(t, t_first, t_break) {
  if_else(t < t_break, base::pmax(t - t_first, 0) / (t_break - t_first), 0)
}

tp_corr_cells <-
  tp_steps |>
  cross_join(tibble(year = yrs)) |>
  mutate(corr = step * built_up(month_index(paste0(year, "-12")), t_first, t_break)) |>
  select(ctr, year, sex, age, corr)

# the extract's totals over sex ("t") and over age ("Total") take the sum of
# their cells' corrections
tp_corr_sex <-
  bind_rows(
    tp_corr_cells,
    tp_corr_cells |> summarise(corr = sum(corr), .by = c(ctr, year, age)) |> mutate(sex = "t")
  )
tp_corr <-
  bind_rows(
    tp_corr_sex,
    tp_corr_sex |> summarise(corr = sum(corr), .by = c(ctr, year, sex)) |> mutate(age = "Total")
  )

stopifnot(nrow(semi_join(tp_corr, dt2_published, by = c("ctr", "year", "sex", "age"))) ==
            nrow(tp_corr))

dt2 <-
  dt2_published |>
  left_join(tp_corr, by = c("ctr", "year", "sex", "age")) |>
  mutate(mix = mix + coalesce(corr, 0)) |>
  select(-corr)

# no cell goes below zero, and from each break on the stock is as published
stopifnot(
  all(dt2$mix >= 0),
  dt2 |>
    inner_join(dt2_published, by = c("ctr", "year", "sex", "age"), suffix = c("", "_pub")) |>
    inner_join(tp_breaks |> select(ctr, break_month), by = "ctr") |>
    filter(year >= as.integer(str_sub(break_month, 1, 4))) |>
    with(isTRUE(all.equal(mix, mix_pub)))
)

print(as.data.frame(
  tp_steps |>
    summarise(step = sum(step), change = sum(change), .by = c(ctr, break_month, sex)) |>
    arrange(ctr, sex)
))

# --- year-end totals and age-sex structure --------------------------------------
eurostat_total <- function(d) {
  d |>
    filter(ctr != eu_row, sex == "t", age == "Total") |>
    summarise(eurostat = sum(mix), .by = year)
}

# unknown sex and age are allocated with the distribution of the known cases
eurostat_structure <- function(d) {
  d |>
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
}

eurostat_tot <- eurostat_total(dt2)
eurostat_bands <- eurostat_structure(dt2)

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
# 1c. CANADA AND USA
# ==============================================================================
# Why not the Data Finder: it counts ~15-22k Ukrainians in the USA (asylum
# seekers and refugees only; U4U parolees are not refugees, and US arrivals
# were ~260k), and holds Canada flat at ~75k in 2022-2024 before jumping to
# 174k in 2025, a reporting change rather than a flow.
#
# Linear interpolation between dated counts, evaluated at 31 December, flat
# after the last one. What the register stock needs is
# PEOPLE PRESENT, so each series ends on a published count of people present
# and uses programme arrivals only for the earlier years, when the programmes
# had just opened and few had yet left, so arrivals still approximate a stock:
#
#   CANADA  298,128 cumulative CUAET entries to 1 Apr 2024 (entry closed
#           31 Mar 2024), then 274,000 people in Canada at May 2025
#           (CES 2026 p. 12). The fall is departure and onward movement, which
#           the earlier arrivals-only series could not show at all.
#   USA     93,928 U4U arrivals to Dec 2022 and over 158,000 to Sep 2023
#           (DHS assessment of U4U, Nov 2024), then 260,000 people present.
#           CES dates that count at Nov 2025, but also reports (p. 22) that
#           U4U "was effectively suspended for new applications in early
#           2025". The anchor is therefore placed at 31 Jan 2025: with intake
#           stopped, the stock had plateaued by then, and dating it at Nov
#           2025 would spread the rise over two years in which it could not
#           have happened.
#
# The two bases are spliced, so the Canadian decline is partly a change of
# definition rather than pure return; see documents/migration_methodology.md
# section 8.
anchors <-
  read_csv("data_input/migration/national_programme_counts.csv",
           show_col_types = FALSE) |>
  mutate(date = as.Date(date))

programme <-
  anchors |>
  group_by(country) |>
  reframe(
    year = yrs,
    present = approx(
      x = as.numeric(date), y = value,
      xout = as.numeric(as.Date(paste0(yrs, "-12-31"))),
      rule = 2
    )$y
  ) |>
  pivot_wider(names_from = country, values_from = present) |>
  rename(canada = CAN, usa = USA)

register_stock <-
  eurostat_tot |>
  left_join(unhcr_other, by = "year") |>
  left_join(programme, by = "year") |>
  mutate(register_stock = eurostat + unhcr_other + canada + usa)

print(as.data.frame(register_stock))

# ==============================================================================
# 2. AGE-SEX NET FLOW BY COHORT DIFFERENCING
# ==============================================================================
# The destinations outside EU + EFTA carry no age-sex detail, so the register
# stock takes Eurostat's structure for that year. Ungrouping is done on the
# stocks, which are positive, and only then differenced.
# The stocks are counts of persons and enter pclm() as they are: it chooses its
# penalty by the BIC of a Poisson likelihood, so the size of the counts is what
# says how far the single ages may depart from a smooth curve within a band.
# Inflating them first (they were once multiplied by 1e5) leaves the fit almost
# unsmoothed and moves men within the band of ages 35 to 64 towards its end.
ung_stock <- function(d) {
  y <- base::pmax(d$mix, 1)
  fit <- suppressWarnings(pclm(x = d$age, y = y, nlast = 36)$fitted)
  tibble(age = 0:100, S = fit)
}

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

# from a year-end age-sex structure and register stock to the flows by single
# age: scaled, ungrouped by year and sex, and differenced by cohort (section 6
# applies it to the stock as published as well)
register_flows <- function(bands, stock) {
  bands |>
    left_join(stock |> select(year, register_stock), by = "year") |>
    mutate(mix = cx * register_stock) |>
    group_by(year, sex) |>
    group_modify(~ ung_stock(.x)) |>
    ungroup() |>
    group_by(sex) |>
    group_modify(~ cohort_diff(.x)) |>
    ungroup() |>
    left_join(mx_exp, by = c("year", "sex", "age")) |>
    mutate(
      deaths_abroad = if_else(year == min(yrs), 0, coalesce(S_prev * mx, 0)),
      D = if_else(age == 0, 0, S - coalesce(S_prev, 0) + deaths_abroad)
    )
}

flows <- register_flows(eurostat_bands, register_stock)

stopifnot(!any(is.na(flows$D)))

# Each year's stock is ungrouped on its own, so where a cohort crosses a
# Eurostat band edge between two year ends, a change in the within-band shape
# shows up as a flow. The flows are smoothed lightly across age by penalised
# least squares on second differences, with each year's total by sex and
# Eurostat band (ages 1-13, 14-17, 18-34, 35-64, 65+) held fixed, so every
# migration total and both readings stay as they are. Age 0 stays at zero.
# At this strength no single age moves by more than 270 people, and few by
# more than a hundred: the flows are already smooth from year to year of age. What remains - the
# step from zero at age 0 to the return of young children at age 1, and the
# broad waves around 35 and 65, where the older women's returns in the 65+
# band are spread by the ungrouping into the ages around it - is set by the
# band totals themselves, which no smoothing that keeps them can remove. The
# flows before smoothing are kept as D_raw.
lambda_D <- 5
smooth_flows <- function(D, age, lambda = lambda_D) {
  keep <- age >= 1
  x <- D[keep]
  n <- length(x)
  band <- cut(age[keep], c(1, 14, 18, 35, 65, Inf), right = FALSE)
  C <- t(model.matrix(~ band - 1))
  P <- diff(diag(n), differences = 2)
  A <- rbind(cbind(diag(n) + lambda * crossprod(P), t(C)),
             cbind(C, matrix(0, nrow(C), nrow(C))))
  out <- D
  out[keep] <- solve(A, c(x, C %*% x))[seq_len(n)]
  out
}
flows <- flows |>
  arrange(year, sex, age) |>
  mutate(D_raw = D, D = smooth_flows(D, age), .by = c(year, sex))
stopifnot(isTRUE(all.equal(
  flows |> summarise(d = sum(D), .by = c(year, sex)) |> pull(d),
  flows |> summarise(d = sum(D_raw), .by = c(year, sex)) |> pull(d)
)))

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
# 6. THE CLEAN-UPS READ AS RETURNS
# ==============================================================================
# The register reading from the stock as published, where each clean-up is a
# fall in the year it happened, beside the reading as used. Both end on the
# same end-2025 stock, so the four-year totals differ only through the deaths
# abroad added back and the children aged 0 left out, which are computed on
# the lower year-end stocks before the breaks.
register_stock_published <-
  register_stock |>
  select(-eurostat, -register_stock) |>
  left_join(eurostat_total(dt2_published), by = "year") |>
  mutate(register_stock = eurostat + unhcr_other + canada + usa)

flows_published <- register_flows(eurostat_structure(dt2_published), register_stock_published)

register_by_sex <-
  flows |>
  summarise(register = sum(D), .by = c(year, sex)) |>
  left_join(
    flows_published |> summarise(register_published = sum(D), .by = c(year, sex)),
    by = c("year", "sex")
  ) |>
  arrange(year, sex)

stock_by_country <- function(d, nm) {
  bind_rows(
    d |>
      filter(ctr %in% tp_breaks$ctr, sex == "t", age == "Total") |>
      select(country = ctr, year, value = mix),
    eurostat_total(d) |> transmute(country = "EU + EFTA", year, value = eurostat)
  ) |>
    rename("{nm}" := value)
}

migration_breaks <- list(
  flagged = tp_flag_steps |>
    select(country = ctr, month = break_month, change, ordinary, step, "break" = is_break),
  breaks = tp_steps |>
    select(country = ctr, month = break_month, sex, age, change, ordinary, step),
  stock = stock_by_country(dt2_published, "published") |>
    left_join(stock_by_country(dt2, "corrected"), by = c("country", "year")) |>
    arrange(country, year),
  register = register_by_sex
)

stopifnot(
  isTRUE(all.equal(register_stock$register_stock[register_stock$year == max(yrs)],
                   register_stock_published$register_stock[register_stock_published$year == max(yrs)])),
  isTRUE(all.equal(sum(register_by_sex$register), sum(readings$register)))
)

write_rds(migration_breaks, "data_inter/ukr_migration_breaks.rds")

print(as.data.frame(migration_breaks$stock))
print(as.data.frame(register_by_sex |>
  bind_rows(register_by_sex |> summarise(across(c(register, register_published), sum),
                                         .by = sex) |> mutate(year = 0L))))

# ==============================================================================
# DIAGNOSTIC FIGURES (exploratory)
# ==============================================================================
nice <- function(d) d |> mutate(sex = if_else(sex == "f", "Females", "Males"))

# register stock and its parts
register_stock |>
  select(year, eurostat, unhcr_other, canada, usa) |>
  pivot_longer(-year, names_to = "part", values_to = "stock") |>
  mutate(part = factor(part, levels = c("usa", "canada",
                                        "unhcr_other", "eurostat"))) |>
  ggplot(aes(factor(year), stock / 1e6, fill = part)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(eurostat = "#0A9396", unhcr_other = "#94D2BD",
                               canada = "#EE9B00", usa = "#CA6702")) +
  labs(x = NULL, y = "Stock at 31 December (millions)", fill = NULL,
       caption = paste("Canada and USA interpolated between published counts.",
                       "Eurostat with the flagged clean-ups back-cast.")) +
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
