# ==============================================================================
# STEP 13r - Sensitivity: the register stock built bottom-up
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# 06 takes the register reading of the western flow from UNHCR's published
# year-end totals of refugees from Ukraine recorded outside Russia. This step
# builds the same stock bottom-up, so that each part is counted in its own
# source's definitions, as the input was specified before UNHCR's totals were
# used, and reports how far the life expectancy loss moves with the choice of
# stock. It reproduces that specification exactly, and its stock and reading
# are asserted against the numbers it gave.
#
#   1. BOTTOM-UP STOCK of people displaced outside Russia and Belarus, at the
#      end of each year:
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
#   2. THE REGISTER READING of that stock, by the same cohort differencing,
#      deaths abroad, age 0 and smoothing as 06 (the mig_* functions of
#      00_setup.R), with Eurostat's age-sex structure with the breaks
#      back-cast; and the same from the stock as published, where each
#      clean-up is a fall in its own year.
#
#   3. THE CLEAN-UPS READ AS RETURNS: the measured steps, the year-end stocks
#      before and after the correction, and the reading by sex from both.
#
#   4. RUSSIA AND BELARUS as one 2022 component, the counterpart of the
#      bottom-up stock, which leaves out both countries: mode = Data Finder
#      stock of the two at end 2022 net of 2021, min 700,000, max 2,000,000;
#      and the allocation at the mode of the bottom-up input.
#
#   5. SCENARIOS. Each one deterministic projection with every conflict input
#      at its central value, as 13c section 5 does: the main input at its mode
#      and at the register end of the bracket; the bottom-up input at its mode
#      and register end; and the bottom-up input with the stock as published.
#      An input with a register reading of its own takes its own age-sex
#      allocation.
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
#          data_inter/ukr_migration_readings.rds, ukr_migration_bounds.rds (06)
#          the same as 13, and data_inter/ukr_migration_decomposition.rds (13)
# OUTPUTS  data_inter/ukr_bottom_up_stock_e0.rds   scenario, year, sex, loss
#          data_inter/ukr_bottom_up_stock.rds      <- used by 15
#            list: stock    year, eurostat (breaks back-cast),
#                           eurostat_published, unhcr_other, canada, usa,
#                           bottom_up, bottom_up_published, unhcr_total,
#                           difference (unhcr_total - bottom_up)
#                  readings year, crossings, register_unhcr (06),
#                           register_bottom_up, register_bottom_up_published
#                  russia   component, min, mode, max: "Russia" (06) and
#                           "Russia and Belarus" (the bottom-up input)
#          data_inter/ukr_migration_breaks.rds     <- the breaks, for the docs
#            list: flagged  country, month, change, ordinary, step, break
#                           (every month Eurostat flags, country totals)
#                  breaks   country, month, sex, age, change, ordinary, step
#                  stock    country, year, published, corrected (all sexes
#                           and ages; "EU + EFTA" is the sum over countries)
#                  register year, sex, register (breaks back-cast),
#                           register_published (the clean-ups differenced as
#                           returns of the year)
#          figures/exploratory/migration_register_stock_components.png
# ==============================================================================

# ==============================================================================
# CANADA AND USA
# ==============================================================================
# The bottom-up stock includes Ukrainians who went to Canada and the USA. The
# UNHCR Data Finder does not count them properly (see section 1c), so they
# come from national programme figures. Neither country publishes a year-end
# count of Ukrainians present, so section 1c interpolates between the counts
# that are published, stored with their sources in
#   data_input/migration/national_programme_counts.csv
#
# A better series, if it appears: CUAET arrivals and holders present at each
# year end (IRCC key-figures page history in the Wayback Machine, or IRCC open
# data), and Uniting for Ukraine parolees by month (DHS Office of Homeland
# Security Statistics). New rows in the CSV need no code change; re-run 13r.
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

yrs <- mig_years

# ==============================================================================
# 1a. EUROSTAT TEMPORARY PROTECTION STOCK, EU + EFTA, 31 DECEMBER
# ==============================================================================
# The extract holds December of each year for EU-27 countries plus Iceland,
# Liechtenstein, Norway and Switzerland (mig_eurostat_extract(), 00_setup.R).
dt2_published <- mig_eurostat_extract(yrs)

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
# on the lowered stocks (section 3).
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

# --- year-end totals and age-sex structure (00_setup.R) -------------------------
eurostat_tot <- mig_eurostat_total(dt2)
eurostat_bands <- mig_eurostat_structure(dt2)

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
# definition rather than pure return; see documents/reports/migration_methodology.md
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

# the bottom-up stock as it was specified
stopifnot(identical(round(register_stock$register_stock), c(4304675, 4846389, 5259885, 5499840)))

# ==============================================================================
# 2. THE REGISTER READING OF THE BOTTOM-UP STOCK
# ==============================================================================
# The destinations outside EU + EFTA carry no age-sex detail, so the register
# stock takes Eurostat's structure for that year, here with the breaks
# back-cast. The rest is 06's construction (00_setup.R).
mx_exp <- mig_mx_exp(yrs)

flows <- mig_register_flows(eurostat_bands, register_stock, mx_exp, yrs)
stopifnot(!any(is.na(flows$D)))
flows <- mig_smoothed(flows)

crossings <- mig_crossings()

readings <-
  mig_flow_totals(flows) |>
  left_join(crossings, by = "year") |>
  left_join(register_stock, by = "year") |>
  left_join(unhcr_ru_by, by = "year")

print(as.data.frame(readings |> select(year, register, crossings, departures, returns,
                                       deaths_abroad, births_abroad_dropped, ru_by)))

# the register reading as it was specified
stopifnot(identical(round(readings$register, 1), c(4227661.1, 488934.5, 369209.8, 203317.9)))

# ==============================================================================
# 3. THE CLEAN-UPS READ AS RETURNS
# ==============================================================================
# The register reading from the stock as published, where each clean-up is a
# fall in the year it happened, beside the reading with the breaks back-cast.
# Both end on the same end-2025 stock, so the four-year totals differ only
# through the deaths abroad added back and the children aged 0 left out, which
# are computed on the lower year-end stocks before the breaks.
register_stock_published <-
  register_stock |>
  select(-eurostat, -register_stock) |>
  left_join(mig_eurostat_total(dt2_published), by = "year") |>
  mutate(register_stock = eurostat + unhcr_other + canada + usa)

flows_published <- mig_register_flows(mig_eurostat_structure(dt2_published),
                                      register_stock_published, mx_exp, yrs)

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
    mig_eurostat_total(d) |> transmute(country = "EU + EFTA", year, value = eurostat)
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
  isTRUE(all.equal(sum(register_by_sex$register), sum(readings$register))),
  # the reading from the stock as published, as it was specified
  identical(
    register_by_sex |>
      filter(paste(year, sex) %in% c("2022 f", "2022 m", "2024 f")) |>
      pull(register_published) |>
      round(2),
    c(2887097.51, 1490206.08, -69050.84)
  )
)

write_rds(migration_breaks, "data_inter/ukr_migration_breaks.rds")

print(as.data.frame(migration_breaks$stock))
print(as.data.frame(register_by_sex |>
  bind_rows(register_by_sex |> summarise(across(c(register, register_published), sum),
                                         .by = sex) |> mutate(year = 0L))))

# ==============================================================================
# 4. RUSSIA AND BELARUS, AND THE ALLOCATION AT THE MODE
# ==============================================================================
# The bottom-up stock leaves out Russia and Belarus, so its counterpart is one
# 2022 component for the two. The mode is UNHCR's end-2022 statistical stock
# in the two, net of 2021. The floor allows for Russian administrative data
# overstating and for onward movement; the ceiling for the rental-contract
# basis missing non-renters, informal renters and the naturalised, and for the
# 1.5-2.0M Russian officials cited in autumn 2022.
ru_by_min <- 700000
ru_by_max <- 2000000

bounds_bottom_up <-
  bind_rows(
    readings |>
      transmute(year, component = "west", min = base::pmin(register, crossings),
                mode = (register + crossings) / 2, max = base::pmax(register, crossings)),
    readings |>
      filter(year == min(yrs)) |>
      transmute(year, component = "ru_by", min = ru_by_min, mode = ru_by, max = ru_by_max)
  )

stopifnot(bounds_bottom_up |> filter(component == "ru_by") |> pull(mode) == 1284402)

alloc_bottom_up <-
  mig_allocate(flows, bounds_bottom_up |> summarise(mode_tot = sum(mode), .by = year))

# the same from the stock as published: its own reading, smoothed as 06
# smooths its flows, and its own allocation at its mode
flows_published_s <- mig_smoothed(flows_published)
register_published <- mig_flow_totals(flows_published_s) |> select(year, register_published = register)
alloc_published <-
  mig_allocate(
    flows_published_s,
    register_published |>
      left_join(crossings, by = "year") |>
      left_join(bounds_bottom_up |> filter(component == "ru_by") |> select(year, ru_mode = mode),
                by = "year") |>
      transmute(year, mode_tot = (register_published + crossings) / 2 + coalesce(ru_mode, 0))
  )

# ==============================================================================
# 5. SCENARIOS
# ==============================================================================
# Every conflict input at its central value, as 13c section 5. The western
# flow takes the path of each scenario, and the Russian component its central
# value, the mean of its PERT (central_values(), 00_setup.R): Russia for the
# input as used, Russia and Belarus for the bottom-up input.
param_table <- read_rds("data_inter/ukr_param_table.rds")

exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") %>%
  filter(year %in% 2022:2025, source == "frcst") %>%
  select(-source)

pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") %>%
  filter(reg == "cnt", year == 2022) %>%
  select(-reg)

asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds")
asfr2 <- asfr %>%
  bind_rows(
    asfr %>% filter(year == 2023) %>% mutate(year = 2024),
    asfr %>% filter(year == 2023) %>% mutate(year = 2025)
  ) %>%
  arrange(year, age)

ohchr2 <- read_rds("data_inter/ukr_ohchr_civilian_casualties.rds") %>%
  select(year, sex, age, prop_cvs = cx)

ual_raw <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")

cmb_profile <- function(keep_status, nm) {
  ual_raw |>
    filter(status == keep_status, year %in% 2022:2025) |>
    select(-status) |>
    complete(year = 2022:2025, sex, age = 0:100, fill = list(dx = 0)) |>
    summarise(dx = sum(dx), .by = c(year, sex, age)) |>
    mutate("{nm}" := dx / sum(dx), .by = year) |>
    select(-dx)
}

profile_props <-
  cmb_profile("dead", "prop_cmb_dead") |>
  left_join(cmb_profile("missing", "prop_cmb_miss"), by = c("year", "sex", "age")) |>
  left_join(ohchr2, by = c("year", "sex", "age"))

# the static inputs of 13c, with the age-sex allocation of the input given
static_from <- function(alloc) {
  migs2 <- alloc |> mutate(ems = -mix) |> select(-mix)
  out <-
    expand_grid(year = 2022:2025, sex = c("f", "m"), age = 0:100) %>%
    left_join(migs2, by = c("year", "sex", "age")) %>%
    left_join(exp_mort2, by = c("year", "sex", "age")) %>%
    left_join(asfr2, by = c("year", "sex", "age")) %>%
    left_join(profile_props, by = c("year", "sex", "age")) %>%
    replace_na(list(ems = 0, w = 0, mx = 0, fx = 0, prop_cvs = 0,
                    prop_cmb_dead = 0, prop_cmb_miss = 0))
  stopifnot(all(out$mx > 0), !any(is.na(out)))
  out
}

static_main <- static_from(read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds"))
static_bottom_up <- static_from(alloc_bottom_up)
static_published <- static_from(alloc_published)

conflict_at_mode <-
  param_table |>
  central_values() |>
  filter(role %in% c("civilians", "combatants")) |>
  select(year, role, mode) |>
  pivot_wider(names_from = role, values_from = mode) |>
  rename(draw_cvs = civilians, draw_cmb = combatants) |>
  left_join(
    param_table |> filter(role == "combatants") |> select(year, conf_cmb = confirmed),
    by = "year"
  )

loss_of <- function(static, draw_mig) {
  sim <-
    run_single_sim(1, conflict_at_mode |> left_join(draw_mig, by = "year") |> mutate(sim_id = 1),
                   static, pop22_ini) |>
    mutate(conflict = civilian + combatant_confirmed + combatant_imputed,
           mx_all = (expected + conflict) / pop, mx_bsn = expected / pop)
  stopifnot(min(sim$pop) >= 0)
  e0 <- function(col) {
    sim |> select(year, sex, age, mx = all_of(col)) |> group_by(year, sex) |>
      do(lifetable(dt_in = .data)) |> ungroup() |> filter(age == 0) |> select(year, sex, ex)
  }
  e0("mx_bsn") |> rename(ex_bsn = ex) |>
    left_join(e0("mx_all") |> rename(ex_all = ex), by = c("year", "sex")) |>
    transmute(year, sex, loss = ex_bsn - ex_all)
}

# the Russian component at its central value: Russia (06), or Russia and
# Belarus with the bounds above
ru_central <- param_table |> filter(role == "mig_ru") |> central_values() |> transmute(year, value = mode)
ru_by_central <-
  bounds_bottom_up |>
  filter(component == "ru_by") |>
  transmute(year, role = "mig_ru", min, mode, max) |>
  central_values() |>
  transmute(year, value = mode)

# a western path plus a Russian component, as the projection's yearly totals
totals <- function(west, ru) bind_rows(west, ru) |> summarise(draw_mig = sum(value), .by = year)

west_main <- param_table |> filter(role == "mig_west")
west_bottom_up <- readings |> select(year, crossings, register)

scenario_order <- c(
  "UNHCR totals: mode (as used)",
  "UNHCR totals: register end",
  "Bottom-up stock: mode",
  "Bottom-up stock: register end",
  "Bottom-up stock, Eurostat as published: mode"
)

bottom_up_e0 <- bind_rows(
  loss_of(static_main,
          totals(west_main |> transmute(year, value = (crossings + register) / 2), ru_central)) |>
    mutate(scenario = scenario_order[1]),
  loss_of(static_main, totals(west_main |> transmute(year, value = register), ru_central)) |>
    mutate(scenario = scenario_order[2]),
  loss_of(static_bottom_up,
          totals(west_bottom_up |> transmute(year, value = (crossings + register) / 2), ru_by_central)) |>
    mutate(scenario = scenario_order[3]),
  loss_of(static_bottom_up, totals(west_bottom_up |> transmute(year, value = register), ru_by_central)) |>
    mutate(scenario = scenario_order[4]),
  loss_of(static_published,
          totals(register_published |> left_join(crossings, by = "year") |>
                   transmute(year, value = (crossings + register_published) / 2), ru_by_central)) |>
    mutate(scenario = scenario_order[5])
) |>
  select(scenario, year, sex, loss)

# the first row must reproduce 13's loss at the mode
stopifnot(isTRUE(all.equal(
  bottom_up_e0 |> filter(scenario == scenario_order[1]) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

# ==============================================================================
# 6. OUTPUTS
# ==============================================================================
readings_main <- read_rds("data_inter/ukr_migration_readings.rds")
bounds_main <- read_rds("data_inter/ukr_migration_bounds.rds")

bottom_up_stock <- list(
  stock = register_stock |>
    left_join(mig_eurostat_total(dt2_published) |> rename(eurostat_published = eurostat), by = "year") |>
    left_join(register_stock_published |> select(year, bottom_up_published = register_stock), by = "year") |>
    left_join(readings_main |> select(year, unhcr_total = register_stock), by = "year") |>
    transmute(year, eurostat, eurostat_published, unhcr_other, canada, usa,
              bottom_up = register_stock, bottom_up_published, unhcr_total,
              difference = unhcr_total - bottom_up),
  readings = readings |>
    select(year, crossings, register_bottom_up = register) |>
    left_join(readings_main |> select(year, register_unhcr = register), by = "year") |>
    left_join(register_published, by = "year") |>
    transmute(year, crossings, register_unhcr, register_bottom_up,
              register_bottom_up_published = register_published),
  russia = bind_rows(
    bounds_main |> filter(component == "ru") |> mutate(component = "Russia"),
    bounds_bottom_up |> filter(component == "ru_by") |> mutate(component = "Russia and Belarus")
  ) |>
    select(component, min, mode, max)
)

write_rds(bottom_up_e0, "data_inter/ukr_bottom_up_stock_e0.rds")
write_rds(bottom_up_stock, "data_inter/ukr_bottom_up_stock.rds")

cat("\n=== REGISTER STOCK AT 31 DECEMBER: UNHCR'S TOTALS AND THE BOTTOM-UP STOCK ===\n")
print(as.data.frame(bottom_up_stock$stock |> mutate(across(-year, round))))
cat("\n=== READINGS OF THE WESTERN FLOW ===\n")
print(as.data.frame(bottom_up_stock$readings |> mutate(across(-year, round))))
cat("\n=== THE RUSSIAN COMPONENT ===\n")
print(as.data.frame(bottom_up_stock$russia))

cat("\n=== LOSS BY STOCK: MEN BY YEAR, AND 2025 WOMEN ===\n")
print(as.data.frame(
  bottom_up_e0 |>
    mutate(col = if_else(sex == "m", paste0("m_", year), paste0("f_", year))) |>
    filter(sex == "m" | year == 2025) |>
    select(scenario, col, loss) |>
    pivot_wider(names_from = col, values_from = loss) |>
    mutate(across(where(is.numeric), \(x) round(x, 3)))
))

# bottom-up stock and its parts
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

message("Done. data_inter/ukr_bottom_up_stock.rds and ukr_bottom_up_stock_e0.rds written for 15.")
