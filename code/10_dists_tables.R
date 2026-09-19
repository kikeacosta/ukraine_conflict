# ==============================================================================
# STEP 10 - Parameter table for the Monte Carlo simulation
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Turns the conflict death estimates from the previous steps into the
# min / mode / max triplets that step 11 feeds to a PERT distribution, one
# per year and role.
#
#   COMBATANTS (ualosses register)
#     The total is the registered dead, completed for registration lag, plus
#     the missing imputed as dead by 09's model. How many of the missing are
#     alive rests on the evidence alive_evidence() (00_setup.R) sets out: the
#     prisoners of war among the missing and the share of the unresolved
#     alive for other reasons. The bounds are the total at the ends and the
#     centre of that evidence, at the model's estimates and 08b's point
#     factors (military_draws(), 00_setup.R):
#     min  = the most alive: every unrecorded prisoner among the missing, the
#            most prisoners held, and the unresolved alive at their ceiling
#     mode = the central values
#     max  = the fewest alive: the unrecorded prisoners among the missing as
#            the returned of every year had been, the fewest held, none of
#            the unresolved alive
#     11 draws the evidence together with the sampling error of the model
#     and of the lag factors, so a draw can fall slightly outside these
#     bounds. Two further columns describe the register rather than the
#     estimate: confirmed (its dead with their late registrations, used to
#     split confirmed from imputed) and listed (confirmed plus everyone still
#     listed as missing).
#
#   CIVILIANS (UCDP georeferenced event data)
#     min / mode / max = the low / best / high bounds UCDP publishes,
#     after the deaths of unknown side have been redistributed (step 07)
#
# INPUTS   data_inter/ukr_ucdp_invals.rds                     (from 07_invals)
#          data_inter/ukr_ualosses_..._sex_age_2022_2025.rds  (from 08)
#          data_inter/ukr_ualosses_..._imputed_...rds         (from 09)
#          data_inter/ukr_military_inputs.rds (09): the evidence on the missing
#          alive, the model and the registered counts by event month
#          data_inter/ukr_migration_bounds.rds                (from 06)
#          data_input/ohchr_civilian_deaths.xlsx (sheet "monthly"),
#          data_input/migration/unhcr/unhcr_border_crossings_2022.csv and
#          data_inter/ukr_registration_completion.rds (08b): the months of
#          2022's events, for the timing within the year
# OUTPUT   data_inter/ukr_param_table.rds
#          data_inter/ukr_timing_absent.rds: the share of each year's net
#          outflow and deaths removed before exposure, read by
#          run_single_sim() (00_setup.R)
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

ucdp <- read_rds("data_inter/ukr_ucdp_invals.rds")
ual <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")
ual_imp <- read_rds(
  "data_inter/ukr_ualosses_conflict_deaths_imputed_miss_sex_age_2022_2025.rds"
)

# complete the age-sex grid so the three bounds are built on the same support
ual2 <-
  ual %>%
  # dead and missing only: prisoners and released prisoners are alive
  filter(status %in% c("dead", "missing"), year %in% 2022:2025) %>%
  summarise(dx = sum(dx), .by = c(status, year, sex, age)) |>
  arrange(status, year, sex, age) |>
  complete(status, year = 2022:2025, sex, age = 0:100, fill = list(dx = 0))

# --- combatants -------------------------------------------------------------
mil <- read_rds("data_inter/ukr_military_inputs.rds")
ev <- mil$evidence
at_evidence <- function(captives, other) military_draws(mil, captives, other)
at_mode <- at_evidence(ev$captives_mode, ev$other_mode)

cmb <-
  at_mode |>
  transmute(
    year,
    role = "combatants",
    mode = military,
    min = at_evidence(ev$captives_max, ev$other_max)$military,
    max = at_evidence(ev$captives_min, ev$other_min)$military,
    confirmed,
    listed = confirmed + missing
  )

# the mode must be the total 09 imputed at the central values; the register
# columns must be 09's counts completed for registration lag (08b), and those
# must rest on the registered counts the profiles are built from
imp_tab <- read_rds("data_inter/ukr_ualosses_imputation_table.rds") |> arrange(year)
stopifnot(
  isTRUE(all.equal(
    cmb$mode,
    ual_imp |> summarise(d = sum(dx), .by = year) |> arrange(year) |> pull(d)
  )),
  isTRUE(all.equal(cmb$confirmed, imp_tab$confirmados_stock)),
  isTRUE(all.equal(cmb$listed, imp_tab$confirmados_stock + imp_tab$missing_stock)),
  isTRUE(all.equal(
    imp_tab$registered,
    ual2 |> filter(status == "dead") |> summarise(d = sum(dx), .by = year) |> arrange(year) |> pull(d)
  ))
)

# --- civilians --------------------------------------------------------------
cvs <-
  ucdp |>
  filter(role == "civilians") |>
  select(year, role, mode = dts, min = dts_l, max = dts_u)

# --- migration --------------------------------------------------------------
# Net emigration enters as TWO components, because two different things are
# unknown about it and they are unknown to different degrees:
#
#   mig_west   how much of the western outflow the sources see. The two
#              readings are CES's net border crossings and the cohort-
#              differenced register stock built in 06; the gap between them
#              is the discrepancy Pozniak (2023) set out - registers that keep
#              people after they return against a crossing balance that
#              missed the peak weeks.
#   mig_ru_by  displacement into Russia and Belarus, which no register sees;
#              one 2022 entry centred on UNHCR's end-2022 statistical stock.
#
# 06 builds the bounds from the input data and allocates its age-sex flows at
# the modes of the same table. Both components are systematic rather than
# year-by-year noise, and 11 draws them accordingly.
#
# The western component is bracketed by SOURCE, not year by year: its two
# readings disagree in opposite directions in different years (the register
# is higher in 2022-2023, the crossings in 2024-2025), so per-year bounds
# would mix the sources and reach four-year totals neither supports. What the
# population denominator depends on is the total abroad, so 11 blends the two
# readings' whole paths. The paths are carried here as `crossings` and
# `register`; min and max stay the per-year range of the two, for display.
migration_bounds <- read_rds("data_inter/ukr_migration_bounds.rds")
migration_readings <- read_rds("data_inter/ukr_migration_readings.rds")

mig <-
  migration_bounds |>
  mutate(role = paste0("mig_", component)) |>
  select(year, role, mode, min, max) |>
  left_join(
    migration_readings |> transmute(year, role = "mig_west", crossings, register),
    by = c("year", "role")
  )

# the western mode is the midpoint of the two paths, and its range their span
west <- mig |> filter(role == "mig_west")
stopifnot(
  isTRUE(all.equal(west$mode, (west$crossings + west$register) / 2)),
  isTRUE(all.equal(west$min, pmin(west$crossings, west$register))),
  isTRUE(all.equal(west$max, pmax(west$crossings, west$register)))
)

param_table <- bind_rows(cmb, cvs, mig)

# rpert() requires min <= mode <= max; catch any ordering problem here rather
# than as an obscure error inside the simulation loop
stopifnot(
  # 4 years x 2 conflict roles, 4 years of western migration, and the single
  # 2022 row for Russia and Belarus
  nrow(param_table) == 13,
  all(param_table$min <= param_table$mode),
  all(param_table$mode <= param_table$max)
)

# 06 allocated its flows at the modes, so a draw at the mode must leave the
# age-sex cells untouched. If these drift apart, run_single_sim() adds the
# difference along the departures profile in every draw, silently shifting
# every year's flows.
stopifnot(
  all.equal(
    param_table |>
      filter(str_starts(role, "mig_")) |>
      summarise(mode = sum(mode), .by = year) |>
      arrange(year) |>
      pull(mode),
    read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds") |>
      summarise(mix = -sum(mix), .by = year) |>
      arrange(year) |>
      pull(mix),
    tolerance = 1e-4
  )
)

print(as.data.frame(param_table))

write_rds(param_table, "data_inter/ukr_param_table.rds")

# totals implied by the table, for reference
param_table |>
  summarise(across(c(min, mode, max), sum), .by = role)

# --- timing within 2022 -------------------------------------------------------
# The projection removes a share of each year's net outflow, civilian and
# military deaths before counting the year's exposure (run_single_sim(),
# 00_setup.R). Half is the mid-year convention, and 2023-2025 keep it. 2022
# was not like that: the war began on 24 February, and most departures fell
# in March and April. A person who leaves or dies at a fraction u of the year
# lives u of it, so the share to remove is 1 - mean(u), with u taken at the
# middle of each month (26 February for the war's first days) and the mean
# weighted by the month's events:
#   civilian deaths   OHCHR's monthly count of civilians killed in 2022
#   military deaths   the register's dead and missing by month of the event
#                     (v19, from 08b)
#   net outflow       UNHCR's border crossings out of Ukraine less those into
#                     it, by month, as a share of the year's net; the western
#                     crossings stand for Russia and Belarus too, for which no
#                     monthly series exists
month_mid <- tibble(month = 2:12) |>
  mutate(mid = if_else(month == 2, as.Date("2022-02-26"), as.Date(sprintf("2022-%02d-15", month))),
         u = as.numeric(mid - as.Date("2022-01-01")) / 365)
mean_u <- function(d) d |> left_join(month_mid, by = "month") |> with(sum(n * u) / sum(n))
timing_civilian <-
  read_xlsx("data_input/ohchr_civilian_deaths.xlsx", sheet = "monthly") |>
  filter(Year == 2022) |>
  transmute(month = match(Month, month.abb), n = Killed)
timing_military <-
  read_rds("data_inter/ukr_registration_completion.rds")$month |>
  filter(year == 2022) |>
  summarise(n = sum(n_v19), .by = m) |>
  transmute(month = month(m), n)
# cumulative since 24 February; the first row, on 1 March, is the war's
# first days, entered as February
timing_crossings <-
  read_csv("data_input/migration/unhcr/unhcr_border_crossings_2022.csv", show_col_types = FALSE) |>
  arrange(date) |>
  mutate(net_cum = crossings_from_ukraine_cumulative - crossings_to_ukraine_cumulative,
         n = net_cum - lag(net_cum, default = 0),
         month = if_else(row_number() == 1, 2, month(date)))
stopifnot(all(timing_crossings$n > 0))
timing_2022 <- tibble(
  component = c("civilian deaths", "military deaths", "net outflow"),
  mean_u = c(mean_u(timing_civilian), mean_u(timing_military), mean_u(timing_crossings))
) |>
  mutate(absent = 1 - mean_u)
timing_absent <-
  tibble(year = 2022:2025, absent_cvs = 0.5, absent_cmb = 0.5, absent_mig = 0.5) |>
  mutate(absent_cvs = if_else(year == 2022, timing_2022$absent[1], absent_cvs),
         absent_cmb = if_else(year == 2022, timing_2022$absent[2], absent_cmb),
         absent_mig = if_else(year == 2022, timing_2022$absent[3], absent_mig))
stopifnot(all(between(unlist(timing_absent[-1]), 0, 1)))
print(as.data.frame(timing_2022 |> mutate(across(where(is.numeric), \(x) round(x, 3)))))
write_rds(list(timing_2022 = timing_2022, absent = timing_absent,
               crossings = timing_crossings |> select(date, month, n)),
          "data_inter/ukr_timing_absent.rds")
