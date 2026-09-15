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
#     min  = individually confirmed deaths only
#     mode = confirmed + the share of the missing imputed as dead (step 09)
#     max  = confirmed + ALL of the missing treated as dead
#
#   CIVILIANS (UCDP georeferenced event data)
#     min / mode / max = the low / best / high bounds UCDP publishes,
#     after the deaths of unknown side have been redistributed (step 07)
#
# INPUTS   data_inter/ukr_ucdp_invals.rds                     (from 07_invals)
#          data_inter/ukr_ualosses_..._sex_age_2022_2025.rds  (from 08)
#          data_inter/ukr_ualosses_..._imputed_...rds         (from 09)
# OUTPUT   data_inter/ukr_param_table.rds
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
  filter(status != "prisoner", year %in% 2022:2025) %>%
  summarise(dx = sum(dx), .by = c(status, year, sex, age)) |>
  arrange(status, year, sex, age) |>
  complete(status, year = 2022:2025, sex, age = 0:100, fill = list(dx = 0))

# --- combatants -------------------------------------------------------------
cmb_m <- ual_imp |> summarise(mode = sum(dx), .by = c(year)) # confirmed + imputed
cmb_l <- ual2 |> filter(status == "dead") |> summarise(min = sum(dx), .by = c(year))
cmb_u <- ual2 |> summarise(max = sum(dx), .by = c(year)) # dead + all missing

cmb <-
  cmb_m |>
  left_join(cmb_l, by = "year") |>
  left_join(cmb_u, by = "year") |>
  mutate(role = "combatants") |>
  select(year, role, mode, min, max)

# --- civilians --------------------------------------------------------------
cvs <-
  ucdp |>
  filter(role == "civilians") |>
  select(year, role, mode = dts, min = dts_l, max = dts_u)

# --- migration --------------------------------------------------------------
# NET EMIGRATION by year (in person-years, cumulative stock change)
# Based on reconciliation of Eurostat (official TP registrations),
# UNHCR global totals, CES refugee survey (5.6M Jan 2026 mode),
# and Pozniak/SBGSU border crossing data sourced from daily Facebook posts.
#
# KEY ASSUMPTIONS:
#   1. CES age-sex distribution (Jan 2026 snapshot) applied to all annual
#      Eurostat flows 2022-2025: assumes stable demographic profile of
#      refugee outflows year-to-year (defensible but not proven; selective
#      returns or multi-wave migrations could alter the structure).
#   2. Pozniak 2022-2023 from published paper; 2024-2026 extrapolated from
#      SBGSU Facebook-sourced tallies — source/method needs confirmation.
#
# Mode from CES estimate of refugee population in EU+EFTA (5.6M)
# scaled by Eurostat annual flow distribution
# Min/Max bounds based on uncertainty in definition and coverage:
#   floor ~2.2-2.5M (lower-bound estimate accounting for returns/measurement)
#   ceiling ~6.2-6.5M (upper-bound estimate accounting for unregistered)
mig <-
  tribble(
    ~year, ~min, ~mode, ~max,
    2022,  1100000, 1900000, 2300000,
    2023,  1200000, 1950000, 2400000,
    2024,  -200000, 1200000, 1800000,
    2025,  -500000, 600000, 1200000
  ) |>
  mutate(role = "migration") |>
  select(year, role, mode, min, max)

param_table <- bind_rows(cmb, cvs, mig)

# rpert() requires min <= mode <= max; catch any ordering problem here rather
# than as an obscure error inside the simulation loop
stopifnot(
  nrow(param_table) == 12,
  all(param_table$min <= param_table$mode),
  all(param_table$mode <= param_table$max)
)

print(as.data.frame(param_table))

write_rds(param_table, "data_inter/ukr_param_table.rds")

# totals implied by the table, for reference
param_table |>
  summarise(across(c(min, mode, max), sum), .by = role)
