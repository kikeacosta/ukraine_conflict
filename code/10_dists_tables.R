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
#     The total is confirmed deaths plus the missing imputed as dead by 09's
#     chain, which is linear in alpha, the share of the never-resolved missing
#     who are alive. The bounds are that total at the ends and the centre of
#     the range of alpha the evidence allows (alpha_evidence(), 00_setup.R):
#     min  = total at alpha_max (the most alive)
#     mode = total at alpha_mode
#     max  = total at alpha_min
#     11 draws alpha itself and computes the total from it; drawing the total
#     from these bounds instead would give the same distribution, because a
#     PERT carried through a linear map is the PERT of the mapped bounds.
#     Two further columns describe the register rather than the estimate:
#     confirmed (its recorded dead, used to split confirmed from imputed) and
#     listed (confirmed plus everyone still listed as missing).
#
#   CIVILIANS (UCDP georeferenced event data)
#     min / mode / max = the low / best / high bounds UCDP publishes,
#     after the deaths of unknown side have been redistributed (step 07)
#
# INPUTS   data_inter/ukr_ucdp_invals.rds                     (from 07_invals)
#          data_inter/ukr_ualosses_..._sex_age_2022_2025.rds  (from 08)
#          data_inter/ukr_ualosses_..._imputed_...rds         (from 09)
#          data_inter/ukr_alpha_missing.rds, ukr_military_alpha_lines.rds (09)
#          data_inter/ukr_migration_bounds.rds                (from 06)
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
  # dead and missing only: prisoners and released prisoners are alive
  filter(status %in% c("dead", "missing"), year %in% 2022:2025) %>%
  summarise(dx = sum(dx), .by = c(status, year, sex, age)) |>
  arrange(status, year, sex, age) |>
  complete(status, year = 2022:2025, sex, age = 0:100, fill = list(dx = 0))

# --- combatants -------------------------------------------------------------
alpha_range <- read_rds("data_inter/ukr_alpha_missing.rds")
alpha_lines <- read_rds("data_inter/ukr_military_alpha_lines.rds")
at_alpha <- function(a) alpha_lines$at_alpha0 - a * alpha_lines$residual

cmb <-
  alpha_lines |>
  transmute(
    year,
    role = "combatants",
    mode = at_alpha(alpha_range$alpha_mode),
    min = at_alpha(alpha_range$alpha_max),
    max = at_alpha(alpha_range$alpha_min),
    confirmed,
    listed = confirmed + missing
  )

# the mode must be the total 09 imputed at the central alpha, and the
# register columns must match the counts the profiles are built from
stopifnot(
  isTRUE(all.equal(
    cmb$mode,
    ual_imp |> summarise(d = sum(dx), .by = year) |> arrange(year) |> pull(d)
  )),
  isTRUE(all.equal(
    cmb$confirmed,
    ual2 |> filter(status == "dead") |> summarise(d = sum(dx), .by = year) |> arrange(year) |> pull(d)
  )),
  isTRUE(all.equal(
    cmb$listed,
    ual2 |> summarise(d = sum(dx), .by = year) |> arrange(year) |> pull(d)
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
migration_bounds <- read_rds("data_inter/ukr_migration_bounds.rds")

mig <-
  migration_bounds |>
  mutate(role = paste0("mig_", component)) |>
  select(year, role, mode, min, max)

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
