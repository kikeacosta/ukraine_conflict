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
# Net emigration enters as TWO components, because two different things are
# unknown about it and they are unknown to different degrees:
#
#   mig_west   how much of the western outflow the sources see. The floor is
#              CES's net border crossings, the ceiling the rescaled Eurostat
#              register, and the gap between them is the discrepancy Pozniak
#              (2023) set out - a register that keeps people after they
#              return against a crossing balance that missed the peak weeks.
#   mig_ru_by  displacement into Russia and Belarus, which no register sees
#              and which UNHCR stopped being able to estimate in May 2025.
#
# The bounds, and the sources behind each number, are in migration_bounds in
# 00_setup.R; 06 scales its age-sex profiles to the modes from the same
# table. Both components are systematic rather than year-by-year noise, and
# 11 draws them accordingly.
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

# 06 scaled its profiles to the modes, so the draw at the mode must leave the
# age-sex cells untouched. If these drift apart the projection silently
# rescales every year (see the guard in run_single_sim).
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
