# ==============================================================================
# STEP 13f - Sensitivity to the timing of deaths and departures within 2022
# ==============================================================================
#
# WHY
# ---
# The projection removes a share of each year's net outflow and conflict
# deaths before the year's exposure is counted (run_single_sim(), 00_setup.R).
# For 2023-2025 that share is a half, the mid-year convention. 2022 was not
# like that: the war began on 24 February, most departures fell in March and
# April, and civilian deaths peaked in March. So for 2022 the share is taken
# from the months of the events, as 10 measures them: OHCHR's monthly count of
# civilians killed, the register's dead and missing by month of the event, and
# UNHCR's net border crossings by month. This step shows how far that moves
# the loss against the mid-year convention, and against an earlier outflow.
#
# SCENARIOS, every other input at its mode
#   months of the events (used)   10's shares for 2022
#   mid-year convention           a half for everything, as in 2023-2025
#   deaths at their months, net outflow mid-year
#   net outflow at 1 April        the year's net outflow at the turn of March
#                                 and April, earlier than the crossings show
#
# INPUTS   data_inter/ukr_timing_absent.rds (10) and the static inputs at the mode
# OUTPUTS  data_inter/ukr_timing_2022.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()
timing <- read_rds("data_inter/ukr_timing_absent.rds")
cat("\n=== MEAN TIMING WITHIN 2022, AS A FRACTION OF THE YEAR ===\n")
print(as.data.frame(timing$timing_2022 |> mutate(across(where(is.numeric), \(x) round(x, 3)))))

a <- set_names(timing$timing_2022$absent, timing$timing_2022$component)
april <- 1 - as.numeric(as.Date("2022-04-01") - as.Date("2022-01-01")) / 365
with_timing <- function(cvs, cmb, mig) {
  mi$draws |>
    mutate(absent_cvs = if_else(year == 2022, cvs, 0.5),
           absent_cmb = if_else(year == 2022, cmb, 0.5),
           absent_mig = if_else(year == 2022, mig, 0.5))
}
scenarios <- list(
  "months of the events (used)" = mi$draws,
  "mid-year convention" = with_timing(0.5, 0.5, 0.5),
  "deaths at their months, net outflow mid-year" =
    with_timing(a[["civilian deaths"]], a[["military deaths"]], 0.5),
  "net outflow at 1 April" =
    with_timing(a[["civilian deaths"]], a[["military deaths"]], april)
)
loss <- imap_dfr(scenarios, \(d, s) loss_at_mode(d, mi$static_inputs, mi$pop22_ini) |> mutate(scenario = s))

# the timing used must reproduce 13's loss at the mode
stopifnot(isTRUE(all.equal(
  loss |> filter(scenario == "months of the events (used)") |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(timing = timing$timing_2022, loss = loss), "data_inter/ukr_timing_2022.rds")

cat("\n=== LOSS BY THE TIMING OF 2022 EVENTS ===\n")
print(as.data.frame(
  loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
message("Done. data_inter/ukr_timing_2022.rds written.")
