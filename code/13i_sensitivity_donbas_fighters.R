# ==============================================================================
# STEP 13i - Residents of occupied Donbas killed in Russian-controlled forces
# ==============================================================================
#
# WHY
# ---
# The population base is SSSU's for continental Ukraine, which carries the
# whole of Donetsk and Luhansk oblasts, the parts occupied since 2014
# included. Their residents who died fighting in the units of the
# self-proclaimed republics are in that denominator but in no numerator: the
# civilians killed there are among UCDP's civilian deaths, and the register
# counts Ukraine's own forces only. UCDP does not separate them either: after
# 24 February 2022 it codes the war between the two states alone, so any such
# deaths it records are Russia's, and its sources for Russia's side exclude
# those units (07).
#
# No source counts them after September 2023. BBC Russian and Mediazona put the
# deaths in the two republics' units at 21,000-23,500 by late September 2023
# (data_input/official_figures.csv). Two scenarios bracket what that implies,
# every other input at its mode:
#   counted       the 21,000 and 23,500 counted by September 2023, none after:
#                 the floor
#   continued     the midpoint to September 2023, then deaths at the same
#                 average monthly count to December 2025
# Within each period the deaths follow the months of Ukraine's own military
# deaths (09), a measure of the intensity of the fighting, and they take the
# age-sex profile of the register's dead of the same year: men mobilised at
# the same ages. They are added to the numerator only; the denominator already
# holds them.
#
# INPUTS   data_input/official_figures.csv, data_inter/ukr_military_by_month.rds
#          (09), the static inputs at the mode
# OUTPUTS  data_inter/ukr_donbas_fighters_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

counted <-
  read_csv("data_input/official_figures.csv", show_col_types = FALSE) |>
  filter(key == "dpr_lpr_killed_estimate")
stopifnot(nrow(counted) == 1)
counted_by <- floor_date(as.Date(counted$date), "month")

months <-
  read_rds("data_inter/ukr_military_by_month.rds") |>
  select(month, year, intensity = military) |>
  mutate(counted_period = month <= counted_by)
n_after <- sum(!months$counted_period)

# the deaths by month for a count by September 2023 and a monthly count after
by_month <- function(level, per_month_after) {
  months |>
    mutate(deaths = if_else(counted_period, level * intensity / sum(intensity[counted_period]), 0),
           deaths = deaths + if_else(!counted_period, per_month_after * intensity / mean(intensity[!counted_period]), 0))
}
mid <- (counted$low + counted$high) / 2
scenarios <- list(
  "Counted by September 2023: 21,000" = by_month(counted$low, 0),
  "Counted by September 2023: 23,500" = by_month(counted$high, 0),
  "Continued to 2025 at the same monthly count" = by_month(mid, mid / sum(months$counted_period))
)
deaths <- imap_dfr(scenarios, \(d, nm) d |> summarise(deaths = sum(deaths), .by = year) |> mutate(scenario = nm))
stopifnot(all(deaths$deaths >= 0))

loss_with <- function(extra) {
  d <- mi$draws |>
    left_join(extra, by = "year") |>
    mutate(deaths = coalesce(deaths, 0),
           draw_cmb = draw_cmb + deaths,
           conf_cmb = conf_cmb + deaths) |>
    select(-deaths)
  loss_at_mode(d, mi$static_inputs, mi$pop22_ini)
}
loss <-
  bind_rows(
    loss_at_mode(mi$draws, mi$static_inputs, mi$pop22_ini) |> mutate(scenario = "None (as estimated)"),
    imap_dfr(scenarios, \(d, nm) loss_with(d |> summarise(deaths = sum(deaths), .by = year)) |> mutate(scenario = nm))
  ) |>
  relocate(scenario)

# with none added, the loss must be 13's at the mode
stopifnot(isTRUE(all.equal(
  loss |> filter(scenario == "None (as estimated)") |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(deaths = deaths, loss = loss), "data_inter/ukr_donbas_fighters_e0.rds")

cat("\n=== DEATHS ADDED BY YEAR ===\n")
print(as.data.frame(deaths |> mutate(deaths = round(deaths)) |> pivot_wider(names_from = year, values_from = deaths)))
cat("\n=== LOSS WITH THE RESIDENTS OF OCCUPIED DONBAS KILLED IN RUSSIAN-CONTROLLED FORCES ===\n")
print(as.data.frame(
  loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
message("Done. data_inter/ukr_donbas_fighters_e0.rds written.")
