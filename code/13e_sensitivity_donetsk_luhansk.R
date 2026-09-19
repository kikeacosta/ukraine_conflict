# ==============================================================================
# STEP 13e - Sensitivity to the Donetsk and Luhansk population base
# ==============================================================================
#
# WHY
# ---
# The projection starts from SSSU's 1 January 2022 population of continental
# Ukraine, which includes Donetsk and Luhansk oblasts: 6.14 of 41.00 million.
# SSSU carried their estimates forward as if registration were complete, but
# only the government-controlled parts registered vital events after 2015, and
# people who left for other oblasts or abroad were largely never deregistered.
# So the base is probably inflated, the rates too low and the loss understated.
# No independent count says by how much. This step measures how much it
# matters: the loss with the two oblasts' 2022 population lower than SSSU's by
# a quarter, by half, and taken out altogether, every other input - conflict
# deaths included - at its mode.
#
# The conflict deaths are not changed. Taking the population out while keeping
# its deaths is the upper end of the bias, not a scenario for the true count.
#
# INPUTS   data_inter/ukr_pop_sssu.rds (01) and the static inputs at the mode
# OUTPUTS  data_inter/ukr_denominator_donetsk_luhansk.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

pop_dl <-
  read_rds("data_inter/ukr_pop_sssu.rds") |>
  filter(reg %in% c("dnk", "luk"), year == 2022) |>
  summarise(pop_dl = sum(pop), .by = c(sex, age))
share_dl <- sum(pop_dl$pop_dl) / sum(mi$pop22_ini$pop)
cat(sprintf("Donetsk and Luhansk: %.2f of %.2f million (%.1f%%)\n",
            sum(pop_dl$pop_dl) / 1e6, sum(mi$pop22_ini$pop) / 1e6, 100 * share_dl))

scenarios <- tibble(
  scenario = c("SSSU (used)", "Donetsk and Luhansk a quarter lower",
               "Donetsk and Luhansk half lower", "Donetsk and Luhansk taken out"),
  keep = c(1, 0.75, 0.5, 0)
)
base_with <- function(keep) {
  mi$pop22_ini |>
    left_join(pop_dl, by = c("sex", "age")) |>
    mutate(pop = pop - (1 - keep) * coalesce(pop_dl, 0)) |>
    select(-pop_dl)
}
loss <-
  scenarios |>
  mutate(res = map(keep, \(k) loss_at_mode(mi$draws, mi$static_inputs, base_with(k)))) |>
  unnest(res)

# the base used must reproduce 13's loss at the mode
stopifnot(isTRUE(all.equal(
  loss |> filter(keep == 1) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(share = share_dl, loss = loss), "data_inter/ukr_denominator_donetsk_luhansk.rds")

cat("\n=== LOSS BY THE POPULATION BASE OF DONETSK AND LUHANSK ===\n")
print(as.data.frame(
  loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
message("Done. data_inter/ukr_denominator_donetsk_luhansk.rds written.")
