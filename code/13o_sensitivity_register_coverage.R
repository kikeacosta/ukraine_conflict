# ==============================================================================
# STEP 13o - Sensitivity: deaths the register never lists
# ==============================================================================
#
# Every correction in the pipeline measures the register against itself: the
# chain ladder completes it by its own growth (08b), the linkage follows it
# across its own releases and the imputation bounds its missing by its own
# resolutions (09). None can see a combatant who dies and is never listed, and
# no source available to this study measures how many there are. So the
# undercount is not estimated here. It is stated as a size: what the combatant
# total and the loss of life expectancy would be if the register missed 10%,
# 25% or 40% of combatant deaths.
#
# The 25% is the one outside figure there is: US officials put Ukrainian
# military deaths at about 70,000 in August 2023, and the register reaches
# about 52,000 by then even with every one of its missing counted dead, so if
# that figure is right the register misses at least a quarter. 10% and 40%
# bracket it. No source supports any of the three as an estimate.
#
# A register that misses a share k of combatant deaths holds 1 - k of them, so
# the true total is the estimate divided by 1 - k. The deaths it misses are
# taken to be like those it lists: the same years, ages and sexes, and the
# same split between the registered and the missing. Everything else is held at
# its central value. The analysis is deterministic, and it is kept out of the
# second interval (14e): it is an unmeasured bias in one direction, not a choice
# between alternatives.
#
# INPUTS   the static inputs and draws at the central values (mode_projection_inputs())
# OUTPUTS  data_inter/ukr_register_coverage_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()
missed <- c(0, 0.10, 0.25, 0.40)

loss <- map_dfr(missed, function(k) {
  d <- mi$draws |> mutate(draw_cmb = draw_cmb / (1 - k), conf_cmb = conf_cmb / (1 - k))
  loss_at_mode(d, mi$static_inputs, mi$pop22_ini) |> mutate(missed = k)
})
military <- map_dfr(missed, function(k) mi$draws |> transmute(year, missed = k, military = draw_cmb / (1 - k)))

# the register as it is must reproduce 13's loss at the central values
stopifnot(isTRUE(all.equal(
  loss |> filter(missed == 0) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))
stopifnot(all(diff(military |> summarise(t = sum(military), .by = missed) |> pull(t)) > 0))

cat("\n=== COMBATANT DEATHS AND THE LOSS IF THE REGISTER MISSES A SHARE OF THEM ===\n")
print(as.data.frame(
  military |> summarise(combatants = round(sum(military)), .by = missed) |>
    left_join(loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
                select(missed, col, loss) |> pivot_wider(names_from = col, values_from = loss),
              by = "missed")))

write_rds(list(loss = loss, military = military), "data_inter/ukr_register_coverage_e0.rds")
cat("\nDone. data_inter/ukr_register_coverage_e0.rds written.\n")
