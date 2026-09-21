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
# 0.5 and 1.0 million, by a quarter, by half, and taken out altogether, every
# other input - conflict deaths included - at its mode.
#
# The conflict deaths are not changed. Taking the population out while keeping
# its deaths is the upper end of the bias, not a scenario for the true count.
#
# The rest of the base has an overcount of its own: people registered in
# Ukraine but living abroad before 2022, mostly labour migrants of working
# age. The 2019 electronic census counted 37.29 million on the territory the
# government controlled; SSSU carries more. So the base outside the two
# oblasts is also made 0.5 and 1.0 million smaller at ages 20-64, men and
# women in proportion to their numbers at those ages.
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

# the two oblasts 0.5 and 1.0 million lower, a quarter and a half lower, and
# taken out altogether
dl_total <- sum(pop_dl$pop_dl)
scenarios <- tibble(
  scenario = c("SSSU (used)", "Donetsk and Luhansk 0.5 million lower",
               "Donetsk and Luhansk 1.0 million lower", "Donetsk and Luhansk a quarter lower",
               "Donetsk and Luhansk half lower", "Donetsk and Luhansk taken out"),
  keep = c(1, 1 - 0.5e6 / dl_total, 1 - 1e6 / dl_total, 0.75, 0.5, 0)
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

# labour migrants: the base outside Donetsk and Luhansk smaller at ages 20-64
working_rest <-
  mi$pop22_ini |>
  left_join(pop_dl, by = c("sex", "age")) |>
  mutate(rest = pop - coalesce(pop_dl, 0), working = age >= 20 & age <= 64)
# The absentees are cut from ages 20-64, either in proportion to the population
# there (male_share NULL) or with a given share of the cut taken from men.
# Pre-war labour migration from Ukraine was predominantly male, and the conflict
# deaths fall on men, so the split matters for the male rates; no source here
# measures it, and the two readings are shown as a range rather than an estimate.
base_without <- function(million, male_share = NULL) {
  w <- working_rest
  if (is.null(male_share)) {
    w <- w |> mutate(cut = if_else(working, million * 1e6 * rest / sum(rest[working]), 0))
  } else {
    by_sex <- c(m = male_share, f = 1 - male_share) * million * 1e6
    w <- w |> mutate(cut = if_else(working, by_sex[sex] * rest / sum(rest[working]), 0), .by = sex)
  }
  stopifnot(all(w$cut <= w$rest))
  w |> mutate(pop = pop - cut) |> select(all_of(names(mi$pop22_ini)))
}
migrants <-
  bind_rows(
    tibble(million = c(0.5, 1.0), male_share = NA_real_,
           scenario = paste0("Outside Donetsk and Luhansk, ", format(c(0.5, 1.0), nsmall = 1),
                             " million fewer aged 20-64")),
    tibble(million = c(0.5, 1.0), male_share = 2 / 3,
           scenario = paste0("Outside Donetsk and Luhansk, ", format(c(0.5, 1.0), nsmall = 1),
                             " million fewer aged 20-64, two thirds of them men"))
  ) |>
  mutate(res = map2(million, male_share,
                    function(m, s) loss_at_mode(mi$draws, mi$static_inputs,
                                                base_without(m, if (is.na(s)) NULL else s)))) |>
  unnest(res) |>
  select(-million, -male_share)

# the base used must reproduce 13's loss at the mode
stopifnot(isTRUE(all.equal(
  loss |> filter(keep == 1) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(share = share_dl, loss = loss, migrants = migrants),
          "data_inter/ukr_denominator_donetsk_luhansk.rds")

cat("\n=== LOSS BY THE POPULATION BASE OF DONETSK AND LUHANSK ===\n")
print(as.data.frame(
  loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
cat("\n=== LOSS WITH THE BASE OUTSIDE DONETSK AND LUHANSK SMALLER AT WORKING AGES ===\n")
print(as.data.frame(
  migrants |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
message("Done. data_inter/ukr_denominator_donetsk_luhansk.rds written.")
