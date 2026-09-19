# ==============================================================================
# STEP 13h - Sensitivity to the age profile of civilian deaths
# ==============================================================================
#
# WHY
# ---
# Civilian deaths take their level from UCDP and their age-sex profile from
# OHCHR (07). OHCHR's ages come from the deaths it verified, mostly in
# government-controlled territory, and it knew the age of only 5,228 of the
# 9,511 civilians it verified as killed in 2022-2023. The deaths UCDP counts
# beyond OHCHR's verified total - most of 2022's, few after - are largely those
# OHCHR could not reach: in occupied territory and in besieged cities such as
# Mariupol. The people who stayed there were older than the population as a
# whole, as in the frontline villages where OHCHR finds people aged 60 or more
# to be over 45% of the civilians killed in 2025 and half in March 2026
# (data_input/official_figures.csv), against about a third in its 2022
# profile. An older death costs fewer years of life expectancy, so the profile
# used may overstate the civilian loss.
#
# This step gives the deaths beyond OHCHR's verified count an older profile,
# every other input at its mode:
#   half over 60     OHCHR's profile of the same year, re-weighted so that people
#                    aged 60 or more are the share OHCHR last reported among
#                    civilians killed near the front line (half), the shape
#                    within the two age bands and the sex ratio at each age kept
#   as pre-war       the counterfactual's own deaths of the year by age and sex,
#   deaths           as if the war had raised everyone's mortality in
#                    proportion: the oldest profile a count of deaths could take
# The deaths OHCHR verified keep its profile, and each year's civilian total is
# unchanged, so only the ages move.
#
# INPUTS   data_inter/ukr_ohchr_annual_totals.rds (07), the static inputs at
#          the mode, data_input/official_figures.csv
# OUTPUTS  data_inter/ukr_civilian_age_profile_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

# the civilian deaths at the mode, and the part beyond OHCHR's verified count
beyond <-
  mi$draws |>
  select(year, civilians = draw_cvs) |>
  left_join(read_rds("data_inter/ukr_ohchr_annual_totals.rds") |> select(year, verified = dts),
            by = "year") |>
  mutate(beyond = civilians - verified, share_beyond = beyond / civilians)
stopifnot(all(beyond$beyond >= 0))

# the share of people aged 60 or more among the civilians killed that OHCHR
# last reported near the front line
older_share <-
  read_csv("data_input/official_figures.csv", show_col_types = FALSE) |>
  filter(key == "older_share_killed_frontline") |>
  slice_max(as.Date(date), n = 1) |>
  pull(value)
stopifnot(length(older_share) == 1, older_share > 0, older_share < 1)
# the scenario is reported as "half over 60" (15, the paper): a later figure
# must change those labels too
stopifnot(older_share == 0.5)

cx <- mi$static_inputs |> select(year, sex, age, prop_cvs)

# OHCHR's profile re-weighted to the older share, within each year
half_over_60 <-
  cx |>
  mutate(old = age >= 60) |>
  mutate(band_share = prop_cvs / sum(prop_cvs), .by = c(year, old)) |>
  transmute(year, sex, age, older = if_else(old, older_share, 1 - older_share) * band_share)

# the counterfactual's own deaths at the mode, by year, sex and age
mode_run <- run_single_sim(1, mi$draws, mi$static_inputs, mi$pop22_ini)
as_prewar <-
  mode_run |>
  mutate(older = expected / sum(expected), .by = year) |>
  select(year, sex, age, older)

# each year's profile: the verified deaths on OHCHR's profile, the rest on the
# older one
profile_with <- function(older) {
  cx |>
    left_join(older, by = c("year", "sex", "age")) |>
    left_join(beyond |> select(year, share_beyond), by = "year") |>
    transmute(year, sex, age, prop_cvs = (1 - share_beyond) * prop_cvs + share_beyond * older)
}
scenarios <- list(
  "OHCHR's profile (used)" = cx,
  "Beyond OHCHR's verified count: half over 60" = profile_with(half_over_60),
  "Beyond OHCHR's verified count: as pre-war deaths" = profile_with(as_prewar)
)

# every profile keeps each year's civilian total, and the older ones are older
# in every year
for (p in scenarios) {
  stopifnot(all(abs(tapply(p$prop_cvs, p$year, sum) - 1) < 1e-9), all(p$prop_cvs >= 0))
}
stats <- imap_dfr(scenarios, \(p, nm) {
  p |>
    summarise(share_60_plus = sum(prop_cvs[age >= 60]),
              mean_age = sum((age + 0.5) * prop_cvs),
              share_female = sum(prop_cvs[sex == "f"]),
              .by = year) |>
    mutate(scenario = nm, .before = 1)
})
stopifnot(stats |>
            mutate(used = share_60_plus[scenario == names(scenarios)[1]], .by = year) |>
            filter(scenario != names(scenarios)[1]) |>
            with(all(share_60_plus > used)))

loss <- imap_dfr(scenarios, \(p, nm) {
  si <- mi$static_inputs |> select(-prop_cvs) |> left_join(p, by = c("year", "sex", "age"))
  loss_at_mode(mi$draws, si, mi$pop22_ini) |> mutate(scenario = nm, .before = 1)
})

# the profile used must reproduce 13's loss at the mode
stopifnot(isTRUE(all.equal(
  loss |> filter(scenario == names(scenarios)[1]) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(beyond = beyond, older_share = older_share, stats = stats, loss = loss),
          "data_inter/ukr_civilian_age_profile_e0.rds")

cat("\n=== CIVILIAN DEATHS BEYOND OHCHR'S VERIFIED COUNT, AT THE MODE ===\n")
print(as.data.frame(beyond |> mutate(across(c(civilians, beyond), round), share_beyond = round(share_beyond, 3))))
cat("\n=== THE CIVILIAN PROFILES ===\n")
print(as.data.frame(stats |> mutate(across(c(share_60_plus, share_female), \(x) round(x, 3)),
                                    mean_age = round(mean_age, 1))))
cat("\n=== LOSS BY THE AGE PROFILE OF CIVILIAN DEATHS ===\n")
print(as.data.frame(
  loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
message("Done. data_inter/ukr_civilian_age_profile_e0.rds written.")
