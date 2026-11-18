# ==============================================================================
# STEP 13p - How the siege of Mariupol is split between civilians and combatants
# ==============================================================================
#
# WHY
# ---
# UCDP records some deaths without saying whose they were, and 07 spreads them
# over the dyad's civilians and combatants in proportion to the deaths of the
# same dyad and year that UCDP does attribute (22% civilian for the Russia-Ukraine
# dyad in 2022). In 2022 nearly all of these deaths are one event, the siege of
# Mariupol (GED event 438282, 2 March to 7 May 2022): 15,996 deaths of unknown
# status, 94% of the year's total, with a low of 12,066 and a high of 64,183. So
# the front-wide mix decides how many of the siege's dead are civilians, and with
# them the level of the 2022 civilian total. Mariupol's own record differs from
# the front's: of the 11,212 deaths UCDP attributes to a side in the city, 87%
# are civilian.
#
# This step splits the event by four rules, every other input at its central
# value, and shows how the 2022 civilian deaths and the loss of life expectancy
# move. Nothing else changes: the other deaths of unknown status in 2022 keep the
# dyad-year share, the combatants and the civilians of 2023-2025 are as estimated,
# and the civilian deaths added or removed take OHCHR's age-sex profile of 2022
# like every other civilian death.
#
# SCENARIOS
# ---------
#   none civilian      none of the event's deaths is civilian: all are combatants
#   dyad-year share    the share of the dyad and year, as 07 splits it (used)
#   own share          the share among the deaths UCDP attributes to a side in
#                      Mariupol raion in 2022
#   all civilian       all of the event's deaths are civilian
#
# The rules move only the mode of the 2022 civilian total. Its low and high are
# scaled as 07 scales them, the deaths of each dyad and year times that year's
# ratios of low and high to best, and the central value is the mean of the
# Beta-PERT, (min + 4 mode + max) / 6, which is what the projection takes under
# civilian_source "pert"; under uncertainUCDP, that model's mean with the event's
# share changed
# (central_values(), 00_setup.R). The scenario used must reproduce 07's civilian
# totals and 13's loss at the central values.
#
# INPUTS   data_inter/ucdp_ged_events_slim.rds (07), ukr_ucdp_invals.rds (07),
#          ukr_migration_decomposition.rds (13), the static inputs at the mode.
#          The raw GED file, data_input/ucdp/GEDEvent_v26_1.csv, is read when it
#          is present, to check the figures for the event and for the city.
# OUTPUTS  data_inter/ukr_mariupol_split_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# 1. THE EVENT AND THE CITY'S OWN SHARE ========================================
# The cached extract carries neither the event's id nor the second administrative
# level, so the event is found by its deaths and bounds, and the city's attributed
# deaths are stated here and checked against the raw file when it is present.
mariupol_event <- list(id = 438282, year = 2022, unknown = 15996, low = 12066, high = 64183)
# deaths UCDP attributes to a side in Mariupol raion in 2022, in every conflict:
# civilians, and combatants of side A (542) and side B (878)
mariupol_attributed <- list(civilians = 9792, combatants = 542 + 878)
share_own <- mariupol_attributed$civilians /
  (mariupol_attributed$civilians + mariupol_attributed$combatants)

raw_ged <- "data_input/ucdp/GEDEvent_v26_1.csv"
if (file.exists(raw_ged)) {
  ged <-
    read_csv(raw_ged, show_col_types = FALSE,
             col_select = c(id, year, country, adm_2, deaths_a, deaths_b, deaths_civilians,
                            deaths_unknown, best, low, high)) |>
    replace_na(list(deaths_a = 0, deaths_b = 0, deaths_civilians = 0, best = 0, low = 0))
  ged_event <- ged |> filter(id == mariupol_event$id)
  stopifnot(
    nrow(ged_event) == 1,
    ged_event$year == mariupol_event$year,
    ged_event$deaths_unknown == mariupol_event$unknown,
    ged_event$best == mariupol_event$unknown,
    ged_event$low == mariupol_event$low,
    ged_event$high == mariupol_event$high,
    ged_event$deaths_a + ged_event$deaths_b + ged_event$deaths_civilians == 0
  )
  # events with no deaths at their best and low estimates are left out, as in 07
  ged_city <- ged |>
    filter(country == "Ukraine", year == mariupol_event$year, adm_2 %in% "Mariupol raion",
           low + best > 0)
  stopifnot(
    sum(ged_city$deaths_civilians) == mariupol_attributed$civilians,
    sum(ged_city$deaths_a + ged_city$deaths_b) == mariupol_attributed$combatants
  )
  rm(ged, ged_event, ged_city)
} else {
  message("The raw GED file is absent: the figures for the event and the city are not re-checked.")
}

# 2. 07's SPLIT OF THE DEATHS OF UNKNOWN STATUS, FOR UKRAINE'S CIVILIANS ========
events <-
  read_rds("data_inter/ucdp_ged_events_slim.rds") |>
  filter(country == "Ukraine") |>
  # Crimea and Sevastopol are outside the population: the civilians killed there
  # come off the event's totals and bounds, as 07 takes them off
  mutate(
    in_crimea = str_detect(coalesce(adm1, ""), regex("crimea|sevastopol", ignore_case = TRUE)),
    t = if_else(in_crimea, t - c, t),
    t_l = if_else(in_crimea, pmax(0, t_l - c), t_l),
    t_u = if_else(in_crimea, pmax(0, t_u - c), t_u),
    c = if_else(in_crimea, 0, c),
    is_event = year == mariupol_event$year & u == mariupol_event$unknown &
      t_l == mariupol_event$low & t_u == mariupol_event$high
  )
stopifnot(sum(events$is_event) == 1, events$a[events$is_event] == 0, events$b[events$is_event] == 0,
          events$c[events$is_event] == 0, !events$in_crimea[events$is_event])

# the share the event takes in 07: the civilians among the deaths attributed in
# the event's dyad and year, all three of whose categories are non-zero
event_row <- events |> filter(is_event)
dyad_year <-
  events |>
  filter(year == event_row$year, type_of_violence == event_row$type_of_violence,
         side_a == event_row$side_a, side_b == event_row$side_b) |>
  summarise(a = sum(a), b = sum(b), c = sum(c), u = sum(u))
stopifnot(dyad_year$a > 0, dyad_year$b > 0, dyad_year$c > 0)
share_dyad <- dyad_year$c / (dyad_year$a + dyad_year$b + dyad_year$c)
cat(sprintf("\nMariupol: %s of the dyad-year's %s deaths of unknown status; the dyad-year share civilian is %.1f%%, the city's own %.1f%%\n",
            format(mariupol_event$unknown, big.mark = ","), format(dyad_year$u, big.mark = ","),
            100 * share_dyad, 100 * share_own))

# The deaths of unknown status go to civilians in proportion to the attributed
# deaths of the same dyad and year, or, where one of the three categories has none,
# of the dyad over all its years (state dyads) or of the country over all its years
# (everything else), as 07 spreads them. The event's deaths leave that pool and
# enter the civilians at the share given, so the other deaths of unknown status
# keep the dyad-year share. Bounds scale by each row's ratios of low and high to
# best, as in 07.
spread_unknown <- function(d) {
  d |>
    mutate(c2 = case_when(
      u == 0 ~ c,
      a != 0 & b != 0 & c != 0 ~ c + u * c / (a + b + c),
      .default = c + u * sum(c) / (sum(a) + sum(b) + sum(c))
    ) + c_event) |>
    ungroup()
}
sum_cols <- c("a", "b", "c", "u", "c_event", "t", "t_l", "t_u")
civilians_of_ukraine <- function(event_share) {
  d <-
    events |>
    mutate(c_event = if_else(is_event, u * event_share, 0),
           u = if_else(is_event, 0, u),
           state_dyad = type_of_violence == 1 & str_detect(coalesce(side_b, ""), "Government"))
  bind_rows(
    d |> filter(state_dyad) |>
      summarise(across(all_of(sum_cols), sum), .by = c(country, side_a, side_b, year)) |>
      group_by(country, side_a, side_b) |>
      spread_unknown(),
    d |> filter(!state_dyad) |>
      summarise(across(all_of(sum_cols), sum), .by = c(country, year)) |>
      group_by(country) |>
      spread_unknown()
  ) |>
    summarise(dts = sum(c2), dts_l = sum(c2 * t_l / t), dts_u = sum(c2 * t_u / t), .by = year) |>
    filter(year >= 2022) |>
    arrange(year)
}

# 3. THE FOUR RULES ============================================================
rules <- tibble(
  rule = c("none civilian", "dyad-year share", "own share", "all civilian"),
  event_share = c(0, share_dyad, share_own, 1),
  scenario = c(
    "Mariupol's unknown-status deaths: none civilian",
    sprintf("Mariupol's unknown-status deaths: dyad-year share, %.0f%% civilian (used)", 100 * share_dyad),
    sprintf("Mariupol's unknown-status deaths: own attributed share, %.0f%% civilian", 100 * share_own),
    "Mariupol's unknown-status deaths: all civilian"
  )
)
used <- rules$scenario[rules$rule == "dyad-year share"]

by_year <-
  rules |>
  mutate(res = map(event_share, civilians_of_ukraine)) |>
  unnest(res) |>
  mutate(central = central_values(tibble(role = "civilians", min = dts_l, mode = dts, max = dts_u))$mode)
# Under uncertainUCDP (civilian_source, 00_setup.R) the central value is that
# model's exact mean (10). Only the event's own civilian share moves between the
# rules, and the mean is linear in each event's share, so the 2022 mean moves by
# the change in the share times the event's own uncertainUCDP mean.
if (civilian_source == "uncertainUCDP") {
  stopifnot(event_row$t == event_row$u)
  event_mean <- uncertainUCDP::mean_uncertainUCDP(event_row$t, tov = c("sb", "ns", "os")[event_row$type_of_violence])
  uu_central <- read_rds("data_inter/ukr_param_table.rds") |> filter(role == "civilians") |> select(year, uu = mode)
  by_year <- by_year |>
    left_join(uu_central, by = "year") |>
    mutate(central = uu + if_else(year == mariupol_event$year, (event_share - share_dyad) * event_mean, 0)) |>
    select(-uu)
}

# the rule used must give 07's totals, bounds included, in every year
ucdp <- read_rds("data_inter/ukr_ucdp_invals.rds") |> filter(role == "civilians") |> arrange(year)
used_res <- by_year |> filter(scenario == used) |> arrange(year)
stopifnot(
  isTRUE(all.equal(as.numeric(used_res$year), as.numeric(ucdp$year))),
  isTRUE(all.equal(used_res$dts, ucdp$dts, tolerance = 1e-9)),
  isTRUE(all.equal(used_res$dts_l, ucdp$dts_l, tolerance = 1e-9)),
  isTRUE(all.equal(used_res$dts_u, ucdp$dts_u, tolerance = 1e-9))
)
civilians <-
  by_year |>
  filter(year == mariupol_event$year) |>
  transmute(scenario, rule, civilians_min = dts_l, civilians_mode = dts, civilians_max = dts_u,
            civilians_central = central, event_share)

# the years after 2022 do not depend on the rule, and the 2022 total rises with
# the share of the event that is civilian
stopifnot(
  by_year |> filter(year > mariupol_event$year) |>
    summarise(spread = max(central) - min(central), .by = year) |> pull(spread) < 1e-6,
  all(diff(civilians$civilians_mode[order(civilians$event_share)]) > 0),
  all(diff(civilians$civilians_central[order(civilians$event_share)]) > 0)
)

# 4. THE LOSS OF LIFE EXPECTANCY AT THE CENTRAL VALUES =========================
mi <- mode_projection_inputs()
stopifnot(isTRUE(all.equal(
  mi$draws |> arrange(year) |> pull(draw_cvs),
  used_res |> pull(central),
  tolerance = 1e-9
)))

loss_with <- function(central_2022) {
  draws <- mi$draws |> mutate(draw_cvs = if_else(year == mariupol_event$year, central_2022, draw_cvs))
  loss_at_mode(draws, mi$static_inputs, mi$pop22_ini) |> select(year, sex, loss)
}
loss <-
  civilians |>
  mutate(res = map(civilians_central, loss_with)) |>
  select(scenario, res) |>
  unnest(res)

# the rule used must reproduce 13's loss at the central values
stopifnot(isTRUE(all.equal(
  loss |> filter(scenario == used) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

# 5. THE FOUR-YEAR CIVILIAN TOTAL ==============================================
conflict_deaths <-
  by_year |>
  summarise(civilians_central_total = sum(central), .by = scenario) |>
  mutate(change = civilians_central_total - civilians_central_total[scenario == used])

write_rds(list(civilians = civilians |> select(scenario, rule, starts_with("civilians_"), event_share),
               loss = loss, conflict_deaths = conflict_deaths),
          "data_inter/ukr_mariupol_split_e0.rds")

cat("\n=== 2022 CIVILIAN DEATHS BY THE RULE THAT SPLITS THE SIEGE OF MARIUPOL ===\n")
print(as.data.frame(
  civilians |>
    mutate(across(starts_with("civilians_"), round), event_share = round(event_share, 3)) |>
    select(-scenario)
))
cat("\n=== FOUR-YEAR CIVILIAN TOTAL AT THE CENTRAL VALUES, AGAINST THE RULE USED ===\n")
print(as.data.frame(conflict_deaths |> mutate(across(where(is.numeric), round))))
cat("\n=== LOSS OF LIFE EXPECTANCY BY THE RULE THAT SPLITS THE SIEGE OF MARIUPOL ===\n")
print(as.data.frame(
  loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
message("Done. data_inter/ukr_mariupol_split_e0.rds written.")
