# ==============================================================================
# STEP 13q - Three figures reported beside the headline estimate
# ==============================================================================
#
# WHY
# ---
# The headline loss rests on two choices a reader may weigh differently: the
# missing are counted dead in the central scenario's share (the register's
# resolved cases applied to the unresolved; 09), and the loss is taken
# over continental Ukraine although the conflict deaths and the counterfactual
# do not all refer to that territory. This step gives three figures that let a
# reader see what each choice carries, every one in the deterministic central
# projection (mode_projection_inputs(), loss_at_mode(), 00_setup.R):
#
#   (a) the loss without the missing counted dead: combatant deaths equal to the
#       confirmed deaths alone - the registered dead with their late
#       registrations - so the missing who are counted dead stay in the
#       population. This is the convention of studies that count identified
#       deaths and leave the missing out. The share of the loss the missing carry
#       is 1 - (the loss without them) / (the loss used).
#   (b) the scale of the share on the missing: the change in the four-year
#       combatant total and in the loss per 10 percentage points of the missing
#       alive, read from 13s's sweep over the share counted dead, between 60%
#       and 70% dead, around the central scenario. No projection is run for it.
#   (c) a territory-consistent second estimate, Ukraine without Crimea,
#       Donetsk and Luhansk, in which numerator, denominator and counterfactual
#       refer to one territory. The headline keeps SSSU's population of
#       Donetsk and Luhansk in the base, counts the civilians killed there and
#       the combatants who come from there, and takes its counterfactual from a
#       pre-war series that already leaves out Crimea, Donetsk and Luhansk
#       (04b). Only 4% of the register's dead and missing come from the two
#       oblasts, against 15% of the base.
#
# METHOD
# ------
#   (a) draw_cmb = conf_cmb in the central draws; everything else as used.
#   (b) the rows of 13s's sweep at 60% and 70% of the missing counted dead
#       (40% and 30% alive), and at 65% to check the sweep is linear between
#       them: the combatant total is linear in the share by construction and
#       the loss nearly so.
#   (c) population  SSSU's 1 January 2022 population of Donetsk and Luhansk
#                   taken out of the base by sex and age, as 13e's "taken out"
#                   scenario takes it out, and the net outflow - the flows by
#                   sex and age and each year's total - in proportion to the
#                   population kept
#       combatants  by event year, the confirmed deaths scaled by the share of
#                   the register's dead, and the imputed dead (draw_cmb -
#                   conf_cmb) by the share of its missing, whose place of origin
#                   is outside the two oblasts (08). Records of unknown origin
#                   are spread pro rata: the shares are taken among records
#                   with a known oblast. Kyiv city, which the register gives as
#                   "Kyiv" alone, is outside the two oblasts.
#       civilians   by year, scaled by the share of UCDP's civilian deaths
#                   located outside the two oblasts, attributed as 07 attributes
#                   them: each event's civilian deaths take its dyad-year's
#                   share of the deaths of unknown side. Events with no
#                   first-level unit are spread pro rata.
#       rates       the counterfactual rates, fertility and the age profiles of
#                   the deaths are those of the headline
#       The changes are also made one after another, each projected, to show
#       how far each moves the loss (`steps`).
#
# CAVEATS (c)
# -----------
#   - The register does not say whether its place of origin is the place of
#     birth or of residence.
#   - Internal displacement. SSSU's population of the two oblasts includes
#     people who left for other oblasts after 2014 or in 2022 and were never
#     deregistered (13e); they are taken out of the base, while their deaths
#     stay in the numerator if the register or UCDP places them outside the
#     two oblasts. Combatants born in the two oblasts and living elsewhere are
#     taken out of the numerator if the field gives the place of birth.
#   - Records from the two oblasts are listed as missing 1.7 times as often as
#     dead, against 1.0 in the register as a whole. The imputed dead are split
#     by the origin of the missing, which assumes the missing of every origin
#     are dead in the same proportion.
#   - The net outflow is taken out in proportion to the population: no source
#     here splits the flows by oblast. The second of the steps shows what
#     this does to the loss.
#
# INPUTS   the static inputs and draws at the central values (mode_projection_inputs())
#          data_inter/ukr_missing_dead_share_sweep.rds (13s)
#          data_inter/ualosses_origin_2022_2025.rds (08)
#          data_inter/ucdp_ged_events_slim.rds (07)
#          data_inter/ukr_pop_sssu.rds (01)
#          data_inter/ukr_denominator_donetsk_luhansk.rds (13e), for a check
# OUTPUTS  data_inter/ukr_headline_companions.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

# 0. THE CENTRAL PROJECTION ===================================================
used <- loss_at_mode(mi$draws, mi$static_inputs, mi$pop22_ini)
# it must reproduce the central scenario's loss of 13s, which 13d and the tables report
sweep13s <- read_rds("data_inter/ukr_missing_dead_share_sweep.rds")
stopifnot(isTRUE(all.equal(
  used |> arrange(year, sex) |> pull(loss),
  sweep13s$loss |> filter(scenario %in% "central", year %in% as.character(2022:2025)) |>
    mutate(year = as.numeric(year)) |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))
loss_used <- used |> select(year, sex, loss_used = loss)

# 1. (a) THE LOSS WITHOUT THE IMPUTED MISSING ==================================
without_imputed <-
  loss_at_mode(mi$draws |> mutate(draw_cmb = conf_cmb), mi$static_inputs, mi$pop22_ini) |>
  select(year, sex, loss) |>
  left_join(loss_used, by = c("year", "sex")) |>
  mutate(share_imputed = 1 - loss / loss_used)
stopifnot(all(without_imputed$loss <= without_imputed$loss_used + 1e-9))

# 2. (b) THE CHANGE PER 10 POINTS OF THE MISSING ALIVE =========================
# 13s sweeps the share of the missing counted dead; "per 10 points of the missing
# alive" is read from 70% dead (30% alive) to 60% dead (40% alive), the grid
# points around the central scenario, with 65% as the check of linearity
sweep_mil <- sweep13s$sweep |> transmute(share_dead, share_alive = 1 - share_dead, combatants = military)
sweep_e0 <-
  sweep13s$loss |>
  filter(year %in% as.character(2022:2025)) |>
  transmute(share_dead, share_alive = 1 - share_dead, year = as.numeric(year), sex, loss)
at_share <- function(d, s) d |> filter(abs(share_alive - s) < 1e-9) |> mutate(share_alive = s)
shares_alive <- c(0.30, 0.35, 0.40)
mil3 <- map_dfr(shares_alive, \(s) at_share(sweep_mil, s))
e03 <- map_dfr(shares_alive, \(s) at_share(sweep_e0, s))
stopifnot(nrow(mil3) == 3, nrow(e03) == 3 * 8)
# the sweep is linear between 30% and 40% alive: the combatant total at 35% is the
# midpoint of the two, and the loss within 1% of the change from it (the women's
# loss of 2023 hardly moves, hence the floor of 1e-6 years)
off_mid <- function(v) v[2] - (v[1] + v[3]) / 2
stopifnot(abs(off_mid(mil3$combatants)) < 1e-6 * abs(diff(mil3$combatants[c(1, 3)])))
e0_lin <-
  e03 |>
  arrange(year, sex, share_alive) |>
  summarise(at_30 = loss[1], at_40 = loss[3], off_mid = off_mid(loss), .by = c(year, sex))
stopifnot(all(abs(e0_lin$off_mid) < 0.01 * abs(e0_lin$at_40 - e0_lin$at_30) + 1e-6))

per_10_points <-
  bind_rows(
    tibble(measure = "combatants, 2022-2025", year = NA_real_, sex = NA_character_,
           at_30 = mil3$combatants[1], at_40 = mil3$combatants[3]),
    e0_lin |> transmute(measure = "loss", year, sex, at_30, at_40)
  ) |>
  mutate(per_10_points = at_40 - at_30)

# 3. (c) UKRAINE WITHOUT CRIMEA, DONETSK AND LUHANSK ===========================
# the combatants: the share of the dead and of the missing from outside the two
# oblasts, among records with a known oblast
origin <- read_rds("data_inter/ualosses_origin_2022_2025.rds")
in_dl <- c("Donetsk", "Luhansk")
cmb_shares <-
  origin |>
  filter(origin != "unknown") |>
  summarise(outside = 1 - sum(n[origin %in% in_dl]) / sum(n), .by = c(year, status)) |>
  pivot_wider(names_from = status, values_from = outside, names_prefix = "cmb_")
# the missing per dead among the two oblasts' records and in the whole register
missing_per_dead <- function(d) sum(d$n[d$status == "missing"]) / sum(d$n[d$status == "dead"])
ratio_dl <- missing_per_dead(origin |> filter(origin %in% in_dl))
ratio_all <- missing_per_dead(origin)

# the civilians: UCDP's civilian deaths in Ukraine as 07 counts them. Crimea and
# Sevastopol are left out, and the deaths of unknown side are spread within each
# dyad-year (inter-state) or country-year (the rest) in proportion to the
# attributed deaths, across the group's years where one category is empty.
crimea <- regex("crimea|sevastopol", ignore_case = TRUE)
ged <-
  read_rds("data_inter/ucdp_ged_events_slim.rds") |>
  mutate(c = if_else(country == "Ukraine" & str_detect(coalesce(adm1, ""), crimea), 0, c),
         group = if_else(type_of_violence == 1 & str_detect(side_b, "Government"),
                         paste(side_a, side_b, sep = " | "), "intra-state"))
# The siege of Mariupol is taken apart, as in 07: its deaths of unknown side leave
# the pool and enter at the mix of the deaths attributed in Mariupol
# (mariupol_siege, 00_setup.R); c_siege is its civilians, all in Donetsk oblast.
by_group <-
  ged |>
  mutate(siege = is_mariupol_siege(ged),
         c_siege = if_else(siege, u * mariupol_siege$share[["c"]], 0),
         u = if_else(siege, 0, u)) |>
  summarise(across(c(a, b, c, u, c_siege), sum), .by = c(country, group, year)) |>
  mutate(c2 = case_when(u == 0 ~ c,
                        a != 0 & b != 0 & c != 0 ~ c + u * c / (a + b + c),
                        .default = c + u * sum(c) / (sum(a) + sum(b) + sum(c))) + c_siege,
         .by = c(country, group)) |>
  filter(country == "Ukraine", year %in% 2022:2025)
siege_by_year <- by_group |> summarise(c_siege = sum(c_siege), .by = year)
stopifnot(sum(siege_by_year$c_siege > 0) == 1)
# the attribution must reproduce UCDP's best civilian counts as 10 records them (its
# ucdp_best column; the mode is uncertainUCDP's mean when that is the source). The
# share outside the two oblasts is UCDP best's share, applied to the civilians as used
stopifnot(isTRUE(all.equal(
  by_group |> summarise(c2 = sum(c2), .by = year) |> arrange(year) |> pull(c2),
  mi$param_table |> filter(role == "civilians") |> arrange(year) |> pull(ucdp_best)
)))
cvs_shares <-
  ged |>
  filter(country == "Ukraine", year %in% 2022:2025, !is.na(adm1)) |>
  # the siege's civilians are not spread over the other events: they are added in
  # Donetsk oblast, where the siege was
  left_join(by_group |> filter(c > 0) |> transmute(group, year, f = (c2 - c_siege) / c), by = c("group", "year")) |>
  filter(c > 0) |>
  summarise(in_dl = sum(c * f * (adm1 %in% c("Donetsk oblast", "Luhansk oblast"))), all = sum(c * f),
            .by = year) |>
  left_join(siege_by_year, by = "year") |>
  transmute(year, cvs = 1 - (in_dl + c_siege) / (all + c_siege))
stopifnot(!any(is.na(cvs_shares$cvs)))
shares <- cmb_shares |> left_join(cvs_shares, by = "year") |> arrange(year)
stopifnot(nrow(shares) == 4, !any(is.na(shares)))

# the population: the two oblasts taken out of the base as 13e takes them out
pop_dl <-
  read_rds("data_inter/ukr_pop_sssu.rds") |>
  filter(reg %in% c("dnk", "luk"), year == 2022) |>
  summarise(pop_dl = sum(pop), .by = c(sex, age))
base_out <-
  mi$pop22_ini |>
  left_join(pop_dl, by = c("sex", "age")) |>
  mutate(pop = pop - coalesce(pop_dl, 0)) |>
  select(-pop_dl)
kept <- sum(base_out$pop) / sum(mi$pop22_ini$pop)
# the net outflow in proportion to the population kept: the flows at the mode
# and each year's drawn total
static_out <- mi$static_inputs |> mutate(ems = ems * kept)

draws_out <- function(outflow = TRUE, civilians = TRUE, combatants = TRUE) {
  mi$draws |>
    left_join(shares, by = "year") |>
    mutate(draw_mig = if (outflow) draw_mig * kept else draw_mig,
           draw_cvs = if (civilians) draw_cvs * cvs else draw_cvs,
           conf_out = if (combatants) conf_cmb * cmb_dead else conf_cmb,
           draw_cmb = if (combatants) conf_out + (draw_cmb - conf_cmb) * cmb_missing else draw_cmb,
           conf_cmb = conf_out) |>
    select(all_of(names(mi$draws)))
}
steps <-
  tribble(
    ~step, ~outflow, ~civilians, ~combatants,
    "Donetsk and Luhansk taken out of the base (13e)", FALSE, FALSE, FALSE,
    "and the net outflow in proportion", TRUE, FALSE, FALSE,
    "and the civilians killed outside the two oblasts", TRUE, TRUE, FALSE,
    "and the combatants from outside them (the estimate)", TRUE, TRUE, TRUE
  ) |>
  mutate(res = pmap(list(outflow, civilians, combatants), function(o, cv, cm) {
    loss_at_mode(draws_out(o, cv, cm), if (o) static_out else mi$static_inputs, base_out)
  })) |>
  select(step, res) |>
  unnest(res) |>
  left_join(loss_used, by = c("year", "sex"))
# the first step is 13e's "taken out" scenario
stopifnot(isTRUE(all.equal(
  steps |> filter(str_detect(step, "13e")) |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_denominator_donetsk_luhansk.rds")$loss |> filter(keep == 0) |>
    arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

territory_loss <-
  steps |>
  filter(str_detect(step, "the estimate")) |>
  select(year, sex, loss, loss_used, e0_bsn, e0_all)
territory_deaths <-
  draws_out() |>
  select(year, confirmed = conf_cmb, combatants = draw_cmb, civilians = draw_cvs) |>
  left_join(mi$draws |> select(year, combatants_used = draw_cmb, civilians_used = draw_cvs),
            by = "year")

write_rds(
  list(without_imputed = without_imputed,
       per_10_points = per_10_points,
       territory = list(shares = shares, deaths = territory_deaths, loss = territory_loss,
                        steps = steps)),
  "data_inter/ukr_headline_companions.rds"
)

# 4. PRINTED ===================================================================
cat("\n=== (a) THE LOSS WITHOUT THE IMPUTED MISSING ===\n")
print(as.data.frame(
  without_imputed |>
    mutate(across(c(loss, loss_used), \(x) round(x, 2)), share_imputed = round(100 * share_imputed, 1)) |>
    arrange(sex, year)
))

cat("\n=== (b) CHANGE PER 10 POINTS OF THE MISSING ALIVE (13s's sweep, 30% to 40% alive) ===\n")
print(as.data.frame(
  per_10_points |> mutate(across(c(at_30, at_40, per_10_points), \(x) round(x, 3)))
))
cat(sprintf("largest departure of the loss at 35%% from the line: %.5f years\n",
            max(abs(e0_lin$off_mid))))

cat("\n=== (c) UKRAINE WITHOUT CRIMEA, DONETSK AND LUHANSK ===\n")
cat("register's dead and missing, 2022-2025, by origin:\n")
print(as.data.frame(origin |> summarise(n = sum(n), .by = c(status, origin)) |>
                      pivot_wider(names_from = status, values_from = n) |>
                      mutate(share = round(100 * (dead + missing) / sum(dead + missing), 2))))
cat(sprintf("missing per dead: Donetsk and Luhansk %.2f, the register %.2f\n", ratio_dl, ratio_all))
cat(sprintf("population kept: %.2f of %.2f million (%.1f%%)\n",
            sum(base_out$pop) / 1e6, sum(mi$pop22_ini$pop) / 1e6, 100 * kept))
cat("shares outside the two oblasts (%):\n")
print(as.data.frame(shares |> mutate(across(-year, \(x) round(100 * x, 2)))))
cat("deaths kept:\n")
print(as.data.frame(territory_deaths |> mutate(across(-year, round))))
cat("loss, step by step:\n")
print(as.data.frame(
  steps |> mutate(col = paste0(sex, "_", year), loss = round(loss, 2)) |>
    select(step, col, loss) |> pivot_wider(names_from = col, values_from = loss) |>
    bind_rows(used |> mutate(step = "as used (continental Ukraine)", col = paste0(sex, "_", year),
                             loss = round(loss, 2)) |>
                select(step, col, loss) |> pivot_wider(names_from = col, values_from = loss))
))
message("Done. data_inter/ukr_headline_companions.rds written.")
