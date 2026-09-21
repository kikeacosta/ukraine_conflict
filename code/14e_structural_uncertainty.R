# ==============================================================================
# STEP 14e - The intervals with the structural choices folded in
# ==============================================================================
#
# WHY
# ---
# The simulation's intervals cover the inputs that have a distribution (11).
# They do not cover the choices the model holds fixed - how far the register's
# growth is carried, the linkage and resolution rules, the age profile of the
# unverified civilian deaths, the counterfactual's window and model, the
# population base, the timing within 2022, the specification of net migration -
# which 13b-13n show one at a time and which move the results further than the
# drawn inputs do. A reader given the first interval alone takes it for the
# uncertainty of the estimate. This step reports a second one beside it.
#
# METHOD - a mixture over the alternatives
# ----------------------------------------
# Each structural choice is a dimension with a handful of alternatives, the one
# used among them. Every alternative has been projected with every other input
# at its central value, so it has a shift from the central projection: in the
# military total and in the loss by year and sex. In each of the simulation's
# draws one alternative per dimension is sampled, independently across
# dimensions and with equal weights, and the shifts are added to the draw. The
# weights are a convention, not a measurement: equal weights say only that no
# alternative shown is thought less defensible than another. Shifts are taken as
# additive, which holds to first order for shifts this size.
#
# What is left out, and why:
#   - alternatives the simulation already draws over (each release window's
#     multipliers alone, which the drawn weights on the windows reach; the share
#     of the unrecorded prisoners at its bounds; the forecast one SD either way;
#     each migration component across its range);
#   - the counterfactual fitted with the pandemic years, which is a different
#     counterfactual and not an uncertainty about this one;
#   - Donetsk and Luhansk a quarter lower, half lower and taken out, the upper
#     end of a bias and not scenarios for the true count (13e), and the base
#     outside them 2.0 million lower (13m);
#   - no net migration;
#   - residents of occupied Donbas killed in Russian-controlled forces, who are
#     outside the estimand and are reported as a second estimand (13i).
#
# INPUTS   the draws (11, 14) and the sensitivity steps' outputs (13b-13n)
# OUTPUTS  data_inter/ukr_structural_uncertainty.rds:
#            dimensions   every alternative's shift, by dimension
#            intervals    the mean and 95% interval, drawn inputs alone and
#                         with the structural choices, for the military total,
#                         all conflict deaths and the loss by year and sex
#            ranges       each dimension's range of shifts, for the tornado figure
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

r <- function(f) read_rds(file.path("data_inter", f))
central_loss <- r("ukr_migration_decomposition.rds") |> select(year, sex, central = loss)

# 1. THE ALTERNATIVES, AS SHIFTS FROM THE CENTRAL PROJECTION ===================
# one row per dimension, alternative, year and sex: d_loss, and d_military (the
# shift in the four-year military total, the same in every row of an alternative)
shifts <- function(d, dimension, alternative, military = NULL) {
  d |>
    transmute(dimension, alternative = as.character(.data[[alternative]]), year, sex, loss,
              military = if (is.null(military)) NA_real_ else .data[[military]]) |>
    left_join(central_loss, by = c("year", "sex")) |>
    mutate(d_loss = loss - central) |>
    select(dimension, alternative, year, sex, d_loss, military)
}
# The additive and the proportional allocation of a migration draw are compared
# at the two ends of the western bracket, not at the mode, so the shift is the
# proportional end against the additive one and is carried to the central loss.
mig_spec <- r("ukr_migration_specification_e0.rds")
mig_spec <-
  bind_rows(
    mig_spec |> filter(!str_detect(scenario, "^A5")),
    mig_spec |>
      filter(str_detect(scenario, "^A5")) |>
      mutate(end = str_extract(scenario, "crossings end|register end"),
             kind = if_else(str_detect(scenario, "proportional"), "proportional", "additive")) |>
      select(-scenario) |>
      pivot_wider(names_from = kind, values_from = loss) |>
      left_join(central_loss, by = c("year", "sex")) |>
      transmute(year, sex, loss = central + proportional - additive,
                scenario = paste0("A5: ", end, ", proportional"))
  )
lag_tail <- r("ukr_registration_lag_tail.rds")
linkage <- r("ukr_linkage_rules_e0.rds")
drawn_already <- "window's multipliers alone|every unrecorded prisoner|as the returned of every year"
dl <- r("ukr_denominator_donetsk_luhansk.rds")

alternatives <- bind_rows(
  lag_tail$loss |>
    left_join(lag_tail$military |> summarise(military = sum(total), .by = scenario), by = "scenario") |>
    shifts("Registration lag beyond four years", "scenario", "military"),
  linkage |>
    filter(!str_detect(design, drawn_already)) |>
    mutate(design = if_else(str_detect(design, "production"), "as used", design)) |>
    shifts("Linkage rules, resolution model and prisoner-of-war evidence", "design", "military"),
  r("ukr_civilian_age_profile_e0.rds")$loss |>
    shifts("Age profile of the unverified civilian deaths", "scenario"),
  bind_rows(r("ukr_baseline_window_e0.rds"), r("ukr_coherent_forecast_e0.rds")$baseline) |>
    filter(!str_detect(window, "COVID|one SD")) |>
    shifts("Counterfactual: window, jump-off and coherent forecast", "window"),
  r("ukr_covid_carryover_e0.rds")$loss |>
    shifts("Pandemic mortality carried into 2022", "scenario"),
  dl$loss |>
    filter(keep > 0.8) |>
    shifts("Population base of Donetsk and Luhansk", "scenario"),
  r("ukr_base_coherent_e0.rds")$loss |>
    filter(million <= 1, million == 0 | counterfactual != "as fitted") |>
    mutate(scenario = if_else(million == 0, "as used", paste(million, "million fewer aged 20-64"))) |>
    shifts("Population base outside Donetsk and Luhansk", "scenario"),
  r("ukr_timing_2022.rds")$loss |>
    shifts("Timing within 2022", "scenario"),
  mig_spec |>
    shifts("Specification of net migration", "scenario")
)
central_military <- lag_tail$military |> filter(scenario == "to 48 months (used)") |> summarise(m = sum(total)) |> pull(m)
alternatives <-
  alternatives |>
  mutate(d_military = coalesce(military - central_military, 0)) |>
  select(-military)

# every dimension must hold the choice used, with no shift
used <- alternatives |>
  summarise(loss = max(abs(d_loss)), military = max(abs(d_military)), .by = c(dimension, alternative)) |>
  filter(loss < 1e-4, military < 1)
stopifnot(setequal(used$dimension, unique(alternatives$dimension)))

cat("\n=== ALTERNATIVES BY DIMENSION: SHIFT IN THE MILITARY TOTAL AND THE 2025 MALE LOSS ===\n")
print(as.data.frame(
  alternatives |> filter(year == 2025, sex == "m") |>
    transmute(dimension = str_trunc(dimension, 45), alternative = str_trunc(alternative, 60),
              d_military = round(d_military), d_loss = round(d_loss, 3))
))

# 2. ONE ALTERNATIVE PER DIMENSION IN EVERY DRAW ==============================
draws <- r("ukr_e0_loss_by_cause_draws_2022_2025.rds") |> as_tibble()
stopifnot(n_distinct(draws$sim_id) == n_sim)
military_draw <-
  r("ukr_sim_param_draws.rds") |>
  filter(role == "combatants") |>
  summarise(military = sum(draw), .by = sim_id) |>
  arrange(sim_id)
civilians_draw <- draws |> summarise(civilians = sum(dx_civilian), .by = sim_id) |> arrange(sim_id)

set.seed(20260920)
dims <- unique(alternatives$dimension)
picked <- map(set_names(dims), function(dm) {
  alts <- unique(alternatives$alternative[alternatives$dimension == dm])
  sample(alts, n_sim, replace = TRUE)
})
shift_military <- reduce(map(dims, function(dm) {
  a <- alternatives |> filter(dimension == dm) |> distinct(alternative, d_military)
  a$d_military[match(picked[[dm]], a$alternative)]
}), `+`)
shift_loss <- map_dfr(dims, function(dm) {
  tibble(sim_id = seq_len(n_sim), alternative = picked[[dm]]) |>
    left_join(alternatives |> filter(dimension == dm) |> select(alternative, year, sex, d_loss),
              by = "alternative", relationship = "many-to-many")
}) |>
  summarise(d_loss = sum(d_loss), .by = c(sim_id, year, sex))

q <- function(x) tibble(mean = mean(x), lo = quantile(x, 0.025), hi = quantile(x, 0.975))
both <- function(x, shift, what, year = NA, sex = NA) {
  bind_rows(q(x) |> mutate(interval = "drawn inputs"),
            q(x + shift) |> mutate(interval = "drawn inputs and structural choices")) |>
    mutate(what = what, year = year, sex = sex, .before = 1)
}
loss_draws <- draws |> select(sim_id, year, sex, loss = loss_total) |>
  left_join(shift_loss, by = c("sim_id", "year", "sex"))
intervals <- bind_rows(
  both(military_draw$military, shift_military, "Military deaths, 2022-2025"),
  both(military_draw$military + civilians_draw$civilians, shift_military, "Conflict deaths, 2022-2025"),
  loss_draws |>
    nest(.by = c(year, sex)) |>
    mutate(res = pmap(list(data, year, sex), \(d, y, s) both(d$loss, d$d_loss, "Loss of life expectancy", y, s))) |>
    pull(res) |> list_rbind()
)

# 3. EACH DIMENSION'S RANGE, FOR THE TORNADO ==================================
ranges <- bind_rows(
  alternatives |>
    summarise(lo = min(d_loss), hi = max(d_loss), .by = c(dimension, year, sex)) |>
    mutate(what = "Loss of life expectancy"),
  alternatives |>
    distinct(dimension, alternative, d_military) |>
    summarise(lo = min(d_military), hi = max(d_military), .by = dimension) |>
    mutate(what = "Military deaths, 2022-2025"),
  intervals |>
    filter(interval == "drawn inputs") |>
    transmute(dimension = "Drawn inputs: 95% interval of the simulation", year, sex, what,
              lo = lo - mean, hi = hi - mean)
)

write_rds(list(dimensions = alternatives, intervals = intervals, ranges = ranges, seed = 20260920),
          "data_inter/ukr_structural_uncertainty.rds")

cat("\n=== THE INTERVALS, DRAWN INPUTS ALONE AND WITH THE STRUCTURAL CHOICES ===\n")
print(as.data.frame(
  intervals |> filter(is.na(sex) | sex == "m" | year == 2022) |>
    mutate(across(c(mean, lo, hi), \(x) if_else(what == "Loss of life expectancy", round(x, 2), round(x))))
))
message("Done. data_inter/ukr_structural_uncertainty.rds written.")
