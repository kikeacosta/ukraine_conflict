# ==============================================================================
# STEP 13s - The estimate as a function of the share of the missing counted dead
# ==============================================================================
#
# WHY
# ---
# About half of the combatant total rests on the register's missing, and no
# source says how many of them are dead. Rather than one number, this step
# gives the whole line: the combatant total and the loss of life expectancy at
# every share of the missing counted dead, from none to all, with the scenarios
# the paper discusses marked on it:
#
#   named deaths only     0%: the register's dead, completed for lag; the
#                         missing all alive
#   floor                 the named dead plus the bodies not yet identified
#                         (data_input/official_figures.csv, unidentified_bodies),
#                         nine in ten assumed combatants listed as missing,
#                         placed across event years by months since the event (09)
#   central               the register's resolved cases applied to the
#                         unresolved: of the missing whose fate a later release
#                         settled, the share recorded dead rather than taken
#                         prisoner or gone from the list beyond list maintenance
#                         (09: composition and captivity_by_window)
#   heavy desertion       20% of the missing alive: about 6% of the 311,327
#                         absence and desertion cases listed as missing and
#                         never resurfacing, the one channel that could hold
#                         many living people the register cannot see
#   upper                 the rule of 09: no evidence of life means dead, the
#                         prisoners from the official counts and half the
#                         register's bound alive outside captivity
#   all dead              100%
#
# The lines are straight: the total is linear in the share, and the loss nearly
# so. Everything else - the registered dead and their late registrations, the
# age-sex profiles, civilians, migration, the counterfactual - is held at its
# central value, as in 13 and 13b.
#
# INPUTS   data_inter/ukr_military_inputs.rds, ukr_ualosses_linkage_checks.rds,
#          ukr_ualosses_imputation_table.rds (09); data_input/official_figures.csv;
#          the static inputs at the mode (mode_projection_inputs())
# OUTPUTS  data_inter/ukr_missing_dead_share_sweep.rds, ukr_missing_dead_share_bars.rds,
#          ukr_missing_dead_share_models.rds. No figure: step 17 draws Figure A3 from them
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")
library(ggrepel)
library(patchwork)

mi <- mode_projection_inputs()
draw_cvs_by_year <- mi$draws |> select(year, draw_cvs)
draw_mig_by_year <- mi$draws |> select(year, draw_mig)
mil <- read_rds("data_inter/ukr_military_inputs.rds")
lc <- read_rds("data_inter/ukr_ualosses_linkage_checks.rds")
imp <- read_rds("data_inter/ukr_ualosses_imputation_table.rds")

# 1. THE MISSING COUNTED DEAD AT A SHARE p ====================================
# captives 0 and alive_other 1 - p: a share p of every year's missing dead
at_dead <- function(p) military_draws(mil, 0, 1 - p) |> transmute(year, draw_cmb = military, conf_cmb = confirmed)
base <- at_dead(1)
missing_total <- sum(base$draw_cmb - base$conf_cmb)
named_total <- sum(base$conf_cmb)

# 2. THE SCENARIOS =============================================================
# the scenarios' shares are 09's (mil$scenarios; 00_setup.R, the three
# scenarios of the missing): the two bounds, none and all of the missing dead;
# the resolution-based (central) from the register's resolved cases; the appendix's bodies and rule scenarios
sc <- mil$scenarios
bodies_unidentified <- sc$bodies_unidentified
p_central <- sc$p_central
scenarios <- tibble(
  scenario = c("bound_low", "bodies", "central", "heavy desertion", "rule", "bound_high"),
  label = c("Lower bound: none of the missing dead",
            "Bodies: the named dead and the\nbodies awaiting identification",
            "Resolution-based: the register's resolved\ncases applied to the unresolved",
            "Heavy desertion: 20% of the missing alive",
            "Rule: no evidence of life means dead",
            "Upper bound: all of the missing dead"),
  share_dead = c(0,
                 sc$p_bodies,
                 p_central,
                 0.80,
                 sum(imp$imputed_dead) / sum(imp$missing_stock),
                 1)
)
stopifnot(all(diff(scenarios$share_dead) > 0))
# the line applies one share to every event year. Two scenarios split their dead
# unevenly over the years (09): the upper, the rule's own table, because the
# prisoners and the bound fall unevenly; the floor, because its bodies are placed
# by months since the event. Their points take that split, at the same total, so
# that the figure agrees with Tables 5, 6 and A19 (they sit off the line,
# the floor's by up to a few tenths of a year). The central applies one share.
own_split <- c("bodies", "rule")
own_mil <- map(set_names(own_split), \(s) scenario_military(mil, s) |> transmute(year, draw_cmb = military, conf_cmb = confirmed))
own_p <- map(own_mil, \(m) (m$draw_cmb - m$conf_cmb) / (base$draw_cmb - base$conf_cmb))
for (s in own_split) {
  stopifnot(identical(own_mil[[s]]$year, base$year), all(own_p[[s]] > 0), all(own_p[[s]] < 1),
            isTRUE(all.equal(sum(own_mil[[s]]$draw_cmb - own_mil[[s]]$conf_cmb) / missing_total,
                             scenarios$share_dead[scenarios$scenario == s], tolerance = 1e-5)))
}
stopifnot(all(own_p$rule > 0.85))
upper_p <- own_p$rule
mil_at <- function(p, scenario) if (scenario %in% own_split) own_mil[[scenario]] else at_dead(p)

# 3. THE LINE, AND THE LOSS AT EVERY POINT ====================================
# the loss in each year, and for the four years as one period: deaths and
# person-years summed over 2022-2025 by sex and age, one life table on the sum
loss_from_sim <- function(sim) {
  stopifnot(min(sim$pop) >= 0)
  pooled <-
    sim |>
    summarise(across(c(expected, civilian, combatant_confirmed, combatant_imputed, pop), sum), .by = c(sex, age)) |>
    mutate(year = "2022–2025")
  e0_loss <- function(d) {
    d <- d |> mutate(mx_all = (expected + civilian + combatant_confirmed + combatant_imputed) / pop, mx_bsn = expected / pop)
    e0 <- function(col) {
      d |> select(year, sex, age, mx = all_of(col)) |> group_by(year, sex) |> do(lifetable(dt_in = .data)) |>
        ungroup() |> filter(age == 0) |> select(year, sex, ex)
    }
    e0("mx_bsn") |> rename(e0_bsn = ex) |> left_join(e0("mx_all") |> rename(e0_all = ex), by = c("year", "sex")) |>
      transmute(year, sex, loss = e0_bsn - e0_all)
  }
  bind_rows(e0_loss(sim |> mutate(year = as.character(year))), e0_loss(pooled))
}
loss_for <- function(mil_by_year) {
  draws_df <-
    draw_cvs_by_year |>
    left_join(draw_mig_by_year, by = "year") |>
    left_join(mil_by_year, by = "year") |>
    mutate(sim_id = 1)
  loss_from_sim(run_single_sim(1, draws_df, mi$static_inputs, mi$pop22_ini))
}
grid <- tibble(share_dead = seq(0, 1, by = 0.05)) |> anti_join(scenarios, by = "share_dead")
sweep <-
  bind_rows(scenarios, grid) |>
  arrange(share_dead) |>
  mutate(mil = map2(share_dead, scenario, mil_at),
         military = map_dbl(mil, \(d) sum(d$draw_cmb)),
         loss = map(mil, loss_for),
         # the line itself applies one share to every year at every point; only the
         # upper scenario's point leaves it
         off_line = scenario %in% own_split,
         mil_line = pmap(list(mil, share_dead, off_line), \(m, p, o) if (o) at_dead(p) else m),
         loss_line = pmap(list(loss, mil_line, off_line), \(l, m, o) if (o) loss_for(m) else l))
stopifnot(isTRUE(all.equal(sweep$military, named_total + sweep$share_dead * missing_total, tolerance = 1e-5)))

# combatant deaths by year along the line, and the four years together: the
# marked points (`*_long`, the scenarios' own values) and the line (`line_*`)
mil_long_of <- function(d, col) {
  bind_rows(
    d |> select(scenario, share_dead, mil = all_of(col)) |> unnest(mil) |>
      transmute(scenario, share_dead, year = as.character(year), military = draw_cmb),
    d |> select(scenario, share_dead, mil = all_of(col)) |>
      transmute(scenario, share_dead, year = "2022–2025", military = map_dbl(mil, \(x) sum(x$draw_cmb)))
  )
}
military_long <- mil_long_of(sweep, "mil")
line_military <- mil_long_of(sweep, "mil_line")
loss_long <- sweep |> select(scenario, label, share_dead, military, loss) |> unnest(loss)
line_loss <- sweep |> select(scenario, share_dead, loss_line) |> unnest(loss_line)
sweep <- sweep |> select(-mil, -loss, -mil_line, -loss_line, -off_line)

# the crossing: the share of the missing counted dead at which 2025 overtakes
# 2022 as the deadliest year for men (the lines are nearly straight, so the
# root is found by linear interpolation on the grid and then projected)
gap <- loss_long |> filter(sex == "m", year %in% c("2022", "2025")) |>
  pivot_wider(id_cols = share_dead, names_from = year, values_from = loss) |>
  arrange(share_dead) |> mutate(d = `2025` - `2022`)
k <- max(which(gap$d < 0))
share_cross <- gap$share_dead[k] + (0 - gap$d[k]) / (gap$d[k + 1] - gap$d[k]) * (gap$share_dead[k + 1] - gap$share_dead[k])
cross <- tibble(scenario = "crossing", label = "Crossing: 2025 overtakes 2022 as the deadliest year for men",
                share_dead = share_cross, mil = list(at_dead(share_cross))) |>
  mutate(military = map_dbl(mil, \(d) sum(d$draw_cmb)), loss = map(mil, loss_for))
sweep <- bind_rows(sweep, cross |> select(scenario, label, share_dead, military)) |> arrange(share_dead)
loss_long <- bind_rows(loss_long, cross |> select(scenario, label, share_dead, military, loss) |> unnest(loss)) |> arrange(share_dead)
military_long <- bind_rows(military_long, mil_long_of(cross, "mil")) |> arrange(share_dead)
# the crossing lies on the line
line_military <- bind_rows(line_military, mil_long_of(cross, "mil")) |> arrange(share_dead)
line_loss <- bind_rows(line_loss, cross |> select(scenario, share_dead, loss) |> unnest(loss)) |> arrange(share_dead)
cat(sprintf("\n2025 overtakes 2022 as the deadliest year for men at %.1f%% of the missing counted dead\n", 100 * share_cross))
write_rds(list(sweep = sweep, loss = loss_long, military = military_long, scenarios = scenarios,
               line_loss = line_loss, line_military = line_military,
               missing_total = missing_total, named_total = named_total,
               bodies_unidentified = bodies_unidentified, bodies_combatant_share = bodies_combatant_share),
          "data_inter/ukr_missing_dead_share_sweep.rds")

cat("\n=== THE SCENARIOS: SHARE OF THE MISSING COUNTED DEAD, COMBATANT DEATHS AND THE MALE LOSS ===\n")
print(as.data.frame(
  loss_long |>
    filter(!is.na(scenario), sex == "m", year %in% c("2022", "2025", "2022–2025")) |>
    pivot_wider(id_cols = c(scenario, share_dead, military), names_from = year, values_from = loss, names_prefix = "male_") |>
    mutate(share_dead = round(100 * share_dead, 1), military = round(military), across(starts_with("male"), \(x) round(x, 2)))
))

# 3b. THE DRAWN INTERVAL AT EACH POINT ==========================================
# At each marked share the missing are fixed and everything else is drawn as
# in 11: each draw's civilians, net migration, Lee-Carter path and lag factors,
# the last through its registered dead and missing, which 11 saves by draw
# (roles confirmed and missing). The interval therefore covers the drawn
# inputs with the share of the missing fixed, and not the evidence on that
# share.
period <- "2022–2025"
pd <- read_rds("data_inter/ukr_sim_param_draws.rds")
draws_wide <-
  pd |>
  filter(role %in% c("confirmed", "missing", "civilians", "lc_f", "lc_m", "mig_west", "mig_ru")) |>
  pivot_wider(id_cols = c(sim_id, year), names_from = role, values_from = draw, values_fill = 0) |>
  arrange(sim_id, year) |>
  mutate(draw_cvs = civilians, draw_mig = mig_west + mig_ru, dk_f = lc_f, dk_m = lc_m, conf_cmb = confirmed) |>
  select(sim_id, year, draw_cvs, draw_mig, dk_f, dk_m, conf_cmb, missing)
stopifnot(all(draws_wide$missing > 0), n_distinct(draws_wide$sim_id) >= 100)
# the recovered missing average close to the point values of 09
cat("\n=== THE MISSING BY YEAR: MEAN OVER THE DRAWS, AND AT THE POINT FACTORS ===\n")
print(as.data.frame(draws_wide |> summarise(missing_draws = round(mean(missing)), .by = year) |>
                      left_join(base |> transmute(year, missing_point = round(draw_cmb - conf_cmb)), by = "year")))
lc_error <- read_rds("data_inter/ukr_lc_forecast_error.rds")
static_full <-
  mi$static_inputs |>
  left_join(lc_error$bx, by = c("sex", "age")) |>
  left_join(lc_error$variance |> select(sex, year, vk), by = c("sex", "year"))
stopifnot(!anyNA(static_full), all(c("bx", "vk") %in% names(static_full)))

bar_points <- bind_rows(scenarios, cross |> select(scenario, label, share_dead)) |>
  # the share by event year: one share, except the upper scenario's own split
  mutate(p_year = map2(share_dead, scenario, \(p, s) if (s %in% own_split) own_p[[s]] else rep(p, length(upper_p))))
sims <- draws_wide |> group_split(sim_id)
stopifnot(all(map_lgl(sims, \(d) all(as.numeric(d$year) == as.numeric(base$year)))))
# The projections below (one per point and draw) take almost all of this step's
# time, so their result is kept with a fingerprint of what it was computed from,
# as 11 keeps the draws: they are run again only when the draws, the points, the
# static inputs, the starting population or the projection's code have changed.
bars_file <- "data_inter/ukr_missing_dead_share_bars.rds"
bars_inputs_file <- sub("[.]rds$", "_inputs.rds", bars_file)
num_sum <- function(...) vapply(list(...), function(d) {
  sum(vapply(Filter(is.numeric, d), function(v) sum(v * seq_along(v), na.rm = TRUE), 0))
}, 0)
code_txt <- utf8ToInt(paste(unlist(lapply(list(run_single_sim, loss_from_sim), deparse)), collapse = "\n"))
bars_inputs_now <- c(num_sum(draws_wide, static_full, mi$pop22_ini, list(unlist(bar_points$p_year))),
                     srb, sum(code_txt * seq_along(code_txt)))
bars_fresh <- file.exists(bars_file) && file.exists(bars_inputs_file) &&
  isTRUE(all.equal(read_rds(bars_inputs_file), bars_inputs_now, tolerance = 1e-12))
if (bars_fresh) {
  message("cache hit  : ", bars_file, " (inputs unchanged; the projections are not run again)")
  bars_loss <- read_rds(bars_file)$loss
  bars_mil <- read_rds(bars_file)$military
  stopifnot(read_rds(bars_file)$n_draws == length(sims))
} else {
  message("cache miss : the drawn interval at ", nrow(bar_points), " points, ", length(sims), " draws each ...")
  bars_raw <- map_dfr(seq_len(nrow(bar_points)), function(k) {
    p <- bar_points$share_dead[k]
    pv <- bar_points$p_year[[k]]
    map_dfr(sims, function(d) {
      d <- d |> mutate(draw_cmb = conf_cmb + pv * missing)
      sim <- run_single_sim(d$sim_id[1], d, static_full, mi$pop22_ini)
      bind_rows(
        # both sexes: Figure 3 of the paper carries women too
        loss_from_sim(sim) |> transmute(sim_id = d$sim_id[1], year, sex, loss),
        d |> transmute(sim_id, year = as.character(year), military = draw_cmb),
        tibble(sim_id = d$sim_id[1], year = period, military = sum(d$draw_cmb))
      )
    }) |> mutate(scenario = bar_points$scenario[k], share_dead = p)
  })
  q <- function(x, p) quantile(x, p, names = FALSE)
  bars_loss <- bars_raw |> filter(!is.na(loss)) |> summarise(lo = q(loss, 0.025), hi = q(loss, 0.975), .by = c(scenario, share_dead, year, sex))
  bars_mil <- bars_raw |> filter(!is.na(military)) |> summarise(lo = q(military, 0.025), hi = q(military, 0.975), .by = c(scenario, share_dead, year))
  write_rds(list(loss = bars_loss, military = bars_mil, n_draws = length(sims)), bars_file)
  write_rds(bars_inputs_now, bars_inputs_file)
}
stopifnot(setequal(bars_loss$sex, c("f", "m")))
# the male loss, for the log
bars_loss_m <- bars_loss |> filter(sex == "m") |> select(-sex)
cat("\n=== THE 95% INTERVAL OF THE DRAWN INPUTS AT EACH POINT, MALE LOSS ===\n")
print(as.data.frame(bars_loss_m |> filter(year %in% c("2022", "2025", period)) |>
                      mutate(share_dead = round(100 * share_dead), across(c(lo, hi), \(x) round(x, 2))) |> arrange(share_dead, year)))

# 4. THE POINTS NAMED IN THE SUPPLEMENT ========================================
# This step saves data only: the figure of the sweep (Figure A3 of the paper) is
# drawn by step 17 from the files written above.

# the supplement: every point named in the text, and the projected models
models <-
  lc$designs |>
  filter(str_detect(design, "duration model|event-year cohorts")) |>
  transmute(design, military, share_dead = (military - named_total) / missing_total)
# every point as a letter, the key in a note under the figure, so that more
# points can be added without crowding the strip
scen_sm <-
  bind_rows(scenarios, cross |> select(scenario, label, share_dead)) |>
  arrange(share_dead) |>
  mutate(letter = LETTERS[row_number()])
models <- models |> mutate(letter = LETTERS[nrow(scen_sm) + 1])
m_at <- function(sc, yr) loss_long |> filter(scenario == sc, sex == "m", year == yr) |> pull(loss)
note_rows <- scen_sm |>
  mutate(military = sweep$military[match(scenario, sweep$scenario)],
         line = sprintf("%s  %s: %.0f%% dead; %s combatant deaths; male loss %.2f years in 2022, %.2f in 2025, %.2f for 2022–2025",
                        letter, str_replace_all(label, "\n", " "), 100 * share_dead, scales::comma(round(military, -2)),
                        map_dbl(scenario, m_at, "2022"), map_dbl(scenario, m_at, "2025"), map_dbl(scenario, m_at, period)))
note_sm <- paste(str_wrap(c(
  note_rows$line,
  sprintf("%s  Resolution models projected in the rule's place: the duration model with one set of hazards and with the events of 2022 at their own rate, and the chain of event-year cohorts; %s dead",
          models$letter[1], paste0(round(100 * range(models$share_dead)), collapse = "–") |> paste0("%")),
  "The line applies one share to every event year; the bodies and rule scenarios' points take their own split by event year, at the same total, so they sit off the line: the bodies are placed by months since the event, and the rule's prisoners and bound fall unevenly.",
  sprintf("Lines: every other input at its central value. Bars: the 95%% interval of the drawn inputs (the counterfactual's path, civilians, net migration and the lag factors) with the share of the missing fixed, from %s draws; they do not cover the register's evidence on the share itself. The loss of male life expectancy at birth by year, and for 2022–2025 as one period.", scales::comma(length(sims)))
), width = 150, exdent = 3), collapse = "\n")
write_rds(list(models = models, scenarios_sm = scen_sm, note = note_sm, share_cross = share_cross),
          "data_inter/ukr_missing_dead_share_models.rds")
# which year is the deadliest for men at each share: the lines of 2024 and 2025
# cross 2022's near the middle of the axis
deadliest <-
  loss_long |>
  filter(sex == "m", year != period) |>
  summarise(deadliest = year[which.max(loss)], .by = c(share_dead, scenario)) |>
  arrange(share_dead)
cat("\n=== THE DEADLIEST YEAR FOR MEN, BY THE SHARE OF THE MISSING COUNTED DEAD ===\n")
print(as.data.frame(deadliest |> mutate(share_dead = round(100 * share_dead, 1))))
message("Done. data_inter/ukr_missing_dead_share_sweep.rds, _bars.rds and _models.rds written.")
