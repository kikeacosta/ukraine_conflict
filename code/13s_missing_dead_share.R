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
# OUTPUTS  data_inter/ukr_missing_dead_share_sweep.rds, ukr_missing_dead_share_models.rds
#          figures/exploratory/missing_dead_share_main.png (the three scenarios)
#          figures/exploratory/missing_dead_share_sm.png (every point named in the text)
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
message("  the drawn interval at ", nrow(bar_points), " points, ", length(sims), " draws each ...")
bars_raw <- map_dfr(seq_len(nrow(bar_points)), function(k) {
  p <- bar_points$share_dead[k]
  pv <- bar_points$p_year[[k]]
  map_dfr(sims, function(d) {
    d <- d |> mutate(draw_cmb = conf_cmb + pv * missing)
    sim <- run_single_sim(d$sim_id[1], d, static_full, mi$pop22_ini)
    bind_rows(
      # both sexes: step 15's figure of the loss by scenario carries women too
      loss_from_sim(sim) |> transmute(sim_id = d$sim_id[1], year, sex, loss),
      d |> transmute(sim_id, year = as.character(year), military = draw_cmb),
      tibble(sim_id = d$sim_id[1], year = period, military = sum(d$draw_cmb))
    )
  }) |> mutate(scenario = bar_points$scenario[k], share_dead = p)
})
q <- function(x, p) quantile(x, p, names = FALSE)
bars_loss <- bars_raw |> filter(!is.na(loss)) |> summarise(lo = q(loss, 0.025), hi = q(loss, 0.975), .by = c(scenario, share_dead, year, sex))
bars_mil <- bars_raw |> filter(!is.na(military)) |> summarise(lo = q(military, 0.025), hi = q(military, 0.975), .by = c(scenario, share_dead, year))
stopifnot(setequal(bars_loss$sex, c("f", "m")))
write_rds(list(loss = bars_loss, military = bars_mil, n_draws = length(sims)), "data_inter/ukr_missing_dead_share_bars.rds")
# the figure below is of the male loss
bars_loss_m <- bars_loss |> filter(sex == "m") |> select(-sex)
cat("\n=== THE 95% INTERVAL OF THE DRAWN INPUTS AT EACH POINT, MALE LOSS ===\n")
print(as.data.frame(bars_loss_m |> filter(year %in% c("2022", "2025", period)) |>
                      mutate(share_dead = round(100 * share_dead), across(c(lo, hi), \(x) round(x, 2))) |> arrange(share_dead, year)))

# 4. THE FIGURE ================================================================
# Four panels on one x-axis. The scenarios are positions on the axis, not
# points of one year, so they are named once in a strip at the top, each with
# a vertical guide running down through the panels below: the loss of male life
# expectancy in each of the four years and for 2022-2025 as one period, the
# combatant deaths by year, and the four-year total. The years take a
# sequential colour scale, the four-year period one colour, thicker.
pts <- loss_long |> filter(!is.na(scenario))
year_levels <- c("2022", "2023", "2024", "2025", "2022–2025")
period <- "2022–2025"
# the lines at one share for every year; the points at the scenarios' own values
line_m <- line_loss |> filter(sex == "m") |> mutate(year = factor(year, year_levels))
pts_m <- pts |> filter(sex == "m") |> mutate(year = factor(year, year_levels))
pal <- c("2022" = "#9ecae1", "2023" = "#4292c6", "2024" = "#08519c", "2025" = "#08306b", "2022–2025" = "#e63946")
# the years thin and partly transparent, so that the four do not crowd the
# period's line, which keeps its full weight
widths <- c("2022" = 0.45, "2023" = 0.45, "2024" = 0.45, "2025" = 0.45, "2022–2025" = 1.6)
sizes <- c("2022" = 1.2, "2023" = 1.2, "2024" = 1.2, "2025" = 1.2, "2022–2025" = 2.8)
alphas <- c("2022" = 0.6, "2023" = 0.6, "2024" = 0.6, "2025" = 0.6, "2022–2025" = 1)
xscale <- scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1), limits = c(-0.01, 1.01), expand = c(0, 0))
# one x-axis for the stacked panels: the lower one keeps its labels
no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
guides <- geom_vline(data = scenarios, aes(xintercept = share_dead), linetype = "dotted", colour = "grey55")
year_scales <- list(
  scale_colour_manual(values = pal, name = "Year", drop = FALSE),
  scale_linewidth_manual(values = widths, name = "Year", drop = FALSE),
  scale_size_manual(values = sizes, name = "Year", drop = FALSE),
  scale_alpha_manual(values = alphas, name = "Year", drop = FALSE)
)
base_theme <- theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.margin = margin(0, 4, 2, 4))

# One figure in two versions: the main text's, with the three scenarios the
# paper discusses, and the supplement's, with every point of the axis that the
# text names - the register's named deaths alone, the heavy-desertion reading,
# all of the missing dead, and the resolution models projected in the rule's
# place (09, Table A19), which fall almost on the upper scenario.
# scen: the scenarios to mark, with the row of the label strip each takes and
# the side its label reads to (hjust 0 right of the guide, 1 left of it);
# models: NULL, or the projected models' shares, drawn as a grey bracket.
# letters: TRUE names each point by a letter in the strip, with the key in a
# note under the figure (`note`), so that many points fit; FALSE writes the
# labels out, for a few points.
main_scenarios <- c("bound_low", "central", "bound_high")
# period_labels: TRUE writes the 2022-2025 values beside the period's points at
# the scenarios, the loss in the top panel and the four-year total in the bottom
# one (the main text's version; the supplement's note gives them instead).
sweep_figure <- function(scen, models = NULL, file, height, letters = FALSE, note = NULL, period_labels = FALSE) {
  # the guides: the three scenarios of the main text thicker than the other points
  scen <- scen |> mutate(main = scenario %in% main_scenarios)
  guides <- list(
    geom_vline(data = scen |> filter(!main), aes(xintercept = share_dead), linetype = "dotted", colour = "grey55", linewidth = 0.4),
    geom_vline(data = scen |> filter(main), aes(xintercept = share_dead), linetype = "dotted", colour = "grey30", linewidth = 1)
  )
  if (letters) {
    # the letter with its share beneath; two points less than 6 points apart
    # take turns at two heights
    # labels at the edges read inwards; the projected models' bracket sits
    # below the letters, its own letter to its left
    strip_labels <- scen |>
      arrange(share_dead) |>
      mutate(tight = coalesce(share_dead - lag(share_dead) < 0.06, FALSE),
             row = if_else(tight, 2.6, 1.5),
             hjust = case_when(share_dead < 0.03 ~ 0, share_dead > 0.97 ~ 1, .default = 0.5),
             x_text = share_dead + case_when(hjust == 0 ~ -0.008, hjust == 1 ~ 0.008, .default = 0),
             letter = sprintf("%s\n(%.0f%%)", letter, 100 * share_dead))
    rows <- 1
    p_strip <-
      ggplot(strip_labels, aes(share_dead, row)) +
      geom_label(aes(x = x_text, label = letter, size = main, hjust = hjust), colour = "grey10", fontface = "bold",
                 label.size = 0, fill = "white", label.padding = unit(0.1, "lines"), lineheight = 0.85) +
      scale_size_manual(values = c(`TRUE` = 3.8, `FALSE` = 3), guide = "none") +
      xscale +
      scale_y_continuous(limits = c(0, 3.4), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.margin = margin(2, 4, 0, 4))
    if (!is.null(models)) {
      p_strip <- p_strip +
        annotate("segment", x = min(models$share_dead), xend = max(models$share_dead), y = 0.4, yend = 0.4, colour = "grey45") +
        annotate("point", x = models$share_dead, y = 0.4, shape = 21, fill = "white", colour = "grey45", size = 1.6) +
        annotate("text", x = min(models$share_dead) - 0.012, y = 0.4, hjust = 1, size = 3, colour = "grey45",
                 fontface = "bold", label = models$letter[1])
    }
  } else {
    strip_labels <-
      scen |>
      mutate(label = str_replace_all(label, "\n", " "),
             label = str_wrap(sprintf("%s (%.0f%%)", label, 100 * share_dead), 34),
             x_text = share_dead + if_else(hjust == 0, 0.008, -0.008))
    rows <- max(strip_labels$row)
    p_strip <-
      ggplot(strip_labels, aes(share_dead, row)) +
      geom_segment(aes(xend = share_dead, y = 0, yend = row - 0.4), colour = "grey55", linetype = "dotted") +
      geom_text(aes(x = x_text, label = label, hjust = hjust), size = 2.9, colour = "grey15",
                lineheight = 0.9, vjust = 0.5) +
      xscale +
      scale_y_continuous(limits = c(0, rows + 0.6), expand = c(0, 0)) +
      theme_void() +
      theme(plot.margin = margin(2, 4, 0, 4))
    if (!is.null(models)) {
      # the projected models as a bracket just above the axis of the strip
      p_strip <- p_strip +
        annotate("segment", x = min(models$share_dead), xend = max(models$share_dead), y = 0.25, yend = 0.25, colour = "grey45") +
        annotate("point", x = models$share_dead, y = 0.25, shape = 21, fill = "white", colour = "grey45", size = 2) +
        annotate("text", x = min(models$share_dead) - 0.012, y = 0.25, hjust = 1, size = 2.7, colour = "grey45", fontface = "italic",
                 label = sprintf("Resolution models projected in the rule's place (%s)",
                                 paste0(round(100 * range(models$share_dead)), collapse = "–") |> paste0("%")))
    }
  }
  # the lines, the points, and at each marked share the 95% interval of the
  # drawn inputs as a bar, the years set a little apart so that the bars can
  # be told apart
  bar_offset <- c("2022" = -0.014, "2023" = -0.007, "2024" = 0, "2025" = 0.007, "2022–2025" = 0)
  lines_of <- function(d, d_pts, y, ylab, b = NULL) {
    g <- ggplot(d, aes(share_dead, {{ y }}, colour = year, linewidth = year, alpha = year)) +
      guides
    if (!is.null(b)) {
      b <- b |> filter(share_dead %in% scen$share_dead) |>
        mutate(year = factor(year, year_levels), x = share_dead + bar_offset[as.character(year)])
      g <- g + geom_linerange(data = b, aes(x = x, ymin = lo, ymax = hi, colour = year, linewidth = year),
                              inherit.aes = FALSE, alpha = 0.45)
    }
    g +
      geom_line() +
      geom_point(data = d_pts, aes(size = year)) +
      year_scales +
      xscale +
      labs(x = NULL, y = ylab) +
      base_theme
  }
  keep <- function(d) d |> filter(share_dead %in% scen$share_dead)
  # the period's value at each of the three scenarios, beside the point: above it
  # and to its left, away from the rising line, on a white ground so that the
  # years' lines pass behind it
  period_label <- function(d, value, fmt, nudge_y) {
    if (!period_labels) return(NULL)
    d <- d |> filter(year == period, scenario %in% main_scenarios) |> mutate(value = {{ value }})
    stopifnot(nrow(d) == length(main_scenarios))
    geom_label(data = d, aes(x = share_dead - 0.012, y = value + nudge_y, label = fmt(value)),
               inherit.aes = FALSE, hjust = 1, vjust = 0, size = 3, colour = pal[[period]], fontface = "bold",
               fill = alpha("white", 0.85), linewidth = 0, label.padding = unit(0.08, "lines"))
  }
  p_loss <-
    lines_of(line_m, keep(pts_m), loss, "Loss of male life expectancy at birth (years)", bars_loss_m) +
    period_label(pts_m, loss, \(v) sprintf("%.2f", v), 0.06) +
    scale_y_continuous(breaks = seq(0, 8, 0.5)) +
    theme(legend.position = c(0.1, 0.8), legend.background = element_rect(fill = "white", colour = NA)) +
    no_x
  # combatant deaths in two panels, the years on an axis of their own so that
  # they can be told apart, the four-year total below
  p_mil_years <-
    lines_of(mil_lines |> filter(year != period), keep(mil_pts) |> filter(year != period), military, "Combatant deaths,\nby year",
             bars_mil |> filter(year != period)) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, 80000, 10000)) +
    theme(legend.position = "none") +
    no_x
  p_mil_total <-
    lines_of(mil_lines |> filter(year == period), keep(mil_pts) |> filter(year == period), military, "Combatant deaths,\n2022–2025",
             bars_mil |> filter(year == period)) +
    period_label(mil_pts, military, \(v) scales::comma(round(v, -2)), 2500) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, 200000, 25000)) +
    labs(x = "Share of the missing counted dead") +
    theme(legend.position = "none")
  fig <- p_strip / p_loss / p_mil_years / p_mil_total +
    plot_layout(heights = c(if (letters) 0.8 else 0.45 * rows + 0.4, 3.2, 1.8, 1.4))
  if (!is.null(note)) {
    fig <- fig + plot_annotation(caption = note,
                                 theme = theme(plot.caption = element_text(hjust = 0, size = 8.3, lineheight = 1.15, colour = "grey15")))
  }
  ggsave(file, fig, width = 8, height = height, dpi = 200, bg = "white")
  invisible(fig)
}

mil_lines <- line_military |> mutate(year = factor(year, year_levels))
mil_pts <- military_long |> filter(!is.na(scenario)) |> mutate(year = factor(year, year_levels))
dir.create("figures/exploratory", showWarnings = FALSE, recursive = TRUE)

# the main text: the three scenarios
scen_main <-
  scenarios |>
  filter(scenario %in% main_scenarios) |>
  mutate(row = c(1, 1, 1), hjust = c(0, 1, 1))
# The main text's figure: four panels, the years on the left (the male loss and
# combatant deaths, each year's value labelled at the marks) and 2022-2025 on
# the right. It marks the two ends of the axis - the floor, none of the missing
# dead (the named dead alone), and the ceiling, all of them - and the central
# scenario, the register's observed resolution, between them. The supplement's
# figure (below) keeps every scenario the text names. Combatant deaths in
# thousands. The loss panels share one scale.
# letters = TRUE is the supplement's version: every mark lettered in the strip
# (the key in `note`, under the figure), the projected models as a bracket, and
# no numbers on the panels, which would crowd seven marks; the main text's three
# marks (main) keep the heavier guide.
four_panel_figure <- function(marks, file, letters = FALSE, models = NULL, note = NULL, height = 7.8) {
  if (!"main" %in% names(marks)) marks$main <- TRUE
  keep4 <- function(d) d |> filter(scenario %in% marks$scenario)
  # light reference lines at a quarter, half and three quarters of the missing
  # dead, unlabelled, under the marks' dotted guides
  guides4 <- list(
    geom_vline(xintercept = c(0.25, 0.5, 0.75), colour = "grey60", alpha = 0.3, linewidth = 0.45),
    geom_vline(data = marks |> filter(!main), aes(xintercept = share_dead), linetype = "dotted", colour = "grey55", linewidth = 0.5),
    geom_vline(data = marks |> filter(main), aes(xintercept = share_dead), linetype = "dotted", colour = "grey30", linewidth = 0.9))
  w4 <- c("2022" = 0.6, "2023" = 0.6, "2024" = 0.6, "2025" = 0.6, "2022–2025" = 1.6)
  s4 <- c("2022" = 1.5, "2023" = 1.5, "2024" = 1.5, "2025" = 1.5, "2022–2025" = 2.8)
  scales4 <- list(
    scale_colour_manual(values = pal, name = "Year", drop = TRUE, aesthetics = c("colour", "segment.colour")),
    scale_linewidth_manual(values = w4, name = "Year", drop = TRUE),
    scale_size_manual(values = s4, name = "Year", drop = TRUE))
  # a strip over each column: its title, and the marks by short name and share,
  # or by letter; letters closer than 6 points apart take a second row, and the
  # projected models' bracket sits beneath the letters, its own letter to its left
  strip4 <- function(title) {
    if (!letters) {
      p <- ggplot(marks, aes(share_dead, 0.55)) +
        geom_text(aes(label = sprintf("%s\n(%.0f%%)", short, 100 * share_dead), hjust = hj),
                  size = 3, colour = "grey15", lineheight = 0.9, fontface = "bold") +
        annotate("text", x = 0, y = 1.55, label = title, hjust = 0, size = 4.2, fontface = "bold")
      ylim <- c(0, 1.9)
    } else {
      lab <- marks |>
        arrange(share_dead) |>
        mutate(tight = coalesce(share_dead - lag(share_dead) < 0.06, FALSE),
               y = if_else(tight, 1.05, 0.5),
               hjl = case_when(share_dead < 0.03 ~ 0, share_dead > 0.97 ~ 1, .default = 0.5),
               label = sprintf("%s\n(%.0f%%)", letter, 100 * share_dead))
      p <- ggplot(lab, aes(share_dead, y)) +
        geom_text(aes(label = label, hjust = hjl, size = main), colour = "grey15", lineheight = 0.85, fontface = "bold") +
        scale_size_manual(values = c(`TRUE` = 3.4, `FALSE` = 2.8), guide = "none") +
        annotate("text", x = 0, y = 1.8, label = title, hjust = 0, size = 4.2, fontface = "bold")
      if (!is.null(models)) {
        p <- p +
          annotate("segment", x = min(models$share_dead), xend = max(models$share_dead), y = 0.08, yend = 0.08, colour = "grey45") +
          annotate("point", x = models$share_dead, y = 0.08, shape = 21, fill = "white", colour = "grey45", size = 1.5) +
          annotate("text", x = min(models$share_dead) - 0.012, y = 0.08, hjust = 1, size = 2.8, colour = "grey45",
                   fontface = "bold", label = models$letter[1])
      }
      ylim <- c(0, 2.1)
    }
    p +
      xscale +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.margin = margin(2, 8, 0, 4))
  }
  # the bars of the years set a little apart; at the ends of the axis, inwards
  b_off <- c("2022" = -0.014, "2023" = -0.007, "2024" = 0, "2025" = 0.007, "2022–2025" = 0)
  lines4 <- function(d, d_pts, y, ylab, b) {
    b <- b |> filter(share_dead %in% marks$share_dead, year %in% as.character(unique(d$year))) |>
      mutate(off = b_off[as.character(year)],
             off = case_when(share_dead < 0.02 ~ off - min(b_off), share_dead > 0.98 ~ off - max(b_off), .default = off),
             year = factor(year, year_levels), x = share_dead + off)
    ggplot(d, aes(share_dead, {{ y }}, colour = year, linewidth = year)) +
      guides4 +
      geom_linerange(data = b, aes(x = x, ymin = lo, ymax = hi, colour = year, linewidth = year),
                     inherit.aes = FALSE, alpha = 0.5) +
      geom_line() +
      geom_point(data = d_pts, aes(size = year)) +
      scales4 +
      xscale +
      labs(x = NULL, y = ylab) +
      theme_minimal(base_size = 10.5) +
      theme(panel.grid.minor = element_blank(), plot.margin = margin(0, 8, 2, 4))
  }
  # each year's value beside its point, bold in the year's colour on a white box
  # with no border (linewidth 0); right of the point except at the ceiling
  year_label4 <- function(d, value, fmt) {
    if (letters) return(NULL)
    d <- keep4(d) |> filter(year != period) |> mutate(value = {{ value }}, right = share_dead < 0.98)
    lab <- function(dd, nudge, hjust) {
      geom_label_repel(data = dd, aes(share_dead, value, label = fmt(value), colour = year, segment.colour = year),
                       inherit.aes = FALSE, nudge_x = nudge, hjust = hjust, direction = "y",
                       size = 3.1, fontface = "bold", fill = "white", linewidth = 0,
                       label.padding = unit(0.06, "lines"), box.padding = 0.15, point.padding = 0.15,
                       min.segment.length = 0, segment.size = 0.35, seed = 1, max.overlaps = Inf, show.legend = FALSE)
    }
    list(lab(d |> filter(right), 0.045, 0), lab(d |> filter(!right), -0.045, 1))
  }
  # the period's value above its point, inwards from the ends of the axis
  period_label4 <- function(d, value, fmt, nudge_y) {
    if (letters) return(NULL)
    d <- keep4(d) |> filter(year == period) |>
      mutate(value = {{ value }}, right = share_dead < 0.02,
             x = share_dead + if_else(right, 0.012, -0.012), hj = if_else(right, 0, 1))
    geom_label(data = d, aes(x = x, y = if_else(right, value - nudge_y, value + nudge_y), label = fmt(value), hjust = hj, vjust = if_else(right, 1, 0)),
               inherit.aes = FALSE, size = 3, colour = pal[[period]], fontface = "bold",
               fill = alpha("white", 0.85), linewidth = 0, label.padding = unit(0.08, "lines"))
  }
  ly <- line_m |> filter(year != period); py <- pts_m |> filter(year != period)
  lp <- line_m |> filter(year == period); pp <- pts_m |> filter(year == period)
  k <- \(d) d |> mutate(military = military / 1000)
  my <- k(mil_lines) |> filter(year != period); mpy <- k(mil_pts) |> filter(year != period)
  mp <- k(mil_lines) |> filter(year == period); mpp <- k(mil_pts) |> filter(year == period)
  bm <- bars_mil |> mutate(across(c(lo, hi), \(x) x / 1000))
  # one scale for the two loss panels
  in_marks <- bars_loss_m |> filter(share_dead %in% marks$share_dead)
  lim <- range(c(line_m$loss, in_marks$lo, in_marks$hi))
  loss_y <- scale_y_continuous(limits = c(floor(lim[1] * 2) / 2, ceiling(lim[2] * 2) / 2), breaks = seq(0, 10, 0.5))
  p_ly <- lines4(ly, keep4(py), loss, "Loss of male life expectancy\nat birth (years)", bars_loss_m) +
    year_label4(py, loss, \(v) sprintf("%.2f", v)) + loss_y +
    theme(legend.position = c(0.12, 0.82), legend.background = element_rect(fill = "white", colour = NA)) + no_x
  p_lp <- lines4(lp, keep4(pp), loss, NULL, bars_loss_m) +
    period_label4(pp, loss, \(v) sprintf("%.2f", v), 0.06) + loss_y + theme(legend.position = "none") + no_x
  p_my <- lines4(my, keep4(mpy), military, "Combatant deaths\n(thousands)", bm) +
    year_label4(mpy, military, \(v) sprintf("%.1f", v)) +
    scale_y_continuous(breaks = seq(0, 80, 10)) + labs(x = "Share of the missing counted dead") +
    theme(legend.position = "none")
  p_mp <- lines4(mp, keep4(mpp), military, NULL, bm) +
    period_label4(mpp, military, \(v) sprintf("%.1f", v), 2.5) +
    scale_y_continuous(breaks = seq(0, 200, 25)) + labs(x = "Share of the missing counted dead") +
    theme(legend.position = "none")
  caption <- if (!is.null(note)) note else str_wrap(paste0(
    "Lower bound (0%): none of the missing dead, the named dead alone. Resolution-based scenario (", sprintf("%.0f%%", 100 * marks$share_dead[marks$scenario == "central"]),
    "), the estimate: the register's resolved cases applied to the unresolved. Upper bound (100%): all of the missing dead. ",
    "Lines: every other input at its central value, one share of the missing for every event year. Bars: the 95% interval ",
    "of the drawn inputs with the share of the missing fixed (", scales::comma(length(sims)), " draws). Numbers: each year's ",
    "value in its colour, and 2022–2025 in red. Deaths in thousands. The bounds are logical limits, not plausible values; the supplement's figure marks the bodies and rule scenarios."), 175)
  fig <- wrap_plots(strip4("By year"), strip4("2022–2025"), p_ly, p_lp, p_my, p_mp,
                    ncol = 2, heights = c(if (letters) 0.75 else 0.55, 3, 2.2)) +
    plot_annotation(caption = caption,
                    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey25", lineheight = 1.1)))
  ggsave(file, fig, width = 12, height = height, dpi = 200, bg = "white")
  invisible(fig)
}
marks_main <-
  scenarios |>
  filter(scenario %in% c("bound_low", "central", "bound_high")) |>
  arrange(share_dead) |>
  mutate(short = c("Lower bound", "Resolution-based", "Upper bound"), hj = c(0, 0.5, 1))
stopifnot(nrow(marks_main) == 3, isTRUE(all.equal(marks_main$share_dead[c(1, 3)], c(0, 1))))
four_panel_figure(marks_main, "figures/exploratory/missing_dead_share_main.png")

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
# the supplement's figure in the main text's four panels, every point lettered
four_panel_figure(scen_sm |> mutate(main = scenario %in% c("bound_low", "central", "bound_high")),
                  "figures/exploratory/missing_dead_share_sm.png", letters = TRUE, models = models, note = note_sm,
                  height = 10.2)
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
message("Done. figures/exploratory/missing_dead_share_main.png and _sm.png, and data_inter/ukr_missing_dead_share_sweep.rds, written.")
