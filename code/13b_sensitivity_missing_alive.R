# ==============================================================================
# STEP 13b - Sensitivity of the missing-combatant imputation
# ==============================================================================
#
# Two things the main Monte Carlo treats in a particular way, quantified here
# instead of being asserted to be small or large:
#   1. how many of the missing are alive. 11 draws it inside the range the
#      evidence allows (alive_evidence(), 00_setup.R); here it is swept from
#      the fewest alive the evidence allows to all but the projected deaths,
#      so the reader can see how far the results move across and beyond that
#      range;
#   2. the linkage rules behind the resolution model, the model itself and
#      the reading of the prisoner-of-war evidence (section 5).
# ==============================================================================
#
# WHY THIS MATTERS
# -----------------
# The military death total is dominated by how many of the missing are
# alive, not by a measurement. The evidence fixes the prisoners of war among
# the missing only within a range, and nothing counts those alive for other
# reasons. So besides drawing both inside the range the evidence allows, the
# pipeline shows what happens across the whole of it and beyond.
#
# THE SWEEP. The missing alive are the prisoners of war among them (captives),
# those the model projects to leave the register, and a share of the rest,
# the unresolved, alive for other reasons (impute_missing(), 00_setup.R). The
# sweep follows one path through the two free quantities, in straight
# segments: from the floor of the evidence (captives at their floor, none of
# the unresolved alive) to its central values (both at their medians), on to
# its ceiling (both at their ceilings), and then to all the unresolved alive
# with the captives at their ceiling. Each point is labelled by the share of
# all the missing alive. Each point is propagated through the same deterministic,
# mode-only projection step 13 uses for its migration check: one scenario per
# point, not a full Monte Carlo re-run at each one, because the question is
# how far the estimate moves.
#
# WHAT DOES NOT CHANGE ALONG THE SWEEP
# ------------------------------------
#   - the registered dead and their late registrations (conf_cmb): observed
#     and completed for registration lag, not imputed
#   - the model's projection of deaths and of people leaving the register
#   - the age-sex profiles (prop_cmb_dead, prop_cmb_miss), civilian deaths,
#     and migration: all held at their mode
#
# The combatant bounds in ukr_param_table.rds are this same imputation at the
# floor, the centre and the ceiling of the evidence, so those points
# reproduce them exactly.
#
# INPUTS   data_inter/ukr_military_inputs.rds, ukr_ualosses_linkage_checks.rds (09)
#          data_inter/ukr_sim_alive_draws.rds (11): the evidence in each draw,
#          for the 95% interval of the missing alive
#          the same static inputs as step 13 (param_table, forecast
#          mortality, population, migration, fertility, age-sex profiles)
# OUTPUTS  data_inter/ukr_alive_sensitivity_military.rds
#          data_inter/ukr_alive_sensitivity_e0.rds
#          data_inter/ukr_linkage_rules_e0.rds (section 5)
#          figures/exploratory/alive_sensitivity_military.png
#          figures/exploratory/alive_sensitivity_e0.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# 1. INPUTS SHARED WITH STEP 13, HELD AT THEIR MODE ===========================
mi <- mode_projection_inputs()
param_table <- mi$param_table
draw_cvs_by_year <- mi$draws |> select(year, draw_cvs)
draw_mig_by_year <- mi$draws |> select(year, draw_mig)

# 2. THE MISSING ALIVE, SWEPT =================================================
mil <- read_rds("data_inter/ukr_military_inputs.rds")
ev <- mil$evidence
at <- function(captives, other) military_draws(mil, captives, other)

# the four-year parts the path moves along, at the model's estimates and the
# point lag factors
parts <- function(captives, other) {
  d <- at(captives, other)
  tibble(missing = sum(d$missing), unlisted = sum(d$imputed_unlisted),
         alive = sum(d$missing - d$imputed_dead))
}
alive_of <- function(captives, other) parts(captives, other)$alive
missing_total <- parts(ev$captives_central, ev$other_central)$missing

# the path's corners: the floor, the centre and the ceiling of the evidence,
# and all the unresolved alive with the captives at their ceiling
nodes <-
  tibble(point = c("floor", "central", "cap", "all but the projected deaths"),
         captives = c(ev$captives_min, ev$captives_central, ev$captives_max, ev$captives_max),
         other = c(ev$other_min, ev$other_central, ev$other_max, 1)) |>
  mutate(alive = map2_dbl(captives, other, alive_of))
stopifnot(all(diff(nodes$alive) > 0))

# the point of the path with a given number of the missing alive
path_at <- function(alive) {
  k <- min(max(findInterval(alive, nodes$alive), 1), nrow(nodes) - 1)
  seg <- function(t) {
    c(captives = nodes$captives[k] + t * (nodes$captives[k + 1] - nodes$captives[k]),
      other = nodes$other[k] + t * (nodes$other[k + 1] - nodes$other[k]))
  }
  if (alive <= nodes$alive[k]) return(seg(0))
  if (alive >= nodes$alive[k + 1]) return(seg(1))
  seg(uniroot(\(t) { p <- seg(t); alive_of(p[["captives"]], p[["other"]]) - alive },
              c(0, 1), tol = 1e-12)$root)
}

# the points of the evidence, and the 95% interval of the draws: each draw's
# missing alive at the model's estimates and the point factors
alive_draws <- read_rds("data_inter/ukr_sim_alive_draws.rds")
alive_in_draws <-
  military_draws(mil, alive_draws$captives, alive_draws$alive_other) |>
  summarise(alive = sum(missing - imputed_dead), .by = sim_id) |>
  pull(alive)
top_alive <- nodes$alive[nrow(nodes)]
between <-
  tibble(point = c("2.5th percentile of draws", "97.5th percentile of draws"),
         alive = quantile(alive_in_draws, c(0.025, 0.975), names = FALSE)) |>
  bind_rows(tibble(alive = missing_total * c(seq(0.10, 0.30, by = 0.05), seq(0.4, 0.8, by = 0.1))) |>
              filter(alive > nodes$alive[1], alive < top_alive)) |>
  mutate(par = map(alive, path_at),
         captives = map_dbl(par, \(v) v[["captives"]]),
         other = map_dbl(par, \(v) v[["other"]])) |>
  select(-par)
sweep_points <-
  bind_rows(nodes, between) |>
  mutate(share_alive = alive / missing_total) |>
  arrange(alive)
# the path must pass through the evidence's own points
stopifnot(
  isTRUE(all.equal(path_at(nodes$alive[2]), c(captives = ev$captives_central, other = ev$other_central))),
  isTRUE(all.equal(path_at(nodes$alive[3]), c(captives = ev$captives_max, other = ev$other_max)))
)

military_by_alive <-
  sweep_points |>
  mutate(res = map2(captives, other, \(c, o) at(c, o) |> select(year, confirmed, imputed_dead, total_military = military))) |>
  unnest(res)

# the points of the evidence must reproduce 10's combatant bounds
cmb_bounds <- param_table |> filter(role == "combatants") |> arrange(year)
at_point <- function(p) {
  military_by_alive |> filter(point == p) |> arrange(year) |> pull(total_military)
}
stopifnot(
  isTRUE(all.equal(at_point("central"), cmb_bounds$mode)),
  isTRUE(all.equal(at_point("cap"), cmb_bounds$min)),
  isTRUE(all.equal(at_point("floor"), cmb_bounds$max))
)

write_rds(military_by_alive, "data_inter/ukr_alive_sensitivity_military.rds")

military_totals <-
  military_by_alive |>
  summarise(total_military = sum(total_military), .by = c(point, alive, share_alive, captives, other)) |>
  arrange(alive)

cat("\n=== TOTAL MILITARY DEATHS, 2022-2025, BY THE SHARE OF THE MISSING ALIVE ===\n")
print(as.data.frame(military_totals |> mutate(across(c(alive, total_military, captives), round),
                                              share_alive = round(share_alive, 4), other = round(other, 4))))

# 3. PROPAGATE EACH POINT TO LIFE EXPECTANCY LOSS (deterministic, at the mode)
# the loss for a military total by year (draw_cmb) with its registered dead
# (conf_cmb), every other input at its mode; shared by the sweep and the
# linkage rules in section 5
loss_for <- function(mil_by_year) {
  draws_df <-
    draw_cvs_by_year |>
    left_join(draw_mig_by_year, by = "year") |>
    left_join(mil_by_year, by = "year") |>
    mutate(sim_id = 1)
  loss_at_mode(draws_df, mi$static_inputs, mi$pop22_ini) |>
    select(year, sex, ex_bsn = e0_bsn, ex_all = e0_all, loss)
}

e0_by_alive <-
  military_by_alive |>
  nest(.by = c(point, alive, share_alive, captives, other)) |>
  mutate(res = map(data, \(d) loss_for(d |> select(year, draw_cmb = total_military, conf_cmb = confirmed)))) |>
  select(-data) |>
  unnest(res)

write_rds(e0_by_alive, "data_inter/ukr_alive_sensitivity_e0.rds")

cat("\n=== LIFE EXPECTANCY LOSS AT THE POINTS OF THE EVIDENCE ===\n")
print(as.data.frame(e0_by_alive |> filter(!is.na(point)) |> arrange(year, sex, alive)))

# 4. DIAGNOSTIC PLOTS (exploratory; the manuscript table and figure are
#    assembled in 15 from the two .rds files written above) ===================
share_of <- function(p) sweep_points$share_alive[sweep_points$point %in% p]
range_band <- function() {
  list(
    annotate("rect", xmin = share_of("floor"), xmax = share_of("cap"), ymin = -Inf, ymax = Inf, alpha = 0.15),
    geom_vline(xintercept = share_of("central"), linetype = "dashed", colour = "grey40")
  )
}

p_mil <-
  military_totals |>
  ggplot(aes(share_alive, total_military)) +
  range_band() +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Share of the missing who are alive",
    y = "Total military deaths, 2022-2025",
    title = "Sensitivity of the military death total to the missing alive",
    caption = "Shaded: the range the evidence allows. Dashed: its central values."
  ) +
  theme_bw()
ggsave("figures/exploratory/alive_sensitivity_military.png", p_mil, w = 7, h = 4.5)

p_e0 <-
  e0_by_alive |>
  mutate(sex = if_else(sex == "f", "Females", "Males")) |>
  ggplot(aes(share_alive, loss, colour = factor(year))) +
  range_band() +
  geom_line(linewidth = 1) +
  facet_wrap(~sex, scales = "free_y") +
  scale_x_continuous(labels = scales::percent) +
  labs(
    x = "Share of the missing who are alive",
    y = "Life expectancy loss (years)",
    colour = "Year",
    title = "Sensitivity of the life expectancy loss to the missing alive",
    caption = "Shaded: the range the evidence allows. Dashed: its central values."
  ) +
  theme_bw()
ggsave("figures/exploratory/alive_sensitivity_e0.png", p_e0, w = 9, h = 4.5)

# ==============================================================================
# 5. THE LINKAGE RULES
# ==============================================================================
# 09 recomputes the military total under alternatives to its linkage rules,
# its resolution model and its reading of the prisoner-of-war evidence, each
# at the central values of the evidence (09 sets them out). Each is projected
# here as the sweep is, every other input at its mode. The rules used must
# reproduce the central point of the sweep.
linkage_e0 <-
  read_rds("data_inter/ukr_ualosses_linkage_checks.rds")$designs |>
  mutate(res = map(by_year, \(d) loss_for(d |> select(year, draw_cmb = total, conf_cmb = confirmed))),
         # summed over rounded years, as the tables build their totals
         military = map_dbl(by_year, \(d) sum(round(d$total)))) |>
  select(design, captives, military, res) |>
  unnest(res)
stopifnot(isTRUE(all.equal(
  linkage_e0 |> filter(str_detect(design, "production")) |> arrange(year, sex) |> pull(loss),
  e0_by_alive |> filter(point %in% "central") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))
write_rds(linkage_e0, "data_inter/ukr_linkage_rules_e0.rds")

cat("\n=== THE LINKAGE RULES: MILITARY TOTAL AND LOSS ===\n")
print(as.data.frame(
  linkage_e0 |>
    mutate(col = paste0(sex, "_", year), loss = round(loss, 3), military = round(military)) |>
    select(design, military, col, loss) |>
    pivot_wider(names_from = col, values_from = loss)
))

message("\nDone. data_inter/ukr_alive_sensitivity_military.rds, ukr_alive_sensitivity_e0.rds ",
        "and ukr_linkage_rules_e0.rds written; the first two are consumed by 15 for the ",
        "manuscript table and figure.")
