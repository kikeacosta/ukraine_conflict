# ==============================================================================
# STEP 13b - Sensitivity of the missing-combatant imputation
# ==============================================================================
#
# Three things the main Monte Carlo treats in a particular way, quantified here
# instead of being asserted to be small or large:
#   1. alpha, the share of the long-term missing who are alive. 11 draws it
#      inside the range the evidence allows; here it is swept over [0, 1], so
#      the reader can see how far the results move across and beyond that
#      range;
#   2. sampling error in the observed transition rates, which are proportions
#      estimated from finite counts and are not propagated (section 5);
#   3. the linkage rules behind those rates (section 6).
# ==============================================================================
#
# WHY THIS MATTERS
# -----------------
# The military death total is dominated by alpha, not by a measurement. Its
# range comes from comparing official prisoner-of-war figures with the
# prisoners and released prisoners the register records (alpha_evidence(),
# 00_setup.R), with one official count of prisoners held behind it. So besides
# drawing alpha inside that range, the pipeline shows what happens across the
# whole of [0, 1].
#
# This re-runs 09's imputation chain (impute_missing(), 00_setup.R) across a
# grid of alpha from 0 (none of the never-resolved are alive) to 1 (all of
# them are), plus the points of the range, and propagates each military
# total through the same deterministic, mode-only projection step 13 uses for
# its migration check: one scenario per grid point, not a full Monte Carlo
# re-run at each one, because the question is how far the estimate moves.
#
# WHAT DOES NOT CHANGE ACROSS THE GRID
# --------------------------------------
#   - the individually confirmed dead (conf_cmb): observed, not imputed
#   - the age-sex profiles (prop_cmb_dead, prop_cmb_miss), civilian deaths,
#     and migration: all held at their mode
#
# The combatant bounds in ukr_param_table.rds are this same chain evaluated
# at the ends and the centre of the range, so the grid points at those three
# values reproduce them exactly.
#
# INPUTS   data_inter/ukr_ualosses_transition_rates.rds            (from 09)
#          data_inter/ukr_ualosses_window_counts.rds               (from 09)
#          data_inter/ukr_alpha_missing.rds, ukr_ualosses_linkage_checks.rds,
#          ukr_ualosses_imputation_table.rds                       (from 09)
#          data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds (08)
#          the same static inputs as step 13 (param_table, forecast
#          mortality, population, migration, fertility, age-sex profiles)
# OUTPUTS  data_inter/ukr_alpha_sensitivity_military.rds
#          data_inter/ukr_alpha_sensitivity_e0.rds
#          data_inter/ukr_transition_rate_sampling.rds
#          data_inter/ukr_linkage_rules_e0.rds (section 6)
#          figures/exploratory/alpha_sensitivity_military.png
#          figures/exploratory/alpha_sensitivity_e0.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# lifetable(), run_single_sim() and impute_missing() all come from 00_setup.R.

# 1. INPUTS SHARED WITH STEP 13, HELD AT THEIR MODE ===========================
param_table <- read_rds("data_inter/ukr_param_table.rds")
alpha_range <- read_rds("data_inter/ukr_alpha_missing.rds")
alpha_mode <- alpha_range$alpha_mode

exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") %>%
  filter(year %in% 2022:2025, source == "frcst") %>%
  select(-source)

pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") %>%
  filter(reg == "cnt", year == 2022) %>%
  select(-reg)

migs2 <- read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds") %>%
  mutate(ems = -mix) |>
  select(-mix)

asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds")
asfr2 <- asfr %>%
  bind_rows(
    asfr %>% filter(year == 2023) %>% mutate(year = 2024),
    asfr %>% filter(year == 2023) %>% mutate(year = 2025)
  ) %>%
  arrange(year, age)

ohchr2 <- read_rds("data_inter/ukr_ohchr_civilian_casualties.rds") %>%
  select(year, sex, age, prop_cvs = cx)

# two profiles, matching 11: registered deaths take the age distribution of
# the register's confirmed dead, imputed deaths that of the missing they
# come from. Neither profile depends on alpha - only the SIZE of the imputed
# total moves across the grid, not its age-sex shape.
ual_raw <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")

cmb_profile <- function(keep_status, nm) {
  ual_raw |>
    filter(status == keep_status, year %in% 2022:2025) |>
    select(-status) |>
    complete(year = 2022:2025, sex, age = 0:100, fill = list(dx = 0)) |>
    summarise(dx = sum(dx), .by = c(year, sex, age)) |>
    mutate("{nm}" := dx / sum(dx), .by = year) |>
    select(-dx)
}

ual2 <-
  cmb_profile("dead", "prop_cmb_dead") |>
  left_join(cmb_profile("missing", "prop_cmb_miss"), by = c("year", "sex", "age"))

profile_props <- ual2 |> left_join(ohchr2, by = c("year", "sex", "age"))

static_inputs <-
  expand_grid(year = 2022:2025, sex = c("f", "m"), age = 0:100) %>%
  left_join(migs2, by = c("year", "sex", "age")) %>%
  left_join(exp_mort2, by = c("year", "sex", "age")) %>%
  left_join(asfr2, by = c("year", "sex", "age")) %>%
  left_join(profile_props, by = c("year", "sex", "age")) %>%
  replace_na(list(ems = 0, w = 0, mx = 0, fx = 0, prop_cvs = 0,
                  prop_cmb_dead = 0, prop_cmb_miss = 0))

stopifnot(all(static_inputs$mx > 0), !any(is.na(static_inputs)))

# civilians and migration held at their mode throughout - only the combatant
# total moves as alpha moves across the grid
draw_cvs_by_year <- param_table |> filter(role == "civilians") |> select(year, draw_cvs = mode)
draw_mig_by_year <- param_table |> filter(str_starts(role, "mig_")) |>
  summarise(draw_mig = sum(mode), .by = year)
conf_cmb_by_year <- param_table |> filter(role == "combatants") |> select(year, conf_cmb = confirmed)

# 2. THE MISSING-COMBATANT IMPUTATION, RE-RUN ACROSS A GRID OF alpha ==========
# tasas_long is 09's own saved output: the twelve-month resolution rates the
# chain is built from. Re-reading it here (rather than re-deriving it from
# the raw registers) means this script cannot silently drift onto a
# different set of empirical rates than 09 actually used.
tasas_long <- read_rds("data_inter/ukr_ualosses_transition_rates.rds")

# the dead and the missing exactly as 09 imputed them: v19 completed for
# registration lag (08b)
imputation_table <- read_rds("data_inter/ukr_ualosses_imputation_table.rds")
stock_missing_2026 <- imputation_table |> select(year, missing_stock)
confirmados_df <- imputation_table |> select(year, confirmados_stock)

# A regular grid, rounded to 2dp so seq()'s floating-point accumulation cannot
# leave 0.30000000000000004 in a manuscript column, plus the three points of
# the range and the 2.5% and 97.5% points of the PERT 11 draws alpha from (the
# 95% interval of the simulated alpha), all marked so 15 can pick them out.
alpha_q <- function(p) {
  qpert(p, min = alpha_range$alpha_min, mode = alpha_range$alpha_mode, max = alpha_range$alpha_max)
}
range_points <- c(
  min = alpha_range$alpha_min,
  mode = alpha_range$alpha_mode,
  max = alpha_range$alpha_max,
  p2.5 = alpha_q(0.025),
  p97.5 = alpha_q(0.975)
)
alpha_grid <- sort(unique(c(round(seq(0, 1, by = 0.05), 2), range_points)))

point_label <- function(a) {
  lab <- names(range_points)[match(a, range_points)]
  if_else(is.na(lab), NA_character_, lab)
}

military_by_alpha <- map_dfr(alpha_grid, function(a) {
  impute_missing(a, tasas_long, stock_missing_2026) |>
    left_join(confirmados_df, by = "year") |>
    mutate(
      alpha = a,
      range_point = point_label(a),
      total_military = confirmados_stock + imputed_dead
    ) |>
    select(
      alpha, range_point, year, confirmados_stock, missing_stock,
      imputed_dead, imputed_alive, imputed_prisoner, total_military
    )
})

# the range points must reproduce 10's combatant bounds
cmb_bounds <- param_table |> filter(role == "combatants") |> arrange(year)
at_point <- function(p) {
  military_by_alpha |> filter(range_point == p) |> arrange(year) |> pull(total_military)
}
stopifnot(
  isTRUE(all.equal(at_point("mode"), cmb_bounds$mode)),
  isTRUE(all.equal(at_point("min"), cmb_bounds$max)),
  isTRUE(all.equal(at_point("max"), cmb_bounds$min))
)

write_rds(military_by_alpha, "data_inter/ukr_alpha_sensitivity_military.rds")

military_totals <-
  military_by_alpha |>
  summarise(total_military = sum(total_military), .by = c(alpha, range_point)) |>
  arrange(alpha)

cat("\n=== TOTAL MILITARY DEATHS, 2022-2025, BY alpha ===\n")
print(as.data.frame(military_totals))

cat(sprintf(
  "\nRange the evidence allows: alpha %.3f / %.3f / %.3f -> %s / %s / %s\n",
  range_points[["min"]], range_points[["mode"]], range_points[["max"]],
  scales::comma(round(at_point("min") |> sum())),
  scales::comma(round(at_point("mode") |> sum())),
  scales::comma(round(at_point("max") |> sum()))
))
cat(sprintf(
  "Whole of alpha in [0, 1]: %s to %s\n",
  scales::comma(round(min(military_totals$total_military))),
  scales::comma(round(max(military_totals$total_military)))
))

# 3. PROPAGATE EACH alpha TO LIFE EXPECTANCY LOSS (deterministic, at the mode) =
e0_of <- function(d, col) {
  d |>
    select(year, sex, age, mx = all_of(col)) |>
    group_by(year, sex) |>
    do(lifetable(dt_in = .data)) |>
    ungroup() |>
    filter(age == 0) |>
    select(year, sex, ex)
}

# the loss for a military total by year (draw_cmb), every other input at its
# mode; shared by the alpha grid and the linkage rules in section 6
loss_for <- function(draw_cmb_by_year) {
  draws_df <-
    draw_cvs_by_year |>
    left_join(draw_mig_by_year, by = "year") |>
    left_join(draw_cmb_by_year, by = "year") |>
    left_join(conf_cmb_by_year, by = "year") |>
    mutate(sim_id = 1)

  sim <-
    run_single_sim(1, draws_df, static_inputs, pop22_ini) |>
    mutate(
      conflict = civilian + combatant_confirmed + combatant_imputed,
      mx_all = (expected + conflict) / pop,
      mx_bsn = expected / pop
    )

  e0_of(sim, "mx_bsn") |>
    rename(ex_bsn = ex) |>
    left_join(e0_of(sim, "mx_all") |> rename(ex_all = ex), by = c("year", "sex")) |>
    mutate(loss = ex_bsn - ex_all)
}

e0_by_alpha <- map_dfr(alpha_grid, function(a) {
  military_by_alpha |>
    filter(alpha == a) |>
    select(year, draw_cmb = total_military) |>
    loss_for() |>
    mutate(alpha = a, range_point = point_label(a))
})

write_rds(e0_by_alpha, "data_inter/ukr_alpha_sensitivity_e0.rds")

cat("\n=== LIFE EXPECTANCY LOSS AT alpha = 0, THE RANGE, AND 1 ===\n")
print(as.data.frame(
  e0_by_alpha |>
    filter(alpha %in% c(0, 1) | !is.na(range_point)) |>
    arrange(year, sex, alpha)
))

# 4. DIAGNOSTIC PLOTS (exploratory; the manuscript table and figure are
#    assembled in 15 from the two .rds files written above) ===================
range_band <- function() {
  list(
    annotate("rect", xmin = range_points[["min"]], xmax = range_points[["max"]],
             ymin = -Inf, ymax = Inf, alpha = 0.15),
    geom_vline(xintercept = range_points[["mode"]], linetype = "dashed", colour = "grey40")
  )
}

p_mil <-
  military_totals |>
  ggplot(aes(alpha, total_military)) +
  range_band() +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "alpha (share of the never-resolved missing who are alive)",
    y = "Total military deaths, 2022-2025",
    title = "Sensitivity of the military death total to alpha",
    caption = "Shaded: the range the evidence allows. Dashed: its central value."
  ) +
  theme_bw()
ggsave("figures/exploratory/alpha_sensitivity_military.png", p_mil, w = 7, h = 4.5)

p_e0 <-
  e0_by_alpha |>
  mutate(sex = if_else(sex == "f", "Females", "Males")) |>
  ggplot(aes(alpha, loss, colour = factor(year))) +
  range_band() +
  geom_line(linewidth = 1) +
  facet_wrap(~sex, scales = "free_y") +
  labs(
    x = "alpha (share of the never-resolved missing who are alive)",
    y = "Life expectancy loss (years)",
    colour = "Year",
    title = "Sensitivity of the life expectancy loss to alpha",
    caption = "Shaded: the range the evidence allows. Dashed: its central value."
  ) +
  theme_bw()
ggsave("figures/exploratory/alpha_sensitivity_e0.png", p_e0, w = 9, h = 4.5)

# ==============================================================================
# 5. SAMPLING ERROR IN THE OBSERVED TRANSITION RATES
# ==============================================================================
# The rates in ukr_ualosses_transition_rates.rds are proportions estimated from
# finite counts, so they carry sampling error that the pipeline does not
# propagate. This measures it rather than asserting it is small.
#
# Each cohort-year and window between releases is treated as an independent
# multinomial sample over the outcomes, with its observed number at risk (09's
# window counts). Draws come from the Dirichlet posterior with a Jeffreys
# prior, counts + 1/2, generated from independent Gammas so no extra package
# is needed; each draw is composed into twelve-month rates by the same
# ual_compose() as 09 and goes through the same impute_missing() chain at the
# central alpha, so the spread is attributable to rate sampling alone.
#
# WHY THIS IS REPORTED RATHER THAN PROPAGATED: the observed transitions carry
# about a sixth of the imputed dead (see the decomposition printed below); the
# terminal alpha assumption carries the rest. The resulting band is about the
# effect of moving alpha by 0.01, small next to the spread alpha produces, so
# folding it in would add machinery without changing any reported figure.
set.seed(42)
B_rates <- 2000

rdirich <- function(a) {
  g <- rgamma(length(a), shape = a, rate = 1)
  g / sum(g)
}

window_counts <-
  read_rds("data_inter/ukr_ualosses_window_counts.rds") |>
  summarise(n = sum(n), .by = c(year, from_release, to)) |>
  complete(nesting(year, from_release), to = ual_outcomes, fill = list(n = 0))

boot_total <- vapply(seq_len(B_rates), function(b) {
  resampled <-
    window_counts |>
    mutate(n = rdirich(n + 0.5), .by = c(year, from_release)) |>
    ual_compose() |>
    mutate(status2 = ual_chain_state(to)) |>
    summarise(prop = sum(prop), .by = c(year, status2))
  imp <- impute_missing(alpha_mode, resampled, stock_missing_2026) |>
    left_join(confirmados_df, by = "year") |>
    mutate(total = confirmados_stock + imputed_dead)
  sum(imp$total)
}, numeric(1))

point_total <- sum(at_point("mode"))
ci_rates <- quantile(boot_total, c(0.025, 0.975))

# where the imputed dead actually come from: observed transitions carry the
# chain terms, the terminal assumption carries the residual
imputed_dead_total <- sum(impute_missing(alpha_mode, tasas_long, stock_missing_2026)$imputed_dead)
terminal_total <- (1 - alpha_mode) * alpha_range$residual
observed_part <- imputed_dead_total - terminal_total

rate_sampling <- tibble(
  alpha = alpha_mode,
  point_total = point_total,
  sd = sd(boot_total),
  lo = ci_rates[[1]],
  hi = ci_rates[[2]],
  width = ci_rates[[2]] - ci_rates[[1]],
  imputed_dead = imputed_dead_total,
  from_observed_transitions = observed_part,
  from_terminal_assumption = terminal_total,
  n_draws = B_rates
)
write_rds(rate_sampling, "data_inter/ukr_transition_rate_sampling.rds")

cat("\n=== SAMPLING ERROR IN THE OBSERVED TRANSITION RATES ===\n")
cat(sprintf("military total at the central alpha: %s\n", scales::comma(round(point_total))))
cat(sprintf("95%% sampling interval             : %s - %s (width %s, %.2f%% of the total)\n",
            scales::comma(round(ci_rates[[1]])), scales::comma(round(ci_rates[[2]])),
            scales::comma(round(ci_rates[[2]] - ci_rates[[1]])),
            100 * (ci_rates[[2]] - ci_rates[[1]]) / point_total))
cat(sprintf("\nof the %s imputed dead:\n", scales::comma(round(imputed_dead_total))))
cat(sprintf("  from observed transitions    : %s (%.1f%%)\n",
            scales::comma(round(observed_part)), 100 * observed_part / imputed_dead_total))
cat(sprintf("  from the terminal assumption : %s (%.1f%%)\n",
            scales::comma(round(terminal_total)), 100 * terminal_total / imputed_dead_total))
cat(sprintf("\nFor scale, moving alpha by 0.01 shifts the total by %s,\n",
            scales::comma(round(0.01 * alpha_range$residual))))
cat("which is several times the entire sampling band above.\n")

# ==============================================================================
# 6. THE LINKAGE RULES
# ==============================================================================
# 09 recomputes the military total under two alternatives to its linkage rules
# - a person no longer listed held as still missing rather than alive, and
# rates from the people listed as missing in v14 only - each at the central
# alpha its own evidence gives. Each is projected here as the alpha grid is,
# every other input at its mode. The rules used must reproduce the central
# value of the grid.
linkage_e0 <-
  read_rds("data_inter/ukr_ualosses_linkage_checks.rds")$designs |>
  mutate(res = map(by_year, \(d) loss_for(d |> select(year, draw_cmb = total)))) |>
  select(design, alpha_mode, military, res) |>
  unnest(res)
stopifnot(isTRUE(all.equal(
  linkage_e0 |> filter(str_detect(design, "production")) |> arrange(year, sex) |> pull(loss),
  e0_by_alpha |> filter(range_point %in% "mode") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))
write_rds(linkage_e0, "data_inter/ukr_linkage_rules_e0.rds")

cat("\n=== THE LINKAGE RULES: MILITARY TOTAL AND LOSS ===\n")
print(as.data.frame(
  linkage_e0 |>
    mutate(col = paste0(sex, "_", year), loss = round(loss, 3), military = round(military)) |>
    select(design, alpha_mode, military, col, loss) |>
    pivot_wider(names_from = col, values_from = loss)
))

message("\nDone. data_inter/ukr_alpha_sensitivity_military.rds, ",
        "ukr_alpha_sensitivity_e0.rds, ukr_transition_rate_sampling.rds and ",
        "ukr_linkage_rules_e0.rds written; the first two are consumed by 15 for ",
        "the manuscript table and figure.")
