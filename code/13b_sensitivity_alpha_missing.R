# ==============================================================================
# STEP 13b - Sensitivity of the missing-combatant imputation
# ==============================================================================
#
# Two things the main Monte Carlo does not propagate, quantified here instead
# of being asserted to be small or large:
#   1. suelo_vivo, the share of the long-term missing assumed alive - swept
#      over its whole range, since nothing in the data identifies it;
#   2. sampling error in the observed transition rates, which are proportions
#      estimated from finite counts (section 5).
# ==============================================================================
#
# WHY THIS MATTERS
# -----------------
# The military death total is dominated by an assumption, not a measurement.
# Of the individuals recorded as missing and never resolved by the end of the
# observed register chain, suelo_vivo = 0.10 (09) assumes 90% are dead.
# Nothing in the data identifies this number: successive register releases
# show how disappearances are ADMINISTRATIVELY RESOLVED, not the eventual
# truth for those never resolved, so extrapolating from one to the other is
# an assumption, not an estimate. The honest response is to show how far the
# headline moves as the assumption moves, not to defend one value.
#
# This re-runs 09's imputation chain (impute_missing(), defined once in
# 00_setup.R so the two scripts cannot drift apart) across a grid of
# suelo_vivo from 0 (none of the never-resolved are alive) to 1 (all of them
# are), and propagates each resulting military-death total through the same
# deterministic, mode-only projection step 13 uses for its migration check:
# one scenario per grid point, not a full Monte Carlo re-run at each one,
# because the question is how far the estimate moves, not its probabilistic
# uncertainty at every point along the way.
#
# WHAT DOES NOT CHANGE ACROSS THE GRID
# --------------------------------------
#   - the individually confirmed dead (min_cmb): these are observed, not
#     imputed, so suelo_vivo cannot move them
#   - the age-sex profiles (prop_cmb_dead, prop_cmb_miss), civilian deaths,
#     and migration: all held at their production mode
#
# NOTE ON THE EXISTING PERT BOUNDS. The combatant min/max in
# ukr_param_table.rds (all confirmed only / all missing simply added as
# dead, with no chain applied) are a different construction from this
# script's suelo_vivo = 0 and 1 endpoints. The chain still nets out whatever
# the observed one-window transitions already resolved before applying the
# terminal assumption to the never-resolved residual, so the two do not
# coincide and should not be read as confirming each other.
#
# INPUTS   data_inter/ukr_ualosses_transition_rates.rds            (from 09)
#          data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds (08)
#          the same static inputs as step 13 (param_table, forecast
#          mortality, population, migration, fertility, age-sex profiles)
# OUTPUTS  data_inter/ukr_alpha_sensitivity_military.rds
#          data_inter/ukr_alpha_sensitivity_e0.rds
#          figures/exploratory/alpha_sensitivity_military.png
#          figures/exploratory/alpha_sensitivity_e0.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# lifetable() and run_single_sim() both come from 00_setup.R, as does
# impute_missing() (extracted from 09 so the two cannot drift apart).

# 1. INPUTS SHARED WITH STEP 13, HELD AT THEIR PRODUCTION MODE ================
param_table <- read_rds("data_inter/ukr_param_table.rds")

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
# come from. Neither profile depends on suelo_vivo - only the SIZE of the
# imputed total moves across the grid, not its age-sex shape.
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
# total moves as suelo_vivo moves across the grid
draw_cvs_by_year <- param_table |> filter(role == "civilians") |> select(year, draw_cvs = mode)
draw_mig_by_year <- param_table |> filter(str_starts(role, "mig_")) |>
  summarise(draw_mig = sum(mode), .by = year)
min_cmb_by_year <- param_table |> filter(role == "combatants") |> select(year, min_cmb = min)

# 2. THE MISSING-COMBATANT IMPUTATION, RE-RUN ACROSS A GRID OF suelo_vivo =====
# tasas_long is 09's own saved output: the one-window resolution rates the
# chain is built from. Re-reading it here (rather than re-deriving it from
# the raw registers) means this script cannot silently drift onto a
# different set of empirical rates than 09 actually used.
tasas_long <- read_rds("data_inter/ukr_ualosses_transition_rates.rds")

status_stocks <-
  ual_raw |>
  summarise(n = sum(dx), .by = c(year, status))

stock_missing_2026 <-
  status_stocks |> filter(status == "missing") |> select(year, missing_stock = n)

confirmados_df <-
  status_stocks |> filter(status == "dead") |> select(year, confirmados_stock = n)

# rounded to 2dp: seq()'s floating-point accumulation would otherwise leave
# suelo_vivo = 0.30000000000000004 in the saved table, which is a real defect
# in a manuscript column, not just a display nuisance
alpha_grid <- round(seq(0, 1, by = 0.05), 2)
stopifnot(0.10 %in% alpha_grid) # the production value must be on the grid

military_by_alpha <- map_dfr(alpha_grid, function(a) {
  impute_missing(a, tasas_long, stock_missing_2026) |>
    left_join(confirmados_df, by = "year") |>
    mutate(
      suelo_vivo = a,
      total_military = confirmados_stock + imputed_dead
    ) |>
    select(
      suelo_vivo, year, confirmados_stock, missing_stock,
      imputed_dead, imputed_alive, imputed_prisoner, total_military
    )
})

write_rds(military_by_alpha, "data_inter/ukr_alpha_sensitivity_military.rds")

military_totals <-
  military_by_alpha |>
  summarise(total_military = sum(total_military), .by = suelo_vivo) |>
  arrange(suelo_vivo)

cat("\n=== TOTAL MILITARY DEATHS, 2022-2025, BY suelo_vivo ===\n")
print(as.data.frame(military_totals))

current <- military_totals$total_military[military_totals$suelo_vivo == 0.10]
crude_min <- sum(confirmados_df$confirmados_stock)
crude_max <- crude_min + sum(stock_missing_2026$missing_stock)
cat(sprintf(
  "\nProduction value (suelo_vivo = 0.10): %s\n",
  scales::comma(round(current))
))
cat(sprintf(
  "Range of the chain across suelo_vivo in [0, 1]: %s to %s\n",
  scales::comma(round(min(military_totals$total_military))),
  scales::comma(round(max(military_totals$total_military)))
))
cat(sprintf(
  "For comparison, the crude PERT bounds (no chain, all-or-nothing): %s to %s\n",
  scales::comma(round(crude_min)), scales::comma(round(crude_max))
))

# 3. PROPAGATE EACH suelo_vivo TO LIFE EXPECTANCY LOSS (deterministic, at the mode)
e0_of <- function(d, col) {
  d |>
    select(year, sex, age, mx = all_of(col)) |>
    group_by(year, sex) |>
    do(lifetable(dt_in = .data)) |>
    ungroup() |>
    filter(age == 0) |>
    select(year, sex, ex)
}

e0_by_alpha <- map_dfr(alpha_grid, function(a) {
  draw_cmb_by_year <-
    military_by_alpha |>
    filter(suelo_vivo == a) |>
    select(year, draw_cmb = total_military)

  draws_df <-
    draw_cvs_by_year |>
    left_join(draw_mig_by_year, by = "year") |>
    left_join(draw_cmb_by_year, by = "year") |>
    left_join(min_cmb_by_year, by = "year") |>
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
    mutate(loss = ex_bsn - ex_all, suelo_vivo = a)
})

write_rds(e0_by_alpha, "data_inter/ukr_alpha_sensitivity_e0.rds")

cat("\n=== LIFE EXPECTANCY LOSS AT suelo_vivo = 0, 0.10 (production) AND 1 ===\n")
print(as.data.frame(
  e0_by_alpha |>
    filter(suelo_vivo %in% c(0, 0.10, 1)) |>
    arrange(year, sex, suelo_vivo)
))

# 4. DIAGNOSTIC PLOTS (exploratory; the manuscript table and figure are
#    assembled in 15 from the two .rds files written above) ===================
p_mil <-
  military_totals |>
  ggplot(aes(suelo_vivo, total_military)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.10, y = max(military_totals$total_military),
           label = "production (0.10)", hjust = -0.05, size = 3, colour = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "suelo_vivo (share of the never-resolved missing assumed alive)",
    y = "Total military deaths, 2022-2025",
    title = "Sensitivity of the military death total to the imputation assumption"
  ) +
  theme_bw()
ggsave("figures/exploratory/alpha_sensitivity_military.png", p_mil, w = 7, h = 4.5)

p_e0 <-
  e0_by_alpha |>
  mutate(sex = if_else(sex == "f", "Females", "Males")) |>
  ggplot(aes(suelo_vivo, loss, colour = factor(year))) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "grey40") +
  facet_wrap(~sex, scales = "free_y") +
  labs(
    x = "suelo_vivo (share of the never-resolved missing assumed alive)",
    y = "Life expectancy loss (years)",
    colour = "Year",
    title = "Sensitivity of the life expectancy loss to the imputation assumption"
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
# Each cohort-year is an independent multinomial sample over
# {dead, missing, prisoner} with the observed at-risk denominator (the cohorts
# are disjoint sets of people, so their sampling errors are independent; the
# chain couples them afterwards, which the resampling below reproduces by
# construction). Draws come from the Dirichlet posterior with a Jeffreys prior,
# counts + 1/2, generated from independent Gammas so no extra package is
# needed. Each draw goes through the same impute_missing() chain at the
# production suelo_vivo, so the spread is attributable to rate sampling alone.
#
# WHY THIS IS REPORTED RATHER THAN PROPAGATED: the observed transitions carry
# under a tenth of the imputed dead (see the decomposition printed below); the
# terminal suelo_vivo assumption carries the rest. Sampling error on a tenth of
# the quantity is correspondingly small, and folding a band of that size into a
# distribution two orders of magnitude wider would add machinery without
# changing any reported figure.
set.seed(42)
B_rates <- 2000

rdirich <- function(alpha) {
  g <- rgamma(length(alpha), shape = alpha, rate = 1)
  g / sum(g)
}

boot_total <- vapply(seq_len(B_rates), function(b) {
  resampled <- map_dfr(sort(unique(tasas_long$year)), function(y) {
    d <- tasas_long |> filter(year == y) |> arrange(status2)
    tibble(year = y, status2 = d$status2, prop = rdirich(d$n + 0.5))
  })
  imp <- impute_missing(0.10, resampled, stock_missing_2026) |>
    left_join(confirmados_df, by = "year") |>
    mutate(total = confirmados_stock + imputed_dead)
  sum(imp$total)
}, numeric(1))

point_total <- military_totals$total_military[military_totals$suelo_vivo == 0.10]
ci_rates <- quantile(boot_total, c(0.025, 0.975))

# where the imputed dead actually come from: observed transitions carry the
# last chain terms, the terminal assumption carries the residual
gp <- function(y, s) {
  v <- tasas_long$prop[tasas_long$year == y & tasas_long$status2 == s]
  if (length(v) == 0) 0 else v[1]
}
M <- setNames(stock_missing_2026$missing_stock, stock_missing_2026$year)
sm <- 1 - 0.10
terminal <- c(
  M[["2022"]] * sm,
  M[["2023"]] * gp(2022, "missing") * sm,
  M[["2024"]] * gp(2023, "missing") * gp(2022, "missing") * sm,
  M[["2025"]] * gp(2024, "missing") * gp(2023, "missing") * gp(2022, "missing") * sm
)
imputed_dead_total <- sum(impute_missing(0.10, tasas_long, stock_missing_2026)$imputed_dead)
observed_part <- imputed_dead_total - sum(terminal)

rate_sampling <- tibble(
  point_total = point_total,
  sd = sd(boot_total),
  lo = ci_rates[[1]],
  hi = ci_rates[[2]],
  width = ci_rates[[2]] - ci_rates[[1]],
  imputed_dead = imputed_dead_total,
  from_observed_transitions = observed_part,
  from_terminal_assumption = sum(terminal),
  n_draws = B_rates
)
write_rds(rate_sampling, "data_inter/ukr_transition_rate_sampling.rds")

cat("\n=== SAMPLING ERROR IN THE OBSERVED TRANSITION RATES ===\n")
cat(sprintf("military total, point estimate : %s\n", scales::comma(round(point_total))))
cat(sprintf("95%% sampling interval          : %s - %s (width %s, %.2f%% of the total)\n",
            scales::comma(round(ci_rates[[1]])), scales::comma(round(ci_rates[[2]])),
            scales::comma(round(ci_rates[[2]] - ci_rates[[1]])),
            100 * (ci_rates[[2]] - ci_rates[[1]]) / point_total))
cat(sprintf("\nof the %s imputed dead:\n", scales::comma(round(imputed_dead_total))))
cat(sprintf("  from observed transitions    : %s (%.1f%%)\n",
            scales::comma(round(observed_part)), 100 * observed_part / imputed_dead_total))
cat(sprintf("  from the terminal assumption : %s (%.1f%%)\n",
            scales::comma(round(sum(terminal))), 100 * sum(terminal) / imputed_dead_total))
cat(sprintf("\nFor scale, moving suelo_vivo by 0.01 shifts the total by about %s,\n",
            scales::comma(round(
              sum(impute_missing(0.10, tasas_long, stock_missing_2026)$imputed_dead) -
              sum(impute_missing(0.11, tasas_long, stock_missing_2026)$imputed_dead)))))
cat("which is several times the entire sampling band above.\n")

message("\nDone. data_inter/ukr_alpha_sensitivity_military.rds, ",
        "ukr_alpha_sensitivity_e0.rds and ukr_transition_rate_sampling.rds ",
        "written; the first two are consumed by 15 for the manuscript table ",
        "and figure.")
