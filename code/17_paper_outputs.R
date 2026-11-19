# ==============================================================================
# STEP 17 - The tables and figures of the paper
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Writes every table and figure that appears in the manuscript, main text and
# Supplement, and nothing else, numbered as the manuscript numbers them:
#
#   outputs/tables/    Tables 1-4 and A1-A11
#   outputs/figures/   Figures 1-5 and A1-A4
#   outputs/manifest.csv   one row per item: its number, title, file, and
#                                the files it is built from
#
# The analysis is done upstream, and the earlier steps save processed data only.
# This script reads what they wrote to data_inter/, lays it out and saves it; it
# computes no estimate of its own, and no other step writes a table or figure
# of the paper. outputs/ holds one version of each item, that of the last run: the
# folder is emptied before it is written.
#
# To change a table or figure of the paper, change it here. To add one, add it
# with out_tab() or out_fig(), which also enter it in the manifest.
#
# INPUTS   data_inter/ (07-14)
# OUTPUTS  outputs/
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

out_dir <- "outputs"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
for (d in c(out_dir, tab_dir, fig_dir)) dir.create(d, showWarnings = FALSE)
# only what this run writes belongs there
old <- list.files(c(tab_dir, fig_dir), full.names = TRUE)
unlink(old)
# a file left open in another program cannot be removed on Windows: the run goes
# on, and says at the end which files it could not replace or remove
locked <- basename(old[file.exists(old)])

manifest <- tibble(item = character(), title = character(), file = character(), built_from = character())
enter <- function(item, title, file, built_from) {
  manifest <<- add_row(manifest, item = item, title = title, file = file,
                       built_from = paste(built_from, collapse = "; "))
  message("  ", item, ": ", file)
}
out_tab <- function(x, item, title, name, built_from) {
  write_csv(x, file.path(tab_dir, name), na = "")
  enter(item, title, file.path("tables", name), built_from)
}
out_fig <- function(plot, item, title, name, w, h, built_from) {
  ggsave(file.path(fig_dir, name), plot = plot, width = w, height = h, dpi = 600, bg = "white")
  enter(item, title, file.path("figures", name), built_from)
}

years <- 2022:2025
yr <- as.character(years)
comma0 <- function(x) scales::comma(round(x), accuracy = 1)
ci <- function(m, lo, hi, f = comma0) sprintf("%s (%s-%s)", f(m), f(lo), f(hi))
dec2 <- function(x) formatC(x, format = "f", digits = 2)
sex_label <- function(s) if_else(s == "f", "Females", "Males")
# a table by year with a last row for the four years
with_total <- function(d) bind_rows(d |> mutate(year = as.character(year)),
                                    d |> summarise(across(-year, sum)) |> mutate(year = "Total"))
# a table of series by year, with the total over the years summed from the
# rounded years, as the other tables of the project build theirs
by_year_wide <- function(d) {
  d |>
    mutate(across(-year, round)) |>
    pivot_longer(-year, names_to = "series") |>
    pivot_wider(names_from = year, values_from = value) |>
    mutate(Total = rowSums(across(all_of(yr))))
}

# ==============================================================================
# INPUTS
# ==============================================================================
f <- list(
  uu = "data_inter/ukr_uncertain_ucdp.rds",
  loss_draws = "data_inter/ukr_e0_loss_by_cause_draws_2022_2025.rds",
  loss_sum = "data_inter/ukr_e0_loss_by_cause_summary_2022_2025.rds",
  rates = "data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds",
  mig = "data_inter/ukr_migration_decomposition.rds",
  mig_age = "data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds",
  ual = "data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds",
  ucdp_events = "data_inter/ucdp_ged_events_slim.rds",
  ucdp_adjusted = "data_inter/ukr_ucdp_invals.rds",
  acled = "data_inter/ukr_acled.rds",
  ohchr = "data_inter/ukr_ohchr_annual_totals.rds",
  life_tables = "data_inter/ukr_life_tables_1989_2021.rds",
  resolution = "data_inter/ukr_ualosses_resolution_model.rds",
  scenarios = "data_inter/ukr_scenarios.rds",
  sweep = "data_inter/ukr_missing_dead_share_sweep.rds",
  bars = "data_inter/ukr_missing_dead_share_bars.rds",
  lag = "data_inter/ukr_registration_lag.rds",
  param_table = "data_inter/ukr_param_table.rds",
  mariupol = "data_inter/ukr_mariupol_split_e0.rds"
)
stopifnot(all(file.exists(unlist(f))))
uu <- read_rds(f$uu)
loss_draws <- as_tibble(read_rds(f$loss_draws))
# life expectancy and its loss by year and sex: the resolution-based scenario
# from the draws (mean and 95% interval), and the loss in the minimum and the maximum scenario of the
# missing from the central projection (13s)
bounds_loss <- read_rds(f$sweep)$loss |>
  filter(scenario %in% c("bound_low", "bound_high"), year %in% as.character(years)) |>
  transmute(year = as.integer(year), sex, scenario = paste0("loss_", scenario), loss) |>
  pivot_wider(names_from = scenario, values_from = loss)
t_e0 <- loss_draws |>
  summarise(e0_expected = mean(e0_bsn), e0_observed = mean(e0_war), loss = mean(loss_total),
            loss_lo = quantile(loss_total, 0.025), loss_hi = quantile(loss_total, 0.975), .by = c(year, sex)) |>
  left_join(bounds_loss, by = c("year", "sex")) |>
  arrange(year, sex)
stopifnot(nrow(t_e0) == 8, !anyNA(t_e0), all(t_e0$loss_bound_low <= t_e0$loss_bound_high + 1e-9))

# ==============================================================================
# TABLE 1 - Conflict deaths as they enter the estimate
# ==============================================================================
t1 <- bind_rows(
  uu$compare |> transmute(year, Low = pert_min, Best = best, High = pert_max) |>
    by_year_wide() |> mutate(group = "Civilians (UCDP adjusted)", .before = 1),
  uu$combatants |> transmute(year, `Minimum scenario` = register_bound_low,
                             `Resolution-based scenario` = register_central,
                             `Maximum scenario` = register_bound_high) |>
    by_year_wide() |> mutate(group = "Combatants (UALosses)", .before = 1)
)
out_tab(t1, "Table 1", "Conflict deaths as they enter the estimate, Ukraine, 2022-2025",
        "table1_conflict_deaths_inputs.csv", f$uu)

# ==============================================================================
# TABLE 2 - The register: registered deaths, late registrations and the missing
# ==============================================================================
t2 <- read_rds(f$scenarios)$by_year |>
  filter(scenario == "central") |>
  transmute(year, registered_deaths = round(registered), late_registrations = round(late_registrations),
            missing = round(missing), missing_counted_dead_resolution_based = round(imputed_dead)) |>
  arrange(year) |>
  with_total()
stopifnot(nrow(t2) == 5, !anyNA(t2))
out_tab(t2, "Table 2", "Ukrainian combatants in the UALosses register by year of the event, 2022-2025",
        "table2_register_combatants.csv", f$scenarios)

# ==============================================================================
# TABLE 3 and TABLE A8 - Conflict deaths by type and sex
# ==============================================================================
# Late registrations are part of the registered deaths: the register's count is
# fixed, and its interval is that of the deaths estimated to be registered late.
types <- c("Civilians", "Registered combatants", "Missing combatants counted as dead")
dx <- loss_draws |>
  transmute(sim_id, year, sex,
            Civilians = dx_civilian,
            `Registered combatants` = dx_cmb_registered + dx_cmb_late,
            `Missing combatants counted as dead` = dx_cmb_imputed) |>
  pivot_longer(all_of(types), names_to = "type", values_to = "dx")
stopifnot(isTRUE(all.equal(sum(dx$dx), sum(loss_draws$dx_civilian + loss_draws$dx_combatant))))
# sums are taken within each draw before the interval, never by adding bounds
cells <- function(by) {
  one <- function(d, col) d |>
    summarise(dx = sum(dx), .by = all_of(c("sim_id", by))) |>
    summarise("{col}" := ci(mean(dx), quantile(dx, .025), quantile(dx, .975)), .by = all_of(by))
  out <- one(filter(dx, sex == "f"), "Females") |>
    bind_cols(one(filter(dx, sex == "m"), "Males") |> select(Males)) |>
    bind_cols(one(dx, "Total") |> select(Total))
  out
}
ord <- function(d) d |> mutate(type = factor(type, types)) |> arrange(type)
t3 <- bind_rows(cells("type") |> ord() |> mutate(type = as.character(type)),
                cells(character(0)) |> mutate(type = "Total", .before = 1))
out_tab(t3, "Table 3", "Conflict deaths by type and sex, Ukraine, 2022-2025",
        "table3_deaths_by_type_sex.csv", f$loss_draws)

tA8 <- bind_rows(
  cells(c("type", "year")) |> ord() |> arrange(type, year) |>
    mutate(type = as.character(type), year = as.character(year)),
  cells("year") |> arrange(year) |> mutate(type = "All conflict deaths", year = as.character(year)),
  cells(character(0)) |> mutate(type = "All conflict deaths", year = "2022-2025")
) |> select(type, year, Females, Males, Total)

# ==============================================================================
# TABLE 4 - Life expectancy at birth and its loss
# ==============================================================================
t4 <- t_e0 |>
  arrange(desc(sex), year) |>
  transmute(sex = sex_label(sex), year,
            expected = dec2(e0_expected), with_war = dec2(e0_observed),
            loss = sprintf("%s (%s-%s)", dec2(loss), dec2(loss_lo), dec2(loss_hi)),
            loss_minimum_scenario = dec2(loss_bound_low), loss_maximum_scenario = dec2(loss_bound_high))
out_tab(t4, "Table 4", "Life expectancy at birth and its loss due to direct conflict deaths, by sex, Ukraine, 2022-2025",
        "table4_life_expectancy_loss.csv", c(f$loss_draws, f$sweep))

# ==============================================================================
# TABLE A1 - The sources, as published and as adjusted
# ==============================================================================
ev <- read_rds(f$ucdp_events) |> mutate(u = coalesce(u, 0)) |> filter(country == "Ukraine", year %in% years)
share_acled <- read_rds(f$acled) |> filter(year %in% years) |>
  pivot_wider(names_from = role, values_from = dts) |>
  transmute(year, civilians, combatants_both_sides = combatants) |> arrange(year) |> with_total() |>
  mutate(source = "ACLED")
share_ucdp <- ev |>
  summarise(civilians = sum(c), combatants_both_sides = sum(a + b),
            combatants_ukrainian_side = sum(b[type_of_violence == 1 & side_b == "Government of Ukraine"]),
            unknown_status = sum(u), .by = year) |>
  arrange(year) |> with_total() |> mutate(source = "UCDP")
share_adj <- read_rds(f$ucdp_adjusted) |> filter(year %in% years) |>
  summarise(dts = sum(dts), .by = c(year, role)) |>
  pivot_wider(names_from = role, values_from = dts) |>
  transmute(year, civilians = round(civilians), combatants_ukrainian_side = round(combatants)) |>
  arrange(year) |> with_total() |> mutate(source = "UCDP adjusted")
share_ohchr <- read_rds(f$ohchr) |> filter(year %in% years) |> transmute(year, civilians = dts) |>
  arrange(year) |> with_total() |> mutate(source = "OHCHR")
share_ual <- read_rds(f$ual) |> filter(year %in% years, status %in% c("dead", "missing")) |>
  summarise(n = round(sum(dx)), .by = c(year, status)) |>
  pivot_wider(names_from = status, values_from = n) |>
  transmute(year, combatants_ukrainian_side = dead, missing) |> arrange(year) |> with_total() |>
  mutate(source = "UALosses")
tA1 <- bind_rows(share_acled, share_ucdp, share_adj, share_ohchr, share_ual) |>
  select(source, year, civilians, combatants_both_sides, combatants_ukrainian_side, missing, unknown_status)
stopifnot(nrow(tA1) == 25)
out_tab(tA1, "Table A1", "Conflict deaths in the sources, as published and as adjusted, Ukraine, 2022-2025",
        "tableA1_sources.csv", c(f$acled, f$ucdp_events, f$ucdp_adjusted, f$ohchr, f$ual))

# ==============================================================================
# TABLE A4 - Life expectancy at birth, 2000-2025
# ==============================================================================
# Observed to 2021. From 2022: expected without war, and with direct conflict
# deaths in the resolution-based scenario and in the minimum and maximum scenarios. The minimum
# of deaths gives the higher life expectancy.
observed <- read_rds(f$life_tables) |>
  filter(age == 0, year >= 2000) |>
  transmute(year, sex = if_else(sex == "females", "Females", "Males"), observed = dec2(ex))
stopifnot(nrow(observed) == 44)
war_years <- t_e0 |>
  transmute(year, sex = sex_label(sex), expected = dec2(e0_expected),
            with_war_resolution_based = dec2(e0_observed),
            with_war_minimum_scenario = dec2(e0_expected - loss_bound_low),
            with_war_maximum_scenario = dec2(e0_expected - loss_bound_high))
tA4 <- bind_rows(observed, war_years) |> arrange(desc(sex), year)
out_tab(tA4, "Table A4", "Life expectancy at birth by sex, Ukraine, 2000-2025",
        "tableA4_life_expectancy_2000_2025.csv", c(f$life_tables, f$loss_draws, f$sweep))

# ==============================================================================
# TABLE A5 - The correction for registration lag
# ==============================================================================
lag <- read_rds(f$lag)
horizons <- c("as registered" = "None: as registered", "to 48 months (used)" = "To 48 months (used)",
              "to 60 months" = "To 60 months", "to 72 months" = "To 72 months")
tA5 <- lag$military |>
  filter(scenario %in% names(horizons)) |>
  summarise(combatant_deaths_resolution_based = sum(round(total)), .by = scenario) |>
  left_join(lag$loss |> filter(sex == "m", scenario %in% names(horizons)) |>
              transmute(scenario, year = paste0("male_loss_", year), loss = dec2(loss)) |>
              pivot_wider(names_from = year, values_from = loss), by = "scenario") |>
  mutate(correction = unname(horizons[scenario]), .before = 1) |>
  arrange(match(scenario, names(horizons))) |>
  select(-scenario)
stopifnot(nrow(tA5) == 4, !anyNA(tA5))
out_tab(tA5, "Table A5", "Combatant deaths and male loss of life expectancy by the correction for registration lag",
        "tableA5_registration_lag.csv", f$lag)

# ==============================================================================
# TABLE A6 - What became of the missing whose status was later settled
# ==============================================================================
res <- read_rds(f$resolution)$resolved
tA6 <- tibble(outcome = c("Recorded as dead", "Found to be a prisoner of war",
                          "No longer listed, beyond list maintenance", "Total settled"),
              persons = round(c(res[["dead"]], res[["prisoner"]], res[["no_longer_listed"]], sum(res)))) |>
  mutate(share = sprintf("%.1f%%", 100 * persons / persons[4]))
stopifnot(abs(res[["dead"]] / sum(res) - 0.696) < 0.001)
out_tab(tA6, "Table A6", "Outcome of the missing whose status was settled in a later release of the register",
        "tableA6_resolution_of_missing.csv", f$resolution)

# ==============================================================================
# TABLE A7 - Net migration: the two estimates and Russia
# ==============================================================================
pt <- read_rds(f$param_table)
tA7 <- pt |> filter(role == "mig_west") |>
  transmute(year, border_crossings = round(crossings), unhcr_totals = round(register), midpoint = round(mode)) |>
  left_join(pt |> filter(role == "mig_ru") |>
              transmute(year, russia_low = round(min), russia_central = round(mode), russia_high = round(max)),
            by = "year") |>
  arrange(year) |>
  mutate(across(starts_with("russia"), function(v) coalesce(v, 0))) |>
  with_total()
out_tab(tA7, "Table A7", "Net migration out of Ukraine by year: the two estimates and the component for Russia",
        "tableA7_net_migration.csv", f$param_table)

out_tab(tA8, "Table A8", "Conflict deaths by year, type and sex, Ukraine, 2022-2025",
        "tableA8_deaths_by_year_type_sex.csv", f$loss_draws)

# ==============================================================================
# TABLE A9 and TABLE A10 - Every death from UCDP
# ==============================================================================
# Two tables, since they hold different things: the deaths by year, which are
# not by sex, and the loss by sex and year.
# the total's mean is the sum of the years' exact means, rounded as in the rows;
# its interval is that of the four years summed within each draw
total_ci <- function(means, draws) ci(sum(round(means)), quantile(colSums(draws), .025), quantile(colSums(draws), .975))
tA9 <- uu$compare |>
  transmute(year = as.character(year), civilians_estimate = comma0(central_projection),
            civilians_ucdp_corrected = ci(uu_mean, uu_lo, uu_hi)) |>
  left_join(uu$combatants |> transmute(year = as.character(year),
                                       combatants_resolution_based = comma0(register_central),
                                       combatants_ucdp_corrected = ci(uu_mean, uu_lo, uu_hi)), by = "year") |>
  add_row(year = "Total",
          civilians_estimate = comma0(sum(round(uu$compare$central_projection))),
          civilians_ucdp_corrected = total_ci(uu$compare$uu_mean, uu$draws),
          combatants_resolution_based = comma0(sum(round(uu$combatants$register_central))),
          combatants_ucdp_corrected = total_ci(uu$combatants$uu_mean, uu$combatant_draws))
out_tab(tA9, "Table A9", "Conflict deaths with every death taken from UCDP, beside the estimate, Ukraine, 2022-2025",
        "tableA9_ucdp_only_deaths.csv", f$uu)
tA10 <- uu$ucdp_only_parts |>
  arrange(desc(sex), year) |>
  transmute(sex = sex_label(sex), year, as_estimated = dec2(loss_used),
            civilians_replaced = dec2(loss_used + `civilians alone`),
            combatants_replaced = dec2(loss_used + `combatants alone`),
            both_replaced = dec2(loss_used + both))
stopifnot(nrow(tA10) == 8, !anyNA(tA10))
out_tab(tA10, "Table A10", "Loss of life expectancy with every death taken from UCDP, by sex, Ukraine, 2022-2025",
        "tableA10_ucdp_only_loss.csv", f$uu)

# ==============================================================================
# TABLES A2 and A3 - the age patterns behind the estimate, by five-year age group
# ==============================================================================
# The estimates are made by single year of age; Figures A1 and A2 draw them so.
age_levels <- c(sprintf("%d-%d", seq(0, 80, 5), seq(4, 84, 5)), "85+", "Total")
age_group <- function(a) {
  lo <- pmin(5 * (a %/% 5), 85)
  if_else(lo == 85, "85+", sprintf("%d-%d", lo, lo + 4))
}
# groups and a row for all ages, by sex; the columns rounded after the sums
by_age_group <- function(d, cols) {
  g <- d |> mutate(age_group = age_group(age)) |> summarise(across(all_of(cols), sum), .by = c(sex, age_group))
  bind_rows(g, g |> summarise(across(all_of(cols), sum), .by = sex) |> mutate(age_group = "Total")) |>
    mutate(sex = sex_label(sex), age_group = factor(age_group, age_levels), across(all_of(cols), round)) |>
    arrange(sex, age_group) |>
    mutate(Total = rowSums(across(all_of(cols))))
}
# TABLE A3 - net migration at the central values: the midpoint of the two estimates
# of the western flow, and Russia's central value in 2022. Negative: net outflow.
tA3 <- read_rds(f$mig_age) |>
  select(year, sex, age, mix) |>
  pivot_wider(names_from = year, values_from = mix) |>
  by_age_group(yr)
stopifnot(nrow(tA3) == 2 * length(age_levels), !anyNA(tA3))
out_tab(tA3, "Table A3", "Net migration by age group and sex, Ukraine, 2022-2025",
        "tableA3_net_migration_by_age.csv", f$mig_age)

# TABLE A2 - conflict deaths of the four years by type: means of the draws in the
# resolution-based scenario, late registrations among the registered combatants
tA2 <- read_rds(f$rates) |>
  filter(cause %in% c("civilian", "combatant_confirmed", "combatant_imputed")) |>
  summarise(deaths = sum(dx_mean), .by = c(sex, age, cause)) |>
  pivot_wider(names_from = cause, values_from = deaths) |>
  rename(civilians = civilian, registered_combatants = combatant_confirmed,
         missing_combatants_counted_dead = combatant_imputed) |>
  by_age_group(c("civilians", "registered_combatants", "missing_combatants_counted_dead"))
stopifnot(nrow(tA2) == 2 * length(age_levels), !anyNA(tA2))
out_tab(tA2, "Table A2", "Conflict deaths by age group, sex and type, Ukraine, 2022-2025",
        "tableA2_conflict_deaths_by_age.csv", f$rates)

# TABLE A11 - the siege of Mariupol: the civilian deaths of 2022 and the loss by
# the rule applied to the event's 15,996 deaths of unknown status (13p), every
# other input at its central value. The rule used is the mix of the deaths
# attributed in Mariupol itself (mariupol_siege, 00_setup.R).
mar <- read_rds(f$mariupol)
tA11 <- mar$civilians |>
  transmute(scenario,
            rule = c(`none civilian` = "None civilian",
                     `dyad-year share` = "Proportions of the conflict and year",
                     `own share` = "Proportions observed in Mariupol (used)",
                     `all civilian` = "All civilian")[rule],
            share_counted_as_civilians = scales::percent(event_share, accuracy = 1),
            civilians_2022_best = round(civilians_mode)) |>
  left_join(mar$conflict_deaths |>
              transmute(scenario, civilians_2022_2025_central = round(civilians_central_total)), by = "scenario") |>
  left_join(mar$loss |> filter(year == 2022) |>
              transmute(scenario, sex = if_else(sex == "f", "loss_2022_females", "loss_2022_males"), loss = dec2(loss)) |>
              pivot_wider(names_from = sex, values_from = loss), by = "scenario") |>
  select(-scenario)
stopifnot(nrow(tA11) == 4, !anyNA(tA11), sum(str_detect(mar$civilians$scenario, "own attributed.*\\(used\\)")) == 1)
out_tab(tA11, "Table A11",
        "Civilian deaths and loss of life expectancy in 2022 by the rule applied to the deaths of unknown status of the siege of Mariupol",
        "tableA11_mariupol_rules.csv", f$mariupol)

# ==============================================================================
# FIGURES
# ==============================================================================
ink <- "#1f2430"
theme_paper <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold", size = base + 1))
}
nice_sex <- function(d) d |> mutate(sex = factor(sex_label(sex), c("Females", "Males")))

# FIGURE 1 - registered dead and missing by year of the event: counts and shares
f1 <- read_rds(f$ual) |>
  filter(status %in% c("dead", "missing"), year %in% years) |>
  summarise(n = sum(dx), .by = c(year, status)) |>
  mutate(pct = scales::percent(n / sum(n), accuracy = 0.1), .by = year) |>
  mutate(status = factor(status, c("missing", "dead"), c("Missing", "Dead")))
col_status <- c(Dead = "#e63946", Missing = "#457b9d")
p1a <- ggplot(f1, aes(factor(year), n, fill = status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = col_status) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Persons", fill = NULL, title = "Counts") + theme_paper()
p1b <- ggplot(f1, aes(factor(year), n, fill = status)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(aes(label = pct), position = position_fill(vjust = 0.5), size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = col_status) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Share", fill = NULL, title = "Proportions") + theme_paper()
out_fig(p1a + p1b + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
        "Figure 1", "Ukrainian combatants registered as dead and as missing by year of the event, 2022-2025",
        "fig1_register_dead_missing.png", 9, 4, f$ual)

# FIGURE 2 - expected and with-war death rates by age
col_rate <- c(all = "#FF4D4D", expected = "#4D4D4D")
lab_rate <- c(all = "With direct conflict deaths", expected = "Expected")
fig2 <- read_rds(f$rates) |>
  filter(cause %in% c("all", "expected")) |>
  nice_sex() |>
  ggplot() +
  geom_ribbon(aes(age, ymin = mx_lower, ymax = mx_upper, fill = cause), alpha = 0.35) +
  geom_line(aes(age, mx_mean, colour = cause), linewidth = 0.5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = col_rate, labels = lab_rate) +
  scale_colour_manual(values = col_rate, labels = lab_rate) +
  facet_grid(sex ~ year) +
  labs(x = "Age", y = "Death rate (log scale)", colour = NULL, fill = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"))
out_fig(fig2, "Figure 2", "Expected and with-war death rates by age and sex, Ukraine, 2022-2025",
        "fig2_death_rates_by_age.png", 9, 5, f$rates)

# FIGURE 3 and FIGURE A4 - the loss in the three scenarios of the missing, without
# and with the scenario that takes every death from UCDP. The minimum and the maximum (keys bound_low and bound_high upstream) are
# scenarios, not the ends of an interval: each is a line of its own, named in the
# legend, with the 95% interval of the draws made with the share of the missing
# counted dead held at that limit (13s). Their points are the central
# projection, the resolution-based scenario's the mean of the draws.
range_dat <- t_e0 |>
  transmute(year, sex, loss, lo = loss_lo, hi = loss_hi, low = loss_bound_low, high = loss_bound_high) |>
  left_join(uu$ucdp_only |> filter(point == "UCDP only, uncertainUCDP mean") |> select(year, sex, ucdp = loss),
            by = c("year", "sex")) |>
  mutate(panel = factor(if_else(sex == "m", "Males", "Females"), c("Males", "Females")))
range_bars <- bind_rows(
  read_rds(f$bars)$loss |>
    filter(scenario %in% c("bound_low", "bound_high"), year %in% yr) |>
    transmute(year = as.integer(year), sex, series = if_else(scenario == "bound_low", "low", "high"), lo, hi),
  range_dat |> transmute(year, sex, series = "loss", lo, hi)
)
# for females the three scenarios coincide, so only the resolution-based one is drawn;
# for males the three are set slightly apart within the year, so their bars do not overlap
range_long <- range_dat |>
  select(year, sex, panel, high, loss, low, ucdp) |>
  pivot_longer(c(high, loss, low, ucdp), names_to = "series", values_to = "y") |>
  filter(sex == "m" | series %in% c("loss", "ucdp")) |>
  left_join(range_bars, by = c("year", "sex", "series")) |>
  mutate(x = year + if_else(sex == "m", c(high = -0.08, loss = 0, low = 0.08, ucdp = 0)[series], 0))
stopifnot(!anyNA(range_long$y),
          with(filter(range_long, series != "ucdp"), !anyNA(lo) & all(lo <= y & y <= hi)))
range_figure <- function(with_ucdp) {
  lab <- c(high = "Maximum scenario: all of the missing are dead",
           loss = "Resolution-based scenario",
           low = "Minimum scenario: none of the missing are dead",
           ucdp = "Every death from UCDP")
  col <- c(high = "#c51b7d", loss = "#22308f", low = "#1b9e77", ucdp = "#c2570c")
  if (!with_ucdp) lab <- lab[names(lab) != "ucdp"]
  d <- range_long |> filter(series %in% names(lab)) |> mutate(series = factor(series, names(lab)))
  # the Supplement's version gives the values of the resolution-based scenario only: the others' are in Figure 3
  val <- function(s) filter(d, series == s, !with_ucdp | (s == "loss" & sex == "m"))
  ggplot(d, aes(x, y)) +
    geom_line(aes(colour = series, linewidth = series, alpha = series)) +
    geom_linerange(aes(ymin = lo, ymax = hi, colour = series, linewidth = series), na.rm = TRUE, show.legend = FALSE) +
    geom_point(aes(colour = series, shape = series, fill = series, size = series), stroke = 0.9) +
    # on a white ground: in 2022 the maximum scenario's line runs behind the value
    geom_label(data = val("loss"), aes(x = x + if_else(sex == "m", 0.13, 0.1), label = dec2(y),
                                       vjust = if_else(sex == "f", -0.1, 0.5)),
               hjust = 0, colour = ink, size = 3.1, fill = "white", linewidth = 0,
               label.padding = unit(0.08, "lines")) +
    geom_text(data = val("high"), aes(y = hi, label = dec2(y)), vjust = -0.6, colour = ink, size = 3.1) +
    geom_text(data = val("low"), aes(y = lo, label = dec2(y)), vjust = 1.6, colour = ink, size = 3.1) +
    facet_wrap(~panel, nrow = 1, scales = "free_x") +
    ggh4x::force_panelsizes(cols = c(1.5, 1)) +
    scale_colour_manual(values = col, labels = lab, name = NULL) +
    scale_fill_manual(values = c(high = "white", loss = col[["loss"]], low = "white", ucdp = "white"),
                      labels = lab, name = NULL) +
    scale_shape_manual(values = c(high = 24, loss = 21, low = 25, ucdp = 23), labels = lab, name = NULL) +
    scale_linewidth_manual(values = c(high = 0.45, loss = 0.9, low = 0.45, ucdp = 0.6), labels = lab, name = NULL) +
    scale_alpha_manual(values = c(high = 0.65, loss = 1, low = 0.65, ucdp = 0.85), labels = lab, name = NULL) +
    scale_size_manual(values = c(high = 2.2, loss = 3, low = 2.2, ucdp = 2.4), labels = lab, name = NULL) +
    scale_x_continuous(breaks = years, expand = expansion(add = 0.45)) +
    scale_y_continuous(limits = c(0, 7.6), breaks = 0:7, expand = c(0, 0)) +
    labs(x = NULL, y = "Loss of life expectancy at birth (years)",
         caption = "Bars: 95% uncertainty interval of the simulation draws within each scenario.") +
    guides(colour = guide_legend(ncol = 1)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#e6e8ee", linewidth = 0.4),
          strip.text = element_text(face = "bold", size = 13, hjust = 0, colour = ink),
          axis.text = element_text(colour = "#4a5060"), axis.title.y = element_text(colour = ink),
          panel.spacing.x = unit(1.2, "lines"),
          legend.position = "bottom", legend.justification = "left",
          legend.key.width = unit(2.2, "lines"), legend.text = element_text(colour = ink, size = 10),
          legend.margin = margin(t = -4),
          plot.caption = element_text(hjust = 0, colour = "#4a5060", size = 9.5),
          plot.caption.position = "plot")
}
out_fig(range_figure(FALSE), "Figure 3",
        "Loss of life expectancy at birth due to direct conflict deaths by sex, Ukraine, 2022-2025",
        "fig3_loss_of_life_expectancy.png", 9, 6.2, c(f$loss_draws, f$sweep, f$bars))

# FIGURE 4 - the loss by conflict mortality and net migration, in shares
f4 <- read_rds(f$mig) |>
  select(year, sex, `Conflict mortality` = loss_nomig, `Net migration` = diff) |>
  pivot_longer(-c(year, sex), names_to = "component", values_to = "loss") |>
  mutate(component = factor(component, c("Conflict mortality", "Net migration"))) |>
  nice_sex() |>
  mutate(pct = loss / sum(loss), .by = c(year, sex))
fig4 <- ggplot(f4, aes(factor(year), loss, fill = component)) +
  geom_col(position = "fill", width = 0.65, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)), position = position_fill(vjust = 0.5),
            size = 3, colour = "white", fontface = "bold") +
  facet_wrap(~sex) +
  scale_fill_manual(values = c("Conflict mortality" = "#555555", "Net migration" = "#FF4D4D")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Share of the loss of life expectancy", fill = NULL) + theme_paper()
out_fig(fig4, "Figure 4",
        "Contribution of conflict mortality and net migration to the loss of life expectancy at birth, by sex, Ukraine, 2022-2025",
        "fig4_loss_mortality_net_migration.png", 7, 4, f$mig)

# FIGURE 5 - the loss by type of death, in shares; late registrations are part
# of the registered combatants
col_type <- c("Civilians" = "#66C2A5", "Registered combatants" = "#FC8D62",
              "Missing combatants counted as dead" = "#8DA0CB")
f5 <- read_rds(f$loss_sum) |>
  mutate(type = case_when(cause == "Civilians" ~ types[1],
                          cause %in% c("Registered combatants", "Late registrations (estimated)") ~ types[2],
                          cause == "Missing combatants (imputed)" ~ types[3])) |>
  filter(!is.na(type)) |>
  summarise(loss = sum(loss_mean), .by = c(year, sex, type)) |>
  mutate(type = factor(type, types)) |>
  nice_sex() |>
  mutate(sh = loss / sum(loss), .by = c(year, sex))
stopifnot(nrow(f5) == 24)
fig5 <- ggplot(f5, aes(factor(year), sh, fill = type)) +
  geom_col(position = "fill", width = 0.68, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = if_else(sh >= 0.03, scales::percent(sh, accuracy = 0.1), "")),
            position = position_fill(vjust = 0.5), size = 3, colour = "black") +
  facet_wrap(~sex) +
  scale_fill_manual(values = col_type) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.03, 0.03))) +
  labs(x = NULL, y = "Share of the loss of life expectancy", fill = NULL) + theme_paper()
out_fig(fig5, "Figure 5",
        "Contribution of civilian and combatant deaths to the loss of life expectancy at birth, by sex, Ukraine, 2022-2025",
        "fig5_loss_by_type_of_death.png", 8, 4.4, f$loss_sum)

# FIGURE A3 - the male loss and combatant deaths at every share of the missing
# counted dead, by year (left) and for 2022-2025 as one period (right). Step 13s
# computes what is drawn: the lines, with every other input at its central value
# and one share of the missing for every event year, and the bars, the 95%
# interval of the draws with the share held fixed. The figure is drawn here.
sw <- read_rds(f$sweep)
sw_bars <- read_rds(f$bars)
period <- "2022–2025"
year_levels <- c(yr, period)
share_pal <- c("2022" = "#9ecae1", "2023" = "#4292c6", "2024" = "#08519c", "2025" = "#08306b", "2022–2025" = "#e63946")
share_marks <- sw$scenarios |>
  filter(scenario %in% c("bound_low", "central", "bound_high")) |>
  arrange(share_dead) |>
  mutate(short = c("Minimum scenario", "Resolution-based", "Maximum scenario"), hj = c(0, 0.5, 1))
stopifnot(nrow(share_marks) == 3, isTRUE(all.equal(share_marks$share_dead[c(1, 3)], c(0, 1))),
          setequal(unique(sw$loss$year), year_levels))
share_year <- function(d) d |> mutate(year = factor(year, year_levels))
share_at_marks <- function(d) d |> filter(scenario %in% share_marks$scenario)
share_line_loss <- sw$line_loss |> filter(sex == "m") |> share_year()
share_pts_loss <- sw$loss |> filter(sex == "m") |> share_at_marks() |> share_year()
# combatant deaths in thousands
share_line_mil <- sw$line_military |> share_year() |> mutate(military = military / 1000)
share_pts_mil <- sw$military |> share_at_marks() |> share_year() |> mutate(military = military / 1000)
share_bars_loss <- sw_bars$loss |> filter(sex == "m") |> select(-sex)
share_bars_mil <- sw_bars$military |> mutate(across(c(lo, hi), \(x) x / 1000))
share_figure <- function() {
  xscale <- scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1), limits = c(-0.01, 1.01), expand = c(0, 0))
  no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  # light reference lines at a quarter, half and three quarters, under the marks' dotted guides
  guides <- list(
    geom_vline(xintercept = c(0.25, 0.5, 0.75), colour = "grey60", alpha = 0.3, linewidth = 0.45),
    geom_vline(data = share_marks, aes(xintercept = share_dead), linetype = "dotted", colour = "grey30", linewidth = 0.9))
  scales_year <- list(
    scale_colour_manual(values = share_pal, name = "Year", drop = TRUE, aesthetics = c("colour", "segment.colour")),
    scale_linewidth_manual(values = c("2022" = 0.6, "2023" = 0.6, "2024" = 0.6, "2025" = 0.6, "2022–2025" = 1.6),
                           name = "Year", drop = TRUE),
    scale_size_manual(values = c("2022" = 1.5, "2023" = 1.5, "2024" = 1.5, "2025" = 1.5, "2022–2025" = 2.8),
                      name = "Year", drop = TRUE))
  # a strip over each column: its title, and the marks by name and share
  strip <- function(title) {
    ggplot(share_marks, aes(share_dead, 0.55)) +
      geom_text(aes(label = sprintf("%s\n(%.0f%%)", short, 100 * share_dead), hjust = hj),
                size = 3, colour = "grey15", lineheight = 0.9, fontface = "bold") +
      annotate("text", x = 0, y = 1.55, label = title, hjust = 0, size = 4.2, fontface = "bold") +
      xscale +
      scale_y_continuous(limits = c(0, 1.9), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.margin = margin(2, 8, 0, 4))
  }
  # the bars of the years set a little apart; at the ends of the axis, inwards
  off <- c("2022" = -0.014, "2023" = -0.007, "2024" = 0, "2025" = 0.007, "2022–2025" = 0)
  lines_of <- function(d, d_pts, y, ylab, b) {
    b <- b |> filter(share_dead %in% share_marks$share_dead, year %in% as.character(unique(d$year))) |>
      mutate(o = off[as.character(year)],
             o = case_when(share_dead < 0.02 ~ o - min(off), share_dead > 0.98 ~ o - max(off), .default = o),
             year = factor(year, year_levels), x = share_dead + o)
    ggplot(d, aes(share_dead, {{ y }}, colour = year, linewidth = year)) +
      guides +
      geom_linerange(data = b, aes(x = x, ymin = lo, ymax = hi, colour = year, linewidth = year),
                     inherit.aes = FALSE, alpha = 0.5) +
      geom_line() +
      geom_point(data = d_pts, aes(size = year)) +
      scales_year +
      xscale +
      labs(x = NULL, y = ylab) +
      theme_minimal(base_size = 10.5) +
      theme(panel.grid.minor = element_blank(), plot.margin = margin(0, 8, 2, 4))
  }
  # each year's value beside its point, in the year's colour on a white ground;
  # right of the point except at the maximum scenario
  year_label <- function(d, value, fmt) {
    d <- d |> filter(year != period) |> mutate(value = {{ value }}, right = share_dead < 0.98)
    lab <- function(dd, nudge, hjust) {
      ggrepel::geom_label_repel(data = dd, aes(share_dead, value, label = fmt(value), colour = year, segment.colour = year),
                                inherit.aes = FALSE, nudge_x = nudge, hjust = hjust, direction = "y",
                                size = 3.1, fontface = "bold", fill = "white", linewidth = 0,
                                label.padding = unit(0.06, "lines"), box.padding = 0.15, point.padding = 0.15,
                                min.segment.length = 0, segment.size = 0.35, seed = 1, max.overlaps = Inf,
                                show.legend = FALSE)
    }
    list(lab(d |> filter(right), 0.045, 0), lab(d |> filter(!right), -0.045, 1))
  }
  # the period's value above its point, inwards from the ends of the axis
  period_label <- function(d, value, fmt, nudge_y) {
    d <- d |> filter(year == period) |>
      mutate(value = {{ value }}, right = share_dead < 0.02,
             x = share_dead + if_else(right, 0.012, -0.012), hj = if_else(right, 0, 1))
    geom_label(data = d, aes(x = x, y = if_else(right, value - nudge_y, value + nudge_y), label = fmt(value),
                             hjust = hj, vjust = if_else(right, 1, 0)),
               inherit.aes = FALSE, size = 3, colour = share_pal[[period]], fontface = "bold",
               fill = alpha("white", 0.85), linewidth = 0, label.padding = unit(0.08, "lines"))
  }
  by_year <- function(d) d |> filter(year != period)
  whole <- function(d) d |> filter(year == period)
  # one scale for the two loss panels
  lim <- range(c(share_line_loss$loss, share_bars_loss$lo[share_bars_loss$share_dead %in% share_marks$share_dead],
                 share_bars_loss$hi[share_bars_loss$share_dead %in% share_marks$share_dead]))
  loss_y <- scale_y_continuous(limits = c(floor(lim[1] * 2) / 2, ceiling(lim[2] * 2) / 2), breaks = seq(0, 10, 0.5))
  p_ly <- lines_of(by_year(share_line_loss), by_year(share_pts_loss), loss, "Loss of male life expectancy\nat birth (years)", share_bars_loss) +
    year_label(share_pts_loss, loss, \(v) sprintf("%.2f", v)) + loss_y +
    theme(legend.position = c(0.12, 0.82), legend.background = element_rect(fill = "white", colour = NA)) + no_x
  p_lp <- lines_of(whole(share_line_loss), whole(share_pts_loss), loss, NULL, share_bars_loss) +
    period_label(share_pts_loss, loss, \(v) sprintf("%.2f", v), 0.06) + loss_y + theme(legend.position = "none") + no_x
  p_my <- lines_of(by_year(share_line_mil), by_year(share_pts_mil), military, "Combatant deaths\n(thousands)", share_bars_mil) +
    year_label(share_pts_mil, military, \(v) sprintf("%.1f", v)) +
    scale_y_continuous(breaks = seq(0, 80, 10)) + labs(x = "Share of the missing counted dead") +
    theme(legend.position = "none")
  p_mp <- lines_of(whole(share_line_mil), whole(share_pts_mil), military, NULL, share_bars_mil) +
    period_label(share_pts_mil, military, \(v) sprintf("%.1f", v), 2.5) +
    scale_y_continuous(breaks = seq(0, 200, 25)) + labs(x = "Share of the missing counted dead") +
    theme(legend.position = "none")
  caption <- str_wrap(paste0(
    "Minimum scenario (0%): none of the missing dead, the named dead alone. Resolution-based scenario (",
    sprintf("%.0f%%", 100 * share_marks$share_dead[share_marks$scenario == "central"]),
    "): the register's resolved cases applied to the unresolved. Maximum scenario (100%): all of the missing dead. ",
    "Lines: every other input at its central value, one share of the missing for every event year. Bars: the 95% interval ",
    "of the drawn inputs with the share of the missing fixed (", scales::comma(sw_bars$n_draws), " draws). Numbers: each year's ",
    "value in its colour, and 2022–2025 in red. Deaths in thousands. The minimum and maximum scenarios are logical limits, not plausible values."), 175)
  wrap_plots(strip("By year"), strip("2022–2025"), p_ly, p_lp, p_my, p_mp, ncol = 2, heights = c(0.55, 3, 2.2)) +
    plot_annotation(caption = caption,
                    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey25", lineheight = 1.1)))
}
out_fig(share_figure(), "Figure A3",
        "Male loss of life expectancy and combatant deaths by the share of the missing counted as dead",
        "figA3_loss_by_share_of_missing.png", 12, 7.8, c(f$sweep, f$bars))

out_fig(range_figure(TRUE), "Figure A4",
        "Loss of life expectancy at birth by sex, with every death taken from UCDP, Ukraine, 2022-2025",
        "figA4_loss_ucdp_only.png", 9, 6.5, c(f$loss_draws, f$sweep, f$bars, f$uu))

# FIGURE A2 - net migration by age and sex, at the central values: the midpoint
# of the two estimates of the western flow, and Russia's central value in 2022.
# The four years share one axis, so that the size of 2022 beside the later years shows.
fA2 <- read_rds(f$mig_age) |>
  # those aged 0 at the end of the year were born in it and are given no net
  # migration: drawn, that zero reads as a drop in the curve, so the lines start at age 1
  filter(age >= 1) |>
  nice_sex()
stopifnot(setequal(fA2$year, years), !anyNA(fA2$mix))
figA2 <- ggplot(fA2, aes(age, mix / 1e3, colour = factor(year))) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~sex) +
  scale_colour_manual(values = c("2022" = "#1f2430", "2023" = "#4292c6", "2024" = "#08519c", "2025" = "#e6842a")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Age at the end of the year", y = "Net migration (thousands, by single year of age)", colour = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"))
out_fig(figA2, "Figure A2", "Net migration by age and sex, Ukraine, 2022-2025",
        "figA2_net_migration_by_age.png", 9, 4.4, f$mig_age)

# FIGURE A1 - conflict deaths by age, sex and type, the four years together:
# means of the draws in the resolution-based scenario, late registrations among
# the registered combatants as in Table 3
fA1 <- read_rds(f$rates) |>
  filter(cause %in% c("civilian", "combatant_confirmed", "combatant_imputed")) |>
  mutate(type = factor(c(civilian = "Civilians", combatant_confirmed = "Registered combatants",
                         combatant_imputed = "Missing combatants counted as dead")[cause], names(col_type))) |>
  summarise(deaths = sum(dx_mean), .by = c(sex, age, type)) |>
  nice_sex()
stopifnot(!anyNA(fA1$type), abs(sum(fA1$deaths) - with(loss_draws, sum(dx_civilian + dx_cmb_confirmed + dx_cmb_imputed)) / n_distinct(loss_draws$sim_id)) < 1)
figA1 <- ggplot(fA1, aes(age, deaths, fill = type)) +
  geom_area(colour = "white", linewidth = 0.15) +
  facet_wrap(~sex, scales = "free_y") +
  scale_fill_manual(values = col_type) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Age", y = "Conflict deaths, 2022–2025 (by single year of age)", fill = NULL) +
  theme_paper() + theme(panel.grid.major.x = element_line(colour = "grey92"))
out_fig(figA1, "Figure A1", "Conflict deaths by age, sex and type, Ukraine, 2022-2025",
        "figA1_conflict_deaths_by_age.png", 9, 4.6, f$rates)

# ==============================================================================
# MANIFEST
# ==============================================================================
# in the order of the paper: tables, then figures; main text before the Supplement
manifest <- manifest |>
  arrange(str_detect(item, "^Figure"), str_detect(item, " A\\d"), parse_number(item)) |>
  mutate(made_by = "code/17_paper_outputs.R", n_draws = n_distinct(loss_draws$sim_id),
         run = format(Sys.time(), "%Y-%m-%d %H:%M"),
         git_commit = tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)[1],
                               error = function(e) NA_character_))
write_csv(manifest, file.path(out_dir, "manifest.csv"))
print(as.data.frame(manifest |> select(item, file)))
stopifnot(all(file.exists(file.path(out_dir, manifest$file))))
extra <- setdiff(list.files(c(tab_dir, fig_dir)), basename(manifest$file))
if (length(c(locked, extra))) {
  warning("open in another program, so not replaced or removed: ", paste(union(locked, extra), collapse = ", "),
          ". Close and run again.", call. = FALSE, immediate. = TRUE)
}
message("Done. ", nrow(manifest), " tables and figures in outputs/.")
