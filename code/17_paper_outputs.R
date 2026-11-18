# ==============================================================================
# STEP 17 - The tables and figures of the paper
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Writes every table and figure that appears in the manuscript, main text and
# Supplement, and nothing else, numbered as the manuscript numbers them:
#
#   outputs/tables/    Tables 1-4 and A1-A8
#   outputs/figures/   Figures 1-5, A1 and A2
#   outputs/manifest.csv   one row per item: its number, title, file, and
#                                the files it is built from
#
# The analysis is done upstream. This script only reads what the earlier steps
# wrote to data_inter/, lays it out and saves it; it computes no estimate of its
# own. outputs/ holds one version of each item, that of the last run: the
# folder is emptied before it is written.
#
# To change a table or figure of the paper, change it here. To add one, add it
# with out_tab() or out_fig(), which also enter it in the manifest.
#
# INPUTS   data_inter/ (07-14), figures/exploratory/ (13s, Figure A1)
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
  ggsave(file.path(fig_dir, name), plot = plot, width = w, height = h, dpi = 300, bg = "white")
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
  ual = "data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds",
  ucdp_events = "data_inter/ucdp_ged_events_slim.rds",
  ucdp_adjusted = "data_inter/ukr_ucdp_invals.rds",
  acled = "data_inter/ukr_acled.rds",
  ohchr = "data_inter/ukr_ohchr_annual_totals.rds",
  life_tables = "data_inter/ukr_life_tables_1989_2021.rds",
  resolution = "data_inter/ukr_ualosses_resolution_model.rds",
  scenarios = "data_inter/ukr_scenarios.rds",
  sweep = "data_inter/ukr_missing_dead_share_sweep.rds",
  lag = "data_inter/ukr_registration_lag.rds",
  param_table = "data_inter/ukr_param_table.rds",
  fig_share = "figures/exploratory/missing_dead_share_main.png"
)
stopifnot(all(file.exists(unlist(f))))
uu <- read_rds(f$uu)
loss_draws <- as_tibble(read_rds(f$loss_draws))
# life expectancy and its loss by year and sex: the resolution-based scenario
# from the draws (mean and 95% interval), and the loss at the two bounds of the
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
  uu$combatants |> transmute(year, `Lower bound` = register_bound_low,
                             `Resolution-based scenario` = register_central,
                             `Upper bound` = register_bound_high) |>
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
# TABLE 3 and TABLE A6 - Conflict deaths by type and sex
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

tA6 <- bind_rows(
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
            loss_lower_bound = dec2(loss_bound_low), loss_upper_bound = dec2(loss_bound_high))
out_tab(t4, "Table 4", "Life expectancy at birth and its loss due to direct conflict deaths, by sex, Ukraine, 2022-2025",
        "table4_life_expectancy_loss.csv", c(f$loss_draws, f$sweep))

# ==============================================================================
# TABLE A1 - The sources, as published and as adjusted
# ==============================================================================
ev <- read_rds(f$ucdp_events) |> mutate(u = coalesce(u, 0)) |> filter(country == "Ukraine", year %in% years)
a1_acled <- read_rds(f$acled) |> filter(year %in% years) |>
  pivot_wider(names_from = role, values_from = dts) |>
  transmute(year, civilians, combatants_both_sides = combatants) |> arrange(year) |> with_total() |>
  mutate(source = "ACLED")
a1_ucdp <- ev |>
  summarise(civilians = sum(c), combatants_both_sides = sum(a + b),
            combatants_ukrainian_side = sum(b[type_of_violence == 1 & side_b == "Government of Ukraine"]),
            unknown_status = sum(u), .by = year) |>
  arrange(year) |> with_total() |> mutate(source = "UCDP")
a1_adj <- read_rds(f$ucdp_adjusted) |> filter(year %in% years) |>
  summarise(dts = sum(dts), .by = c(year, role)) |>
  pivot_wider(names_from = role, values_from = dts) |>
  transmute(year, civilians = round(civilians), combatants_ukrainian_side = round(combatants)) |>
  arrange(year) |> with_total() |> mutate(source = "UCDP adjusted")
a1_ohchr <- read_rds(f$ohchr) |> filter(year %in% years) |> transmute(year, civilians = dts) |>
  arrange(year) |> with_total() |> mutate(source = "OHCHR")
a1_ual <- read_rds(f$ual) |> filter(year %in% years, status %in% c("dead", "missing")) |>
  summarise(n = round(sum(dx)), .by = c(year, status)) |>
  pivot_wider(names_from = status, values_from = n) |>
  transmute(year, combatants_ukrainian_side = dead, missing) |> arrange(year) |> with_total() |>
  mutate(source = "UALosses")
tA1 <- bind_rows(a1_acled, a1_ucdp, a1_adj, a1_ohchr, a1_ual) |>
  select(source, year, civilians, combatants_both_sides, combatants_ukrainian_side, missing, unknown_status)
stopifnot(nrow(tA1) == 25)
out_tab(tA1, "Table A1", "Conflict deaths in the sources, as published and as adjusted, Ukraine, 2022-2025",
        "tableA1_sources.csv", c(f$acled, f$ucdp_events, f$ucdp_adjusted, f$ohchr, f$ual))

# ==============================================================================
# TABLE A2 - Life expectancy at birth, 2000-2025
# ==============================================================================
# Observed to 2021. From 2022: expected without war, and with direct conflict
# deaths in the resolution-based scenario and at the two bounds. The lower bound
# of deaths gives the higher life expectancy.
observed <- read_rds(f$life_tables) |>
  filter(age == 0, year >= 2000) |>
  transmute(year, sex = if_else(sex == "females", "Females", "Males"), observed = dec2(ex))
stopifnot(nrow(observed) == 44)
war_years <- t_e0 |>
  transmute(year, sex = sex_label(sex), expected = dec2(e0_expected),
            with_war_resolution_based = dec2(e0_observed),
            with_war_lower_bound = dec2(e0_expected - loss_bound_low),
            with_war_upper_bound = dec2(e0_expected - loss_bound_high))
tA2 <- bind_rows(observed, war_years) |> arrange(desc(sex), year)
out_tab(tA2, "Table A2", "Life expectancy at birth by sex, Ukraine, 2000-2025",
        "tableA2_life_expectancy_2000_2025.csv", c(f$life_tables, f$loss_draws, f$sweep))

# ==============================================================================
# TABLE A3 - The correction for registration lag
# ==============================================================================
lag <- read_rds(f$lag)
horizons <- c("as registered" = "None: as registered", "to 48 months (used)" = "To 48 months (used)",
              "to 60 months" = "To 60 months", "to 72 months" = "To 72 months")
tA3 <- lag$military |>
  filter(scenario %in% names(horizons)) |>
  summarise(combatant_deaths_resolution_based = sum(round(total)), .by = scenario) |>
  left_join(lag$loss |> filter(sex == "m", scenario %in% names(horizons)) |>
              transmute(scenario, year = paste0("male_loss_", year), loss = dec2(loss)) |>
              pivot_wider(names_from = year, values_from = loss), by = "scenario") |>
  mutate(correction = unname(horizons[scenario]), .before = 1) |>
  arrange(match(scenario, names(horizons))) |>
  select(-scenario)
stopifnot(nrow(tA3) == 4, !anyNA(tA3))
out_tab(tA3, "Table A3", "Combatant deaths and male loss of life expectancy by the correction for registration lag",
        "tableA3_registration_lag.csv", f$lag)

# ==============================================================================
# TABLE A4 - What became of the missing whose status was later settled
# ==============================================================================
res <- read_rds(f$resolution)$resolved
tA4 <- tibble(outcome = c("Recorded as dead", "Found to be a prisoner of war",
                          "No longer listed, beyond list maintenance", "Total settled"),
              persons = round(c(res[["dead"]], res[["prisoner"]], res[["no_longer_listed"]], sum(res)))) |>
  mutate(share = sprintf("%.1f%%", 100 * persons / persons[4]))
stopifnot(abs(res[["dead"]] / sum(res) - 0.696) < 0.001)
out_tab(tA4, "Table A4", "Outcome of the missing whose status was settled in a later release of the register",
        "tableA4_resolution_of_missing.csv", f$resolution)

# ==============================================================================
# TABLE A5 - Net migration: the two estimates and Russia
# ==============================================================================
pt <- read_rds(f$param_table)
tA5 <- pt |> filter(role == "mig_west") |>
  transmute(year, border_crossings = round(crossings), unhcr_totals = round(register), midpoint = round(mode)) |>
  left_join(pt |> filter(role == "mig_ru") |>
              transmute(year, russia_low = round(min), russia_central = round(mode), russia_high = round(max)),
            by = "year") |>
  arrange(year) |>
  mutate(across(starts_with("russia"), function(v) coalesce(v, 0))) |>
  with_total()
out_tab(tA5, "Table A5", "Net migration out of Ukraine by year: the two estimates and the component for Russia",
        "tableA5_net_migration.csv", f$param_table)

out_tab(tA6, "Table A6", "Conflict deaths by year, type and sex, Ukraine, 2022-2025",
        "tableA6_deaths_by_year_type_sex.csv", f$loss_draws)

# ==============================================================================
# TABLE A7 and TABLE A8 - Every death from UCDP
# ==============================================================================
# Two tables, since they hold different things: the deaths by year, which are
# not by sex, and the loss by sex and year.
# the total's mean is the sum of the years' exact means, rounded as in the rows;
# its interval is that of the four years summed within each draw
total_ci <- function(means, draws) ci(sum(round(means)), quantile(colSums(draws), .025), quantile(colSums(draws), .975))
tA7 <- uu$compare |>
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
out_tab(tA7, "Table A7", "Conflict deaths with every death taken from UCDP, beside the estimate, Ukraine, 2022-2025",
        "tableA7_ucdp_only_deaths.csv", f$uu)
tA8 <- uu$ucdp_only_parts |>
  arrange(desc(sex), year) |>
  transmute(sex = sex_label(sex), year, as_estimated = dec2(loss_used),
            civilians_replaced = dec2(loss_used + `civilians alone`),
            combatants_replaced = dec2(loss_used + `combatants alone`),
            both_replaced = dec2(loss_used + both))
stopifnot(nrow(tA8) == 8, !anyNA(tA8))
out_tab(tA8, "Table A8", "Loss of life expectancy with every death taken from UCDP, by sex, Ukraine, 2022-2025",
        "tableA8_ucdp_only_loss.csv", f$uu)

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

# FIGURE 3 and FIGURE A2 - the loss between the two bounds, without and with the
# scenario that takes every death from UCDP
range_dat <- t_e0 |>
  transmute(year, sex, loss, lo = loss_lo, hi = loss_hi, low = loss_bound_low, high = loss_bound_high) |>
  left_join(uu$ucdp_only |> filter(point == "UCDP only, uncertainUCDP mean") |> select(year, sex, ucdp = loss),
            by = c("year", "sex")) |>
  mutate(panel = factor(if_else(sex == "m", "Males", "Females"), c("Males", "Females")))
range_vals <- range_dat |> pivot_longer(c(high, loss, low), names_to = "series", values_to = "y") |>
  filter(sex == "m" | series == "loss")
range_figure <- function(with_ucdp) {
  band <- "#c9ced8"; main <- "#22308f"; alt <- "#c2570c"
  ends <- range_dat |> filter(year == 2025, sex == "m") |>
    pivot_longer(c(high, loss, low, ucdp), names_to = "series", values_to = "y") |>
    mutate(label = c(high = "Upper bound (all of the missing dead)", loss = "Resolution-based scenario",
                     low = "Lower bound (none of the missing dead)", ucdp = "Every death from UCDP")[series],
           x = if_else(series == "loss", 2025.5, 2025.15),
           y = y + c(high = 0, loss = 0, low = -0.12, ucdp = 0.12)[series])
  p <- ggplot(range_dat, aes(year)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = band, alpha = 0.7) +
    geom_line(aes(y = low), colour = "#8b93a3", linewidth = 0.4) +
    geom_line(aes(y = high), colour = "#8b93a3", linewidth = 0.4)
  if (with_ucdp) p <- p +
    geom_line(aes(y = ucdp), colour = alt, linewidth = 0.7, linetype = "22") +
    geom_point(aes(y = ucdp), colour = alt, fill = "white", shape = 23, size = 2.4, stroke = 0.9)
  p +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = main, linewidth = 0.9) +
    geom_line(aes(y = loss), colour = main, linewidth = 0.8) +
    geom_point(aes(y = loss), colour = "white", size = 4.2) +
    geom_point(aes(y = loss), colour = main, size = 3) +
    # the Supplement's version labels the scenario for males only: the bounds' values are in Figure 3
    geom_text(data = if (with_ucdp) filter(range_vals, sex == "m", series == "loss") else range_vals,
              aes(y = y, label = dec2(y),
                  vjust = case_when(series == "high" ~ -0.5, series == "low" ~ 1.4, sex == "f" ~ -0.9, .default = 0.4),
                  hjust = case_when(sex == "f" ~ -0.3, series == "loss" ~ -0.45, .default = 1.15)),
              colour = ink, size = 3.1) +
    # the series are named on the figure in the Supplement only; in the main text the note names them
    geom_text(data = if (with_ucdp) ends else ends[0, ], aes(x = x, y = y, label = label),
              colour = ink, size = 3.2, hjust = 0) +
    geom_blank(data = tibble(year = if (with_ucdp) 2027.9 else 2025, panel = factor("Males", c("Males", "Females")))) +
    facet_wrap(~panel, nrow = 1, scales = "free_x") +
    ggh4x::force_panelsizes(cols = if (with_ucdp) c(2.4, 1) else c(1.5, 1)) +
    scale_x_continuous(breaks = years, expand = expansion(add = 0.45)) +
    scale_y_continuous(limits = c(0, 7.4), breaks = 0:7, expand = c(0, 0)) +
    labs(x = NULL, y = "Loss of life expectancy at birth (years)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#e6e8ee", linewidth = 0.4),
          strip.text = element_text(face = "bold", size = 13, hjust = 0, colour = ink),
          axis.text = element_text(colour = "#4a5060"), axis.title.y = element_text(colour = ink),
          panel.spacing.x = unit(1.2, "lines"))
}
out_fig(range_figure(FALSE), "Figure 3",
        "Loss of life expectancy at birth due to direct conflict deaths by sex, Ukraine, 2022-2025",
        "fig3_loss_of_life_expectancy.png", 9, 5.2, c(f$loss_draws, f$sweep))

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

# FIGURE A1 - the loss and combatant deaths at every share of the missing counted
# dead. Step 13s draws it with the sweep it computes; it is copied as drawn.
stopifnot(file.copy(f$fig_share, file.path(fig_dir, "figA1_loss_by_share_of_missing.png"), overwrite = TRUE))
enter("Figure A1", "Male loss of life expectancy and combatant deaths by the share of the missing counted as dead",
      "figures/figA1_loss_by_share_of_missing.png", c(f$fig_share, "drawn by code/13s_missing_dead_share.R"))

out_fig(range_figure(TRUE), "Figure A2",
        "Loss of life expectancy at birth by sex, with every death taken from UCDP, Ukraine, 2022-2025",
        "figA2_loss_ucdp_only.png", 11, 5.6, c(f$loss_draws, f$sweep, f$uu))

# ==============================================================================
# MANIFEST
# ==============================================================================
manifest <- manifest |>
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
