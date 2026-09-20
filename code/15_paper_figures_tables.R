# ==============================================================================
# STEP 15 - Manuscript figures and tables
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Assembles every figure and table that appears in the paper, in one place, so
# the manuscript can be regenerated in a single run and the styling stays
# consistent across all of them.
#
# The earlier scripts still draw their own diagnostic plots, but those go to
# figures/exploratory/. Everything written to figures/ by THIS script is a
# manuscript deliverable, numbered as it is numbered in the draft.
#
# FIGURES                                            source step
#   fig1_combatant_deaths_missing.png                08, 09
#   fig2_mortality_rates_by_age.png                  11
#   fig3_life_expectancy_loss.png                    14
#   fig4_decomposition_migration_mortality.png       13
#   fig5_decomposition_by_cause.png                  14   <- new
#   fig6_uncertainty_shares.png                      11, 14 <- new
#   figA1_pert_draw_distributions.png                10, 11
#   figA2_pert_draw_distributions_migration.png      10, 11 <- new
#   figA3_cumulative_migration_draws.png             11   <- new
#   figA6_duration_lexis.png                         09
#
# TABLES (written to tables/ as .csv)
#   table1_source_totals.csv                         07_*, 08
#   table3_missing_imputation.csv                    09
#   table2_pert_input_bounds.csv                     10
#   tableA26_missing_alive_inputs.csv                09
#   table4_conflict_deaths_by_cause.csv              14   <- new
#   table7_life_expectancy_loss.csv                  14
#   tableA1_source_reconciliation.csv                07_ucdp, 07_acled
#   tableA2_status_transitions.csv                   09
#   table5_totals_by_cause.csv                       14   <- new
#   table6_totals_by_year.csv                        14   <- new
#   tableA3_e0_loss_by_cause.csv                     14
#   tableA4_uncertainty_shares.csv                   11, 14 <- new
#   tableA5_military_reconciliation.csv              07_ucdp, 08, 10, 11
#   tableA6_missing_alive_sensitivity.csv            09, 13b
#   figA4_missing_alive_sensitivity.png              09, 13b
#   tableA7_migration_sensitivity.csv                13, 13c
#   figA5_migration_sensitivity.png                  13, 13c
#   tableA8_migration_specification.csv              13c
#   tableA9_linkage_and_chain.csv                    09, 13b
#   tableA10_registration_lag.csv                    08b, 13d
#   tableA11_population_base_and_timing.csv          13e, 13f
#   tableA12_counterfactual_window.csv               13g, 13j
#   tableA13_pert_shape.csv                          13g
#   tableA14_years_of_life_lost.csv                  14b
#   tableA15_adult_mortality_45q15.csv               14b
#   tableA25_civilian_age_profile.csv                13h
#   tableA27_donbas_fighters.csv                     13i
#   tableA28_out_of_sample.csv                       09i (where it has run)
#   table8_structural_sensitivity.csv                13, 13b-13j
#   tableA16_missing_alive_by_lag_horizon.csv        13d
#   tableA17_military_triangulation.csv              09, data_input/official_figures.csv
#   tableA18_returned_prisoners_prior_status.csv     09
#   tableA19_register_dropout.csv                    09
#   tableA20_resolution_hazards.csv                  09
#   tableA21_civil_register_check.csv                11, data_input/official_figures.csv
#
# INPUTS   the .rds products of steps 07-14
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

fig_dir <- "figures"
tab_dir <- "tables"
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

# The plot is passed explicitly rather than relying on ggsave()'s default of
# last_plot(): patchwork compositions only register there when they are
# printed, which does not happen under source(print.eval = FALSE).
save_fig <- function(plot, name, w, h) {
  ggsave(file.path(fig_dir, name), plot = plot, width = w, height = h, dpi = 300)
  message("  figure: ", name)
}
# A .csv left open in Excel (or mid-sync on a cloud drive) is locked for
# writing on Windows. Rather than abort the whole run on one file, fall back
# to a sibling ".new.csv" and report it, so the remaining tables still build.
save_tab <- function(x, name) {
  path <- file.path(tab_dir, name)
  ok <- tryCatch({
    write_csv(x, path)
    TRUE
  }, error = function(e) FALSE)

  if (ok) {
    message("  table : ", name)
  } else {
    alt <- file.path(tab_dir, sub("[.]csv$", ".new.csv", name))
    write_csv(x, alt)
    warning(
      "could not overwrite ", name,
      " (file open elsewhere?). Wrote ", basename(alt), " instead.",
      call. = FALSE, immediate. = TRUE
    )
  }
}

# ==============================================================================
# SHARED STYLING
# ==============================================================================
COL_STATUS <- c("dead" = "#e63946", "missing" = "#457b9d")
COL_CAUSE <- c(
  "Civilians" = "#66C2A5",
  "Registered combatants" = "#FC8D62",
  "Late registrations (estimated)" = "#E5C494",
  "Missing combatants (imputed)" = "#8DA0CB"
)
COL_EFFECT <- c("Mortality" = "#555555", "Migration" = "#FF4D4D")
COL_RATE <- c("all" = "#FF4D4D", "expected" = "#4D4D4D")

theme_paper <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold", size = base + 1),
      plot.caption = element_text(colour = "grey35", size = base - 3)
    )
}
nice_sex <- function(d) {
  d |> mutate(sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males"))
}

# ==============================================================================
# INPUTS
# ==============================================================================
ual <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")
imput <- read_rds("data_inter/ukr_ualosses_imputation_table.rds")
param_table <- read_rds("data_inter/ukr_param_table.rds")
sce <- read_rds("data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds")
loss_draws <- read_rds("data_inter/ukr_e0_loss_by_cause_draws_2022_2025.rds")
loss_sum <- read_rds("data_inter/ukr_e0_loss_by_cause_summary_2022_2025.rds")
deaths_sum <- read_rds("data_inter/ukr_conflict_deaths_by_cause_summary_2022_2025.rds")
mig <- read_rds("data_inter/ukr_migration_decomposition.rds")
alive_mil <- read_rds("data_inter/ukr_alive_sensitivity_military.rds")
alive_e0 <- read_rds("data_inter/ukr_alive_sensitivity_e0.rds")
evidence <- read_rds("data_inter/ukr_alive_missing.rds")

# The draws file is named for the simulation size, so 15 reads the one matching
# the n_sim it was configured with rather than whatever the last run left behind.
draws_file <- sprintf("data_inter/ukr_sim_draws_2022_2025_n%d.rds", n_sim)
if (!file.exists(draws_file)) {
  stop("Run step 11 first (at n_sim = ", n_sim, "): ", draws_file, call. = FALSE)
}
sim_wide <- as_tibble(readRDS(draws_file))

# and it must actually hold that many draws: if it does not, the file has been
# renamed or truncated, and every interval below would be reported at a size it
# was not computed at
n_found <- length(unique(sim_wide$sim_id))
if (n_found != n_sim) {
  stop("draws file holds ", n_found, " simulations but n_sim is ", n_sim,
       ": ", draws_file, call. = FALSE)
}
message("simulation draws available: ", n_sim)

# the drawn input parameters themselves, one row per simulation, year and role
param_draws <- read_rds("data_inter/ukr_sim_param_draws.rds")

MIG_ROLES <- c(mig_west = "Western displacement", mig_ru_by = "Russia / Belarus")
COL_MIG <- c("Western displacement" = "#0A9396", "Russia / Belarus" = "#EE9B00")

# ==============================================================================
# FIGURE 1 - Death and disappearance of Ukrainian combatants
# ==============================================================================
# Two panels: absolute counts, and the same as yearly proportions, so the
# growing share of unresolved disappearances is visible.
f1_dat <-
  ual |>
  filter(status %in% c("dead", "missing")) |>
  summarise(dts = sum(dx), .by = c(year, status)) |>
  mutate(
    prop = dts / sum(dts),
    pct_label = scales::percent(prop, accuracy = 0.1),
    .by = year
  ) |>
  mutate(status = factor(status, levels = c("missing", "dead")))

p1a <-
  f1_dat |>
  ggplot(aes(x = factor(year), y = dts, fill = status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = COL_STATUS) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Counts", fill = "Status", title = "Counts") +
  theme_paper()

p1b <-
  f1_dat |>
  ggplot(aes(x = factor(year), y = dts, fill = status)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(
    aes(label = pct_label),
    position = position_fill(vjust = 0.5),
    size = 3, colour = "black"
  ) +
  scale_fill_manual(values = COL_STATUS) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year", y = "Percentage", fill = "Status", title = "Proportions") +
  theme_paper()

fig1 <- p1a + p1b + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_fig(fig1, "fig1_combatant_deaths_missing.png", 9, 4)

# ==============================================================================
# FIGURE 2 - All-cause and non-conflict mortality rates by age
# ==============================================================================
fig2 <-
  sce |>
  filter(cause %in% c("all", "expected")) |>
  nice_sex() |>
  ggplot() +
  geom_ribbon(
    aes(age, ymin = mx_lower, ymax = mx_upper, fill = cause),
    alpha = 0.35
  ) +
  geom_line(aes(age, mx_median, colour = cause), linewidth = 0.5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = COL_RATE, labels = c("All-cause", "Expected")) +
  scale_colour_manual(values = COL_RATE, labels = c("All-cause", "Expected")) +
  facet_grid(sex ~ year) +
  labs(
    x = "Age", y = "Death rate (log scale)",
    colour = "Cause", fill = "Cause",
    caption = paste0(
      "Median and 95% uncertainty interval across ", n_sim, " draws."
    )
  ) +
  theme_paper() +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(fig2, "fig2_mortality_rates_by_age.png", 9, 5)

# ==============================================================================
# FIGURE 3 - Life expectancy loss
# ==============================================================================
# Plotted as the CHANGE in life expectancy (war minus counterfactual), so the
# values are negative and the bars read as a fall - which is what the figure
# is about. Everywhere else in the pipeline the same quantity is carried as a
# positive "loss"; only the sign of the display is flipped here.
#
# Each box is annotated with its median (bold) and its quartiles, so the
# spread across draws is readable without going back to table 5.
f3_dat <-
  loss_draws |>
  as_tibble() |>
  nice_sex() |>
  mutate(chg = -loss_total)

f3_sum <-
  f3_dat |>
  summarise(
    Q1 = quantile(chg, 0.25),
    Median = median(chg),
    Q3 = quantile(chg, 0.75),
    .by = c(year, sex)
  )

fig3 <-
  f3_dat |>
  ggplot(aes(x = sex, y = chg)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, fill = "white") +
  # seeded, so the figure is the same in every run
  geom_point(position = position_jitter(width = 0.2, seed = 1), alpha = 0.04, size = 0.1, colour = "black") +
  geom_text(
    data = f3_sum,
    aes(y = Median, label = sprintf("%.2f", Median)),
    hjust = -1.5, size = 2.5, colour = "black", fontface = "bold"
  ) +
  geom_text(
    data = f3_sum,
    aes(y = Q1, label = sprintf("%.2f", Q1)),
    vjust = 1.5, hjust = -1.8, size = 1.8, colour = "black", alpha = 0.7
  ) +
  geom_text(
    data = f3_sum,
    aes(y = Q3, label = sprintf("%.2f", Q3)),
    vjust = -1 / 1.5, hjust = -1.8, size = 1.8, colour = "black", alpha = 0.7
  ) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.9))) +
  facet_grid(~year) +
  labs(
    y = "Change in life expectancy at birth (years)", x = NULL,
    caption = paste0(
      "Distribution across ", n_sim,
      " simulation draws. Bold: median. Small: 25th and 75th percentiles."
    )
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(size = 8))
save_fig(fig3, "fig3_life_expectancy_loss.png", 8, 3.8)

# ==============================================================================
# FIGURE 4 - Mortality vs migration decomposition
# ==============================================================================
f4_dat <-
  mig |>
  select(year, sex, Mortality = loss_nomig, Migration = diff) |>
  pivot_longer(c(Mortality, Migration), names_to = "Component", values_to = "Loss") |>
  mutate(Component = factor(Component, levels = c("Mortality", "Migration"))) |>
  nice_sex() |>
  mutate(pct = Loss / sum(Loss), .by = c(year, sex))

# Proportions only. The absolute panel was dominated by the male bars, which
# left the female contribution unreadable and added nothing the shares do not
# already convey. Labels are white here: the fills are dark grey and red, so
# white holds contrast far better than black.
fig4 <-
  f4_dat |>
  ggplot(aes(x = factor(year), y = Loss, fill = Component)) +
  geom_col(position = "fill", width = 0.65, colour = "white", linewidth = 0.2) +
  geom_text(
    aes(label = scales::percent(pct, accuracy = 0.1)),
    position = position_fill(vjust = 0.5),
    size = 3, colour = "white", fontface = "bold"
  ) +
  facet_wrap(~sex) +
  scale_fill_manual(values = COL_EFFECT) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Year",
    y = "Contribution to the life expectancy loss",
    fill = NULL
  ) +
  theme_paper()
save_fig(fig4, "fig4_decomposition_migration_mortality.png", 7, 4)

# ==============================================================================
# FIGURE 5 - Decomposition by cause of conflict death  (NEW)
# ==============================================================================
three_way <- c(
  "Civilians", "Registered combatants", "Late registrations (estimated)",
  "Missing combatants (imputed)"
)

f5_share <-
  loss_sum |>
  filter(cause %in% three_way) |>
  nice_sex() |>
  mutate(sh = loss_median / sum(loss_median), .by = c(year, sex))

# Proportions only, for the same reason as figure 4: on a shared axis the
# female bars are too small to read against the male ones, and the absolute
# magnitudes are already reported in table 3.
#
# Labels stay BLACK here (unlike figure 4): these fills are the light Set2
# pastels, on which black is the higher-contrast choice.
fig5 <-
  f5_share |>
  ggplot(aes(factor(year), sh, fill = cause)) +
  geom_col(position = "fill", width = 0.68, colour = "white", linewidth = 0.2) +
  geom_text(
    aes(label = if_else(sh >= 0.03, scales::percent(sh, accuracy = 0.1), "")),
    position = position_fill(vjust = 0.5), size = 3, colour = "black"
  ) +
  facet_wrap(~sex) +
  scale_fill_manual(values = COL_CAUSE, breaks = three_way) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    x = "Year",
    y = "Share of the life expectancy loss",
    fill = "Cause of loss"
  ) +
  theme_paper()
save_fig(fig5, "fig5_decomposition_by_cause.png", 8, 4.4)

# ==============================================================================
# FIGURE 6 - Where the uncertainty in the loss comes from
# ==============================================================================
# fig4 asks how much of the LOSS is the migration denominator. This asks how
# much of the SPREAD in that loss is migration - a different question, and the
# one that decides which input is worth more data.
#
# The inputs are drawn independently of one another, so the variance of the
# loss partitions additively across them: each share is the variance the
# input carries through its own first-order effect. What that does not
# account for is reported as its own category rather than distributed over
# the others. Civilians and combatants are kept apart rather than summed into
# one "conflict" term. Combatant draws are an order of magnitude larger, so a
# combined term is dominated by combatant variation - which is nearly
# irrelevant to the female loss, where civilians drive almost all of it. The
# combatants are split by what moves them: the evidence on the missing alive,
# the resolution model's estimates and the registration-lag factors, each
# summarised by the military total it gives with the other two at their point
# values (simulation_draws(), 00_setup.R). The counterfactual enters through
# the draw's forecast index for the sex and year.
UNC_LAB <- c(
  cvs   = "Civilian deaths",
  alive = "Missing alive",
  model = "Resolution model",
  lag   = "Registration lag",
  lc    = "Counterfactual forecast",
  mig_w = "Migration: western",
  mig_r = "Migration: Russia / Belarus",
  resid = "Interaction"
)
COL_UNC <- c(
  "Civilian deaths"             = "#AE2012",
  "Missing alive"               = "#333333",
  "Resolution model"            = "#7F7F7F",
  "Registration lag"            = "#B8B8B8",
  "Counterfactual forecast"     = "#94D2BD",
  "Migration: western"          = "#0A9396",
  "Migration: Russia / Belarus" = "#EE9B00",
  "Interaction"                 = "#E9E9E9"
)

# The loss in year y is that year's death RATES: its own deaths over a
# population that is still missing everyone who left in any earlier year. So
# conflict enters as the year's own draw and migration as the cumulative one.
# Using cumulative deaths instead pushes most of the variance into the
# residual, because it is not what drives the rate.
unc_inputs <-
  param_draws |>
  mutate(bucket = case_when(
    role == "civilians" ~ "cvs",
    role == "mil_alive" ~ "alive",
    role == "mil_model" ~ "model",
    role == "mil_lag"   ~ "lag",
    role == "lc_f"      ~ "lc_f",
    role == "lc_m"      ~ "lc_m",
    role == "mig_west"  ~ "mig_w",
    role == "mig_ru_by" ~ "mig_r"
  )) |>
  filter(!is.na(bucket)) |>
  summarise(draw = sum(draw), .by = c(sim_id, year, bucket)) |>
  # mig_ru_by is entered once, in 2022, but the people it removes are still
  # absent in 2025, so the grid is filled with zeros before accumulating
  complete(sim_id, year, bucket, fill = list(draw = 0)) |>
  arrange(sim_id, bucket, year) |>
  mutate(
    val = if_else(bucket %in% c("mig_w", "mig_r"), cumsum(draw), draw),
    .by = c(sim_id, bucket)
  ) |>
  select(sim_id, year, bucket, val) |>
  pivot_wider(names_from = bucket, values_from = val)

# First-order variance index, Var(E[Y | X]) / Var(Y), estimated by binning X
# on its own quantiles. This is model-free on purpose: the loss is a ratio of
# deaths to a population that migration shrinks, so it is not linear in its
# inputs and a regression slope attributes real structure to the residual.
# What is left over here is genuine interaction, not a failure to fit.
first_order <- function(y, x, nbin = 40) {
  br <- unique(quantile(x, probs = seq(0, 1, length.out = nbin + 1)))
  if (length(br) < 3) return(0)
  b <- cut(x, breaks = br, include.lowest = TRUE, labels = FALSE)
  m <- tapply(y, b, mean)
  n <- tapply(y, b, length)
  N <- length(y)
  k <- length(m)
  ss_b <- sum(n * (m - mean(y))^2)
  ss_w <- sum((y - m[as.character(b)])^2)
  # Each bin mean carries its own sampling noise, which inflates the between-bin
  # sum of squares and would otherwise push the shares past 100% and leave a
  # spuriously negative remainder. This is the standard ANOVA correction.
  max(ss_b - (k - 1) * ss_w / (N - k), 0) / N / var(y)
}

unc_shares <-
  loss_draws |>
  select(sim_id, year, sex, loss_total) |>
  left_join(unc_inputs, by = c("sim_id", "year")) |>
  nest(.by = c(year, sex)) |>
  mutate(sh = map2(data, sex, function(d, sx) {
    s <- c(
      cvs   = first_order(d$loss_total, d$cvs),
      alive = first_order(d$loss_total, d$alive),
      model = first_order(d$loss_total, d$model),
      lag   = first_order(d$loss_total, d$lag),
      lc    = first_order(d$loss_total, if (sx == "f") d$lc_f else d$lc_m),
      mig_w = first_order(d$loss_total, d$mig_w),
      mig_r = first_order(d$loss_total, d$mig_r)
    )
    tibble(source = c(names(s), "resid"), share = c(s, 1 - sum(s)))
  })) |>
  select(-data) |>
  unnest(sh) |>
  mutate(source = factor(UNC_LAB[source], levels = unname(UNC_LAB))) |>
  nice_sex()

fig6 <-
  unc_shares |>
  ggplot(aes(factor(year), share, fill = source)) +
  geom_col(width = 0.7) +
  facet_wrap(~sex) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = COL_UNC) +
  labs(
    x = NULL, y = "Share of the variance in the life expectancy loss",
    fill = NULL,
    caption = paste0(
      "Variance of the loss across ", scales::comma(n_sim),
      " draws, partitioned over inputs that are drawn independently."
    )
  ) +
  theme_paper() +
  guides(fill = guide_legend(nrow = 2))
save_fig(fig6, "fig6_uncertainty_shares.png", 9, 4.8)

save_tab(
  unc_shares |>
    mutate(share = scales::percent(share, accuracy = 0.1)) |>
    pivot_wider(names_from = source, values_from = share),
  "tableA4_uncertainty_shares.csv"
)

# ==============================================================================
# FIGURE A1 - Distributions of the simulated conflict death totals
# ==============================================================================
# param_table also carries the two migration components now, and those belong
# in A2, not here.
param_cnf <- param_table |> filter(role %in% c("combatants", "civilians"))

drawn_totals <-
  sim_wide |>
  summarise(
    civilians = sum(civilian),
    combatants = sum(combatant_confirmed + combatant_imputed),
    .by = c(sim_id, year)
  ) |>
  pivot_longer(c(civilians, combatants), names_to = "role", values_to = "dts")

figA1 <-
  drawn_totals |>
  ggplot(aes(dts, colour = role, fill = role)) +
  geom_density(alpha = 0.55) +
  geom_vline(data = param_cnf, aes(xintercept = mode, colour = role),
             linetype = "dashed", show.legend = FALSE) +
  geom_vline(data = param_cnf, aes(xintercept = min, colour = role),
             linetype = "dotted", show.legend = FALSE) +
  geom_vline(data = param_cnf, aes(xintercept = max, colour = role),
             linetype = "dotted", show.legend = FALSE) +
  geom_text(
    data = param_cnf,
    aes(x = -Inf, y = Inf, label = paste0(
      scales::comma(mode, accuracy = 1), "\n(",
      scales::comma(min, accuracy = 1), " - ",
      scales::comma(max, accuracy = 1), ")"
    )),
    colour = "black", hjust = -0.08, vjust = 1.2, size = 2.4,
    show.legend = FALSE, lineheight = 0.95
  ) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_manual(values = c(civilians = "#AE2012", combatants = "#005F73")) +
  scale_colour_manual(values = c(civilians = "#AE2012", combatants = "#005F73")) +
  facet_nested_wrap(role ~ year, scales = "free", ncol = 4) +
  labs(
    x = "Death counts", y = "Density", colour = "Role", fill = "Role",
    caption = paste0(
      "Draws (n = ", scales::comma(n_sim),
      "). Dashed line: mode. Dotted lines: min and max, for combatants the evidence on the missing alive ",
      "at its ends.\nCombatant totals follow from one draw of that evidence per simulation, so their four ",
      "years move together;\nthe resolution model and the registration-lag factors are drawn too, so a draw ",
      "can pass the dotted lines."
    )
  ) +
  theme_paper() +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(figA1, "figA1_pert_draw_distributions.png", 10, 5)

# ==============================================================================
# FIGURE A2 - Distributions of the simulated net migration totals
# ==============================================================================
# The counterpart of A1 for the migration side. It reads differently from A1
# in one respect: a simulation takes ONE quantile per component and holds it
# across all four years, so within a row the panels move together and the
# spread is a coverage scenario rather than four independent accidents. For
# the western component the dotted lines are that year's two readings, which
# every draw blends with one weight across the four years.
# west first, so its four years fill the top row and the single Russia and
# Belarus panel sits on its own below rather than breaking the row
mig_component <- function(r) factor(MIG_ROLES[r], levels = unname(MIG_ROLES))

mig_bounds <-
  param_table |>
  filter(role %in% names(MIG_ROLES)) |>
  mutate(
    component = mig_component(role),
    across(c(min, mode, max), \(x) x / 1e6)
  )

mig_drawn <-
  param_draws |>
  filter(role %in% names(MIG_ROLES)) |>
  mutate(component = mig_component(role), draw = draw / 1e6)

figA2 <-
  mig_drawn |>
  ggplot(aes(draw, colour = component, fill = component)) +
  geom_density(alpha = 0.55) +
  geom_vline(data = mig_bounds, aes(xintercept = mode, colour = component),
             linetype = "dashed", show.legend = FALSE) +
  geom_vline(data = mig_bounds, aes(xintercept = min, colour = component),
             linetype = "dotted", show.legend = FALSE) +
  geom_vline(data = mig_bounds, aes(xintercept = max, colour = component),
             linetype = "dotted", show.legend = FALSE) +
  geom_text(
    data = mig_bounds,
    aes(x = -Inf, y = Inf, label = paste0(
      scales::number(mode, accuracy = 0.01), "\n(",
      scales::number(min, accuracy = 0.01), " - ",
      scales::number(max, accuracy = 0.01), ")"
    )),
    colour = "black", hjust = -0.08, vjust = 1.2, size = 2.4,
    show.legend = FALSE, lineheight = 0.95
  ) +
  scale_fill_manual(values = COL_MIG) +
  scale_colour_manual(values = COL_MIG) +
  facet_nested_wrap(component ~ year, scales = "free", ncol = 4) +
  labs(
    x = "Net migration, net outflow (millions)", y = "Density",
    colour = "Component", fill = "Component",
    caption = paste0(
      "Draws (n = ", scales::comma(n_sim), "). Dashed line: mode. Dotted lines: ",
      "western, the year's two readings, blended with one Beta-PERT weight per ",
      "simulation;\nRussia and Belarus, its PERT minimum and maximum. Each ",
      "component is drawn once per simulation and held across all four years."
    )
  ) +
  theme_paper() +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(figA2, "figA2_pert_draw_distributions_migration.png", 10, 5)

# ==============================================================================
# FIGURE A3 - Cumulative net migration implied by the draws
# ==============================================================================
# A2 shows the components year by year; this shows what they add up to, which
# is the quantity the published sources actually disagree about.
mig_cum <-
  param_draws |>
  filter(role %in% names(MIG_ROLES)) |>
  summarise(total = sum(draw) / 1e6, .by = sim_id)

# both reference values are computed from the saved sources, not typed in
ces_fig <- read_csv("data_input/migration/ces_2026_figures.csv", show_col_types = FALSE)
unhcr_2025 <-
  read_csv("data_input/migration/unhcr/unhcr_population_coo_UKR_2021_2025.csv",
           show_col_types = FALSE) |>
  filter(year == 2025) |>
  summarise(v = sum(coalesce(refugees, 0) + coalesce(asylum_seekers, 0))) |>
  pull(v)

mig_refs <- tribble(
  ~label,                                                ~value,
  "CES 2026 total (west + via RU/BY + in RU/BY)",
  ces_fig$value[ces_fig$key == "refugees_total"] / 1e6,
  "UNHCR Data Finder, end 2025 (excl. US parolees)",   unhcr_2025 / 1e6
)

figA3 <-
  mig_cum |>
  ggplot(aes(total)) +
  geom_density(fill = "#0A9396", colour = "#0A9396", alpha = 0.55) +
  geom_vline(data = mig_refs, aes(xintercept = value, linetype = label),
             colour = "grey20") +
  scale_linetype_manual(values = c("dashed", "dotdash")) +
  guides(linetype = guide_legend(nrow = 2)) +
  labs(
    x = "Cumulative net migration 2022–2025 (millions, net outflow)", y = "Density",
    linetype = "Published estimate",
    caption = paste0(
      "Draws (n = ", scales::comma(n_sim), ").\nThe published ",
      "figures count different populations: the Data Finder\nomits US parolees ",
      "and most of Canada, CES pairs crossings with 1.3M in Russia/Belarus."
    )
  ) +
  theme_paper() +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(figA3, "figA3_cumulative_migration_draws.png", 7, 4)

# ==============================================================================
# TABLE 1 - What each source reports, cumulative 2022-2025
# ==============================================================================
# One row per data source, with its civilian and combatant totals and, where
# the source provides them, low and high bounds. NA means the source does not
# report that quantity at all - it is not a zero.
#
#   OHCHR         verified civilian deaths only. A single documented count,
#                 with no bounds. Combatants are outside its mandate.
#   UCDP          all deaths recorded in events located in Ukraine, in UCDP's
#                 own categories, BEFORE nationality is assigned. Combatants
#                 here therefore include Russian soldiers killed in Ukraine.
#   UCDP adjusted after combatant deaths are reassigned to the combatant's own
#                 nationality and the deaths of unknown side are redistributed
#                 across civilians and combatants.
#   ACLED         point estimates only; ACLED publishes no bounds.
#   UALosses      an individually named register of Ukrainian military
#                 personnel, so it has no civilian counterpart and no natural
#                 mid-point. Low = deaths recorded in the register;
#                 high = those deaths plus everyone still listed as missing.
#
# The last column records whether the source is disaggregated by age and sex,
# which is what determines whether it can supply a mortality PROFILE rather
# than just a total.
ohchr_tot <- read_rds("data_inter/ukr_ohchr_annual_totals.rds")
acled <- read_rds("data_inter/ukr_acled.rds")
ucdp_cmp_full <- read_rds("data_inter/ukr_rus_ucdp_source_comparison.rds")
ucdp_ukr <- ucdp_cmp_full |> filter(country == "Ukraine")

# The register's own low and high, taken from the parameter table's register
# columns and summed over rounded years, the way the other tables build their
# Total rows. The counts are fractional (08 redistributes records of unknown
# age or year), and rounding a four-year sum once can differ by one from the
# sum of rounded years.
cmb_bounds <- param_table |> filter(role == "combatants")
ual_low <- sum(round(cmb_bounds$confirmed))
ual_high <- sum(round(cmb_bounds$listed))

acled_tot <- acled |> summarise(dts = sum(dts), .by = role)

tab_src <- tribble(
  ~source, ~civ_mid, ~civ_low, ~civ_high, ~cmb_mid, ~cmb_low, ~cmb_high, ~age_sex,

  "OHCHR",
  sum(ohchr_tot$dts), NA, NA,
  NA, NA, NA,
  "Yes",

  "ACLED",
  acled_tot$dts[acled_tot$role == "civilians"], NA, NA,
  acled_tot$dts[acled_tot$role == "combatants"], NA, NA,
  "No",
  
  "UCDP",
  ucdp_ukr$civilians[ucdp_ukr$type == "unadjusted"],
  ucdp_ukr$civilians_l[ucdp_ukr$type == "unadjusted"],
  ucdp_ukr$civilians_u[ucdp_ukr$type == "unadjusted"],
  ucdp_ukr$combatants[ucdp_ukr$type == "unadjusted"],
  ucdp_ukr$combatants_l[ucdp_ukr$type == "unadjusted"],
  ucdp_ukr$combatants_u[ucdp_ukr$type == "unadjusted"],
  "No",

  "UCDP adjusted",
  ucdp_ukr$civilians[ucdp_ukr$type == "adjusted"],
  ucdp_ukr$civilians_l[ucdp_ukr$type == "adjusted"],
  ucdp_ukr$civilians_u[ucdp_ukr$type == "adjusted"],
  ucdp_ukr$combatants[ucdp_ukr$type == "adjusted"],
  ucdp_ukr$combatants_l[ucdp_ukr$type == "adjusted"],
  ucdp_ukr$combatants_u[ucdp_ukr$type == "adjusted"],
  "No",

  "UALosses",
  NA, NA, NA,
  NA, ual_low, ual_high,
  "Yes"
) |>
  # plain integers rather than formatted strings: these import as numbers in
  # Excel, and thousands separators would force every cell to be quoted.
  # NA here means the source does not report that quantity, not zero.
  mutate(
    across(
      c(civ_mid, civ_low, civ_high, cmb_mid, cmb_low, cmb_high),
      ~ round(.x)
    )
  ) |>
  rename(
    `Civilians: mid` = civ_mid,
    `Civilians: low` = civ_low,
    `Civilians: high` = civ_high,
    `Combatants: mid` = cmb_mid,
    `Combatants: low` = cmb_low,
    `Combatants: high` = cmb_high,
    `Age-sex disaggregation` = age_sex,
    Source = source
  )

save_tab(tab_src, "table1_source_totals.csv")
print(tab_src)

# ==============================================================================
# TABLE 2 - Imputation of missing combatants
# ==============================================================================
# The missing imputed alive in their three parts: the prisoners of war among
# them, set by the official figures; those projected to leave the register
# (found alive); and those alive for other reasons, none at the central
# values.
tab_imput <-
  imput |>
  mutate(
    pct_dead = imputed_dead / missing_stock,
    pct_alive = imputed_alive / missing_stock
  ) |>
  select(
    year,
    registered_deaths = registered,
    late_registrations,
    missing = missing_stock,
    imputed_dead,
    pct_dead,
    prisoners_of_war = captives,
    projected_no_longer_listed = imputed_unlisted,
    alive_other = alive_other,
    pct_alive,
    total_deaths = total_estimado
  )
stopifnot(isTRUE(all.equal(imput$imputed_dead + imput$imputed_alive, imput$missing_stock)),
          isTRUE(all.equal(imput$captives + imput$imputed_unlisted + imput$alive_other, imput$imputed_alive)))

# The Total row sums the rounded years, as tables 1 and 3 do, so the tables
# print the same totals; rounding the four-year sum once can differ by one.
count_cols <- c("registered_deaths", "late_registrations", "missing", "imputed_dead",
                "prisoners_of_war", "projected_no_longer_listed", "alive_other",
                "total_deaths")
tab_imput <- bind_rows(
  tab_imput,
  tab_imput |>
    summarise(across(all_of(count_cols), ~ sum(round(.x)))) |>
    mutate(
      year = NA_integer_,
      pct_dead = imputed_dead / missing,
      pct_alive = (prisoners_of_war + projected_no_longer_listed + alive_other) / missing
    )
) |>
  mutate(
    across(all_of(count_cols), ~round(.x)),
    across(c(pct_dead, pct_alive), ~scales::percent(.x, accuracy = 0.1)),
    year = if_else(is.na(year), "Total", as.character(year))
  )

# beside the imputation at the central values, the military total the
# simulation gives: its median and 95% interval, by year and in total
mil_sim <- param_draws |> filter(role == "combatants")
fmt_ui_count <- \(m, lo, hi) sprintf("%s (%s-%s)", scales::comma(round(m)), scales::comma(round(lo)), scales::comma(round(hi)))
mil_ui <- bind_rows(
  mil_sim |>
    summarise(m = median(draw), lo = quantile(draw, 0.025, names = FALSE), hi = quantile(draw, 0.975, names = FALSE),
              .by = year) |>
    mutate(year = as.character(year)),
  mil_sim |>
    summarise(t = sum(draw), .by = sim_id) |>
    summarise(m = median(t), lo = quantile(t, 0.025, names = FALSE), hi = quantile(t, 0.975, names = FALSE)) |>
    mutate(year = "Total")
) |>
  transmute(year, military_deaths_simulated = fmt_ui_count(m, lo, hi))
tab_imput <- tab_imput |> left_join(mil_ui, by = "year")
stopifnot(!anyNA(tab_imput$military_deaths_simulated))
save_tab(tab_imput, "table3_missing_imputation.csv")
print(tab_imput)

# ==============================================================================
# TABLE 3 - PERT simulation input bounds
# ==============================================================================
# The western component is bracketed by source, so its columns are the two
# readings' paths (crossings, register) either side of the mode, not per-year
# minima and maxima: a draw blends the paths, and the Total row is then the
# range of the four-year total.
tab_pert <-
  param_table |>
  select(year, role, mode, min, max, crossings, register) |>
  mutate(across(c(mode, min, max, crossings, register), ~round(.x))) |>
  pivot_wider(
    names_from = role,
    values_from = c(mode, min, max, crossings, register),
    names_glue = "{role}_{.value}"
  ) |>
  select(
    year,
    combatants_mode, combatants_min, combatants_max,
    civilians_mode, civilians_min, civilians_max,
    mig_west_mode, mig_west_crossings, mig_west_register,
    mig_ru_by_mode, mig_ru_by_min, mig_ru_by_max
  )

tab_pert <- bind_rows(
  tab_pert |> mutate(year = as.character(year)),
  # mig_ru_by has a 2022 row only, so the other years are empty, not zero
  tab_pert |> summarise(across(-year, \(x) sum(x, na.rm = TRUE))) |> mutate(year = "Total")
)

save_tab(tab_pert, "table2_pert_input_bounds.csv")
print(tab_pert)

# The inputs on the missing alive behind the combatant bounds, each drawn once
# per simulation from a Beta-PERT (alive_evidence(), 00_setup.R), with the two
# counts they give: the prisoners of war the register does not record as such,
# and those among its missing.
tab_alive <-
  tibble(
    input = c("Share of the unrecorded prisoners of war among the missing",
              "Prisoners of war held, February 2026",
              "Prisoners of war the register does not record (held + returned − recorded)",
              "Prisoners of war among the missing",
              "Share of the unresolved missing alive for other reasons"),
    min = with(evidence, c(s_min, held_min, unrecorded_min, captives_min, other_min)),
    mode = with(evidence, c(s_mode, held_mode, unrecorded_mode, captives_mode, other_mode)),
    central = with(evidence, c(s_central, held_central, unrecorded_central, captives_central, other_central)),
    max = with(evidence, c(s_max, held_max, unrecorded_max, captives_max, other_max)),
    basis = c(
      "Returned prisoners listed as missing before their return: every event year (min), 2024–2025 events (mode), all (max)",
      "\"About 7,000\" (President of Ukraine, 14 February 2026), rounded to the thousand",
      sprintf("Held plus %s military personnel returned, less the %s the register records as prisoners or released",
              scales::comma(evidence$returned_military), scales::comma(round(evidence$register_alive))),
      "Share among the missing times the unrecorded prisoners",
      "None (min and mode); the share of the register's first resolutions that leave it beyond list maintenance (max)"
    )
  ) |>
  mutate(across(c(min, mode, central, max), \(x) if_else(x < 1, round(x, 3), round(x))))
save_tab(tab_alive, "tableA26_missing_alive_inputs.csv")
print(tab_alive)

# ==============================================================================
# TABLE 4 - Conflict deaths by cause, sex and year  (NEW)
# ==============================================================================
# Rows are year x cause, sexes are columns, and every cell carries a 95%
# uncertainty interval.
#
# The Total column and the Total row are computed by summing WITHIN each draw
# and then taking quantiles of that sum - not by adding the two sexes' bounds
# together, which would overstate the interval.
d4_draws <-
  loss_draws |>
  as_tibble() |>
  select(
    sim_id, year, sex,
    Civilians = dx_civilian,
    `Registered combatants` = dx_cmb_registered,
    `Late registrations (estimated)` = dx_cmb_late,
    `Missing combatants (imputed)` = dx_cmb_imputed
  ) |>
  pivot_longer(
    -c(sim_id, year, sex),
    names_to = "role",
    values_to = "dx"
  )

fmt_ci <- function(m, lo, hi) {
  sprintf(
    "%s (%s-%s)",
    scales::comma(round(m), accuracy = 1),
    scales::comma(round(lo), accuracy = 1),
    scales::comma(round(hi), accuracy = 1)
  )
}

# one cell per year x cause x sex
cell_sex <-
  d4_draws |>
  summarise(
    m = median(dx), lo = quantile(dx, 0.025), hi = quantile(dx, 0.975),
    .by = c(year, role, sex)
  ) |>
  mutate(
    val = fmt_ci(m, lo, hi),
    sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males")
  ) |>
  select(year, role, sex, val) |>
  pivot_wider(names_from = sex, values_from = val)

# both sexes combined, per year x cause
cell_tot <-
  d4_draws |>
  summarise(dx = sum(dx), .by = c(sim_id, year, role)) |>
  summarise(
    m = median(dx), lo = quantile(dx, 0.025), hi = quantile(dx, 0.975),
    .by = c(year, role)
  ) |>
  mutate(Total = fmt_ci(m, lo, hi)) |>
  select(year, role, Total)

# bottom row: everything summed, per sex and overall
grand_sex <-
  d4_draws |>
  summarise(dx = sum(dx), .by = c(sim_id, sex)) |>
  summarise(
    m = median(dx), lo = quantile(dx, 0.025), hi = quantile(dx, 0.975),
    .by = sex
  ) |>
  mutate(
    val = fmt_ci(m, lo, hi),
    sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males")
  ) |>
  select(sex, val) |>
  pivot_wider(names_from = sex, values_from = val)

grand_tot <-
  d4_draws |>
  summarise(dx = sum(dx), .by = sim_id) |>
  summarise(
    m = median(dx), lo = quantile(dx, 0.025), hi = quantile(dx, 0.975)
  ) |>
  mutate(Total = fmt_ci(m, lo, hi)) |>
  select(Total)

tab_cause <-
  cell_sex |>
  left_join(cell_tot, by = c("year", "role")) |>
  mutate(role = factor(role, levels = three_way)) |>
  arrange(role, year) |>
  mutate(year = as.character(year), role = as.character(role)) |>
  bind_rows(bind_cols(tibble(year = "", role = "Total"), grand_sex, grand_tot)) |>
  select(year, role, Females, Males, Total)

save_tab(tab_cause, "table4_conflict_deaths_by_cause.csv")
print(tab_cause)

# ==============================================================================
# TABLES 6 and 7 - the two margins of table 4, for the main text
# ==============================================================================
# Table 4 is the full year x cause grid. These collapse it one way or the
# other: totals by cause (summed over years) and totals by year (summed over
# causes). Both are built from the draws in the same way, so the intervals
# stay correct.
totals_by <- function(gv) {
  by_sex <-
    d4_draws |>
    summarise(dx = sum(dx), .by = all_of(c("sim_id", gv, "sex"))) |>
    summarise(
      m = median(dx), lo = quantile(dx, 0.025), hi = quantile(dx, 0.975),
      .by = all_of(c(gv, "sex"))
    ) |>
    mutate(
      val = fmt_ci(m, lo, hi),
      sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males")
    ) |>
    select(all_of(c(gv, "sex")), val) |>
    pivot_wider(names_from = sex, values_from = val)

  both <-
    d4_draws |>
    summarise(dx = sum(dx), .by = all_of(c("sim_id", gv))) |>
    summarise(
      m = median(dx), lo = quantile(dx, 0.025), hi = quantile(dx, 0.975),
      .by = all_of(gv)
    ) |>
    mutate(Total = fmt_ci(m, lo, hi)) |>
    select(all_of(gv), Total)

  by_sex |> left_join(both, by = gv)
}

grand_row <- bind_cols(grand_sex, grand_tot)

tab_by_cause <-
  totals_by("role") |>
  mutate(role = factor(role, levels = three_way)) |>
  arrange(role) |>
  mutate(role = as.character(role)) |>
  bind_rows(bind_cols(tibble(role = "Total"), grand_row)) |>
  select(role, Females, Males, Total)

save_tab(tab_by_cause, "table5_totals_by_cause.csv")
print(tab_by_cause)

tab_by_year <-
  totals_by("year") |>
  arrange(year) |>
  mutate(year = as.character(year)) |>
  bind_rows(bind_cols(tibble(year = "Total"), grand_row)) |>
  select(year, Females, Males, Total)

save_tab(tab_by_year, "table6_totals_by_year.csv")
print(tab_by_year)

# ==============================================================================
# TABLE 5 - Life expectancy loss
# ==============================================================================
tab_e0 <-
  loss_draws |>
  as_tibble() |>
  summarise(
    e0_expected = median(e0_bsn),
    e0_observed = median(e0_war),
    loss_median = median(loss_total),
    loss_lo = quantile(loss_total, 0.025),
    loss_hi = quantile(loss_total, 0.975),
    .by = c(year, sex)
  ) |>
  arrange(year, sex) |>
  mutate(
    across(c(e0_expected, e0_observed), ~round(.x, 2)),
    loss = sprintf("%.2f (%.2f-%.2f)", loss_median, loss_lo, loss_hi)
  ) |>
  select(year, sex, e0_expected, e0_observed, loss)
# the loss with the counterfactual at its point forecast (14c): the
# uncertainty of the war's toll alone
tiers <- read_rds("data_inter/ukr_fixed_counterfactual_e0.rds")
tab_e0 <-
  tab_e0 |>
  left_join(tiers |>
              filter(counterfactual == "fixed at its point forecast") |>
              transmute(year, sex, loss_counterfactual_fixed = sprintf("%.2f (%.2f-%.2f)", median, lo, hi)),
            by = c("year", "sex"))
stopifnot(!anyNA(tab_e0$loss_counterfactual_fixed))

save_tab(tab_e0, "table7_life_expectancy_loss.csv")
print(tab_e0)

# ==============================================================================
# TABLE A3 - Life expectancy loss attributed to each cause
# ==============================================================================
# Table 4 reports death COUNTS by cause; this keeps the corresponding years of
# life expectancy lost, which the previous combined table carried and which is
# otherwise only visible as shares in figure 5. Drop it if the paper does not
# need it.
tab_e0_cause <-
  loss_sum |>
  filter(cause %in% three_way) |>
  mutate(
    loss = sprintf("%.3f (%.3f-%.3f)", loss_median, loss_lo, loss_hi),
    sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males")
  ) |>
  select(year, role = cause, sex, loss) |>
  pivot_wider(names_from = sex, values_from = loss) |>
  arrange(role, year)

save_tab(tab_e0_cause, "tableA3_e0_loss_by_cause.csv")

# ==============================================================================
# TABLE A1 - Source reconciliation, UCDP vs ACLED
# ==============================================================================
ucdp_cmp <- read_rds("data_inter/ukr_rus_ucdp_source_comparison.rds")
acled_cmp <- read_rds("data_inter/ukr_rus_acled_source_comparison.rds")

tA1 <-
  bind_rows(ucdp_cmp, acled_cmp) |>
  mutate(
    source = if_else(type == "adjusted", paste(source, "adjusted"), source)
  ) |>
  select(source, country, civilians, combatants, unknown) |>
  mutate(across(c(civilians, combatants, unknown), ~round(.x))) |>
  arrange(country, source)

save_tab(tA1, "tableA1_source_reconciliation.csv")
print(tA1)

# ==============================================================================
# TABLE A5 - Ukrainian military deaths: UCDP against the register and this study
# ==============================================================================
# UCDP is the only source in Table 1 with an independent military count and
# bounds, and its best estimate for Ukraine sits BELOW the register's
# individually named dead - so the gap between it and this study is not only a
# disagreement about the missing. This table sets the three side by side and
# splits the gap into three parts: UCDP short of the named dead; the deaths
# not yet registered (08b), which the named count will still grow to include;
# and the imputed dead among the missing, whom an event-based count does not
# see.
#
# Two columns test the obvious alternative explanations:
#   ukr_share   the Ukrainian share of the military deaths UCDP attributes to a
#               side. A falling share is what one-sided thinning of reporting
#               looks like: Ukraine does not disclose its military losses, while
#               Russian losses are widely claimed.
#   even_split  UCDP best with the unknown-side deaths split evenly between the
#               belligerents, instead of in proportion to each side's reported
#               losses as 07 does. Anchored to 07's own figure, so the column
#               differs from it only in that rule.
#
# UCDP deaths in the Russia-Ukraine state dyad are assigned by side, not by
# where the event happened, so Ukrainian losses in Kursk and Belgorod count.
ged <- read_rds("data_inter/ucdp_ged_events_slim.rds")

ucdp_ukr <-
  read_rds("data_inter/ukr_ucdp_invals.rds") |>
  filter(role == "combatants") |>
  select(year, ucdp = dts, ucdp_lo = dts_l, ucdp_hi = dts_u)

dyad <-
  ged |>
  filter(
    year %in% 2022:2025, type_of_violence == 1,
    str_detect(side_a, "Russia"), str_detect(side_b, "Ukraine")
  ) |>
  summarise(
    ukr = sum(b), rus = sum(a), civ = sum(c), unk = sum(u, na.rm = TRUE),
    .by = year
  ) |>
  mutate(
    ukr_prop = ukr + unk * ukr / (ukr + rus + civ),
    ukr_even = ukr + unk * (1 - civ / (civ + ukr + rus)) / 2
  )

cmb_draws <- param_draws |> filter(role == "combatants")
q3 <- function(x) tibble(med = median(x), lo = quantile(x, 0.025), hi = quantile(x, 0.975))

recon_year <-
  dyad |>
  left_join(ucdp_ukr, by = "year") |>
  left_join(
    imput |> select(year, registered, late = late_registrations),
    by = "year"
  ) |>
  left_join(cmb_draws |> reframe(q3(draw), .by = year), by = "year") |>
  mutate(even_split = ucdp + (ukr_even - ukr_prop)) |>
  arrange(year)

recon_total <-
  recon_year |>
  summarise(across(c(ukr, rus, ucdp, ucdp_lo, ucdp_hi, even_split, registered, late), sum)) |>
  bind_cols(cmb_draws |> summarise(draw = sum(draw), .by = sim_id) |> reframe(q3(draw)))

tA5 <-
  bind_rows(
    recon_year |> mutate(year = as.character(year)),
    recon_total |> mutate(year = "Total")
  ) |>
  transmute(
    year,
    ucdp_ukrainian_side_reported = round(ukr),
    ucdp_best = fmt_ci(ucdp, ucdp_lo, ucdp_hi),
    ucdp_unknowns_split_evenly = round(even_split),
    ukrainian_share_of_ucdp_military = scales::percent(ukr / (ukr + rus), accuracy = 0.1),
    ualosses_registered_dead = round(registered),
    this_study = fmt_ci(med, lo, hi),
    ucdp_as_pct_of_registered = scales::percent(ucdp / registered, accuracy = 0.1),
    gap_to_ucdp = round(med - ucdp),
    of_which_ucdp_below_registered = round(registered - ucdp),
    of_which_late_registrations = round(late),
    of_which_imputed_missing = round(med - registered - late)
  )

save_tab(tA5, "tableA5_military_reconciliation.csv")
print(tA5)

# ==============================================================================
# TABLE A6 / FIGURE A4 - Sensitivity to the missing-combatant imputation
# ==============================================================================
# The military death total is dominated by how many of the missing are
# alive. The evidence (09, alive_evidence() in 00_setup.R) sets the prisoners
# of war among them within a range, and allows a share of the rest to be alive
# for other reasons; 11 draws both inside that range. 13b sweeps the missing
# alive along one path, from the fewest alive the evidence allows to all but
# the projected deaths, re-running the imputation and the downstream
# projection with everything else at its central value, and this assembles the result
# into a manuscript table and figure with the range marked on it.
tA6 <-
  alive_mil |>
  # summed over rounded years, as table 3 builds its totals, so the two agree
  summarise(total_military_deaths = sum(round(total_military)),
            .by = c(point, alive, share_alive, captives, other)) |>
  left_join(
    alive_e0 |>
      summarise(
        male_loss_2025 = loss[year == 2025 & sex == "m"],
        female_loss_2025 = loss[year == 2025 & sex == "f"],
        .by = alive
      ),
    by = "alive"
  ) |>
  arrange(alive) |>
  transmute(
    share_of_missing_alive = round(share_alive, 3),
    evidence = point,
    missing_alive = round(alive),
    prisoners_of_war_among_missing = round(captives),
    share_of_unresolved_alive_otherwise = round(other, 3),
    total_military_deaths,
    male_e0_loss_2025 = round(male_loss_2025, 2),
    female_e0_loss_2025 = round(female_loss_2025, 3)
  )

save_tab(tA6, "tableA6_missing_alive_sensitivity.csv")
print(tA6)

# "central (-down / +up)": the value at the central input and how far it moves
# at the two ends of the 95% interval of the simulated input. Used by A4 and A5.
pm_label <- function(central, a, b, digits) {
  lo <- pmin(a, b)
  hi <- pmax(a, b)
  negligible <- pmax(central - lo, hi - central) < 0.5 * 10^-digits
  if_else(
    negligible,
    sprintf("%.*f (±<%s)", digits, central, format(10^-digits, scientific = FALSE)),
    sprintf("%.*f (−%.*f / +%.*f)", digits, central, digits, central - lo, digits, hi - central)
  )
}

# text beside one end of each line, spread apart where lines end close
# together; the x axis is widened on that side to make room for it
end_labels <- function(d, side = c("right", "left")) {
  side <- match.arg(side)
  ggrepel::geom_text_repel(
    data = d, aes(x = x, y = y, label = label, colour = colour_key),
    inherit.aes = FALSE, hjust = if (side == "right") 0 else 1,
    direction = "y", nudge_x = 0,
    # labels move only vertically, so each keeps its own x in every panel
    xlim = c(-Inf, Inf),
    size = 2.7, segment.colour = NA,
    box.padding = 0.15, min.segment.length = Inf, show.legend = FALSE,
    # seeded, so the labels land in the same place in every run
    seed = 1
  )
}
# breaks only over the data when the axis is widened for labels: ggplot
# computes breaks on the widened range, which would put ticks under the text
data_breaks <- function(left, right) {
  function(lims) {
    r <- diff(lims) / (1 + left + right)
    lo <- lims[1] + left * r
    hi <- lo + r
    b <- scales::breaks_pretty()(c(lo, hi))
    b[b >= lo - 1e-9 & b <= hi + 1e-9]
  }
}

P_LO <- "2.5th percentile of draws"
P_HI <- "97.5th percentile of draws"
share_at <- function(p) alive_mil$share_alive[which(alive_mil$point == p)[1]]

# the range as a light band, the 95% interval of the draws as a darker one
# inside it, the central values as a dashed line
alive_band <- list(
  annotate(
    "rect",
    xmin = share_at("floor"), xmax = share_at("cap"),
    ymin = -Inf, ymax = Inf, fill = "#0A9396", alpha = 0.14
  ),
  annotate(
    "rect",
    xmin = share_at(P_LO), xmax = share_at(P_HI),
    ymin = -Inf, ymax = Inf, fill = "#0A9396", alpha = 0.18
  ),
  geom_vline(xintercept = share_at("central"), linetype = "dashed", colour = "grey30")
)
alive_x <- scale_x_continuous(
  "Share of the missing who are alive",
  labels = scales::percent, breaks = seq(0, 0.8, 0.2),
  expand = expansion(mult = c(0.85, 0.02))
)

mil_curve <-
  alive_mil |>
  summarise(total_military = sum(round(total_military)), .by = c(point, share_alive))
mil_at <- function(p) mil_curve$total_military[which(mil_curve$point == p)[1]]
mil_lo <- min(mil_at(P_LO), mil_at(P_HI))
mil_hi <- max(mil_at(P_LO), mil_at(P_HI))
# labels at the floor, the left end, where the lines are furthest apart
lab4a <- tibble(
  x = share_at("floor"),
  y = mil_at("floor"),
  colour_key = "military",
  label = sprintf(
    "%s (−%s / +%s)", scales::comma(mil_at("central")),
    scales::comma(mil_at("central") - mil_lo), scales::comma(mil_hi - mil_at("central"))
  )
)

figA4a <-
  mil_curve |>
  ggplot(aes(share_alive, total_military)) +
  alive_band +
  geom_line(linewidth = 1, colour = "#B23A48") +
  end_labels(lab4a, "left") +
  scale_colour_manual(values = c(military = "#B23A48"), guide = "none") +
  alive_x +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Total military deaths, 2022–2025") +
  theme_paper()

lab4b <-
  alive_e0 |>
  filter(point %in% c("central", P_LO, P_HI)) |>
  select(year, sex, point, loss) |>
  pivot_wider(names_from = point, values_from = loss) |>
  left_join(alive_e0 |> filter(point %in% "floor") |> select(year, sex, y = loss), by = c("year", "sex")) |>
  mutate(
    x = share_at("floor"),
    colour_key = factor(year),
    label = pm_label(central, .data[[P_LO]], .data[[P_HI]], if_else(sex == "m", 2, 3))
  ) |>
  nice_sex()

figA4b <-
  alive_e0 |>
  nice_sex() |>
  ggplot(aes(share_alive, loss, colour = factor(year))) +
  alive_band +
  geom_line(linewidth = 1) +
  end_labels(lab4b, "left") +
  facet_wrap(~sex, scales = "free_y") +
  alive_x +
  labs(y = "Life expectancy loss (years)", colour = "Year") +
  theme_paper()

figA4 <- figA4a + figA4b +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "Sensitivity to the share of the missing who are alive",
    caption = sprintf(
      paste0(
        "Light band: the range the evidence allows (%s to %s of the missing alive), within which ",
        "every simulation draws the evidence. Dark band: the 95%% interval of the draws\n",
        "(%s to %s). Dashed: the central values, %s: the %s prisoners of war the official figures ",
        "place among the missing, and those projected to leave the register.\n",
        "Labels: the value at the central values, and how far it moves at the two ends of ",
        "the dark band. Deterministic, at the central value of every other input - not a Monte ",
        "Carlo re-run at each point."
      ),
      scales::percent(share_at("floor"), 0.1), scales::percent(share_at("cap"), 0.1),
      scales::percent(share_at(P_LO), 0.1), scales::percent(share_at(P_HI), 0.1),
      scales::percent(share_at("central"), 0.1), scales::comma(round(evidence$captives_mode))
    )
  )
save_fig(figA4, "figA4_missing_alive_sensitivity.png", 14, 5)

# ==============================================================================
# TABLE A7 / FIGURE A5 - Sensitivity to net migration
# ==============================================================================
# The counterpart of A4 for migration. 13c sweeps each migration component
# across its range, along the quantile 11 draws, with the other component
# and every conflict input at their mode. The band is the 95% interval of the
# component's simulated four-year total, the dashed line its mode. Figure 4
# gives the other end of the question: the loss with no migration at all.
mig_sens <- read_rds("data_inter/ukr_migration_sensitivity_e0.rds")
mig_decomp_2025 <- mig |> filter(year == 2025)

mig_labels <- c(mig_west = "Western", mig_ru_by = "Russia and Belarus")
point_order <- c("min", "p2.5", "mode", "p97.5", "max")

tA7 <-
  bind_rows(
    mig_sens |>
      filter(!is.na(point), year == 2025) |>
      select(component, point, cumulative_net_migration, sex, loss),
    # no migration at all (13), the reference for both components
    tibble(component = "none", point = "no migration", cumulative_net_migration = 0) |>
      cross_join(mig_decomp_2025 |> select(sex, loss = loss_nomig))
  ) |>
  pivot_wider(names_from = sex, values_from = loss) |>
  mutate(
    component = factor(component, levels = c("mig_west", "mig_ru_by", "none"),
                       labels = c("Western", "Russia and Belarus", "None")),
    point = factor(point, levels = c(point_order, "no migration"))
  ) |>
  arrange(component, point) |>
  transmute(
    component,
    point,
    cumulative_net_migration_2022_2025 = round(cumulative_net_migration),
    male_e0_loss_2025 = round(m, 2),
    female_e0_loss_2025 = round(f, 3)
  )

save_tab(tA7, "tableA7_migration_sensitivity.csv")
print(tA7)

# one strip per panel, "Western · Males", rather than two stacked strips
mig_panels <- as.vector(t(outer(mig_labels, c("Females", "Males"), paste, sep = " · ")))
mig_panel <- function(component, sex) {
  factor(paste(mig_labels[component], sex, sep = " · "), levels = mig_panels)
}

mig_curves <-
  mig_sens |>
  filter(!is.na(u)) |>
  nice_sex() |>
  mutate(panel = mig_panel(component, sex))

mig_marks <-
  mig_sens |>
  filter(point %in% c("p2.5", "mode", "p97.5")) |>
  distinct(component, point, cumulative_net_migration) |>
  pivot_wider(names_from = point, values_from = cumulative_net_migration) |>
  cross_join(tibble(sex = c("Females", "Males"))) |>
  mutate(panel = mig_panel(component, sex))

# labels beside the right-hand end of each line: the loss at the mode and how
# far it moves at the two ends of the band (pm_label() and end_labels() are
# defined with A4)
lab5 <-
  mig_sens |>
  filter(point %in% c("mode", "p2.5", "p97.5")) |>
  select(component, year, sex, point, loss) |>
  pivot_wider(names_from = point, values_from = loss) |>
  left_join(
    mig_sens |> filter(point == "max") |> select(component, year, sex, x = cumulative_net_migration, y = loss),
    by = c("component", "year", "sex")
  ) |>
  mutate(
    x = x / 1e6,
    colour_key = factor(year),
    label = pm_label(mode, p2.5, p97.5, if_else(sex == "m", 2, 3))
  ) |>
  nice_sex() |>
  mutate(panel = mig_panel(component, sex))

figA5 <-
  mig_curves |>
  ggplot(aes(cumulative_net_migration / 1e6, loss, colour = factor(year))) +
  geom_rect(
    data = mig_marks, inherit.aes = FALSE,
    aes(xmin = p2.5 / 1e6, xmax = p97.5 / 1e6, ymin = -Inf, ymax = Inf),
    fill = "#0A9396", alpha = 0.18
  ) +
  geom_vline(data = mig_marks, aes(xintercept = mode / 1e6),
             linetype = "dashed", colour = "grey30") +
  geom_line(linewidth = 1) +
  end_labels(lab5) +
  facet_wrap(~panel, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = data_breaks(0.02, 0.7), expand = expansion(mult = c(0.02, 0.7))) +
  labs(
    x = "Net migration of the component, cumulative 2022–2025 (millions, net outflow)",
    y = "Life expectancy loss (years)",
    colour = "Year",
    title = "Sensitivity to net migration",
    caption = sprintf(
      paste0(
        "Each component is swept across its range with the other and every ",
        "conflict input at their mode. Shaded: the 95%% interval of the\n",
        "component's simulated total. Dashed: its mode. Labels: the loss at the mode, ",
        "and how far it moves at the two ends of the shaded interval.\n",
        "With no migration at all the 2025 loss would be %.2f years for men and ",
        "%.2f for women (figure 4)."
      ),
      mig_decomp_2025$loss_nomig[mig_decomp_2025$sex == "m"],
      mig_decomp_2025$loss_nomig[mig_decomp_2025$sex == "f"]
    )
  ) +
  theme_paper()
save_fig(figA5, "figA5_migration_sensitivity.png", 12, 7)

# ==============================================================================
# TABLE A8 - Migration specification checks
# ==============================================================================
# Four assumptions of the migration input changed one at a time (13c,
# section 5): where the western mode sits in the bracket (A2), the age-sex
# profile of the Russia and Belarus flow (A9), additive versus proportional
# allocation of a draw away from the mode (A5, at both ends of the western
# bracket), and the Canada and USA interpolation 30% lower and higher (A10).
# Loss in years, every conflict input at its mode.
spec <- read_rds("data_inter/ukr_migration_specification_e0.rds")
spec_order <- unique(spec$scenario)

tA8 <-
  spec |>
  filter(sex == "m" | year == 2025) |>
  mutate(col = if_else(sex == "m", paste0("male_", year), paste0("female_", year))) |>
  select(scenario, col, loss) |>
  pivot_wider(names_from = col, values_from = loss) |>
  mutate(
    scenario = factor(scenario, levels = spec_order),
    across(starts_with("male_"), \(x) round(x, 2)),
    across(starts_with("female_"), \(x) round(x, 3))
  ) |>
  arrange(scenario)

save_tab(tA8, "tableA8_migration_specification.csv")
print(tA8)

# ==============================================================================
# TABLES A9-A15 - The remaining sensitivity analyses and summary measures
# ==============================================================================
# One row per scenario, the loss in years at the central value of every other input:
# men in each year and women in 2025, the cells that move. A13 is the PERT
# shape, and A14 and A15 the summary measures from the draws, with intervals.
loss_wide <- function(d, keep_cols) {
  d |>
    filter(sex == "m" | year == 2025) |>
    mutate(col = if_else(sex == "m", paste0("male_", year), paste0("female_", year))) |>
    select(all_of(keep_cols), col, loss) |>
    pivot_wider(names_from = col, values_from = loss) |>
    mutate(across(starts_with("male_"), \(x) round(x, 2)),
           across(starts_with("female_"), \(x) round(x, 3)))
}

# A9: the linkage rules and the design of the chain (09, 13b)
tA9 <-
  read_rds("data_inter/ukr_linkage_rules_e0.rds") |>
  mutate(prisoners_of_war_among_missing = round(captives), military_deaths = round(military)) |>
  loss_wide(c("design", "prisoners_of_war_among_missing", "military_deaths"))
save_tab(tA9, "tableA9_linkage_and_chain.csv")
print(tA9)

# A10: the horizon of the registration-lag correction (08b, 13d)
lag <- read_rds("data_inter/ukr_registration_lag.rds")
tA10 <-
  lag$loss |>
  left_join(lag$military |> summarise(military_deaths = sum(round(total)), .by = scenario),
            by = "scenario") |>
  mutate(scenario = factor(scenario, levels = unique(lag$military$scenario))) |>
  loss_wide(c("scenario", "military_deaths")) |>
  arrange(scenario)
save_tab(tA10, "tableA10_registration_lag.csv")
print(tA10)

# A11: the population base of Donetsk and Luhansk (13e) and the timing of
# 2022's events (13f)
dl <- read_rds("data_inter/ukr_denominator_donetsk_luhansk.rds")
timing <- read_rds("data_inter/ukr_timing_2022.rds")
tA11 <-
  bind_rows(
    dl$loss |> mutate(check = "Population base of Donetsk and Luhansk"),
    dl$migrants |> mutate(check = "Population base outside Donetsk and Luhansk"),
    timing$loss |> mutate(check = "Timing of 2022's deaths and departures")
  ) |>
  mutate(scenario = factor(scenario, levels = unique(scenario))) |>
  loss_wide(c("check", "scenario")) |>
  arrange(check, scenario)
save_tab(tA11, "tableA11_population_base_and_timing.csv")
print(tA11)

# A12: the Lee-Carter window (13g), with the counterfactual e0 in 2025
# the windows (13g), with the coherent forecast (13j) before the rows one SD
# either way
baseline <- read_rds("data_inter/ukr_baseline_window_e0.rds")
coherent <- read_rds("data_inter/ukr_coherent_forecast_e0.rds")
sd_rows <- str_detect(baseline$window, "one SD")
baseline <- bind_rows(baseline[!sd_rows, ], coherent$baseline, baseline[sd_rows, ])
tA12 <-
  baseline |>
  mutate(window = factor(window, levels = unique(window))) |>
  loss_wide("window") |>
  left_join(baseline |> filter(year == 2025) |>
              transmute(window = factor(window, levels = unique(baseline$window)), sex,
                        e0 = round(e0_bsn, 2)) |>
              pivot_wider(names_from = sex, values_from = e0, names_prefix = "counterfactual_e0_2025_"),
            by = "window") |>
  arrange(window)
spread <-
  bind_rows(read_rds("data_inter/ukr_forecast_spread_by_window.rds"), coherent$spread) |>
  filter(year == 2025, sex == "m")
tA12 <-
  tA12 |>
  left_join(spread |>
              summarise(band_counterfactual_e0_2025_m = sprintf("%.2f-%.2f", min(e0_bsn), max(e0_bsn)),
                        band_male_loss_2025 = sprintf("%.2f-%.2f", min(loss), max(loss)),
                        .by = window) |>
              mutate(window = factor(window, levels = levels(tA12$window))),
            by = "window")
save_tab(tA12, "tableA12_counterfactual_window.csv")
print(tA12)

# A13: the PERT shape (13g): median and 95% interval of the loss
tA13 <-
  read_rds("data_inter/ukr_pert_shape_e0.rds") |>
  mutate(sex = if_else(sex == "f", "Females", "Males"),
         loss = sprintf("%.2f (%.2f-%.2f)", median, lo, hi)) |>
  select(year, sex, shape, loss) |>
  pivot_wider(names_from = shape, values_from = loss, names_prefix = "pert_shape_") |>
  arrange(sex, year)
save_tab(tA13, "tableA13_pert_shape.csv")
print(tA13)

# A25: the age profile of civilian deaths (13h). The deaths beyond OHCHR's
# verified count given an older profile, with the share aged 60 or more and
# the mean age of all civilian deaths in 2022, and the loss in 2022, when
# nearly all of those deaths fall, and in 2025
civ_age <- read_rds("data_inter/ukr_civilian_age_profile_e0.rds")
tA25 <-
  civ_age$loss |>
  filter(year %in% c(2022, 2025)) |>
  mutate(col = paste0(if_else(sex == "m", "male_", "female_"), year),
         loss = round(loss, if_else(sex == "f" & year == 2025, 3, 2))) |>
  select(scenario, col, loss) |>
  pivot_wider(names_from = col, values_from = loss) |>
  left_join(civ_age$stats |>
              filter(year == 2022) |>
              transmute(scenario, civilian_60_plus_2022 = round(share_60_plus, 3),
                        civilian_mean_age_2022 = round(mean_age, 1)),
            by = "scenario") |>
  mutate(scenario = factor(scenario, levels = unique(civ_age$loss$scenario))) |>
  arrange(scenario) |>
  select(scenario, civilian_60_plus_2022, civilian_mean_age_2022, female_2022, male_2022, female_2025, male_2025)
save_tab(tA25, "tableA25_civilian_age_profile.csv")
print(tA25)

# A27: residents of occupied Donbas killed in Russian-controlled forces (13i)
donbas <- read_rds("data_inter/ukr_donbas_fighters_e0.rds")
tA27 <-
  donbas$loss |>
  filter(year %in% c(2022, 2023, 2025)) |>
  mutate(col = paste0(if_else(sex == "m", "male_", "female_"), year),
         loss = round(loss, if_else(sex == "f", 3, 2))) |>
  select(scenario, col, loss) |>
  pivot_wider(names_from = col, values_from = loss) |>
  left_join(donbas$deaths |> summarise(deaths_added = round(sum(deaths)), .by = scenario), by = "scenario") |>
  mutate(deaths_added = coalesce(deaths_added, 0),
         scenario = factor(scenario, levels = unique(donbas$loss$scenario))) |>
  arrange(scenario) |>
  select(scenario, deaths_added, male_2022, male_2023, male_2025, female_2025)
save_tab(tA27, "tableA27_donbas_fighters.csv")
print(tA27)

# A28: the resolution model against the releases held out of the fit (09i).
# Written only where 09i has run: it needs the two releases the estimates do
# not use, v15 and v17.
oos_file <- "data_inter/ukr_ualosses_out_of_sample.rds"
if (file.exists(oos_file)) {
  oos <- read_rds(oos_file)
  nxt <- set_names(oos$releases$release[-1], head(oos$releases$release, -1))
  tA28 <-
    oos$by_window |>
    transmute(window = paste0(from_release, " to ", nxt[from_release]),
              boundary = if_else(held_out, "held out of the fit", "fitted"),
              missing_at_start = round(at_risk),
              dead_observed = round(dead), dead_predicted = round(pred_dead),
              resolved_observed = round(resolved), resolved_predicted = round(pred_resolved),
              predicted_over_observed = round(ratio, 2))
  save_tab(tA28, "tableA28_out_of_sample.csv")
  print(tA28)
} else {
  message("  no out-of-sample test (09i has not run): table A28 is not written")
}

# A14: years of life lost and 45q15 (14b), medians and 95% intervals
yq <- read_rds("data_inter/ukr_yll_45q15_summary.rds")
fmt_ui <- function(m, lo, hi, digits = 0) {
  f <- if (digits == 0) \(x) scales::comma(round(x)) else \(x) formatC(x, format = "f", digits = digits)
  sprintf("%s (%s-%s)", f(m), f(lo), f(hi))
}
yll_labels <- c(yll_civilian = "Civilians", yll_registered = "Registered combatants",
                yll_late = "Late registrations (estimated)", yll_imputed = "Missing combatants (imputed)",
                yll_total = "Total")
tA14_yll <-
  bind_rows(yq$both_sexes |> mutate(year = as.character(year)),
            yq$all_years_both |> mutate(year = "2022–2025")) |>
  filter(measure %in% names(yll_labels)) |>
  mutate(cause = factor(yll_labels[as.character(measure)], levels = yll_labels),
         yll = fmt_ui(median, lo, hi)) |>
  select(year, cause, yll) |>
  pivot_wider(names_from = year, values_from = yll) |>
  arrange(cause)
q_labels <- c(q4515_bsn = "Counterfactual", q4515_war = "With conflict deaths",
              q4515_diff = "Difference")
tA14_q <-
  yq$by_year_sex |>
  filter(measure %in% names(q_labels)) |>
  mutate(sex = if_else(sex == "f", "Females", "Males"),
         measure = factor(q_labels[as.character(measure)], levels = q_labels),
         q = fmt_ui(1000 * median, 1000 * lo, 1000 * hi, digits = 1)) |>
  select(sex, measure, year, q) |>
  pivot_wider(names_from = year, values_from = q) |>
  arrange(sex, measure)
save_tab(tA14_yll, "tableA14_years_of_life_lost.csv")
save_tab(tA14_q, "tableA15_adult_mortality_45q15.csv")
print(tA14_yll)
print(tA14_q)
# ==============================================================================
# TABLE A2 - Resolution of the missing over twelve months
# ==============================================================================
# From 09: everyone listed as missing in v14, v16 or v18, followed to v19, the
# three windows between the four releases composed into the twelve months
# from v14 to v19. A return from captivity is a resolution to captivity; the
# people whose first resolution it was are counted apart. Shares of the
# cohort's missing.
res12 <- read_rds("data_inter/ukr_ualosses_resolution_12m.rds")
returned_first <-
  read_rds("data_inter/ukr_ualosses_linkage_checks.rds")$released_first |>
  summarise(returned_from_captivity = sum(n), .by = year)
tA2 <-
  res12 |>
  left_join(returned_first, by = "year") |>
  transmute(
    year,
    missing_v14 = at_risk_v14,
    first_listed_v16_v18 = later_entrants,
    returned_from_captivity,
    across(c(dead, prisoner, no_longer_listed, still_missing = missing),
           \(x) scales::percent(x, accuracy = 0.1))
  ) |>
  arrange(year)

save_tab(tA2, "tableA2_status_transitions.csv")
print(tA2)

# ==============================================================================
# TABLE 8 - How far the results move with each structural choice
# ==============================================================================
# The simulation's intervals cover the drawn inputs only. Each row here is a
# set of deterministic projections, every input not under study at its central
# value (the mode, or the median for the inputs on the missing alive),
# and gives the range its alternatives span: the military total, and the male
# and female loss in 2022 and 2025. The rows behind each range are
# in the supplementary table named.
# Military totals are summed over rounded years, as tables 3 and A6 build theirs.
mil_total <- function(d) d |> summarise(military = sum(round(total_military)), .by = c(point, alive))
mode_military <- alive_mil |> filter(point %in% "central") |> summarise(m = sum(round(total_military))) |> pull(m)
alive_rows <- function(points) {
  alive_e0 |>
    filter(point %in% points) |>
    left_join(mil_total(alive_mil), by = c("point", "alive")) |>
    select(military, year, sex, loss)
}
mig_sens_t8 <- read_rds("data_inter/ukr_migration_sensitivity_e0.rds")
spec_t8 <- read_rds("data_inter/ukr_migration_specification_e0.rds")
lag_t8 <- read_rds("data_inter/ukr_registration_lag.rds")
alternatives <- bind_rows(
  alive_rows(c("2.5th percentile of draws", "97.5th percentile of draws")) |>
    mutate(analysis = "The missing alive: 95% interval of the draws", table = "A6"),
  alive_rows(c("floor", "cap")) |>
    mutate(analysis = "The missing alive: floor to cap", table = "A6"),
  lag_t8$loss |>
    filter(scenario != "to 48 months (used)") |>
    left_join(lag_t8$military |> summarise(military = sum(round(total)), .by = scenario), by = "scenario") |>
    select(military, year, sex, loss) |>
    mutate(analysis = "Registration lag: no correction, or corrected to 60 or 72 months", table = "A10, A16"),
  read_rds("data_inter/ukr_linkage_rules_e0.rds") |>
    filter(!str_detect(design, "production")) |>
    mutate(analysis = paste0("Linkage rules, resolution model and prisoner-of-war evidence: ",
                             c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
                               "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen",
                               "seventeen", "eighteen", "nineteen", "twenty")[n_distinct(design)], " alternatives"),
           table = "A9") |>
    select(military, year, sex, loss, analysis, table),
  read_rds("data_inter/ukr_denominator_donetsk_luhansk.rds")$loss |>
    filter(scenario != "SSSU (used)") |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Population base of Donetsk and Luhansk: 0.5 million lower to taken out",
              table = "A11"),
  read_rds("data_inter/ukr_timing_2022.rds")$loss |>
    filter(scenario != "months of the events (used)") |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Timing within 2022: mid-year convention, or net outflow by 1 April", table = "A11"),
  read_rds("data_inter/ukr_denominator_donetsk_luhansk.rds")$migrants |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Population base outside Donetsk and Luhansk: 0.5 or 1.0 million fewer aged 20-64",
              table = "A11"),
  donbas$loss |>
    filter(scenario != "None (as estimated)") |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Residents of occupied Donbas killed in Russian-controlled forces: counted to September 2023, or continued to 2025",
              table = "A27"),
  civ_age$loss |>
    filter(scenario != "OHCHR's profile (used)") |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Civilian age profile: deaths beyond OHCHR's verified count with half over 60, or aged as pre-war deaths",
              table = "A25"),
  bind_rows(read_rds("data_inter/ukr_baseline_window_e0.rds"),
            read_rds("data_inter/ukr_coherent_forecast_e0.rds")$baseline) |>
    filter(window != "2000-2019 (used)") |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Counterfactual: Lee–Carter window, a coherent forecast with eight neighbours, and the forecast one SD lower or higher",
              table = "A12"),
  mig_sens_t8 |>
    filter(point %in% c("min", "max")) |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Net migration: each component across its range", table = "A7"),
  spec_t8 |>
    filter(scenario != "Mode") |>
    transmute(military = mode_military, year, sex, loss,
              analysis = "Specification of the net migration input", table = "A8"),
  mig |>
    transmute(military = mode_military, year, sex, loss = loss_nomig,
              analysis = "No net migration", table = "A7")
)
fmt_range <- function(x, f) {
  lo <- f(min(x)); hi <- f(max(x))
  if (lo == hi) lo else paste0(lo, " - ", hi)
}
f_mil <- \(x) scales::comma(round(x))
f_m <- \(x) formatC(x, format = "f", digits = 2)
f_f <- \(x) formatC(x, format = "f", digits = 3)
cells_t8 <- function(d) {
  tibble(
    military_deaths = fmt_range(d$military, f_mil),
    male_loss_2022 = fmt_range(d$loss[d$sex == "m" & d$year == 2022], f_m),
    female_loss_2022 = fmt_range(d$loss[d$sex == "f" & d$year == 2022], f_m),
    male_loss_2025 = fmt_range(d$loss[d$sex == "m" & d$year == 2025], f_m),
    female_loss_2025 = fmt_range(d$loss[d$sex == "f" & d$year == 2025], f_f)
  )
}
mode_rows <- alive_e0 |> filter(point %in% "central") |> mutate(military = mode_military)
tab8 <-
  bind_rows(
    cells_t8(mode_rows) |> mutate(analysis = "Every input at its central value", table = "", .before = 1),
    alternatives |>
      mutate(analysis = factor(analysis, levels = unique(analysis))) |>
      nest(.by = c(analysis, table)) |>
      mutate(cells = map(data, cells_t8)) |>
      select(-data) |>
      unnest(cells) |>
      mutate(analysis = as.character(analysis))
  ) |>
  rename(supplementary_table = table)
save_tab(tab8, "table8_structural_sensitivity.csv")
print(tab8)

# ==============================================================================
# TABLE A16 - The missing alive and the registration-lag horizon together
# ==============================================================================
tA16 <-
  lag_t8$grid |>
  filter(year == 2025, sex == "m") |>
  transmute(evidence = factor(point, levels = unique(lag_t8$grid$point)),
            share_of_missing_alive = round(share_alive, 3),
            scenario = factor(scenario, levels = unique(lag_t8$grid$scenario)),
            cell = sprintf("%s / %.2f", scales::comma(round(military)), loss)) |>
  arrange(scenario) |>
  pivot_wider(names_from = scenario, values_from = cell) |>
  arrange(evidence)
save_tab(tA16, "tableA16_missing_alive_by_lag_horizon.csv")
print(tA16)

# ==============================================================================
# TABLE A17 - The military estimate against official statements and other
# estimates, at the dates they were made
# ==============================================================================
# Each statement counts the deaths known at its date. The register (v19, as
# registered) and this study (at the mode) count the deaths among events up
# to that date, as the register knew them in July 2026. As in every other table, complete
# years enter as their rounded totals, so a date at the end of 2025 or later
# gives the four-year total of Tables 3 and 3.
official <- read_csv("data_input/official_figures.csv", show_col_types = FALSE)
mil_month <- read_rds("data_inter/ukr_military_by_month.rds")
through <- function(d, col) {
  d <- as.Date(d)
  years_done <- mil_month |>
    filter(year(month) < year(d)) |>
    summarise(n = sum(.data[[col]]), .by = year)
  this_year <- mil_month |> filter(year(month) == year(d))
  part <- sum(this_year[[col]][this_year$month < floor_date(d, "month")]) +
    sum(this_year[[col]][this_year$month == floor_date(d, "month")]) * day(d) / days_in_month(d)
  sum(round(years_done$n)) + round(part)
}
# the same in every draw: the draw's complete years, and its current year in
# the proportion of that year's deaths the months to the date hold at the
# central values
mil_year_draws <- param_draws |> filter(role == "combatants") |> select(sim_id, year, draw)
through_draws <- function(d) {
  d <- as.Date(d)
  this_year <- mil_month |> filter(year(month) == year(d))
  share <- if (nrow(this_year) == 0) 0 else
    (sum(this_year$military[this_year$month < floor_date(d, "month")]) +
       sum(this_year$military[this_year$month == floor_date(d, "month")]) * day(d) / days_in_month(d)) /
    sum(this_year$military)
  mil_year_draws |>
    summarise(v = sum(draw[year < year(d)]) + share * sum(draw[year == year(d)]), .by = sim_id) |>
    pull(v)
}
tA17 <-
  official |>
  filter(key %in% c("military_killed_official", "military_killed_estimate")) |>
  mutate(
    date = as.Date(date),
    kind = if_else(key == "military_killed_official", "Official Ukrainian statement", "Other estimate"),
    stated = if_else(is.na(value), sprintf("%s–%s", scales::comma(low), scales::comma(high)),
                     scales::comma(value)),
    register_named_dead = round(map_dbl(date, \(d) through(d, "registered_dead"))),
    this_study = round(map_dbl(date, \(d) through(d, "military"))),
    ui = map(date, \(d) {
      v <- through_draws(d)
      tibble(this_study_median = round(median(v)), this_study_lo = round(quantile(v, 0.025, names = FALSE)),
             this_study_hi = round(quantile(v, 0.975, names = FALSE)))
    })
  ) |>
  unnest(ui) |>
  arrange(date) |>
  select(date, kind, source, stated, register_named_dead, this_study, this_study_median, this_study_lo, this_study_hi)
save_tab(tA17, "tableA17_military_triangulation.csv")
print(tA17)

# ==============================================================================
# TABLE A18 - Where the returned prisoners had been listed before their return
# ==============================================================================
linkage_checks <- read_rds("data_inter/ukr_ualosses_linkage_checks.rds")
tA18 <-
  linkage_checks$released_prior |>
  mutate(prior = factor(prior, levels = c("prisoner", "missing", "dead", "not listed"),
                        labels = c("as prisoner", "as missing", "as dead", "in no earlier release"))) |>
  pivot_wider(names_from = prior, values_from = n, values_fill = 0) |>
  arrange(year) |>
  mutate(year = as.character(year))
tA18 <-
  bind_rows(tA18, tA18 |> summarise(across(-year, sum)) |> mutate(year = "Total")) |>
  mutate(total = rowSums(across(-year)),
         listed_as_missing_of_those_not_recorded_as_prisoners =
           scales::percent(`as missing` / (`as missing` + `in no earlier release`), accuracy = 0.1))
save_tab(tA18, "tableA18_returned_prisoners_prior_status.csv")
print(tA18)

# ==============================================================================
# TABLE A19 - How often the missing and the dead leave the register
# ==============================================================================
# Window by window, the share of those at risk at the window's start who are
# in no later release: the missing no longer listed, and the dead, who cannot
# have been found alive and leave only through list maintenance.
windows_cache <- read_rds("data_inter/ualosses_window_transitions.rds")
missing_dropout <-
  windows_cache |>
  filter(table == "windows", rule == "production") |>
  summarise(missing_at_risk = sum(n), missing_no_longer_listed = sum(n[to == "no_longer_listed"]),
            .by = c(year, from_release))
dead_dropout_all <-
  windows_cache |>
  filter(table == "dead_windows") |>
  summarise(dead_at_risk = sum(n), dead_dropped = sum(n[to == "no_longer_listed"]), .by = c(year, from_release))
tA19 <-
  missing_dropout |>
  left_join(dead_dropout_all, by = c("year", "from_release")) |>
  mutate(window = paste0(from_release, "-", ual_releases$release[match(from_release, ual_releases$release) + 1]),
         missing_pct = round(100 * missing_no_longer_listed / missing_at_risk, 2),
         dead_pct = round(100 * dead_dropped / dead_at_risk, 2)) |>
  arrange(year, match(from_release, ual_releases$release)) |>
  select(cohort = year, window, missing_at_risk, missing_no_longer_listed, missing_pct,
         dead_at_risk, dead_dropped, dead_pct)
stopifnot(nrow(tA19) > 0, !anyNA(tA19))
save_tab(tA19, "tableA19_register_dropout.csv")
print(tA19)

# ==============================================================================
# TABLE A20 - Hazards of resolution by months since disappearance
# ==============================================================================
# Monthly hazards (%) of leaving "missing" for each outcome, as projected: the
# first window's hazards times each outcome's average multiplier over the
# three windows. Below them, each window's multiplier against the first.
model_09 <- read_rds("data_inter/ukr_ualosses_resolution_model.rds")
mult <- model_09$window_multipliers
colnames(mult) <- ual_resolutions
tA20 <-
  bind_rows(
    read_rds("data_inter/ukr_ualosses_resolution_hazards.rds") |>
      mutate(across(-band, \(x) round(100 * x, 3)), row = "monthly hazard (%)", .before = 1) |>
      rename(months_since_event = band),
    as_tibble(round(mult, 3)) |>
      mutate(row = "multiplier against the first window",
             months_since_event = c("v14-v16", "v16-v18", "v18-v19"), .before = 1)
  )
save_tab(tA20, "tableA20_resolution_hazards.csv")
print(tA20)

# ==============================================================================
# TABLE A21 - The projection against the civil register
# ==============================================================================
# The Ministry of Justice counts the deaths and births registered on the
# territory the government controls. The projection covers continental
# Ukraine, occupied parts included. Dividing the registrations by the
# projection's crude rates gives the population they imply at those rates, a
# check on the denominator rather than on the deaths themselves.
fx_by_age <-
  read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds") |>
  select(year, age, fx) |>
  bind_rows(read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds") |> filter(year == 2023) |>
              select(age, fx) |> expand_grid(year = c(2024, 2025)))
vital <-
  sim_wide |>
  left_join(fx_by_age |> mutate(sex = "f"), by = c("year", "sex", "age")) |>
  summarise(population = sum(pop),
            expected = sum(expected),
            conflict = sum(civilian + combatant_confirmed + combatant_imputed),
            births = sum(coalesce(fx, 0) * pop),
            .by = c(sim_id, year)) |>
  mutate(deaths = expected + conflict) |>
  summarise(across(c(population, expected, conflict, deaths, births), median), .by = year)
registered <-
  official |>
  filter(key %in% c("registered_deaths", "registered_births")) |>
  mutate(year = year(as.Date(date))) |>
  select(year, key, value) |>
  pivot_wider(names_from = key, values_from = value)
tA21 <-
  vital |>
  left_join(registered, by = "year") |>
  mutate(population_implied_by_deaths = registered_deaths / (deaths / population),
         population_implied_by_births = registered_births / (births / population)) |>
  transmute(year,
            population_millions = round(population / 1e6, 2),
            projected_deaths = round(deaths), projected_conflict_deaths = round(conflict),
            projected_births = round(births),
            registered_deaths, registered_births,
            population_implied_by_deaths_millions = round(population_implied_by_deaths / 1e6, 2),
            population_implied_by_births_millions = round(population_implied_by_births / 1e6, 2))
save_tab(tA21, "tableA21_civil_register_check.csv")
print(tA21)

# ==============================================================================
# FIGURE A6 - The durations each event month is observed at
# ==============================================================================
# One vertical segment per event month and window between releases, from the
# months since the event at the window's start to those at its end; the
# projection carries each month from its duration in v19 to the horizon.
cells_lexis <-
  model_09$fitted |>
  distinct(month, from_release, d0, d1) |>
  mutate(window = paste0(from_release, "-", ual_releases$release[match(from_release, ual_releases$release) + 1]))
projection_lexis <-
  cells_lexis |>
  distinct(month) |>
  mutate(d19 = ual_months_since(ual_releases$date[nrow(ual_releases)], month)) |>
  filter(d19 < model_09$horizon)
figA6 <-
  ggplot() +
  geom_segment(data = projection_lexis,
               aes(x = month, xend = month, y = d19, yend = model_09$horizon),
               colour = "grey65", linewidth = 0.5, linetype = "22") +
  geom_segment(data = cells_lexis,
               aes(x = month, xend = month, y = d0, yend = d1, colour = window),
               linewidth = 1.1) +
  geom_hline(yintercept = model_09$horizon, linewidth = 0.3) +
  scale_colour_manual(values = c("v14-v16" = "#0A9396", "v16-v18" = "#EE9B00", "v18-v19" = "#AE2012")) +
  scale_y_continuous(breaks = seq(0, 60, 6)) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y") +
  labs(x = "Month of the event", y = "Months since the event", colour = "Window",
       caption = paste0("Solid: the durations observed in each window between register releases. Dashed: ",
                        "the projection from each month's duration in v19 to ", model_09$horizon,
                        " months (line).")) +
  theme_paper()
save_fig(figA6, "figA6_duration_lexis.png", 9, 4.6)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# provenance stamp for the whole table set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The figures carry the draw count in their captions; a .csv cannot, so it is
# recorded once here and covers the whole set. The leading underscore keeps it
# out of the manuscript tables - everything else directly in tables/ is a
# deliverable.
run_provenance <- tibble(
  field = c("n_draws", "seed", "draws_file", "run_finished", "git_commit"),
  value = c(
    format(n_sim),
    "42",
    basename(draws_file),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    tryCatch(
      suppressWarnings(
        system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)
      )[1],
      error = function(e) NA_character_
    )
  )
)
save_tab(run_provenance, "_run_provenance.csv")
print(run_provenance)

message("\nDone. Figures in figures/, tables in tables/.")
