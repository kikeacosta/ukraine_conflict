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
#
# TABLES (written to tables/ as .csv)
#   table1_source_totals.csv                         07_*, 08
#   table2_missing_imputation.csv                    09
#   table3_pert_input_bounds.csv                     10
#   table4_conflict_deaths_by_cause.csv              14   <- new
#   table5_life_expectancy_loss.csv                  14
#   tableA1_source_reconciliation.csv                07_ucdp, 07_acled
#   tableA2_status_transitions.csv                   09
#   table6_totals_by_cause.csv                       14   <- new
#   table7_totals_by_year.csv                        14   <- new
#   tableA3_e0_loss_by_cause.csv                     14
#   tableA4_uncertainty_shares.csv                   11, 14 <- new
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
trans <- read_rds("data_inter/ukr_ualosses_transition_rates.rds")
param_table <- read_rds("data_inter/ukr_param_table.rds")
sce <- read_rds("data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds")
loss_draws <- read_rds("data_inter/ukr_e0_loss_by_cause_draws_2022_2025.rds")
loss_sum <- read_rds("data_inter/ukr_e0_loss_by_cause_summary_2022_2025.rds")
deaths_sum <- read_rds("data_inter/ukr_conflict_deaths_by_cause_summary_2022_2025.rds")
mig <- read_rds("data_inter/ukr_migration_decomposition.rds")

draws_file <- "data_inter/ukr_sim_draws_2022_2025.rds"
if (!file.exists(draws_file)) {
  stop("Run step 11 first: ", draws_file, call. = FALSE)
}
sim_wide <- as_tibble(readRDS(draws_file))
n_sim <- length(unique(sim_wide$sim_id))
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
  filter(status != "prisoner") |>
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
  geom_jitter(width = 0.2, alpha = 0.04, size = 0.1, colour = "black") +
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
  "Civilians", "Registered combatants", "Missing combatants (imputed)"
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
# The three inputs are drawn independently of one another, so the variance of
# the loss partitions additively across them: each share is the variance the
# input carries through its own slope. What linear terms do not account for is
# reported as its own category rather than distributed over the others.
# Civilians and combatants are kept apart rather than summed into one
# "conflict" term. Combatant draws are an order of magnitude larger, so a
# combined term is dominated by combatant variation - which is nearly
# irrelevant to the female loss, where civilians drive almost all of it.
UNC_LAB <- c(
  cvs   = "Civilian deaths",
  cmb   = "Combatant deaths",
  mig_w = "Migration: western",
  mig_r = "Migration: Russia / Belarus",
  resid = "Interaction"
)
COL_UNC <- c(
  "Civilian deaths"             = "#AE2012",
  "Combatant deaths"            = "#333333",
  "Migration: western"          = "#0A9396",
  "Migration: Russia / Belarus" = "#EE9B00",
  "Interaction"                 = "#CCCCCC"
)

# The loss in year y is that year's death RATES: its own deaths over a
# population that is still missing everyone who left in any earlier year. So
# conflict enters as the year's own draw and migration as the cumulative one.
# Using cumulative deaths instead pushes most of the variance into the
# residual, because it is not what drives the rate.
unc_inputs <-
  param_draws |>
  mutate(bucket = case_when(
    role == "civilians"  ~ "cvs",
    role == "combatants" ~ "cmb",
    role == "mig_west"   ~ "mig_w",
    role == "mig_ru_by"  ~ "mig_r"
  )) |>
  summarise(draw = sum(draw), .by = c(sim_id, year, bucket)) |>
  # mig_ru_by is entered once, in 2022, but the people it removes are still
  # absent in 2025, so the grid is filled with zeros before accumulating
  complete(sim_id, year, bucket, fill = list(draw = 0)) |>
  arrange(sim_id, bucket, year) |>
  mutate(
    val = if_else(bucket %in% c("cvs", "cmb"), draw, cumsum(draw)),
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
  mutate(sh = map(data, function(d) {
    s <- c(
      cvs   = first_order(d$loss_total, d$cvs),
      cmb   = first_order(d$loss_total, d$cmb),
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
save_fig(fig6, "fig6_uncertainty_shares.png", 8, 4.6)

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
      "Beta-PERT draws (n = ", scales::comma(n_sim),
      "). Dashed line: mode. Dotted lines: min and max."
    )
  ) +
  theme_paper() +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(figA1, "figA1_pert_draw_distributions.png", 10, 5)

# ==============================================================================
# FIGURE A2 - Distributions of the simulated net emigration totals
# ==============================================================================
# The counterpart of A1 for the migration side. It reads differently from A1
# in one respect: a simulation takes ONE quantile per component and holds it
# across all four years, so within a row the panels move together and the
# spread is a coverage scenario rather than four independent accidents.
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
    x = "Net emigration (millions)", y = "Density",
    colour = "Component", fill = "Component",
    caption = paste0(
      "Beta-PERT draws (n = ", scales::comma(n_sim),
      "). Dashed line: mode. Dotted lines: min and max. Each component is ",
      "drawn once per simulation and held across all four years."
    )
  ) +
  theme_paper() +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(figA2, "figA2_pert_draw_distributions_migration.png", 10, 5)

# ==============================================================================
# FIGURE A3 - Cumulative net emigration implied by the draws
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
  labs(
    x = "Cumulative net emigration 2022-2025 (millions)", y = "Density",
    linetype = "Published estimate",
    caption = paste0(
      "Beta-PERT draws (n = ", scales::comma(n_sim), ").\nThe published ",
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

ual_tot <-
  ual |>
  filter(status != "prisoner") |>
  summarise(n = sum(dx), .by = status)

ual_low <- ual_tot$n[ual_tot$status == "dead"]
ual_high <- sum(ual_tot$n)

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
tab_imput <-
  imput |>
  mutate(
    imputed_alive_total = imputed_alive + imputed_prisoner,
    pct_dead = imputed_dead / missing_stock,
    pct_alive = imputed_alive_total / missing_stock
  ) |>
  select(
    year,
    registered_deaths = confirmados_stock,
    missing = missing_stock,
    imputed_dead,
    pct_dead,
    imputed_alive = imputed_alive_total,
    pct_alive,
    total_deaths = total_estimado
  )

tab_imput <- bind_rows(
  tab_imput,
  tab_imput |>
    summarise(across(c(registered_deaths, missing, imputed_dead,
                       imputed_alive, total_deaths), sum)) |>
    mutate(
      year = NA_integer_,
      pct_dead = imputed_dead / missing,
      pct_alive = imputed_alive / missing
    )
) |>
  mutate(
    across(c(registered_deaths, missing, imputed_dead, imputed_alive,
             total_deaths), ~round(.x)),
    across(c(pct_dead, pct_alive), ~scales::percent(.x, accuracy = 0.1)),
    year = if_else(is.na(year), "Total", as.character(year))
  )

save_tab(tab_imput, "table2_missing_imputation.csv")
print(tab_imput)

# ==============================================================================
# TABLE 3 - PERT simulation input bounds
# ==============================================================================
tab_pert <-
  param_table |>
  mutate(across(c(mode, min, max), ~round(.x))) |>
  pivot_wider(
    names_from = role,
    values_from = c(mode, min, max),
    names_glue = "{role}_{.value}"
  ) |>
  select(
    year,
    combatants_mode, combatants_min, combatants_max,
    civilians_mode, civilians_min, civilians_max,
    mig_west_mode, mig_west_min, mig_west_max,
    mig_ru_by_mode, mig_ru_by_min, mig_ru_by_max
  )

tab_pert <- bind_rows(
  tab_pert |> mutate(year = as.character(year)),
  # mig_ru_by has a 2022 row only, so the other years are empty, not zero
  tab_pert |> summarise(across(-year, \(x) sum(x, na.rm = TRUE))) |> mutate(year = "Total")
)

save_tab(tab_pert, "table3_pert_input_bounds.csv")
print(tab_pert)

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
    `Registered combatants` = dx_cmb_confirmed,
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

save_tab(tab_by_cause, "table6_totals_by_cause.csv")
print(tab_by_cause)

tab_by_year <-
  totals_by("year") |>
  arrange(year) |>
  mutate(year = as.character(year)) |>
  bind_rows(bind_cols(tibble(year = "Total"), grand_row)) |>
  select(year, Females, Males, Total)

save_tab(tab_by_year, "table7_totals_by_year.csv")
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

save_tab(tab_e0, "table5_life_expectancy_loss.csv")
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
# TABLE A2 - Empirical status transitions of the missing
# ==============================================================================
tA2 <-
  trans |>
  filter(status2 != "missing") |>
  mutate(pct_of_changed = n / sum(n), .by = year) |>
  left_join(
    trans |>
      summarise(
        still_missing = sum(n[status2 == "missing"]),
        changed = sum(n[status2 != "missing"]),
        .by = year
      ),
    by = "year"
  ) |>
  mutate(
    missing_start = still_missing + changed,
    pct_of_changed = scales::percent(pct_of_changed, accuracy = 0.1)
  ) |>
  select(year, missing_start, changed, to = status2, counts = n, pct_of_changed) |>
  arrange(year, to)

save_tab(tA2, "tableA2_status_transitions.csv")
print(tA2)

message("\nDone. Figures in figures/, tables in tables/.")
