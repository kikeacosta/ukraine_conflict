# ==============================================================================
# STEP 12 - Figures for the mortality estimates
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Draws the four main descriptive figures from the simulation output of 11:
#
#   1. observed vs counterfactual death rates by age, with credible bands
#   2. the same, as a spaghetti plot of individual draws
#   3. the distribution of the life expectancy loss across draws
#   4. the drawn conflict death totals, against their PERT parameters
#
# The cause-of-death decomposition lives in 14, not here.
#
# INPUTS   data_inter/ukr_sim_draws_2022_2025_n<n_sim>.rds              (from 11)
#          data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds (from 11)
#          data_inter/ukr_param_table.rds                          (from 10)
# OUTPUTS  figures/mort_rates_*.png, figures/losses_boxplot.png,
#          figures/dists_combatants_civilians_draws.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# named for the simulation size set in 00_setup.R, so this reads the draws that
# match the configured n_sim rather than whatever a previous run left behind
draws_file <- sprintf("data_inter/ukr_sim_draws_2022_2025_n%d.rds", n_sim)
if (!file.exists(draws_file)) {
  stop(
    "Simulation draws not found at n_sim = ", n_sim, ". Run step 11 first:\n  ",
    draws_file,
    call. = FALSE
  )
}

# components, one row per draw x year x sex x age
sim_wide <- as_tibble(readRDS(draws_file))

# expand to the long cause format used by the plots
sim_output_raw <-
  sim_wide %>%
  mutate(
    conflict = civilian + combatant_confirmed + combatant_imputed,
    all = expected + conflict
  ) %>%
  pivot_longer(
    c(all, conflict, expected, civilian,
      combatant_confirmed, combatant_imputed),
    names_to = "cause",
    values_to = "dx"
  ) %>%
  mutate(mx = dx / pop)

sce_all_probabilistic <- read_rds(
  "data_inter/ukr_probabilistic_deaths_rates_2022_2025.rds"
)

param_table <- read_rds("data_inter/ukr_param_table.rds")

# the drawn yearly totals are recoverable from the components, so they do not
# need to be stored separately
draws_df <-
  sim_wide %>%
  summarise(
    draw_cvs = sum(civilian),
    draw_cmb = sum(combatant_confirmed + combatant_imputed),
    .by = c(sim_id, year)
  )

# ==============================================================================
# 1. DEATH RATES BY AGE, WITH CREDIBLE BANDS
# ==============================================================================
sce_all_probabilistic |>
  filter(cause %in% c("all", "expected")) |>
  mutate(
    sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males")
  ) |>
  ggplot() +
  geom_ribbon(
    aes(age, ymin = mx_lower, ymax = mx_upper, fill = cause),
    alpha = 0.4
  ) +
  geom_line(aes(age, mx_mean, color = cause)) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_fill_manual(values = c("red", "black")) +
  scale_color_manual(values = c("red", "black")) +
  facet_grid(sex ~ year) +
  labs(y = "Death rates", col = "Cause", fill = "Cause", x = "Age") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave("figures/exploratory/mort_rates_2022_2025_invals.png", w = 8, h = 5)

# ==============================================================================
# 2. THE SAME, AS INDIVIDUAL DRAWS
# ==============================================================================
# thousands of overplotted lines per panel is slow to render and visually identical
# to a subsample, so the spaghetti uses the first 300 draws.
n_show <- 300

sim_output_raw |>
  filter(cause == "all", sim_id <= n_show) |>
  ggplot() +
  geom_line(
    aes(age, mx, group = sim_id, color = "all"),
    lwd = 0.01,
    alpha = 0.05
  ) +
  geom_line(
    data = sce_all_probabilistic |>
      filter(cause == "all") |>
      rename(mx = mx_mean),
    aes(age, mx, color = "all"),
    lwd = 0.5
  ) +
  geom_line(
    data = sim_output_raw |> filter(cause == "expected", sim_id == 1),
    aes(age, mx, color = "expected"),
    lwd = 0.5
  ) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_color_manual(
    name = "cause",
    values = c("all" = "#FF4D4D", "expected" = "#4D4D4D")
  ) +
  facet_grid(sex ~ year) +
  labs(y = "Death rates", x = "Age") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 5, alpha = 1)))

ggsave("figures/exploratory/mort_rates_2022_2025_invals2.png", w = 10, h = 7, dpi = 300)

# ==============================================================================
# 3. DISTRIBUTION OF THE LIFE EXPECTANCY LOSS
# ==============================================================================
# One life table per draw x year x sex, computed column-wise via lt_cols()
# from 00_setup.R. e0 is row 1 of the ex matrix.
n_age <- 101
d <- as.data.table(sim_wide)
setorder(d, sim_id, year, sex, age)
key <- d[seq(1L, .N, by = n_age), .(sim_id, year, sex)]
as_mat <- function(v) matrix(v, nrow = n_age)

POP <- as_mat(d$pop)
DEX <- as_mat(d$expected)
CNF <- as_mat(d$civilian + d$combatant_confirmed + d$combatant_imputed)

e0_war <- lt_cols((DEX + CNF) / POP, key$sex)$ex[1, ]
e0_bsn <- lt_cols(DEX / POP, key$sex)$ex[1, ]

lt3 <-
  as_tibble(key) |>
  mutate(
    ex = e0_war,
    ex_bsn = e0_bsn,
    loss = ex_bsn - ex, # positive = years lost
    sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males")
  )

loss_sum <-
  lt3 |>
  summarise(
    lss_m = mean(loss),
    lss_l = quantile(loss, 0.025),
    lss_u = quantile(loss, 0.975),
    ex_m = mean(ex),
    ex_bsn = mean(ex_bsn),
    .by = c(year, sex)
  ) |>
  arrange(year, sex)

print(loss_sum)

lt3_summary <-
  lt3 |>
  summarise(
    Q1 = quantile(loss, 0.25),
    Median = median(loss),
    Q3 = quantile(loss, 0.75),
    .by = c(year, sex)
  )

lt3 |>
  ggplot(aes(x = sex, y = loss)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, fill = "white") +
  # seeded, so the figure is the same in every run
  geom_point(position = position_jitter(width = 0.2, seed = 1), alpha = 0.05, size = 0.1, color = "black") +
  geom_text(
    data = lt3_summary,
    aes(y = Median, label = sprintf("%.2f", Median)),
    hjust = -1.5, size = 2.5, color = "black", fontface = "bold"
  ) +
  geom_text(
    data = lt3_summary,
    aes(y = Q1, label = sprintf("%.2f", Q1)),
    vjust = 1.5, hjust = -1.8, size = 1.8, color = "black", alpha = 0.7
  ) +
  geom_text(
    data = lt3_summary,
    aes(y = Q3, label = sprintf("%.2f", Q3)),
    vjust = -1 / 1.5, hjust = -1.8, size = 1.8, color = "black", alpha = 0.7
  ) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.9))) +
  facet_grid(~year) +
  labs(y = "Life expectancy loss (years)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave("figures/exploratory/losses_boxplot.png", w = 8, h = 3.5)

# ==============================================================================
# 4. THE DRAWN CONFLICT DEATH TOTALS
# ==============================================================================
# Densities of the draws, with the min / mode / max that defined each
# PERT distribution marked as dashed lines.
cols <- c("#AE2012", "#005F73")

# param_table also carries the migration components; this plot is the conflict
# side only, and 15 draws the migration counterpart as figA2.
param_cnf <- param_table |> filter(role %in% c("combatants", "civilians"))

draws_df |>
  rename(civilians = draw_cvs, combatants = draw_cmb) |>
  pivot_longer(
    c(civilians, combatants),
    names_to = "role",
    values_to = "dts"
  ) |>
  ggplot(aes(x = dts, color = role, fill = role)) +
  geom_density(alpha = 0.6) +
  geom_text(
    data = param_cnf,
    aes(
      x = -Inf,
      y = Inf,
      label = paste0(
        scales::comma(mode, accuracy = 1),
        " (",
        scales::comma(min, accuracy = 1),
        "-",
        scales::comma(max, accuracy = 1),
        ")"
      )
    ),
    col = "black",
    fontface = "bold",
    hjust = -0.1,
    vjust = 1.5,
    size = 2.5,
    show.legend = FALSE
  ) +
  geom_vline(
    data = param_cnf,
    aes(xintercept = mode, color = role),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_vline(
    data = param_cnf,
    aes(xintercept = min, color = role),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_vline(
    data = param_cnf,
    aes(xintercept = max, color = role),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_nested_wrap(role ~ year, scales = "free", ncol = 4) +
  labs(x = "Death counts", y = "Probability", col = "Role", fill = "Role") +
  theme_bw()

ggsave("figures/exploratory/dists_combatants_civilians_draws.png", w = 10, h = 5)

message("Done.")
