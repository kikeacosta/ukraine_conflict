# ==============================================================================
# STEP 09E - MULTISTATE OUTPUT PLOT: Stacked State-Occupation Probabilities
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# The standard way to display a competing-risks / multistate model's fitted
# output is a stacked state-occupation probability plot: at each time t,
# P(missing) + P(dead) + P(not dead) = 1 by construction, so stacking them
# as areas makes that constraint visible instead of asserted. A single CIF
# line (e.g. 09d's mstate_cif_by_cohort.png, which only plots P(dead)) loses
# the competing-risks story - you can't see what the "not dead" outcome is
# doing, or that the "missing" band is literally where the probability mass
# for "dead"/"not dead" is coming FROM as time passes.
#
# NOTE ON STATES: 09d fits dead vs. "alive" as the two competing outcomes,
# where "alive" already absorbs the raw register's prisoner/released-prisoner
# categories (see 09d section 1 for why: the finer 3-outcome split caused
# severe separation in the Cox fit and isn't needed for the dead/not-dead
# imputation this feeds into). This plot therefore shows 3 stacked bands,
# not 4.
#
# INPUT   data_inter/ukr_ualosses_mstate_cifs.rds   (from 09d)
# OUTPUT  figures/exploratory/mstate_stacked_state_occupation.png
#         figures/exploratory/mstate_horizon_comparison.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

if (!file.exists("data_inter/ukr_ualosses_mstate_cifs.rds")) {
  stop("Run 09d_competing_risks_mstate.R first.", call. = FALSE)
}

cifs <- read_rds("data_inter/ukr_ualosses_mstate_cifs.rds")

# Long format for geom_area(position = "stack"). Order matters: the level
# listed FIRST is stacked at the BOTTOM.
plot_long <- cifs %>%
  select(cohort, time_months, p_dead, p_alive, p_missing) %>%
  pivot_longer(
    cols = starts_with("p_"),
    names_to = "state", names_prefix = "p_",
    values_to = "probability"
  ) %>%
  mutate(
    state = recode(state,
                   dead = "Confirmed dead",
                   alive = "Resolved alive (incl. POW/released)",
                   missing = "Still missing"),
    state = factor(state, levels = c("Confirmed dead",
                                     "Resolved alive (incl. POW/released)",
                                     "Still missing")),
    cohort_label = paste0(cohort, " cohort")
  )

# Palette consistent with the rest of the pipeline's existing red-for-deaths
# convention (see 09_ualosses_missing_analysis.R's cols <- c("grey30",
# "#e2606b", "#e63946")): dark red for the outcome that matters most for
# mortality estimation, grey for the still-unresolved majority, muted green
# for the resolved-alive outcome.
state_colors <- c(
  "Confirmed dead"                       = "#8b1a1a",
  "Resolved alive (incl. POW/released)"  = "#4a7c59",
  "Still missing"                        = "#c9c9c9"
)

p <- ggplot(plot_long, aes(x = time_months, y = probability, fill = state)) +
  geom_area(position = "stack", color = "white", linewidth = 0.15, alpha = 0.95) +
  facet_wrap(~cohort_label, nrow = 1) +
  scale_fill_manual(values = state_colors, name = NULL) +
  # No hard upper limit: stacking sums the 3 bands to their cumulative
  # height, and a hard limits=c(0,1) censors (to NA, silently) any point
  # that floating-point noise pushes fractionally above 1.0 - which then
  # leaves a geom_area() group empty at that x and crashes the grob
  # renderer entirely rather than just clipping a pixel. The CIF sum-to-1
  # check in 09d already guarantees correctness to 1e-6, so the limit was
  # redundant protection that carried real crash risk.
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Months since disappearance",
    y = "State-occupation probability",
    title = "Multistate Model: Where the Missing-Persons Cohort Goes Over Time",
    subtitle = paste0(
      "Competing-risks cumulative incidence (mstate cause-specific hazards). ",
      "Bands sum to 100% at every timepoint by construction."
    ),
    caption = "Source: ualosses registers v14/v18/v19, cause-specific Cox model stratified by transition (code/09d_competing_risks_mstate.R)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(size = 7, color = "grey50")
  )

ggsave("figures/exploratory/mstate_stacked_state_occupation.png", p, w = 11, h = 4.5, dpi = 150)
message("✓ Saved: figures/exploratory/mstate_stacked_state_occupation.png")

# ==============================================================================
# SECONDARY PLOT: transition probabilities at a fixed horizon, by cohort
# ==============================================================================
# Complements the stacked plot with a direct "read the number off" view at
# one common horizon, useful for a text callout or supplementary table figure.

horizon_months <- 12
at_horizon <- cifs %>%
  mutate(dist_to_horizon = abs(time_months - horizon_months)) %>%
  slice_min(dist_to_horizon, by = cohort, n = 1) %>%
  select(cohort, p_dead, p_dead_lo, p_dead_hi, p_alive, p_missing)

message(sprintf("\n=== State occupation at ~%d months since disappearance ===", horizon_months))
print(at_horizon)

p2 <- at_horizon %>%
  pivot_longer(c(p_dead, p_alive, p_missing),
              names_to = "state", names_prefix = "p_", values_to = "probability") %>%
  mutate(
    state = recode(state, dead = "Confirmed dead",
                   alive = "Resolved alive (incl. POW/released)",
                   missing = "Still missing"),
    state = factor(state, levels = c("Still missing", "Confirmed dead",
                                     "Resolved alive (incl. POW/released)"))
  ) %>%
  ggplot(aes(x = factor(cohort), y = probability, fill = state)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = state_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Disappearance cohort", y = "Probability",
    title = sprintf("State Occupation at %d Months Since Disappearance, by Cohort", horizon_months)
  ) +
  theme_bw() + theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave("figures/exploratory/mstate_horizon_comparison.png", p2, w = 7, h = 5, dpi = 150)
message("✓ Saved: figures/exploratory/mstate_horizon_comparison.png")
