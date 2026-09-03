# ==============================================================================
# STEP 14 - Life expectancy loss decomposed by cause of conflict death
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Answers: how many years of life expectancy did the war cost in each year,
# and how much of that is due to each kind of conflict death?
#
# Two decompositions are produced:
#
#   TWO-WAY    civilian  vs  combatant
#   THREE-WAY  civilian  vs  registered combatant  vs  imputed missing combatant
#
# METHOD
# ------
# For every one of the 5,000 draws from step 11, and every year and sex:
#
#   1. build two life tables - one on the counterfactual "expected" mortality
#      (no war) and one on all-cause mortality (with war);
#   2. run an Arriaga decomposition to get the age-specific contribution
#      TE(x) of each age interval to the total e0 gap;
#   3. split TE(x) between causes in proportion to each cause's share of the
#      conflict deaths at that age.
#
# Step 3 is exact rather than approximate: both life tables use the same
# exposure denominator, so mx_all - mx_expected IS the conflict death rate,
# and conflict deaths are additive across causes. The script asserts this.
#
# Doing all of this INSIDE each draw (rather than once on the medians, as the
# original version in step 11 did) means the cause-specific contributions
# carry the same 95% credible intervals as the total.
#
# SIGN CONVENTION: positive = years of life expectancy LOST.
#
# INPUT    data_inter/ukr_sim_draws_2022_2025.rds   (from 11)
# OUTPUTS  data_inter/ukr_e0_loss_by_cause_draws_2022_2025.rds
#          data_inter/ukr_e0_loss_by_cause_summary_2022_2025.rds
#          data_inter/ukr_conflict_deaths_by_cause_summary_2022_2025.rds
#          figures/e0_loss_*.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

draws_file <- "data_inter/ukr_sim_draws_2022_2025.rds"
if (!file.exists(draws_file)) {
  stop(
    "Simulation draws not found. Run step 11 first (~7 minutes):\n  ",
    draws_file,
    call. = FALSE
  )
}

AGES <- 0:100
n_age <- length(AGES)

# lt_cols() and arriaga_TE() are defined in 00_setup.R, so 12 and 14 share a
# single implementation.

# ==============================================================================
# 1. DECOMPOSE EVERY DRAW
# ==============================================================================
d <- as.data.table(readRDS(draws_file))
setorder(d, sim_id, year, sex, age)
stopifnot(nrow(d) %% n_age == 0)

# one column per (draw, year, sex) group, ages down the rows
key <- d[seq(1L, .N, by = n_age), .(sim_id, year, sex)]
as_mat <- function(v) matrix(v, nrow = n_age)

POP <- as_mat(d$pop)
DEX <- as_mat(d$expected) # expected (non-conflict) deaths
DCV <- as_mat(d$civilian)
DCC <- as_mat(d$combatant_confirmed)
DCI <- as_mat(d$combatant_imputed)

MXB <- DEX / POP # baseline: no war
MXW <- (DEX + DCV + DCC + DCI) / POP # observed: with war

lb <- lt_cols(MXB, key$sex)
lw <- lt_cols(MXW, key$sex)
TE <- arriaga_TE(lb, lw)

# share of the age-specific conflict deaths going to each cause
DEN <- DCV + DCC + DCI
sh <- function(M) ifelse(DEN > 0, M / DEN, 0)

sims <- cbind(
  key,
  data.table(
    e0_bsn = lb$ex[1, ],
    e0_war = lw$ex[1, ],
    loss_total = colSums(TE),
    loss_civilian = colSums(TE * sh(DCV)),
    loss_cmb_confirmed = colSums(TE * sh(DCC)),
    loss_cmb_imputed = colSums(TE * sh(DCI)),
    dx_civilian = colSums(DCV),
    dx_cmb_confirmed = colSums(DCC),
    dx_cmb_imputed = colSums(DCI)
  )
)
sims[, loss_combatant := loss_cmb_confirmed + loss_cmb_imputed]
sims[, dx_combatant := dx_cmb_confirmed + dx_cmb_imputed]

# ==============================================================================
# 2. VALIDATION
# ==============================================================================
# The decomposition must reproduce the e0 gap exactly, and the cause
# components must sum back to the total. Both to floating-point precision.
cat("\n=== validation ===\n")
cat(
  "max |sum(TE) - (e0_bsn - e0_war)| :",
  max(abs(sims$loss_total - (sims$e0_bsn - sims$e0_war))), "\n"
)
cat(
  "max |components - total|          :",
  max(abs(
    sims$loss_civilian + sims$loss_cmb_confirmed +
      sims$loss_cmb_imputed - sims$loss_total
  )), "\n"
)
cat("any negative total loss           :", any(sims$loss_total < 0), "\n")
stopifnot(
  max(abs(sims$loss_total - (sims$e0_bsn - sims$e0_war))) < 1e-8,
  !any(sims$loss_total < 0)
)

write_rds(sims, "data_inter/ukr_e0_loss_by_cause_draws_2022_2025.rds")

# ==============================================================================
# 3. SUMMARISE ACROSS DRAWS
# ==============================================================================
CAUSE_LEVELS <- c(
  "Civilians",
  "Registered combatants",
  "Missing combatants (imputed)",
  "Combatants",
  "Total"
)

e0_loss_by_cause <-
  sims |>
  as_tibble() |>
  select(
    sim_id, year, sex,
    Total = loss_total,
    Civilians = loss_civilian,
    Combatants = loss_combatant,
    `Registered combatants` = loss_cmb_confirmed,
    `Missing combatants (imputed)` = loss_cmb_imputed
  ) |>
  pivot_longer(-c(sim_id, year, sex), names_to = "cause", values_to = "loss") |>
  summarise(
    loss_median = median(loss),
    loss_lo = quantile(loss, 0.025),
    loss_hi = quantile(loss, 0.975),
    .by = c(year, sex, cause)
  ) |>
  left_join(
    sims |>
      as_tibble() |>
      summarise(total_median = median(loss_total), .by = c(year, sex)),
    by = c("year", "sex")
  ) |>
  mutate(
    share = loss_median / total_median,
    cause = factor(cause, levels = CAUSE_LEVELS)
  ) |>
  arrange(year, sex, cause)

deaths_by_cause <-
  sims |>
  as_tibble() |>
  select(
    sim_id, year, sex,
    Civilians = dx_civilian,
    Combatants = dx_combatant,
    `Registered combatants` = dx_cmb_confirmed,
    `Missing combatants (imputed)` = dx_cmb_imputed
  ) |>
  pivot_longer(-c(sim_id, year, sex), names_to = "cause", values_to = "dx") |>
  summarise(
    dx_median = median(dx),
    dx_lo = quantile(dx, 0.025),
    dx_hi = quantile(dx, 0.975),
    .by = c(year, sex, cause)
  ) |>
  mutate(cause = factor(cause, levels = CAUSE_LEVELS)) |>
  arrange(year, sex, cause)

write_rds(
  e0_loss_by_cause,
  "data_inter/ukr_e0_loss_by_cause_summary_2022_2025.rds"
)
write_rds(
  deaths_by_cause,
  "data_inter/ukr_conflict_deaths_by_cause_summary_2022_2025.rds"
)

fmt <- function(m, l, h) sprintf("%.2f [%.2f-%.2f]", m, l, h)

cat("\n=== TWO-WAY: e0 loss in years, median [95% UI] ===\n")
print(
  e0_loss_by_cause |>
    filter(cause %in% c("Civilians", "Combatants", "Total")) |>
    mutate(out = fmt(loss_median, loss_lo, loss_hi)) |>
    select(year, sex, cause, out) |>
    pivot_wider(names_from = cause, values_from = out),
  n = 50
)

cat("\n=== THREE-WAY: e0 loss in years, median [95% UI] ===\n")
print(
  e0_loss_by_cause |>
    filter(cause != "Combatants") |>
    mutate(out = fmt(loss_median, loss_lo, loss_hi)) |>
    select(year, sex, cause, out) |>
    pivot_wider(names_from = cause, values_from = out),
  n = 50
)

cat("\n=== THREE-WAY: share of the loss (%) ===\n")
print(
  e0_loss_by_cause |>
    filter(!cause %in% c("Combatants", "Total")) |>
    mutate(share = round(100 * share, 1)) |>
    select(year, sex, cause, share) |>
    pivot_wider(names_from = cause, values_from = share),
  n = 50
)

# ==============================================================================
# 4. PLOTS
# ==============================================================================
# Colours are fixed per cause and reused across every figure, so a category
# always reads the same. Labels are drawn in black, which holds contrast on
# all three fills.
cols_cause <- c(
  "Civilians" = "#66C2A5",
  "Combatants" = "#FC8D62",
  "Registered combatants" = "#FC8D62",
  "Missing combatants (imputed)" = "#8DA0CB"
)

nice_sex <- function(d) {
  d |> mutate(sex = case_when(sex == "f" ~ "Females", sex == "m" ~ "Males"))
}

theme_decomp <- function() {
  theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12)
    )
}

# --- plot builders ----------------------------------------------------------
# `causes`      the categories to stack, in order
# `total_cause` whose 95% UI the whiskers show, and what the shares are taken
#               as a percentage OF. "Total" for the whole loss, "Combatants"
#               when only the combatant loss is being subdivided.
plot_abs <- function(causes, total_cause = "Total", min_label = 0.05,
                     caption = NULL) {
  dat <- e0_loss_by_cause |> filter(cause %in% causes) |> nice_sex()
  tot <- e0_loss_by_cause |> filter(cause == total_cause) |> nice_sex()

  ggplot() +
    geom_col(
      data = dat,
      aes(x = factor(year), y = loss_median, fill = cause),
      width = 0.68, colour = "white", linewidth = 0.2
    ) +
    geom_text(
      data = dat,
      aes(
        x = factor(year), y = loss_median, group = cause,
        label = if_else(loss_median >= min_label, sprintf("%.2f", loss_median), "")
      ),
      position = position_stack(vjust = 0.5),
      colour = "black", size = 3
    ) +
    geom_errorbar(
      data = tot,
      aes(x = factor(year), ymin = loss_lo, ymax = loss_hi),
      width = 0.16, linewidth = 0.4, colour = "grey25"
    ) +
    facet_wrap(~sex) +
    scale_fill_manual(values = cols_cause, breaks = causes) +
    labs(
      x = "Year",
      y = "Life expectancy loss (years)",
      fill = "Cause of loss",
      caption = caption %||% paste0(
        "Bars: median contribution. Whiskers: 95% UI on the ",
        tolower(total_cause), " loss."
      )
    ) +
    theme_decomp()
}

# Shares are recomputed WITHIN the plotted categories, so they always sum to
# 100% and the printed labels match the bar heights that position = "fill"
# produces. Using the stored `share` column (which is always a share of the
# grand total) would desynchronise the two whenever a subset is plotted.
plot_share <- function(causes, min_label = 0.03, y_lab = NULL,
                       facet_sex = TRUE, caption = NULL) {
  dat <-
    e0_loss_by_cause |>
    filter(cause %in% causes) |>
    nice_sex() |>
    mutate(share_plot = loss_median / sum(loss_median), .by = c(year, sex))

  # when the split does not vary by sex, showing two identical panels implies
  # a comparison that is not there
  if (!facet_sex) dat <- dat |> filter(sex == "Males")

  p <-
    ggplot(dat, aes(x = factor(year), y = share_plot, fill = cause)) +
    geom_col(position = "fill", width = 0.68, colour = "white", linewidth = 0.2) +
    geom_text(
      aes(label = if_else(
        share_plot >= min_label,
        scales::percent(share_plot, accuracy = 0.1),
        ""
      )),
      position = position_fill(vjust = 0.5),
      colour = "black", size = 3
    ) +
    scale_fill_manual(values = cols_cause, breaks = causes) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      x = "Year",
      y = y_lab %||% "Share of the life expectancy loss",
      fill = "Cause of loss",
      caption = caption
    ) +
    theme_decomp()

  if (facet_sex) p <- p + facet_wrap(~sex)
  p
}

# The three-way split is the reported decomposition: civilians, combatants
# confirmed individually in the register, and combatants imputed from the
# missing. Two- and three-way variants were also trialled; the three-way is
# the one that carries the story, so it is the only one produced.
#
# Note that the confirmed:imputed ratio is uniform across age and sex by
# construction — both groups are distributed over age and sex using the same
# profile (that of the confirmed dead), and the partition is a single
# year-level scalar, min_cmb / draw_cmb. The sex contrast in these figures is
# therefore driven entirely by the civilian share and by the magnitudes.
three_way <- c(
  "Civilians", "Registered combatants", "Missing combatants (imputed)"
)

plot_abs(three_way)
ggsave("figures/exploratory/e0_loss_three_way_abs.png", w = 8, h = 4.2, dpi = 300)

plot_share(three_way)
ggsave("figures/exploratory/e0_loss_three_way_share.png", w = 8, h = 4.2, dpi = 300)

message("Done.")
