# ==============================================================================
# STEP 09j - The projection's own rule against a window it was not fitted to
# ==============================================================================
#
# WHY
# ---
# S8.6 sets out why the two interior releases cannot test the duration model: a
# half-window predicted with the multiplier of the window it sits inside
# reproduces a fitted total by construction. This is the test that can be run
# with the same data.
#
# The projection carries each event month's missing to 48 months with the
# LENGTH-WEIGHTED AVERAGE of the fitted windows' multipliers (ual_theta_hazards(),
# the `projection` element). Whether that rule forecasts a window it has not seen
# is a question the three windows can answer: fit on two of them and predict the
# third.
#
#   case A   fit v14>v16 and v16>v18, predict v18>v19
#   case B   fit v16>v18 and v18>v19, predict v14>v16
#   case C   fit v14>v16 and v18>v19, predict v16>v18
#
# Reported by cause, because the causes do not matter equally. By the identity of
# S3.5 a death the model misses returns to the unresolved pool and is imputed dead
# at 1 - a anyway, so the military total barely moves; a resolution ALIVE that the
# model misses does change it. The death channel is what carries the estimate, and
# the captivity channel is expected to miss badly in case A, because the last
# release recorded a batch of past returns that no earlier window anticipates.
# That is a property of the register worth reporting, not a failure of the model.
#
# INPUTS   data_inter/ualosses_window_transitions.rds (09's cached linkage),
#          data_inter/ukr_registration_completion.rds (08b)
# OUTPUTS  data_inter/ukr_ualosses_held_out_window.rds:
#            by_cause   observed and predicted, per case and outcome, with the
#                       interval the simulation's two draws of the model give
#            by_band    the same by months since the event
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

feb22 <- as.Date("2022-02-01")
completion <- read_rds("data_inter/ukr_registration_completion.rds")
linkage <- read_rds("data_inter/ualosses_window_transitions.rds")

# 1. THE SAME COUNTS THE PRODUCTION MODEL IS FITTED TO ==========================
# List maintenance, exactly as 09 applies it: only drop-out beyond the rate at
# which the dead leave the register counts as a resolution alive.
dead_rate <-
  linkage |>
  filter(table == "dead_windows", month != feb22) |>
  summarise(at_risk = sum(n), dropped = sum(n[to == "no_longer_listed"]),
            .by = c(year, from_release)) |>
  mutate(rate = dropped / at_risk)

maintained <-
  linkage |>
  filter(table == "windows", rule == "production", month != feb22) |>
  select(year, month, entry, from_release, to, n) |>
  left_join(dead_rate |> select(year, from_release, rate), by = c("year", "from_release")) |>
  mutate(maint_share = pmin(1, rate * sum(n) / sum(n[to == "no_longer_listed"])),
         .by = c(year, from_release)) |>
  mutate(n_maint = if_else(to == "no_longer_listed", n * maint_share, 0))
maintained <-
  bind_rows(maintained |> mutate(n = n - n_maint),
            maintained |> filter(n_maint > 0) |> mutate(to = "missing", n = n_maint)) |>
  summarise(n = sum(n), .by = c(year, month, entry, from_release, to))

# Durations come from the full four-release calendar, so the held-out window's
# cells keep the durations the projection would meet.
cells_all <- ual_duration_cells(maintained)
rel_all <- ual_releases
months_all <- ual_window_months

# 2. FIT ON TWO WINDOWS, PREDICT THE THIRD ======================================
# ual_fit_durations() reads ual_releases and ual_window_months, so the two are
# restricted around the fit and restored after it.
fit_without <- function(keep_rows, fit_windows) {
  ual_releases <<- rel_all[keep_rows, ]
  ual_window_months <<- months_all[head(keep_rows, -1)]
  on.exit({ ual_releases <<- rel_all; ual_window_months <<- months_all })
  ual_fit_durations(cells_all |> filter(from_release %in% fit_windows),
                    horizon = completion$L_ref)
}

predict_window <- function(m, w) {
  cells <- cells_all |> filter(from_release == w)
  p <- ual_resolve(m$h, cells$d0, cells$d1, m$breaks)   # m$h is the projection rule
  at_risk <- rowSums(cells[, c("missing", ual_resolutions)])
  bind_cols(cells |> select(month, from_release, d0, d1),
            tibble(at_risk = at_risk),
            as_tibble(p[, ual_resolutions, drop = FALSE] * at_risk) |>
              rename_with(\(x) paste0("pred_", x)))
}

cases <- tribble(
  ~case, ~keep,       ~fit_windows,        ~held_out,
  "A",   1:3,         c("v14", "v16"),     "v18",
  "B",   2:4,         c("v16", "v18"),     "v14",
  "C",   c(1, 3, 4),  c("v14", "v18"),     "v16"
)

models <- pmap(cases, function(case, keep, fit_windows, held_out) {
  message("  case ", case, ": fitting on ", paste(fit_windows, collapse = " and "),
          ", to predict ", held_out)
  fit_without(keep, fit_windows)
}) |> set_names(cases$case)

runs <- pmap(cases, function(case, keep, fit_windows, held_out) {
  obs <- cells_all |> filter(from_release == held_out) |> select(month, all_of(ual_resolutions))
  predict_window(models[[case]], held_out) |>
    left_join(obs, by = "month") |>
    mutate(case = case, fitted_on = paste(fit_windows, collapse = " + "), .before = 1)
})
cmp <- list_rbind(runs)

# 3. BY CAUSE AND BY DURATION ===================================================
by_cause <-
  cmp |>
  summarise(at_risk = sum(at_risk),
            across(c(all_of(ual_resolutions), starts_with("pred_")), sum),
            .by = c(case, fitted_on, from_release)) |>
  pivot_longer(c(all_of(ual_resolutions), starts_with("pred_")),
               names_to = "k", values_to = "v") |>
  mutate(what = if_else(str_starts(k, "pred_"), "predicted", "observed"),
         outcome = str_remove(k, "^pred_")) |>
  select(-k) |>
  pivot_wider(names_from = what, values_from = v) |>
  mutate(ratio = predicted / observed)

band_of <- function(d) cut(d, ual_duration_breaks, right = FALSE,
                           labels = paste0(head(ual_duration_breaks, -1), "-",
                                           ual_duration_breaks[-1]))
by_band <-
  cmp |>
  mutate(band = band_of(d0)) |>
  summarise(at_risk = sum(at_risk),
            across(c(all_of(ual_resolutions), starts_with("pred_")), sum),
            .by = c(case, band)) |>
  arrange(case, band)

# 4. DOES THE SPREAD THE SIMULATION CARRIES COVER THE ERROR? ====================
# In every simulation the resolution model varies in two ways (11): its estimates
# are drawn from their sampling distribution, and the projection averages the
# fitted windows' multipliers under weights drawn around the windows' lengths
# (window_weight_draws(), 00_setup.R). The same two draws, made from each case's
# fit, give the interval the simulation would have put around the held-out
# window's resolutions. An observed count outside it is an error the
# simulation's interval does not carry: the interval is for the average of many
# future windows, and one window can sit far from it.
set.seed(20260920)
n_pred <- 2000
predictive <- pmap_dfr(cases, function(case, keep, fit_windows, held_out) {
  m <- models[[case]]
  lens <- months_all[head(keep, -1)]
  n_theta <- length(m$theta)
  theta <- sweep(matrix(rnorm(n_pred * n_theta), n_pred) %*% chol(m$vcov + diag(1e-12, n_theta)),
                 2, m$theta, "+")
  w <- window_weight_draws(n_pred, lens)
  cells <- cells_all |> filter(from_release == held_out)
  at_risk <- rowSums(cells[, c("missing", ual_resolutions)])
  pred <- vapply(seq_len(n_pred), function(i) {
    hz <- ual_theta_hazards(theta[i, ], m$nb, m$k, windows = lens)
    p <- ual_resolve(sweep(hz$h, 2, colSums(hz$mult * w[i, ]), "*"), cells$d0, cells$d1, m$breaks)
    colSums(p[, ual_resolutions, drop = FALSE] * at_risk)
  }, numeric(length(ual_resolutions)))
  tibble(case = case, outcome = ual_resolutions,
         pred_lo = apply(pred, 1, quantile, 0.025), pred_hi = apply(pred, 1, quantile, 0.975))
})
by_cause <-
  by_cause |>
  left_join(predictive, by = c("case", "outcome")) |>
  mutate(covered = observed >= pred_lo & observed <= pred_hi)

cat("\n=== THE PROJECTION'S RULE AGAINST A WINDOW IT WAS NOT FITTED TO ===\n")
print(as.data.frame(by_cause |> mutate(across(where(is.double), \(x) round(x, 3)))))
cat("\n=== THE SAME, BY MONTHS SINCE THE EVENT ===\n")
print(as.data.frame(by_band |> mutate(across(where(is.double), \(x) round(x, 1)))))

write_rds(list(by_cause = by_cause, by_band = by_band, cells = cmp,
               dispersion = map_dbl(models, "dispersion")),
          "data_inter/ukr_ualosses_held_out_window.rds")
cat("\nWritten: data_inter/ukr_ualosses_held_out_window.rds\n")
