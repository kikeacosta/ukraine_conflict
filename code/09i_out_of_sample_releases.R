# ==============================================================================
# STEP 09i - The resolution model against the releases it was not fitted to
# ==============================================================================
#
# WHY
# ---
# The imputation of the missing rests on one assumption that the data used to
# fit it cannot test: that resolution depends on the months since the
# disappearance and not on the calendar, so that the cohorts of 2024-2025,
# observed only at short durations, will resolve as those of 2022 did at the
# same durations (Methods 3.2). Two further releases of the register, v15 of
# 2 November 2025 and v17 of 6 February 2026, fall INSIDE the first two windows
# the model is fitted to and are held out of every fit. They cut those windows
# in two and say where the register stood part of the way through, which is
# what the model claims to know.
#
# THE TEST
# --------
# The missing and the dead are followed through all six releases by the rules of
# 09, giving five sub-windows. Four of them - v14>v15, v15>v16, v16>v17 and
# v17>v18 - have a boundary the fit never saw. The production model predicts
# each one from the duration each event month had reached at its start, with the
# multiplier of the fitted window it falls in, and the prediction is compared
# with what the register did.
#
# What a difference means. The multiplier of a fitted window is that window's
# average: a release that records a batch of resolutions at once puts them in
# one of the two halves, so an excess in one sub-window that a shortfall in the
# other cancels is evidence about when a release recorded, not about the
# duration model. A difference that does NOT cancel across the pair, or that
# runs the same way at every duration, is evidence about the model.
#
# INPUTS   the six releases in data_input/ualosses_hubert_datasets/,
#          data_inter/ukr_ualosses_resolution_model.rds (09)
# OUTPUTS  data_inter/ukr_ualosses_out_of_sample.rds:
#            by_window  observed and predicted resolutions per sub-window
#            by_band    the same by duration band, over the held-out windows
#            summary    the totals and the ratio predicted / observed
#          data_inter/ualosses_window_transitions_six.rds (cached linkage)
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# The four releases the estimates use, with v15 and v17 between them. The file
# names carry the date each file was saved; the dates are Kaggle's versions.
releases_6 <- tibble(
  release = c("v14", "v15", "v16", "v17", "v18", "v19"),
  date = as.Date(c("2025-09-16", "2025-11-02", "2025-12-04",
                   "2026-02-06", "2026-04-23", "2026-07-21")),
  path = file.path("data_input/ualosses_hubert_datasets", c(
    "250916_UKR_ualosses_Personnel_v14.xlsx", "251102_UKR_ualosses_Personnel_v15.xlsx",
    "251204_UKR_ualosses_Personnel_v16.xlsx", "260206_UKR_ualosses_Personnel_v17.xlsx",
    "260423_UKR_ualosses_Personnel_v18.xlsx", "260919_UKR_ualosses_Personnel_v19.xlsx"))
)
stopifnot(all(file.exists(releases_6$path)))

# ual_windows() and ual_duration_cells() read `ual_releases`; the six-release
# linkage needs all six. Nothing is fitted here - the model comes from 09 - so
# the window months of the production fit are left as they are.
ual_releases <- releases_6

model <- read_rds("data_inter/ukr_ualosses_resolution_model.rds")
prod_rel <- c("v14", "v16", "v18")                  # the fitted windows, by their first release
parent <- c(v14 = "v14", v15 = "v14", v16 = "v16", v17 = "v16", v18 = "v18")
held_out <- c("v14", "v15", "v16", "v17")           # sub-windows with an unseen boundary
feb22 <- as.Date("2022-02-01")

# 1. THE SIX-RELEASE LINKAGE ===================================================
# The same two searches, the same rules as 09, over six releases instead of four.
linkage6 <- cache_rds(
  "data_inter/ualosses_window_transitions_six.rds",
  {
    regs <- map(set_names(ual_releases$path, ual_releases$release), ual_read_release)
    hist <- ual_follow_missing(regs)
    # the dead followed the same way, the two labels swapped, so that a dead
    # record leaving the register shows as "no longer listed"
    swap <- \(r) mutate(r, status = case_when(status == "dead" ~ "missing",
                                              status == "missing" ~ "dead",
                                              .default = status))
    hist_dead <- ual_follow_missing(map(regs, swap))
    bind_rows(
      ual_windows(hist, "alive", released = "alive") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "windows"),
      ual_windows(hist_dead, "alive", released = "exclude") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "dead_windows")
    )
  }
)

# 2. LIST MAINTENANCE, AS IN 09 ================================================
# Only drop-out beyond the rate at which the dead leave the register counts as
# a resolution alive; the rest is held as still missing.
dead_rate6 <-
  linkage6 |>
  filter(table == "dead_windows", month != feb22) |>
  summarise(at_risk = sum(n), dropped = sum(n[to == "no_longer_listed"]),
            .by = c(year, from_release)) |>
  mutate(rate = dropped / at_risk)

maintain6 <- function(w) {
  w <-
    w |>
    left_join(dead_rate6 |> select(year, from_release, rate), by = c("year", "from_release")) |>
    mutate(maint_share = pmin(1, rate * sum(n) / sum(n[to == "no_longer_listed"])),
           .by = c(year, from_release)) |>
    mutate(n_maint = if_else(to == "no_longer_listed", n * maint_share, 0))
  bind_rows(w |> mutate(n = n - n_maint),
            w |> filter(n_maint > 0) |> mutate(to = "missing", n = n_maint)) |>
    summarise(n = sum(n), .by = c(year, month, entry, from_release, to))
}

observed <-
  linkage6 |>
  filter(table == "windows", month != feb22) |>
  select(year, month, entry, from_release, to, n) |>
  maintain6()

# 3. WHAT THE MODEL PREDICTS FOR EACH SUB-WINDOW ===============================
# Each cell is an event month over one sub-window: missing at its start, the
# durations d0 to d1 it spans, and where the register put them by its end.
cells <- ual_duration_cells(observed)
at_risk <- rowSums(cells[, c("missing", ual_resolutions)])

pred <- matrix(0, nrow(cells), length(ual_resolutions) + 1,
               dimnames = list(NULL, c("missing", ual_resolutions)))
for (w in prod_rel) {
  i <- parent[cells$from_release] == w
  if (!any(i)) next
  h <- sweep(model$h_reference, 2, model$window_multipliers[match(w, prod_rel), ], "*")
  pred[i, ] <- ual_resolve(h, cells$d0[i], cells$d1[i], model$breaks)
}

cmp <-
  cells |>
  select(month, from_release, d0, d1, all_of(c("missing", ual_resolutions))) |>
  mutate(at_risk = at_risk, held_out = from_release %in% held_out) |>
  bind_cols(as_tibble(pred * at_risk) |> rename_with(\(x) paste0("pred_", x)))

# 4. BY WINDOW, BY BAND, AND IN TOTAL ==========================================
by_window <-
  cmp |>
  summarise(across(c(at_risk, all_of(ual_resolutions), starts_with("pred_")), sum),
            .by = c(from_release, held_out)) |>
  mutate(resolved = dead + prisoner + no_longer_listed,
         pred_resolved = pred_dead + pred_prisoner + pred_no_longer_listed,
         ratio = pred_resolved / resolved) |>
  arrange(match(from_release, ual_releases$release))

band_of <- function(d) cut(d, model$breaks, labels = paste0(head(model$breaks, -1), "-",
                                                            model$breaks[-1]), right = FALSE)
by_band <-
  cmp |>
  filter(held_out) |>
  mutate(band = band_of(d0)) |>
  summarise(across(c(at_risk, all_of(ual_resolutions), starts_with("pred_")), sum), .by = band) |>
  mutate(resolved = dead + prisoner + no_longer_listed,
         pred_resolved = pred_dead + pred_prisoner + pred_no_longer_listed,
         ratio = pred_resolved / resolved) |>
  arrange(band)

# the pair of sub-windows of each fitted window: a batch recorded at one moment
# lands in one half, so the pair together is the fairer comparison
by_pair <-
  cmp |>
  filter(held_out) |>
  mutate(fitted_window = parent[from_release]) |>
  summarise(across(c(at_risk, all_of(ual_resolutions), starts_with("pred_")), sum),
            .by = fitted_window) |>
  mutate(resolved = dead + prisoner + no_longer_listed,
         pred_resolved = pred_dead + pred_prisoner + pred_no_longer_listed,
         ratio = pred_resolved / resolved)

summary_all <-
  cmp |>
  filter(held_out) |>
  summarise(across(c(at_risk, all_of(ual_resolutions), starts_with("pred_")), sum)) |>
  mutate(resolved = dead + prisoner + no_longer_listed,
         pred_resolved = pred_dead + pred_prisoner + pred_no_longer_listed,
         ratio = pred_resolved / resolved)

cat("\n=== THE REGISTER'S SIX RELEASES ===\n")
print(as.data.frame(releases_6 |> select(release, date)))
cat("\n=== OBSERVED AND PREDICTED RESOLUTIONS BY SUB-WINDOW ===\n")
cat("(held_out = the fit never saw this boundary; v18 is a fitted window)\n")
print(as.data.frame(by_window |> mutate(across(where(is.double), \(x) round(x, 3)))))
cat("\n=== THE HELD-OUT SUB-WINDOWS, BY MONTHS SINCE THE EVENT ===\n")
print(as.data.frame(by_band |> mutate(across(where(is.double), \(x) round(x, 3)))))
cat("\n=== THE TWO HALVES OF EACH FITTED WINDOW TOGETHER ===\n")
print(as.data.frame(by_pair |> mutate(across(where(is.double), \(x) round(x, 3)))))
cat("\n=== THE HELD-OUT WINDOWS IN TOTAL ===\n")
print(as.data.frame(summary_all |> mutate(across(where(is.double), \(x) round(x, 3)))))

write_rds(list(releases = releases_6, by_window = by_window, by_band = by_band,
               by_pair = by_pair, summary = summary_all, cells = cmp),
          "data_inter/ukr_ualosses_out_of_sample.rds")
cat("\nWritten: data_inter/ukr_ualosses_out_of_sample.rds\n")
