# ==============================================================================
# STEP 09l - How stable is the bound on the missing alive?
# ==============================================================================
#
# WHY
# ---
# The imputation rests on one ratio by event year: of the missing whose fate a
# later release settles outside captivity, the share that leave the register
# beyond list maintenance rather than being recorded dead (composition_bound(),
# 00_setup.R). The drop-outs come in clean-ups of the list, in some releases and
# not others, while the resolutions to death arrive in every window. So the
# ratio is a running one: it falls as windows without a clean-up are added and
# jumps when one comes. This step says how far it moves, and why the events of
# 2022 take the bound of the later cohorts.
#
# WHAT IT REPORTS
# ---------------
#   spans      the bound on each cumulative span of releases, v14-v15 to
#              v14-v19, by event year: each cohort's own, and the one used,
#              in which 2022 takes the pooled bound of 2023-2025
#   windows    by event year and window between releases: first resolutions to
#              death and drop-outs beyond list maintenance, the missing at risk
#              at the window's start, the person-months they give (at risk x the
#              window's length) and both as rates per 1,000 person-months
#   cohorts    the same over the five windows, with the share of each cohort's
#              missing at risk in v14 who left or were recorded dead
#   left_out   the bound used with each window left out in turn
#   resampled  the bound used under the window weights the simulation draws
#              (window_weight_draws()), its median and 95% interval
#
# A rate per person-month does not replace the ratio. Exits and deaths compete
# for the same people over the same months, so the ratio of their rates in a
# window is the ratio of their counts; summed over windows it is the bound
# itself. What the rates add is the comparison across cohorts, which the ratio
# hides: a high bound can come from many exits or from few deaths.
#
# INPUTS   data_inter/ukr_military_inputs.rds       (09: the composition)
#          data_inter/ualosses_window_transitions.rds (09's anonymous linkage
#          counts: the missing at risk by event year and window)
# OUTPUTS  data_inter/ukr_bound_stability.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mil <- read_rds("data_inter/ukr_military_inputs.rds")
comp <- mil$composition
wins <- ual_releases$release[-nrow(ual_releases)]
span_name <- \(k) paste0("v14-", ual_releases$release[k + 1])
years <- sort(unique(comp$year))

# 1. THE BOUND ON EACH CUMULATIVE SPAN =========================================
spans <- map_dfr(seq_along(wins), function(k) {
  d <- comp |> filter(from_release %in% wins[seq_len(k)])
  # a span in which no later cohort has resolved at all has no bound to lend
  used <- if (sum(d$out[d$year != 2022] + d$dead[d$year != 2022]) > 0) composition_bound(d)[1, ] else rep(NA_real_, length(years))
  tibble(span = span_name(k), windows = k, year = years,
         own = unname(composition_bound(d, borrow = NULL)[1, ]), used = unname(used),
         pooled = sum(d$out) / sum(d$out + d$dead))
})
stopifnot(isTRUE(all.equal(spans$used[spans$windows == length(wins)], unname(composition_bound(comp)[1, ]))))

# 2. BY EVENT YEAR AND WINDOW, WITH THE MISSING AT RISK =========================
linkage <- read_rds("data_inter/ualosses_window_transitions.rds")
at_risk <-
  linkage |>
  filter(table == "windows", rule == "production", month != as.Date("2022-02-01")) |>
  summarise(at_risk = sum(n), .by = c(year, from_release))
windows <-
  comp |>
  left_join(at_risk, by = c("year", "from_release")) |>
  mutate(months = ual_window_months[match(from_release, wins)],
         person_months = at_risk * months,
         out_rate = 1000 * out / person_months,
         dead_rate = 1000 * dead / person_months)
stopifnot(!anyNA(windows))

cohorts <-
  windows |>
  summarise(dead = sum(dead), out = sum(out), at_risk_v14 = at_risk[from_release == wins[1]],
            person_months = sum(person_months), .by = year) |>
  mutate(out_share = out / at_risk_v14, dead_share = dead / at_risk_v14,
         out_rate = 1000 * out / person_months, dead_rate = 1000 * dead / person_months,
         own = unname(composition_bound(comp, borrow = NULL)[1, ]), used = unname(composition_bound(comp)[1, ]))

# 3. EACH WINDOW LEFT OUT, AND THE WINDOWS RESAMPLED ============================
left_out <- map_dfr(seq_along(wins), function(k) {
  w <- matrix(1, 1, length(wins))
  w[1, k] <- 0
  tibble(left_out = paste(wins[k], ual_releases$release[k + 1], sep = "-"), year = years,
         own = unname(composition_bound(comp, w, borrow = NULL)[1, ]), used = unname(composition_bound(comp, w)[1, ]))
})
set.seed(20260922)
b <- composition_bound(comp, window_weight_draws(20000))
resampled <- tibble(year = years,
                    median = apply(b, 2, median),
                    lo = apply(b, 2, quantile, 0.025, names = FALSE),
                    hi = apply(b, 2, quantile, 0.975, names = FALSE))

cat("\n=== THE BOUND ON EACH CUMULATIVE SPAN (%): OWN, AND USED ===\n")
print(as.data.frame(spans |> mutate(across(c(own, used, pooled), \(v) round(100 * v, 1))) |>
                      pivot_wider(id_cols = c(span, pooled), names_from = year, values_from = c(own, used))))
cat("\n=== DROP-OUTS AND DEATHS PER 1,000 PERSON-MONTHS, BY EVENT YEAR AND WINDOW ===\n")
print(as.data.frame(windows |> mutate(across(c(out, out_rate, dead_rate), \(v) round(v, 2)), person_months = round(person_months))))
cat("\n=== OVER THE FIVE WINDOWS ===\n")
print(as.data.frame(cohorts |> mutate(across(c(out_share, dead_share, own, used), \(v) round(100 * v, 2)),
                                      across(c(out, out_rate, dead_rate), \(v) round(v, 2)))))
cat("\n=== THE BOUND USED WITH EACH WINDOW LEFT OUT, AND RESAMPLED (%) ===\n")
print(as.data.frame(left_out |> mutate(across(c(own, used), \(v) round(100 * v, 1)))))
print(as.data.frame(resampled |> mutate(across(-year, \(v) round(100 * v, 1)))))

write_rds(list(spans = spans, windows = windows, cohorts = cohorts, left_out = left_out, resampled = resampled),
          "data_inter/ukr_bound_stability.rds")
message("Done. data_inter/ukr_bound_stability.rds written.")
