# ==============================================================================
# STEP 13t - Civilian deaths: uncertainUCDP (Vesco et al. 2026) and the Beta-PERT
# ==============================================================================
#
# WHY
# ---
# Vesco, Randahl, Hegre, Hogbladh and Yilmaz (2026, Journal of Conflict
# Resolution 70(7-8)) show that UCDP's low and high are not an uncertainty
# interval - they come from coding rules about sources - and that its best
# estimate undercounts, most of all for small events. From a survey of UCDP's
# own coders they fit, for any reported count and type of violence, a
# distribution of plausible deaths, released as the R package uncertainUCDP.
# The estimate draws each year's civilian deaths from it (civilian_source,
# 00_setup.R; drawn in 10). This step sets it beside the method it replaced, a
# Beta-PERT on UCDP's low, best and high values, and projects the loss under
# whichever of the two is not used, for Table 6 and Table A49.
#
# METHOD
# ------
#   1. Each event's civilian deaths at UCDP's best, as 07 counts them
#      (ucdp_civilian_events(), 00_setup.R), by the reported size of the event.
#   2. uncertainUCDP: 10's draws (2.5th and 97.5th percentiles) and the exact
#      mean. The Beta-PERT: its mean, (low + 4 best + high) / 6, and its 2.5th
#      and 97.5th percentiles.
#   3. The loss at the central projection with the civilians at each method's
#      mean and at the ends of each, every other input as used.
#
# CAVEATS
#   - The model is fitted on all of an event's deaths; the civilian share is
#     applied to the draws, so the civilians are corrected as the event is.
#   - It covers documented events only: deaths that never reach a source (much
#     of Mariupol) stay outside it, as they stay outside UCDP (13p).
#   - Draws are independent between events, so the yearly totals' intervals
#     are narrow: the coders' uncertainty is per event, and it largely cancels
#     in sums of thousands of events.
#   - It is not calibrated to Ukraine: after 2022 most civilian deaths are in
#     events of one or two reported deaths, which it nearly doubles.
#
# INPUTS   data_inter/ukr_civilians_uncertain_ucdp.rds (10), ukr_ucdp_invals.rds (07),
#          data_inter/ucdp_ged_events_slim.rds (07), ukr_ohchr_annual_totals.rds,
#          the static inputs and draws at the central values (mode_projection_inputs())
# OUTPUT   data_inter/ukr_uncertain_ucdp.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

years <- 2022:2025
q <- \(x, p) quantile(x, p, names = FALSE)

# 1. THE EVENTS ================================================================
ev <- ucdp_civilian_events()
by_size <- ev |>
  filter(s > 0) |>
  mutate(size = cut(t, c(0, 1, 2, 5, 10, 100, Inf), labels = c("1", "2", "3-5", "6-10", "11-100", ">100"))) |>
  summarise(events = n(), civilians_best = sum(civ_best), .by = c(year, size)) |>
  mutate(share = civilians_best / sum(civilians_best), .by = year) |>
  arrange(year, size)

# 2. THE TWO METHODS ===========================================================
uu <- read_rds("data_inter/ukr_civilians_uncertain_ucdp.rds")
draws <- uu$draws
stopifnot(identical(as.integer(rownames(draws)), years))
ucdp <- read_rds("data_inter/ukr_ucdp_invals.rds") |> filter(role == "civilians", year %in% years) |> arrange(year)
ohchr <- read_rds("data_inter/ukr_ohchr_annual_totals.rds") |> arrange(year)
compare <-
  tibble(year = years, best = ucdp$dts, pert_min = ucdp$dts_l, pert_mode = ucdp$dts, pert_max = ucdp$dts_u) |>
  mutate(pert_mean = (pert_min + 4 * pert_mode + pert_max) / 6,
         pert_lo = qpert(0.025, pert_min, pert_mode, pert_max, shape = 4),
         pert_hi = qpert(0.975, pert_min, pert_mode, pert_max, shape = 4),
         uu_mean = uu$mean$mean[match(year, uu$mean$year)],
         uu_draw_mean = rowMeans(draws), uu_median = apply(draws, 1, median),
         uu_lo = apply(draws, 1, q, .025), uu_hi = apply(draws, 1, q, .975),
         ratio_mean = uu_mean / pert_mean,
         ohchr = ohchr$dts[match(year, ohchr$year)])
cat("\n=== CIVILIAN DEATHS: THE BETA-PERT ON UCDP'S BOUNDS AGAINST uncertainUCDP ===\n")
print(as.data.frame(compare |> mutate(across(-c(year, ratio_mean), round), ratio_mean = round(ratio_mean, 3))))
print(as.data.frame(by_size), digits = 3)

# 3. THE LOSS ==================================================================
mi <- mode_projection_inputs()
central_civ <- mi$draws |> arrange(year) |> pull(draw_cvs)
used_mean <- if (civilian_source == "uncertainUCDP") compare$uu_mean else compare$pert_mean
stopifnot(isTRUE(all.equal(central_civ, used_mean, tolerance = 1e-9)))
compare <- compare |> mutate(central_projection = central_civ)
loss_with <- function(civ) {
  d <- mi$draws |> mutate(draw_cvs = civ[match(year, years)])
  loss_at_mode(d, mi$static_inputs, mi$pop22_ini) |> select(year, sex, loss)
}
used <- loss_at_mode(mi$draws, mi$static_inputs, mi$pop22_ini) |> select(year, sex, loss_used = loss)
loss <-
  bind_rows(
    loss_with(compare$uu_mean) |> mutate(point = "uncertainUCDP mean"),
    loss_with(compare$uu_lo) |> mutate(point = "uncertainUCDP 2.5th percentile"),
    loss_with(compare$uu_hi) |> mutate(point = "uncertainUCDP 97.5th percentile"),
    loss_with(compare$pert_mean) |> mutate(point = "Beta-PERT mean"),
    loss_with(compare$pert_min) |> mutate(point = "Beta-PERT minimum (UCDP low)"),
    loss_with(compare$pert_max) |> mutate(point = "Beta-PERT maximum (UCDP high)")
  ) |>
  left_join(used, by = c("year", "sex")) |>
  mutate(shift = loss - loss_used)
# the method in use must reproduce the loss as used
used_point <- if (civilian_source == "uncertainUCDP") "uncertainUCDP mean" else "Beta-PERT mean"
alternative <- setdiff(c("uncertainUCDP mean", "Beta-PERT mean"), used_point)
stopifnot(all(abs(loss$shift[loss$point == used_point]) < 1e-9))
cat("\n=== THE LOSS AT THE CENTRAL PROJECTION, CIVILIANS SET AS SHOWN ===\n")
print(as.data.frame(loss |> mutate(across(c(loss, loss_used, shift), \(x) round(x, 3))) |> arrange(point, sex, year)))

write_rds(list(compare = compare, draws = draws, by_size = by_size, loss = loss,
               used = used_point, alternative = alternative, civilian_source = civilian_source,
               n_draw = uu$n_draw, package_version = uu$package_version),
          "data_inter/ukr_uncertain_ucdp.rds")
