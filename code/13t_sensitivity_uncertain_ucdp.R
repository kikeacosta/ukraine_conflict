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
# whichever of the two is not used, for Table 6 and Table A14.
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
#   4. The loss with every death from UCDP: the combatants too drawn from
#      uncertainUCDP, in place of the register (section 4).
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
# OUTPUT   data_inter/ukr_uncertain_ucdp.rds: compare, draws, by_size and loss for
#          the civilians; combatants, combatant_draws and ucdp_only for the
#          scenario with every death from UCDP, and ucdp_only_parts for its
#          shift by part (civilians alone, combatants alone, both)
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

# 4. EVERY DEATH FROM UCDP =====================================================
# The estimate takes the civilians from UCDP and the combatants from the
# register. Here both come from UCDP, each event's deaths drawn from
# uncertainUCDP's own distributions: what a study with that one source would
# find. Ukraine's combatants are the deaths 07 assigns to them, event by event
# (ucdp_combatant_events(), 00_setup.R), drawn as the civilians are.
#   - UCDP gives no ages, so the combatants keep the register's: in each year the
#     same mix of the profiles of its dead and of its missing as in the central
#     projection. The scenario changes the totals and nothing else.
#   - The share of the missing counted dead and the registration lag are
#     properties of the register and have no part in it.
#   - The two groups are drawn separately; the rows at the percentiles set both
#     at the same percentile, which is wider than their joint 95% interval.
#   - A second row leaves uncertainUCDP out: the combatants at UCDP's best value
#     as 07 adjusts it, the civilians as the estimate has them.
cmb_ev <- ucdp_combatant_events()
cmb_draws <- uncertain_ucdp_draws(cmb_ev, uu$n_draw)
cmb_mean <- uncertain_ucdp_mean(cmb_ev)
stopifnot(identical(as.integer(rownames(cmb_draws)), years), identical(as.integer(cmb_mean$year), years),
          all(abs(rowMeans(cmb_draws) / cmb_mean$mean - 1) < 0.01))
ucdp_cmb <- read_rds("data_inter/ukr_ucdp_invals.rds") |> filter(role == "combatants", year %in% years) |> arrange(year)
register <- mi$param_table |> filter(role == "combatants") |> arrange(year)
central_cmb <- mi$draws |> arrange(year)
stopifnot(identical(as.integer(register$year), years), identical(as.integer(central_cmb$year), years))
combatants <-
  tibble(year = years, ucdp_best = ucdp_cmb$dts, ucdp_low = ucdp_cmb$dts_l, ucdp_high = ucdp_cmb$dts_u,
         uu_mean = cmb_mean$mean, uu_lo = apply(cmb_draws, 1, q, .025), uu_hi = apply(cmb_draws, 1, q, .975),
         register_bound_low = register$min, register_central = central_cmb$draw_cmb,
         register_bound_high = register$max) |>
  mutate(ratio_mean = uu_mean / ucdp_best)
cat("\n=== UKRAINE'S COMBATANTS: UCDP, uncertainUCDP AND THE REGISTER ===\n")
print(as.data.frame(combatants |> mutate(across(-c(year, ratio_mean), round), ratio_mean = round(ratio_mean, 3))))

# the combatants replaced, with the confirmed part scaled so that the two age
# profiles keep the central projection's weights
loss_ucdp_only <- function(civ, cmb) {
  d <- mi$draws |>
    mutate(draw_cvs = civ[match(year, years)],
           conf_cmb = conf_cmb / draw_cmb * cmb[match(year, years)],
           draw_cmb = cmb[match(year, years)])
  loss_at_mode(d, mi$static_inputs, mi$pop22_ini) |> select(year, sex, loss)
}
ucdp_only_points <- c("UCDP only, uncertainUCDP mean", "UCDP only, uncertainUCDP 2.5th percentile",
                      "UCDP only, uncertainUCDP 97.5th percentile", "UCDP only, combatants at UCDP best")
ucdp_only <-
  bind_rows(
    loss_ucdp_only(compare$uu_mean, combatants$uu_mean) |>
      mutate(point = ucdp_only_points[1], combatants = sum(combatants$uu_mean)),
    loss_ucdp_only(compare$uu_lo, combatants$uu_lo) |>
      mutate(point = ucdp_only_points[2], combatants = sum(combatants$uu_lo)),
    loss_ucdp_only(compare$uu_hi, combatants$uu_hi) |>
      mutate(point = ucdp_only_points[3], combatants = sum(combatants$uu_hi)),
    loss_ucdp_only(central_civ, combatants$ucdp_best) |>
      mutate(point = ucdp_only_points[4], combatants = sum(combatants$ucdp_best))
  ) |>
  left_join(used, by = c("year", "sex")) |>
  mutate(shift = loss - loss_used, combatants_used = sum(central_cmb$draw_cmb))
# with the register's own totals put back, the function must give the loss as used
stopifnot(max(abs(loss_ucdp_only(central_civ, central_cmb$draw_cmb)$loss - used$loss_used)) < 1e-9)
cat("\n=== THE LOSS WITH EVERY DEATH FROM UCDP ===\n")
print(as.data.frame(ucdp_only |> mutate(across(c(loss, loss_used, shift), \(x) round(x, 3)),
                                        across(starts_with("combatants"), round)) |> arrange(point, sex, year)))

# The scenario in its two parts, each group changed alone at uncertainUCDP's
# mean with the other as the estimate has it, and the two together. The
# civilians alone are section 3's row; what the parts leave of the whole is
# their interaction, through the population the two share.
ucdp_only_parts <-
  bind_rows(
    loss |> filter(point == "uncertainUCDP mean") |> transmute(year, sex, loss, part = "civilians alone"),
    loss_ucdp_only(central_civ, combatants$uu_mean) |> mutate(part = "combatants alone"),
    ucdp_only |> filter(point == ucdp_only_points[1]) |> transmute(year, sex, loss, part = "both")
  ) |>
  left_join(used, by = c("year", "sex")) |>
  mutate(shift = loss - loss_used) |>
  select(year, sex, part, shift) |>
  pivot_wider(names_from = part, values_from = shift) |>
  mutate(interaction = both - `civilians alone` - `combatants alone`) |>
  left_join(used, by = c("year", "sex")) |>
  arrange(sex, year)
# the parts must account for the whole to within a hundredth of a year
stopifnot(all(abs(ucdp_only_parts$interaction) < 0.01))
cat("\n=== EVERY DEATH FROM UCDP: THE SHIFT IN THE LOSS BY PART ===\n")
print(as.data.frame(ucdp_only_parts |> mutate(across(-c(year, sex), \(x) round(x, 3)))))

write_rds(list(compare = compare, draws = draws, by_size = by_size, loss = loss,
               combatants = combatants, combatant_draws = cmb_draws, ucdp_only = ucdp_only,
               ucdp_only_parts = ucdp_only_parts,
               used = used_point, alternative = alternative, civilian_source = civilian_source,
               n_draw = uu$n_draw, package_version = uu$package_version),
          "data_inter/ukr_uncertain_ucdp.rds")
