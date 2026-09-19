rm(list = ls())
gc()
source("code/00_setup.R")

# data compiled by Olivier Hubert
# in https://www.kaggle.com/datasets/ol4ubert/confirmed-ukrainian-military-personnel-losses
#
# How do the personnel recorded as "missing" resolve? Estimated by following
# everyone listed as missing through the four register releases - v14 (16 Sep
# 2025), v16 (4 Dec 2025), v18 (23 Apr 2026) and v19 (21 Jul 2026) - and
# fitting, to the three windows between them, how fast the missing leave that
# status by months since the event. Each event month's missing in v19 are then
# carried forward from the duration they have reached.
#
# The releases are large individual-level files that are not tracked in git.
# Only the anonymous count summaries below cross the cache boundary, so the
# repository carries no names or dates of birth.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KEY ASSUMPTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How many of the missing are alive drives the combatant total. The model
# projects who will leave the register beyond list maintenance (counted alive
# by the rules below); the prisoners of war among the missing are set by the
# official figures, compared with the prisoners the register records; and a
# share of the rest may be alive for other reasons. alive_evidence() in
# 00_setup.R sets out the evidence, its ranges and the sources. This script
# imputes at the central values - the medians of the evidence's inputs - and
# records the ranges, 10 turns them into the combatant bounds, and 11 draws
# them together with the sampling error of the model and of the
# registration-lag factors (military_draws()).

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# individual-level linkage -> anonymous window counts (cached)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The rules live in 00_setup.R (ual_follow_missing, ual_windows):
#   - a person enters at their first listing as missing, event in 2022-2025
#     and Ukrainian, in v14, v16 or v18, and is looked up in every later
#     release searched in full, whatever the event year or nationality field;
#   - by exact full name and date of birth; if absent, under a corrected key:
#     a key new in that release with the same surname and first name and the
#     same date of birth or event month;
#   - a person no longer listed - absent from a release and from every later
#     one - leaves the register at the first absence. Absent from one release
#     but listed again later, they were still missing. Below, drop-out at the
#     rate at which the dead leave the register is list maintenance and held
#     as still missing; only the excess is resolved alive;
#   - a history stops at its first resolution;
#   - a return from captivity is a resolution to captivity: the person was a
#     prisoner the register had not recorded, and is alive.
linkage <- cache_rds(
  "data_inter/ualosses_window_transitions.rds",
  {
    regs <- map(set_names(ual_releases$path, ual_releases$release), ual_read_release)
    hist <- ual_follow_missing(regs)

    # The dead drop out of the register too, and cannot have been found alive:
    # v14's dead with no trace in v19, by the same two searches, measure how
    # often a record disappears through list maintenance.
    lu14 <- ual_one_per_key(regs$v14)
    lu19 <- ual_one_per_key(regs$v19)
    dead14 <-
      regs$v14 |>
      filter(status == "dead", year %in% 2022:2025, ukrainian) |>
      distinct(key, .keep_all = TRUE)
    dead_absent <- dead14 |> anti_join(lu19, by = "key")
    dead_corrected <- ual_corrected_key(dead_absent, lu19, lu14)

    # The same measured window by window: the dead followed through the
    # releases exactly as the missing are, the two labels swapped, so a dead
    # record that leaves the register shows as "no longer listed".
    swap <- \(r) mutate(r, status = case_when(status == "dead" ~ "missing",
                                               status == "missing" ~ "dead", .default = status))
    hist_dead <- ual_follow_missing(map(regs, swap))

    # Where the people v19 records as returned from captivity were listed
    # before their return: as prisoners, as missing, as dead, or nowhere. By
    # exact key in v14, v16 or v18, else by a corrected key among the records
    # v19 no longer holds.
    lu <- map(regs, ual_one_per_key)
    returned <- lu$v19 |>
      filter(status == "released_prisoner", year %in% 2022:2025, ukrainian) |>
      select(key, name, dob, date_evnt, year)
    earlier <- bind_rows(lu[c("v14", "v16", "v18")])
    exact <- returned |> select(key) |> inner_join(earlier |> select(key, status), by = "key")
    gone <- earlier |> anti_join(lu$v19 |> select(key), by = "key")
    corrected <-
      returned |>
      anti_join(exact, by = "key") |>
      transmute(key, name, dob, month = floor_date(date_evnt, "month")) |>
      inner_join(gone |> transmute(name, dob_e = dob, month_e = floor_date(date_evnt, "month"), status),
                 by = "name", relationship = "many-to-many") |>
      filter((dob == dob_e) %in% TRUE | (month == month_e) %in% TRUE) |>
      select(key, status)
    prior <-
      bind_rows(exact, corrected) |>
      summarise(to = case_when(any(status == "prisoner") ~ "prisoner", any(status == "missing") ~ "missing",
                               any(status == "dead") ~ "dead", .default = "other"), .by = key)

    # the rules: production, and the two alternatives 09 compares with it
    bind_rows(
      ual_windows(hist, "alive", released = "alive") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "windows", rule = "production"),
      ual_windows(hist, "missing", released = "alive") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "windows", rule = "dropouts_missing"),
      ual_windows(hist, "alive", released = "exclude") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "windows", rule = "released_excluded"),
      ual_windows(hist_dead, "alive", released = "exclude") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "dead_windows"),
      returned |>
        left_join(prior, by = "key") |>
        mutate(to = coalesce(to, "not listed")) |>
        count(year, to) |>
        mutate(table = "released_prior"),
      hist |> count(year, entry, release, to = status, found) |> mutate(table = "found"),
      # each person's sequence of statuses, for the documentation's counts of
      # people absent from one release and listed again
      hist |>
        arrange(pid, match(release, ual_releases$release)) |>
        summarise(to = paste(coalesce(status, "absent"), collapse = " > "), .by = c(pid, year, entry)) |>
        count(year, entry, to) |>
        mutate(table = "patterns"),
      dead14 |>
        mutate(found = case_when(!key %in% dead_absent$key ~ "exact",
                                 key %in% dead_corrected$key ~ "corrected",
                                 .default = "absent")) |>
        count(year, found) |>
        mutate(table = "dead_v14_in_v19")
    )
  }
)

windows <- linkage |> filter(table == "windows") |> select(rule, year, month, entry, from_release, to, n)
production_windows <- windows |> filter(rule == "production")

# February 2022 is left out of every fit, as it is out of the registration-lag
# rates (08b): v19 deleted a batch of its records at once, which a duration
# model would read as drop-out at the durations those records had reached.
# Its missing are still projected with everyone else's.
feb22 <- as.Date("2022-02-01")

# List maintenance. The dead cannot have been found alive, yet they leave the
# register too: their drop-out rate, by cohort and window, measures how often
# a record disappears through list maintenance. The missing who leave at that
# rate are held as still missing; only the excess over it is resolved alive.
dead_rates <- function(drop_feb22 = TRUE) {
  linkage |>
    filter(table == "dead_windows", !drop_feb22 | month != feb22) |>
    summarise(at_risk = sum(n), dropped = sum(n[to == "no_longer_listed"]), .by = c(year, from_release)) |>
    mutate(rate = dropped / at_risk)
}
dead_dropout <- dead_rates()
maintain <- function(w, drop_feb22 = TRUE) {
  w <-
    w |>
    filter(!drop_feb22 | month != feb22) |>
    left_join(dead_rates(drop_feb22) |> select(year, from_release, rate), by = c("year", "from_release")) |>
    mutate(maint_share = pmin(1, rate * sum(n) / sum(n[to == "no_longer_listed"])), .by = c(year, from_release)) |>
    mutate(n_maint = if_else(to == "no_longer_listed", n * maint_share, 0))
  bind_rows(w |> mutate(n = n - n_maint),
            w |> filter(n_maint > 0) |> mutate(to = "missing", n = n_maint)) |>
    summarise(n = sum(n), .by = c(rule, year, month, entry, from_release, to))
}

at_risk_v14 <-
  production_windows |>
  filter(from_release == "v14") |>
  summarise(at_risk = sum(n), .by = year)

# twelve-month probabilities in the chain's states; n is prop times the
# cohort's missing in v14
chain_rates <- function(w) {
  ual_compose(w) |>
    mutate(status2 = ual_chain_state(to)) |>
    summarise(prop = sum(prop), .by = c(year, status2)) |>
    left_join(at_risk_v14, by = "year") |>
    transmute(year, status2, n = prop * at_risk, prop, from = "missing")
}

# Stocks by event year, taken from the OUTPUT OF STEP 08 rather than re-read
# from the raw register.
#
# This matters: step 08 redistributes the ~3,300 records whose age or event
# year is not recorded. Counting the raw file again here would silently drop
# them, so the confirmed-death total used in this script would fall short of
# the one used everywhere downstream by those records, and the manuscript
# tables built in 15 would disagree with each other.
stocks_registered <-
  read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") |>
  summarise(n = sum(dx), .by = c(year, status)) |>
  mutate(status = as.character(status)) |>
  arrange(year, status)

# By event month, and completed for registration lag (08b). Each year's
# registered count is spread over its months as v19 lists them, and each month
# is brought to the completeness events reach at four years by its own factor.
# The deaths added are late registrations: people not yet in the register at
# all. Prisoners and released prisoners are left as recorded. `row` points
# into 08b's factor tables, so 11 can apply its resampled factors.
completion <- read_rds("data_inter/ukr_registration_completion.rds")
stock_month <-
  completion$month |>
  mutate(row = row_number()) |>
  mutate(share = n_v19 / sum(n_v19), .by = c(year, status)) |>
  left_join(stocks_registered |> rename(n_year = n), by = c("year", "status")) |>
  transmute(month = m, year, status, row, registered = n_year * share, completed = registered * factor)
stopifnot(isTRUE(all.equal(
  stock_month |> summarise(r = sum(registered), .by = c(year, status)) |> arrange(year, status) |> pull(r),
  stocks_registered |> filter(status %in% c("dead", "missing")) |> arrange(year, status) |> pull(n)
)))
stocks <-
  bind_rows(
    stocks_registered |> filter(!status %in% c("dead", "missing")),
    stock_month |> summarise(n = sum(completed), .by = c(year, status))
  ) |>
  arrange(year, status)

print(stocks |> pivot_wider(names_from = status, values_from = n))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# resolution of the missing over twelve months, by cohort year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# every outcome kept apart, with the people at risk: those listed as missing
# in v14, and those first listed in v16 or v18, who join the later windows
later_entrants <-
  production_windows |>
  filter(entry != "v14", from_release == entry) |>
  summarise(later_entrants = sum(n), .by = year)
resolution_12m <-
  ual_compose(production_windows) |>
  pivot_wider(names_from = to, values_from = prop) |>
  left_join(at_risk_v14, by = "year") |>
  left_join(later_entrants, by = "year") |>
  select(year, at_risk_v14 = at_risk, later_entrants, all_of(ual_outcomes))
print(as.data.frame(resolution_12m |> mutate(across(all_of(ual_outcomes), \(x) round(100 * x, 2)))))

# the same composed by event-year cohort, in the states of the earlier chain,
# kept for the comparison of designs below; list maintenance held as missing,
# as in the model
tasas_long <- chain_rates(maintain(production_windows))
stopifnot(tasas_long |> summarise(p = sum(prop), .by = year) |> with(all(abs(p - 1) < 1e-9)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# resolution by duration since disappearance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hazards of leaving "missing" by months since the event (00_setup.R,
# ual_fit_durations), fitted to one cell per event month and window. The
# projection runs to the horizon 08b settled for registration lag, 48 months:
# the longest duration measured on at least a year of event months. Beyond it
# the data rest on a few months of 2022 events and on v19's one-off deletion
# of a batch of February 2022 records. The model is fitted to the windows with
# list maintenance held as missing and February 2022 left out (above).
fit_model <- function(w, drop_feb22 = TRUE) {
  if (drop_feb22) w <- w |> filter(month != feb22)
  ual_fit_durations(ual_duration_cells(w), horizon = completion$L_ref)
}
model <- fit_model(maintain(production_windows))

hazards <-
  as_tibble(model$h) |>
  mutate(band = paste0(head(model$breaks, -1), "-", model$breaks[-1]), .before = 1)
cat("\n=== MONTHLY HAZARDS OF RESOLUTION BY MONTHS SINCE THE EVENT (%) ===\n")
print(as.data.frame(hazards |> mutate(across(-band, \(x) round(100 * x, 3)))))

# the fit against the counts, by event year and window
fit_check <-
  model$fitted |>
  mutate(year = year(month), at_risk = missing + dead + prisoner + no_longer_listed,
         across(c(p_dead, p_prisoner, p_no_longer_listed), \(p) p * at_risk, .names = "fit_{.col}")) |>
  summarise(across(c(at_risk, dead, prisoner, no_longer_listed, starts_with("fit_")), sum),
            .by = c(year, from_release))
cat("\n=== OBSERVED AND FITTED RESOLUTIONS BY COHORT AND WINDOW ===\n")
print(as.data.frame(fit_check |> mutate(across(where(is.double), round))))

# moves into captivity by window: v19 records past returns in one batch, the
# reason captivity is set by the official figures and not projected
captivity_by_window <-
  production_windows |>
  summarise(n = sum(n[to == "prisoner"]), .by = from_release)
print(captivity_by_window)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the missing who are alive, and the imputation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# impute_missing() (00_setup.R) holds the projection itself, so this call,
# the draws in 11 and the sensitivity steps share one formula.
stock_missing_month <-
  stock_month |>
  filter(status == "missing") |>
  select(month, year, missing_stock = completed)
stock_missing_2026 <-
  stock_missing_month |>
  summarise(missing_stock = sum(missing_stock), .by = year)

# the evidence (00_setup.R sets out the sources). Prisoners are alive whether
# still held or released, so both statuses count as the prisoners the
# register knows of. Where v19's released prisoners had been listed before
# their return gives the share of the unrecorded prisoners among the missing.
register_alive <- sum(stocks$n[stocks$status %in% c("prisoner", "released_prisoner")])
released_prior <- linkage |> filter(table == "released_prior") |> select(year, prior = to, n)
evidence <- alive_evidence(model$resolved, register_alive, released_prior)
print(as.data.frame(evidence))
write_rds(evidence, "data_inter/ukr_alive_missing.rds")

impute <- function(captives = evidence$captives_central, other = evidence$other_central, m = model,
                   stock = stock_missing_month, by = "year") {
  impute_missing(m, stock, captives, other, by)
}
imputation_final <- impute()
print(imputation_final)

# the dead, registered and completed for registration lag; downstream the
# completed count is the "confirmed" military deaths, spread over the
# registered dead's age-sex profile, and 14 and 15 report the late
# registrations within it apart
confirmados_df <-
  stocks |>
  filter(status == "dead") |>
  select(year, confirmados_stock = n) |>
  left_join(
    stocks_registered |> filter(status == "dead") |> select(year, registered = n),
    by = "year"
  ) |>
  mutate(late_registrations = confirmados_stock - registered)

combined_losses <- confirmados_df %>%
  left_join(imputation_final, by = "year") %>%
  mutate(
    total_estimado = confirmados_stock + imputed_dead,
    year = factor(year)
  )

print("=== BALANCE DE PERDIDAS MORTALES TOTALES ESTIMADAS ===")
print(combined_losses)

# What 10 and 11 need to rebuild the military total in any draw
# (military_draws(), 00_setup.R): the registered dead and missing by event
# month, 08b's point and resampled completion factors, the model and the
# evidence. At the point values it must give this script's totals.
military_inputs <- list(
  month = stock_month |> select(month, year, status, row, registered),
  factor_point = completion$month$factor,
  factor_draws = completion$factor_draws,
  model = model,
  evidence = evidence
)
check_point <- military_draws(military_inputs, evidence$captives_central, evidence$other_central)
stopifnot(isTRUE(all.equal(check_point$military, combined_losses$total_estimado)))
write_rds(military_inputs, "data_inter/ukr_military_inputs.rds")

# saved for the manuscript tables assembled in 15
write_rds(
  combined_losses |> mutate(year = as.integer(as.character(year))),
  "data_inter/ukr_ualosses_imputation_table.rds"
)
# the fitted model and the missing by event month, for 13b and 13d
write_rds(model, "data_inter/ukr_ualosses_resolution_model.rds")
write_rds(stock_missing_month, "data_inter/ukr_ualosses_missing_by_month.rds")
write_rds(hazards, "data_inter/ukr_ualosses_resolution_hazards.rds")
# the twelve-month rates by event-year cohort, every outcome kept apart, for
# table A2, and in the earlier chain's states
write_rds(resolution_12m, "data_inter/ukr_ualosses_resolution_12m.rds")
write_rds(tasas_long, "data_inter/ukr_ualosses_transition_rates.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how much the linkage rules, the model and the evidence matter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The military total under alternatives, each at the central values of the
# evidence (the medians of its inputs), each differing from the rules used in
# one respect:
#   - every person no longer listed found alive, list maintenance included;
#   - a person no longer listed held as still missing;
#   - February 2022 events kept in the fit;
#   - a return from captivity leaving the person out of the population at
#     risk, as someone held rather than disappeared;
#   - rates from the people listed as missing in v14 only, without those first
#     listed later;
#   - the projection run to the longest duration observed in the fit
#     instead of 48;
#   - each release window's multipliers alone instead of their average;
#   - every unrecorded prisoner among the missing;
#   - the unrecorded prisoners among the missing in the proportion the
#     returned prisoners of every event year had been;
#   - the register's prisoners with no recorded event year left out, which
#     step 08 spreads over the event years, 2014-2021 included;
#   - captivity as the model projects it, instead of from the official
#     figures;
#   - the captives spread over event years in proportion to the prisoners the
#     register records (held and released) rather than as the model spreads
#     captivity: the same total, a different split by year;
#   - the earlier chain, in which event-year cohorts stand in for duration.
# And how often records vanish: v14's missing and v14's dead with no trace in
# v19 after both searches.
design_total <- function(imp) {
  by_year <-
    confirmados_df |>
    left_join(imp |> select(year, imputed_dead), by = "year") |>
    transmute(year, confirmed = confirmados_stock, total = confirmados_stock + imputed_dead)
  tibble(military = sum(by_year$total), by_year = list(by_year))
}
model_all_alive <- fit_model(production_windows)
model_dropouts_missing <- fit_model(windows |> filter(rule == "dropouts_missing"))
model_with_feb22 <- fit_model(maintain(production_windows, drop_feb22 = FALSE), drop_feb22 = FALSE)
model_released_excluded <- fit_model(maintain(windows |> filter(rule == "released_excluded")))
model_v14_only <- fit_model(maintain(production_windows |> filter(entry == "v14")))
# each window's multipliers alone in place of their average over the year
model_window <- function(w) {
  modifyList(model, list(h = sweep(model$h_reference, 2, model$window_multipliers[w, ], "*")))
}
window_label <- paste0(ual_releases$release[-nrow(ual_releases)], "-", ual_releases$release[-1])

# the register's prisoners and released prisoners with a recorded event year
# in 2022-2025, before 08 spreads the records without one
register_alive_known <-
  read_rds("data_inter/ualosses_counts_2022_2025.rds") |>
  filter(year %in% 2022:2025, status %in% c("prisoner", "released_prisoner")) |>
  pull(dx) |>
  sum()
captives_known_year <- evidence$s_central * (evidence$held_central + returned_military - register_alive_known)

# the share listed as missing among the returned of every event year
share_listed_missing <- evidence$s_min

# captivity as the model projects it: the captives are the projected
# captivity itself
projected_captivity <- sum(imputation_final$projected_captivity)

# the captives by event year in proportion to the register's prisoners, each
# year capped at its projected captivity and never-resolved
prisoners_by_year <-
  stocks |>
  filter(status %in% c("prisoner", "released_prisoner")) |>
  summarise(prisoners = sum(n), .by = year)
by_year_captives <-
  imputation_final |>
  left_join(prisoners_by_year, by = "year") |>
  left_join(confirmados_df |> select(year, confirmados_stock), by = "year") |>
  mutate(captives_year = pmin(projected_captivity + never_resolved,
                              evidence$captives_central * prisoners / sum(prisoners)),
         unresolved_year = projected_captivity + never_resolved - captives_year,
         imputed_dead = projected_dead + (1 - evidence$other_central) * unresolved_year,
         total = confirmados_stock + imputed_dead)

# the earlier chain: its imputed alive are those who leave the register beyond
# list maintenance and the central share of those still missing at its end,
# and the captives are taken out of its dead
cohort_chain <-
  impute_missing_cohorts(evidence$other_central, tasas_long, stock_missing_2026) |>
  mutate(imputed_dead = missing_stock - imputed_alive -
           evidence$captives_central * missing_stock / sum(missing_stock))

cc <- evidence$captives_central
linkage_designs <-
  bind_rows(
    design_total(impute()) |>
      mutate(design = "duration model, drop-outs beyond list maintenance found alive, captives from the official figures (production)",
             captives = cc),
    design_total(impute(m = model_all_alive)) |>
      mutate(design = "every person no longer listed found alive", captives = cc),
    design_total(impute(m = model_dropouts_missing)) |>
      mutate(design = "no longer listed = still missing", captives = cc),
    design_total(impute(m = model_with_feb22)) |>
      mutate(design = "February 2022 events kept in the fit", captives = cc),
    design_total(impute(m = model_released_excluded)) |>
      mutate(design = "returns from captivity left out of the population at risk", captives = cc),
    design_total(impute(m = model_v14_only)) |>
      mutate(design = "people listed as missing in v14 only", captives = cc),
    design_total(impute(m = modifyList(model, list(horizon = max(model$fitted$d1))))) |>
      mutate(design = "duration model to the longest duration observed", captives = cc),
    map_dfr(seq_along(window_label), \(w) {
      design_total(impute(m = model_window(w))) |>
        mutate(design = paste0("the ", window_label[w], " window's multipliers alone"), captives = cc)
    }),
    design_total(impute(captives = evidence$unrecorded_central)) |>
      mutate(design = "every unrecorded prisoner among the missing", captives = evidence$unrecorded_central),
    design_total(impute(captives = share_listed_missing * evidence$unrecorded_central)) |>
      mutate(design = "unrecorded prisoners among the missing as the returned of every year had been",
             captives = share_listed_missing * evidence$unrecorded_central),
    design_total(impute(captives = captives_known_year)) |>
      mutate(design = "the register's prisoners with no recorded event year left out",
             captives = captives_known_year),
    design_total(impute(captives = projected_captivity)) |>
      mutate(design = "captivity as the model projects it", captives = projected_captivity),
    tibble(military = sum(by_year_captives$total),
           by_year = list(by_year_captives |> transmute(year, confirmed = confirmados_stock, total)),
           design = "captives spread over event years by the register's prisoners",
           captives = cc),
    design_total(cohort_chain) |>
      mutate(design = "event-year cohorts standing in for duration", captives = cc)
  ) |>
  relocate(design)
stopifnot(isTRUE(all.equal(linkage_designs$military[1], sum(combined_losses$total_estimado))))

no_trace <-
  bind_rows(
    linkage |> filter(table == "found", entry == "v14", release == "v19") |> mutate(who = "missing"),
    linkage |> filter(table == "dead_v14_in_v19") |> mutate(who = "dead")
  ) |>
  summarise(listed_v14 = sum(n), corrected_key = sum(n[found == "corrected"]),
            no_trace = sum(n[found == "absent"]), .by = c(who, year)) |>
  mutate(pct_no_trace = 100 * no_trace / listed_v14)

patterns <-
  linkage |>
  filter(table == "patterns") |>
  select(year, entry, pattern = to, n) |>
  mutate(first_resolution = map_chr(strsplit(pattern, " > "), function(s) {
    r <- s[!s %in% c("missing", "absent")]
    if (length(r) > 0) r[1] else if (tail(s, 1) == "absent") "no_longer_listed" else "missing"
  }))
# people whose first resolution is a return from captivity, by cohort and
# list of entry
released_first <-
  patterns |>
  filter(first_resolution == "released_prisoner") |>
  summarise(n = sum(n), .by = c(year, entry))

print(as.data.frame(linkage_designs |> select(-by_year)))
print(as.data.frame(no_trace))
print(as.data.frame(released_first |> pivot_wider(names_from = entry, values_from = n, values_fill = 0)))
cat("\n=== THE DEAD'S DROP-OUT BY WINDOW (%) ===\n")
print(as.data.frame(dead_dropout |> mutate(rate = round(100 * rate, 2))))
cat("\n=== WHERE THE RETURNED PRISONERS WERE LISTED BEFORE THEIR RETURN ===\n")
print(as.data.frame(released_prior |> pivot_wider(names_from = prior, values_from = n, values_fill = 0)))
cat(sprintf("share listed as missing, of those not recorded as prisoners: %.3f (all years), %.3f (2024-2025)\n",
            evidence$s_min, evidence$s_mode))
cat(sprintf("register prisoners and released, 2022-2025: %s with the unknown-year records spread, %s without\n",
            scales::comma(round(register_alive)), scales::comma(register_alive_known)))
write_rds(list(designs = linkage_designs, no_trace = no_trace, patterns = patterns,
               released_first = released_first, dead_dropout = dead_dropout,
               released_prior = released_prior, captivity_by_window = captivity_by_window,
               by_year_captives = by_year_captives, register_alive_known = register_alive_known,
               window_multipliers = model$window_multipliers),
          "data_inter/ukr_ualosses_linkage_checks.rds")

# military deaths by event month at the central values, for 15's comparison
# with official statements made at given dates: the registered dead, the dead
# completed for registration lag, and the missing imputed as dead
military_by_month <-
  stock_month |>
  filter(status == "dead") |>
  select(month, year, registered_dead = registered, completed_dead = completed) |>
  full_join(impute(by = c("month", "year")) |> select(month, year, missing = missing_stock, imputed_dead),
            by = c("month", "year")) |>
  mutate(across(c(registered_dead, completed_dead, missing, imputed_dead), \(x) coalesce(x, 0)),
         military = completed_dead + imputed_dead) |>
  arrange(month)
stopifnot(isTRUE(all.equal(sum(military_by_month$military), sum(combined_losses$total_estimado))))
write_rds(military_by_month, "data_inter/ukr_military_by_month.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rescale the confirmed-dead age-sex profile up to the imputed total
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ual <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")

adj_fct <-
  ual |>
  summarise(dts = sum(dx), .by = c(year, status)) |>
  filter(status == "dead") |>
  left_join(
    combined_losses |>
      select(year, total_estimado) |>
      mutate(year = as.character(year) |> as.integer()),
    by = "year"
  ) |>
  mutate(adj = total_estimado / dts) |>
  select(year, adj)

print(adj_fct)

ual2 <-
  ual |>
  filter(status == "dead") |>
  left_join(adj_fct, by = "year") |>
  mutate(dx = dx * adj) |>
  select(-adj)

ual2 |>
  summarise(dts = sum(dx), .by = c(year, status))

write_rds(
  ual2,
  "data_inter/ukr_ualosses_conflict_deaths_imputed_miss_sex_age_2022_2025.rds"
)

copy_this(
  combined_losses |>
    select(
      year,
      confirmed_deaths = confirmados_stock,
      missing = missing_stock,
      imputed_dead,
      imputed_alive,
      total_deaths = total_estimado
    )
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt_plot <-
  combined_losses |>
  select(
    year,
    confirmed_dead = confirmados_stock,
    imputed_dead,
    imputed_alive
  ) |>
  gather(-year, key = "status", value = "dts") |>
  mutate(
    status = factor(
      status,
      levels = c("imputed_alive", "imputed_dead", "confirmed_dead")
    )
  )

cols <- c("grey30", "#e2606b", "#e63946")

dt_plot |>
  ggplot() +
  geom_bar(
    aes(fill = status, y = dts, x = year),
    position = "stack",
    stat = "identity"
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma, breaks = seq(0, 70000, 10000)) +
  labs(y = "Counts", x = "Year", fill = "Status") +
  theme_minimal()
ggsave("figures/exploratory/labtalk/missing_time_qt_imputed.png", w = 6, h = 3)

dt_plot |>
  mutate(
    prop = dts / sum(dts),
    pct_label = scales::percent(prop, accuracy = 0.1),
    .by = year,
  ) |>
  ggplot(aes(x = factor(year), y = dts, fill = status)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(
    aes(label = pct_label),
    position = position_fill(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage", x = "Year", fill = "Status") +
  theme_minimal()
ggsave("figures/exploratory/labtalk/missing_time_imputed.png", w = 6, h = 3)

dt_plot |>
  filter(status != "confirmed_dead") |>
  mutate(
    prop = dts / sum(dts),
    pct_label = scales::percent(prop, accuracy = 0.1),
    .by = year,
  ) |>
  ggplot(aes(x = factor(year), y = dts, fill = status)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(
    aes(label = pct_label),
    position = position_fill(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage", x = "Year", fill = "Status") +
  theme_minimal()
ggsave("figures/exploratory/labtalk/missing_time_imputed_v2.png", w = 6, h = 3)
