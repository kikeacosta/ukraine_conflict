# ==============================================================================
# STEP 08b - Registration lag: how complete is each event month in v19?
# ==============================================================================
#
# WHY
# ---
# The register keeps adding deaths and disappearances long after the event:
# for 2025 events, v14 listed 9,001 dead and v19, nine months on, 16,457, and
# even 2022 events keep growing. The imputation in 09 only reclassifies people
# already listed, so a loss not yet registered would not be counted. This step
# measures how far each event month's v19 count falls short of the completeness
# older events have reached; 09 scales the v19 dead and missing up by that
# factor before the imputation.
#
# METHOD - a chain-ladder on the six releases
# --------------------------------------------
# For each event month m, each status s (dead, missing) and each consecutive
# pair of releases r -> r+1, the growth of the register:
#   dead      records listed as dead in r+1 whose key is absent from r in any
#             status, less those that are a missing person of r under a
#             corrected key, less records listed as dead in r whose key is
#             absent from r+1 in any status
#   missing   new persons listed as missing in r+1, a record whose key is
#             absent from every release up to r and is not a corrected key of
#             a record those releases held
#
# The two statuses are treated differently on purpose. The dead keep the
# SYMMETRIC net, additions less drops, so that a record which leaves and returns
# cancels; only a missing person recoded as dead under a corrected key is taken
# out of the additions, since that is a resolution and not a new death. The
# missing count new persons only, because a missing person who drops off the
# register is a resolution that 09 counts as found alive, and netting the drops
# here would take them out twice.
#
# Applying the corrected-key test to the dead's additions against ALL departed
# records breaks this: a re-keyed death is then dropped without being re-added
# and nets to -1, which drives the growth rates negative and the completion
# factors below one - a register that loses deaths as it ages. The guards below
# stop the asymmetry recurring.
# A key present in both releases is left out whatever its status, so a missing
# person who is found dead, taken prisoner or released is a resolution, not a
# registration or a de-registration. The dead leave the register only through
# list maintenance - corrections and merged duplicates - so their drops are
# netted out. The missing are different: a missing person who drops off the
# register is a resolution in 09's model, counted there as found alive, so
# netting the drops here as well would take them out twice. Their growth is
# therefore new persons only; a record re-keyed, re-spelt or listed again after
# a gap is not new. The monthly growth rate
#   lambda = log(1 + growth / N_r) / (months between the releases)
# is attached to the lag since the event (midpoint of the pair), and pooled
# by lag band, weighting by N_r x months. Lags and months run between the dates
# of the releases, each dated by the data it holds (ual_releases, 00_setup.R).
#
# Chaining those rates from each month's lag at v19 up to 48 months - the
# longest lag measured on at least twelve event months - gives the factor that
# brings the month to the completeness events reach at four years. Beyond 48
# months the rates rest on three event months of a single pair of releases,
# where v19 deleted a batch of February 2022 missing records in one clean-up:
# a one-off, not a lag. The register still grows after four years, so the
# horizon makes the factor a lower bound on the lag; 13d shows how far a longer
# horizon would move the results.
#
# WHAT THE CHAIN ASSUMES, AND ITS ERROR
# -------------------------------------
# The chain gives a recent event month the growth older months showed at the
# lags it has yet to reach: the register is taken to keep growing by lag as it
# did over these releases. That is a model, measured on five windows spanning
# nine months, and its error is the error of resting on five windows: what a
# release adds differs far more from one window to the next - a batch arrives
# or it does not - than between event months within a window. So the factors
# are rebuilt under weights on the windows (lag_factor_draws(), 00_setup.R),
# and 11 draws the weights once per simulation for the factors and for the
# bound on the missing alive together (window_weight_draws()).
#
# TWO ALTERNATIVES, KEPT FOR 13d
# ------------------------------
#   exposure   The rates above are per record of the status's own stock. The
#              missing deplete as they are resolved and the dead accrue, so an
#              older cohort's rate per missing person is measured on a stock
#              that has shrunk, and a recent cohort, most of whose listed
#              losses are still missing, is given too much growth of the
#              missing and too little of the dead. The alternative takes both
#              statuses' growth per person the register lists for the event
#              month, dead or missing, and chains the listed stock.
#   dating     The releases dated by their Kaggle versions, not by the data
#              they hold.
#
# INPUT   the six register releases (gitignored), through a cache of
#         aggregate counts: data_inter/ualosses_registration_by_month.rds
# OUTPUT  data_inter/ukr_registration_completion.rds: the growth rates by lag
#         band and by window, the horizon, v19's dead and missing by event
#         month with each month's completion factor, and the factors under the
#         two alternatives
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

releases <- ual_releases

# 1. AGGREGATE COUNTS FROM THE RELEASES (cached; no names leave this block) ====
reg_counts <- cache_rds("data_inter/ualosses_registration_by_month.rds", {
  regs <- map(set_names(releases$path, releases$release), ual_read_release)
  with_np <- \(r) mutate(r, m = floor_date(date_evnt, "month"),
                         np = str_squish(str_remove(key, "\\s+\\S+$")))
  rd <- function(r) {
    with_np(r) |>
      filter(ukrainian, !is.na(m), m >= as.Date("2022-02-01"), m <= as.Date("2025-12-01")) |>
      ual_one_per_key() |>
      select(key, np, name, surname, first, dob, date_evnt, m, status)
  }
  d <- map(regs, rd)
  losses <- \(x) filter(x, status %in% c("dead", "missing"))

  # one long table: stocks by release, and additions and drops between
  # consecutive releases. A record counts as added or dropped only if its key
  # is absent from the other release in ANY status, so a missing person who
  # becomes a prisoner or is released is a resolution, not a de-registration.
  # added_new: the missing who are new persons (see METHOD). The corrected-key
  # test looks among the earlier records whose key the later release no
  # longer holds, as the linkage in 09 does.
  stock <- imap_dfr(d, \(x, r) count(losses(x), m, status, name = "n") |>
                      mutate(kind = "stock", from = r, to = NA_character_))
  flows <- map_dfr(seq_len(length(d) - 1), function(i) {
    a <- d[[i]]; b <- d[[i + 1]]
    earlier <- bind_rows(regs[1:i]) |>
      with_np() |>
      select(key, np, name, surname, first, dob, date_evnt, m, status) |>
      distinct()
    gone <- earlier |> anti_join(regs[[i + 1]] |> distinct(key), by = "key")
    # The corrected-key test (ual_match_keys(), 00_setup.R), applied to each
    # status by what it is for.
    #
    # For the MISSING it removes anyone already in the register under another
    # key, so that "new persons" are new.
    #
    # For the DEAD it must remove only the missing who became dead under a
    # corrected key: those are resolutions, not newly registered deaths. It must
    # NOT remove a dead record that was merely re-keyed, because that record's
    # old key is counted among the drops, and excluding the new key as well
    # would net it to -1 - which is what makes the growth rates go negative and
    # the completion factors fall below one. A register cannot lose deaths as it
    # ages.
    # The missing keep the linkage's own test against every departed record.
    # The dead look only among departed MISSING records and, for the rule that
    # compares names, use the patronymic as well, because in a register of this
    # size a surname, first name and event month coincide often: matched loosely
    # against all departures that rule removes 848, 105, 266, 92 and 494 of the
    # new dead over the five pairs, against one in all once the patronymic is
    # required - so the loose version is matching namesakes, not the same
    # person.
    drop_keys <- function(x, g, k) x |> anti_join(ual_match_keys(x, g, k), by = "key")
    new_missing <- drop_keys(losses(b) |> filter(status == "missing") |>
                               anti_join(earlier, by = "key"), gone, "name")
    # The dead keep the symmetric net, added - dropped, so that a record which
    # leaves and returns cancels. Only the missing who became dead under a
    # corrected key are taken out of the additions, counted here so that the
    # growth can subtract them.
    #
    # A record is added only if its key is absent from the WHOLE earlier
    # release, and dropped only if absent from the whole later one - not from
    # `a` and `b`, which hold the dated records of 2022-2025 alone. The register
    # lists thousands of deaths with no event date (4,199 in v14, 3,151 in v19)
    # and dates them as it goes: tested against `a`, a death undated in one
    # release and dated in the next counts as newly registered (110, 103, 527,
    # 186 and 1,161 of the additions over the five pairs), although 08 has already spread
    # the undated deaths over the event years, so the base the factors multiply
    # holds it. What this leaves out is the death that ENTERS the register
    # undated, which has no event month to be counted in: 421, 266, 538, 228 and 529 such
    # records against 725, 63, 247, 140 and 112 undated ones dropped, a net 695 over the
    # same releases beside the 5,901 measured on them, so the factors remain a
    # slight lower bound on that account.
    in_a <- regs[[i]] |> distinct(key)
    in_b <- regs[[i + 1]] |> distinct(key)
    added_dead <- losses(b) |> filter(status == "dead") |> anti_join(in_a, by = "key")
    mtd <- added_dead |>
      anti_join(drop_keys(added_dead, gone |> filter(status == "missing"), "np"), by = "key")
    bind_rows(
      losses(b) |> anti_join(in_a, by = "key") |> count(m, status, name = "n") |> mutate(kind = "added"),
      losses(a) |> anti_join(in_b, by = "key") |> count(m, status, name = "n") |> mutate(kind = "dropped"),
      new_missing |> count(m, status, name = "n") |> mutate(kind = "added_new"),
      mtd |> count(m, status, name = "n") |> mutate(kind = "added_recoded")
    ) |> mutate(from = releases$release[i], to = releases$release[i + 1])
  })
  # the latest event each release lists and its number of records, to check the
  # dates the releases are given (ual_releases)
  latest <- imap_dfr(regs, \(x, r) tibble(m = max(x$date_evnt, na.rm = TRUE), status = NA_character_,
                                           n = nrow(x), kind = "latest_event", from = r, to = NA_character_))
  bind_rows(stock, flows, latest)
})

stock <- reg_counts |> filter(kind == "stock") |> transmute(m, status, n, release = from)
flows <- reg_counts |> filter(!kind %in% c("stock", "latest_event")) |> select(m, status, n, kind, from, to)

# the releases are dated from the latest event each lists: the files must agree
latest <- reg_counts |> filter(kind == "latest_event")
stopifnot(identical(latest$m[match(releases$release, latest$from)], releases$latest_event))
cat("\n=== THE RELEASES: LATEST EVENT, DATE GIVEN AND KAGGLE VERSION ===\n")
print(as.data.frame(releases |> transmute(release, records = latest$n[match(release, latest$from)],
                                          latest_event, date, version,
                                          days_to_version = as.integer(version - latest_event))))

# 2. GROWTH RATES BY LAG ======================================================
lag_months <- function(release_date, m) as.numeric(release_date - (m + 14)) / 30.44
rel_date <- set_names(releases$date, releases$release)
last_release <- releases$release[nrow(releases)]
windows <- tibble(from = head(releases$release, -1), to = releases$release[-1])
bands <- c(1, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48, 56)

# One row per event month, status and window: the stock at the window's start
# and the growth over it, with the lags the dates of the releases give. Every
# month and status the earlier release lists has a row, a window that brought
# it nothing included: its growth is zero, not missing. n_listed and
# growth_listed are the two statuses together, the exposure of the alternative.
growth_cells <- function(dates) {
  stock |>
    rename(from = release, n_from = n) |>
    inner_join(windows, by = "from") |>
    left_join(flows |> pivot_wider(names_from = kind, values_from = n, values_fill = 0),
              by = c("m", "status", "from", "to")) |>
    mutate(
      across(c(added, dropped, added_new, added_recoded), \(x) coalesce(x, 0)),
      L0 = lag_months(dates[from], m), L1 = lag_months(dates[to], m),
      dL = L1 - L0, Lmid = (L0 + L1) / 2,
      growth = if_else(status == "missing", added_new, added - added_recoded - dropped)
    ) |>
    mutate(n_listed = sum(n_from), growth_listed = sum(growth), .by = c(m, from))
}

# The monthly rate of each cell and its weight, by lag band. exposure "own":
# growth per record of the status's own stock. "listed": growth per person
# listed for the event month, dead or missing - the listed stock's own rate,
# log(1 + all growth / listed) / months, split between the statuses as the
# growth is.
# February 2022 holds five days of events and takes the register's batch
# operations on the start of the war - v18 added 221 missing to it at once,
# v19 deleted 335 - so it is left out of the rates. At its lag it is past
# the horizon and needs no factor itself.
band_cells <- function(g, exposure = c("own", "listed")) {
  exposure <- match.arg(exposure)
  g <-
    if (exposure == "own") {
      g |> mutate(n_exp = n_from, lambda = log1p(growth / n_from) / dL)
    } else {
      g |> mutate(n_exp = n_listed, x = growth_listed / n_listed,
                  lambda = growth / n_listed / dL * if_else(abs(x) < 1e-12, 1, log1p(x) / x))
    }
  g |>
    filter(L0 >= 1, n_exp >= 50, m != as.Date("2022-02-01")) |>
    mutate(band = cut(Lmid, bands, right = TRUE), wt = n_exp * dL) |>
    filter(!is.na(band))
}
band_rates <- function(gb) {
  gb |>
    summarise(lambda = weighted.mean(lambda, wt), n_obs = n(), .by = c(status, band)) |>
    mutate(lo = bands[as.integer(band)], hi = bands[as.integer(band) + 1]) |>
    arrange(status, lo)
}

growth_band <- band_cells(growth_cells(rel_date))
lambda_band <- band_rates(growth_band)

# Guard against the asymmetry that broke this step once: the register cannot
# lose deaths as it ages, so no band rate may be negative. The tolerance is for
# a band where additions and drops balance exactly and the weighted mean lands a
# few parts in ten million below zero.
stopifnot(all(lambda_band$lambda >= -1e-6))

cat("\n=== MONTHLY GROWTH OF THE REGISTER BY LAG SINCE THE EVENT ===\n")
print(as.data.frame(lambda_band |> mutate(lambda_pct = round(100 * lambda, 3)) |>
                      select(status, band, lambda_pct, n_obs)))

# 3. CHAIN TO THE LONGEST WELL-MEASURED LAG ===================================
# the last band with at least a year's worth of event months (12) behind it
L_ref <- lambda_band |>
  summarise(ok = all(n_obs >= 12), .by = c(lo, hi)) |>
  filter(ok) |>
  pull(hi) |>
  max()
stopifnot(L_ref == 48)

# growth from lag L_from up to lag L_to, chaining the band rates
chain_factor <- function(L_from, s, L_to = L_ref, lb = lambda_band) {
  lb <- lb |> filter(status == s)
  exp(sum(pmax(0, pmin(lb$hi, L_to) - pmax(lb$lo, L_from)) * lb$lambda))
}

completion_month <-
  stock |>
  filter(release == last_release) |>
  mutate(year = year(m), lag = lag_months(rel_date[[last_release]], m),
         factor = map2_dbl(lag, status, chain_factor)) |>
  select(m, year, status, lag, n_v19 = n, factor)
stopifnot(all(completion_month$factor >= 1 - 1e-9))   # a lag correction can only add

by_year <- function(month, f) {
  month |>
    mutate(f = f) |>
    summarise(n_completed = sum(n_v19 * f), n_v19 = sum(n_v19), .by = c(year, status)) |>
    mutate(ratio = n_completed / n_v19) |>
    arrange(status, year)
}
completion_year <- by_year(completion_month, completion_month$factor)

cat("\n=== V19 COUNTS BROUGHT TO THE COMPLETENESS EVENTS REACH AT FOUR YEARS ===\n")
print(as.data.frame(completion_year |> mutate(across(c(n_v19, n_completed), round),
                                              ratio = round(ratio, 4))))

# 4. THE ERROR OF THE FACTORS: THE RELEASE WINDOWS ============================
# Each band's rate is a sum over the windows of rate x weight over a sum of
# weight. Kept by window, the sums rebuild the rates, and with them the factors,
# under any weights on the windows (lag_factor_draws(), 00_setup.R): one set of
# weights per simulation in 11, the same that weigh the windows behind the
# bound on the missing alive. Here, what the windows say on their own: each
# left out in turn, and the spread under the weights the simulation draws.
window_cells <-
  growth_band |>
  summarise(num = sum(lambda * wt), den = sum(wt), .by = c(status, band, from)) |>
  mutate(lo = bands[as.integer(band)], hi = bands[as.integer(band) + 1]) |>
  arrange(status, lo, match(from, releases$release))
lag_inputs <- list(cells = window_cells, month = completion_month |> select(status, lag), L_ref = L_ref)
stopifnot(isTRUE(all.equal(as.vector(lag_factor_draws(lag_inputs)), completion_month$factor)))

completed_under <- function(f) {
  s <- completion_month$status
  tibble(dead = colSums(completion_month$n_v19 * f * (s == "dead")),
         missing = colSums(completion_month$n_v19 * f * (s == "missing")))
}
nw <- nrow(windows)
cat("\n=== COMPLETED COUNTS WITH EACH WINDOW LEFT OUT ===\n")
print(as.data.frame(
  bind_cols(tibble(left_out = c("none", paste(windows$from, windows$to, sep = "-"))),
            completed_under(lag_factor_draws(lag_inputs, rbind(1, 1 - diag(nw))))) |>
    mutate(across(c(dead, missing), round))
))
set.seed(20260919)
B_lag <- 2000
window_spread <- completed_under(lag_factor_draws(lag_inputs, window_weight_draws(B_lag)))
cat("\n=== COMPLETED COUNTS UNDER", B_lag, "DRAWS OF THE WINDOW WEIGHTS ===\n")
print(as.data.frame(
  window_spread |>
    pivot_longer(everything(), names_to = "status", values_to = "completed") |>
    summarise(sd = sd(completed), lo = quantile(completed, 0.025), mid = median(completed),
              hi = quantile(completed, 0.975), .by = status) |>
    mutate(across(-status, round))
))

# 5. THE ALTERNATIVES =========================================================
# exposure: both statuses' growth per person listed for the event month. The
# listed stock N of a month grows at the two rates together; over a band it
# grows by exp(rate x months), and each status takes its rate's share of what
# is added. The factor is the completed count of the status over v19's.
gb_listed <- band_cells(growth_cells(rel_date), "listed")
lambda_listed <- band_rates(gb_listed)
stopifnot(all(lambda_listed$lambda >= -1e-6))
factor_listed <- function(month, lb, L_to = L_ref) {
  bd <- lb |>
    filter(hi <= L_to) |>
    select(status, lo, hi, lambda) |>
    pivot_wider(names_from = status, values_from = lambda, values_fill = 0) |>
    arrange(lo)
  added <- map2_dbl(month$lag, month$status, function(l, s) {
    N <- 1
    add <- 0
    for (b in seq_len(nrow(bd))) {
      d <- max(0, min(bd$hi[b], L_to) - max(bd$lo[b], l))
      total <- bd$dead[b] + bd$missing[b]
      add <- add + N * bd[[s]][b] * (if (abs(total) < 1e-12) d else expm1(total * d) / total)
      N <- N * exp(total * d)
    }
    add
  })
  month |>
    mutate(listed = sum(n_v19), .by = m) |>
    mutate(f = 1 + added * listed / n_v19) |>
    pull(f)
}

# dating: the same chain with the releases dated by their Kaggle versions
rel_version <- set_names(releases$version, releases$release)
lambda_version <- band_rates(band_cells(growth_cells(rel_version)))
factor_version <- map2_dbl(lag_months(rel_version[[last_release]], completion_month$m), completion_month$status,
                           \(l, s) chain_factor(l, s, lb = lambda_version))

alternatives <- bind_rows(
  completion_month |> transmute(m, status, alternative = "exposure", f = factor_listed(completion_month, lambda_listed)),
  completion_month |> transmute(m, status, alternative = "dating", f = factor_version)
)
stopifnot(all(alternatives$f >= 1 - 1e-9))

cat("\n=== MONTHLY GROWTH PER PERSON LISTED, DEAD OR MISSING ===\n")
print(as.data.frame(lambda_listed |> mutate(lambda_pct = round(100 * lambda, 3)) |>
                      select(status, band, lambda_pct, n_obs)))
cat("\n=== COMPLETED COUNTS: AS USED, AND UNDER THE TWO ALTERNATIVES ===\n")
print(as.data.frame(
  bind_rows(
    completion_year |> mutate(alternative = "as used"),
    alternatives |> nest(.by = alternative) |>
      mutate(res = map(data, \(a) by_year(completion_month, a$f))) |> select(-data) |> unnest(res)
  ) |>
    summarise(n_v19 = sum(n_v19), n_completed = round(sum(n_completed)), .by = c(alternative, status))
))

# growth: the rates the bands pool, one row per event month, status and pair of
# releases, for 13l to fit a curve through
write_rds(
  list(lambda_band = lambda_band, L_ref = L_ref, month = completion_month, year = completion_year,
       cells = window_cells, alternatives = alternatives,
       lambda_listed = lambda_listed, lambda_version = lambda_version,
       growth = growth_band |> select(m, status, from, to, Lmid, band, lambda, wt)),
  "data_inter/ukr_registration_completion.rds"
)
message("Done. data_inter/ukr_registration_completion.rds written for 09, 11, 13d and 13l.")
