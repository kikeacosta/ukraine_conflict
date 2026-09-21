# ==============================================================================
# STEP 08b - Registration lag: how complete is each event month in v19?
# ==============================================================================
#
# WHY
# ---
# The register keeps adding deaths and disappearances long after the event:
# for 2025 events, v14 (16 Sep 2025) listed 9,001 dead, v19 (21 Jul 2026)
# 16,457, and even 2022 events keep growing. The imputation in 09 only
# reclassifies people already listed, so a loss not yet registered would not
# be counted. This step measures how far each event month's v19 count falls
# short of the completeness older events have reached; 09 scales the v19 dead
# and missing up by that factor before the imputation.
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
# factors below one - a register that loses deaths as it ages. Measured, the
# loose test removes 934, 345 and 494 records between v14, v16, v18 and v19 against 1, 0 and
# 0 for the strict one, so nearly all of it was matching namesakes. The guards
# below stop the asymmetry recurring.
# A key present in both releases is left out whatever its status, so a missing
# person who is found dead, taken prisoner or released is a resolution, not a
# registration or a de-registration. The dead leave the register only through
# list maintenance - corrections and merged duplicates - so their drops are
# netted out. The missing are different: a missing person who drops off the
# register is a resolution in 09's model, counted there as found alive, so
# netting the drops here as well would take them out twice. Their growth is
# therefore new persons only; a record re-keyed or listed again after a gap is
# not new. The monthly growth rate
#   lambda = log(1 + growth / N_r) / (months between the releases)
# is attached to the lag since the event (midpoint of the pair), and pooled
# by lag band, weighting by N_r x months.
#
# Chaining those rates from each month's lag at v19 up to 48 months - the
# longest lag measured on at least twelve event months - gives the factor that
# brings the month to the completeness events reach at four years. Beyond 48
# months the rates rest on five event months of a single pair of releases,
# where v19 deleted a batch of February 2022 missing records in one clean-up:
# a one-off, not a lag. The register still grows after four years, so the
# factor is a lower bound on the lag; 13d shows how far a longer horizon would
# move the results.
#
# INPUT   the six register releases (gitignored), through a cache of
#         aggregate counts: data_inter/ualosses_registration_by_month.rds
# OUTPUT  data_inter/ukr_registration_completion.rds: the growth rates by lag
#         band, the horizon, and v19's dead and missing by event month with
#         each month's completion factor
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

releases <- ual_releases

# 1. AGGREGATE COUNTS FROM THE RELEASES (cached; no names leave this block) ====
reg_counts <- cache_rds("data_inter/ualosses_registration_by_month.rds", {
  regs <- map(set_names(releases$path, releases$release), ual_read_release)
  rd <- function(r) {
    r |>
      mutate(m = floor_date(date_evnt, "month"),
             np = str_squish(str_remove(key, "\\s+\\S+$"))) |>
      filter(ukrainian, !is.na(m), m >= as.Date("2022-02-01"), m <= as.Date("2025-12-01")) |>
      ual_one_per_key() |>
      select(key, np, name, dob, m, status)
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
      transmute(key, np = str_squish(str_remove(key, "\\s+\\S+$")), name, dob, status,
                m = floor_date(date_evnt, "month")) |>
      distinct()
    gone <- earlier |> anti_join(regs[[i + 1]] |> distinct(key), by = "key")
    # The corrected-key test, applied to each status by what it is for.
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
    # The missing keep the linkage's own test: surname and first name, and
    # either the date of birth or the event month, against every departed
    # record. The dead use the patronymic as well and look only among departed
    # MISSING records, because in a register of this size a surname, first name
    # and event month coincide often: matched loosely against all departures
    # the test removes 934, 345 and 494 of the new dead between v14, v16, v18 and v19,
    # against 1, 0 and 0 once the patronymic is required - so the loose version
    # is matching namesakes, not the same person.
    drop_keys <- function(x, g, k) {
      x |>
        anti_join(g |> filter(!is.na(dob)) |> distinct(across(all_of(c(k, "dob")))),
                  by = c(k, "dob")) |>
        anti_join(g |> filter(!is.na(m)) |> distinct(across(all_of(c(k, "m")))),
                  by = c(k, "m"))
    }
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
    # same releases beside the 5,981 measured on them, so the factors remain a
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
  bind_rows(stock, flows)
})

stock <- reg_counts |> filter(kind == "stock") |> transmute(m, status, n, release = from)
flows <- reg_counts |> filter(kind != "stock") |> select(m, status, n, kind, from, to)

# 2. GROWTH RATES BY LAG ======================================================
lag_months <- function(release_date, m) as.numeric(release_date - (m + 14)) / 30.44
rel_date <- set_names(releases$date, releases$release)

growth <-
  flows |>
  pivot_wider(names_from = kind, values_from = n, values_fill = 0) |>
  left_join(stock |> rename(from = release, n_from = n), by = c("m", "status", "from")) |>
  mutate(
    L0 = lag_months(rel_date[from], m), L1 = lag_months(rel_date[to], m),
    dL = L1 - L0, Lmid = (L0 + L1) / 2,
    growth = if_else(status == "missing", added_new, added - added_recoded - dropped),
    lambda = log1p(growth / n_from) / dL
  ) |>
  # February 2022 holds five days of events and takes the register's batch
  # operations on the start of the war - v18 added 221 missing to it at once,
  # v19 deleted 335 - so it is left out of the rates. At its lag it is past
  # the horizon and needs no factor itself.
  filter(L0 >= 1, n_from >= 50, m != as.Date("2022-02-01"))

bands <- c(1, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48, 56)
lambda_band <-
  growth |>
  mutate(band = cut(Lmid, bands, right = TRUE)) |>
  filter(!is.na(band)) |>
  summarise(lambda = weighted.mean(lambda, n_from * dL), n_obs = n(), .by = c(status, band)) |>
  mutate(lo = bands[as.integer(band)], hi = bands[as.integer(band) + 1]) |>
  arrange(status, lo)

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
  filter(release == "v19") |>
  mutate(year = year(m), lag = lag_months(rel_date[["v19"]], m),
         factor = map2_dbl(lag, status, chain_factor)) |>
  select(m, year, status, lag, n_v19 = n, factor)
stopifnot(all(completion_month$factor >= 1 - 1e-9))   # a lag correction can only add

completion_year <-
  completion_month |>
  summarise(n_completed = sum(n_v19 * factor), n_v19 = sum(n_v19), .by = c(year, status)) |>
  mutate(ratio = n_completed / n_v19) |>
  arrange(status, year)

cat("\n=== V19 COUNTS BROUGHT TO THE COMPLETENESS EVENTS REACH AT FOUR YEARS ===\n")
print(as.data.frame(completion_year |> mutate(across(c(n_v19, n_completed), round),
                                              ratio = round(ratio, 4))))

# 4. SAMPLING ERROR OF THE FACTORS ============================================
# The band rates are averages over event months, and the months differ: a
# release can add a batch to some months and not others. The error of the
# factors is measured by resampling event months with replacement - each
# month keeps all its release pairs, dead and missing together - and
# re-estimating the band rates and the chained factors, B times. A band that a
# resample happens to leave empty keeps its estimate. 11 draws one replicate
# per simulation, so the late registrations and the completed missing carry
# this error into the results.
set.seed(20260919)
B_lag <- 2000
months_rates <- unique(growth$m)
growth_band <-
  growth |>
  mutate(band = cut(Lmid, bands, right = TRUE)) |>
  filter(!is.na(band)) |>
  mutate(wt = n_from * dL)
lag_replicate <- function(ms) {
  lb <-
    tibble(m = ms) |>
    count(m, name = "times") |>
    inner_join(growth_band, by = "m", relationship = "one-to-many") |>
    summarise(lambda = sum(lambda * wt * times) / sum(wt * times), .by = c(status, band))
  lb <-
    lambda_band |>
    select(status, band, lo, hi, lambda0 = lambda) |>
    left_join(lb, by = c("status", "band")) |>
    mutate(lambda = coalesce(lambda, lambda0))
  map2_dbl(completion_month$lag, completion_month$status, \(l, s) chain_factor(l, s, lb = lb))
}
lag_factors <- sapply(seq_len(B_lag), \(b) lag_replicate(sample(months_rates, replace = TRUE)))
stopifnot(all(lag_factors > 0), nrow(lag_factors) == nrow(completion_month))

cat("\n=== SAMPLING ERROR OF THE COMPLETED COUNTS (", B_lag, "resamples of event months) ===\n")
print(as.data.frame(
  completion_month |>
    select(year, status, n_v19) |>
    bind_cols(as_tibble(lag_factors, .name_repair = \(x) paste0("b", seq_along(x)))) |>
    summarise(across(starts_with("b"), \(f) sum(n_v19 * f)), n_v19 = sum(n_v19), .by = c(year, status)) |>
    pivot_longer(starts_with("b"), values_to = "completed") |>
    summarise(n_v19 = first(n_v19), lo = quantile(completed, 0.025), mid = median(completed),
              hi = quantile(completed, 0.975), .by = c(year, status)) |>
    mutate(across(c(lo, mid, hi), round)) |>
    arrange(status, year)
))

# growth: the rates the bands pool, one row per event month, status and pair of
# releases, for 13l to fit a curve through
write_rds(
  list(lambda_band = lambda_band, L_ref = L_ref, month = completion_month, year = completion_year,
       factor_draws = lag_factors,
       growth = growth_band |> select(m, status, from, to, Lmid, band, lambda, wt)),
  "data_inter/ukr_registration_completion.rds"
)
message("Done. data_inter/ukr_registration_completion.rds written for 09, 11, 13d and 13l.")
