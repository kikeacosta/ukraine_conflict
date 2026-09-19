# ==============================================================================
# STEP 08b - Registration lag: how complete is each event month in v19?
# ==============================================================================
#
# WHY
# ---
# The register keeps adding deaths and disappearances long after the event:
# for 2025 events, v14 (16 Sep 2025) listed 9,001 dead, v19 (19 Sep 2026)
# 16,457, and even 2022 events keep growing. The imputation in 09 only
# reclassifies people already listed, so a loss not yet registered would not
# be counted. This step measures how far each event month's v19 count falls
# short of the completeness older events have reached; 09 scales the v19 dead
# and missing up by that factor before the imputation.
#
# METHOD - a chain-ladder on the four releases
# --------------------------------------------
# For each event month m, each status s (dead, missing) and each consecutive
# pair of releases r -> r+1:
#   net additions = records listed with status s in r+1 whose key is absent
#                   from r in any status
#                   - records with status s in r whose key is absent from r+1
#                   in any status
# A key present in both releases is left out whatever its status, so a missing
# person who is found dead, taken prisoner or released is a resolution, not a
# registration or a de-registration. The monthly growth rate
#   lambda = log(1 + net additions / N_r) / (months between the releases)
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
# INPUT   the four register releases (gitignored), through a cache of
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
  rd <- function(p) {
    ual_read_release(p) |>
      mutate(m = floor_date(date_evnt, "month")) |>
      filter(ukrainian, !is.na(m), m >= as.Date("2022-02-01"), m <= as.Date("2025-12-01")) |>
      ual_one_per_key() |>
      select(key, m, status)
  }
  d <- set_names(map(releases$path, rd), releases$release)
  losses <- \(x) filter(x, status %in% c("dead", "missing"))

  # one long table: stocks by release, and additions and drops between
  # consecutive releases. A record counts as added or dropped only if its key
  # is absent from the other release in ANY status, so a missing person who
  # becomes a prisoner or is released is a resolution, not a de-registration.
  stock <- imap_dfr(d, \(x, r) count(losses(x), m, status, name = "n") |>
                      mutate(kind = "stock", from = r, to = NA_character_))
  flows <- map_dfr(1:3, function(i) {
    a <- d[[i]]; b <- d[[i + 1]]
    bind_rows(
      losses(b) |> anti_join(a, by = "key") |> count(m, status, name = "n") |> mutate(kind = "added"),
      losses(a) |> anti_join(b, by = "key") |> count(m, status, name = "n") |> mutate(kind = "dropped")
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
    lambda = log1p((added - dropped) / n_from) / dL
  ) |>
  filter(L0 >= 1, n_from >= 50)

bands <- c(1, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48, 56)
lambda_band <-
  growth |>
  mutate(band = cut(Lmid, bands, right = TRUE)) |>
  filter(!is.na(band)) |>
  summarise(lambda = weighted.mean(lambda, n_from * dL), n_obs = n(), .by = c(status, band)) |>
  mutate(lo = bands[as.integer(band)], hi = bands[as.integer(band) + 1]) |>
  arrange(status, lo)

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
stopifnot(all(completion_month$factor > 0))

completion_year <-
  completion_month |>
  summarise(n_completed = sum(n_v19 * factor), n_v19 = sum(n_v19), .by = c(year, status)) |>
  mutate(ratio = n_completed / n_v19) |>
  arrange(status, year)

cat("\n=== V19 COUNTS BROUGHT TO THE COMPLETENESS EVENTS REACH AT FOUR YEARS ===\n")
print(as.data.frame(completion_year |> mutate(across(c(n_v19, n_completed), round),
                                              ratio = round(ratio, 4))))

write_rds(
  list(lambda_band = lambda_band, L_ref = L_ref, month = completion_month, year = completion_year),
  "data_inter/ukr_registration_completion.rds"
)
message("Done. data_inter/ukr_registration_completion.rds written for 09 and 13d.")
