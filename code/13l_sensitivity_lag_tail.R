# ==============================================================================
# STEP 13l - The registration lag beyond four years: a modelled tail
# ==============================================================================
#
# WHY
# ---
# 08b brings every event month to the completeness events reach at four years,
# the longest lag measured on a year of event months, so the correction is
# complete at 48 months by construction. 13d shows what a longer horizon would
# do by carrying the rate of the 42-48 month band on, flat, to 60 and 72 months.
# Neither is a model of the tail: the first truncates it and the second holds a
# rate constant that has been falling since the event.
#
# This step fits a curve to the monthly growth rates at lags of 12 to 48 months
# and extrapolates it, the way a reserving actuary closes a development triangle
# with a tail factor. Two curves, the two in common use:
#   exponential      lambda(L) = exp(a - b L): growth dies out, and the register
#                    has a limit
#   inverse power    lambda(L) = exp(a) L^-c: growth slows but, for c <= 1, never
#                    sums to a limit, so it is reported to a horizon only
# The curves are fitted to the month-level rates 08b pools into bands, by
# least squares weighted as 08b weights them (records x months), on the rates
# themselves since a month's net growth can be zero or negative. Below 48 months
# the measured band rates are kept; the curve replaces the truncation beyond.
#
# The band at 48-56 months is measured on three event months of one pair of
# releases and is left out of the fit, as it is out of the correction (08b). For
# the missing, who carry nearly all of the tail, it lies ABOVE both curves, so
# the curves are not an upper bound: what they show is how much the four-year
# horizon leaves out if growth keeps slowing as it has.
#
# Each horizon goes through production's imputation at the central values of
# the evidence on the missing alive and through the deterministic projection,
# every other input at its mode, as in 13d. The error of the tail is measured
# as 08b measures the factors': by resampling event months.
#
# INPUTS   data_inter/ukr_registration_completion.rds (08b)
#          data_inter/ukr_ualosses_resolution_model.rds, ukr_alive_missing.rds (09)
#          data_inter/ukr_registration_lag.rds (13d), the static inputs at the mode
# OUTPUTS  data_inter/ukr_registration_lag_tail.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

completion <- read_rds("data_inter/ukr_registration_completion.rds")
lambda_band <- completion$lambda_band
L_ref <- completion$L_ref
growth <- completion$growth
fit_from <- 12     # the first lag the curve is fitted to: past the first year's catch-up

# 1. THE CURVES ================================================================
curve_x <- list(exponential = identity, `inverse power` = log)
fit_curve <- function(g, form) {
  g <- g |> filter(Lmid > fit_from, Lmid <= L_ref)
  x <- curve_x[[form]](g$Lmid)
  # started from the log-linear fit to the positive rates
  pos <- g$lambda > 0
  start <- coef(lm(log(g$lambda[pos]) ~ x[pos], weights = g$wt[pos]))
  sse <- function(p) sum(g$wt * (g$lambda - exp(p[1] + p[2] * x))^2) / sum(g$wt)
  o <- optim(start, sse, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  c(a = o$par[[1]], slope = o$par[[2]])
}
# growth of the register from lag L0 to lag H under a fitted curve: the integral
# of lambda. H may be Inf for a curve with a limit.
tail_growth <- function(p, form, L0, H) {
  a <- p[["a"]]; s <- p[["slope"]]
  if (form == "exponential") {
    if (abs(s) < 1e-9) return(exp(a) * (H - L0))
    return(exp(a) * (exp(s * H) - exp(s * L0)) / s)
  }
  if (abs(s + 1) < 1e-9) return(exp(a) * (log(H) - log(L0)))
  exp(a) * (H^(s + 1) - L0^(s + 1)) / (s + 1)
}
has_limit <- function(p, form) if (form == "exponential") p[["slope"]] < 0 else p[["slope"]] < -1

# the measured band rates to L_ref, as 08b chains them
chain_factor <- function(L_from, s, lb) {
  lb <- lb |> filter(status == s, hi <= L_ref)
  exp(sum(pmax(0, pmin(lb$hi, L_ref) - pmax(lb$lo, L_from)) * lb$lambda))
}
# each event month's completion factor: measured to L_ref, the curve from there
# (or from the month's own lag, if it is past L_ref already) to the horizon
factors_under <- function(form, horizon, lb = lambda_band, fits = NULL) {
  completion$month |>
    mutate(f = map2_dbl(lag, status, \(l, s) chain_factor(l, s, lb)),
           f = if (is.na(form) || horizon <= L_ref) f else
             f * exp(map2_dbl(lag, status, \(l, s) {
               if (l >= horizon) 0 else tail_growth(fits[[s]][[form]], form, max(l, L_ref), horizon)
             })))
}

fits <- map(set_names(c("dead", "missing")), \(s)
            map(set_names(names(curve_x)), \(form) fit_curve(growth |> filter(status == s), form)))

cat("\n=== THE FITTED CURVES, MONTHLY GROWTH (%) AT THE BAND MIDPOINTS ===\n")
print(as.data.frame(
  lambda_band |>
    filter(lo >= fit_from) |>
    mutate(mid = (lo + hi) / 2,
           exponential = map2_dbl(status, mid, \(s, L) exp(fits[[s]]$exponential[["a"]] + fits[[s]]$exponential[["slope"]] * L)),
           `inverse power` = map2_dbl(status, mid, \(s, L) exp(fits[[s]]$`inverse power`[["a"]] + fits[[s]]$`inverse power`[["slope"]] * log(L))),
           fitted_to = hi <= L_ref) |>
    transmute(status, band, n_obs, fitted_to, measured = round(100 * lambda, 3),
              exponential = round(100 * exponential, 3), `inverse power` = round(100 * `inverse power`, 3))
))

scenarios <-
  bind_rows(
    tibble(scenario = "to 48 months (used)", form = NA_character_, horizon = L_ref),
    expand_grid(form = names(curve_x), horizon = c(72, 96, 120)) |>
      mutate(scenario = paste0(form, " tail, to ", horizon, " months")),
    tibble(form = names(curve_x), horizon = Inf) |>
      filter(map_lgl(form, \(fm) all(map_lgl(fits, \(f) has_limit(f[[fm]], fm))))) |>
      mutate(scenario = paste0(form, " tail, to its limit"))
  )

# 2. THROUGH THE IMPUTATION ====================================================
stocks_registered <-
  read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") |>
  summarise(n = sum(dx), .by = c(year, status)) |>
  mutate(status = as.character(status))
model <- read_rds("data_inter/ukr_ualosses_resolution_model.rds")
evidence <- read_rds("data_inter/ukr_alive_missing.rds")

# as 13d: each year's registered count spread over its months as v19 lists
# them, each month completed by its factor, and the missing imputed
military_with <- function(fac) {
  st <- fac |>
    mutate(share = n_v19 / sum(n_v19), .by = c(year, status)) |>
    left_join(stocks_registered |> rename(n_year = n), by = c("year", "status")) |>
    mutate(completed = n_year * share * f)
  dead <- st |> filter(status == "dead") |> summarise(confirmed = sum(completed), .by = year)
  miss <- st |> filter(status == "missing") |> select(month = m, year, missing_stock = completed)
  impute_missing(model, miss, evidence$captives_central, evidence$other_central) |>
    left_join(dead, by = "year") |>
    transmute(year, confirmed, missing = missing_stock, imputed_dead, total = confirmed + imputed_dead)
}
military <-
  scenarios |>
  mutate(res = map2(form, horizon, \(fm, h) military_with(factors_under(fm, h, fits = fits)))) |>
  unnest(res)

# the horizon used must reproduce 13d, and through it production
lag13d <- read_rds("data_inter/ukr_registration_lag.rds")
stopifnot(isTRUE(all.equal(
  military |> filter(scenario == "to 48 months (used)") |> arrange(year) |> pull(total),
  lag13d$military |> filter(scenario == "to 48 months (used)") |> arrange(year) |> pull(total)
)))

# 3. THE DETERMINISTIC PROJECTION ==============================================
mi <- mode_projection_inputs()
loss_with <- function(mil) {
  d <- mi$draws |>
    select(-draw_cmb, -conf_cmb) |>
    left_join(mil |> select(year, draw_cmb = total, conf_cmb = confirmed), by = "year")
  loss_at_mode(d, mi$static_inputs, mi$pop22_ini)
}
loss <-
  military |>
  nest(.by = c(scenario, form, horizon)) |>
  mutate(res = map(data, loss_with)) |>
  select(-data) |>
  unnest(res)
stopifnot(isTRUE(all.equal(
  loss |> filter(scenario == "to 48 months (used)") |> arrange(year, sex) |> pull(loss),
  lag13d$loss |> filter(scenario == "to 48 months (used)") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

# 4. THE ERROR OF THE TAIL =====================================================
# Event months resampled with replacement, each keeping all its release pairs,
# dead and missing together, as in 08b: the band rates and the curves are
# re-estimated on each resample and the military total recomputed.
set.seed(20260920)
B_tail <- 500
months_rates <- unique(growth$m)
replicate_total <- function(ms) {
  g <- tibble(m = ms) |> count(m, name = "times") |>
    inner_join(growth, by = "m", relationship = "one-to-many") |>
    mutate(wt = wt * times)
  lb <- lambda_band |>
    select(status, band, lo, hi, lambda0 = lambda) |>
    left_join(g |> summarise(lambda = sum(lambda * wt) / sum(wt), .by = c(status, band)),
              by = c("status", "band")) |>
    mutate(lambda = coalesce(lambda, lambda0))
  fb <- map(set_names(c("dead", "missing")), \(s)
            map(set_names(names(curve_x)), \(form) fit_curve(g |> filter(status == s), form)))
  scenarios |>
    filter(is.finite(horizon)) |>
    mutate(total = map2_dbl(form, horizon, \(fm, h) sum(military_with(factors_under(fm, h, lb, fb))$total)))
}
boot <- map_dfr(seq_len(B_tail), \(b) replicate_total(sample(months_rates, replace = TRUE)) |> mutate(b = b))
boot_summary <-
  boot |>
  summarise(lo = quantile(total, 0.025), mid = median(total), hi = quantile(total, 0.975), .by = scenario)

write_rds(list(fits = fits, fit_from = fit_from, scenarios = scenarios, military = military,
               loss = loss, boot = boot, boot_summary = boot_summary),
          "data_inter/ukr_registration_lag_tail.rds")

cat("\n=== MILITARY TOTAL BY TAIL, AT THE CENTRAL VALUES, WITH THE RESAMPLED 95% INTERVAL ===\n")
print(as.data.frame(
  military |>
    summarise(across(c(confirmed, missing, imputed_dead, total), \(x) round(sum(x))), .by = scenario) |>
    left_join(boot_summary |> mutate(across(c(lo, mid, hi), round)), by = "scenario")
))
cat("\n=== LOSS, MEN ===\n")
print(as.data.frame(loss |> filter(sex == "m") |> mutate(loss = round(loss, 2)) |>
                      select(scenario, year, loss) |> pivot_wider(names_from = year, values_from = loss)))
message("Done. data_inter/ukr_registration_lag_tail.rds written.")
