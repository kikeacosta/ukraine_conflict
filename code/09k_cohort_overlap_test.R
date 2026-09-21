# ==============================================================================
# STEP 09k - Do later cohorts resolve as earlier ones did? The overlap test
# ==============================================================================
#
# WHY
# ---
# The resolution model takes the hazards of leaving "missing" to depend on the
# months since the event and on nothing else: one pattern for every event month
# (ual_fit_durations(), 00_setup.R). That is what lets a hazard measured at 36-56
# months on the events of 2022 be carried onto the missing of 2024 and 2025, who
# have not reached those durations. It cannot be tested where it matters most,
# because no release sees a 2024-25 cohort beyond three years. It can be tested
# where cohorts OVERLAP: within a six-month band of durations, in one and the
# same window between releases, the event months of two adjacent years are both
# at risk - 2024 and 2025 events at 9-18 months, 2023 and 2024 at 18-30, 2022 and
# 2023 at 30-42 - so the window's batch effects are the same for both and what
# differs is the cohort.
#
# THE TEST
# --------
# The model is fitted again with one multiplier per event year and outcome, the
# events of 2022 as the reference: hazard = band x window x cohort. The
# multipliers are identified by the overlaps alone, as a chain of adjacent
# contrasts (2023 against 2022, 2024 against 2023, 2025 against 2024), so they
# are reported as those ratios. Two things limit what they say, and both are
# stated with the results:
#   - within a band the hazard is taken as constant, so a hazard that falls or
#     rises with duration INSIDE a band is read as a cohort effect: the later
#     cohort sits at the band's shorter durations;
#   - the counts are overdispersed (a release records resolutions in batches),
#     so the test is a quasi-likelihood F test on the Pearson dispersion, and
#     the intervals use the covariance scaled by it, as production's do.
# This is not the multistate model with cohort effects that the project tested
# and rejected (09g): there the cohort effects were hazard ratios over the whole
# follow-up, identified by almost nothing, since the 2022 and 2025 cohorts share
# no time at risk. Here each ratio rests only on the bands two cohorts share.
#
# WHAT IT IS FOR
# --------------
# It is the reason the duration model is not the imputation. A model with one
# pattern for every cohort, carried to 48 months, takes the rate at which the
# missing leave the register beyond list maintenance - observed almost only on
# the events of 2022, at long durations, in the releases that cleaned the
# register up - and applies it to the missing of 2023-2025, who show no such
# excess where they can be seen. The last section gives the people that channel
# would put alive under the pooled hazards, under each cohort's own multiplier
# (which cannot be carried beyond the durations a cohort has reached, as the
# result shows), and with the channel kept to the events of 2022. The estimate
# uses none of them: it counts the missing as dead but for the prisoners and a
# share bounded by each cohort's own resolutions (military_draws(), 00_setup.R).
#
# INPUTS   data_inter/ualosses_window_transitions.rds (09's cached linkage),
#          data_inter/ukr_registration_completion.rds (08b),
#          data_inter/ukr_military_inputs.rds (09)
# OUTPUTS  data_inter/ukr_ualosses_cohort_overlap.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

feb22 <- as.Date("2022-02-01")
completion <- read_rds("data_inter/ukr_registration_completion.rds")
linkage <- read_rds("data_inter/ualosses_window_transitions.rds")

# 1. THE SAME COUNTS THE PRODUCTION MODEL IS FITTED TO ==========================
# list maintenance exactly as 09 and 09j apply it
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
cells <- ual_duration_cells(maintained) |> mutate(cohort = year(month))
cohorts <- sort(unique(cells$cohort))
breaks <- ual_duration_breaks
nb <- length(breaks) - 1
k <- length(ual_resolutions)
nw <- length(ual_window_months)
nc <- length(cohorts)

# 2. WHAT THE OVERLAPS HOLD =====================================================
# the observed monthly rate of each outcome by band, window and cohort, where a
# band and window hold two cohorts: events over the person-months at risk
band_of <- function(d) cut(d, breaks, right = FALSE,
                           labels = paste0(head(breaks, -1), "-", breaks[-1]))
overlap <-
  cells |>
  mutate(band = band_of(d0), at_risk = missing + dead + prisoner + no_longer_listed,
         months = at_risk * (d1 - d0)) |>
  summarise(at_risk = sum(at_risk), months = sum(months),
            across(all_of(ual_resolutions), sum), .by = c(band, from_release, cohort)) |>
  filter(n() > 1, .by = c(band, from_release)) |>
  mutate(across(all_of(ual_resolutions), \(x) 100 * x / months, .names = "rate_{.col}")) |>
  arrange(band, from_release, cohort)

cat("\n=== MONTHLY RATES (%) WHERE TWO COHORTS SHARE A BAND AND A WINDOW ===\n")
print(as.data.frame(overlap |> transmute(band, window = from_release, cohort, at_risk = round(at_risk),
                                         dead = round(rate_dead, 3), prisoner = round(rate_prisoner, 3),
                                         out_of_register = round(rate_no_longer_listed, 3))))

# 3. THE MODEL WITH COHORT MULTIPLIERS ==========================================
# theta: log-hazards by band and outcome, log-multipliers of the later windows,
# then log-multipliers of the later cohorts for the outcomes in `free`
window <- match(cells$from_release, ual_releases$release)
cohort <- match(cells$cohort, cohorts)
y <- as.matrix(cells[, c("missing", ual_resolutions)])
fit_with <- function(free) {
  nf <- length(free)
  n_base <- nb * k + (nw - 1) * k
  unpack <- function(theta) {
    hz <- ual_theta_hazards(theta[seq_len(n_base)], nb, k)
    gam <- matrix(1, nc, k, dimnames = list(NULL, ual_resolutions))
    if (nf) gam[-1, free] <- exp(theta[n_base + seq_len((nc - 1) * nf)])
    list(h = hz$h, mult = hz$mult, gam = gam)
  }
  cell_probs <- function(theta) {
    u <- unpack(theta)
    p <- matrix(0, nrow(cells), k + 1, dimnames = list(NULL, c("missing", ual_resolutions)))
    for (w in seq_len(nw)) for (cc in seq_len(nc)) {
      i <- window == w & cohort == cc
      if (any(i)) p[i, ] <- ual_resolve(sweep(u$h, 2, u$mult[w, ] * u$gam[cc, ], "*"),
                                        cells$d0[i], cells$d1[i], breaks)
    }
    p
  }
  nll <- function(theta) -sum(y * log(pmax(cell_probs(theta)[, colnames(y)], 1e-300)))
  n_par <- n_base + (nc - 1) * nf
  start <- c(rep(log(0.003), nb * k), rep(0, (nw - 1) * k), rep(0, (nc - 1) * nf))
  fit <- optim(start, nll, method = "L-BFGS-B",
               lower = c(rep(-20, nb * k), rep(-10, n_par - nb * k)),
               upper = c(rep(0, nb * k), rep(10, n_par - nb * k)),
               hessian = TRUE, control = list(maxit = 10000, factr = 1e5))
  stopifnot(fit$convergence == 0)
  rest <- fit$par[-seq_len(nb * k)]
  at_bound <- c(fit$par[seq_len(nb * k)] <= -12, abs(rest) >= 9.9)
  expected <- cell_probs(fit$par)[, colnames(y)] * rowSums(y)
  pearson <- sum(((y - expected)^2 / expected)[expected > 1])
  df <- nrow(y) * k - sum(!at_bound)
  vcov <- matrix(0, n_par, n_par)
  vcov[!at_bound, !at_bound] <- solve(fit$hessian[!at_bound, !at_bound])
  list(par = fit$par, nll = fit$value, n_free = sum(!at_bound), df = df,
       dispersion = max(1, pearson / df), vcov = vcov * max(1, pearson / df),
       at_bound = at_bound, unpack = unpack, free = free, n_base = n_base)
}

message("  fitting the pooled model and the models with cohort multipliers ...")
pooled <- fit_with(character())
by_outcome <- map(set_names(ual_resolutions), fit_with)
all_three <- fit_with(ual_resolutions)

# the pooled model is production's
production <- read_rds("data_inter/ukr_ualosses_resolution_model.rds")
stopifnot(isTRUE(all.equal(unname(pooled$par), unname(production$theta), tolerance = 1e-3)))

# quasi-likelihood F test of the cohort multipliers, on the larger model's dispersion
f_test <- function(big, small = pooled) {
  q <- big$n_free - small$n_free
  f <- (2 * (small$nll - big$nll) / q) / big$dispersion
  tibble(multipliers = q, deviance = 2 * (small$nll - big$nll), dispersion = big$dispersion,
         F = f, p = pf(f, q, big$df, lower.tail = FALSE))
}
tests <- bind_rows(
  imap_dfr(by_outcome, \(m, o) f_test(m) |> mutate(outcome = o, .before = 1)),
  f_test(all_three) |> mutate(outcome = "all three", .before = 1)
)

# each cohort against the one before it, with a 95% interval from the scaled covariance
adjacent <- imap_dfr(by_outcome, function(m, o) {
  idx <- m$n_base + seq_len(nc - 1)
  lg <- c(0, m$par[idx])                      # log-multipliers against 2022
  V <- rbind(0, cbind(0, m$vcov[idx, idx]))
  bound <- c(FALSE, m$at_bound[idx])
  map_dfr(2:nc, function(j) {
    d <- lg[j] - lg[j - 1]
    se <- sqrt(V[j, j] + V[j - 1, j - 1] - 2 * V[j, j - 1])
    tibble(outcome = o, cohort = cohorts[j], against = cohorts[j - 1], ratio = exp(d),
           lo = exp(d - 1.96 * se), hi = exp(d + 1.96 * se),
           at_bound = bound[j] | bound[j - 1])
  })
})

cat("\n=== DO THE HAZARDS DIFFER BY COHORT? QUASI-LIKELIHOOD F TESTS ===\n")
print(as.data.frame(tests |> mutate(across(c(deviance, dispersion, F), \(x) round(x, 2)), p = signif(p, 3))))
cat("\n=== EACH COHORT'S HAZARD AGAINST THE COHORT BEFORE IT, IN THE BANDS THEY SHARE ===\n")
print(as.data.frame(adjacent |> mutate(across(c(ratio, lo, hi), \(x) signif(x, 3)))))

# 4. WHAT IT WOULD DO TO THE PEOPLE PUT ALIVE OUT OF THE REGISTER ===============
# v19's missing by event month, completed for registration lag, as production
# carries them; each month projected to the horizon under the pooled hazards,
# under its own cohort's multipliers (all three outcomes), and with the channel
# out of the register kept to the events of 2022
mil <- read_rds("data_inter/ukr_military_inputs.rds")
stock <- mil$month |> filter(status == "missing") |>
  transmute(month, year, missing_stock = registered * mil$factor_point[row])
d_now <- pmin(ual_months_since(ual_releases$date[nrow(ual_releases)], stock$month), production$horizon)
leaves <- function(h_of_year) {
  map_dbl(seq_len(nrow(stock)), \(i)
          stock$missing_stock[i] *
            ual_resolve(h_of_year(stock$year[i]), d_now[i], production$horizon, breaks)[, "no_longer_listed"])
}
u_all <- all_three$unpack(all_three$par)
proj_mult <- colSums(u_all$mult * ual_window_months) / sum(ual_window_months)
h_2022_only <- production$h
h_2022_only[, "no_longer_listed"] <- 0
out_of_register <-
  stock |>
  mutate(pooled = leaves(\(yr) production$h),
         own_cohort = leaves(\(yr) sweep(u_all$h, 2, proj_mult * u_all$gam[match(yr, cohorts), ], "*")),
         events_2022_only = leaves(\(yr) if (yr == 2022) production$h else h_2022_only)) |>
  summarise(missing = sum(missing_stock), across(c(pooled, own_cohort, events_2022_only), sum), .by = year)
stopifnot(abs(sum(out_of_register$pooled) -
                sum(read_rds("data_inter/ukr_ualosses_projection_table.rds")$imputed_unlisted)) < 1)

cat("\n=== PEOPLE THE CHANNEL OUT OF THE REGISTER PUTS ALIVE, BY EVENT YEAR ===\n")
print(as.data.frame(out_of_register |> mutate(across(-year, round)) |>
                      bind_rows(out_of_register |> summarise(across(-year, \(x) round(sum(x)))) |> mutate(year = NA))))

write_rds(list(overlap = overlap, tests = tests, adjacent = adjacent, out_of_register = out_of_register,
               dispersion = c(pooled = pooled$dispersion, all_three = all_three$dispersion)),
          "data_inter/ukr_ualosses_cohort_overlap.rds")
message("Done. data_inter/ukr_ualosses_cohort_overlap.rds written.")
