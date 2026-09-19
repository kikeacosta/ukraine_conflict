options(scipen = 9999)
try(dev.off(), silent = T)

# Under Rscript, printing a plot opens the default graphics device, which
# writes a stray "Rplots.pdf" into the working directory. Every figure here is
# saved explicitly with ggsave(), so that file is never wanted: send the
# default device to a null device when running non-interactively.
if (!interactive()) {
  options(device = function(...) grDevices::pdf(NULL))
}

# installing and loading required packages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# install pacman to streamline further package installation
if (!require("pacman", character.only = TRUE)) {
  install.packages("pacman", dep = TRUE)
  if (!require("pacman", character.only = TRUE)) {
    stop("Package pacman not found")
  }
}

library(pacman)
packages_CRAN <- c(
  "tidyverse",
  "countrycode",
  "lubridate",
  "readxl",
  "ungroup",
  "mgcv",
  "data.table",
  "MortalityLaws",
  "purrr",
  "demography",
  "patchwork",
  "fields",
  "forecast",
  "ggrepel",
  "ggh4x",
  "R.utils",
  "vital",
  "fst",
  "mc2d"
)

# Install required CRAN packages if not available yet
if (sum(!p_isinstalled(packages_CRAN)) > 0) {
  install_these <- packages_CRAN[!p_isinstalled(packages_CRAN)]
  for (i in 1:length(install_these)) {
    install.packages(install_these[i], dependencies = "Depends")
  }
}

# Load the required CRAN/github packages
p_load(packages_CRAN, character.only = TRUE)

copy_this <- function(x, row.names = FALSE, col.names = TRUE, ...) {
  write.table(
    x,
    file = paste0("clipboard-", object.size(x)),
    sep = "\t",
    row.names = row.names,
    col.names = col.names,
    ...
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# caching of the heavy raw sources ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Several raw files in data_input/ are far too large for version control
# (UCDP GED is 239 MB, the WPP fertility workbook 78 MB, the ualosses
# registers 27-46 MB). They are therefore NOT tracked in git - see
# data_input/README.md for the download locations and .gitignore for the
# exclusion rules.
#
# Every script that touches one of those files wraps the expensive part in
# cache_rds(). The cached extract is small, is tracked in git, and is read
# back on subsequent runs, so a fresh clone can run the whole pipeline end to
# end WITHOUT the raw downloads. `expr` is a lazily evaluated promise: it is
# only forced when the cache is missing.
#
#   dt <- cache_rds("data_inter/x.rds", { ...heavy processing... })
#
# Pass refresh = TRUE (or delete the .rds) to rebuild from the raw file.
#
# FORMAT: plain gzipped .rds, like every other intermediate in this pipeline.
# An earlier version used parquet, on the assumption that its ability to
# split a table across several files would be needed to keep everything under
# GitHub's 100 MB limit. It is not - the only file that large is the
# simulation output, which is derived, reproducible from a fixed seed, and
# therefore excluded from version control rather than committed. Without that
# need parquet earned nothing here: measured on the simulation draws it was
# LARGER than gzipped .rds (104 vs 98 MB) and slower to read, while adding
# the heavy `arrow` dependency to every clone.
cache_rds <- function(path, expr, refresh = FALSE) {
  if (!refresh && file.exists(path)) {
    message("cache hit  : ", path)
    return(as_tibble(readRDS(path)))
  }
  message("cache miss : rebuilding ", path, " from data_input/ ...")
  out <- expr
  saveRDS(out, path, compress = "gzip")
  message("cache built: ", path, " (", round(file.size(path) / 1e6, 1), " MB)")
  as_tibble(out)
}

# Fail early, and with a useful message, when a raw source is absent because
# the user cloned the repository without downloading the large inputs.
require_raw <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Raw input not found:\n  ",
      path,
      "\nThis file is too large for git. Either download it (see ",
      "data_input/README.md)\nor keep the cached .rds extract in ",
      "data_inter/ so this step can be skipped.",
      call. = FALSE
    )
  }
  path
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# life table ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Single life table from a data frame with columns age, mx, sex.
# Single years of age with an open interval at the last age.
lifetable <- function(dt_in) {
  x <- dt_in$age
  mx <- dt_in$mx
  sex <- unique(dt_in$sex)

  m <- length(x)
  n <- c(diff(x), NA)
  ax <- rep(0, m)
  if (x[1] != 0 | x[2] != 1) {
    ax <- n / 2
    ax[m] <- 1 / mx[m]
  } else {
    if (sex == "f") {
      if (mx[1] < 0.01724) {
        ax[1] <- 0.14903 - 2.05527 * mx[1]
      } else if (mx[1] >= 0.01724 & mx[1] < 0.06891) {
        ax[1] <- 0.04667 + 3.88089 * mx[1]
      } else {
        ax[1] <- 0.31411
      }
    }
    if (sex == "m") {
      if (mx[1] < 0.02300) {
        ax[1] <- 0.14929 - 1.99545 * mx[1]
      } else if (mx[1] >= 0.02300 & mx[1] < 0.08307) {
        ax[1] <- 0.02832 + 3.26021 * mx[1]
      } else {
        ax[1] <- 0.29915
      }
    }
    ax[-1] <- n[-1] / 2
    ax[m] <- 1 / mx[m]
  }
  qx <- n * mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px <- 1 - qx
  lx <- cumprod(c(1, px)) * 100000
  dx <- -diff(lx)
  Lx <- n * lx[-1] + ax * dx
  lx <- lx[-(m + 1)]
  Lx[m] <- lx[m] / mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx
  dt_out <- data.frame(age = x, mx, ex)
  return(dt_out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# vectorised life table ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Same arithmetic as lifetable() above, but computed column-wise over many
# groups at once, so the n_sim x 4 years x 2 sexes life tables of a full run
# take seconds instead of minutes. Used by 12 and 14.
#
# MX : n_age x G matrix of single-year death rates, open interval at the last
#      age. sx : length-G vector of sex ("f"/"m"), selecting the a0 rule.
# Returns lx, Lx and ex as n_age x G matrices.
lt_cols <- function(MX, sx) {
  m <- nrow(MX)
  G <- ncol(MX)
  ax <- matrix(0.5, m, G)
  m0 <- MX[1, ]
  ax[1, ] <- ifelse(
    sx == "f",
    ifelse(
      m0 < 0.01724,
      0.14903 - 2.05527 * m0,
      ifelse(m0 < 0.06891, 0.04667 + 3.88089 * m0, 0.31411)
    ),
    ifelse(
      m0 < 0.02300,
      0.14929 - 1.99545 * m0,
      ifelse(m0 < 0.08307, 0.02832 + 3.26021 * m0, 0.29915)
    )
  )
  ax[m, ] <- 1 / MX[m, ]

  qx <- MX / (1 + (1 - ax) * MX)
  qx[m, ] <- 1 # everyone dies in the open interval
  lx <- rbind(1, matrix(apply(1 - qx, 2, cumprod), nrow = m)) * 1e5
  dx <- -diff(lx)
  Lx <- lx[-1, , drop = FALSE] + ax * dx # n = 1 in every closed interval
  lx <- lx[-(m + 1), , drop = FALSE]
  Lx[m, ] <- lx[m, ] / MX[m, ] # open interval
  Lx[!is.finite(Lx)] <- 0
  Tx <- matrix(apply(Lx, 2, function(z) rev(cumsum(rev(z)))), nrow = m)
  list(lx = lx, Lx = Lx, ex = Tx / lx)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Arriaga decomposition ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-specific contribution of the mortality difference in each age interval
# to the total gap in life expectancy at birth, as the sum of
#   DE - person-years lost within the interval itself, and
#   IE - person-years lost later, because fewer people survive the interval.
#
# lb / lw are lt_cols() outputs for the baseline and the "war" schedule.
# BOTH terms are written in the same direction (baseline minus war), so both
# are positive when the war schedule costs life expectancy and sum(TE)
# reproduces e0_baseline - e0_war exactly. Writing one term in each direction
# makes the components cancel and the total come out negative.
arriaga_TE <- function(lb, lw) {
  m <- nrow(lb$lx)
  l0 <- lb$lx[1, ]
  nxt <- function(M) rbind(M[-1, , drop = FALSE], NA_real_)

  DE <- sweep(lb$lx, 2, l0, "/") * (lb$Lx / lb$lx - lw$Lx / lw$lx)
  IE <- (sweep(nxt(lb$lx), 2, l0, "/") -
    sweep(lb$lx, 2, l0, "/") * (nxt(lw$lx) / lw$lx)) * nxt(lw$ex)

  # open interval: nobody survives past it, so there is no indirect effect
  DE[m, ] <- (lb$lx[m, ] / l0) * (lb$ex[m, ] - lw$ex[m, ])
  IE[m, ] <- 0
  DE + IE
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# missing-combatant imputation ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UALosses nationality. Some releases hold the place of origin in the
# Nationality field (about 1,170 records in v18, 16 in v19) ("Kyiv, None",
# "Zaporizhzhja, Zaporizka urban community", ...), with From = "Unknown": the
# field is shifted, and the places are Ukrainian. Filtering on
# Nationality == "Ukraine" silently dropped them. A value that is not a
# country name is therefore read as Ukrainian, as is one the country matcher
# reads as Ukraine (the village of Ukrayinka); genuine foreign nationals keep
# their country and are still excluded.
ual_is_ukrainian <- function(nationality) {
  iso <- countrycode(nationality, "country.name", "iso3c", warn = FALSE)
  !is.na(nationality) & (nationality == "Ukraine" | is.na(iso) | iso %in% "UKR")
}

# UALosses dates are Excel serial numbers. Excel's day 1 is 1 January 1900
# but its calendar includes a 29 February 1900 that never existed, so for
# every date after February 1900 the serial counts from 30 December 1899.
# Counting from 1 January 1900 instead put every date two days late - the
# busiest day of disappearance came out as 26 February 2022 rather than the
# 24th - which moved events of 30 and 31 December into the following year.
excel_date <- function(x) {
  as.Date(suppressWarnings(as.numeric(x)), origin = "1899-12-30")
}

# Synthetic-cohort Markov imputation of the long-term missing. Defined once
# so that 09 (the central value of alpha), 13b (the sweep over alpha) and the
# range computed by alpha_evidence() below all use the SAME formula.
#
# tasas_long     observed one-window resolution rates by event-year cohort
#                and destination status: transitions |> filter(from ==
#                "missing") |> mutate(prop = n / sum(n), .by = year), as
#                built in 09.
# stock_missing  tibble(year, missing_stock) - the missing stock by event
#                year, from step 08's output.
# alpha          share of the NEVER-RESOLVED long-term missing who are
#                alive; 1 - alpha are dead. The one free parameter: its range
#                comes from alpha_evidence(), 11 draws it, 13b sweeps it.
#
# A "step" in the chain is the window between the two register releases,
# twelve months from v14 to v19; event-year cohorts, one year apart, stand in
# for duration since disappearance. Returns stock_missing with imputed_dead,
# imputed_alive and imputed_prisoner added. Every output is linear in alpha.
impute_missing <- function(alpha, tasas_long, stock_missing) {
  suelo_vivo <- alpha
  suelo_muerto <- 1 - alpha

  extraer_tasa <- function(target_year, target_status) {
    valor <- tasas_long$prop[
      tasas_long$year == target_year & tasas_long$status2 == target_status
    ]
    if (length(valor) == 0) return(0) else return(valor[1])
  }

  p_a22 <- extraer_tasa(2022, "alive")
  p_d22 <- extraer_tasa(2022, "dead")
  p_p22 <- extraer_tasa(2022, "prisoner")
  p_m22 <- extraer_tasa(2022, "missing")

  p_a23 <- extraer_tasa(2023, "alive")
  p_d23 <- extraer_tasa(2023, "dead")
  p_p23 <- extraer_tasa(2023, "prisoner")
  p_m23 <- extraer_tasa(2023, "missing")

  p_a24 <- extraer_tasa(2024, "alive")
  p_d24 <- extraer_tasa(2024, "dead")
  p_p24 <- extraer_tasa(2024, "prisoner")
  p_m24 <- extraer_tasa(2024, "missing")

  stock_missing %>%
    mutate(
      imputed_dead = case_when(
        year == 2022 ~ missing_stock * suelo_muerto,

        year == 2023 ~ (missing_stock * p_d22) +
          (missing_stock * p_m22 * suelo_muerto),

        year == 2024 ~ (missing_stock * p_d23) +
          (missing_stock * p_m23 * p_d22) +
          (missing_stock * p_m23 * p_m22 * suelo_muerto),

        year == 2025 ~ (missing_stock * p_d24) +
          (missing_stock * p_m24 * p_d23) +
          (missing_stock * p_m24 * p_m23 * p_d22) +
          (missing_stock * p_m24 * p_m23 * p_m22 * suelo_muerto)
      ),

      imputed_alive = case_when(
        year == 2022 ~ missing_stock * suelo_vivo,

        year == 2023 ~ (missing_stock * p_a22) +
          (missing_stock * p_m22 * suelo_vivo),

        year == 2024 ~ (missing_stock * p_a23) +
          (missing_stock * p_m23 * p_a22) +
          (missing_stock * p_m23 * p_m22 * suelo_vivo),

        year == 2025 ~ (missing_stock * p_a24) +
          (missing_stock * p_m24 * p_a23) +
          (missing_stock * p_m24 * p_m23 * p_a22) +
          (missing_stock * p_m24 * p_m23 * p_m22 * suelo_vivo)
      ),

      imputed_prisoner = case_when(
        year == 2022 ~ 0,

        year == 2023 ~ (missing_stock * p_p22),

        year == 2024 ~ (missing_stock * p_p23) +
          (missing_stock * p_m23 * p_p22),

        year == 2025 ~ (missing_stock * p_p24) +
          (missing_stock * p_m24 * p_p23) +
          (missing_stock * p_m24 * p_m23 * p_p22)
      )
    )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the range of alpha allowed by the evidence ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# alpha is the share of the never-resolved missing who are alive. Prisoners
# are alive whether still held or released, so the evidence on alpha is about
# prisoners of war, compared across two sources:
#
#   the register  v19 records prisoners and, as released_prisoner, returns
#                 from captivity; together they are the prisoners it knows of
#   official      about 7,000 Ukrainian prisoners of war held, plus 9,606
#   figures       people returned through exchanges: everyone ever taken
#                 prisoner, whether held or released since
#
# Prisoners the register does not record as such are either among its missing
# or not in the register at all.
#
#   min   0 - every unrecorded prisoner is outside the register (it never
#         listed 3,857 of the people it now records as returned)
#   mode  (pow_held + returned_from_captivity - register_alive) / residual -
#         every unrecorded prisoner is among the never-resolved missing
#   max   the register's own share of resolutions that are alive (to prisoner
#         or released), assuming the never-resolved are no more often alive
#         than those resolved: a living prisoner is listed or exchanged more
#         readily than a body is recovered
#
# The mode leans high: the returned figure includes civilians, whom the
# register does not list, and prisoners held in February and released by June
# appear in both figures. Captures after February are missing from it.
#
# residual is the number of missing never resolved at the end of the chain,
# the quantity alpha applies to: impute_missing() is linear in alpha, so it is
# the imputed dead at alpha = 0 minus those at alpha = 1.
#
# SOURCES
#   pow_held                 "about 7,000 Ukrainian prisoners of war" held by
#                            Russia - President Zelensky, Munich, 14 Feb 2026
#                            (Ukrinform, 14 Feb 2026).
#   returned_from_captivity  9,606 military personnel and civilians returned
#                            through exchanges by late June 2026 (Euromaidan
#                            Press, 24 Aug 2026).
pow_held <- 7000
returned_from_captivity <- 9606

# register_alive: the register's prisoners plus released prisoners, events
# 2022-2025, from 08's redistributed counts
alpha_evidence <- function(tasas_long, stock_missing, register_alive) {
  imp0 <- impute_missing(0, tasas_long, stock_missing)
  imp1 <- impute_missing(1, tasas_long, stock_missing)
  residual <- sum(imp0$imputed_dead - imp1$imputed_dead)

  resolved <- tasas_long |> filter(status2 != "missing")
  resolved_alive <- sum(resolved$n[resolved$status2 %in% c("alive", "prisoner")])
  resolved_all <- sum(resolved$n)
  unrecorded <- pow_held + returned_from_captivity - register_alive

  out <- tibble(
    alpha_min = 0,
    alpha_mode = unrecorded / residual,
    alpha_max = resolved_alive / resolved_all,
    residual = residual,
    unrecorded_prisoners = unrecorded,
    register_alive = register_alive,
    resolved_alive = resolved_alive,
    resolved_dead = resolved_all - resolved_alive,
    pow_held = pow_held,
    returned_from_captivity = returned_from_captivity
  )
  stopifnot(
    out$alpha_mode > out$alpha_min,
    out$alpha_mode < out$alpha_max,
    out$alpha_max < 1
  )
  out
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# demographic assumptions shared by the projection scripts ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sex ratio at birth (male births per female birth). Used by 11 and 13 to
# split the projected births between sexes. Ukraine has run at about 1.06
# for decades; the pipeline previously assumed 1.00 (an even split).
srb <- 1.06
prop_male_birth <- srb / (1 + srb)
prop_female_birth <- 1 / (1 + srb)

# Net emigration inputs and PERT bounds are built from data in 06 and written
# to data_inter/ (ukr_migration_bounds.rds); see documents/migration_methodology.md.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulation size ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THE number of Monte Carlo draws for the whole pipeline. Change it here and
# everything follows: 11 draws this many and names its cache after it, 15
# reads that cache and stamps the number into the tables, and the figures
# carry it in their captions.
#
#   1000   CONSTRUCTION (current). ~1.5 min in 11. Monte Carlo error of the
#          reported quantiles is about 0.08 years in the worst cell, the 2.5%
#          quantile of the male loss in 2025, so individual bounds move in the
#          second decimal between runs. Medians stay good to roughly 0.05.
#          Because the seed is fixed, two runs of the SAME size draw the same
#          PERT quantiles, so differences between pipeline states are far more
#          precise than either state on its own - which is what makes this
#          size usable for comparing changes.
#
#   20000  PRODUCTION, for the estimates that get reported. ~28 min in 11.
#          Brings that worst cell to about 0.018 years, inside the two decimal
#          places the intervals are printed to.
#
# SET THIS TO 20000 BEFORE GENERATING THE FINAL ESTIMATES. Nothing silently
# depends on remembering: the cache file name carries the size, 15 refuses a
# cache that does not match, and tables/_run_provenance.csv records what every
# table was built from.
n_sim <- 1000

# One-off override that does not need the file edited, e.g. a quick check at
# another size: UKR_N_SIM=20000 Rscript code/11_...R
if (nzchar(Sys.getenv("UKR_N_SIM"))) {
  n_sim <- as.integer(Sys.getenv("UKR_N_SIM"))
}
stopifnot(!is.na(n_sim), n_sim >= 100)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cohort-component projection for a single simulation draw ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Projects the population from 1 January 2022 to 31 December 2025, one year
# at a time, and returns the death components for each year-sex-age cell.
#
# Shared by 11 (probabilistic estimation) and 13 (migration sensitivity) so
# the two cannot drift apart.
#
# ARGUMENTS
#   sim_id         identifier carried through to the output
#   draws_this_sim one row per year: draw_cvs, draw_cmb, draw_mig, conf_cmb
#                  (conf_cmb = the register's confirmed dead that year)
#   static_inputs  year-sex-age grid with ems, w, mx, fx, prop_cvs,
#                  prop_cmb_dead, prop_cmb_miss
#   pop22_ini      population by sex and completed age on 1 January 2022
#
# AGE CONVENTION: row x in year t is ONE birth cohort, the one born in year
# t - x, i.e. completed age x on 31 December of year t. Newborns are row 0.
# This is the convention the migration input from 06 is built in. The SSSU
# population on 1 January 2022 at age x is the cohort born in 2021 - x, so it
# enters 2022 in row x + 1, exactly as if it had been aged on from 31 December
# 2021; nothing is dropped. (An earlier version joined it at row x and then
# overwrote row 0 with 2022 births, which lost the ~271k children born in
# 2021 and left rows from 2023 onward in two different conventions.)
#
# Period rates and age-at-death profiles (mx, fx, OHCHR, ualosses) are by
# completed age at the event, while row x spends year t aged x - 1 and then x.
# They are applied to row x unchanged, which is a half-year age offset applied
# uniformly to every row; see documents/migration_methodology.md, 11.12.
#
# TIMING: mid-year convention. Exposure is the average of the cohort's stock at
# the start and end of the year. For existing cohorts that means half the
# emigration and half the conflict deaths are removed before expected deaths
# are computed, and half the expected deaths after. Newborns start the year at
# zero and arrive through it, so their exposure is half of what survives.
run_single_sim <- function(sim_id, draws_this_sim, static_inputs, pop22_ini) {
  # carries a stock at 31 December into next year's rows: everyone moves up one
  # row, 100+ stays open, and row 0 is left empty for next year's births
  age_on <- function(d, next_year) {
    d %>%
      mutate(age = pmin(age + 1, 100)) %>%
      summarise(pop = sum(pop), .by = c(sex, age)) %>%
      bind_rows(tibble(sex = c("f", "m"), age = 0, pop = 0)) %>%
      mutate(year = next_year) %>%
      arrange(sex, age)
  }

  # ems is net emigration at the mode (POSITIVE = leaving, negative = return)
  # and w is the departures profile, summing to 1 within a year. The draw
  # changes only the year's total, and the difference from the mode is ADDED
  # along the departures profile, so the return cells stay as the registers
  # show them. Scaling the signed profile proportionally instead would inflate
  # returns whenever a draw raises emigration, and would blow up in years
  # whose net total is small next to its departures and returns. 13's
  # no-migration scenario passes ems = 0 and w = 0.
  mig_base_by_year <- static_inputs %>%
    group_by(year) %>%
    summarise(mig_base = sum(ems), .groups = "drop")

  df_sim <- static_inputs %>%
    left_join(draws_this_sim, by = "year") %>%
    left_join(mig_base_by_year, by = "year") %>%
    mutate(
      ems = ems + if_else(is.na(draw_mig), 0, draw_mig - mig_base) * w,

      # spread the drawn yearly totals over age and sex
      cvs = draw_cvs * prop_cvs,

      # partition combatants: everything up to the individually confirmed
      # count is "confirmed", the excess is attributed to the imputed missing
      ratio_confirmed = if_else(draw_cmb > 0, conf_cmb / draw_cmb, 1),
      ratio_confirmed = pmin(pmax(ratio_confirmed, 0), 1),

      # the two parts carry DIFFERENT age profiles. The register's missing are
      # older than its confirmed dead, so spreading both over the confirmed
      # profile placed imputed deaths at ages they did not occur - and 14's
      # decomposition is age-weighted, so that biases the loss they account
      # for. Registered deaths take the dead profile, imputed deaths take the
      # profile of the missing they come from.
      cmb_confirmed = draw_cmb * ratio_confirmed * prop_cmb_dead,
      cmb_imputed = draw_cmb * (1 - ratio_confirmed) * prop_cmb_miss,
      cmb = cmb_confirmed + cmb_imputed,

      cnf = cvs + cmb
    )

  results_list <- list()
  # 1 January 2022 by completed age is 31 December 2021 by completed age
  current_pop_ini <- age_on(pop22_ini %>% select(sex, age, pop), 2022)

  for (yr in 2022:2025) {
    yr_data <- df_sim %>%
      filter(year == yr) %>%
      left_join(current_pop_ini, by = c("year", "sex", "age")) %>%
      replace_na(list(pop = 0))

    yr_processed <- yr_data %>%
      mutate(
        pop2 = pop - 0.5 * ems - 0.5 * cnf, # population at risk
        # expected (non-conflict) deaths. mx is a rate per PERSON-YEAR (04
        # fits it to deaths over the mean of consecutive 1 January stocks),
        # so deaths must equal mx x exposure. exposure = pop2 - noc / 2, which
        # solves to noc = pop2 * mx / (1 + mx / 2). Using pop2 * mx instead
        # would apply mx / (1 - mx / 2): 6% too high at 80, 21% at 100.
        noc = pop2 * mx / (1 + 0.5 * mx),
        exposure = pop2 - 0.5 * noc, # person-years lived

        # row 0 is the cohort born during this year: it starts empty and
        # receives births = sum(ASFR x female exposure), split by the sex
        # ratio at birth. Births never touch the other rows.
        births = sum(fx * exposure),
        pop = case_when(
          age == 0 & sex == "f" ~ prop_female_birth * births,
          age == 0 & sex == "m" ~ prop_male_birth * births,
          .default = pop
        ),
        # same bookkeeping for the newborns, but they start the year at zero:
        # the average of the start (0) and end stock is half of what survives
        pop2 = ifelse(age == 0, 0.5 * (pop - ems - cnf), pop2),
        noc = ifelse(age == 0, pop2 * mx / (1 + 0.5 * mx), noc),
        exposure = ifelse(age == 0, pop2 - 0.5 * noc, exposure),

        dx = noc + cnf, # all-cause deaths
        pop_end = pop - ems - dx # survivors at 31 December
      )

    # keep the components only; every rate downstream is derived from these
    results_list[[as.character(yr)]] <- yr_processed %>%
      select(
        year, sex, age,
        pop = exposure,
        expected = noc,
        civilian = cvs,
        combatant_confirmed = cmb_confirmed,
        combatant_imputed = cmb_imputed
      )

    # age the survivors on into next year's starting population
    if (yr < 2025) {
      current_pop_ini <- age_on(yr_processed %>% select(sex, age, pop = pop_end), yr + 1)
    }
  }

  bind_rows(results_list) %>% mutate(sim_id = sim_id)
}
