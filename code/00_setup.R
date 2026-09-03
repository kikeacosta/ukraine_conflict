options(scipen = 9999)
try(dev.off(), silent = T)

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
  "arrow",
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
# cache_parquet(). The cached extract is small, is tracked in git, and is
# read back on subsequent runs, so a fresh clone can run the whole pipeline
# end to end WITHOUT the raw downloads. `expr` is a lazily evaluated promise:
# it is only forced when the cache is missing.
#
#   dt <- cache_parquet("data_inter/x.parquet", { ...heavy processing... })
#
# Pass refresh = TRUE (or delete the .parquet) to rebuild from the raw file.
cache_parquet <- function(path, expr, refresh = FALSE) {
  if (!refresh && file.exists(path)) {
    message("cache hit  : ", path)
    return(as_tibble(arrow::read_parquet(path)))
  }
  message("cache miss : rebuilding ", path, " from data_input/ ...")
  out <- expr
  arrow::write_parquet(out, path, compression = "zstd")
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
      "data_input/README.md)\nor keep the cached .parquet extract in ",
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
# groups at once, so the 40,000 life tables implied by 5,000 draws x 4 years
# x 2 sexes take seconds instead of minutes. Used by 12 and 14.
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
# demographic assumptions shared by the projection scripts ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sex ratio at birth (male births per female birth). Used by 11 and 13 to
# split the projected births between sexes. Ukraine has run at about 1.06
# for decades; the pipeline previously assumed 1.00 (an even split).
srb <- 1.06
prop_male_birth <- srb / (1 + srb)
prop_female_birth <- 1 / (1 + srb)

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
#   draws_this_sim one row per year: draw_cvs, draw_cmb, min_cmb
#   static_inputs  year-sex-age grid with ems, mx, fx, prop_cvs, prop_cmb
#   pop22_ini      population by sex and age on 1 January 2022
#
# TIMING: mid-year convention. Half the emigration and half the conflict
# deaths are removed before exposure is computed, half the expected deaths
# after, so exposure approximates person-years lived.
run_single_sim <- function(sim_id, draws_this_sim, static_inputs, pop22_ini) {
  df_sim <- static_inputs %>%
    left_join(draws_this_sim, by = "year") %>%
    mutate(
      # spread the drawn yearly totals over age and sex
      cvs = draw_cvs * prop_cvs,
      cmb = draw_cmb * prop_cmb,

      # partition combatants: everything up to the individually confirmed
      # count is "confirmed", the excess is attributed to the imputed missing
      ratio_confirmed = if_else(draw_cmb > 0, min_cmb / draw_cmb, 1),
      ratio_confirmed = pmin(pmax(ratio_confirmed, 0), 1),
      cmb_confirmed = cmb * ratio_confirmed,
      cmb_imputed = cmb * (1 - ratio_confirmed),

      cnf = cvs + cmb
    )

  results_list <- list()
  current_pop_ini <- pop22_ini %>% mutate(year = 2022)

  for (yr in 2022:2025) {
    yr_data <- df_sim %>%
      filter(year == yr) %>%
      left_join(current_pop_ini, by = c("year", "sex", "age")) %>%
      replace_na(list(pop = 0))

    yr_processed <- yr_data %>%
      mutate(
        pop2 = pop - 0.5 * ems - 0.5 * cnf, # population at risk
        noc = pop2 * mx, # expected (non-conflict) deaths
        exposure = pop2 - 0.5 * noc, # person-years lived

        # age 0 is not carried in from the previous year: it is born during
        # this one. Births = sum(ASFR x female exposure), split by the sex
        # ratio at birth.
        births = sum(fx * exposure),
        pop = case_when(
          age == 0 & sex == "f" ~ prop_female_birth * births,
          age == 0 & sex == "m" ~ prop_male_birth * births,
          .default = pop
        ),
        # then re-run the same bookkeeping for the newborn cohort
        pop2 = ifelse(age == 0, pop - 0.5 * ems - 0.5 * cnf, pop2),
        noc = ifelse(age == 0, pop2 * mx, noc),
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
      current_pop_ini <- yr_processed %>%
        select(year, sex, age, pop = pop_end) %>%
        mutate(
          age = ifelse(age == 100, 100, age + 1), # 100+ is open
          year = year + 1
        ) %>%
        summarise(pop = sum(pop), .by = c(year, sex, age)) %>%
        bind_rows(tibble(
          age = 0,
          sex = c("f", "m"),
          pop = 0,
          year = yr + 1
        )) %>%
        arrange(sex, age)
    }
  }

  bind_rows(results_list) %>% mutate(sim_id = sim_id)
}
