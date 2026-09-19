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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# following the missing across register releases ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Four releases of the register, each a full individual-level list with
# names and dates of birth (gitignored). 09 turns them into anonymous
# transition counts, 09f into per-person histories without names, 13d into
# counts by event month. Nothing these functions return may be written out
# with the key, name or date-of-birth columns still in it.
ual_releases <- tibble(
  release = c("v14", "v16", "v18", "v19"),
  date = as.Date(c("2025-09-16", "2025-12-04", "2026-04-23", "2026-09-19")),
  path = file.path("data_input/ualosses_hubert_datasets", c(
    "250916_UKR_ualosses_Personnel_v14.xlsx", "251204_UKR_ualosses_Personnel_v16.xlsx",
    "260423_UKR_ualosses_Personnel_v18.xlsx", "260919_UKR_ualosses_Personnel_v19.xlsx"
  ))
)
ual_status_order <- c("dead", "released_prisoner", "prisoner", "missing")

# What linkage needs from one release. A person's key is full name plus date
# of birth; `name` (surname and first name) serves the corrected-key search.
ual_read_release <- function(path) {
  read_xlsx(require_raw(path), sheet = "Database", col_types = "text") |>
    transmute(
      dob = excel_date(DateBirth),
      key = paste(LastName, FirstName, Patronym, dob),
      name = paste(LastName, FirstName),
      date_evnt = excel_date(DateEvent),
      year = year(date_evnt),
      ukrainian = ual_is_ukrainian(Nationality),
      status = Status
    )
}

# One row per key. Where a key appears under two statuses the most resolved
# is kept: a death record, or a return from captivity, is newer information
# than a missing one.
ual_one_per_key <- function(r) {
  r |> arrange(match(status, ual_status_order)) |> distinct(key, .keep_all = TRUE)
}

# A person not found under their key may be listed under a corrected name or
# date of birth. The candidates are the later release's NEW keys - absent from
# `before`, the release the search starts from - sharing surname and first
# name and either the date of birth or the month of the event. Where several
# qualify, the most resolved status is kept.
ual_corrected_key <- function(p, later, before) {
  candidates <-
    later |>
    anti_join(before, by = "key") |>
    transmute(name, dob_new = dob, month_new = floor_date(date_evnt, "month"),
              key_new = key, status_new = status)
  p |>
    transmute(key, name, dob, month = floor_date(date_evnt, "month")) |>
    inner_join(candidates, by = "name", relationship = "many-to-many") |>
    filter((dob == dob_new) %in% TRUE | (month == month_new) %in% TRUE) |>
    arrange(match(status_new, ual_status_order)) |>
    distinct(key, .keep_all = TRUE) |>
    select(key, key_new, status_new)
}

# Everyone listed as missing - event in 2022-2025, Ukrainian - in any release
# but the last, entering at their first listing, and their status in every
# later release: `found` says whether under the same key, a corrected key, or
# not at all. regs: the releases in order, named, as read by
# ual_read_release(). One row per person and later release; the person is
# identified by a sequential `pid`, and no key, name or date of birth is kept.
ual_follow_missing <- function(regs) {
  rel <- names(regs)
  lu <- map(regs, ual_one_per_key)
  seen <- character(0)    # keys listed in an earlier release
  claimed <- character(0) # corrected keys of people already followed
  out <- list()
  for (i in seq_len(length(rel) - 1)) {
    entrants <-
      regs[[i]] |>
      filter(status == "missing", year %in% 2022:2025, ukrainian) |>
      distinct(key, .keep_all = TRUE) |>
      filter(!key %in% seen, !key %in% claimed) |>
      mutate(pid = paste0(rel[i], "-", row_number()))
    cur <- entrants |> select(pid, key, name, dob, date_evnt)
    rows <- list()
    for (j in (i + 1):length(rel)) {
      at_j <- cur |> left_join(lu[[j]] |> select(key, status), by = "key")
      ck <- ual_corrected_key(at_j |> filter(is.na(status)), lu[[j]], lu[[j - 1]])
      at_j <-
        at_j |>
        left_join(ck, by = "key") |>
        mutate(found = case_when(!is.na(status) ~ "exact", !is.na(status_new) ~ "corrected",
                                 .default = "absent"),
               status = coalesce(status, status_new),
               key = coalesce(key_new, key))
      claimed <- c(claimed, at_j$key[at_j$found == "corrected"])
      rows[[j]] <- at_j |> transmute(pid, release = rel[j], status, found)
      cur <- at_j |> select(pid, key, name, dob, date_evnt)
    }
    out[[i]] <-
      entrants |>
      transmute(pid, year, date_evnt, entry = rel[i]) |>
      left_join(bind_rows(rows), by = "pid", relationship = "one-to-many")
    seen <- union(seen, lu[[i]]$key)
  }
  bind_rows(out)
}

# The windows between consecutive releases in which each person is at risk
# (missing at the window's start), up to their first resolution, with the
# outcome at the window's end. A person absent from a release but listed again
# in a later one is carried as still missing. A person absent from a release
# and from every later one is no longer listed, and counts as resolved alive
# at the first absence (absent_as = "alive", production) or, for the
# sensitivity, as still missing to the end (absent_as = "missing").
#
# A released prisoner was held, not disappeared. v19 is the first release to
# record returns from captivity, so a person listed as missing whose first
# resolution is a release was a prisoner the register had not recorded. Such
# a person is left out of the population at risk altogether, in every window
# (released = "exclude", production). The alternative, released = "alive",
# counts the return as a resolution to captivity; 09 compares the two.
ual_outcomes <- c("dead", "prisoner", "no_longer_listed", "missing")
ual_windows <- function(hist, absent_as = c("alive", "missing"), released = c("exclude", "alive")) {
  absent_as <- match.arg(absent_as)
  released <- match.arg(released)
  rel <- ual_releases$release
  w <-
    hist |>
    arrange(pid, match(release, rel)) |>
    mutate(listed_later = rev(cumsum(rev(!is.na(status)))) > 0, .by = pid) |>
    mutate(
      to = case_when(!is.na(status) ~ status,
                     listed_later | absent_as == "missing" ~ "missing",
                     .default = "no_longer_listed"),
      from_release = rel[match(release, rel) - 1]
    ) |>
    filter(cumsum(lag(to != "missing", default = FALSE)) == 0, .by = pid) |>
    select(pid, year, entry, date_evnt, from_release, release, to)
  if (released == "alive") return(w |> mutate(to = if_else(to == "released_prisoner", "prisoner", to)))
  w |> anti_join(w |> filter(to == "released_prisoner") |> distinct(pid), by = "pid")
}

# Twelve-month probabilities of moving from missing to each outcome, by
# event-year cohort, from window counts (year, from_release, to, n): a chain
# through the windows, P(outcome) = sum over windows of P(still missing at
# the window's start) x the window's share moving to that outcome.
ual_compose <- function(counts) {
  rel <- ual_releases$release
  counts |>
    summarise(n = sum(n), .by = c(year, from_release, to)) |>
    complete(nesting(year, from_release), to = ual_outcomes, fill = list(n = 0)) |>
    mutate(p = n / sum(n), .by = c(year, from_release)) |>
    select(-n) |>
    pivot_wider(names_from = to, values_from = p) |>
    arrange(year, match(from_release, rel)) |>
    mutate(s_before = lag(cumprod(missing), default = 1), .by = year) |>
    summarise(across(all_of(setdiff(ual_outcomes, "missing")), \(x) sum(s_before * x)),
              missing = prod(missing), .by = year) |>
    pivot_longer(-year, names_to = "to", values_to = "prop")
}

# The chain's states: a person no longer listed is alive; a prisoner is alive
# too, but kept apart.
ual_chain_state <- function(to) {
  case_when(to == "no_longer_listed" ~ "alive", .default = to)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# resolution of the missing by duration since disappearance ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How fast the missing are resolved depends on how long they have been
# missing. The cause-specific hazards of leaving "missing" - to death, to
# captivity, out of the register - are taken as constant within bands of
# months since the event and the same for every event month: duration, not
# cohort, carries the pattern. They are fitted by maximum likelihood to the
# windows between releases, one cell per event month and window, and each
# event month's missing are then carried forward from the duration they have
# reached in v19.
ual_duration_breaks <- c(0, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48, 60)
ual_resolutions <- c("dead", "prisoner", "no_longer_listed")

# months from the middle of an event month to a date
ual_months_since <- function(date, month) as.numeric(date - (month + 14)) / 30.44

# For someone missing at duration d0, the probability of each outcome by d1:
# still missing, or first resolved to each of ual_resolutions, under hazards
# h[band, outcome] per month (competing risks, piecewise constant).
ual_resolve <- function(h, d0, d1, breaks = ual_duration_breaks) {
  still <- rep(1, length(d0))
  p <- matrix(0, length(d0), ncol(h), dimnames = list(NULL, colnames(h)))
  for (b in seq_len(nrow(h))) {
    exposure <- pmax(0, pmin(d1, breaks[b + 1]) - pmax(d0, breaks[b]))
    total <- sum(h[b, ])
    leave <- 1 - exp(-exposure * total)
    if (total > 0) p <- p + outer(still * leave, h[b, ] / total)
    still <- still * (1 - leave)
  }
  cbind(missing = still, p)
}

# One row per event month and window between consecutive releases: the
# durations at its start and end, and how many of those missing at its start
# are in each outcome at its end. counts: month, from_release, to, n.
ual_duration_cells <- function(counts) {
  rel_date <- set_names(ual_releases$date, ual_releases$release)
  next_rel <- set_names(ual_releases$release[-1], ual_releases$release[-nrow(ual_releases)])
  counts |>
    summarise(n = sum(n), .by = c(month, from_release, to)) |>
    complete(nesting(month, from_release), to = c("missing", ual_resolutions), fill = list(n = 0)) |>
    pivot_wider(names_from = to, values_from = n) |>
    mutate(d0 = pmax(0, ual_months_since(rel_date[from_release], month)),
           d1 = ual_months_since(rel_date[next_rel[from_release]], month)) |>
    filter(d1 > d0)
}

# Maximum-likelihood hazards from the cells: each cell's counts are
# multinomial over still missing and the three resolutions, with the
# probabilities ual_resolve() gives for its durations.
#
# A release can record one kind of resolution in a batch - v19 added
# prisoners and removed records wholesale - which a duration pattern alone
# would read as belonging to whatever durations were at risk in that window.
# So each window has its own multiplier per outcome (the first window is the
# reference), and the projection uses their average over the year the windows
# span, weighted by each window's length: the hazard at duration d is
# h(d) x that average.
#
# theta holds the log-hazards by band and outcome, then the log-multipliers of
# the later windows by outcome. Log-hazards are bounded below, so an outcome
# never seen in a band gets a hazard of effectively zero and no variance.
# horizon: the duration the projection runs to.
ual_window_months <- diff(ual_releases$date) |> as.numeric() / 30.44
ual_theta_hazards <- function(theta, nb, k, windows = ual_window_months) {
  nw <- length(windows)
  h <- matrix(exp(theta[seq_len(nb * k)]), nb, k, dimnames = list(NULL, ual_resolutions))
  mult <- rbind(1, matrix(exp(theta[nb * k + seq_len((nw - 1) * k)]), nw - 1, k))
  list(h = h, mult = mult,
       projection = sweep(h, 2, colSums(mult * windows) / sum(windows), "*"))
}
ual_fit_durations <- function(cells, horizon, breaks = ual_duration_breaks) {
  nb <- length(breaks) - 1
  k <- length(ual_resolutions)
  nw <- length(ual_window_months)
  window <- match(cells$from_release, ual_releases$release)
  y <- as.matrix(cells[, c("missing", ual_resolutions)])
  cell_probs <- function(theta) {
    hz <- ual_theta_hazards(theta, nb, k)
    p <- matrix(0, nrow(cells), k + 1, dimnames = list(NULL, c("missing", ual_resolutions)))
    for (w in seq_len(nw)) {
      i <- window == w
      p[i, ] <- ual_resolve(sweep(hz$h, 2, hz$mult[w, ], "*"), cells$d0[i], cells$d1[i], breaks)
    }
    p
  }
  nll <- function(theta) -sum(y * log(pmax(cell_probs(theta)[, colnames(y)], 1e-300)))
  n_par <- nb * k + (nw - 1) * k
  fit <- optim(c(rep(log(0.003), nb * k), rep(0, (nw - 1) * k)), nll, method = "L-BFGS-B",
               lower = c(rep(-20, nb * k), rep(-10, (nw - 1) * k)),
               upper = c(rep(0, nb * k), rep(10, (nw - 1) * k)),
               hessian = TRUE, control = list(maxit = 5000, factr = 1e5))
  stopifnot(fit$convergence == 0)
  free <- c(fit$par[seq_len(nb * k)] > -19.9, rep(TRUE, (nw - 1) * k))
  vcov <- matrix(0, n_par, n_par)
  vcov[free, free] <- solve(fit$hessian[free, free])
  hz <- ual_theta_hazards(fit$par, nb, k)
  list(breaks = breaks, h = hz$projection, h_reference = hz$h, window_multipliers = hz$mult,
       theta = fit$par, vcov = vcov, nb = nb, k = k, horizon = horizon,
       resolved = colSums(y[, ual_resolutions]),
       fitted = cells |> bind_cols(as_tibble(cell_probs(fit$par), .name_repair = \(x) paste0("p_", x))))
}

# Imputation of the long-term missing. Defined once so that 09 (the central
# value of alpha), 13b (the sweep over alpha, and rate sampling), 13d and the
# range computed by alpha_evidence() below all use the SAME formula.
#
# model          the fitted durations (ual_fit_durations()): hazards by band
#                and the horizon
# stock_missing  tibble(month, year, missing_stock): v19's missing by event
#                month, completed for registration lag (09)
# alpha          share of the NEVER-RESOLVED long-term missing who are alive;
#                1 - alpha are dead. The one free parameter: its range comes
#                from alpha_evidence(), 11 draws it, 13b sweeps it.
#
# Each event month's missing move from the duration they have reached in v19
# to the horizon, resolving along the way to death, to captivity or out of the
# register; those still missing at the horizon are the never-resolved, and a
# share alpha of them is alive. A month already past the horizon is not moved.
# Returns, by event year (or by = c("month", "year")), missing_stock and its
# parts: projected_dead, imputed_prisoner (projected captivity),
# imputed_unlisted (projected to leave the register), never_resolved, and of
# these alive_never_resolved (alpha x never_resolved); imputed_dead is
# projected_dead plus the dead share of the never-resolved, imputed_alive the
# unlisted plus the alive share. Every output is linear in alpha.
impute_missing <- function(alpha, model, stock_missing, by = "year") {
  d_v19 <- ual_months_since(ual_releases$date[nrow(ual_releases)], stock_missing$month)
  p <- ual_resolve(model$h, pmin(d_v19, model$horizon), model$horizon, model$breaks)
  stock_missing |>
    mutate(
      projected_dead = missing_stock * p[, "dead"],
      imputed_prisoner = missing_stock * p[, "prisoner"],
      imputed_unlisted = missing_stock * p[, "no_longer_listed"],
      never_resolved = missing_stock * p[, "missing"],
      alive_never_resolved = never_resolved * alpha,
      imputed_dead = projected_dead + never_resolved - alive_never_resolved,
      imputed_alive = imputed_unlisted + alive_never_resolved
    ) |>
    summarise(across(c(missing_stock, imputed_dead, imputed_alive, imputed_prisoner, imputed_unlisted,
                       projected_dead, never_resolved, alive_never_resolved), sum),
              .by = all_of(by))
}

# The earlier form of the chain, kept for comparison (09): twelve-month steps,
# with event-year cohorts, one year apart, standing in for duration since
# disappearance. tasas_long holds the twelve-month probabilities by cohort
# and state (alive, dead, prisoner, missing) in `prop`, and stock_missing the
# missing by event year.
impute_missing_cohorts <- function(alpha, tasas_long, stock_missing) {
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
#   official      the prisoners of war held in February 2026 plus the
#   figures       military personnel returned from captivity by then:
#                 everyone taken prisoner, whether held or released since.
#                 Both counts are of the military at the same date, so no one
#                 is counted twice and no civilian is counted
#
# Prisoners the register does not record as such are either among its missing
# or not in the register at all.
#
#   min   0 - every unrecorded prisoner is outside the register (most of the
#         people it records as returned were in no earlier release; 09)
#   mode  (share_among_missing x (pow_held + returned_military - register_alive)
#          - projected captivity) / residual, and no lower than 0. With
#         share_among_missing = 1 (production), every unrecorded prisoner is
#         among the missing, and those the chain itself projects into
#         captivity are not among the never-resolved
#   max   the register's own share of resolutions that are alive (to
#         prisoner or no longer listed), assuming the never-resolved are no
#         more often alive than those resolved: a living prisoner is listed
#         or exchanged more readily than a body is recovered
#
# Each unrecorded prisoner is counted once. Those among today's missing whom
# the chain projects into captivity are resolved alive by the chain, so they
# are taken out of the unrecorded prisoners before the rest are set against
# the never-resolved. The chain resolves nobody to a release - released
# prisoners are outside the population at risk (ual_windows) - so a prisoner
# later exchanged is left to alpha. The mode leans high: of the people the
# register records as returned who were not recorded as prisoners before,
# about a quarter had been listed as missing and the rest were in no earlier
# release, and 09 shows the mode with the unrecorded prisoners among the
# missing in that proportion. The returned figure also holds the few taken
# prisoner before 2022, outside the register's event years.
#
# residual is the number of missing never resolved at the end of the chain,
# the quantity alpha applies to: impute_missing() is linear in alpha, so it is
# the imputed dead at alpha = 0 minus those at alpha = 1.
#
# SOURCES
#   pow_held           "about 7,000 Ukrainian prisoners of war" held by
#                      Russia - President Zelensky, Munich, 14 Feb 2026
#                      (Ukrinform, 14 Feb 2026).
#   returned_military  7,291 military personnel (235 women, 7,056 men)
#                      returned from captivity since 24 Feb 2022, in the
#                      Coordination Headquarters for the Treatment of Prisoners
#                      of War's report for February 2026
#                      (koordshtab.gov.ua/report, read 19 Sep 2026). The 876
#                      civilians returned by then are left out: the register
#                      lists military personnel only.
pow_held <- 7000
returned_military <- 7291

# impute          function(alpha) giving impute_missing()'s output for the
#                 chain in use
# resolved        the first resolutions the chain rests on, as counts named
#                 dead and by the alive outcomes (prisoner, no longer listed)
# register_alive  the register's prisoners plus released prisoners, events
#                 2022-2025, from 08's redistributed counts
# net_projected_captivity  whether the chain's own projected transitions to
#                 captivity are taken out of the unrecorded prisoners before
#                 they are set against the never-resolved
# share_among_missing  the share of the unrecorded prisoners taken to be among
#                 the register's missing rather than outside it
alpha_evidence <- function(impute, resolved, register_alive, net_projected_captivity = TRUE,
                           share_among_missing = 1) {
  imp0 <- impute(0)
  imp1 <- impute(1)
  residual <- sum(imp0$imputed_dead - imp1$imputed_dead)
  projected_captivity <- sum(imp0$imputed_prisoner)

  resolved_alive <- sum(resolved[names(resolved) != "dead"])
  resolved_all <- sum(resolved)
  unrecorded <- pow_held + returned_military - register_alive
  among_never_resolved <- share_among_missing * unrecorded -
    if (net_projected_captivity) projected_captivity else 0

  out <- tibble(
    alpha_min = 0,
    alpha_mode = max(0, among_never_resolved) / residual,
    alpha_max = resolved_alive / resolved_all,
    residual = residual,
    unrecorded_prisoners = unrecorded,
    projected_captivity = projected_captivity,
    net_projected_captivity = net_projected_captivity,
    share_among_missing = share_among_missing,
    register_alive = register_alive,
    resolved_alive = resolved_alive,
    resolved_dead = resolved_all - resolved_alive,
    pow_held = pow_held,
    returned_military = returned_military
  )
  stopifnot(
    out$alpha_mode >= out$alpha_min,
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
#   1000   CONSTRUCTION. ~2 min in 11. The bootstrap standard error of the
#          reported quantiles is about 0.02 years in the worst cell, the upper
#          bound of the 2022 male loss, so individual bounds move in the
#          second decimal between runs. Because the seed is fixed, two runs of
#          the SAME size draw the same PERT quantiles, so differences between
#          pipeline states are far more precise than either state on its own -
#          which is what makes this size usable for comparing changes.
#
#   20000  PRODUCTION (current), for the estimates that get reported. About
#          half an hour in 11. Brings that worst cell to about 0.004 years, so
#          bounds are reproducible to about +/-0.01 years; strict
#          reproducibility of the second decimal would need about 60,000.
#
# Nothing silently depends on remembering the size: the cache file name
# carries it, 15 refuses a cache that does not match, and
# tables/_run_provenance.csv records what every table was built from.
n_sim <- 20000

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
# PERIOD INPUTS ON COHORT ROWS. Rates and age-at-death profiles (mx, fx,
# OHCHR, ualosses) are by completed age at the event, while row x spends year
# t aged x - 1 before its birthday and x after. So each is carried onto the
# rows as the cohort lives it: row x (x >= 1) takes the mean of the rates at
# ages x - 1 and x, and half of the deaths at each of those ages; row 0, born
# during the year, takes the age-0 rate and half of the age-0 deaths; the open
# row 100 keeps what would fall to row 101. The output is returned by
# completed age again: each row's exposure is split half to each of its two
# ages (row 0 wholly to age 0), expected deaths are the rate at each age times
# that exposure, and conflict deaths are the counts by age at death, so the
# life tables of 12 and 14 are period tables that reproduce the inputs.
#
# TIMING: exposure is the average of the cohort's stock at the start and end
# of the year. For existing cohorts, a share of the year's net outflow and of
# its conflict deaths is removed before expected deaths are computed, and half
# the expected deaths after. The shares - absent_mig, absent_cvs, absent_cmb -
# are 1/2 (the mid-year convention) in 2023-2025 and, in 2022, the timing 10
# measures from the months of the events (data_inter/ukr_timing_absent.rds);
# draws_this_sim may carry its own, as 13f does for the mid-year convention.
# Newborns start the year at zero and arrive through it, so their exposure is
# half of what survives.
#
# The three conversions between completed age (0-100, 100 open) and cohort
# rows, for one year and sex ordered by age:
rates_to_rows <- function(v) c(v[1], (head(v, -1) + tail(v, -1)) / 2)
counts_to_rows <- function(v) {
  n <- length(v)
  r <- c(v / 2, 0) + c(0, v / 2)
  r[n] <- r[n] + r[n + 1]
  r[seq_len(n)]
}
rows_to_ages <- function(e) {
  n <- length(e)
  a <- c(e[1], e[-1] / 2)
  a[-n] <- a[-n] + e[-1] / 2
  a
}
# The timing shares by year (see TIMING), read once from 10's output.
absent_timing <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      f <- "data_inter/ukr_timing_absent.rds"
      if (!file.exists(f)) stop("Run step 10 first: ", f, " is missing.", call. = FALSE)
      cache <<- read_rds(f)$absent
    }
    cache
  }
})
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

  timing_cols <- setdiff(c("absent_mig", "absent_cvs", "absent_cmb"), names(draws_this_sim))
  if (length(timing_cols) > 0) {
    draws_this_sim <- draws_this_sim |>
      left_join(absent_timing() |> select(year, all_of(timing_cols)), by = "year")
  }

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
      cmb_imputed = draw_cmb * (1 - ratio_confirmed) * prop_cmb_miss
    ) %>%
    # the period inputs carried onto the cohort rows (see AGE CONVENTION)
    arrange(year, sex, age) %>%
    mutate(
      mx_row = rates_to_rows(mx),
      fx_row = rates_to_rows(fx),
      cvs_row = counts_to_rows(cvs),
      cmb_row = counts_to_rows(cmb_confirmed + cmb_imputed),
      cnf_row = cvs_row + cmb_row,
      .by = c(year, sex)
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
        # population at risk
        pop2 = pop - absent_mig * ems - absent_cvs * cvs_row - absent_cmb * cmb_row,
        # expected (non-conflict) deaths. mx is a rate per PERSON-YEAR (04
        # fits it to deaths over the mean of consecutive 1 January stocks),
        # so deaths must equal mx x exposure. exposure = pop2 - noc / 2, which
        # solves to noc = pop2 * mx / (1 + mx / 2). Using pop2 * mx instead
        # would apply mx / (1 - mx / 2): 6% too high at 80, 21% at 100.
        noc = pop2 * mx_row / (1 + 0.5 * mx_row),
        exposure = pop2 - 0.5 * noc, # person-years lived

        # row 0 is the cohort born during this year: it starts empty and
        # receives births = sum(ASFR x female exposure), split by the sex
        # ratio at birth. Births never touch the other rows.
        births = sum(fx_row * exposure),
        pop = case_when(
          age == 0 & sex == "f" ~ prop_female_birth * births,
          age == 0 & sex == "m" ~ prop_male_birth * births,
          .default = pop
        ),
        # same bookkeeping for the newborns, but they start the year at zero:
        # the average of the start (0) and end stock is half of what survives
        pop2 = ifelse(age == 0, 0.5 * (pop - ems - cnf_row), pop2),
        noc = ifelse(age == 0, pop2 * mx_row / (1 + 0.5 * mx_row), noc),
        exposure = ifelse(age == 0, pop2 - 0.5 * noc, exposure),

        dx = noc + cnf_row, # all-cause deaths
        pop_end = pop - ems - dx # survivors at 31 December
      )

    # keep the components only, back by completed age (see AGE CONVENTION);
    # every rate downstream is derived from these
    results_list[[as.character(yr)]] <- yr_processed %>%
      arrange(sex, age) %>%
      mutate(exposure_age = rows_to_ages(exposure), .by = sex) %>%
      select(
        year, sex, age,
        pop = exposure_age,
        expected = mx,
        civilian = cvs,
        combatant_confirmed = cmb_confirmed,
        combatant_imputed = cmb_imputed
      ) %>%
      mutate(expected = expected * pop)

    # age the survivors on into next year's starting population
    if (yr < 2025) {
      current_pop_ini <- age_on(yr_processed %>% select(sex, age, pop = pop_end), yr + 1)
    }
  }

  bind_rows(results_list) %>% mutate(sim_id = sim_id)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the simulation's draws ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Every disputed input of one set of n simulations, drawn in a fixed order from
# one seed so that 11 and the PERT-shape check in 13g use the same random
# numbers. Each draw is a Beta-PERT quantile of a uniform, with the given
# shape (4 is the standard PERT):
#   civilians   one uniform per year and simulation: each year's count comes
#               from that year's own event reporting
#   combatants  alpha once per simulation, which sets every year's total on
#               09's line: at_alpha0 - alpha x residual
#   migration   the western blend weight once per simulation, applied to every
#               year's two readings; Russia and Belarus once
# Returns the draws by simulation and year (draws_df), the same by role
# (param_draws) and alpha by simulation (alpha_draws).
simulation_draws <- function(param_table, alpha_range, alpha_lines, n, shape = 4, seed = 42) {
  set.seed(seed)
  q <- function(u, a, m, b) qpert(u, min = a, mode = m, max = b, shape = shape)

  cvs <- param_table |> filter(role == "civilians") |> arrange(year)
  draws_cvs_long <- map_dfr(seq_len(nrow(cvs)), function(i) {
    tibble(role = "civilians", year = cvs$year[i], sim_id = seq_len(n),
           draw = q(runif(n), cvs$min[i], cvs$mode[i], cvs$max[i]))
  })

  alpha_draws <- tibble(
    sim_id = seq_len(n),
    alpha = q(runif(n), alpha_range$alpha_min, alpha_range$alpha_mode, alpha_range$alpha_max)
  )
  draws_cmb_long <-
    expand_grid(alpha_draws, alpha_lines |> select(year, at_alpha0, residual)) |>
    transmute(role = "combatants", year, sim_id, draw = at_alpha0 - alpha * residual)

  west <- param_table |> filter(role == "mig_west")
  w <- q(runif(n), 0, 0.5, 1)
  draws_west_long <-
    expand_grid(sim_id = seq_len(n), west |> select(year, crossings, register)) |>
    mutate(role = "mig_west", draw = crossings + w[sim_id] * (register - crossings)) |>
    select(role, year, sim_id, draw)
  ru <- param_table |> filter(role == "mig_ru_by")
  u_ru <- runif(n)
  draws_ru_long <-
    expand_grid(sim_id = seq_len(n), ru |> select(year, min, mode, max)) |>
    mutate(role = "mig_ru_by", draw = q(u_ru[sim_id], min, mode, max)) |>
    select(role, year, sim_id, draw)

  param_draws <- bind_rows(draws_cvs_long, draws_cmb_long, draws_west_long, draws_ru_long)
  draws_df <-
    param_draws |>
    mutate(role = case_when(role == "civilians" ~ "draw_cvs", role == "combatants" ~ "draw_cmb",
                            .default = "draw_mig")) |>
    summarise(draw = sum(draw), .by = c(sim_id, year, role)) |>
    pivot_wider(names_from = role, values_from = draw) |>
    left_join(param_table |> filter(role == "combatants") |> select(year, conf_cmb = confirmed),
              by = "year") |>
    arrange(sim_id, year)
  list(draws_df = draws_df, param_draws = param_draws, alpha_draws = alpha_draws)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the deterministic projection at the mode, for the sensitivity steps ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Every input at its mode, set up as 13 sets it up: the conflict and migration
# totals from the parameter table, the counterfactual, fertility carried
# forward from 2023, and the age-sex profiles of civilian, registered and
# imputed deaths. A sensitivity step changes one of these and projects again.
mode_projection_inputs <- function() {
  param_table <- read_rds("data_inter/ukr_param_table.rds")
  exp_mort <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") |>
    filter(year %in% 2022:2025, source == "frcst") |>
    select(-source)
  pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") |>
    filter(reg == "cnt", year == 2022) |>
    select(-reg)
  migs <- read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds") |>
    mutate(ems = -mix) |>
    select(-mix)
  asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds")
  asfr <- bind_rows(asfr, asfr |> filter(year == 2023) |> mutate(year = 2024),
                    asfr |> filter(year == 2023) |> mutate(year = 2025))
  ohchr <- read_rds("data_inter/ukr_ohchr_civilian_casualties.rds") |>
    select(year, sex, age, prop_cvs = cx)
  ual <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")
  profile <- function(keep, nm) {
    ual |>
      filter(status == keep, year %in% 2022:2025) |>
      select(-status) |>
      complete(year = 2022:2025, sex, age = 0:100, fill = list(dx = 0)) |>
      summarise(dx = sum(dx), .by = c(year, sex, age)) |>
      mutate("{nm}" := dx / sum(dx), .by = year) |>
      select(-dx)
  }
  static_inputs <-
    expand_grid(year = 2022:2025, sex = c("f", "m"), age = 0:100) |>
    left_join(migs, by = c("year", "sex", "age")) |>
    left_join(exp_mort, by = c("year", "sex", "age")) |>
    left_join(asfr, by = c("year", "sex", "age")) |>
    left_join(profile("dead", "prop_cmb_dead"), by = c("year", "sex", "age")) |>
    left_join(profile("missing", "prop_cmb_miss"), by = c("year", "sex", "age")) |>
    left_join(ohchr, by = c("year", "sex", "age")) |>
    replace_na(list(ems = 0, w = 0, mx = 0, fx = 0, prop_cvs = 0,
                    prop_cmb_dead = 0, prop_cmb_miss = 0))
  stopifnot(all(static_inputs$mx > 0), !any(is.na(static_inputs)))
  draws <-
    param_table |>
    mutate(role = if_else(str_starts(role, "mig_"), "migration", role)) |>
    summarise(mode = sum(mode), .by = c(year, role)) |>
    pivot_wider(names_from = role, values_from = mode) |>
    transmute(year, draw_cmb = combatants, draw_cvs = civilians, draw_mig = migration) |>
    left_join(param_table |> filter(role == "combatants") |> select(year, conf_cmb = confirmed),
              by = "year") |>
    mutate(sim_id = 1)
  list(static_inputs = static_inputs, pop22_ini = pop22_ini, draws = draws,
       param_table = param_table)
}

# The loss of life expectancy at birth by year and sex for one set of draws,
# with the counterfactual and the with-war life expectancy it rests on.
loss_at_mode <- function(draws, static_inputs, pop22_ini) {
  sim <-
    run_single_sim(1, draws, static_inputs, pop22_ini) |>
    mutate(mx_all = (expected + civilian + combatant_confirmed + combatant_imputed) / pop,
           mx_bsn = expected / pop)
  stopifnot(min(sim$pop) >= 0)
  e0 <- function(col) {
    sim |>
      select(year, sex, age, mx = all_of(col)) |>
      group_by(year, sex) |>
      do(lifetable(dt_in = .data)) |>
      ungroup() |>
      filter(age == 0) |>
      select(year, sex, ex)
  }
  e0("mx_bsn") |>
    rename(e0_bsn = ex) |>
    left_join(e0("mx_all") |> rename(e0_all = ex), by = c("year", "sex")) |>
    mutate(loss = e0_bsn - e0_all)
}
