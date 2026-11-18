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
  "mc2d",
  "uncertainUCDP"
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
# Six releases of the register, every one between the first and the last the
# project holds, each a full individual-level list with names and dates of
# birth (gitignored). 09 turns them into anonymous transition counts, 09f into
# per-person histories without names, 08b into counts by event month. Nothing
# these functions return may be written out with the key, name or
# date-of-birth columns still in it.
# The releases are Kaggle versions 14 to 19 of Olivier Hubert's "Confirmed
# Ukrainian military personnel losses". The file names carry the date each file
# was saved.
#
# Each release is dated by the data it holds, not by the day it was published.
# What the lag of an event and the length of a window need is the day the
# register was copied, which is not recorded. Two dates bracket it: the latest
# event the release lists (latest_event) and its publication on Kaggle
# (version). The interval between the two runs from 4 days (v16) to 30 (v17)
# and 33 (v19), while in every release the daily counts fall away over the last
# one or two weeks before the latest event, as they do in a list copied a few
# days after it. A release published a month after its latest event was copied
# long before it was published: v19 lists nothing after 18 June 2026, and the
# register itself, read on 1 October 2026, stops at 21 June with about a
# hundred more people than v19 holds. So every release is dated by its latest
# event plus the shortest interval to publication observed, 4 days: the latest
# day on which all six can have been copied if the register lists events as
# promptly in each. Dating the releases by their Kaggle versions instead runs
# the clock of v19 a month long; 08b measures what that does to the lag
# correction (13d, Table A20). 08b also checks the latest events against the
# files.
# Six releases give five windows. What they show that four could not: the
# register's clean-ups come in some releases and not others - the dead and the
# missing both leave it in v15, v17 and v19 and hardly at all in between - so a
# rate measured over a window that spans two releases averages a clean-up with
# a lull.
ual_releases <- tibble(
  release = c("v14", "v15", "v16", "v17", "v18", "v19"),
  version = as.Date(c("2025-09-16", "2025-11-02", "2025-12-04",
                      "2026-02-06", "2026-04-23", "2026-07-21")),
  latest_event = as.Date(c("2025-09-09", "2025-10-16", "2025-11-30",
                           "2026-01-07", "2026-04-10", "2026-06-18")),
  path = file.path("data_input/ualosses_hubert_datasets", c(
    "250916_UKR_ualosses_Personnel_v14.xlsx", "251102_UKR_ualosses_Personnel_v15.xlsx",
    "251204_UKR_ualosses_Personnel_v16.xlsx", "260206_UKR_ualosses_Personnel_v17.xlsx",
    "260423_UKR_ualosses_Personnel_v18.xlsx", "260919_UKR_ualosses_Personnel_v19.xlsx"
  ))
) |>
  mutate(date = latest_event + min(version - latest_event), .after = release)
ual_status_order <- c("dead", "released_prisoner", "prisoner", "missing")

# What linkage needs from one release. A person's key is full name plus date
# of birth; `name` (surname and first name), and the two apart, serve the
# corrected-key search.
ual_read_release <- function(path) {
  read_xlsx(require_raw(path), sheet = "Database", col_types = "text") |>
    transmute(
      dob = excel_date(DateBirth),
      key = paste(LastName, FirstName, Patronym, dob),
      name = paste(LastName, FirstName),
      surname = LastName,
      first = FirstName,
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

# The two rules by which a record listed under another key is the same person:
#   name    the same surname and first name, and the same date of birth or the
#           same month of the event: a corrected patronymic, date of birth or
#           event date
#   dates   the same date of birth AND the same event date, and the same
#           surname or the same first name: a re-spelt name. The register
#           re-spells names in batches - of the records a release drops and the
#           first rule cannot follow, this one finds 72, 11, 259, 10 and 357 of
#           the missing and 53, 10, 427, 12 and 470 of the dead over the five
#           windows, nearly all of them a surname spelt another way. A date of
#           birth and an event date coincide by chance too seldom to matter:
#           with the event date moved by one or two weeks the rule finds no one.
# x: the records looked for; y: the candidates, with their status. Both carry
# key, surname, first, dob, date_evnt and the column `name` names, the name the
# first rule compares (08b adds the patronymic to it for the dead). Returns one
# row per matching pair: key (of x), key_new and status_new (of y) and the rule.
ual_match_keys <- function(x, y, name = "name") {
  x <- x |> transmute(key, nm = .data[[name]], surname, first, dob, date_evnt)
  y <- y |> transmute(nm = .data[[name]], surname_new = surname, first_new = first, dob_new = dob,
                      date_new = date_evnt, key_new = key, status_new = status)
  month <- \(d) floor_date(d, "month")
  bind_rows(
    inner_join(x, y, by = "nm", relationship = "many-to-many") |>
      filter((dob == dob_new) %in% TRUE | (month(date_evnt) == month(date_new)) %in% TRUE) |>
      mutate(rule = "name"),
    inner_join(x |> filter(!is.na(dob), !is.na(date_evnt)), y,
               by = c("dob" = "dob_new", "date_evnt" = "date_new"), relationship = "many-to-many") |>
      filter((surname == surname_new | first == first_new) %in% TRUE) |>
      mutate(rule = "dates")
  ) |>
    select(key, key_new, status_new, rule)
}

# A person not found under their key may be listed under a corrected one. The
# candidates are the later release's NEW keys - absent from `before`, the
# release the search starts from - that either rule of ual_match_keys() accepts.
# Where several qualify, the most resolved status is kept.
ual_corrected_key <- function(p, later, before) {
  ual_match_keys(p, later |> anti_join(before, by = "key")) |>
    arrange(match(status_new, ual_status_order), match(rule, c("name", "dates"))) |>
    distinct(key, .keep_all = TRUE)
}

# Everyone listed as missing - event in 2022-2025, Ukrainian - in any release
# but the last, entering at their first listing, and their status in every
# later release: `found` says whether under the same key, a corrected key
# (`rule`: by which rule of ual_match_keys()), or not at all. regs: the
# releases in order, named, as read by
# ual_read_release(). One row per person and later release; the person is
# identified by a sequential `pid`, and no key, name or date of birth is kept.
# keep_keys = TRUE adds the key at entry and the key found at each release,
# for the clerical review only (09h), whose output stays local.
ual_follow_missing <- function(regs, keep_keys = FALSE) {
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
    cur <- entrants |> select(pid, key, name, surname, first, dob, date_evnt)
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
      rows[[j]] <- at_j |> transmute(pid, release = rel[j], status, found, rule,
                                      key = if_else(found == "absent", NA_character_, key))
      cur <- at_j |> select(pid, key, name, surname, first, dob, date_evnt)
    }
    out[[i]] <-
      entrants |>
      transmute(pid, year, date_evnt, entry = rel[i], entry_key = key) |>
      left_join(bind_rows(rows), by = "pid", relationship = "one-to-many")
    seen <- union(seen, lu[[i]]$key)
  }
  out <- bind_rows(out)
  if (keep_keys) out else out |> select(-entry_key, -key)
}

# The windows between consecutive releases in which each person is at risk
# (missing at the window's start), up to their first resolution, with the
# outcome at the window's end. A person absent from a release but listed again
# in a later one is carried as still missing. A person absent from a release
# and from every later one is no longer listed, and counts as resolved alive
# at the first absence (absent_as = "alive", production) or, for the
# sensitivity, as still missing to the end (absent_as = "missing").
#
# A missing person whose first resolution is a return from captivity was a
# prisoner the register had not recorded, and is alive: the return counts as a
# resolution to captivity (released = "alive", production). v19 is the first
# release to record returns, so these resolutions fall in its window. The
# alternative, released = "exclude", leaves such a person out of the
# population at risk altogether, as someone held rather than disappeared; 09
# compares the two.
ual_outcomes <- c("dead", "prisoner", "no_longer_listed", "missing")
ual_windows <- function(hist, absent_as = c("alive", "missing"), released = c("alive", "exclude")) {
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
  # A parameter estimated at or near its bound has no usable normal
  # approximation, so it is held at its estimate, with no variance: a
  # log-hazard below -12 (a monthly hazard under six in a million: an outcome
  # effectively never seen in that band), or a window multiplier at either bound
  # (an outcome the window records almost never, or in one batch).
  log_mult <- fit$par[nb * k + seq_len((nw - 1) * k)]
  free <- c(fit$par[seq_len(nb * k)] > -12, log_mult > -9.9 & log_mult < 9.9)
  vcov <- matrix(0, n_par, n_par)
  vcov[free, free] <- solve(fit$hessian[free, free])
  # The cells are fitted as multinomial, but event months differ by more than
  # that allows - a release records resolutions in batches, some months and not
  # others - so the counts are overdispersed and the inverse Hessian alone
  # understates the error. The covariance is scaled by the Pearson dispersion,
  # as a quasi-likelihood fit would scale it, and never below one. Cells expected
  # to hold less than one person are left out of the statistic, where it is
  # unstable.
  expected <- cell_probs(fit$par)[, colnames(y)] * rowSums(y)
  pearson <- ((y - expected)^2 / expected)[expected > 1]
  dispersion <- max(1, sum(pearson) / (nrow(y) * k - sum(free)))
  vcov <- vcov * dispersion
  hz <- ual_theta_hazards(fit$par, nb, k)
  list(breaks = breaks, h = hz$projection, h_reference = hz$h, window_multipliers = hz$mult,
       theta = fit$par, vcov = vcov, dispersion = dispersion, nb = nb, k = k, horizon = horizon,
       resolved = colSums(y[, ual_resolutions]),
       fitted = cells |> bind_cols(as_tibble(cell_probs(fit$par), .name_repair = \(x) paste0("p_", x))))
}

# Imputation of the missing. Defined once so that 09 (the central values),
# 10 (the bounds), 11 (every draw, through military_draws() below), 13b (the
# sweep over the missing alive) and 13d all use the SAME formula.
#
# model          the fitted durations (ual_fit_durations()): hazards by band
#                and the horizon
# stock_missing  tibble(month, year, missing_stock): v19's missing by event
#                month, completed for registration lag (09)
# captives       the prisoners of war among the missing, alive: set by the
#                official figures (alive_evidence() below), not projected
# alive_other    the share of the unresolved who are alive for any other
#                reason, as far as the rules can tell
#
# Each event month's missing move from the duration they have reached in v19
# to the horizon, resolving along the way to death, to captivity or out of the
# register (no longer listed, counted as found alive), or staying missing. A
# month already past the horizon is not moved. The projection of captivity
# rests on one release's catch-up of past returns (alive_evidence()), so
# captivity is not taken from it: the captives are spread over event months
# as the projection spreads captivity and scaled to the number given. The rest
# of the projected captivity and those still missing at the horizon are the
# unresolved; alive_other of them are alive and the rest dead.
#
# Returns, by event year (or by = c("month", "year")), missing_stock and its
# parts: projected_dead, projected_captivity, captives, imputed_unlisted
# (projected to leave the register), never_resolved, unresolved (projected
# captivity plus never-resolved, less the captives), alive_other
# (alive_other x unresolved), imputed_dead (projected deaths plus the dead
# among the unresolved) and imputed_alive (captives, no longer listed and
# alive_other).
impute_missing <- function(model, stock_missing, captives, alive_other = 0, by = "year") {
  d_v19 <- ual_months_since(ual_releases$date[nrow(ual_releases)], stock_missing$month)
  p <- ual_resolve(model$h, pmin(d_v19, model$horizon), model$horizon, model$breaks)
  proj <-
    stock_missing |>
    mutate(
      projected_dead = missing_stock * p[, "dead"],
      projected_captivity = missing_stock * p[, "prisoner"],
      imputed_unlisted = missing_stock * p[, "no_longer_listed"],
      never_resolved = missing_stock * p[, "missing"]
    )
  scale <- captives / sum(proj$projected_captivity)
  out <-
    proj |>
    mutate(
      captives = projected_captivity * scale,
      unresolved = projected_captivity + never_resolved - captives,
      alive_other = unresolved * alive_other,
      imputed_dead = projected_dead + unresolved - alive_other,
      imputed_alive = captives + imputed_unlisted + alive_other
    )
  stopifnot(all(out$unresolved >= 0))
  out |>
    summarise(across(c(missing_stock, imputed_dead, imputed_alive, projected_dead, projected_captivity,
                       captives, imputed_unlisted, never_resolved, unresolved, alive_other), sum),
              .by = all_of(by))
}

# The earlier form of the chain, kept for comparison (09): twelve-month steps,
# with event-year cohorts, one year apart, standing in for duration since
# disappearance. tasas_long holds the twelve-month probabilities by cohort
# and state (alive, dead, prisoner, missing) in `prop`, and stock_missing the
# missing by event year. alive_other: the share of those still missing at the
# end of the chain who are alive, as in impute_missing().
impute_missing_cohorts <- function(alive_other, tasas_long, stock_missing) {
  suelo_vivo <- alive_other
  suelo_muerto <- 1 - alive_other

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
# the missing who are alive: the evidence ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# By the rules of the linkage, a missing person is alive if they turn out to
# be a prisoner of war, held or since released, or if they leave the register
# (no longer listed). The model projects the second. The first is set by the
# official prisoner-of-war figures, compared with the register:
#
#   official   the prisoners of war held in February 2026 plus the military
#   figures    personnel returned from captivity by then: everyone taken
#              prisoner, held or released since. Both counts are of the
#              military at the same date, so no one is counted twice and no
#              civilian is counted
#   register   the prisoners and released prisoners v19 records
#
# The difference, the unrecorded prisoners, are either among the register's
# missing or not in it at all. The share among the missing comes from the
# register's own history of returns: of the people v19 records as returned who
# had not been listed as prisoners, the share listed as missing in an earlier
# release. For events of 2024-2025 the releases from September 2025 on cover
# the time from capture to return, and about half had been listed as missing:
# the mode. Over all event years the share is lower, because most returns from
# 2022 events came before the first release that can be followed: the floor.
# The ceiling assumes every unrecorded prisoner is among the missing.
#
# Why the model's projection is not used for captivity: v19 is the first
# release to record returns from captivity, and it recorded past returns in
# one batch - nearly all the moves from missing to captivity the model is
# fitted to fall in the window before it (09). Carried forward as a rate, that
# batch projects more future prisoners among today's missing than the
# official figures leave unrecorded at all.
#
# Everyone else among the missing is counted as dead unless the register's own
# record of resolutions gives a reason not to. What it gives is a bound. Of the
# missing of an event year whose fate a later release settles outside
# captivity, some are recorded dead and some leave the register in excess of
# list maintenance (09), the only way it has of recording a person found alive.
# If the unresolved resembled the resolved of their own cohort, that share -
# excess drop-outs over excess drop-outs plus deaths - would be alive. It is an
# UPPER bound on one condition: that the living outside captivity leave the
# list at least as fast as the dead are identified. Then they are
# over-represented among the resolved, since the living resurface faster than
# bodies are recovered and identified, and a drop-out need not be a person
# found alive, so the share overstates the living twice over. A living person
# whom no one reports, who stays listed as missing for good, is what would
# break it: the compiler says the Interior Ministry's database excludes
# unauthorised absence, and a search ends only by a return from captivity, a
# whereabouts established or a body identified, so the condition is likely to
# hold, but nothing measures it. Every event year takes the pooled bound of
# the events of 2023-2025 (composition_bound(): why not each cohort's own, and
# why not 2022's). So in every simulation the share of a cohort's missing,
# prisoners apart, that is alive is drawn uniformly between none and that
# bound, one draw for all cohorts, and the bound itself carries the error of
# resting on five windows (composition_bound()). No source counts these
# people; the bound is what the register can say about them.
#
# SOURCES
#   pow_held           "about 7,000 Ukrainian prisoners of war" held by
#                      Russia - President Zelensky, Munich, 14 Feb 2026
#                      (Ukrinform, 14 Feb 2026). Rounded to the thousand, so
#                      drawn between 6,500 and 7,500.
#   returned_military  7,291 military personnel (235 women, 7,056 men)
#                      returned from captivity since 24 Feb 2022, in the
#                      Coordination Headquarters for the Treatment of Prisoners
#                      of War's report for February 2026
#                      (koordshtab.gov.ua/report, read 19 Sep 2026). The 876
#                      civilians returned by then are left out: the register
#                      lists military personnel only.
pow_held <- c(min = 6500, mode = 7000, max = 7500)
returned_military <- 7291

# composition     the first resolutions outside captivity by event year and
#                 window: year, from_release, dead, out (out: no longer listed
#                 beyond list maintenance), from 09
# register_alive  the register's prisoners plus released prisoners, events
#                 2022-2025, from 08's redistributed counts
# released_prior  where v19's released prisoners had been listed before their
#                 return, by event year: year, prior (prisoner, missing, dead,
#                 not listed), n
alive_evidence <- function(composition, register_alive, released_prior) {
  share_missing <- function(d) {
    d <- d |> summarise(n = sum(n), .by = prior)
    d$n[d$prior == "missing"] / sum(d$n[d$prior %in% c("missing", "not listed")])
  }
  unrecorded <- pow_held + returned_military - register_alive
  out <- tibble(
    s_min = share_missing(released_prior),
    s_mode = share_missing(released_prior |> filter(year >= 2024)),
    s_max = 1,
    held_min = pow_held[["min"]],
    held_mode = pow_held[["mode"]],
    held_max = pow_held[["max"]],
    returned_military = returned_military,
    register_alive = register_alive,
    unrecorded_min = unrecorded[["min"]],
    unrecorded_mode = unrecorded[["mode"]],
    unrecorded_max = unrecorded[["max"]],
    captives_min = s_min * unrecorded_min,
    captives_mode = s_mode * unrecorded_mode,
    captives_max = s_max * unrecorded_max,
    # the bound on the alive outside captivity, all cohorts' resolutions pooled:
    # an alternative in 09; the imputation uses a bound by cohort
    # (composition_bound())
    other_min = 0,
    other_max = sum(composition$out) / sum(composition$out + composition$dead)
  )
  # The central values every deterministic analysis uses: the MEAN of each
  # input's Beta-PERT (shape 4), (min + 4 mode + max) / 6, so the analyses sit
  # where the estimate does - the estimate is the mean of the draws (11), and the
  # military total is linear in these inputs. The share alive outside captivity
  # is uniform between none and its bound, so its mean is half the bound. The
  # two inputs behind the captives are drawn independently, so the captives at
  # the centre are the means' product.
  pert_mean <- function(min, mode, max, shape = 4) (min + shape * mode + max) / (shape + 2)
  out <- out |>
    mutate(
      s_central = pert_mean(s_min, s_mode, s_max),
      held_central = pert_mean(held_min, held_mode, held_max),
      unrecorded_central = held_central + returned_military - register_alive,
      captives_central = s_central * unrecorded_central,
      other_central = (other_min + other_max) / 2
    )
  stopifnot(
    out$s_min <= out$s_mode, out$s_mode <= out$s_max,
    out$unrecorded_min > 0, out$other_max < 1,
    out$captives_min <= out$captives_central, out$captives_central <= out$captives_max,
    out$other_min <= out$other_central, out$other_central <= out$other_max
  )
  out
}
# the prisoners among the missing, for a share s among the missing and a
# count of prisoners held
captives_for <- function(ev, s, held) s * (held + ev$returned_military - ev$register_alive)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the military total in every draw ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The imputation is an accounting identity, by event year c:
#
#   dead among the missing  =  (M_c - C_c) x (1 - a_c)
#
# M_c the register's missing completed for registration lag, C_c the prisoners
# of war among them (the official counts, alive_evidence()), and a_c the share
# of the rest who are alive. Everyone with no evidence of life is counted dead;
# the evidence is the prisoners and the bound on a_c that the composition of the
# register's own resolutions gives (composition_bound()). Nothing is projected:
# the duration model of ual_fit_durations() describes HOW the missing resolve and
# is not part of the estimate, because a missing person it resolves to death and
# one it leaves unresolved are both counted dead.
#
# composition_bound(): the share of a cohort's first resolutions outside
# captivity that leave the register beyond list maintenance, r_c = out / (out +
# dead), summed over the windows between releases. comp: year, from_release,
# dead, out. weights: NULL for the estimate, or a matrix with one row per draw
# and one column per window (in the order of ual_releases), under which the
# windows are summed - resampling the windows, since the clean-ups that produce
# the drop-outs come in some releases and not others. borrow: the event years
# whose own resolutions are not trusted to speak for them; pool: TRUE gives
# every event year the bound of the other years' resolutions pooled under the
# same weights (the estimate), FALSE gives each year its own bound and the
# borrowers the others' pooled. Returns a matrix, draws x event years.
#
# The events of 2022 are not trusted on their own. Of their 318 drop-outs
# beyond list maintenance, 185 fall in one window, v14 to v15, at 18 per 1,000
# person-months of the missing at risk, where no other cohort exceeds 1.1, and
# another 97 in the last, v18 to v19, at 5.4 where no other cohort exceeds 1.8;
# over the five windows they leave at 4.3 per 1,000 person-months against
# 0.3-0.5 for the later cohorts, while their resolutions to death run at about
# half the later cohorts' rate (09l). That is a clean-up of one cohort's
# records, in the months when repatriated remains ran ahead of the register's
# resolutions to death, and not a rate at which the living resurface: a bound
# of 51% read from it would say more about the list than about the people.
#
# The later cohorts are pooled. Of the exits beyond list maintenance behind
# their own bounds, 72% (2023), 91% (2024) and 100% (2025) fall in the last
# window, v18 to v19, so each cohort's own bound - 6.6%, 1.9% and 6.2% - divides
# one clean-up of the list by that cohort's resolutions to death, and the
# differences between the three are as likely the pace at which each cohort's
# dead are identified as a difference in who is alive. One bound for the three,
# 4.0%, is what the register can tell apart; it moves the total by about 300
# against the own bounds. Each cohort's own bound, and 2022's own, are
# alternatives in 09 (Table A19).
bound_borrowers <- 2022
bound_pool <- TRUE
composition_bound <- function(comp, weights = NULL, borrow = bound_borrowers, pool = bound_pool) {
  years <- sort(unique(comp$year))
  wins <- ual_releases$release[-nrow(ual_releases)]
  wide <- function(col) {
    m <- matrix(0, length(wins), length(years), dimnames = list(wins, years))
    m[cbind(match(comp$from_release, wins), match(comp$year, years))] <- comp[[col]]
    m
  }
  if (is.null(weights)) weights <- matrix(1, 1, length(wins))
  stopifnot(ncol(weights) == length(wins), all(borrow %in% years), !all(years %in% borrow))
  out <- weights %*% wide("out")
  all <- out + weights %*% wide("dead")
  lend <- !years %in% borrow
  take <- if (pool) rep(TRUE, length(years)) else !lend
  out[, take] <- rowSums(out[, lend, drop = FALSE])
  all[, take] <- rowSums(all[, lend, drop = FALSE])
  r <- ifelse(all > 0, out / all, 0)
  dimnames(r) <- list(NULL, years)
  r
}

# The identity once, by event year, for 09's tables and the sensitivity steps
# that change the stock of the missing (13d, 13l). stock_missing: month, year,
# missing_stock. captives: the prisoners among the missing, split over event
# years by captive_share (year, share). alive_other: the share of each year's
# missing, prisoners apart, who are alive - one number for every year, or one
# per year in the order of the years.
impute_composition <- function(stock_missing, captives, alive_other, captive_share) {
  stock_missing |>
    summarise(missing_stock = sum(missing_stock), .by = year) |>
    arrange(year) |>
    mutate(captives = captives * captive_share$share[match(year, captive_share$year)],
           unresolved = missing_stock - captives,
           share_alive_other = unname(alive_other),
           alive_other = unresolved * share_alive_other,
           imputed_dead = unresolved - alive_other,
           imputed_alive = captives + alive_other)
}

# impute by the identity, for many draws at once, adding the registered dead
# with their late registrations. What varies between draws:
#   the prisoners        captives: from the share among the missing and the
#                        prisoners held (alive_evidence())
#   the alive otherwise  alive_other: the share a_c of each cohort's missing,
#                        prisoners apart, who are alive
#   registration lag     the completion factors under the draw's weights on the
#                        release windows (lag_factor_draws())
# mil holds 09's inputs (data_inter/ukr_military_inputs.rds): the registered
# dead and missing by event month (month, year, status, registered, row: the
# row of 08b's factor tables), the point factors and what rebuilds them under
# other window weights (lag), the evidence, the composition by event year and
# window, and the split of the prisoners over event years (captive_share).
# alive_other: NULL for the central values, half of each cohort's bound; a
# vector, one share per draw for every cohort alike (the sweeps); or a matrix,
# draws x event years. lag: the completion factors of each draw, rows of 08b's
# factor table x draws (lag_factor_draws()), or NULL for the point factors.
# Returns one row per draw and year: confirmed (registered and late), missing,
# captives, alive_other, imputed_dead and military.
military_draws <- function(mil, captives, alive_other = NULL, lag = NULL) {
  n <- length(captives)
  years <- sort(unique(mil$month$year))
  if (is.null(alive_other)) alive_other <- composition_bound(mil$composition)[rep(1, n), , drop = FALSE] / 2
  if (!is.matrix(alive_other)) alive_other <- matrix(alive_other, n, length(years))
  stopifnot(nrow(alive_other) == n, ncol(alive_other) == length(years),
            all(alive_other >= 0 & alive_other <= 1),
            is.null(lag) || (ncol(lag) == n && nrow(lag) == length(mil$factor_point)))
  dead <- mil$month |> filter(status == "dead")
  miss <- mil$month |> filter(status == "missing")
  factors <- function(rows) {
    if (is.null(lag)) matrix(mil$factor_point[rows], length(rows), n)
    else lag[rows, , drop = FALSE]
  }
  by_year <- function(x, yr) (outer(years, yr, "==") * 1) %*% x
  confirmed <- by_year(dead$registered * factors(dead$row), dead$year)
  missing <- by_year(miss$registered * factors(miss$row), miss$year)
  share <- mil$captive_share$share[match(years, mil$captive_share$year)]
  captives_y <- outer(share, captives)
  unresolved <- missing - captives_y
  stopifnot(all(unresolved >= 0))
  other <- unresolved * t(alive_other)
  imputed_dead <- unresolved - other
  long <- function(x) as.vector(x)
  tibble(
    sim_id = rep(seq_len(n), each = length(years)),
    year = rep(years, n),
    confirmed = long(confirmed),
    missing = long(missing),
    captives = long(captives_y),
    alive_other = long(other),
    imputed_dead = long(imputed_dead),
    military = long(confirmed + imputed_dead)
  )
}

# The register's growth rates by lag band, held non-increasing by lag within
# each status: each band's rate is the minimum over it and the bands before it.
# A lag process cannot speed up with time since the event, yet the missing's
# pooled rates rise again after 24 months, because the bands beyond are
# measured on the 2022-23 cohorts alone, whose missing lists were rebuilt when
# the register imported the Interior Ministry's database from February 2025
# (08b). lambda: a vector, or a matrix with one row per band and one column per
# set of weights; status: the status of each band, the bands ordered by status
# and lag. The dead's rates already decay and are unchanged.
lag_rates_monotone <- function(lambda, status) {
  lambda <- as.matrix(lambda)
  for (s in unique(status)) {
    i <- which(status == s)
    lambda[i, ] <- apply(lambda[i, , drop = FALSE], 2, cummin)
  }
  lambda
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the scenarios of the missing ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No source says how many of the register's missing are dead. The main text
# gives the estimate - the resolution-based scenario - between the BOUNDS of the
# share of the missing counted dead, which assume nothing about the missing
# (scenario_names); two evidence-based scenarios go to the appendix
# (appendix_scenarios). 09 computes the shares and keeps them in mil$scenarios.
#
#   bound_low   none of the missing dead: the named dead alone, completed for
#               listing lag. A logical limit, not a plausible value.
#   central     the RESOLUTION-BASED SCENARIO, the estimate (described below).
#   bound_high  all of the missing dead. A logical limit, not a plausible value.
#
# In the appendix:
#   bodies   the named dead plus the bodies not yet identified, about 17,000
#            (data_input/official_figures.csv, key unidentified_bodies; the
#            source does not say whether they were repatriated or found in
#            Ukraine, nor whether they are combatants). Nine in ten of them
#            are assumed to be Ukrainian combatants listed as missing
#            (bodies_combatant_share), a judgement stated here and nowhere
#            else. Placed across event years as the register now resolves its
#            missing to death by months since the event (09), and applied as
#            each year's share of its missing (p_bodies_year), so that a draw's
#            missing carry it. p_bodies is the share over the four years.
#   rule     no evidence of life means dead. The prisoners of war among the
#            missing from the official counts (alive_evidence()), and a share
#            alive outside captivity drawn under the cap the register's own
#            resolutions give (composition_bound()).
#
# The resolution-based scenario (central): the register's resolved cases
# applied to the unresolved: of the missing whose fate a later release settled,
# the share recorded dead rather than taken prisoner or gone from the list
# beyond list maintenance (p_central, all event years and windows pooled,
# February 2022 left out). It is conservative by construction: the living
# resurface sooner than bodies are identified, so the resolved overstate the
# living, and its captivity share implies more prisoners than the official
# counts admit.
#
# scenario_military(): the combatant total of a scenario in every draw,
# through military_draws(). lag: the completion factors of each draw, or NULL
# for the point factors; n: the number of draws (one, or ncol(lag)). Every
# scenario but the rule fixes its share; the rule takes its evidence at its
# central values (the draws of its inputs are simulation_draws()'s business).
bodies_combatant_share <- 0.9
scenario_names <- c("bound_low", "central", "bound_high")
appendix_scenarios <- c("bodies", "rule")
scenario_military <- function(mil, scenario = c(scenario_names, appendix_scenarios), lag = NULL, n = NULL) {
  scenario <- match.arg(scenario)
  if (is.null(n)) n <- if (is.null(lag)) 1 else ncol(lag)
  sc <- mil$scenarios
  switch(scenario,
    bound_low = military_draws(mil, rep(0, n), 1, lag),
    central = military_draws(mil, rep(0, n), 1 - sc$p_central, lag),
    bound_high = military_draws(mil, rep(0, n), 0, lag),
    # one share per event year, the same in every draw: a row per draw
    bodies = military_draws(mil, rep(0, n),
                            matrix(1 - sc$p_bodies_year[as.character(sort(unique(mil$month$year)))],
                                   n, length(sc$p_bodies_year), byrow = TRUE), lag),
    rule = military_draws(mil, rep(mil$evidence$captives_central, n), NULL, lag)
  )
}

# The registration-lag factors under weights on the release windows. 08b pools
# the register's monthly growth by lag band over the five windows; what the
# register adds differs far more between windows than between event months - a
# release brings a batch or it does not - so the error of the factors is the
# error of resting on five windows, and it is carried the way the bound's is
# (composition_bound()): by weighting the windows. lag: what 08b keeps for this
# (cells: each band's sum of rate x weight and of weight by window; month: the
# lag of each event month and status at the last release; L_ref: the horizon).
# weights: NULL for the estimate, or a matrix, draws x windows. The rates are
# held non-increasing by lag under every set of weights, as 08b holds the
# estimate's (lag_rates_monotone()). Returns the factors, rows of 08b's factor
# table x draws.
lag_factor_draws <- function(lag, weights = NULL) {
  wins <- ual_releases$release[-nrow(ual_releases)]
  if (is.null(weights)) weights <- matrix(1, 1, length(wins))
  stopifnot(ncol(weights) == length(wins))
  bands <- lag$cells |> filter(lo < lag$L_ref) |> distinct(status, lo, hi) |> arrange(status, lo)
  cells <- lag$cells |> filter(lo < lag$L_ref)
  wide <- function(col) {
    m <- matrix(0, nrow(bands), length(wins))
    m[cbind(match(paste(cells$status, cells$lo), paste(bands$status, bands$lo)),
            match(cells$from, wins))] <- cells[[col]]
    m
  }
  num <- wide("num")
  den <- wide("den")
  lambda <- (num %*% t(weights)) / (den %*% t(weights))
  # a band measured only in windows the weights leave out keeps its estimate
  empty <- !is.finite(lambda)
  lambda[empty] <- (rowSums(num) / rowSums(den))[row(lambda)[empty]]
  lambda <- lag_rates_monotone(lambda, bands$status)
  # the months each row of the factor table spends in each band on its way from
  # its lag at the last release to the horizon
  months <- outer(seq_len(nrow(lag$month)), seq_len(nrow(bands)), function(i, j) {
    (lag$month$status[i] == bands$status[j]) *
      pmax(0, pmin(bands$hi[j], lag$L_ref) - pmax(bands$lo[j], lag$month$lag[i]))
  })
  exp(months %*% lambda)
}

# The weights a draw puts on the windows between releases when it sums what
# they record: the resolutions behind the bound (composition_bound()) and the
# register's growth behind the lag factors (lag_factor_draws()), the same
# weights for both. The estimate sums the windows as they
# are. Five windows say little about what a sixth would record - the register's
# clean-ups fall in some releases and not others - and how little is what the
# draws carry: the weights come from a Dirichlet distribution centred on
# `windows` and as concentrated as there are windows, the Bayesian bootstrap
# (Rubin 1981) of a weighted sum. Counts already grow with a window's length,
# so the default weighs the windows alike. n draws, one row each, summing to
# the number of windows so that the counts keep their scale.
window_weight_draws <- function(n, windows = rep(1, length(ual_window_months))) {
  nw <- length(windows)
  g <- matrix(rgamma(n * nw, shape = rep(nw * windows / sum(windows), each = n)), n, nw)
  nw * g / rowSums(g)
}

# The late registrations' share of the confirmed military deaths in every
# draw of 11: the draw's registered-and-late total against the register's own
# count. Late registrations take the registered dead's age-sex profile, so 14
# and 14b split the confirmed deaths, and what they cost, by this share.
late_share_draws <- function() {
  registered <-
    read_rds("data_inter/ukr_ualosses_imputation_table.rds") |>
    transmute(year = as.numeric(year), registered)
  read_rds("data_inter/ukr_sim_param_draws.rds") |>
    filter(role == "confirmed") |>
    transmute(sim_id, year = as.numeric(year), confirmed = draw) |>
    left_join(registered, by = "year") |>
    transmute(sim_id, year, late_share = 1 - registered / confirmed)
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

# Net migration inputs and PERT bounds are built from data in 06 and written
# to data_inter/ (ukr_migration_bounds.rds); see documents/reports/migration_methodology.md.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# net migration: the register reading and the allocation at the mode ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Shared by 06, which builds the input from UNHCR's year-end totals, and 13r,
# which builds it from the bottom-up stock, so that both go from a stock to
# the flows by age and sex the same way. Functions only: each reads its file
# when it is called.
mig_years <- 2022:2025

# Eurostat temporary protection stock at 31 December, EU-27 countries plus
# Iceland, Liechtenstein, Norway and Switzerland, as published. The EU-27
# aggregate row is kept in the table and dropped by the two functions below,
# so nothing is counted twice.
mig_eu_row <- "European Union - 27 countries (from 2020)"
mig_eurostat_extract <- function(yrs = mig_years) {
  cts_rmv <- c("Special value", ":", "Observation flags:", "dp", "d")
  cols_names <-
    expand_grid(date = yrs, sex = c("t", "m", "f", "u"), dat = c("data", "flag")) |>
    mutate(nms = paste(date, sex, dat, sep = "_")) |>
    pull(nms)
  read_xlsx(
    "data_input/refugees_eurostat/migr_asytpsm__custom_21093689_page_spreadsheet.xlsx",
    skip = 10
  ) |>
    set_names(c("ctr", "age", cols_names)) |>
    filter(!ctr %in% cts_rmv) |>
    drop_na(ctr) |>
    select(-ends_with("flag")) |>
    pivot_longer(-c(ctr, age), names_to = "vars", values_to = "mix") |>
    mutate(
      year = as.integer(str_sub(vars, 1, 4)),
      sex = str_sub(vars, 6, 6),
      mix = suppressWarnings(as.double(mix))
    ) |>
    replace_na(list(mix = 0)) |>
    select(ctr, year, sex, age, mix)
}

# year-end totals over the countries
mig_eurostat_total <- function(d) {
  d |>
    filter(ctr != mig_eu_row, sex == "t", age == "Total") |>
    summarise(eurostat = sum(mix), .by = year)
}

# the age-sex structure; unknown sex and age are allocated with the
# distribution of the known cases. France is left out of it: its figures hold
# almost no children (none under 14 and 35 to 115 aged 14-17 among 50,000 to
# 66,000 people in each year), so its adults would enter without the children
# who came with them. Its people stay in the totals and take the structure of
# the other countries. Hungary gives no ages before 2025 and so is in the
# structure of that year only; it is under 1% of the stock.
mig_structure_excluded <- "France"
mig_eurostat_structure <- function(d) {
  stopifnot(all(mig_structure_excluded %in% d$ctr))
  d |>
    filter(ctr != mig_eu_row, !ctr %in% mig_structure_excluded,
           sex %in% c("m", "f"), !age %in% c("Total", "Unknown")) |>
    summarise(mix = sum(mix), .by = c(year, sex, age)) |>
    mutate(cx = mix / sum(mix), .by = year) |>
    mutate(age = case_when(
      age == "Less than 14 years" ~ 0,
      age == "From 14 to 17 years" ~ 14,
      age == "From 18 to 34 years" ~ 18,
      age == "From 35 to 64 years" ~ 35,
      age == "65 years or over" ~ 65
    )) |>
    select(year, sex, age, cx)
}

# the forecast death rates, for the deaths abroad
mig_mx_exp <- function(yrs = mig_years) {
  read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") |>
    filter(source == "frcst", year %in% yrs) |>
    select(year, sex, age, mx)
}

# Ungrouping is done on the stocks, which are positive, and only then
# differenced. The stocks are counts of persons and enter pclm() as they are:
# it chooses its penalty by the BIC of a Poisson likelihood, so the size of the
# counts is what says how far the single ages may depart from a smooth curve
# within a band. Inflating them first (they were once multiplied by 1e5) leaves
# the fit almost unsmoothed and moves men within the band of ages 35 to 64
# towards its end.
mig_ung_stock <- function(d) {
  y <- base::pmax(d$mix, 1)
  fit <- suppressWarnings(pclm(x = d$age, y = y, nlast = 36)$fitted)
  tibble(age = 0:100, S = fit)
}

# Rows are indexed by age at the END of year t, the cohort born in t - x, which
# is how 11 carries its population. The open interval 100+ collects the
# previous year's 99 and 100+.
mig_cohort_diff <- function(d, yrs = mig_years) {
  map_dfr(yrs, function(t) {
    now <- d |> filter(year == t) |> arrange(age) |> pull(S)
    if (t == min(yrs)) {
      return(tibble(year = t, age = 0:100, S_prev = 0, S = now))
    }
    b <- d |> filter(year == t - 1) |> arrange(age) |> pull(S)
    tibble(year = t, age = 0:100, S_prev = c(NA, b[1:99], b[100] + b[101]), S = now)
  })
}

# From a year-end age-sex structure (bands) and a register stock (stock: year,
# register_stock) to the flows by single age: the structure scaled to the
# stock, ungrouped by year and sex, and differenced by cohort,
#   D_t(x) = S_t(x) - S_{t-1}(x-1) + deaths abroad
# Deaths abroad are added back, otherwise a cohort that merely died abroad
# would read as having returned. Age 0 is set to zero: it is mostly children
# born abroad, who were never in Ukraine's population. D is negative where a
# cohort shrank abroad, i.e. net return.
mig_register_flows <- function(bands, stock, mx_exp, yrs = mig_years) {
  bands |>
    left_join(stock |> select(year, register_stock), by = "year") |>
    mutate(mix = cx * register_stock) |>
    group_by(year, sex) |>
    group_modify(~ mig_ung_stock(.x)) |>
    ungroup() |>
    group_by(sex) |>
    group_modify(~ mig_cohort_diff(.x, yrs)) |>
    ungroup() |>
    left_join(mx_exp, by = c("year", "sex", "age")) |>
    mutate(
      deaths_abroad = if_else(year == min(yrs), 0, coalesce(S_prev * mx, 0)),
      D = if_else(age == 0, 0, S - coalesce(S_prev, 0) + deaths_abroad)
    )
}

# Each year's stock is ungrouped on its own, so where a cohort crosses a
# Eurostat band edge between two year ends, a change in the within-band shape
# shows up as a flow. The flows are smoothed lightly across age by penalised
# least squares on second differences, with each year's total by sex and
# Eurostat band (ages 1-13, 14-17, 18-34, 35-64, 65+) held fixed, so every
# migration total and both readings stay as they are. Age 0 stays at zero.
# At this strength no single age moves by more than 300 people (291 with
# UNHCR's totals, 269 with the bottom-up stock of 13r), and few by more than a
# hundred: the flows are already smooth from year to year of age. What remains - the step from zero at age 0 to the return of young
# children at age 1, and the broad waves around 35 and 65, where the older
# women's returns in the 65+ band are spread by the ungrouping into the ages
# around it - is set by the band totals themselves, which no smoothing that
# keeps them can remove.
mig_lambda_D <- 5
mig_smooth_flows <- function(D, age, lambda = mig_lambda_D) {
  keep <- age >= 1
  x <- D[keep]
  n <- length(x)
  band <- cut(age[keep], c(1, 14, 18, 35, 65, Inf), right = FALSE)
  C <- t(model.matrix(~ band - 1))
  P <- diff(diag(n), differences = 2)
  A <- rbind(cbind(diag(n) + lambda * crossprod(P), t(C)),
             cbind(C, matrix(0, nrow(C), nrow(C))))
  out <- D
  out[keep] <- solve(A, c(x, C %*% x))[seq_len(n)]
  out
}

# the smoothing applied to a table of flows by year and sex, with the flows
# before smoothing kept as D_raw
mig_smoothed <- function(flows) {
  out <-
    flows |>
    arrange(year, sex, age) |>
    mutate(D_raw = D, D = mig_smooth_flows(D, age), .by = c(year, sex))
  stopifnot(isTRUE(all.equal(
    out |> summarise(d = sum(D), .by = c(year, sex)) |> pull(d),
    out |> summarise(d = sum(D_raw), .by = c(year, sex)) |> pull(d)
  )))
  out
}

# each year's register reading and its parts
mig_flow_totals <- function(flows) {
  flows |>
    summarise(
      register = sum(D),
      departures = sum(base::pmax(D, 0)),
      returns = sum(base::pmin(D, 0)),
      deaths_abroad = sum(deaths_abroad),
      births_abroad_dropped = sum(S[age == 0]),
      .by = year
    )
}

# The crossings reading: CES (2026) table 1 net crossings of citizens, plus the
# 277k who left through Russia or Belarus and are now elsewhere (they are in
# the registers but crossed no western border).
mig_crossings <- function() {
  ces <- read_csv("data_input/migration/ces_2026_figures.csv", show_col_types = FALSE)
  ces |>
    filter(key == "net_outflow_citizens") |>
    select(year, ces_net = value) |>
    left_join(
      ces |>
        filter(key == "left_via_russia_belarus_now_elsewhere") |>
        select(year, via_ru_by = value),
      by = "year"
    ) |>
    mutate(year = as.integer(year), via_ru_by = coalesce(via_ru_by, 0),
           crossings = ces_net + via_ru_by)
}

# The allocation at the mode. The mode total (mode_tot: year, mode_tot)
# differs from the register reading, and the difference is spread over the
# DEPARTURES profile (positive part of D), leaving the return cells as the
# register has them:  ems(x) = D(x) + (mode - sum D) * w(x).  11 applies the
# same additive rule to every draw. Returns year, sex, age, mix (negative =
# leaving, the convention 11 and 13 read) and w.
mig_allocate <- function(flows, mode_tot) {
  out <-
    flows |>
    mutate(w = base::pmax(D, 0) / sum(base::pmax(D, 0)), .by = year) |>
    left_join(flows |> summarise(register = sum(D), .by = year), by = "year") |>
    left_join(mode_tot, by = "year") |>
    mutate(
      ems = D + (mode_tot - register) * w,
      mix = -ems
    ) |>
    select(year, sex, age, mix, w) |>
    arrange(year, sex, age)
  stopifnot(
    isTRUE(all.equal(
      out |> summarise(v = -sum(mix), .by = year) |> arrange(year) |> pull(v),
      mode_tot |> arrange(year) |> pull(mode_tot),
      tolerance = 1e-8
    )),
    all(abs(out |> summarise(s = sum(w), .by = year) |> pull(s) - 1) < 1e-9)
  )
  out
}

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
# outputs/manifest.csv records what every table and figure of the paper was built from.
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
#                  (conf_cmb = the register's dead that year, with their late
#                  registrations), and optionally dk_f and dk_m, the draw's
#                  deviation of each sex's Lee-Carter index from its mean
#                  path
#   static_inputs  year-sex-age grid with ems, w, mx, fx, prop_cvs,
#                  prop_cmb_dead, prop_cmb_miss, and for dk the index loading
#                  bx and the forecast variance vk (04). A draw's rates are
#                  mx x exp(bx dk - bx^2 vk / 2), whose mean over draws is mx
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

  has_dk <- all(c("dk_f", "dk_m") %in% names(draws_this_sim)) &&
    all(c("bx", "vk") %in% names(static_inputs))

  # the confirmed deaths are part of the drawn combatant total in every draw and
  # every sensitivity analysis (the imputed dead are never negative), so the
  # clamp on their ratio below never binds; this makes sure of it
  stopifnot(all(draws_this_sim$draw_cmb >= draws_this_sim$conf_cmb * (1 - 1e-9)))

  df_sim <- static_inputs %>%
    left_join(draws_this_sim, by = "year") %>%
    left_join(mig_base_by_year, by = "year") %>%
    mutate(
      # the draw's counterfactual (see ARGUMENTS)
      mx = if (has_dk) mx * exp(bx * if_else(sex == "f", dk_f, dk_m) - bx^2 * vk / 2) else mx,
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
# civilian deaths from uncertainUCDP ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UCDP's low and high values come from its coding rules about sources, not from
# an estimate of measurement error, and its best value undercounts, most of all
# for small events (Vesco, Randahl, Hegre, Hogbladh and Yilmaz 2026, Journal of
# Conflict Resolution 70(7-8)). From a survey of UCDP's coders they fit, for any
# reported count and type of violence, a distribution of plausible deaths,
# released as the R package uncertainUCDP. Under it, each year's civilian deaths
# are every event's deaths given UCDP's best count and its type of violence,
# times the event's civilian share as 07 computes it, summed by year. 10 draws
# them (uncertain_ucdp_n sets of the four years) and takes their exact mean.
#
# The estimate draws the civilians from a Beta-PERT on UCDP's low, best and high
# values ("pert"); uncertainUCDP is the alternative, which 13t reports beside it
# (Table 6, Table A14). With "uncertainUCDP" here the two change places: 11 then
# samples 10's draws through the same uniform numbers the Beta-PERT takes, so
# every other input keeps its draws.
civilian_source <- "pert"
uncertain_ucdp_n <- 4000

# Each event's civilian deaths at UCDP's best, as 07 counts them: the deaths of
# civilians plus the event's deaths of unknown side times the civilians' share of
# the attributed deaths in its dyad-year, with 07's fallback to the dyad's
# all-year shares, Crimea and Sevastopol left out. Returns, for 2022-2025, each
# event's year, best count t, type of violence tov, civilian share s and
# civilian deaths civ_best; their sums must reproduce 07's civilian totals.
ucdp_civilian_events <- function() {
  ev <-
    read_rds("data_inter/ucdp_ged_events_slim.rds") |>
    filter(country == "Ukraine") |>
    mutate(u = coalesce(u, 0))
  in_crimea <- str_detect(coalesce(ev$adm1, ""), regex("crimea|sevastopol", ignore_case = TRUE))
  ev$t[in_crimea] <- ev$t[in_crimea] - ev$c[in_crimea]
  ev$c[in_crimea] <- 0
  ev <- ev |>
    mutate(inter = type_of_violence == 1 & str_detect(side_b, "Government"),
           grp = if_else(inter, paste(side_a, side_b, sep = " | "), "intra"))
  gy <- ev |> summarise(a = sum(a), b = sum(b), c = sum(c), u = sum(u), .by = c(grp, year))
  ga <- gy |> summarise(pi_all = sum(c) / sum(a + b + c), .by = grp)
  gy <- gy |>
    left_join(ga, by = "grp") |>
    mutate(pi_c = case_when(u == 0 ~ 0,
                            a != 0 & b != 0 & c != 0 ~ c / (a + b + c),
                            .default = pi_all)) |>
    select(grp, year, pi_c)
  ev <- ev |>
    left_join(gy, by = c("grp", "year")) |>
    mutate(civ_best = c + u * pi_c,
           s = if_else(t > 0, civ_best / t, 0),
           tov = c("sb", "ns", "os")[type_of_violence]) |>
    filter(year %in% 2022:2025) |>
    select(year, t, tov, s, civ_best, type_of_violence, a, b, c, u)
  stopifnot(all(ev$s >= 0 & ev$s <= 1 + 1e-9), !anyNA(ev$s))
  check07 <- read_rds("data_inter/ukr_ucdp_invals.rds") |>
    filter(role == "civilians", year %in% 2022:2025) |>
    summarise(dts = sum(dts), .by = year) |>
    left_join(ev |> summarise(ours = sum(civ_best), .by = year), by = "year")
  stopifnot(isTRUE(all.equal(check07$dts, check07$ours, tolerance = 1e-6)))
  ev
}

# The same for Ukraine's combatants, for the scenario that takes every death from
# UCDP (13t): each event's deaths of Ukraine's side as 07 counts them - in the
# state dyad the deaths of the side that is the Government of Ukraine, wherever
# the event was, and in the other events located in Ukraine both sides - with the
# event's deaths of unknown side shared out as 07 does, within the dyad, country
# and year. Returns the columns uncertain_ucdp_draws() and uncertain_ucdp_mean()
# take: year, best count t, type of violence tov, the combatants' share s and
# their deaths cmb_best, whose sums must reproduce 07's combatant totals.
ucdp_combatant_events <- function() {
  ev <- read_rds("data_inter/ucdp_ged_events_slim.rds") |> mutate(u = coalesce(u, 0))
  in_crimea <- ev$country == "Ukraine" &
    str_detect(coalesce(ev$adm1, ""), regex("crimea|sevastopol", ignore_case = TRUE))
  ev$t[in_crimea] <- ev$t[in_crimea] - ev$c[in_crimea]
  ev$c[in_crimea] <- 0
  ev <- ev |>
    mutate(inter = type_of_violence == 1 & str_detect(side_b, "Government"),
           grp = if_else(inter, paste(country, side_a, side_b, sep = " | "), paste("intra", country)))
  gy <- ev |> summarise(a = sum(a), b = sum(b), c = sum(c), u = sum(u), .by = c(grp, year))
  ga <- gy |> summarise(pa = sum(a) / sum(a + b + c), pb = sum(b) / sum(a + b + c), .by = grp)
  gy <- gy |>
    left_join(ga, by = "grp") |>
    mutate(full = a != 0 & b != 0 & c != 0,
           pi_a = case_when(u == 0 ~ 0, full ~ a / (a + b + c), .default = pa),
           pi_b = case_when(u == 0 ~ 0, full ~ b / (a + b + c), .default = pb)) |>
    select(grp, year, pi_a, pi_b)
  ev <- ev |>
    left_join(gy, by = c("grp", "year")) |>
    mutate(a2 = a + u * pi_a, b2 = b + u * pi_b,
           cmb_best = case_when(inter & side_a == "Government of Ukraine" ~ a2,
                                inter & side_b == "Government of Ukraine" ~ b2,
                                !inter & country == "Ukraine" ~ a2 + b2,
                                .default = 0),
           s = if_else(t > 0, cmb_best / t, 0),
           tov = c("sb", "ns", "os")[type_of_violence]) |>
    filter(year %in% 2022:2025, cmb_best > 0) |>
    select(year, t, tov, s, cmb_best, type_of_violence)
  stopifnot(all(ev$s >= 0 & ev$s <= 1 + 1e-9), !anyNA(ev$s))
  check07 <- read_rds("data_inter/ukr_ucdp_invals.rds") |>
    filter(role == "combatants", year %in% 2022:2025) |>
    summarise(dts = sum(dts), .by = year) |>
    left_join(ev |> summarise(ours = sum(cmb_best), .by = year), by = "year")
  stopifnot(isTRUE(all.equal(check07$dts, check07$ours, tolerance = 1e-6)))
  ev
}

# n_draw sets of the four years' civilian deaths: a matrix, years by draws
uncertain_ucdp_draws <- function(ev, n_draw, seed = 20261018) {
  civ_ev <- ev |> filter(s > 0)
  years <- 2022:2025
  set.seed(seed)
  draws <- matrix(0, length(years), n_draw, dimnames = list(years, NULL))
  yr_idx <- match(civ_ev$year, years)
  for (k in seq_len(n_draw)) {
    d <- numeric(nrow(civ_ev))
    for (tv in unique(civ_ev$tov)) {
      i <- which(civ_ev$tov == tv)
      d[i] <- uncertainUCDP::runcertainUCDP(length(i), fatalities = civ_ev$t[i], tov = tv)
    }
    draws[, k] <- rowsum(d * civ_ev$s, yr_idx, reorder = TRUE)[, 1]
  }
  stopifnot(all(is.finite(draws)), all(draws > 0))
  draws
}

# the exact mean of each year's civilian deaths: each event's mean deaths given
# its best count, times its civilian share, summed
uncertain_ucdp_mean <- function(ev) {
  m <- numeric(nrow(ev))
  for (tv in unique(ev$tov)) {
    i <- which(ev$tov == tv & ev$t > 0)
    m[i] <- uncertainUCDP::mean_uncertainUCDP(ev$t[i], tov = tv)
  }
  ev |> mutate(m = m * s) |> summarise(mean = sum(m), .by = year) |> arrange(year)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the simulation's draws ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Every uncertain input of one set of n simulations, drawn in a fixed order
# from one seed, so that 11 and the PERT-shape check in 13g use the same
# random numbers. The disputed inputs are Beta-PERT quantiles of uniforms,
# with the given shape (4 is the standard PERT); the estimated ones are drawn
# from their sampling distributions:
#   civilians        one uniform per year and simulation: each year's count
#                    comes from that year's own event reporting
#   the missing      the share of the unrecorded prisoners among the missing,
#   alive            the prisoners held and the share of the unresolved alive
#                    for other reasons, once per simulation each: they are
#                    properties of the missing, not of a year
#   the release      the weights on the windows between releases, once per
#   windows          simulation from a Dirichlet (window_weight_draws()): the
#                    bound on the alive and the registration-lag factors are
#                    both rebuilt under them, so a draw that leans on a window
#                    leans on it for both
#   counterfactual  the Lee-Carter index of each sex, a random walk with
#                    drift: the drift's error and each year's innovation,
#                    correlated between the sexes as 04 measured
#   migration        the western blend weight once per simulation, uniform
#                    between the two readings and applied to every year's
#                    two; Russia once, from its Beta-PERT
# The military total follows from the three inputs on the missing alive and
# the window weights (military_draws()). The two groups are also run alone, the
# other at its central values or point estimates, so that 15 can attribute the
# spread of the results.
# Returns the draws by simulation and year (draws_df), the same by role
# (param_draws) and the inputs on the missing by simulation (alive_draws).
simulation_draws <- function(param_table, mil, lc_error, n, shape = 4, seed = 42) {
  set.seed(seed)
  q <- function(u, a, m, b) qpert(u, min = a, mode = m, max = b, shape = shape)
  ev <- mil$evidence

  cvs <- param_table |> filter(role == "civilians") |> arrange(year)
  # uncertainUCDP's draws (10) through the uniform the Beta-PERT would take
  uu <- if (civilian_source == "uncertainUCDP") read_rds("data_inter/ukr_civilians_uncertain_ucdp.rds")$draws
  draws_cvs_long <- map_dfr(seq_len(nrow(cvs)), function(i) {
    tibble(role = "civilians", year = cvs$year[i], sim_id = seq_len(n), u = runif(n)) |>
      mutate(draw = if (is.null(uu)) q(u, cvs$min[i], cvs$mode[i], cvs$max[i])
                    else quantile(uu[as.character(cvs$year[i]), ], u, names = FALSE, type = 7)) |>
      select(-u)
  })

  # The missing alive. The prisoners among them follow from two drawn inputs,
  # the share of the unrecorded prisoners who are among the missing and the
  # prisoners held. The alive outside captivity are, in each cohort, a share of
  # the rest between none and the cohort's bound: ONE uniform draw for all
  # cohorts, since what is unknown - how much faster the living resurface than
  # the dead are identified - is one thing and not four, and the bound itself
  # from the windows between releases resampled (composition_bound()).
  alive_draws <-
    tibble(
      sim_id = seq_len(n),
      s = q(runif(n), ev$s_min, ev$s_mode, ev$s_max),
      held = q(runif(n), ev$held_min, ev$held_mode, ev$held_max),
      u_other = runif(n)
    ) |>
    mutate(captives = captives_for(ev, s, held))
  # The release windows: one set of weights per simulation, under which both
  # the bound and the registration-lag factors are rebuilt. A release that adds
  # many records tends to drop many as well - the last brought both, the one
  # before it little of either - and more growth raises the total while more
  # exits lower it, so the two move in opposite directions under the same
  # weights, and drawing them apart would add a spread no set of releases could
  # produce.
  window <- window_weight_draws(n)
  bound <- composition_bound(mil$composition, window)
  bound_point <- composition_bound(mil$composition)[rep(1, n), , drop = FALSE]
  lag <- lag_factor_draws(mil$lag, window)
  a_draw <- alive_draws$u_other * bound

  # THE ESTIMATE is the central scenario (scenario_military()): the register's
  # resolved cases applied to the unresolved, the share fixed, so that between
  # draws only the registered dead and the missing move, through the lag
  # factors under the draw's window weights. The bounds need no draw of their
  # own: in every draw they are the confirmed dead, and the confirmed plus all
  # the missing (roles confirmed and missing). The bodies scenario is drawn the
  # same way as the estimate (mil_bodies). The RULE scenario draws its own
  # evidence: the prisoners and where between none and the cap the rest sit, and
  # the cap itself under the window weights (mil_rule). Then, for the shares of
  # the variance (15), the rule's evidence alone (mil_alive), and the windows
  # through the cap alone (mil_model), through the lag alone (mil_lag), and
  # through both (mil_window), each with the rest at its centre.
  mil_all <- scenario_military(mil, "central", lag, n)
  mil_bodies <- scenario_military(mil, "bodies", lag, n)
  mil_upper <- military_draws(mil, alive_draws$captives, a_draw, lag)
  captives_c <- rep(ev$captives_central, n)
  mil_alive <- military_draws(mil, alive_draws$captives, alive_draws$u_other * bound_point)
  mil_window <- military_draws(mil, captives_c, bound / 2, lag)
  mil_model <- military_draws(mil, captives_c, bound / 2)
  mil_lag <- military_draws(mil, captives_c, lag = lag)
  alive_draws <-
    alive_draws |>
    left_join(mil_upper |> summarise(across(c(missing, alive_other, imputed_dead), sum),
                                     .by = sim_id) |>
                rename(alive_other_n = alive_other),
              by = "sim_id")

  # the counterfactual: the deviation of each sex's forecast index from its
  # mean path, h = years after the last fitted year
  walk <- lc_error$walk
  h_max <- 2025 - max(walk$last_year)
  rho <- lc_error$rho
  pair <- function(k) {
    z1 <- matrix(rnorm(n * k), n)
    z2 <- matrix(rnorm(n * k), n)
    list(f = z1, m = rho * z1 + sqrt(1 - rho^2) * z2)
  }
  e_drift <- pair(1)
  e_step <- pair(h_max)
  dk_long <- map_dfr(c("f", "m"), function(sx) {
    wk <- walk |> filter(sex == sx)
    drift_err <- e_drift[[sx]][, 1] * wk$sigma / sqrt(wk$n - 1)
    path <- t(apply(e_step[[sx]] * wk$sigma, 1, cumsum)) + outer(drift_err, seq_len(h_max))
    expand_grid(sim_id = seq_len(n), year = 2022:2025) |>
      mutate(role = paste0("lc_", sx), draw = path[cbind(sim_id, year - wk$last_year)])
  })

  west <- param_table |> filter(role == "mig_west")
  # uniform between the two readings: nothing says the truth lies nearer their
  # midpoint than either end, so no weight is favoured. The same uniform number
  # as before, so every other input keeps its draws.
  w <- runif(n)
  draws_west_long <-
    expand_grid(sim_id = seq_len(n), west |> select(year, crossings, register)) |>
    mutate(role = "mig_west", draw = crossings + w[sim_id] * (register - crossings)) |>
    select(role, year, sim_id, draw)
  ru <- param_table |> filter(role == "mig_ru")
  u_ru <- runif(n)
  draws_ru_long <-
    expand_grid(sim_id = seq_len(n), ru |> select(year, min, mode, max)) |>
    mutate(role = "mig_ru", draw = q(u_ru[sim_id], min, mode, max)) |>
    select(role, year, sim_id, draw)

  as_role <- function(d, col, role) d |> transmute(role = role, year, sim_id, draw = .data[[col]])
  param_draws <- bind_rows(
    draws_cvs_long,
    as_role(mil_all, "military", "combatants"),
    as_role(mil_all, "confirmed", "confirmed"),
    as_role(mil_all, "missing", "missing"),
    as_role(mil_bodies, "military", "mil_bodies"),
    as_role(mil_upper, "military", "mil_rule"),
    as_role(mil_alive, "military", "mil_alive"),
    as_role(mil_window, "military", "mil_window"),
    as_role(mil_model, "military", "mil_model"),
    as_role(mil_lag, "military", "mil_lag"),
    dk_long,
    draws_west_long,
    draws_ru_long
  )
  draws_df <-
    param_draws |>
    filter(role %in% c("civilians", "combatants", "confirmed", "mig_west", "mig_ru", "lc_f", "lc_m")) |>
    mutate(role = case_when(role == "civilians" ~ "draw_cvs", role == "combatants" ~ "draw_cmb",
                            role == "confirmed" ~ "conf_cmb", role == "lc_f" ~ "dk_f",
                            role == "lc_m" ~ "dk_m", .default = "draw_mig")) |>
    summarise(draw = sum(draw), .by = c(sim_id, year, role)) |>
    pivot_wider(names_from = role, values_from = draw) |>
    arrange(sim_id, year)
  list(draws_df = draws_df, param_draws = param_draws, alive_draws = alive_draws)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the deterministic projection at the mode, for the sensitivity steps ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Every input at its central value, set up as 13 sets it up: the conflict and
# migration totals from the parameter table, the counterfactual, fertility
# carried forward from 2023, and the age-sex profiles of civilian, registered
# and imputed deaths. A sensitivity step changes one of these and projects
# again. "At the mode" is shorthand, here and in the sensitivity steps, for
# these central values: every drawn input at the MEAN of its distribution, so
# that the central projection stands where the mean of the draws does (11) -
# the three inputs on how many of the missing are alive (alive_evidence()),
# which set the military total (the "mode" column of the parameter table's
# combatant rows, 10), and the civilians and Russia (central_values()).
#
# central_values(): a parameter table, or the migration bounds, with the "mode"
# of the two inputs drawn from a Beta-PERT of its own - civilians and Russia -
# replaced by that distribution's mean, (min + 4 mode + max) / 6.
# The civilian range is longer above its mode than below, so the mode sits
# below the mean of the draws. The western net outflow is a blend of two
# readings whose mode is their midpoint and its mean already.
central_values <- function(p) {
  # the civilians are a Beta-PERT only when civilian_source is "pert"; uncertainUCDP
  # gives 10 their exact mean as the mode already
  pert_roles <- c(if (civilian_source == "pert") "civilians", "mig_ru")
  p |> mutate(mode = if_else(role %in% pert_roles, (min + 4 * mode + max) / 6, mode))
}
mode_projection_inputs <- function() {
  # returned as it is: 14c and 14d draw the simulation's inputs from it again,
  # and those take the modes. Only the central draws below take the means.
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
    central_values() |>
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
