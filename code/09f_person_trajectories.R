# ==============================================================================
# STEP 09F - PER-PERSON TRAJECTORIES ACROSS ALL REGISTER VERSIONS
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Fixes a real independence violation in 09d/the Phase 1 pairwise-window
# design: a person missing in v14, still missing in v16, then dead in v18
# currently contributes TWO separate "person-window" rows - (v14->v16:
# missing->missing) and (v16->v18: missing->dead) - each treated as an
# independent Surv(0, time) observation. That's the same person's one
# continuous history, artificially cut into two renewal observations that
# both start the clock at 0 instead of respecting that the second interval
# begins where the first one left off.
#
# This script instead builds ONE row per person per resolved interval, using
# proper counting-process (Tstart, Tstop] structure with REAL elapsed time
# since disappearance (date_evnt), which coxph()/mstate already support
# natively (Surv(Tstart, Tstop, status)) - no new package needed.
#
# ENTRY-POINT DESIGN (avoids double-counting a person missing in multiple
# consecutive registers): a person is only added to the panel at their
# FIRST appearance as "missing" (v14 first, else v16, else v18 - v19 is
# endpoint-only, nothing comes after it to track forward). From their entry
# register onward, we look up their status in every later register via the
# same fuzzy-matching logic already validated in 09.
#
# INPUT   raw registers v14, v16, v18, v19 (same reads as 09)
# OUTPUT  data_inter/ukr_ualosses_person_trajectories.rds
#         data_inter/ukr_ualosses_intervals_countingprocess.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

if (!require("stringdist", quietly = TRUE)) {
  install.packages("stringdist")
  library(stringdist)
}

# ==============================================================================
# 1. READ REGISTERS (same as 09)
# ==============================================================================

read_reg <- function(path) {
  read_xlsx(require_raw(path), sheet = "Database") |>
    mutate(
      name_last = as.character(LastName),
      name_first = as.character(FirstName),
      name_patr = as.character(Patronym),
      name_full = paste(name_last, name_first, name_patr, sep = " "),
      date_bth = ymd("1900-01-01") + as.numeric(DateBirth),
      date_evnt = ymd("1900-01-01") + as.numeric(DateEvent),
      year = year(date_evnt)
    ) |>
    filter(year %in% 2022:2025, Nationality == "Ukraine") |>
    select(name_last, name_first, name_patr, name_full, date_bth, year,
           date_evnt, status = Status)
}

message("Reading registers...")
v14 <- read_reg("data_input/ualosses_hubert_datasets/250916_UKR_ualosses_Personnel_v14.xlsx")
v16 <- read_reg("data_input/ualosses_hubert_datasets/251204_UKR_ualosses_Personnel_v16.xlsx")
v18 <- read_reg("data_input/ualosses_hubert_datasets/260530_UKR_ualosses_Personnel_v18.xlsx")
v19 <- read_reg("data_input/ualosses_hubert_datasets/260919_UKR_ualosses_Personnel_v19.xlsx")

register_dates <- c(
  v14 = as.Date("2025-09-16"), v16 = as.Date("2025-12-04"),
  v18 = as.Date("2026-05-30"), v19 = as.Date("2026-09-19")
)

# ==============================================================================
# 2. LOOKUP FUNCTION: find a person's status in a target register
# ==============================================================================
# Same exact-then-fuzzy logic validated in 09 (compound name+DOB key primary,
# vectorized data.table non-equi join + stringdist for the residual), reused
# here as a lookup (source person -> status in target) rather than a full
# transition extraction, since we need to look each anchor person up in
# potentially 3 different later registers independently.

lookup_status <- function(anchor, target, dob_window = 30, threshold = 0.90) {
  anchor <- anchor %>% mutate(.row_id = row_number())

  target_dedup <- target %>%
    distinct(name_full, date_bth, .keep_all = TRUE) %>%
    select(name_full, date_bth, status_found = status)

  matched <- anchor %>%
    left_join(target_dedup, by = c("name_full", "date_bth"), relationship = "many-to-one")

  unmatched <- matched %>% filter(is.na(status_found))
  if (nrow(unmatched) > 0) {
    src_dt <- as.data.table(unmatched %>% select(.row_id, name_full, date_bth))
    tgt_dt <- as.data.table(target %>% select(name_full_target = name_full,
                                              status_found = status, date_bth))
    src_dt[, `:=`(dob_lo = date_bth - dob_window, dob_hi = date_bth + dob_window)]

    pairs <- tgt_dt[src_dt,
                    on = .(date_bth >= dob_lo, date_bth <= dob_hi),
                    allow.cartesian = TRUE,
                    .(.row_id, name_full, name_full_target, status_found)]
    pairs <- pairs[!is.na(name_full_target)]

    if (nrow(pairs) > 0) {
      pairs[, name_dist := stringdist(name_full, name_full_target, method = "lv")]
      pairs[, name_sim := 1 - (name_dist / pmax(nchar(name_full), nchar(name_full_target)))]
      pairs <- pairs[name_sim >= threshold]
      setorder(pairs, .row_id, -name_sim, name_dist)
      best <- pairs[!duplicated(.row_id)]

      matched <- matched %>%
        left_join(best %>% select(.row_id, status_found_fuzzy = status_found), by = ".row_id") %>%
        mutate(status_found = coalesce(status_found, status_found_fuzzy)) %>%
        select(-status_found_fuzzy)
    }
  }

  # NOT folding released_prisoner -> alive here; that normalization is
  # applied once, later, after all lookups are combined (section 4) so it
  # stays in exactly one place regardless of which register a status came
  # from - same principle as 09's normalize_status2(), applied at the right
  # layer for a multi-register panel instead of a single pairwise window.
  matched %>% select(-.row_id) %>% pull(status_found)
}

# ==============================================================================
# 3. BUILD ENTRY COHORTS: each person enters the panel ONCE, at their FIRST
#    appearance as "missing", to avoid double-counting someone missing in
#    consecutive registers as two different people.
# ==============================================================================

message("\nBuilding entry cohorts (avoiding double-counting across registers)...")

entry_v14 <- v14 %>% filter(status == "missing") %>%
  mutate(entry_register = "v14", entry_date = register_dates[["v14"]])

# v16 entrants: missing in v16, but NOT already captured via v14 (i.e., not
# found in v14 at all - genuinely new to the panel, not just "still missing
# since v14" which entry_v14 already covers).
v16_missing <- v16 %>% filter(status == "missing")
already_in_v14 <- lookup_status(v16_missing, v14)
entry_v16 <- v16_missing %>%
  mutate(.found_in_v14 = already_in_v14) %>%
  filter(is.na(.found_in_v14)) %>%
  select(-.found_in_v14) %>%
  mutate(entry_register = "v16", entry_date = register_dates[["v16"]])

message(sprintf("  v14 entrants: %d | v16 new entrants: %d (of %d missing-in-v16)",
                nrow(entry_v14), nrow(entry_v16), nrow(v16_missing)))

# v18 entrants: missing in v18, not found in EITHER v14 or v16.
v18_missing <- v18 %>% filter(status == "missing")
found_v14 <- lookup_status(v18_missing, v14)
found_v16 <- lookup_status(v18_missing, v16)
entry_v18 <- v18_missing %>%
  mutate(.f14 = found_v14, .f16 = found_v16) %>%
  filter(is.na(.f14) & is.na(.f16)) %>%
  select(-.f14, -.f16) %>%
  mutate(entry_register = "v18", entry_date = register_dates[["v18"]])

message(sprintf("  v18 new entrants: %d (of %d missing-in-v18)",
                nrow(entry_v18), nrow(v18_missing)))

# v19 is endpoint-only: nothing comes after it, so a person newly missing at
# v19 contributes no resolved OR censored forward interval - skip.

entrants <- bind_rows(entry_v14, entry_v16, entry_v18) %>%
  mutate(person_id = row_number())

message(sprintf("\nTotal panel entrants: %d", nrow(entrants)))

# ==============================================================================
# 4. FOR EACH ENTRANT, LOOK UP STATUS AT EVERY LATER REGISTER
# ==============================================================================

message("\nLooking up each entrant's status at every later register...")

normalize_status <- function(x) ifelse(x == "released_prisoner", "alive", x)

traj <- entrants %>% select(person_id, name_full, date_bth, year, date_evnt,
                            entry_register, entry_date)

if (nrow(entry_v14) > 0 || nrow(entry_v16) > 0) {
  # v14 and v16 entrants both need a v16 reading if they entered at v14
  # (v16 entrants ARE v16's reading, trivially "missing" - no lookup needed).
  need_v16_lookup <- traj %>% filter(entry_register == "v14")
  status_at_v16 <- lookup_status(need_v16_lookup, v16)
  traj <- traj %>%
    left_join(
      bind_rows(
        need_v16_lookup %>% mutate(status_v16 = normalize_status(status_at_v16)) %>%
          select(person_id, status_v16),
        entrants %>% filter(entry_register == "v16") %>%
          mutate(status_v16 = "missing") %>% select(person_id, status_v16)
      ),
      by = "person_id"
    )
} else {
  traj$status_v16 <- NA_character_
}

need_v18_lookup <- traj %>% filter(entry_register %in% c("v14", "v16"))
status_at_v18 <- lookup_status(need_v18_lookup, v18)
traj <- traj %>%
  left_join(
    bind_rows(
      need_v18_lookup %>% mutate(status_v18 = normalize_status(status_at_v18)) %>%
        select(person_id, status_v18),
      entrants %>% filter(entry_register == "v18") %>%
        mutate(status_v18 = "missing") %>% select(person_id, status_v18)
    ),
    by = "person_id"
  )

need_v19_lookup <- traj  # everyone needs a v19 reading, it's the final register
status_at_v19 <- lookup_status(need_v19_lookup, v19)
traj$status_v19 <- normalize_status(status_at_v19)
traj$status_v19 <- ifelse(is.na(traj$status_v19), "alive", traj$status_v19)  # not found anywhere later = resolved alive

message("✓ Trajectories built")
print(traj %>% count(entry_register))

write_rds(traj, "data_inter/ukr_ualosses_person_trajectories.rds")

# ==============================================================================
# 5. BUILD COUNTING-PROCESS (Tstart, Tstop, status] INTERVALS
# ==============================================================================
# One or more rows per person, each interval running from their last
# CONFIRMED "still missing" observation to their next observation, with
# Tstart/Tstop as REAL elapsed days since date_evnt (not window-relative).
# This is what fixes the independence violation: a person's second interval
# starts where their first one ended, not back at 0.

message("\nBuilding counting-process intervals...")

build_intervals <- function(row) {
  # Sequence of (register, date, status) for this person, starting at their
  # entry (always "missing" by construction) and including every later
  # register - v16/v18/v19 as applicable, skipping any NA (not yet resolved
  # observation gap that isn't actually informative, e.g. a v14 entrant has
  # no "v14" column to re-observe, entry itself IS the first observation).
  obs <- list()
  obs[[1]] <- list(date = row$entry_date, status = "missing")

  add_if_present <- function(reg_name, status_val) {
    if (!is.na(status_val)) {
      obs[[length(obs) + 1]] <<- list(date = register_dates[[reg_name]], status = status_val)
    }
  }
  if (row$entry_register == "v14") add_if_present("v16", row$status_v16)
  if (row$entry_register %in% c("v14", "v16")) add_if_present("v18", row$status_v18)
  add_if_present("v19", row$status_v19)

  obs_df <- bind_rows(lapply(obs, as_tibble))

  # Build intervals between consecutive observations
  n <- nrow(obs_df)
  if (n < 2) return(NULL)

  purrr::map_dfr(1:(n - 1), function(i) {
    tibble(
      person_id = row$person_id,
      year = row$year,
      Tstart = as.numeric(obs_df$date[i] - row$date_evnt),
      Tstop  = as.numeric(obs_df$date[i + 1] - row$date_evnt),
      status_from = obs_df$status[i],
      status_to = obs_df$status[i + 1]
    )
  })
}

intervals <- traj %>%
  split(seq_len(nrow(.))) %>%
  purrr::map_dfr(build_intervals)

# Only intervals starting from "missing" are at-risk periods for the
# dead/alive competing-risks question; an interval that starts already
# resolved (shouldn't occur by construction, since we stop extending once
# resolved - guarded here anyway) is dropped.
intervals <- intervals %>%
  filter(status_from == "missing", Tstart < Tstop) %>%
  mutate(
    to_dead = as.integer(status_to == "dead"),
    to_alive = as.integer(status_to %in% c("alive", "prisoner")),
    # once a person resolves, no further interval should exist for them;
    # if status_to is itself "missing" (still unresolved at that reading),
    # both event indicators are 0 (censored at Tstop, matching the
    # single-window design's censoring convention)
    cohort = factor(year)
  )

message(sprintf("Counting-process intervals built: %d (from %d people)",
                nrow(intervals), n_distinct(intervals$person_id)))
message(sprintf("Events: %d dead, %d alive/prisoner, %d censored (still missing)",
                sum(intervals$to_dead), sum(intervals$to_alive),
                sum(intervals$to_dead == 0 & intervals$to_alive == 0)))

# Sanity check: a person should never appear as "still at risk" in a LATER
# interval after already resolving in an earlier one within their own chain
resolved_then_more <- intervals %>%
  arrange(person_id, Tstart) %>%
  group_by(person_id) %>%
  mutate(prior_resolved = lag(cumsum(to_dead + to_alive) > 0, default = FALSE)) %>%
  ungroup() %>%
  filter(prior_resolved)
message(sprintf("QC: intervals continuing after a prior resolution (should be 0): %d",
                nrow(resolved_then_more)))
stopifnot(nrow(resolved_then_more) == 0)

write_rds(intervals, "data_inter/ukr_ualosses_intervals_countingprocess.rds")

message("\n✓ Person trajectories and counting-process intervals saved.")
message("  data_inter/ukr_ualosses_person_trajectories.rds")
message("  data_inter/ukr_ualosses_intervals_countingprocess.rds")
