# ==============================================================================
# STEP 09F - PER-PERSON TRAJECTORIES ACROSS ALL REGISTER RELEASES
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Builds, for the multistate analysis in 09g, one history per person recorded
# as missing, followed across the four register releases:
#   v14 (16 Sep 2025), v16 (4 Dec 2025), v18 (23 Apr 2026), v19 (19 Sep 2026).
# Each history becomes counting-process intervals (Tstart, Tstop] in real days
# since the date of the event, so a person's later interval starts where the
# earlier one ended instead of restarting the clock.
#
# LINKAGE AND OUTCOMES: the same as production step 09, through the same
# functions in 00_setup.R (ual_follow_missing, ual_windows):
#   - a person enters at their first listing as missing, event in 2022-2025
#     and Ukrainian, in v14, v16 or v18; v19 admits no entrants, since nothing
#     follows it;
#   - each later release is searched in full, by exact full name and date of
#     birth, then under a corrected key;
#   - a person no longer listed - absent from a release and from every later
#     one - is resolved alive at the first absence; absent from one release
#     but listed again later, they were still missing;
#   - a person whose first resolution is a return from captivity was held,
#     not disappeared, and is left out altogether;
#   - outcomes kept apart: dead, prisoner, no longer listed. A history stops
#     at its first resolution.
#
# PERSONAL DATA. The registers hold names and dates of birth. Nothing written
# here carries either: persons are identified by a sequential number, and the
# files are gitignored (individual-level rows, not aggregate counts).
#
# INPUT   the four register releases (gitignored)
# OUTPUT  data_inter/ukr_ualosses_person_trajectories.rds      (gitignored)
#         data_inter/ukr_ualosses_intervals_countingprocess.rds (gitignored)
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

message("Reading registers...")
regs <- map(set_names(ual_releases$path, ual_releases$release), ual_read_release)
hist <- ual_follow_missing(regs)
rm(regs)

message(sprintf("  entrants: %s", paste(
  hist |> distinct(pid, entry) |> count(entry) |> with(paste(entry, n)), collapse = " | "
)))

# ==============================================================================
# COUNTING-PROCESS INTERVALS
# ==============================================================================
rel_date <- set_names(ual_releases$date, ual_releases$release)
ids <- hist |> distinct(pid) |> mutate(person_id = row_number())
intervals <-
  ual_windows(hist, "alive") |>
  left_join(ids, by = "pid") |>
  transmute(
    person_id, year,
    Tstart = as.numeric(rel_date[from_release] - date_evnt),
    Tstop = as.numeric(rel_date[release] - date_evnt),
    status_to = to,
    to_dead = as.integer(to == "dead"),
    to_prisoner = as.integer(to == "prisoner"),
    to_unlisted = as.integer(to == "no_longer_listed"),
    cohort = factor(year)
  ) |>
  filter(Tstart < Tstop)

# a person must never be at risk again after resolving
stopifnot(
  intervals |> summarise(r = sum(to_dead + to_prisoner + to_unlisted), .by = person_id) |>
    with(all(r <= 1))
)

message(sprintf(
  "Intervals: %d from %d people | dead %d, prisoner %d, no longer listed %d, still missing %d",
  nrow(intervals), n_distinct(intervals$person_id),
  sum(intervals$to_dead), sum(intervals$to_prisoner),
  sum(intervals$to_unlisted),
  sum(intervals$to_dead + intervals$to_prisoner + intervals$to_unlisted == 0)
))

# no names or dates of birth leave this script
write_rds(hist |> left_join(ids, by = "pid") |> select(-pid),
          "data_inter/ukr_ualosses_person_trajectories.rds")
write_rds(intervals, "data_inter/ukr_ualosses_intervals_countingprocess.rds")

message("Done: trajectories and intervals written (gitignored).")
