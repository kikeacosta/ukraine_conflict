# ==============================================================================
# STEP 09h - A sample of links and drop-outs for clerical review
# ==============================================================================
#
# WHY
# ---
# Two rules of the linkage rest on assumptions no count can check: that a
# corrected key - same surname and first name and the same date of birth or
# event month, or the same date of birth and event date and the same surname or
# first name (ual_match_keys()) - is the same person (M5), and that a missing
# person who leaves the register beyond list maintenance is alive (M3). A
# clerical review settles both: two people check each case by hand against the
# register's source pages (the person's page on ualosses.org, the Ministry of
# Internal Affairs' database of missing persons, the obituaries the register
# cites).
#
# THE SAMPLE (simple random within each kind, fixed seed)
#   200 exact-key links     a missing person found in a later release under the
#                           same full name and date of birth
#   200 corrected-key links found under a corrected key: 100 by each of its two
#                           rules, the same name and the re-spelt name (all of
#                           them if fewer), so that each rule has its own rate
#                           of false links
#   200 drop-outs           absent from a release and from every later one
# Each row gives the record at entry and, for a link, the record it was linked
# to, with empty columns for the reviewers' verdicts.
#
# PERSONAL DATA. The output holds names and dates of birth. It is written to
# documents/team/clerical_review/, which is gitignored with the rest of documents/,
# and must never be committed or shared outside the team.
#
# Besides the sample, the full list of the 2022 cohort's drop-outs in the first
# window (February 2022 left out), whose reading as clean-ups rather than
# survivals is the one judgement of the imputation that moves the 2022 loss.
#
# INPUTS   the six register releases (ual_releases, 00_setup.R; gitignored)
# OUTPUTS  documents/team/clerical_review/clerical_review_sample.csv (local only)
#          documents/team/clerical_review/removals_2022_v14_v15.csv (local only)
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

out_dir <- "documents/team/clerical_review"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# what a reviewer needs to find the person: the linkage fields, and place of
# origin, unit, rank and place of the event
review_fields <- function(path, release) {
  read_xlsx(require_raw(path), sheet = "Database", col_types = "text") |>
    transmute(
      dob = excel_date(DateBirth),
      key = paste(LastName, FirstName, Patronym, dob),
      last_name = LastName, first_name = FirstName, patronym = Patronym,
      date_of_birth = dob, date_event = excel_date(DateEvent),
      place_of_origin = From, unit = Unit, rank = Rank, place_event = LocationEvent,
      status = Status
    ) |>
    arrange(match(status, ual_status_order)) |>
    distinct(key, .keep_all = TRUE) |>
    mutate(release = release) |>
    select(-dob)
}
regs <- set_names(map(ual_releases$path, ual_read_release), ual_releases$release)
hist <- ual_follow_missing(regs, keep_keys = TRUE)
records <- map2_dfr(ual_releases$path, ual_releases$release, review_fields)

set.seed(2026)
take <- function(d, n) d |> slice_sample(n = min(n, nrow(d)))
rule_label <- c(name = "same name", dates = "re-spelt name")
links <- function(kind, n = 200, by_rule = NULL) {
  hist |>
    filter(found == kind, is.null(by_rule) | rule %in% by_rule) |>
    take(n) |>
    transmute(case = if (is.null(by_rule)) kind else paste0(kind, ", ", rule_label[[by_rule]]),
              pid, entry, release, entry_key, linked_key = key)
}
dropouts <-
  ual_windows(hist) |>
  filter(to == "no_longer_listed") |>
  take(200) |>
  left_join(hist |> distinct(pid, entry_key), by = "pid") |>
  transmute(case = "drop-out", pid, entry, release, entry_key, linked_key = NA_character_)

sample_cases <-
  bind_rows(links("exact"), links("corrected", 100, "name"), links("corrected", 100, "dates"), dropouts) |>
  mutate(case_id = sprintf("C%03d", row_number()), .before = 1) |>
  left_join(records |> rename_with(\(x) paste0("entry_", x), -key), by = c("entry_key" = "key", "entry" = "entry_release")) |>
  left_join(records |> rename_with(\(x) paste0("linked_", x), -key), by = c("linked_key" = "key", "release" = "linked_release")) |>
  select(-pid, -entry_key, -linked_key) |>
  mutate(same_person = "", fate_found = "", source_checked = "", reviewer = "", notes = "")
stopifnot(nrow(sample_cases) <= 600, !anyNA(sample_cases$entry_last_name))

write_excel_csv(sample_cases, file.path(out_dir, "clerical_review_sample.csv"))
cat("clerical review sample:", nrow(sample_cases), "cases\n")
print(count(sample_cases, case))

# THE 2022 COHORT'S FIRST-WINDOW REMOVALS, IN FULL. The estimate reads the
# drop-outs of the events of 2022 as clean-ups of their records, not as people
# found alive (composition_bound(), 00_setup.R); read as survivals they would
# lower the 2022 male loss by about 0.2 years (Table A22). Most of them fall in
# the first window, v14 to v15, and they are not re-spellings: the corrected
# keys find almost none. They are few enough to check one by one against the
# register's source pages, which is the only test of that reading available.
# February 2022 events are left out, as they are of the bound.
removals_2022 <-
  ual_windows(hist) |>
  filter(to == "no_longer_listed", year == 2022, from_release == ual_releases$release[1],
         floor_date(date_evnt, "month") != as.Date("2022-02-01")) |>
  left_join(hist |> distinct(pid, entry_key), by = "pid") |>
  transmute(case = "removal, 2022 events, v14 to v15", pid, entry, release, entry_key) |>
  mutate(case_id = sprintf("R%03d", row_number()), .before = 1) |>
  left_join(records |> rename_with(\(x) paste0("entry_", x), -key), by = c("entry_key" = "key", "entry" = "entry_release")) |>
  select(-pid, -entry_key) |>
  mutate(fate_found = "", source_checked = "", reviewer = "", notes = "")
stopifnot(!anyNA(removals_2022$entry_last_name))
write_excel_csv(removals_2022, file.path(out_dir, "removals_2022_v14_v15.csv"))
cat("removals of 2022 events in the first window, February 2022 left out:", nrow(removals_2022), "\n")
