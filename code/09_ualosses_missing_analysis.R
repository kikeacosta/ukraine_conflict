rm(list = ls())
gc()
source("code/00_setup.R")

# data compiled by Olivier Hubert
# in https://www.kaggle.com/datasets/ol4ubert/confirmed-ukrainian-military-personnel-losses
#
# How do the personnel recorded as "missing" resolve? Estimated by following
# everyone listed as missing through the four register releases - v14 (16 Sep
# 2025), v16 (4 Dec 2025), v18 (23 Apr 2026) and v19 (19 Sep 2026) - and
# fitting, to the three windows between them, how fast the missing leave that
# status by months since the event. Each event month's missing in v19 are then
# carried forward from the duration they have reached.
#
# The releases are large individual-level files that are not tracked in git.
# Only the anonymous count summaries below cross the cache boundary, so the
# repository carries no names or dates of birth.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KEY ASSUMPTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# alpha, the share of the long-term missing (those still unresolved at the end
# of the chain) who are alive, drives the combatant total. Its range comes from
# alpha_evidence() in 00_setup.R, which compares the prisoners the register
# records with official prisoner-of-war figures and sets out the sources. This script imputes
# at the central value and records the range, 10 turns the range into the
# combatant bounds, and 11 draws alpha within it.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# individual-level linkage -> anonymous window counts (cached)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The rules live in 00_setup.R (ual_follow_missing, ual_windows):
#   - a person enters at their first listing as missing, event in 2022-2025
#     and Ukrainian, in v14, v16 or v18, and is looked up in every later
#     release searched in full, whatever the event year or nationality field;
#   - by exact full name and date of birth; if absent, under a corrected key:
#     a key new in that release with the same surname and first name and the
#     same date of birth or event month;
#   - a person no longer listed - absent from a release and from every later
#     one - is resolved alive at the first absence. Absent from one release
#     but listed again later, they were still missing;
#   - a history stops at its first resolution;
#   - a person whose first resolution is a return from captivity was held,
#     not disappeared, and is left out of the population at risk altogether.
linkage <- cache_rds(
  "data_inter/ualosses_window_transitions.rds",
  {
    regs <- map(set_names(ual_releases$path, ual_releases$release), ual_read_release)
    hist <- ual_follow_missing(regs)

    # The dead drop out of the register too, and cannot have been found alive:
    # v14's dead with no trace in v19, by the same two searches, measure how
    # often a record disappears through list maintenance.
    lu14 <- ual_one_per_key(regs$v14)
    lu19 <- ual_one_per_key(regs$v19)
    dead14 <-
      regs$v14 |>
      filter(status == "dead", year %in% 2022:2025, ukrainian) |>
      distinct(key, .keep_all = TRUE)
    dead_absent <- dead14 |> anti_join(lu19, by = "key")
    dead_corrected <- ual_corrected_key(dead_absent, lu19, lu14)

    bind_rows(
      ual_windows(hist, "alive") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "windows", rule = "alive"),
      ual_windows(hist, "missing") |>
        count(year, month = floor_date(date_evnt, "month"), entry, from_release, to) |>
        mutate(table = "windows", rule = "missing"),
      hist |> count(year, entry, release, to = status, found) |> mutate(table = "found"),
      # each person's sequence of statuses, for the documentation's counts of
      # people absent from one release and listed again
      hist |>
        arrange(pid, match(release, ual_releases$release)) |>
        summarise(to = paste(coalesce(status, "absent"), collapse = " > "), .by = c(pid, year, entry)) |>
        count(year, entry, to) |>
        mutate(table = "patterns"),
      dead14 |>
        mutate(found = case_when(!key %in% dead_absent$key ~ "exact",
                                 key %in% dead_corrected$key ~ "corrected",
                                 .default = "absent")) |>
        count(year, found) |>
        mutate(table = "dead_v14_in_v19")
    )
  }
)

windows <- linkage |> filter(table == "windows") |> select(rule, year, month, entry, from_release, to, n)
at_risk_v14 <-
  windows |>
  filter(rule == "alive", from_release == "v14") |>
  summarise(at_risk = sum(n), .by = year)

# twelve-month probabilities in the chain's states; n is prop times the
# cohort's missing in v14, the weight alpha_evidence() gives each cohort
chain_rates <- function(w) {
  ual_compose(w) |>
    mutate(status2 = ual_chain_state(to)) |>
    summarise(prop = sum(prop), .by = c(year, status2)) |>
    left_join(at_risk_v14, by = "year") |>
    transmute(year, status2, n = prop * at_risk, prop, from = "missing")
}

# Stocks by event year, taken from the OUTPUT OF STEP 08 rather than re-read
# from the raw register.
#
# This matters: step 08 redistributes the ~3,300 records whose age or event
# year is not recorded. Counting the raw file again here would silently drop
# them, so the confirmed-death total used in this script would fall short of
# the one used everywhere downstream by those records, and the manuscript
# tables built in 15 would disagree with each other.
stocks_registered <-
  read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") |>
  summarise(n = sum(dx), .by = c(year, status)) |>
  mutate(status = as.character(status)) |>
  arrange(year, status)

# By event month, and completed for registration lag (08b). Each year's
# registered count is spread over its months as v19 lists them, and each month
# is brought to the completeness events reach at four years by its own factor.
# The deaths added are late registrations: people not yet in the register at
# all. Prisoners and released prisoners are left as recorded.
completion <- read_rds("data_inter/ukr_registration_completion.rds")
stock_month <-
  completion$month |>
  mutate(share = n_v19 / sum(n_v19), .by = c(year, status)) |>
  left_join(stocks_registered |> rename(n_year = n), by = c("year", "status")) |>
  transmute(month = m, year, status, registered = n_year * share, completed = registered * factor)
stopifnot(isTRUE(all.equal(
  stock_month |> summarise(r = sum(registered), .by = c(year, status)) |> arrange(year, status) |> pull(r),
  stocks_registered |> filter(status %in% c("dead", "missing")) |> arrange(year, status) |> pull(n)
)))
stocks <-
  bind_rows(
    stocks_registered |> filter(!status %in% c("dead", "missing")),
    stock_month |> summarise(n = sum(completed), .by = c(year, status))
  ) |>
  arrange(year, status)

print(stocks |> pivot_wider(names_from = status, values_from = n))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# resolution of the missing over twelve months, by cohort year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# every outcome kept apart, with the people at risk: those listed as missing
# in v14, and those first listed in v16 or v18, who join the later windows
later_entrants <-
  windows |>
  filter(rule == "alive", entry != "v14", from_release == entry) |>
  summarise(later_entrants = sum(n), .by = year)
resolution_12m <-
  ual_compose(windows |> filter(rule == "alive")) |>
  pivot_wider(names_from = to, values_from = prop) |>
  left_join(at_risk_v14, by = "year") |>
  left_join(later_entrants, by = "year") |>
  select(year, at_risk_v14 = at_risk, later_entrants, all_of(ual_outcomes))
print(as.data.frame(resolution_12m |> mutate(across(all_of(ual_outcomes), \(x) round(100 * x, 2)))))

# the same composed by event-year cohort, in the states of the earlier chain,
# kept for the comparison of designs below
tasas_long <- chain_rates(windows |> filter(rule == "alive"))
stopifnot(tasas_long |> summarise(p = sum(prop), .by = year) |> with(all(abs(p - 1) < 1e-9)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# resolution by duration since disappearance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hazards of leaving "missing" by months since the event (00_setup.R,
# ual_fit_durations), fitted to one cell per event month and window. The
# projection runs to the horizon 08b settled for registration lag, 48 months:
# the longest duration measured on at least a year of event months. Beyond it
# the data rest on a few months of 2022 events and on v19's one-off deletion
# of a batch of February 2022 records.
fit_model <- function(w) {
  ual_fit_durations(ual_duration_cells(w), horizon = completion$L_ref)
}
model <- fit_model(windows |> filter(rule == "alive"))

hazards <-
  as_tibble(model$h) |>
  mutate(band = paste0(head(model$breaks, -1), "-", model$breaks[-1]), .before = 1)
cat("\n=== MONTHLY HAZARDS OF RESOLUTION BY MONTHS SINCE THE EVENT (%) ===\n")
print(as.data.frame(hazards |> mutate(across(-band, \(x) round(100 * x, 3)))))

# the fit against the counts, by event year and window
fit_check <-
  model$fitted |>
  mutate(year = year(month), at_risk = missing + dead + prisoner + no_longer_listed,
         across(c(p_dead, p_prisoner, p_no_longer_listed), \(p) p * at_risk, .names = "fit_{.col}")) |>
  summarise(across(c(at_risk, dead, prisoner, no_longer_listed, starts_with("fit_")), sum),
            .by = c(year, from_release))
cat("\n=== OBSERVED AND FITTED RESOLUTIONS BY COHORT AND WINDOW ===\n")
print(as.data.frame(fit_check |> mutate(across(where(is.double), round))))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# imputation of the missing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# impute_missing() (00_setup.R) holds the projection itself, so this call,
# the alpha sensitivity in 13b and the lag sensitivity in 13d share one
# formula instead of copies that could drift apart.
stock_missing_month <-
  stock_month |>
  filter(status == "missing") |>
  select(month, year, missing_stock = completed)
stock_missing_2026 <-
  stock_missing_month |>
  summarise(missing_stock = sum(missing_stock), .by = year)
impute <- function(a, m = model) impute_missing(a, m, stock_missing_month)

# the range of alpha the evidence allows (00_setup.R sets out the sources).
# Prisoners are alive whether still held or released, so both statuses count
# as the prisoners the register knows of.
register_alive <- sum(stocks$n[stocks$status %in% c("prisoner", "released_prisoner")])
alpha_range <- alpha_evidence(impute, model$resolved, register_alive)
print(as.data.frame(alpha_range))
write_rds(alpha_range, "data_inter/ukr_alpha_missing.rds")

imputation_final <- impute(alpha_range$alpha_mode)

print(imputation_final)

# the dead, registered and completed for registration lag; downstream the
# completed count is the "confirmed" military deaths, spread over the
# registered dead's age-sex profile, and 14 and 15 report the late
# registrations within it apart
confirmados_df <-
  stocks |>
  filter(status == "dead") |>
  select(year, confirmados_stock = n) |>
  left_join(
    stocks_registered |> filter(status == "dead") |> select(year, registered = n),
    by = "year"
  ) |>
  mutate(late_registrations = confirmados_stock - registered)

# The military total in each year is linear in alpha: at_alpha0 - alpha x
# residual. 10 builds the combatant bounds from these two columns and 11
# draws alpha, so both use exactly the chain run here.
military_alpha_lines <-
  confirmados_df |>
  left_join(stock_missing_2026, by = "year") |>
  left_join(impute(0) |> select(year, dead0 = imputed_dead), by = "year") |>
  left_join(impute(1) |> select(year, dead1 = imputed_dead), by = "year") |>
  transmute(
    year,
    confirmed = confirmados_stock,
    missing = missing_stock,
    at_alpha0 = confirmados_stock + dead0,
    residual = dead0 - dead1
  )

stopifnot(isTRUE(all.equal(sum(military_alpha_lines$residual), alpha_range$residual)))
print(military_alpha_lines)
write_rds(military_alpha_lines, "data_inter/ukr_military_alpha_lines.rds")

combined_losses <- confirmados_df %>%
  left_join(imputation_final, by = "year") %>%
  mutate(
    total_estimado = confirmados_stock + imputed_dead,
    year = factor(year)
  )

print("=== BALANCE DE PERDIDAS MORTALES TOTALES ESTIMADAS ===")
print(combined_losses)

# saved for the manuscript tables assembled in 15
write_rds(
  combined_losses |> mutate(year = as.integer(as.character(year))),
  "data_inter/ukr_ualosses_imputation_table.rds"
)
# the fitted model and the missing by event month, for 13b and 13d
write_rds(model, "data_inter/ukr_ualosses_resolution_model.rds")
write_rds(stock_missing_month, "data_inter/ukr_ualosses_missing_by_month.rds")
write_rds(hazards, "data_inter/ukr_ualosses_resolution_hazards.rds")
# the twelve-month rates by event-year cohort, every outcome kept apart, for
# table A2, and in the earlier chain's states
write_rds(resolution_12m, "data_inter/ukr_ualosses_resolution_12m.rds")
write_rds(tasas_long, "data_inter/ukr_ualosses_transition_rates.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how much the linkage rules and the chain matter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The military total under four alternatives, each at the central alpha its
# own evidence gives: a person no longer listed held as still missing rather
# than alive; rates from the people listed as missing in v14 only, without
# those first listed later; the projection run to the longest duration
# observed (about 55 months) instead of 48; and the earlier chain, in which
# event-year cohorts stand in for duration. And how often records vanish: v14's
# missing and v14's dead with no trace in v19 after both searches.
military_under <- function(impute_fn, resolved) {
  a <- alpha_evidence(impute_fn, resolved, register_alive)
  by_year <-
    confirmados_df |>
    left_join(impute_fn(a$alpha_mode) |> select(year, imputed_dead), by = "year") |>
    transmute(year, confirmed = confirmados_stock, total = confirmados_stock + imputed_dead)
  tibble(alpha_mode = a$alpha_mode, alpha_max = a$alpha_max, residual = a$residual,
         military = sum(by_year$total), by_year = list(by_year))
}
model_missing_rule <- fit_model(windows |> filter(rule == "missing"))
model_v14_only <- fit_model(windows |> filter(rule == "alive", entry == "v14"))
cohort_resolved <- with(tasas_long |> filter(status2 != "missing") |> summarise(n = sum(n), .by = status2),
                        set_names(n, status2))
linkage_designs <-
  bind_rows(
    military_under(impute, model$resolved) |>
      mutate(design = "duration model, no longer listed = alive (production)"),
    military_under(\(a) impute(a, model_missing_rule), model_missing_rule$resolved) |>
      mutate(design = "no longer listed = still missing"),
    military_under(\(a) impute(a, model_v14_only), model_v14_only$resolved) |>
      mutate(design = "people listed as missing in v14 only"),
    military_under(\(a) impute(a, modifyList(model, list(horizon = max(model$fitted$d1)))),
                   model$resolved) |>
      mutate(design = "duration model to the longest duration observed"),
    military_under(\(a) impute_missing_cohorts(a, tasas_long, stock_missing_2026), cohort_resolved) |>
      mutate(design = "event-year cohorts standing in for duration")
  ) |>
  relocate(design)
stopifnot(isTRUE(all.equal(linkage_designs$military[1], sum(combined_losses$total_estimado))))

no_trace <-
  bind_rows(
    linkage |> filter(table == "found", entry == "v14", release == "v19") |> mutate(who = "missing"),
    linkage |> filter(table == "dead_v14_in_v19") |> mutate(who = "dead")
  ) |>
  summarise(listed_v14 = sum(n), corrected_key = sum(n[found == "corrected"]),
            no_trace = sum(n[found == "absent"]), .by = c(who, year)) |>
  mutate(pct_no_trace = 100 * no_trace / listed_v14)

patterns <-
  linkage |>
  filter(table == "patterns") |>
  select(year, entry, pattern = to, n) |>
  mutate(first_resolution = map_chr(strsplit(pattern, " > "), function(s) {
    r <- s[!s %in% c("missing", "absent")]
    if (length(r) > 0) r[1] else if (tail(s, 1) == "absent") "no_longer_listed" else "missing"
  }))
# people left out as released prisoners, by cohort and list of entry
released_excluded <-
  patterns |>
  filter(first_resolution == "released_prisoner") |>
  summarise(n = sum(n), .by = c(year, entry))

print(as.data.frame(linkage_designs |> select(-by_year)))
print(as.data.frame(no_trace))
print(as.data.frame(released_excluded |> pivot_wider(names_from = entry, values_from = n, values_fill = 0)))
write_rds(list(designs = linkage_designs, no_trace = no_trace, patterns = patterns,
               released_excluded = released_excluded),
          "data_inter/ukr_ualosses_linkage_checks.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rescale the confirmed-dead age-sex profile up to the imputed total
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ual <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")

adj_fct <-
  ual |>
  summarise(dts = sum(dx), .by = c(year, status)) |>
  filter(status == "dead") |>
  left_join(
    combined_losses |>
      select(year, total_estimado) |>
      mutate(year = as.character(year) |> as.integer()),
    by = "year"
  ) |>
  mutate(adj = total_estimado / dts) |>
  select(year, adj)

print(adj_fct)

ual2 <-
  ual |>
  filter(status == "dead") |>
  left_join(adj_fct, by = "year") |>
  mutate(dx = dx * adj) |>
  select(-adj)

ual2 |>
  summarise(dts = sum(dx), .by = c(year, status))

write_rds(
  ual2,
  "data_inter/ukr_ualosses_conflict_deaths_imputed_miss_sex_age_2022_2025.rds"
)

copy_this(
  combined_losses |>
    mutate(imputed_alive = imputed_alive + imputed_prisoner) |>
    select(
      year,
      confirmed_deaths = confirmados_stock,
      missing = missing_stock,
      imputed_dead,
      imputed_alive,
      total_deaths = total_estimado
    )
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt_plot <-
  combined_losses |>
  mutate(imputed_alive = imputed_alive + imputed_prisoner) |>
  select(
    year,
    confirmed_dead = confirmados_stock,
    imputed_dead,
    imputed_alive
  ) |>
  gather(-year, key = "status", value = "dts") |>
  mutate(
    status = factor(
      status,
      levels = c("imputed_alive", "imputed_dead", "confirmed_dead")
    )
  )

cols <- c("grey30", "#e2606b", "#e63946")

dt_plot |>
  ggplot() +
  geom_bar(
    aes(fill = status, y = dts, x = year),
    position = "stack",
    stat = "identity"
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::comma, breaks = seq(0, 70000, 10000)) +
  labs(y = "Counts", x = "Year", fill = "Status") +
  theme_minimal()
ggsave("figures/exploratory/labtalk/missing_time_qt_imputed.png", w = 6, h = 3)

dt_plot |>
  mutate(
    prop = dts / sum(dts),
    pct_label = scales::percent(prop, accuracy = 0.1),
    .by = year,
  ) |>
  ggplot(aes(x = factor(year), y = dts, fill = status)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(
    aes(label = pct_label),
    position = position_fill(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage", x = "Year", fill = "Status") +
  theme_minimal()
ggsave("figures/exploratory/labtalk/missing_time_imputed.png", w = 6, h = 3)

dt_plot |>
  filter(status != "confirmed_dead") |>
  mutate(
    prop = dts / sum(dts),
    pct_label = scales::percent(prop, accuracy = 0.1),
    .by = year,
  ) |>
  ggplot(aes(x = factor(year), y = dts, fill = status)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(
    aes(label = pct_label),
    position = position_fill(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage", x = "Year", fill = "Status") +
  theme_minimal()
ggsave("figures/exploratory/labtalk/missing_time_imputed_v2.png", w = 6, h = 3)
