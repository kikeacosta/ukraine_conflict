rm(list = ls())
gc()
source("code/00_setup.R")

# data compiled by Olivier Hubert
# in https://www.kaggle.com/datasets/ol4ubert/confirmed-ukrainian-military-personnel-losses
#
# How do the personnel recorded as "missing" resolve? Estimated by following
# individuals from register v14 (16 Sep 2025) into register v19 (19 Sep 2026),
# twelve months later - the same one-year spacing as the event-year cohorts the
# chain below uses to stand in for duration since disappearance.
#
# Both registers are large individual-level files (44 MB and 31 MB) that are
# not tracked in git. Only the anonymous count summaries below cross the cache
# boundary, so the repository carries no names or dates of birth.

file_v14 <- "data_input/ualosses_hubert_datasets/250916_UKR_ualosses_Personnel_v14.xlsx"
file_v19 <- "data_input/ualosses_hubert_datasets/260919_UKR_ualosses_Personnel_v19.xlsx"

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
# individual-level matching -> anonymous transition counts (cached)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
transitions <- cache_rds(
  "data_inter/ualosses_status_transitions.rds",
  {
    read_reg <- function(path) {
      read_xlsx(require_raw(path), sheet = "Database") |>
        mutate(
          name2 = paste(LastName, FirstName, Patronym, sep = " "),
          date_bth = excel_date(DateBirth),
          date_evnt = excel_date(DateEvent),
          year = year(date_evnt),
          ukrainian = ual_is_ukrainian(Nationality)
        ) |>
        select(name2, date_bth, year, date_evnt, ukrainian, status = Status)
    }

    # The population at risk is defined on the EARLIER register only: people
    # listed there as missing (or prisoner) with an event in 2022-2025. The
    # later register is searched in full. Filtering it on event year and
    # nationality too lost every person whose later record had a corrected
    # event date or a shifted nationality field - about 460 people, some 360
    # of whom had been found dead - and counted them as still missing.
    v14 <- read_reg(file_v14) |> filter(year %in% 2022:2025, ukrainian)
    v19 <- read_reg(file_v19)

    # one outcome per person: name + date of birth can collide, and a
    # many-to-many join would count those individuals more than once. Where
    # the same key appears under two statuses the most resolved one is kept:
    # a death record, or a return from captivity, is newer information than a
    # missing one.
    v19_lu <-
      v19 |>
      arrange(match(status, c("dead", "released_prisoner", "prisoner", "missing"))) |>
      distinct(name2, date_bth, .keep_all = TRUE)
    message(
      "  v19 rows: ", nrow(v19),
      " | unique name+dob keys: ", nrow(v19_lu),
      " (", nrow(v19) - nrow(v19_lu), " collapsed)"
    )

    cmp <-
      v14 |>
      select(name2, date_bth, date_evnt, status) |>
      left_join(
        v19_lu |> select(name2, date_bth, status2 = status),
        by = c("name2", "date_bth")
      ) |>
      # A name+DOB not found in v19 is treated as CENSORED (still missing),
      # not resolved alive. Absence is not a sign of survival: the dead
      # drop out of the register at least as often as the missing do (4.0%
      # of v14's dead and 3.1% of its missing are not found in v19 under the
      # same name and date of birth), and many of those not found are still
      # listed under a corrected name or date of birth.
      # v19 records returns from captivity as released_prisoner: an alive
      # resolution, which the chain calls "alive".
      mutate(
        status2 = ifelse(is.na(status2), "missing", status2),
        status2 = ifelse(status2 == "released_prisoner", "alive", status2),
        year = year(date_evnt)
      )

    bind_rows(
      cmp |> filter(status == "missing") |> count(year, status2) |>
        mutate(from = "missing"),
      cmp |> filter(status == "prisoner") |> count(year, status2) |>
        mutate(from = "prisoner")
    )
  }
)

# Stocks by event year, taken from the OUTPUT OF STEP 08 rather than re-read
# from the raw register.
#
# This matters: step 08 redistributes the ~3,300 records whose age or event
# year is not recorded. Counting the raw file again here would silently drop
# them, so the confirmed-death total used in this script would fall short of
# the one used everywhere downstream by those records, and the manuscript
# tables built in 15 would disagree with each other.
stocks <-
  read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") |>
  summarise(n = sum(dx), .by = c(year, status)) |>
  mutate(status = as.character(status)) |>
  arrange(year, status)

print(stocks |> pivot_wider(names_from = status, values_from = n))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# resolution rates of the missing, by cohort year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tasas_long <-
  transitions |>
  filter(from == "missing") |>
  mutate(prop = n / sum(n), .by = year)

print(tasas_long |> arrange(year, status2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# synthetic-cohort Markov imputation of the missing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# impute_missing() (00_setup.R) holds the chain itself, extracted here so
# this call and the alpha sensitivity in 13b share one formula instead of
# two copies that could drift apart. A "step" in the chain is the twelve-month
# window between v14 and v19, and event-year cohorts, one year apart, stand in
# for duration since disappearance.
stock_missing_2026 <-
  stocks |>
  filter(status == "missing") |>
  select(year, missing_stock = n)

# the range of alpha the evidence allows (00_setup.R sets out the sources).
# Prisoners are alive whether still held or released, so both statuses count
# as the prisoners the register knows of.
register_alive <- sum(stocks$n[stocks$status %in% c("prisoner", "released_prisoner")])
alpha_range <- alpha_evidence(tasas_long, stock_missing_2026, register_alive)
print(as.data.frame(alpha_range))
write_rds(alpha_range, "data_inter/ukr_alpha_missing.rds")

imputation_final <- impute_missing(alpha_range$alpha_mode, tasas_long, stock_missing_2026)

print(imputation_final)

confirmados_df <-
  stocks |>
  filter(status == "dead") |>
  select(year, confirmados_stock = n)

# The military total in each year is linear in alpha: at_alpha0 - alpha x
# residual. 10 builds the combatant bounds from these two columns and 11
# draws alpha, so both use exactly the chain run here.
military_alpha_lines <-
  confirmados_df |>
  left_join(stock_missing_2026, by = "year") |>
  left_join(
    impute_missing(0, tasas_long, stock_missing_2026) |> select(year, dead0 = imputed_dead),
    by = "year"
  ) |>
  left_join(
    impute_missing(1, tasas_long, stock_missing_2026) |> select(year, dead1 = imputed_dead),
    by = "year"
  ) |>
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
write_rds(tasas_long, "data_inter/ukr_ualosses_transition_rates.rds")

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
