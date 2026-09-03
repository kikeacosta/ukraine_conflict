rm(list = ls())
gc()
source("code/00_setup.R")

# data compiled by Olivier Hubert
# in https://www.kaggle.com/datasets/ol4ubert/confirmed-ukrainian-military-personnel-losses
#
# What fraction of the personnel recorded as "missing" is eventually resolved
# to "dead"? Estimated by following individuals from register v14 (Sep 2025)
# into register v18 (May 2026).
#
# Both registers are large individual-level files (44 MB and 27 MB) that are
# not tracked in git. Only the anonymous count summaries below cross the cache
# boundary, so the repository carries no names or dates of birth.

file_v14 <- "data_input/ualosses_hubert_datasets/250916_UKR_ualosses_Personnel_v14.xlsx"
file_v18 <- "data_input/260514_UKR_ualosses_Personnel.xlsx"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KEY ASSUMPTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Share of the long-term missing (those still unresolved at the end of the
# chain) assumed to be alive. This single number drives the combatant death
# "mode" in 10 and therefore the headline estimate, so it is worth a
# sensitivity check before quoting the results.
suelo_vivo <- 0.1
suelo_muerto <- 1 - suelo_vivo

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# individual-level matching -> anonymous transition counts (cached)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
transitions <- cache_parquet(
  "data_inter/ualosses_status_transitions.parquet",
  {
    read_reg <- function(path) {
      read_xlsx(require_raw(path), sheet = "Database") |>
        mutate(
          name2 = paste(LastName, FirstName, Patronym, sep = " "),
          date_bth = ymd("1900-01-01") + as.numeric(DateBirth),
          date_evnt = ymd("1900-01-01") + as.numeric(DateEvent),
          year = year(date_evnt)
        ) |>
        filter(year %in% 2022:2025, Nationality == "Ukraine") |>
        select(name2, date_bth, year, date_evnt, status = Status)
    }

    v14 <- read_reg(file_v14)
    v18 <- read_reg(file_v18)

    # one outcome per person: name + date of birth can collide, and a
    # many-to-many join would count those individuals more than once
    v18_lu <- v18 |> distinct(name2, date_bth, .keep_all = TRUE)
    message(
      "  v18 rows: ", nrow(v18),
      " | unique name+dob keys: ", nrow(v18_lu),
      " (", nrow(v18) - nrow(v18_lu), " collapsed)"
    )

    cmp <-
      v14 |>
      select(name2, date_bth, date_evnt, status) |>
      left_join(
        v18_lu |> select(name2, date_bth, status2 = status),
        by = c("name2", "date_bth")
      ) |>
      # absent from the later register = resurfaced alive
      mutate(
        status2 = ifelse(is.na(status2), "alive", status2),
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
# them, so the confirmed-death total used in this script (83,182) would not
# match the one used everywhere downstream (86,526) — and the manuscript
# tables built in 16 would disagree with each other.
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

# intensity lookup
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# synthetic-cohort Markov imputation of the missing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE: a "step" in this chain is the ~8-month window between the two
# registers, not a calendar year; the event-year cohorts stand in for
# duration since disappearance.
stock_missing_2026 <-
  stocks |>
  filter(status == "missing") |>
  select(year, missing_stock = n)

imputation_final <- stock_missing_2026 %>%
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

print(imputation_final)

confirmados_df <-
  stocks |>
  filter(status == "dead") |>
  select(year, confirmados_stock = n)

combined_losses <- confirmados_df %>%
  left_join(imputation_final, by = "year") %>%
  mutate(
    total_estimado = confirmados_stock + imputed_dead,
    year = factor(year)
  )

print("=== BALANCE DE PERDIDAS MORTALES TOTALES ESTIMADAS ===")
print(combined_losses)

# saved for the manuscript tables assembled in 16
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
