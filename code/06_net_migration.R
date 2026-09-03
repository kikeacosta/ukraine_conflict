# ==============================================================================
# STEP 06 - Net emigration by age and sex, 2022-2025
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Builds the refugee outflow that the projection subtracts from exposure.
#
#   1. Eurostat temporary-protection registrations give the STOCK of Ukrainian
#      refugees in the EU by broad age group and sex, for each year.
#   2. Registrations of unknown sex or unknown age are set aside, and the
#      distribution observed among the known cases is used to allocate each
#      country's total across age and sex.
#   3. The broad age groups are ungrouped to single years of age with a
#      penalised composite link model (pclm).
#   4. Totals are rescaled so the 2025 stock matches the UNHCR global figure
#      of 5,923,870 Ukrainian refugees.
#   5. Year-on-year differences in the stock give the annual NET FLOW, which
#      is what step 11 consumes. A negative flow (as in 2024) means net return.
#
# INPUTS   data_input/refugees_eurostat/migr_asytpsm_*.xlsx
#          data_input/unhcr_refugees_annual.csv
# OUTPUT   data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds
#          <- used by 11 and 13
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# according to UNHCR, as of February 2026, there are 5.923.870 refugees from Ukraine reported globally
# https://data.unhcr.org/en/situations/ukraine
lst <- 5923870

# eurostat 2025
dt <-
  read_xlsx(
    "data_input/refugees_eurostat/migr_asytpsm__custom_21093689_page_spreadsheet.xlsx",
    skip = 10
  )

unique(dt$`GEO (Labels)`)
cts_rmv <- c("Special value", ":", "Observation flags:", "dp", "d")

cols_names <-
  expand_grid(
    date = 2022:2025,
    sex = c("t", "m", "f", "u"),
    dat = c("data", "flag")
  ) %>%
  mutate(nms = paste(date, sex, dat, sep = "_")) %>%
  pull(nms)

all_names <- c("ctr", "age", cols_names)

dt2 <-
  dt |>
  set_names(all_names) |>
  filter(!ctr %in% cts_rmv) |>
  drop_na(ctr) |>
  select(-ends_with("flag")) |>
  pivot_longer(-c(ctr, age), names_to = "vars", values_to = "mix") |>
  mutate(year = str_sub(vars, 1, 4), sex = str_sub(vars, 6, 6)) |>
  select(ctr, year, sex, age, mix)

# total EU
tot_eu <-
  dt2 |>
  filter(
    ctr == "European Union - 27 countries (from 2020)",
    sex == "t",
    age == "Total"
  )

tot_cts <-
  dt2 |>
  filter(
    ctr != "European Union - 27 countries (from 2020)",
    sex == "t",
    age == "Total"
  ) |>
  mutate(mix = mix |> as.double()) |>
  replace_na(list(mix = 0)) |>
  summarise(mix_tot = sum(mix), .by = c(year))

# age and sex distribution among known variables
dst_sex_age <-
  dt2 |>
  filter(
    ctr != "European Union - 27 countries (from 2020)",
    !sex %in% c("t", "u"),
    !age %in% c("Total", "Unknown")
  ) |>
  mutate(mix = mix |> as.double()) |>
  replace_na(list(mix = 0)) |>
  summarise(mix = sum(mix), .by = c(sex, age, year)) |>
  arrange(year, sex, age) |>
  group_by(year) |>
  mutate(cx = mix / sum(mix)) |>
  ungroup() |>
  mutate(
    age = case_when(
      age == "Less than 14 years" ~ 0,
      age == "From 14 to 17 years" ~ 14,
      age == "From 18 to 34 years" ~ 18,
      age == "From 35 to 64 years" ~ 35,
      age == "65 years or over" ~ 65
    )
  ) |>
  arrange(year, sex, age) |>
  select(-mix)

# imputing age and sex
dt3 <-
  dst_sex_age |>
  left_join(tot_cts) |>
  mutate(mix = cx * mix_tot, year = as.integer(year)) |>
  select(-cx, -mix_tot)

dt3 |>
  summarise(mix = sum(mix), .by = c(year))

copy_this(
  dt3 |>
    summarise(mix = sum(mix), .by = c(year))
)

# ungrouping in sigle-year ages
chunk <-
  dt3 %>%
  filter(year == 2022, sex == "f")

ung_age_mig <- function(chunk) {
  dt_in <-
    tibble(age = chunk$age, mix = chunk$mix, mix_mt = mix * 1e5) %>%
    mutate(mix_mt = ifelse(mix_mt == 0, 1, mix_mt))
  nl <- 36
  mix <- pclm(x = dt_in$age, y = dt_in$mix_mt, nlast = nl)$fitted

  fit <- tibble(age = 0:100, mix = round(mix / 1e5))

  out <-
    chunk %>%
    select(year, sex) %>%
    unique() %>%
    left_join(fit, by = character())
  return(out)
}

dt4 <-
  dt3 %>%
  group_by(year, sex) %>%
  do(ung_age_mig(chunk = .data)) %>%
  ungroup()

# testing ungrouping consistency
dt3 |>
  summarise(mix = sum(mix), .by = c(year))
dt4 |>
  summarise(mix = sum(mix), .by = c(year))

dt4 |>
  mutate(mix = ifelse(sex == "m", -mix, mix)) |>
  ggplot() +
  geom_line(aes(age, mix, col = year, group = interaction(sex, year))) +
  geom_hline(yintercept = 0, lty = "dashed") +
  coord_flip() +
  theme_bw()

dt4 |>
  summarise(mix = sum(mix), .by = c(year)) |>
  filter(year == 2025) |>
  mutate(mix_tot = lst) |>
  mutate(adj = mix_tot / mix)

# plotting both

# Standardize Grouped Data
dt3_standardized <- dt3 %>%
  arrange(sex, year, age) %>%
  group_by(sex, year) %>%
  mutate(
    next_age = lead(age),
    next_age = if_else(is.na(next_age), 100, next_age),
    interval_width = next_age - age,
    standardized_mix = mix / interval_width,
    plot_mix = if_else(sex == "m", -standardized_mix, standardized_mix)
  ) %>%
  ungroup() |>
  mutate(
    sex = case_when(
      sex == "f" ~ "Females",
      sex == "m" ~ "Males"
    )
  )

# Single-Year data
dt4_plot <-
  dt4 %>%
  mutate(
    plot_mix = if_else(sex == "m", -mix, mix)
  ) |>
  mutate(
    sex = case_when(
      sex == "f" ~ "Females",
      sex == "m" ~ "Males"
    )
  )


# Superposed Pyramid plot

ggplot() +
  facet_grid(. ~ year) +
  geom_rect(
    data = dt3_standardized,
    aes(
      xmin = age, # Starts at the interval base
      xmax = next_age, # Ends at the next interval start
      ymin = 0,
      ymax = plot_mix,
      fill = sex
    ),
    alpha = 0.35,
    color = "white", # Adds a subtle border to distinguish the brackets
    linewidth = 0.2
  ) +
  geom_line(
    data = dt4_plot,
    aes(x = age + 0.5, y = plot_mix, color = sex), # mid-interval adjustment for single years
    linewidth = 0.8
  ) +
  coord_flip() +
  scale_x_continuous(
    labels = abs,
    breaks = seq(0, 100, 10)
  ) +
  scale_fill_manual(values = c("Females" = "#cf6a87", "Males" = "#574b90")) +
  scale_color_manual(values = c("Females" = "#b83b5e", "Males" = "#3f3280")) +
  labs(
    x = "Age",
    y = "Ukraininan refugees in EU",
    fill = "Sex",
    color = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("figures/exploratory/migs_pyramids_stock.png", w = 6, h = 3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# loading data on total refugees
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
an <- read_csv("data_input/unhcr_refugees_annual.csv")

an2 <-
  an |>
  select(-5:-2) |>
  rename(year = 1, mix_r = 2, mix_a = 3) |>
  mutate(mix_tot = mix_r + mix_a) |>
  select(-mix_r, -mix_a) |>
  arrange(year)

copy_this(an2)


# adjusting the total to match UNHCR data
# there is UNHCR data by year and a total as of feb 2026, but they do not match,
# so I will adjust the total to match the UNHCR total and keep the distribution by year as is
adj_fct <-
  dt4 |>
  summarise(mix = sum(mix), .by = c(year)) |>
  filter(year == 2025) |>
  mutate(adj = lst / mix) |>
  pull(adj)

# looking at adjustments to match UNHCR total by year
dt5_old <-
  dt4 |>
  mutate(mix_sum = sum(mix), .by = c(year)) |>
  left_join(an2) |>
  mutate(mix = mix * mix_tot / mix_sum) |>
  select(-mix_sum, -mix_tot) |>
  spread(year, mix) |>
  mutate(
    e2022 = -`2022`,
    e2023 = `2022` - `2023`,
    e2024 = `2023` - `2024`,
    e2025 = `2024` - `2025`
  ) |>
  select(sex, age, e2022, e2023, e2024, e2025) |>
  pivot_longer(-c(sex, age), names_to = "year", values_to = "mix") |>
  mutate(year = str_sub(year, 2, 5) |> as.integer(), mix = round(mix, 0))

dt5_old |>
  summarise(mix = sum(mix), .by = c(year))

# adjusting to match UNHCR total final
dt5 <-
  dt4 |>
  mutate(mix = mix * adj_fct) |>
  spread(year, mix) |>
  mutate(
    e2022 = -`2022`,
    e2023 = `2022` - `2023`,
    e2024 = `2023` - `2024`,
    e2025 = `2024` - `2025`
  ) |>
  select(sex, age, e2022, e2023, e2024, e2025) |>
  pivot_longer(-c(sex, age), names_to = "year", values_to = "mix") |>
  mutate(year = str_sub(year, 2, 5) |> as.integer(), mix = round(mix, 0))

dt5 |>
  mutate(mix = ifelse(sex == "m", -mix, mix)) |>
  ggplot() +
  geom_line(aes(age, mix, col = year, group = interaction(sex, year))) +
  geom_hline(yintercept = 0, lty = "dashed") +
  coord_flip() +
  theme_bw()

dt5 |>
  ggplot() +
  geom_line(aes(age, mix, col = factor(year), group = year)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  facet_grid(~sex) +
  theme_bw()

# saving estimates
write_rds(dt5, "data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds")

dt5 <- read_rds(
  "data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds"
)

dt5 |>
  mutate(mix = ifelse(sex == "m", -mix, mix)) |>
  ggplot() +
  geom_line(aes(age, mix, col = year, group = interaction(sex, year))) +
  geom_hline(yintercept = 0, lty = "dashed") +
  coord_flip() +
  theme_bw()

cols <- rev(c("grey60", "grey40", "grey20", "black"))
dt5 |>
  mutate(
    sex = case_when(sex == "m" ~ "Male", sex == "f" ~ "Female"),
    sex = factor(sex, levels = c("Male", "Female"))
  ) |>
  ggplot() +
  geom_line(aes(age, mix, col = factor(year), group = year)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  facet_grid(~sex) +
  coord_cartesian(ylim = c(-60000, 10000)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(
    breaks = seq(-60000, 10000, 10000),
    labels = scales::comma
  ) +
  scale_color_manual(values = cols) +
  labs(color = "Year") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(
  "figures/exploratory/labtalk/ukr_migrants_unchr_eurostat_sex_age_2022_2025.png",
  w = 7,
  h = 3.5
)


# testing net migration with grouped and ungrouped data ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt3 |>
  spread(year, mix) |>
  mutate(
    e2022 = -`2022`,
    e2023 = `2022` - `2023`,
    e2024 = `2023` - `2024`,
    e2025 = `2024` - `2025`
  ) |>
  select(sex, age, e2022, e2023, e2024, e2025)

dt4 |>
  mutate(
    age = case_when(
      age %in% 0:13 ~ 0,
      age %in% 14:17 ~ 14,
      age %in% 18:34 ~ 18,
      age %in% 35:64 ~ 35,
      age >= 65 ~ 65
    )
  ) |>
  summarise(mix = sum(mix), .by = c(year, sex, age)) |>
  spread(year, mix) |>
  mutate(
    e2022 = -`2022`,
    e2023 = `2022` - `2023`,
    e2024 = `2023` - `2024`,
    e2025 = `2024` - `2025`
  ) |>
  select(sex, age, e2022, e2023, e2024, e2025)
