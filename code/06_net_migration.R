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
#   2. Registrations of unknown sex or unknown age are set aside; the distribution
#      observed among the known cases is used to allocate each country's total.
#   3. The broad age groups are ungrouped to single years of age using pclm
#      smoothing, SEPARATELY FOR EACH YEAR AND SEX, so every year carries its
#      own composition rather than a single profile imposed on all of them.
#   4. Each year's smoothed stock is normalised to a profile and scaled to the
#      reconciled total for that year from migration_bounds in 00_setup.R.
#      Eurostat gives the shape; it cannot give the size, because it covers
#      the EU only and its stock is flat after 2022 for administrative
#      reasons rather than demographic ones.
#
# The output is the annual NET FLOW by single year of age and sex, negative
# for people leaving, which is what 11 and 13 consume.
#
# The CES survey is read below and plotted against the Eurostat profile as a
# cross-check. It is NOT used to build the output: its own age-sex structure
# describes refugees outside Russia and Belarus only (the survey excluded
# them by design), while its value here is the annual totals in its table 1,
# which enter through migration_bounds.
#
# INPUT    data_input/refugees_eurostat/migr_asytpsm_*.xlsx
# OUTPUT   data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds
#          <- used by 11 and 13
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

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

# ============================================================================
# READ CES AGE-SEX FRACTIONATION (observed refugee demographics from surveys)
# ============================================================================
# CES data is given as % of stock within each broad age group, by sex
# We use this to distribute within Eurostat bands 0-14, 14-18, 18-35, 35-65
# For ages 65+, we fall back to pclm smoothing (limited CES data quality there)

ces_raw <- read_csv("data_input/ukr_centre_for_economic_strategy_age_sex_migrants.csv",
                    col_types = cols(.default = "c")) %>%
  rename(grp = 1) %>%
  mutate(m = abs(as.numeric(str_replace(Men, ",", "."))),
         f = as.numeric(str_replace(Women, ",", "."))) %>%
  select(grp, m, f)

# Map CES groups to age ranges
ces_raw <- ces_raw %>%
  mutate(lo = c(0, 6, 10, 14, 18, 25, 35, 45, 55, 65, 75),
         hi = c(5, 9, 13, 17, 24, 34, 44, 54, 64, 74, 89))

# Expand to single ages and normalize within each group
ces_shape <- ces_raw %>%
  pivot_longer(c(m, f), names_to = "sex", values_to = "pct") %>%
  mutate(age = map2(lo, hi, seq)) %>%
  unnest(age) %>%
  group_by(grp, sex) %>%
  mutate(dens = pct / (hi - lo + 1)) %>%
  ungroup() %>%
  select(sex, age, dens) %>%
  complete(sex, age = 0:100, fill = list(dens = 0))

# Ungrouping to single ages using pclm
# NOTE: CES age-sex fractionation (ces_shape, loaded above) is available for future
# refinement of within-band distributions, but currently using pclm for simplicity.
ung_age_mig <- function(chunk) {
  # chunk has columns: age, mix, and grouping info (year, sex)
  year_val <- chunk$year[1]
  sex_val <- chunk$sex[1]

  dt_in <- tibble(age = chunk$age, mix = chunk$mix, mix_mt = mix * 1e5) %>%
    mutate(mix_mt = ifelse(mix_mt == 0, 1, mix_mt))
  suppressWarnings({
    mix <- pclm(x = dt_in$age, y = dt_in$mix_mt, nlast = 36)$fitted
  })
  tibble(year = year_val, sex = sex_val, age = 0:100, mix = mix / 1e5)
}

dt4 <-
  dt3 %>%
  group_by(year, sex) %>%
  do(ung_age_mig(.)) %>%
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
# Eurostat supplies the SHAPE, 00_setup.R supplies the SIZE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The registrations describe the age-sex composition of the outflow well, and
# they describe its size badly: they cover the EU only, they exclude Russia
# and Belarus entirely, and from 2023 the stock is flat because purges of
# inactive status offset new grants rather than because migration stopped.
# The totals therefore come from migration_bounds, which reconciles the
# western and Russia/Belarus components; see 00_setup.R for the sources.
#
# The profile is the year's smoothed STOCK, not the year-on-year difference
# in it. The difference cannot serve: in 2024 it is the near-cancellation of
# two 4.3M stocks - about 1,200 people - whose age-sex shape is noise, and
# scaling that up to a quarter of a million would amplify it 200-fold.
# Taking the profile from the stock keeps what 06 is good for, which is that
# each year is smoothed on its own and so carries its own composition: the
# share of men rises over the period and the share above 45 falls, and both
# survive into the output.
mig_tot <-
  migration_bounds |>
  summarise(tot = sum(mode), .by = year)

dt5 <-
  dt4 |>
  mutate(cx = mix / sum(mix), .by = year) |>
  left_join(mig_tot, by = "year") |>
  # negative = people leaving, the convention 11 and 13 read
  mutate(mix = -round(cx * tot, 0)) |>
  select(sex, age, year, mix) |>
  arrange(year, sex, age)

# the profile must not have changed the totals it was scaled to
stopifnot(
  all.equal(
    dt5 |> summarise(mix = -sum(mix), .by = year) |> arrange(year) |> pull(mix),
    mig_tot |> arrange(year) |> pull(tot),
    tolerance = 1e-4
  )
)

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
