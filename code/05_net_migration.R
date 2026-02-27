rm(list = ls())
gc()
source("code/00_setup.R")

# load refugee data
out_refugee_file <- 'data_input/unhcr_refugees_out_daily.csv'
in_refugee_file <- 'data_input/unhcr_refugees_in_daily.csv'

# in_refugee_file <- 'data_inter/data_pdr/data/refugees_in_daily.csv'
# refugee_file <- 'data_inter/data_pdr/data/refugees_daily.csv'

# download border crossing data
download.file(
  url = 'https://data.unhcr.org/population/get/timeseries?export=csv&widget_id=324084&sv_id=54&population_group=5460&frequency=day&fromDate=1900-01-01',
  destfile = out_refugee_file
)

download.file(
  url = 'https://data.unhcr.org/population/get/timeseries?export=csv&widget_id=324085&sv_id=54&population_group=5472&frequency=day&fromDate=1900-01-01',
  destfile = in_refugee_file
)

# load border crossing data
dat_refugees_out <-
  read.csv(
    out_refugee_file,
    stringsAsFactors = F,
    skip = 1,
    skipNul = T,
    sep = ';'
  )

# |>
#   rename(date = 1, individuals = 2) |>
#   as_tibble()

dat_refugees_in <-
  read.csv(
    in_refugee_file,
    stringsAsFactors = F,
    skip = 1,
    skipNul = T,
    sep = ';'
  )

#---- net border crossings ----#

# calculate
dat_refugees <-
  dat_refugees_out |>
  rename(date = 1, outs = 2) |>
  as_tibble() |>
  left_join(
    dat_refugees_in |>
      rename(date = 1, ins = 2)
  ) %>%
  mutate(net = ins - outs, date = ymd(date)) |>
  drop_na(date)

# 18% subregistration according to PDR
net_day <- dat_refugees |>
  mutate(net2 = net * 1.18)


dates2 <- c("2022-12-31", "2023-12-31", "2024-12-31")

emi_unhcr <-
  net_day %>%
  filter(date %in% dates2) %>%
  mutate(ems = -net, year = year(date)) %>%
  select(year, ems) %>%
  spread(year, ems) %>%
  mutate(em_2023 = `2023` - `2022`, em_2024 = `2024` - `2023`) %>%
  select(`2022`, `2023` = em_2023, `2024` = em_2024) %>%
  gather(`2022`, `2023`, `2024`, key = year, value = mix)

# loading data from EUROSTAT

# Ukranian refuggies in EU countries

mths <- format(
  seq(as.Date("2022-03-01"), as.Date("2025-12-01"), by = "month"),
  "%Y-%m"
)

shts <- tibble(
  sht = 1:24,
  sex = c(rep("t", 8), rep("m", 8), rep("f", 8)),
  age = c(rep(
    c(
      "Total",
      "Less than 14 years",
      "From 14 to 17 years",
      "Less than 18 years",
      "From 18 to 34 years",
      "From 35 to 64 years",
      "65 years or over",
      "Unknown"
    ),
    3
  ))
)

# loading data from Eurostat
mixs <- tibble()
for (i in 1:24) {
  sx <-
    shts %>%
    filter(sht == i) %>%
    pull(sex)

  ag <-
    shts %>%
    filter(sht == i) %>%
    pull(age)

  tmp <-
    read_xlsx(
      "data_input/refugees_eurostat/migr_asytpsm__sex_age_month.xlsx",
      sheet = paste0("Sheet ", i),
      skip = 10
    ) %>%
    rename(state = 1) %>%
    gather(-state, key = month, value = mix) %>%
    filter(
      month %in% mths,
      !state %in% c("European Union - 27 countries (from 2020)", "GEO (Labels)")
    ) %>%
    mutate(mix = mix %>% as.double()) %>%
    drop_na(mix) %>%
    mutate(sex = sx, age = ag)

  mixs <-
    bind_rows(mixs, tmp)
}

unique(mixs$age)
unique(mixs$state)

ags2 <- c(
  "Less than 14 years",
  "From 14 to 17 years",
  "From 18 to 34 years",
  "From 35 to 64 years",
  "65 years or over"
)

mixs2 <-
  mixs %>%
  filter(str_sub(month, 6, 7) == 12) %>%
  summarise(mix = sum(mix), .by = c(month, sex, age))

mixs3 <-
  mixs2 %>%
  filter(age %in% ags2, str_sub(month, 6, 7) == 12) %>%
  mutate(
    age = case_when(
      age == "Less than 14 years" ~ 0,
      age == "From 14 to 17 years" ~ 14,
      age == "From 18 to 34 years" ~ 18,
      age == "From 35 to 64 years" ~ 35,
      age == "65 years or over" ~ 65
    )
  ) %>%
  group_by(month, sex) %>%
  mutate(cx = mix / sum(mix)) %>%
  ungroup() %>%
  left_join(
    mixs2 %>%
      filter(age == "Total") %>%
      select(-age) %>%
      rename(mix_tot = mix)
  ) %>%
  mutate(mix = cx * mix_tot) %>%
  select(month, sex, age, mix)

mixs4 <-
  mixs3 %>%
  filter(sex != "t") %>%
  group_by(month, age) %>%
  mutate(cx = mix / sum(mix)) %>%
  ungroup() %>%
  left_join(
    mixs3 %>%
      filter(sex == "t") %>%
      select(-sex) %>%
      rename(mix_tot = mix)
  ) %>%
  mutate(mix = cx * mix_tot) %>%
  mutate(date = ymd(paste0(month, "-15"))) %>%
  select(date, sex, age, mix) %>%
  mutate(year = year(date)) %>%
  summarise(mix = sum(mix), .by = c(year, sex, age))

mixs5 <-
  mixs4 %>%
  group_by(year) %>%
  mutate(migs = sum(mix)) |>
  ungroup()

mixs4 %>%
  filter(age != "Total", sex != "t") %>%
  summarise(mix = sum(mix), .by = c(year))

mixs %>%
  mutate(year = str_sub(month, 1, 4)) %>%
  filter(age == "Total", sex != "t", str_sub(month, 6, 7) == 12) %>%
  summarise(mix = sum(mix), .by = c(year))

mixs5 |>
  filter(age != "Total", sex != "t") %>%
  summarise(mix = sum(mix), .by = c(year))


unique(mixs4$age)

mixs5 <-
  mixs4 %>%
  spread(year, mix) %>%
  mutate(em_2023 = `2023` - `2022`, em_2024 = `2024` - `2023`) %>%
  select(sex, age, `2022`, `2023` = em_2023, `2024` = em_2024) %>%
  gather(`2022`, `2023`, `2024`, key = year, value = mix) %>%
  mutate(ems = ifelse(mix > 0, mix, 0), ims = ifelse(mix < 0, -mix, 0))

chunk <-
  mixs5 %>%
  filter(year == 2022, sex == "m")

ung_age_mig <- function(chunk) {
  dt_in <-
    tibble(age = chunk$age, ems = chunk$ems, ems_mt = ems * 1e5) %>%
    mutate(ems_mt = ifelse(ems_mt == 0, 1, ems_mt))
  nl <- 36
  ems <- pclm(x = dt_in$age, y = dt_in$ems_mt, nlast = nl)$fitted

  fit <- tibble(age = 0:100, ems = round(ems / 1e5))

  out <-
    chunk %>%
    select(year, sex) %>%
    unique() %>%
    left_join(fit, by = character())
  return(out)
}

mixs6_ems <-
  mixs5 %>%
  group_by(year, sex) %>%
  do(ung_age_mig(chunk = .data)) %>%
  ungroup()

mixs6_ims <-
  mixs5 %>%
  select(sex, age, year, ems = ims) %>%
  group_by(year, sex) %>%
  do(ung_age_mig(chunk = .data)) %>%
  ungroup() %>%
  rename(ims = ems)


# combining UNHCR and EUROSTAT data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mixs6 <-
  mixs6_ems %>%
  left_join(mixs6_ims) %>%
  mutate(mix = ims - ems) %>%
  group_by(year) %>%
  mutate(mix_sum = sum(-mix)) %>%
  ungroup() %>%
  left_join(
    emi_unhcr %>%
      rename(mix_tot = mix)
  ) %>%
  mutate(mix2 = round(mix * (mix_tot / mix_sum)))

mixs_out <-
  mixs6 %>%
  select(year, sex, age, mix = mix2) %>%
  mutate(year = as.double(year))

mixs_out %>%
  summarise(mix = sum(mix), .by = c(year))

write_rds(
  mixs_out,
  "data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2024.rds"
)
