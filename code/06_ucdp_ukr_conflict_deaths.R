rm(list = ls())
source("code/00_setup.R")

# UCDP Georeferenced Event Dataset (GED)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1989 - 2022
dt5 <- read_csv("data_input/ucdp/GEDEvent_v25_1.csv")
# january - June 2023
dt6 <- read_csv("data_input/ucdp/GEDEvent_v25_01_25_06.csv")
# September 2023
dt7 <- read_csv("data_input/ucdp/GEDEvent_v25_0_7.csv")

in1 <-
  dt5 %>%
  select(
    year,
    type_of_violence,
    conflict_name,
    side_a,
    side_b,
    country,
    deaths_a,
    deaths_b,
    deaths_civilians,
    deaths_unknown,
    best
  )

in2 <-
  dt6 %>%
  select(
    year,
    type_of_violence,
    conflict_name,
    side_a,
    side_b,
    country,
    deaths_a,
    deaths_b,
    deaths_civilians,
    deaths_unknown,
    best
  )

in3 <-
  dt7 %>%
  select(
    year,
    type_of_violence,
    conflict_name,
    side_a,
    side_b,
    country,
    deaths_a,
    deaths_b,
    deaths_civilians,
    deaths_unknown,
    best
  )

all <-
  bind_rows(
    in1,
    in2,
    in3
  )

# imputing unknown status (comb a, comb b, or civilians)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all2 <-
  all %>%
  select(
    year,
    country,
    type_of_violence,
    side_a,
    side_b,
    a = deaths_a,
    b = deaths_b,
    c = deaths_civilians,
    u = deaths_unknown,
    t = best
  ) %>%
  replace_na(list(a = 0, b = 0, c = 0, t = 0)) %>%
  filter(t > 0)


# identifying civilians and foreign and local soldiers in inter-state conflicts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inter_state <-
  all2 %>%
  filter(type_of_violence == 1 & str_detect(side_b, "Government"))

inter_state2 <-
  inter_state %>%
  group_by(country, side_a, side_b, year) %>%
  summarise(a = sum(a), b = sum(b), c = sum(c), u = sum(u), t = sum(t)) %>%
  ungroup() %>%
  group_by(country, side_a, side_b) %>%
  mutate(
    a2 = case_when(
      u == 0 ~ a,
      u != 0 & a != 0 & b != 0 & c != 0 ~ a + u * (a / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ a +
        u * sum(a) / (sum(a) + sum(b) + sum(c))
    ),
    b2 = case_when(
      u == 0 ~ b,
      u != 0 & a != 0 & b != 0 & c != 0 ~ b + u * (b / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ b +
        u * sum(b) / (sum(a) + sum(b) + sum(c))
    ),
    c2 = case_when(
      u == 0 ~ c,
      u != 0 & a != 0 & b != 0 & c != 0 ~ c + u * (c / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ c +
        u * sum(c) / (sum(a) + sum(b) + sum(c))
    ),
    t2 = a2 + b2 + c2,
    diff = t - t2
  ) %>%
  ungroup() %>%
  select(country, side_a, side_b, year, a = a2, b = b2, c = c2, t)

# identifying civilians and foreign and local soldiers in intra-state conflicts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
intra_state <-
  all2 %>%
  anti_join(inter_state) %>%
  group_by(country, year) %>%
  summarise(a = sum(a), b = sum(b), c = sum(c), u = sum(u), t = sum(t)) %>%
  # group_by(country) %>%
  # filter(sum(t) > 1000) %>%
  ungroup() %>%
  group_by(country) %>%
  mutate(
    a2 = case_when(
      u == 0 ~ a,
      u != 0 & a != 0 & b != 0 & c != 0 ~ a + u * (a / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ a +
        u * sum(a) / (sum(a) + sum(b) + sum(c))
    ),
    b2 = case_when(
      u == 0 ~ b,
      u != 0 & a != 0 & b != 0 & c != 0 ~ b + u * (b / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ b +
        u * sum(b) / (sum(a) + sum(b) + sum(c))
    ),
    c2 = case_when(
      u == 0 ~ c,
      u != 0 & a != 0 & b != 0 & c != 0 ~ c + u * (c / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ c +
        u * sum(c) / (sum(a) + sum(b) + sum(c))
    ),
    t2 = a2 + b2 + c2,
    diff = t - t2
  )

intra_state2 <-
  intra_state %>%
  select(country, year, a = a2, b = b2, c = c2, t)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identifying combatants nationalities and civilians
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combatants killed in inter_state conflict
comb_inter <-
  inter_state2 %>%
  mutate(
    country_a = str_replace(side_a, "Government of ", ""),
    country_b = str_replace(side_b, "Government of ", "")
  ) %>%
  group_by(country_a, country_b, year) %>%
  summarise(dts_a = sum(a), dts_b = sum(b)) %>%
  ungroup()

comb_inter2 <-
  bind_rows(
    comb_inter %>%
      select(country = country_a, year, dts = dts_a),
    comb_inter %>%
      select(country = country_b, year, dts = dts_b)
  ) %>%
  mutate(role = "combatants")

# combatants killed in intra_state conflict
comb_intra <-
  intra_state2 %>%
  group_by(country, year) %>%
  summarise(dts = sum(a + b)) %>%
  ungroup() %>%
  mutate(role = "combatants")

civils <-
  bind_rows(inter_state2, intra_state2) %>%
  group_by(country, year) %>%
  summarise(dts = sum(c)) %>%
  ungroup() %>%
  mutate(role = "civilians")

# testing consistency before and after imputation ~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dts <-
  bind_rows(comb_inter2, comb_intra, civils) %>%
  group_by(country, year, role) %>%
  summarise(dts = sum(dts)) %>%
  ungroup() %>%
  mutate(
    country = ifelse(country == "Yemen (North Yemen)", "Yemen", country),
    code = countrycode(country, origin = "country.name", destination = "iso3c")
  ) %>%
  drop_na(dts)

dts %>%
  summarise(dts = sum(dts))

all2 %>%
  summarise(dts = sum(t))

ukr <-
  dts %>%
  filter(country == "Ukraine", year >= 2022)

copy_this(
  ukr %>%
    summarise(dts = sum(dts), .by = c(year, role))
)

write_rds(ukr, "data_inter/ukr_ucdp.rds")
