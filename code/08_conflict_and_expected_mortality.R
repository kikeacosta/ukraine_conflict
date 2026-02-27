rm(list = ls())
gc()
source("code/00_setup.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# loading the data ====
# ~~~~~~~~~~~~~~~~~~~~~~~
# UN crisis mortality patterns
as <- read_rds("data_input/un_mort_crises_age_str.rds")

# conflict deaths scenarios
ucdp <- read_rds("data_inter/ukr_ucdp.rds") %>% select(-country, -code)
acled <- read_rds("data_inter/ukr_acled.rds")
ual <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")

# excluding prisioners from UALOSSES and 2025
ual2 <-
  ual %>%
  filter(status != "prisioner", year %in% 2022:2024) %>%
  summarise(cnf = sum(dx), .by = c(year, sex, age)) |>
  arrange(year, sex, age)

# expected (nonconflict) mortality
exp_mort <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds")

# exposures
# pop <- read_rds("data_inter/ukr_pop_wpp_age_single_january.rds")
pop <- read_rds("data_inter/ukr_pop_sssu.rds")

# net migration
migs <- read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2024.rds")

# age-specific fertility rates
asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2024.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# analyzing the data ====
# ~~~~~~~~~~~~~~~~~~~~~~~

# comparing conflict mortality by source / scenario
ucdp %>%
  filter(year %in% 2022:2024) %>%
  summarise(ucdp = sum(dts), .by = year) |>
  left_join(
    acled %>%
      filter(year %in% 2022:2024) %>%
      summarise(acled = sum(dts), .by = year),
    by = "year"
  ) |>
  left_join(
    bind_rows(
      ual2 %>%
        summarise(dts = sum(cnf), .by = year),
      ucdp %>%
        filter(year %in% 2022:2024, role == "civilians")
    ) |>
      summarise(uals = sum(dts), .by = year)
  )

# adding the same fertility to 2024 as in 2023
asfr2 <-
  asfr %>%
  bind_rows(
    asfr %>%
      filter(year == 2023) %>%
      mutate(year = 2024)
  )

# net migration
# for now, not immigration, only emigration
migs2 <-
  migs %>%
  mutate(ems = ifelse(mix < 0, -mix, 0)) %>%
  select(-mix) %>%
  mutate(ems = ifelse(year == 2024, 0, ems))

# baseline population (January 1st, 2022)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop22_ini <-
  pop %>%
  filter(reg == "cnt") |>
  select(-reg) |>
  filter(year == 2022)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scenario 1: UCDP civilian and combatant deaths ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2022 ====
# ~~~~~~~~~
exposures_22 <- est_exposures(pop22_ini)
cnf_y22 <- est_age_sex_conf(exposures_22, ucdp, "all")
# population at the end of the year and expected mortality
pop22_end <- est_pop_end(exposures_22, cnf_y22)
sum_y22 <- est_summary(pop22_end)

# 2023 ====
# ~~~~~~~~~
pop23_ini <- end_to_ini(pop22_end)
exposures_23 <- est_exposures(pop23_ini)
cnf_y23 <- est_age_sex_conf(exposures_23, ucdp, "all")
# population at the end of the year and expected mortality
pop23_end <- est_pop_end(exposures_23, cnf_y23)
sum_y23 <- est_summary(pop23_end)

# 2024 ====
# ~~~~~~~~~
pop24_ini <- end_to_ini(pop23_end)
exposures_24 <- est_exposures(pop24_ini)
cnf_y24 <-
  est_age_sex_conf(exposures_24, ucdp, "all")
# population at the end of the year and expected mortality
pop24_end <- est_pop_end(exposures_24, cnf_y24)
sum_y24 <- est_summary(pop24_end)

# all together ====
# ~~~~~~~~~~~~~~~~~

ys_ucdp <-
  bind_rows(
    sum_y22,
    sum_y23,
    sum_y24
  ) %>%
  mutate(sce = "ucdp")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scenario 2: ACLED civilian and combatant deaths ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2022 ====
# ~~~~~~~~~
exposures_22 <- est_exposures(pop22_ini)
cnf_y22 <- est_age_sex_conf(exposures_22, acled, "all")
# population at the end of the year and expected mortality
pop22_end <- est_pop_end(exposures_22, cnf_y22)
sum_y22 <- est_summary(pop22_end)

# 2023 ====
# ~~~~~~~~~
pop23_ini <- end_to_ini(pop22_end)
exposures_23 <- est_exposures(pop23_ini)
cnf_y23 <- est_age_sex_conf(exposures_23, acled, "all")
# population at the end of the year and expected mortality
pop23_end <- est_pop_end(exposures_23, cnf_y23)
sum_y23 <- est_summary(pop23_end)

# 2024 ====
# ~~~~~~~~~
pop24_ini <- end_to_ini(pop23_end)
exposures_24 <- est_exposures(pop24_ini)
cnf_y24 <-
  est_age_sex_conf(exposures_24, acled, "all")
# population at the end of the year and expected mortality
pop24_end <- est_pop_end(exposures_24, cnf_y24)
sum_y24 <- est_summary(pop24_end)

# all together ====
# ~~~~~~~~~~~~~~~~~

ys_acled <-
  bind_rows(
    sum_y22,
    sum_y23,
    sum_y24
  ) %>%
  mutate(sce = "acled")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scenario UALOSSES combatant and UCDP civilian deaths ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exposures_22_sc3 <- exposures_22

# ungrouping civilians from UCDP
cvs_y22 <- est_age_sex_conf(exposures_22_sc3, ucdp, "civilians")
cnf_y22_sc3 <-
  ual2 |>
  filter(year == 2022) |>
  bind_rows(cvs_y22) |>
  summarise(cnf = sum(cnf), .by = c(year, sex, age))

# population at the end of the year and expected mortality
pop22_end_sc3 <- est_pop_end(exposures_22_sc3, cnf_y22_sc3)
sum_y22_sc3 <- est_summary(pop22_end_sc3)

# 2023 ====
# ~~~~~~~~~
pop23_ini_sc3 <- end_to_ini(pop22_end_sc3)
exposures_23_sc3 <- est_exposures(pop23_ini_sc3)
# ungrouping civilians from UCDP
cvs_y23 <- est_age_sex_conf(exposures_23_sc3, ucdp, "civilians")
cnf_y23_sc3 <-
  ual2 |>
  filter(year == 2023) |>
  bind_rows(cvs_y23) |>
  summarise(cnf = sum(cnf), .by = c(year, sex, age))
# population at the end of the year and expected mortality
pop23_end_sc3 <- est_pop_end(exposures_23_sc3, cnf_y23_sc3)
sum_y23_sc3 <- est_summary(pop23_end_sc3)

# 2024 ====
# ~~~~~~~~~
pop24_ini_sc3 <- end_to_ini(pop23_end_sc3)
exposures_24_sc3 <- est_exposures(pop24_ini_sc3)
# ungrouping civilians from UCDP
cvs_y24 <- est_age_sex_conf(exposures_24_sc3, ucdp, "civilians")
# adding combatant deaths from UALOSSES
cnf_y24_sc3 <-
  ual2 |>
  filter(year == 2024) |>
  bind_rows(cvs_y24) |>
  summarise(cnf = sum(cnf), .by = c(year, sex, age))
# population at the end of the year and expected mortality
pop24_end_sc3 <- est_pop_end(exposures_24_sc3, cnf_y24_sc3)
sum_y24_sc3 <- est_summary(pop24_end_sc3)

# merging the three years together
ys_ual <-
  bind_rows(
    sum_y22_sc3,
    sum_y23_sc3,
    sum_y24_sc3
  ) %>%
  mutate(sce = "ualosses_ucdp")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merging all scenarios together and plotting ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out1 <-
  bind_rows(
    ys_ucdp,
    ys_acled,
    ys_ual
  )

out2 <-
  out1 %>%
  bind_rows(
    out1 %>%
      summarise(
        dx = sum(dx),
        pop = sum(pop),
        .by = c(sce, year, age, cause)
      ) %>%
      mutate(sex = "t", mx = dx / pop)
  )

unique(out2$cause)
unique(out2$sex)
table(out2$cause)
table(out2$sex)
table(out2$sce)
table(out2$year, out2$age)
table(out2$sce, out2$sex)

write_rds(
  out2,
  "data_inter/ukr_expected_and_conflict_deaths_2022_2024.rds"
)

out2 %>%
  ggplot() +
  geom_line(aes(age, mx, col = sex)) +
  facet_grid(cause + year ~ sce) +
  scale_y_log10() +
  theme_bw()

out2 %>%
  ggplot() +
  geom_line(aes(age, dx, col = sex)) +
  facet_grid(sce + cause ~ year) +
  # scale_y_log10()+
  theme_bw()

out2 %>%
  filter(cause != "conflict", sex != "t") %>%
  ggplot() +
  geom_line(aes(age, mx, col = cause)) +
  facet_grid(sce ~ sex + year) +
  scale_y_log10() +
  theme_bw()
