rm(list = ls())
gc()
source("code/00_setup.R")

# population estimates
pop1 <- fread(file.path(
  "data_input",
  "WPP2024_Population1JanuaryBySingleAgeSex_Medium_1950-2023.csv.gz"
))
pop2 <- fread(file.path(
  "data_input",
  "WPP2024_Population1JanuaryBySingleAgeSex_Medium_2024-2100.csv.gz"
))

pop <-
  bind_rows(pop1, pop2) %>%
  filter(Location == "Ukraine") %>%
  select(year = Time, age = AgeGrpStart, m = PopMale, f = PopFemale) %>%
  gather("m", "f", key = sex, value = pop) %>%
  mutate(pop = 1e3 * pop) %>%
  as_tibble()

write_rds(pop, "data_inter/ukr_pop_wpp_age_single_january.rds")


asfr_y1 <-
  read_xlsx(
    "data_input/WPP2024_FERT_F01_FERTILITY_RATES_BY_SINGLE_AGE_OF_MOTHER.xlsx",
    skip = 16
  ) %>%
  filter(`Location code` == 804, Year %in% 2022:2024) %>%
  select(11:46) %>%
  gather(-Year, key = age, value = fx) %>%
  mutate(fx = as.double(fx) / 1000, sex = "f", age = age %>% as.double()) %>%
  rename(year = Year)

write_rds(asfr_y1, "data_inter/ukr_asfr_wpp_2022_2024.rds")
