rm(list = ls())
source("code/00_setup.R")

# loading death counts and population estimates from SSSU
dts_cnt <- read_xlsx(
  "data_input/sssu_ukr_data.xlsx",
  sheet = "DeathsContinental"
)
dts_dnk <- read_xlsx(
  "data_input/sssu_ukr_data.xlsx",
  sheet = "DeathsDonetskReg"
)
dts_luk <- read_xlsx(
  "data_input/sssu_ukr_data.xlsx",
  sheet = "DeathsLuhanskReg"
)
pop_cnt <- read_xlsx("data_input/sssu_ukr_data.xlsx", sheet = "PopContinental")
pop_dnk <- read_xlsx("data_input/sssu_ukr_data.xlsx", sheet = "PopDonetskReg")
pop_luk <- read_xlsx("data_input/sssu_ukr_data.xlsx", sheet = "PopLuhanskReg")

# putting all together
# deaths
dts_cnt2 <-
  dts_cnt %>%
  rename(age = 1) %>%
  mutate(
    sex = ifelse(is.na(as.numeric(age)), str_to_lower(str_sub(age, 1, 1)), NA)
  ) %>%
  fill(sex) %>%
  drop_na() %>%
  filter(age != "Males", age != "Females") %>%
  gather(-sex, -age, key = year, value = dts) %>%
  mutate(reg = "cnt", age = age %>% as.integer(), year = year %>% as.integer())

dts_dnk2 <-
  dts_dnk %>%
  rename(age = 1) %>%
  mutate(
    sex = ifelse(is.na(as.numeric(age)), str_to_lower(str_sub(age, 1, 1)), NA)
  ) %>%
  fill(sex) %>%
  drop_na() %>%
  filter(age != "Males", age != "Females") %>%
  gather(-sex, -age, key = year, value = dts) %>%
  mutate(reg = "dnk", age = age %>% as.integer(), year = year %>% as.integer())

dts_luk2 <-
  dts_luk %>%
  rename(age = 1) %>%
  mutate(
    sex = ifelse(is.na(as.numeric(age)), str_to_lower(str_sub(age, 1, 1)), NA)
  ) %>%
  fill(sex) %>%
  drop_na() %>%
  filter(age != "Males", age != "Females") %>%
  gather(-sex, -age, key = year, value = dts) %>%
  mutate(reg = "luk", age = age %>% as.integer(), year = year %>% as.integer())

dts_all <-
  bind_rows(dts_cnt2, dts_dnk2, dts_luk2)

# population
pop_cnt2 <-
  pop_cnt %>%
  rename(age = 1) %>%
  mutate(
    sex = ifelse(is.na(as.numeric(age)), str_to_lower(str_sub(age, 1, 1)), NA)
  ) %>%
  fill(sex) %>%
  drop_na() %>%
  filter(age != "Males", age != "Females") %>%
  gather(-sex, -age, key = year, value = pop) %>%
  mutate(reg = "cnt", age = age %>% as.integer(), year = year %>% as.integer())

pop_dnk2 <-
  pop_dnk %>%
  rename(age = 1) %>%
  mutate(
    sex = ifelse(is.na(as.numeric(age)), str_to_lower(str_sub(age, 1, 1)), NA)
  ) %>%
  fill(sex) %>%
  drop_na() %>%
  filter(age != "Males", age != "Females") %>%
  gather(-sex, -age, key = year, value = pop) %>%
  mutate(reg = "dnk", age = age %>% as.integer(), year = year %>% as.integer())

pop_luk2 <-
  pop_luk %>%
  rename(age = 1) %>%
  mutate(
    sex = ifelse(is.na(as.numeric(age)), str_to_lower(str_sub(age, 1, 1)), NA)
  ) %>%
  fill(sex) %>%
  drop_na() %>%
  filter(age != "Males", age != "Females") %>%
  gather(-sex, -age, key = year, value = pop) %>%
  mutate(reg = "luk", age = age %>% as.integer(), year = year %>% as.integer())

pop_all <-
  bind_rows(pop_cnt2, pop_dnk2, pop_luk2)

write_rds(pop_all, "data_inter/ukr_pop_sssu.rds")

# deaths + population
dt <-
  dts_all %>%
  left_join(pop_all) %>%
  filter(age != "tot") %>%
  mutate(age = age %>% as.integer(), year = year %>% as.integer())


# ~~~~~~~~~~~~~~~~~~~~~
# from life tables ====
# ~~~~~~~~~~~~~~~~~~~~~

exs <-
  read_xlsx("data_input/LifeTables.xlsx", sheet = "LT")

unique(exs$Sex)

mxs <-
  exs %>%
  filter(LT_function == "qx") %>%
  select(-LT_function) %>%
  rename_with(tolower) %>%
  gather(-sex, -year, key = age, value = qx) %>%
  mutate(
    age = str_remove(age, "age_") %>% as.integer(),
    year = year %>% as.integer(),
    sex = str_sub(sex, 1, 1),
    reg = "cnt"
  ) %>%
  filter(year >= 2002) %>%
  mutate(
    ax = case_when(age == 0 ~ 0.13, age == 100 ~ 1.25, TRUE ~ 0.5),
    mx = qx / (1 - (1 - ax) * qx)
  ) %>%
  left_join(
    dt %>%
      select(-dts)
  ) %>%
  mutate(dts = mx * pop, source = "grig") %>%
  select(-ax, -qx)

dt2 <-
  dt %>%
  mutate(mx = dts / pop, source = "sssu") %>%
  bind_rows(mxs) %>%
  rename(region = reg)

write_rds(dt2, "data_inter/ukr_dts_pop_sssu_grig.rds")

# # life tables
# ltb_f1 <- fread(file.path(
#   "data_input",
#   "WPP2024_Life_Table_Complete_Medium_Female_1950-2023.csv.gz"
# ))
# ltb_f2 <- fread(file.path(
#   "data_input",
#   "WPP2024_Life_Table_Complete_Medium_Female_2024-2100.csv.gz"
# ))
# ltb_m1 <- fread(file.path(
#   "data_input",
#   "WPP2024_Life_Table_Complete_Medium_Male_1950-2023.csv.gz"
# ))
# ltb_m2 <- fread(file.path(
#   "data_input",
#   "WPP2024_Life_Table_Complete_Medium_Male_2024-2100.csv.gz"
# ))

# # population estimates
# pop1 <- fread(file.path(
#   "data_input",
#   "WPP2024_PopulationBySingleAgeSex_Medium_1950-2023.csv.gz"
# ))
# pop2 <- fread(file.path(
#   "data_input",
#   "WPP2024_PopulationBySingleAgeSex_Medium_2024-2100.csv.gz"
# ))
