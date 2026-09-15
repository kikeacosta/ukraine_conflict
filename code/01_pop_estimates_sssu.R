# ==============================================================================
# STEP 01 - Population and death counts from the State Statistics Service
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Reads the SSSU workbook (population and deaths by single year of age, sex
# and region) and reshapes it into tidy long form. Three regions are kept
# separately:
#     cnt = continental Ukraine, i.e. the whole country WITHOUT Crimea and
#           Sevastopol (sheets "*Continental"). It INCLUDES Donetsk and Luhansk.
#     dnk = Donetsk region   - a subset of cnt
#     luk = Luhansk region   - a subset of cnt
# dnk and luk are shipped separately because their registration is incomplete:
# since 2015 only their government-controlled parts register events, but SSSU
# kept estimating their population as if coverage were complete (see the
# workbook's ReadMe sheet). Because they are nested inside cnt, rows must never
# be summed across regions.
# The projection in step 11 starts from "cnt" only.
#
# The population is what the pipeline consumes. The death counts are tidied
# here too, but no later step reads them: they previously fed a comparison of
# SSSU death rates against an independent life-table source (LifeTables.xlsx),
# written to ukr_dts_pop_sssu_grig.rds. Nothing ever consumed that output, so
# the comparison and its input have been retired.
#
# INPUT    data_input/sssu_ukr_data.xlsx
# OUTPUT   data_inter/ukr_pop_sssu.rds   <- used by 11 and 13
# ==============================================================================

rm(list = ls())
source("code/00_setup.R")

# loading death counts and population estimates from SSSU
sssu_file <- "data_input/sssu_ukr_data.xlsx"

dts_cnt <- read_xlsx(sssu_file, sheet = "DeathsContinental")
dts_dnk <- read_xlsx(sssu_file, sheet = "DeathsDonetskReg")
dts_luk <- read_xlsx(sssu_file, sheet = "DeathsLuhanskReg")
pop_cnt <- read_xlsx(sssu_file, sheet = "PopContinental")
pop_dnk <- read_xlsx(sssu_file, sheet = "PopDonetskReg")
pop_luk <- read_xlsx(sssu_file, sheet = "PopLuhanskReg")

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

# NOTE: the deaths sheets above are read and tidied into dts_all, but nothing
# downstream consumes them. They previously fed a comparison of SSSU death
# rates against an independent life-table source (LifeTables.xlsx), written
# to ukr_dts_pop_sssu_grig.rds. That output was never read by any later step,
# so the comparison has been retired along with it.
