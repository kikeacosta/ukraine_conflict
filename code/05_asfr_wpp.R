rm(list = ls())
gc()
source("code/00_setup.R")

# Age-specific fertility rates for Ukraine, WPP 2024.
#
# The source workbook is ~78 MB and is not tracked in git. The extract below
# is a few kB, so once it has been built once the raw file is never needed
# again (see cache_parquet() in 00_setup.R).
#
# NOTE: WPP2024 only publishes observed fertility up to 2023. Scripts 11 and
# 13 carry the 2023 schedule forward to 2024 and 2025.
asfr_y1 <- cache_parquet("data_inter/ukr_asfr_wpp_2022_2025.parquet", {
  read_xlsx(
    require_raw(paste0(
      "data_input/",
      "WPP2024_FERT_F01_FERTILITY_RATES_BY_SINGLE_AGE_OF_MOTHER.xlsx"
    )),
    skip = 16
  ) %>%
    filter(`Location code` == 804, Year %in% 2022:2025) %>%
    select(11:46) %>%
    gather(-Year, key = age, value = fx) %>%
    mutate(fx = as.double(fx) / 1000, sex = "f", age = age %>% as.double()) %>%
    rename(year = Year)
})

# the hard-coded column range above is fragile, so check what came out
stopifnot(
  all(asfr_y1$age %in% 15:49),
  nrow(asfr_y1) == length(unique(asfr_y1$year)) * 35,
  all(asfr_y1$fx >= 0 & asfr_y1$fx < 0.5)
)

asfr_y1 %>%
  summarise(TFR = sum(fx), .by = year)

write_rds(asfr_y1, "data_inter/ukr_asfr_wpp_2022_2025.rds")


