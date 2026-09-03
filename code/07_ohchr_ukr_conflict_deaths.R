# ==============================================================================
# STEP 07 (OHCHR) - Age-sex PROFILE of civilian conflict deaths
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Produces the age-sex distribution of civilian conflict deaths. The yearly
# TOTALS come from UCDP (see 07_ucdp_ukr_conflict_deaths.R); this script only says how
# those totals are spread across age and sex.
#
#   1. OHCHR civilian casualties by five-year age group and sex, obtained
#      through UNICEF-IGME.
#   2. Deaths of unknown sex are redistributed proportionally over the known
#      male/female split within each age group.
#   3. Each year is rescaled to the OHCHR published annual total.
#   4. Five-year groups are ungrouped to single years of age with pclm.
#
# The saved column `cx` is the share of that year's civilian deaths falling
# in each age-sex cell, and sums to 1 within each year.
#
# INPUTS   data_input/ohchr_civilian_deaths.xlsx   (annual totals, sheet 1)
#          data_input/UKR_2022-2025_OHCHR.xlsx     (age-sex detail)
# OUTPUT   data_inter/ukr_ohchr_civilian_casualties.rds   <- used by 11, 13
# ==============================================================================

rm(list = ls())
source("code/00_setup.R")

# annual civilian death totals used to rescale the age-sex profile below.
# (sheets 2 and 3 of this workbook hold monthly and broad-age breakdowns that
# the pipeline does not use, so they are no longer read in)
yr <- read_xlsx("data_input/ohchr_civilian_deaths.xlsx", sheet = 1)

# OHCHR publishes a single documented count with no uncertainty bounds; saved
# here for the source comparison table assembled in 15.
write_rds(
  yr |> filter(year %in% 2022:2025),
  "data_inter/ukr_ohchr_annual_totals.rds"
)

# data obtained throuh UNICEF-IGME
dt <- read_xlsx("data_input/UKR_2022-2025_OHCHR.xlsx", skip = 1)

dt2 <-
  dt |>
  rename(age = 1) |>
  drop_na() |>
  gather(-age, key = sex_year, value = dx) |>
  separate(sex_year, c("sex", "year"), sep = " ") |>
  mutate(
    year = as.numeric(year),
    sex = str_sub(sex, 1, 1) |> str_to_lower(),
    age = ifelse(age == c("100-104"), "100", str_sub(age, 1, 2)),
    age = str_remove(age, "-"),
    age = as.numeric(age)
  ) |>
  spread(sex, dx) |>
  mutate(
    tot = f + m + u,
    tot_sum = m + f,
    f = ifelse(tot_sum == 0, 0, tot / tot_sum * f),
    m = ifelse(tot_sum == 0, 0, tot / tot_sum * m)
  ) |>
  select(age, year, f, m) |>
  gather(f, m, key = sex, value = dx) |>
  left_join(yr, by = join_by(year)) |>
  group_by(year) |>
  mutate(tot_sum = sum(dx)) |>
  ungroup() |>
  mutate(dx = dx * dts / tot_sum) |>
  select(-c(dts, tot_sum)) |>
  mutate(dx = ifelse(dx == 0, 1, dx))


chunk <-
  dt2 |>
  filter(year == 2025, sex == "f")

ung_age_ohchr <- function(chunk) {
  nl <- 1
  dxs <- pclm(x = chunk$age, y = chunk$dx, nlast = nl)$fitted
  fit <-
    tibble(age = 0:100, dx = dxs) |>
    mutate(sex = unique(chunk$sex))
  return(fit)
}

dt3 <-
  dt2 %>%
  group_by(year, sex) %>%
  do(ung_age_ohchr(chunk = .data)) %>%
  ungroup() |>
  select(year, sex, age, dx) |>
  arrange(year, sex, age) |>
  group_by(year) |>
  mutate(cx = dx / sum(dx)) |>
  ungroup()

write_rds(dt3, "data_inter/ukr_ohchr_civilian_casualties.rds")


dt3 |>
  mutate(cx = ifelse(sex == "m", -cx, cx)) |>
  ggplot() +
  geom_line(aes(age, cx, col = factor(year), group = interaction(sex, year))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw()

ggsave("figures/exploratory/ukr_unchr_civilians_2022_2025_cx.png")

dt3 |>
  mutate(dx = ifelse(sex == "m", -dx, dx)) |>
  ggplot() +
  geom_line(aes(age, dx, col = factor(year), group = interaction(sex, year))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw()

ggsave("figures/exploratory/ukr_unchr_civilians_2022_2025_dx.png")
