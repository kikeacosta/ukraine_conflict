rm(list = ls(all = TRUE))

## loading LC functions
source("code/00_setup.R")

# loading data
mx_all <- read_rds("data_inter/ukr_life_tables_1989_2021.rds")

unique(mx_all$region)
unique(mx_all$year)

dt <-
  mx_all %>%
  select(year, sex, age, mx) %>%
  mutate(sex = str_sub(sex, 1, 1))

unique(dt$year)

dt %>%
  filter(year %in% c(2010, 2015, 2021)) %>%
  ggplot() +
  facet_wrap(~sex, scales = "free_y") +
  geom_line(aes(x = age, mx, col = factor(year), group = interaction(year))) +
  scale_y_log10() +
  theme_bw()

# ggsave("figures/mx_grigoriev_rates.png",width = 8,height = 4)

# exposures
pop_cnt <- read.csv("data_input/DataDxEx.csv", header = TRUE)

pop_cnt2 <-
  pop_cnt %>%
  rename_with(tolower) %>%
  filter(data == "Ex") %>%
  gather(starts_with("age"), key = "age", value = "pop") |>
  mutate(age = as.numeric(gsub("age", "", age)), sex = str_sub(sex, 1, 1)) |>
  select(year, sex, age, pop)

# write_rds(pop_cnt2, "data_inter/ukr_pop_sssu.rds")

data <-
  dt %>%
  left_join(pop_cnt2) %>%
  drop_na(pop) %>%
  mutate(dts = mx * pop)

# LC forecasting
tst_19 <-
  data %>%
  group_by(sex) %>%
  arrange(year, sex, age) %>%
  do(frcst(data = .data, t = 2010:2019)) %>%
  ungroup() %>%
  mutate(mx = mx / 1e5, source = "lc_from_2019")

tst_21 <-
  data %>%
  group_by(sex) %>%
  arrange(year, sex, age) %>%
  do(frcst(data = .data, t = 2012:2021)) %>%
  ungroup() %>%
  mutate(mx = mx / 1e5, source = "lc_from_2021")

all_mxs <-
  data %>%
  select(year, sex, age, mx) %>%
  mutate(source = "obs") %>%
  bind_rows(tst_19, tst_21)

## plotting e0 outcomes
all_mxs %>%
  ## remove 2023 from wpp
  filter(!(source %in% c("obs") & year %in% c(2023, 2024))) %>%
  group_by(year, sex, source) %>%
  summarise(e0 = e0.mx(x = age, mx = mx)) %>%
  ggplot(aes(x = year, y = e0, color = source)) +
  facet_wrap(~sex, scales = "free_y") +
  geom_vline(xintercept = 2019, linetype = "dotted") +
  geom_vline(xintercept = 2022, linetype = "dotted") +
  geom_point() +
  geom_line() +
  coord_cartesian(ylim = c(60, 79)) +
  theme_bw()

# ggsave("figures/frcst_rates_e0_pavlo.png", width = 8, height = 6)

all_mxs %>%
  filter(
    !(source %in% c("obs") & year %in% c(2023, 2024)),
    year %in% 2020:2021,
    !(source %in% c("lc_from_2021"))
  ) %>%
  ggplot(aes(x = age, mx, lty = source)) +
  facet_wrap(year ~ sex, scales = "free_y") +
  scale_y_log10() +
  geom_line() +
  theme_bw()

unique(all_mxs$source)

all_mxs2 <-
  all_mxs %>%
  filter(source %in% c("lc_from_2019", "obs"))

## saving data
write_rds(
  all_mxs2,
  file = "data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds"
)
