rm(list = ls())
gc()
source("code/00_setup.R")

dt <- read_rds("data_inter/ukr_expected_and_conflict_deaths_2022_2024.rds")

unique(dt$cause)
unique(dt$sce)

dt2 <-
  dt %>%
  select(sce, year, sex, age, cause, dx, mx, exposure = pop) %>%
  mutate(
    cause = recode(
      cause,
      "all" = "observed",
      "expected" = "expected",
      "conflict" = "conflict"
    ),
    sex = recode(sex, "f" = "female", "m" = "male", "t" = "total")
  )

tst_all <-
  dt2 %>%
  filter(
    sce == "ualosses_ucdp",
    sex == "male",
    year == 2022,
    cause == "observed"
  )

tst_exp <-
  dt2 %>%
  filter(
    sce == "ualosses_ucdp",
    sex == "male",
    year == 2022,
    cause == "expected"
  )

lifetable(tst_all) %>% filter(age == 0)
lifetable(tst_exp) %>% filter(age == 0)

lt2 <-
  dt2 %>%
  filter(cause != "conflict") %>%
  group_by(sce, cause, year, sex) %>%
  do(lifetable(dt_in = .data)) %>%
  ungroup() %>%
  filter(age == 0)

lt3 <-
  lt2 %>%
  select(cause, year, sex, sce, ex) %>%
  spread(cause, ex) %>%
  mutate(loss = observed - expected)

lt3 %>%
  ggplot() +
  geom_point(aes(factor(year), loss, col = sce)) +
  facet_grid(~sex) +
  theme_bw()

lt3 %>%
  select(-loss) |>
  gather(expected, observed, key = "cause", value = "e0") %>%
  ggplot() +
  geom_point(aes(factor(year), e0, col = cause)) +
  facet_grid(sex ~ sce) +
  coord_cartesian(ylim = c(60, 87)) +
  theme_bw()

unique(dt$cause)

# pyramids population
dt2 %>%
  filter(sex != "total") %>%
  mutate(exposure = ifelse(sex == "female", -exposure, exposure)) %>%
  ggplot() +
  geom_line(aes(
    age,
    exposure,
    col = factor(year),
    group = interaction(sex, year)
  )) +
  geom_hline(yintercept = 0) +
  geom_text(aes(100, 0), label = "females", hjust = 1.5) +
  geom_text(aes(100, 0), label = "males", hjust = -.5) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_color_manual(values = c("black", "grey30", "grey50")) +
  facet_grid(~sce) +
  coord_flip() +
  theme_bw()
ggsave("figures/pops.png", w = 10, h = 5)

# pyramids conflict deaths
dt %>%
  filter(cause == "conflict", sex != "total") %>%
  mutate(dx = ifelse(sex == "female", -dx, dx)) %>%
  ggplot() +
  geom_line(aes(age, dx, col = factor(year), group = interaction(sex, year))) +
  geom_hline(yintercept = 0) +
  geom_text(aes(100, 0), label = "females", hjust = 1.5) +
  geom_text(aes(100, 0), label = "males", hjust = -.5) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_color_manual(values = c("black", "grey30", "grey50")) +
  facet_grid(~sce) +
  labs(title = "Conflict Deaths") +
  coord_flip() +
  theme_bw()
ggsave("figures/conf_deaths.png", w = 10, h = 5)

dt2 %>%
  filter(cause != "conflict") %>%
  ggplot() +
  geom_line(aes(age, mx, col = cause, group = interaction(sex, cause))) +
  geom_hline(yintercept = 0) +
  geom_text(
    data = lt3,
    aes(x = 80, y = .005, label = paste0("e0: ", round(expected, 2))),
    col = "black",
    size = 3
  ) +
  geom_text(
    data = lt3,
    aes(x = 80, y = .001, label = paste0("e0: ", round(observed, 2))),
    col = "blue",
    size = 3
  ) +
  geom_text(
    data = lt3,
    aes(x = 80, y = .0003, label = paste0("loss: ", round(loss, 2))),
    col = "red",
    size = 3
  ) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  facet_nested(sex + sce ~ year) +
  scale_color_manual(values = c("black", "blue")) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("figures/ex_scenarios.png", w = 15, h = 10)

dt2 %>%
  filter(cause == "conflict") %>%
  ggplot() +
  geom_line(aes(age, mx, col = sce, group = interaction(sex, sce))) +
  geom_hline(yintercept = 0) +
  scale_y_log10() +
  facet_grid(sex ~ year) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_color_manual(values = c("black", "blue", "red")) +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave("figures/mx_conflict.png", w = 10, h = 5)


unique(dt_out2$sex)

dt <-
  dt2 %>%
  summarise(
    dx = sum(dx),
    exposure = sum(exposure),
    .by = c(cause, year, age, sce)
  ) %>%
  mutate(sex = "total") %>%
  bind_rows(dt2) %>%
  mutate(
    age = case_when(
      age == 0 ~ 0,
      age %in% 1:4 ~ 1,
      age %in% 5:79 ~ age - age %% 5,
      age >= 80 ~ 80
    )
  ) %>%
  summarise(
    dx = sum(dx),
    pop = sum(exposure),
    .by = c(cause, year, sex, age, sce)
  ) %>%
  mutate(mx = dx / pop)

bsns <-
  dt %>%
  filter(cause == "expected") %>%
  select(everything(), -cause, -pop, -dx, bsn = mx)

psc_age_sex <-
  dt %>%
  filter(cause == "observed") %>%
  left_join(bsns) %>%
  mutate(psc = mx / bsn) %>%
  select(-bsn)

# plotting sex-age-specific mortality relative risks

cols <- c("#ae2012", "#03071e", "#00a6fb")

psc_age_sex %>%
  # filter(sex == "conflict") %>%
  mutate(age_ad = ifelse(age == 80, 85, (age + lead(age)) / 2)) %>%
  ggplot() +
  # geom_segment(aes(x = age_ad, y = psc_l, xend = age_ad, yend = psc_u, col = sex),
  #              alpha = 0.5)+
  geom_point(aes(age_ad, psc, col = sex, group = interaction(sex))) +
  geom_hline(yintercept = 1, lty = "dashed") +
  scale_y_log10(breaks = c(1, 2, 3, 4)) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  scale_color_manual(values = cols) +
  facet_grid(sce ~ year) +
  labs(
    x = "Age",
    y = "Mortality Relative Risks (times)",
    col = "Sex",
    fill = "Sex",
    lty = "Mortality"
  ) +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("figures/ukr_igme2025_rr_all_ages.png", w = 7, h = 4)
