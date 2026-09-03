# ==============================================================================
# STEP 13 - Sensitivity: how much of the life expectancy loss is a
#           denominator effect of emigration?
# ==============================================================================
#
# WHY MIGRATION AFFECTS THE ESTIMATED LIFE EXPECTANCY LOSS
# --------------------------------------------------------
# Migration does not kill anyone, and in a closed accounting it would not
# change life expectancy at all. It matters here because of how the conflict
# deaths enter the calculation:
#
#   * conflict deaths are ABSOLUTE COUNTS taken from external sources (UCDP
#     for civilians, the ualosses register for combatants). They are not
#     modelled as a rate, so they do not shrink when the population shrinks.
#
#   * roughly 5.9 million people left Ukraine between 2022 and 2025. That
#     emigration is subtracted from the exposure, i.e. from the DENOMINATOR
#     of every death rate.
#
# The same numerator over a smaller denominator gives higher death rates, and
# therefore a larger estimated loss of life expectancy. This script quantifies
# how much of the estimated loss comes from that denominator effect, by
# re-running the identical projection with emigration switched off (ems = 0)
# and comparing.
#
# The implicit assumption is that the published conflict death counts refer to
# the population that REMAINED. If some of the deaths counted were of people
# who had already emigrated, the full-migration scenario overstates the rates,
# and the truth lies between the two lines below.
#
# NOTE: this is a deterministic sensitivity analysis. It runs the projection
# once at the MODE of each conflict-death distribution rather than over the
# 5,000 draws, because the question is about the direction and size of the
# migration effect, not about its uncertainty.
#
# INPUTS   the same as step 11
# OUTPUTS  figures/loss_decomp_migration_mortality*.png
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# lifetable() and run_single_sim() both come from 00_setup.R.

# 1. INPUTS ====================================================================
param_table <- read_rds("data_inter/ukr_param_table.rds")

# expected (non-conflict) mortality
exp_mort2 <- read_rds("data_inter/ukr_mxs_obs_plus_frcst_1989_2025.rds") %>%
  filter(year %in% 2022:2025, source == "frcst") %>%
  select(-source)

# baseline population, 1 January 2022, government-controlled mainland
pop22_ini <- read_rds("data_inter/ukr_pop_sssu.rds") %>%
  filter(reg == "cnt", year == 2022) %>%
  select(-reg)

# net emigration (positive = people leaving)
migs2 <- read_rds("data_inter/ukr_migrants_unchr_eurostat_sex_age_2022_2025.rds") %>%
  mutate(ems = -mix) |>
  select(-mix)

migs2 |> summarise(ems = sum(ems), .by = year)

# fertility: WPP2024 stops at 2023, so 2023 is carried forward
asfr <- read_rds("data_inter/ukr_asfr_wpp_2022_2025.rds")
asfr2 <- asfr %>%
  bind_rows(
    asfr %>% filter(year == 2023) %>% mutate(year = 2024),
    asfr %>% filter(year == 2023) %>% mutate(year = 2025)
  ) %>%
  arrange(year, age)

# 2. AGE-SEX PROFILES OF CONFLICT DEATHS =======================================
ohchr2 <- read_rds("data_inter/ukr_ohchr_civilian_casualties.rds") %>%
  select(year, sex, age, prop_cvs = cx)

ual2 <- read_rds("data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds") %>%
  filter(status == "dead", year %in% 2022:2025) |>
  select(-status) |>
  complete(year = 2022:2025, sex, age = 0:100, fill = list(dx = 0)) |>
  mutate(prop_cmb = dx / sum(dx), .by = year) |>
  select(-dx)

profile_props <- ual2 |> left_join(ohchr2, by = c("year", "sex", "age"))

# 3. DETERMINISTIC "DRAW" AT THE MODE ==========================================
draws_df <-
  param_table |>
  select(year, role, mode) |>
  spread(role, mode) %>%
  rename(draw_cmb = combatants, draw_cvs = civilians) |>
  left_join(
    param_table |> filter(role == "combatants") |> select(year, min_cmb = min),
    by = "year"
  ) |>
  mutate(sim_id = 1)

static_inputs <-
  expand_grid(
    year = 2022:2025,
    sex = c("f", "m"),
    age = 0:100
  ) %>%
  left_join(migs2, by = c("year", "sex", "age")) %>%
  left_join(exp_mort2, by = c("year", "sex", "age")) %>%
  left_join(asfr2, by = c("year", "sex", "age")) %>%
  left_join(profile_props, by = c("year", "sex", "age")) %>%
  replace_na(list(ems = 0, mx = 0, fx = 0, prop_cvs = 0, prop_cmb = 0))

stopifnot(all(static_inputs$mx > 0), !any(is.na(static_inputs)))

# 4. TWO SCENARIOS =============================================================
# full     - observed emigration
# no_mig   - identical, but nobody leaves
sim_full <- run_single_sim(1, draws_df, static_inputs, pop22_ini)

sim_nomig <- run_single_sim(
  1,
  draws_df,
  static_inputs |> mutate(ems = 0),
  pop22_ini
)

add_rates <- function(d) {
  d |>
    mutate(
      conflict = civilian + combatant_confirmed + combatant_imputed,
      all = expected + conflict,
      mx_all = all / pop,
      mx_bsn = expected / pop
    )
}

sim_full <- add_rates(sim_full)
sim_nomig <- add_rates(sim_nomig)

# how much larger is the exposure when nobody leaves?
bind_rows(
  sim_full |> summarise(pop = sum(pop), .by = year) |> mutate(sce = "full"),
  sim_nomig |> summarise(pop = sum(pop), .by = year) |> mutate(sce = "no migration")
) |>
  pivot_wider(names_from = sce, values_from = pop) |>
  mutate(pct_larger = round(100 * (`no migration` / full - 1), 1)) |>
  print()

# 5. LIFE EXPECTANCY UNDER EACH SCENARIO =======================================
e0_of <- function(d, col) {
  d |>
    select(year, sex, age, mx = all_of(col)) |>
    group_by(year, sex) |>
    do(lifetable(dt_in = .data)) |>
    ungroup() |>
    filter(age == 0) |>
    select(year, sex, ex)
}

lts <-
  e0_of(sim_full, "mx_bsn") |>
  rename(ex_bsn = ex) |>
  left_join(e0_of(sim_full, "mx_all") |> rename(ex_all = ex), by = c("year", "sex")) |>
  left_join(e0_of(sim_nomig, "mx_all") |> rename(ex_nomig = ex), by = c("year", "sex")) |>
  mutate(
    loss = ex_bsn - ex_all, # positive = years lost, with migration
    loss_nomig = ex_bsn - ex_nomig, # years lost from mortality alone
    diff = loss - loss_nomig, # the denominator effect
    cnt_mig = diff / loss # its share of the total loss
  )

print(lts)

# saved for the manuscript figure assembled in 16
write_rds(lts, "data_inter/ukr_migration_decomposition.rds")

# 6. PLOTS =====================================================================
plot_dat <-
  lts |>
  select(
    year,
    sex,
    Mortality = loss_nomig,
    Migration = diff
  ) |>
  pivot_longer(
    c(Mortality, Migration),
    names_to = "Component",
    values_to = "Loss"
  ) |>
  mutate(Component = factor(Component, levels = c("Mortality", "Migration")))

cols_comp <- c("Mortality" = "#555555", "Migration" = "#FF4D4D")

# absolute
plot_dat |>
  ggplot(aes(x = factor(year), y = Loss, fill = Component)) +
  geom_col(position = "stack", width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = if_else(Loss >= 0.05, sprintf("%.2f", Loss), "")),
    position = position_stack(vjust = 0.5),
    colour = "black", size = 3
  ) +
  facet_wrap(
    ~sex,
    labeller = labeller(sex = c("f" = "Females", "m" = "Males"))
  ) +
  scale_fill_manual(values = cols_comp) +
  labs(
    x = "Year",
    y = "Life expectancy loss (years)",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("figures/exploratory/loss_decomp_migration_mortality_abs.png", w = 7, h = 4)

# relative
plot_dat |>
  mutate(pct = Loss / sum(Loss), .by = c(year, sex)) |>
  ggplot(aes(x = factor(year), y = Loss, fill = Component)) +
  geom_col(position = "fill", width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = scales::percent(pct, accuracy = 0.1)),
    position = position_fill(vjust = 0.5),
    colour = "black", fontface = "bold", size = 3
  ) +
  facet_wrap(
    ~sex,
    labeller = labeller(sex = c("f" = "Females", "m" = "Males"))
  ) +
  scale_fill_manual(values = cols_comp) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Year",
    y = "Percentage contribution to loss",
    fill = "Effect"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("figures/exploratory/loss_decomp_migration_mortality.png", w = 7, h = 4)

message("Done.")
