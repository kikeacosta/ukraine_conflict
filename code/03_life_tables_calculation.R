# ==============================================================================
# STEP 03 - Period life tables, 1989-2021
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Builds full period life tables from the adjusted rates of step 02.
#   - infancy uses the a0 separation factors supplied in a0.csv
#   - ages 1-99 use the standard ax = 0.5 assumption
#   - the open interval at 100+ uses ax = 1/mx
# lx monotonicity is checked before the tables are saved.
#
# INPUTS   data_inter/ukr_mx_1989_2021_adj.rds   (from 02)
#          data_input/a0.csv
# OUTPUT   data_inter/ukr_life_tables_1989_2021.rds   <- used by 04
# ==============================================================================

# Created 2026-02-03 to calculate life tables
rm(list = ls())
source("code/00_setup.R")

# DATA MANIPULATION ===========================================================
# Read a0 values
a0 <- read.csv("data_input/a0.csv", header = TRUE)

# Load mx data
mx <- read_rds("data_inter/ukr_mx_1989_2021_adj.rds")

# Join a0 with mx_LT
mxax <-
  mx |>
  left_join(
    a0 |> select(Year, Region, Sex, Age, ax_input = ax),
    by = c("Year", "Region", "Sex", "Age")
  ) |>
  rename_with(tolower) |>
  mutate(
    ax = case_when(
      age == 0 ~ ax_input, # Use value from a0 data frame
      age == 100 ~ 1 / mx, # Calculate 1/mx for age 100
      TRUE ~ 0.5 # Default assumption for other ages (1-99)
    )
  ) |>
  # 3. Clean up the temporary column
  select(-ax_input)

# LIFE TABLE ==================================================================
# qx calculation
lt <-
  mxax |>
  mutate(qx = (mx) / (1 + (1 - ax) * mx)) |> # qx calculation
  arrange(year, sex, age) |> # lx calculation
  group_by(year, region, sex) |>
  mutate(lx = cumprod(c(1, 1 - qx[-length(qx)]))) |>
  mutate(
    dx = lx - dplyr::lead(lx, default = 0),
    Lx = lx - (1 - ax) * dx,
    Tx = rev(cumsum(rev(Lx))),
    ex = Tx / lx
  ) |>
  ungroup()

lt |>
  group_by(year, region, sex) |>
  summarise(n_viol = sum(diff(lx) > 0, na.rm = TRUE), .groups = "drop") |>
  filter(n_viol > 0)

lt |>
  group_by(year, region, sex) |>
  summarise(n_viol = sum(diff(lx) > 0, na.rm = TRUE), .groups = "drop") |>
  summarise(max_found = max(n_viol)) # the highest number of violations found in any group

# Save for later use
# save(LT_UA, file = "LT_UA.RData")

# saving life tables
# ~~~~~~~~~~~~~~~~~~
write_rds(lt, "data_inter/ukr_life_tables_1989_2021.rds")


# VISUALISATION ===============================================================
# e0
lt |>
  filter(age == 0) |>
  ggplot(aes(x = year, y = ex, color = sex)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(60, 78)) +
  labs(
    title = "Life expectancy at birth",
    y = "Years",
    x = "Year",
    color = "Sex"
  ) +
  theme_minimal()

# qx
qx_plot <- lt |>
  filter(
    sex == "males", # swap females/males
    age >= 0,
    age <= 99
  )

ggplot(qx_plot, aes(x = age, y = qx, color = year, group = year)) +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_viridis_c(option = "turbo", name = "Year", direction = -1) +
  labs(
    x = "Age",
    y = "qx (log scale)",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "gray40", hjust = 0),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    legend.position = c(0.15, 0.8),
    legend.background = element_rect(
      fill = "white",
      color = "gray80",
      linewidth = 0.3
    ),
    legend.key.width = unit(1.5, "cm"),
    plot.caption = element_text(color = "gray50", size = 9, hjust = 1)
  ) +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))
