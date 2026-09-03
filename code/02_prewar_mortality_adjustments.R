# ==============================================================================
# STEP 02 - Old-age mortality: Kannisto extrapolation above age 90
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Observed death rates get noisy at the oldest ages, where exposures are tiny.
# A Kannisto model is fitted to ages 75-94 and extrapolated to 95-100.
#
# Observed and modelled rates are then BLENDED with a linear weight that runs
# from fully observed at age 80 to fully modelled at age 90, so the two join
# smoothly instead of jumping at a cut point.
#
# INPUT    data_input/DataDxEx.csv   (deaths and exposures, 1989-2021)
# OUTPUT   data_inter/ukr_mx_1989_2021_adj.rds     <- used by 03
# ==============================================================================

# Created 2026-02-02 to calculate the final tails for life tables

rm(list = ls())
source("code/00_setup.R")

# DATA MANIPULATION ===========================================================
# Read data
DataForKannisto <- read.csv("data_input/DataDxEx.csv", header = TRUE)

# Rename variables
dt2 <-
  DataForKannisto |>
  rename_with(
    .fn = ~ gsub("^age", "", .x),
    .cols = starts_with("age")
  ) |>
  pivot_longer(
    cols = matches("^\\d+$"), # all age columns
    names_to = "Age",
    values_to = "Value"
  ) |>
  mutate(Age = as.numeric(Age)) |>
  pivot_wider(
    names_from = Data,
    values_from = Value
  )

unique(dt2$Year)
unique(dt2$Region)

# NOTE: an exposures extract (Ex_txt -> ukr_pop_1989_2021_pavlo.rds) used to
# be built here. Step 04 reads DataDxEx.csv directly for the same exposures,
# so nothing ever consumed that file and it has been retired.


# We've agreed with the fitting interval
Adult <-
  dt2 |>
  filter(Age >= 75, Age <= 94)

# NESTED TIBBLE ===============================================================
# Store all age-specific mortality data for each group in a single row
Nested <-
  Adult |>
  group_by(Year, Region, Sex) |>
  nest()

# FIT KANNISTO MODEL WITHOUT LOOPS ============================================
Fits <- Nested |>
  mutate(
    fit = map(
      data,
      ~ MortalityLaw(
        x = .x$Age,
        Dx = .x$Dx,
        Ex = .x$Ex,
        law = "kannisto"
      )
    ),
    coef = map(fit, ~ .x$coefficients)
  )

# EXTRAPOLATION ===============================================================
Extrapolated <- Fits |>
  mutate(
    pred = map(
      fit,
      ~ tibble(
        Age = 95:100,
        mx = predict(.x, x = 95:100)
      )
    )
  ) |>
  unnest(pred)

# Drop modelling columns to avoid clutter for later use
Extrapolated_mx <- Extrapolated |>
  select(Year, Region, Sex, Age, mx)

# Check mx >= 1
any(Extrapolated_mx$mx >= 1, na.rm = TRUE)

# Run if you've got TRUE
Extrapolated_mx |> # Identify where
  ungroup() |>
  filter(mx >= 1)

which(Extrapolated_mx$mx >= 1) # this gives row numbers in the data frame

Extrapolated_mx |> # which (Year, Sex) combinations are problematic
  ungroup() |>
  filter(mx >= 1) |>
  count(Year, Sex)

# OBSERVED mx CALCULATION =====================================================
Observed_mx <-
  dt2 |>
  mutate(mx = Dx / Ex) |>
  select(Year, Region, Sex, Age, mx)

# SMOOTHING THE TRANSITION FROM OBSERVED TO EXTRAPOLATED mx ===================
# Fitted mx
Kannisto_fitted_mx <- Fits |>
  mutate(
    fitted = map(
      fit,
      ~ tibble(
        Age = as.numeric(names(.x$fitted.values)),
        mx_theoretical = as.numeric(.x$fitted.values)
      )
    )
  ) |>
  unnest(fitted) |>
  filter(Age >= 75, Age <= 94) |>
  select(Year, Region, Sex, Age, mx_theoretical)

# Extrapolated mx
Kannisto_extrapolated_mx <- Extrapolated |>
  select(
    Year,
    Region,
    Sex,
    Age,
    mx_theoretical = mx
  ) |>
  filter(Age >= 95, Age <= 100)

# Slap them together
Theoretical_mx <- bind_rows(
  Kannisto_fitted_mx,
  Kannisto_extrapolated_mx
) |>
  arrange(Year, Region, Sex, Age)

# Build mx_all
mx_all <- Observed_mx |>
  full_join(
    Theoretical_mx,
    by = c("Year", "Region", "Sex", "Age")
  ) |>
  arrange(Year, Region, Sex, Age)

# Define the weighting variable
mx_all <- mx_all |>
  mutate(
    w_obs = case_when(
      Age <= 80 ~ 1,
      Age >= 90 ~ 0,
      TRUE ~ (90 - Age) / (90 - 80)
    ),
    w_theoretical = 1 - w_obs
  )

# FINAL BLENDED mx ============================================================
mx_all <- mx_all |>
  mutate(
    mx_theoretical = coalesce(mx_theoretical, 0), # Replace NA in mx_theoretical with 0
    mx_final = (w_obs * mx) + (w_theoretical * mx_theoretical)
  )

# Check max value
mx_all |>
  filter(Age == 100) |>
  group_by(Sex) |>
  summarise(max_mx_100 = max(mx_final, na.rm = TRUE))

# # Save the data

# Take mx_final for LC model
mx <-
  mx_all |>
  select(-c(w_obs, w_theoretical, mx_theoretical, mx)) |>
  rename(mx = mx_final)

# saving adjusted mortality rates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_rds(mx, "data_inter/ukr_mx_1989_2021_adj.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# VISUALISATION ===============================================================
## Semi-automatic facet plot ====
plot_mortality_batch <- function(data, start_year, sex_label = "females") {
  end_year <- start_year + 8

  plot_data <- data |>
    filter(
      Sex == sex_label,
      Year >= start_year,
      Year <= end_year,
      Age >= 75,
      Age <= 100
    ) |>
    # Ensure we only look at one region if multiple exist
    filter(Region == Region[1]) |>
    select(
      Year,
      Age,
      Observed = mx,
      Theoretical = mx_theoretical,
      Blended = mx_final
    ) |>
    pivot_longer(
      cols = c(Observed, Theoretical, Blended),
      names_to = "Type",
      values_to = "Value"
    )

  ggplot(plot_data, aes(x = Age, y = Value, color = Type, linetype = Type)) +
    annotate(
      "rect",
      xmin = 80,
      xmax = 90,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.1,
      fill = "grey50"
    ) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~Year, ncol = 3, scales = "free_y") +
    scale_color_manual(
      values = c(
        "Blended" = "red",
        "Observed" = "black",
        "Theoretical" = "blue"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Blended" = "solid",
        "Observed" = "dashed",
        "Theoretical" = "dotted"
      )
    ) +
    labs(
      title = paste("Mortality:", sex_label, start_year, "-", end_year),
      subtitle = "Red line (Blended) should bridge the gap between Black (Observed) and Blue (Kannisto)",
      y = "mx"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(
        face = "bold"
      )
    )
}

# Usage:
plot_mortality_batch(mx_all, 1989)
plot_mortality_batch(mx_all, 1997)
plot_mortality_batch(mx_all, 2005)
plot_mortality_batch(mx_all, 2013)

# to save
plot_mortality_batch(mx_all, 2016)
ggsave("figures/exploratory/labtalk/kannisto_adj.png", w = 8, h = 4)

## Single year plot ====
plot_mortality_year <- function(data, year, sex_label = "females") {
  plot_data <- data |>
    filter(
      Sex == sex_label,
      Year == year,
      Age >= 75,
      Age <= 100
    ) |>
    # Ensure we only look at one region if multiple exist
    filter(Region == Region[1]) |>
    select(
      Year,
      Age,
      Observed = mx,
      Theoretical = mx_theoretical,
      Blended = mx_final
    ) |>
    pivot_longer(
      cols = c(Observed, Theoretical, Blended),
      names_to = "Type",
      values_to = "Value"
    )

  ggplot(plot_data, aes(x = Age, y = Value, color = Type, linetype = Type)) +
    annotate(
      "rect",
      xmin = 80,
      xmax = 90,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.1,
      fill = "grey50"
    ) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      values = c(
        "Blended" = "red",
        "Observed" = "black",
        "Theoretical" = "blue"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Blended" = "solid",
        "Observed" = "dashed",
        "Theoretical" = "dotted"
      )
    ) +
    labs(
      title = paste("Mortality:", sex_label, year),
      subtitle = "Red (Blended) bridges Black (Observed) and Blue (Kannisto)",
      y = "mx"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    )
}

# Usage
plot_mortality_year(mx_all, 1999)
