# ==============================================================================
# STEP 13f - Sensitivity to the timing of deaths and departures within 2022
# ==============================================================================
#
# WHY
# ---
# The projection spreads each year's net outflow and conflict deaths evenly
# over the year: half of them are taken off before the year's exposure is
# counted (run_single_sim(), the mid-year convention). 2022 was not like that.
# The war began on 24 February, most departures fell in March and April, and
# civilian deaths peaked in March. Leaving earlier in the year means fewer
# person-years lived in Ukraine, so the mid-year convention overstates 2022's
# exposure and understates its rates.
#
# METHOD
# ------
# A person who leaves or dies at a fraction u of the year lives u of it, so
# the share of the year's events to take off before exposure is 1 - mean(u).
# The mean timing comes from the months of the events:
#   civilian deaths   OHCHR's monthly count of civilians killed in 2022
#   military deaths   the register's dead and missing by month of the event
#                     (v19, from 08b)
#   net outflow       no monthly series covers it (the border data are annual
#                     and Eurostat's monthly decisions start in August), so a
#                     scenario: the year's net outflow at the turn of March
#                     and April, where the reports place the bulk of it
# Only 2022 changes; later years keep the mid-year convention. Every input is
# at its mode.
#
# INPUTS   data_input/ohchr_civilian_deaths.xlsx (sheet "monthly")
#          data_inter/ukr_registration_completion.rds (08b)
#          the static inputs at the mode
# OUTPUTS  data_inter/ukr_timing_2022.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()

# the fraction of 2022 elapsed at the middle of each month of the war; February
# counts from 24 February, so its middle is 26 February
month_mid <- tibble(month = 2:12) |>
  mutate(mid = if_else(month == 2, as.Date("2022-02-26"),
                       as.Date(sprintf("2022-%02d-15", month))),
         u = as.numeric(mid - as.Date("2022-01-01")) / 365)

civilian <-
  read_xlsx("data_input/ohchr_civilian_deaths.xlsx", sheet = "monthly") |>
  filter(Year == 2022) |>
  transmute(month = match(Month, month.abb), n = Killed)
military <-
  read_rds("data_inter/ukr_registration_completion.rds")$month |>
  filter(year == 2022) |>
  summarise(n = sum(n_v19), .by = m) |>
  transmute(month = month(m), n)

mean_u <- function(d) d |> left_join(month_mid, by = "month") |> with(sum(n * u) / sum(n))
timing <- tibble(
  component = c("civilian deaths", "military deaths", "net outflow"),
  mean_u = c(mean_u(civilian), mean_u(military),
             as.numeric(as.Date("2022-04-01") - as.Date("2022-01-01")) / 365)
) |>
  mutate(absent = 1 - mean_u)
cat("\n=== MEAN TIMING WITHIN 2022, AS A FRACTION OF THE YEAR ===\n")
print(as.data.frame(timing |> mutate(across(where(is.numeric), \(x) round(x, 3)))))

with_timing <- function(cvs, cmb, mig) {
  mi$draws |>
    mutate(absent_cvs = if_else(year == 2022, cvs, 0.5),
           absent_cmb = if_else(year == 2022, cmb, 0.5),
           absent_mig = if_else(year == 2022, mig, 0.5))
}
a <- set_names(timing$absent, timing$component)
scenarios <- list(
  "mid-year (used)" = with_timing(0.5, 0.5, 0.5),
  "deaths at their month" = with_timing(a[["civilian deaths"]], a[["military deaths"]], 0.5),
  "deaths at their month, net outflow at 1 April" =
    with_timing(a[["civilian deaths"]], a[["military deaths"]], a[["net outflow"]])
)
loss <- imap_dfr(scenarios, \(d, s) loss_at_mode(d, mi$static_inputs, mi$pop22_ini) |> mutate(scenario = s))

# the convention used must reproduce 13's loss at the mode
stopifnot(isTRUE(all.equal(
  loss |> filter(scenario == "mid-year (used)") |> arrange(year, sex) |> pull(loss),
  read_rds("data_inter/ukr_migration_decomposition.rds") |> arrange(year, sex) |> pull(loss),
  tolerance = 1e-8
)))

write_rds(list(timing = timing, loss = loss), "data_inter/ukr_timing_2022.rds")

cat("\n=== LOSS BY THE TIMING OF 2022 EVENTS ===\n")
print(as.data.frame(
  loss |> mutate(col = paste0(sex, "_", year), loss = round(loss, 3)) |>
    select(scenario, col, loss) |> pivot_wider(names_from = col, values_from = loss)
))
message("Done. data_inter/ukr_timing_2022.rds written.")
