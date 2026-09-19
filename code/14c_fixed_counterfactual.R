# ==============================================================================
# STEP 14c - The loss with the counterfactual fixed at its point forecast
# ==============================================================================
#
# WHY
# ---
# The simulation's interval for the loss answers two questions at once: how
# uncertain the war's toll is, and how uncertain the counterfactual it is
# measured against is. From 2023 the forecast of the counterfactual carries
# about nine-tenths of the variance of the male loss (table A4). This step
# reruns the same draws with the counterfactual held at its point forecast,
# so the loss can be reported in two tiers: given the counterfactual, with
# every conflict and migration input drawn; and with the counterfactual drawn
# too (11).
#
# The draws are 11's own (simulation_draws(), same seed and size); only the
# forecast index is not applied, because the static inputs carry no Lee-Carter
# loadings (run_single_sim(), 00_setup.R). The projection is run in chunks and
# only the life expectancies are kept.
#
# INPUTS   as 11, and 14's draws of the loss with the counterfactual drawn
# OUTPUTS  data_inter/ukr_fixed_counterfactual_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()
stopifnot(!any(c("bx", "vk") %in% names(mi$static_inputs)))
mil <- read_rds("data_inter/ukr_military_inputs.rds")
lc_error <- read_rds("data_inter/ukr_lc_forecast_error.rds")
draws_list <- simulation_draws(mi$param_table, mil, lc_error, n = n_sim)$draws_df |> group_split(sim_id)

n_age <- 101L
e0_chunk <- function(ids) {
  d <- as.data.table(map_dfr(ids, \(i) run_single_sim(i, draws_list[[i]], mi$static_inputs, mi$pop22_ini)))
  setorder(d, sim_id, year, sex, age)
  key <- d[seq(1L, .N, by = n_age), .(sim_id, year, sex)]
  as_mat <- function(v) matrix(v, nrow = n_age)
  POP <- as_mat(d$pop)
  DEX <- as_mat(d$expected)
  DCF <- as_mat(d$civilian + d$combatant_confirmed + d$combatant_imputed)
  cbind(key, data.table(e0_bsn = lt_cols(DEX / POP, key$sex)$ex[1, ],
                        e0_war = lt_cols((DEX + DCF) / POP, key$sex)$ex[1, ]))
}
message("  running ", n_sim, " simulations with the counterfactual at its point forecast ...")
chunks <- split(seq_len(n_sim), ceiling(seq_len(n_sim) / 1000))
sims <- rbindlist(map(chunks, e0_chunk))
sims[, loss := e0_bsn - e0_war]

# the counterfactual is the same in every draw
stopifnot(sims[, sd(e0_bsn), by = .(year, sex)][, all(V1 < 1e-9)], nrow(sims) == n_sim * 8)

fixed <-
  as_tibble(sims) |>
  summarise(median = median(loss), lo = quantile(loss, 0.025, names = FALSE),
            hi = quantile(loss, 0.975, names = FALSE), sd = sd(loss), e0_bsn = e0_bsn[1],
            .by = c(year, sex)) |>
  mutate(counterfactual = "fixed at its point forecast")
drawn <-
  read_rds("data_inter/ukr_e0_loss_by_cause_draws_2022_2025.rds") |>
  as_tibble() |>
  summarise(median = median(loss_total), lo = quantile(loss_total, 0.025, names = FALSE),
            hi = quantile(loss_total, 0.975, names = FALSE), sd = sd(loss_total),
            e0_bsn = median(e0_bsn), .by = c(year, sex)) |>
  mutate(counterfactual = "drawn")
two_tiers <- bind_rows(fixed, drawn) |> arrange(year, sex, counterfactual)
write_rds(two_tiers, "data_inter/ukr_fixed_counterfactual_e0.rds")

cat("\n=== THE LOSS WITH THE COUNTERFACTUAL FIXED AND DRAWN ===\n")
print(as.data.frame(two_tiers |> mutate(across(c(median, lo, hi, sd, e0_bsn), \(x) round(x, 3)))))
message("Done. data_inter/ukr_fixed_counterfactual_e0.rds written.")
