# ==============================================================================
# STEP 14d - The loss with the population base of Donetsk and Luhansk drawn
# ==============================================================================
#
# WHY
# ---
# The base holds the State Statistics Service's estimates for Donetsk and
# Luhansk, which assume the two oblasts were fully registered in 2022. They are
# 6.14 of the 41.0 million the projection starts from, and no source can check
# them. Held at one value at a time, the base moves the loss further than any
# input the simulation draws (table A11): the 2025 male loss runs from 6.63
# years as used to 7.74 with the two oblasts taken out, against a drawn
# interval of 6.10-7.15.
#
# So the loss is reported in two tiers, as it is for the counterfactual (14c):
#   tier 1  the base at the published estimates, every other input drawn - the
#           estimate of the paper, resting on a source that can be named;
#   tier 2  the base of the two oblasts drawn as well, PERT(0, 0.5, 1.0 million
#           fewer, shape 4), one draw per simulation, everything else as in 11.
#
# THE RANGE IS ONE-SIDED, AND IT IS A JUDGEMENT
# ---------------------------------------------
# Every indication runs the same way - the government's July 2024 figure of
# about 1.8 million in the occupied areas outside Crimea, the 2019 electronic
# census, the registered deaths against the projection's rates (Discussion D4) -
# the base is too large, never too small. A symmetric draw would therefore be
# wrong, and no source fixes the floor, the mode or the ceiling: the three come
# from the team's reading of that evidence. Because the range is one-sided, the
# MEDIAN moves too, and tier 2's central value is not tier 1's. That is why this
# is reported beside the estimate and not inside it.
#
# INPUTS   as 11 (draws, static inputs, the forecast loadings) and
#          data_inter/ukr_pop_sssu.rds; 14's draws of the loss for tier 1
# OUTPUTS  data_inter/ukr_drawn_base_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

base_cut <- c(min = 0, mode = 0.5e6, max = 1.0e6)  # people fewer in the two oblasts
base_seed <- 4114                                  # its own stream: 11's draws are untouched

mi <- mode_projection_inputs()
mil <- read_rds("data_inter/ukr_military_inputs.rds")
lc_error <- read_rds("data_inter/ukr_lc_forecast_error.rds")
draws_list <- simulation_draws(mi$param_table, mil, lc_error, n = n_sim)$draws_df |> group_split(sim_id)

# The static inputs of 11: those of the deterministic steps, with the loading
# and the forecast variance that carry the counterfactual's own uncertainty.
static_inputs <-
  mi$static_inputs |>
  left_join(lc_error$bx, by = c("sex", "age")) |>
  left_join(lc_error$variance |> select(sex, year, vk), by = c("sex", "year"))
stopifnot(all(c("bx", "vk") %in% names(static_inputs)), !any(is.na(static_inputs)),
          nrow(static_inputs) == 4 * 2 * 101)

# 1. THE DRAWN BASE ============================================================
pop_dl <-
  read_rds("data_inter/ukr_pop_sssu.rds") |>
  filter(reg %in% c("dnk", "luk"), year == 2022) |>
  summarise(pop_dl = sum(pop), .by = c(sex, age))
dl_total <- sum(pop_dl$pop_dl)
ini <-
  mi$pop22_ini |>
  left_join(pop_dl, by = c("sex", "age")) |>
  mutate(pop_dl = coalesce(pop_dl, 0))
stopifnot(all(ini$pop_dl <= ini$pop))

set.seed(base_seed)
cut_draw <- qpert(runif(n_sim), min = base_cut[["min"]], mode = base_cut[["mode"]],
                  max = base_cut[["max"]], shape = 4)
keep <- 1 - cut_draw / dl_total
cat(sprintf("\nDonetsk and Luhansk: %.2f million of %.2f; the draw takes out %.0f (%.0f-%.0f) thousand\n",
            dl_total / 1e6, sum(mi$pop22_ini$pop) / 1e6, median(cut_draw) / 1e3,
            quantile(cut_draw, 0.025) / 1e3, quantile(cut_draw, 0.975) / 1e3))

base_for <- function(i) ini |> mutate(pop = pop - (1 - keep[i]) * pop_dl) |> select(-pop_dl)

# 2. THE PROJECTION, AS 14c RUNS IT ============================================
n_age <- 101L
e0_chunk <- function(ids, bases) {
  d <- as.data.table(map_dfr(ids, \(i) run_single_sim(i, draws_list[[i]], static_inputs, bases(i))))
  setorder(d, sim_id, year, sex, age)
  key <- d[seq(1L, .N, by = n_age), .(sim_id, year, sex)]
  as_mat <- function(v) matrix(v, nrow = n_age)
  POP <- as_mat(d$pop)
  DEX <- as_mat(d$expected)
  DCF <- as_mat(d$civilian + d$combatant_confirmed + d$combatant_imputed)
  cbind(key, data.table(e0_bsn = lt_cols(DEX / POP, key$sex)$ex[1, ],
                        e0_war = lt_cols((DEX + DCF) / POP, key$sex)$ex[1, ]))
}

# The check: with nothing taken out, this must reproduce 11's own loss. The
# draws only match 11's when the run has its size: simulation_draws() consumes
# one stream, so n draws are not the first n of a longer run.
ref_all <- read_rds("data_inter/ukr_e0_loss_by_cause_draws_2022_2025.rds") |> as_tibble()
if (n_distinct(ref_all$sim_id) == n_sim) {
  check_ids <- seq_len(min(100L, n_sim))
  chk <- e0_chunk(check_ids, \(i) mi$pop22_ini)
  chk[, loss := e0_bsn - e0_war]
  cmp <-
    as_tibble(chk) |>
    inner_join(ref_all |> filter(sim_id %in% check_ids) |> select(sim_id, year, sex, loss_total),
               by = c("sim_id", "year", "sex"))
  stopifnot(nrow(cmp) == length(check_ids) * 8, max(abs(cmp$loss - cmp$loss_total)) < 1e-9)
  cat("  the base unchanged reproduces 11's loss in", length(check_ids), "simulations\n")
} else {
  message("  11's draws are of a different size (", n_distinct(ref_all$sim_id),
          " against ", n_sim, "), so the reproduction check is skipped")
}

message("  running ", n_sim, " simulations with the population base drawn ...")
chunks <- split(seq_len(n_sim), ceiling(seq_len(n_sim) / 1000))
sims <- rbindlist(map(chunks, \(ids) e0_chunk(ids, base_for)))
sims[, loss := e0_bsn - e0_war]
stopifnot(nrow(sims) == n_sim * 8)

# 3. THE TWO TIERS =============================================================
drawn_base <-
  as_tibble(sims) |>
  summarise(median = median(loss), lo = quantile(loss, 0.025, names = FALSE),
            hi = quantile(loss, 0.975, names = FALSE), sd = sd(loss),
            e0_bsn = median(e0_bsn), .by = c(year, sex)) |>
  mutate(base = "drawn")
fixed_base <-
  ref_all |>
  summarise(median = median(loss_total), lo = quantile(loss_total, 0.025, names = FALSE),
            hi = quantile(loss_total, 0.975, names = FALSE), sd = sd(loss_total),
            e0_bsn = median(e0_bsn), .by = c(year, sex)) |>
  mutate(base = "at the published estimates")
two_tiers <- bind_rows(fixed_base, drawn_base) |> arrange(year, sex, base)

write_rds(list(tiers = two_tiers, cut = base_cut, seed = base_seed,
               dl_total = dl_total, cut_draw = cut_draw),
          "data_inter/ukr_drawn_base_e0.rds")

cat("\n=== THE LOSS WITH THE POPULATION BASE FIXED AND DRAWN ===\n")
print(as.data.frame(two_tiers |> mutate(across(c(median, lo, hi, sd, e0_bsn), \(x) round(x, 3)))))
message("Done. data_inter/ukr_drawn_base_e0.rds written.")
