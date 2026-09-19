# ==============================================================================
# STEP 14b - Years of life lost and adult mortality (45q15)
# ==============================================================================
#
# Two summary measures besides the loss of life expectancy, computed inside
# every draw as 14 does, so they carry the same uncertainty intervals:
#
#   45q15  the probability that someone aged 15 dies before 60, under the
#          counterfactual and with conflict deaths added, and the difference
#   YLL    years of life lost: each conflict death weighted by the life
#          expectancy it cut short, the counterfactual's at the age of death.
#          A death at completed age x falls on average at exact age x + 1/2,
#          so the weight is the mean of the counterfactual e(x) and e(x + 1).
#          By cause, with the late registrations split from the confirmed
#          deaths as in 14.
#
# INPUTS   data_inter/ukr_sim_draws_2022_2025_n<n_sim>.rds   (from 11)
#          data_inter/ukr_ualosses_imputation_table.rds      (from 09)
#          data_inter/ukr_sim_param_draws.rds               (from 11)
# OUTPUTS  data_inter/ukr_yll_45q15_summary.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

draws_file <- sprintf("data_inter/ukr_sim_draws_2022_2025_n%d.rds", n_sim)
if (!file.exists(draws_file)) stop("Run step 11 first (at n_sim = ", n_sim, ").", call. = FALSE)

n_age <- 101
d <- as.data.table(readRDS(draws_file))
setorder(d, sim_id, year, sex, age)
key <- d[seq(1L, .N, by = n_age), .(sim_id, year, sex)]
as_mat <- function(v) matrix(v, nrow = n_age)

POP <- as_mat(d$pop)
DEX <- as_mat(d$expected)
DCV <- as_mat(d$civilian)
DCC <- as_mat(d$combatant_confirmed)
DCI <- as_mat(d$combatant_imputed)
rm(d)
gc()

lb <- lt_cols(DEX / POP, key$sex)
lw <- lt_cols((DEX + DCV + DCC + DCI) / POP, key$sex)

# 45q15: rows are ages 0-100, so age 15 is row 16 and age 60 row 61
key[, `:=`(q4515_bsn = 1 - lb$lx[61, ] / lb$lx[16, ],
           q4515_war = 1 - lw$lx[61, ] / lw$lx[16, ])]
key[, q4515_diff := q4515_war - q4515_bsn]

# the counterfactual life expectancy at the mean age of death in each interval
e_mid <- (lb$ex + rbind(lb$ex[-1, , drop = FALSE], lb$ex[n_age, ])) / 2
# each draw's share of late registrations in its confirmed deaths, as in 14
late_share <- late_share_draws()
ls_vec <- late_share$late_share[match(paste(key$sim_id, key$year), paste(late_share$sim_id, late_share$year))]
stopifnot(!anyNA(ls_vec))
yll_civ <- colSums(DCV * e_mid)
yll_conf <- colSums(DCC * e_mid)
key[, `:=`(yll_civilian = yll_civ,
           yll_registered = yll_conf * (1 - ls_vec),
           yll_late = yll_conf * ls_vec,
           yll_imputed = colSums(DCI * e_mid),
           deaths = colSums(DCV + DCC + DCI))]
key[, yll_total := yll_civilian + yll_registered + yll_late + yll_imputed]

q3 <- function(x) list(median = median(x), lo = quantile(x, 0.025), hi = quantile(x, 0.975))
long <- melt(key, id.vars = c("sim_id", "year", "sex"), variable.name = "measure")

by_year_sex <- long[, q3(value), by = .(year, sex, measure)]
# both sexes, and all four years, summed within each draw
yll_cols <- c("yll_civilian", "yll_registered", "yll_late", "yll_imputed", "yll_total", "deaths")
both_sexes <- long[measure %in% yll_cols, .(value = sum(value)), by = .(sim_id, year, measure)][
  , q3(value), by = .(year, measure)]
all_years <- long[measure %in% yll_cols, .(value = sum(value)), by = .(sim_id, sex, measure)][
  , q3(value), by = .(sex, measure)]
all_years_both <- long[measure %in% yll_cols, .(value = sum(value)), by = .(sim_id, measure)][
  , q3(value), by = .(measure)]

out <- list(
  by_year_sex = as_tibble(by_year_sex),
  both_sexes = as_tibble(both_sexes),
  all_years = as_tibble(all_years),
  all_years_both = as_tibble(all_years_both),
  n_sim = n_sim
)
write_rds(out, "data_inter/ukr_yll_45q15_summary.rds")

cat("\n=== 45q15 (per 1,000): COUNTERFACTUAL, WITH CONFLICT DEATHS, DIFFERENCE ===\n")
print(as.data.frame(
  out$by_year_sex |> filter(str_starts(measure, "q4515")) |>
    mutate(across(c(median, lo, hi), \(x) round(1000 * x, 1))) |>
    arrange(sex, measure, year)
))
cat("\n=== YEARS OF LIFE LOST, BOTH SEXES, ALL FOUR YEARS ===\n")
print(as.data.frame(out$all_years_both |> mutate(across(c(median, lo, hi), round))))
message("Done. data_inter/ukr_yll_45q15_summary.rds written.")
