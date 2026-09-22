# Run full pipeline steps 01-15
# Configure CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("========================================\n")
cat("Ukraine Conflict Mortality Pipeline\n")
cat("Running all steps in sequence\n")
cat("========================================\n\n")

start_time <- Sys.time()

# Execute each step
steps_to_run <- list(
  list(name = "00", file = "code/00_setup.R"),
  list(name = "01", file = "code/01_pop_estimates_sssu.R"),
  list(name = "02", file = "code/02_prewar_mortality_adjustments.R"),
  list(name = "03", file = "code/03_life_tables_calculation.R"),
  list(name = "04", file = "code/04_mort_forecast_lee_carter_miller.R"),
  list(name = "05", file = "code/05_asfr_wpp.R"),
  list(name = "06", file = "code/06_net_migration.R"),
  list(name = "07a", file = "code/07_ucdp_ukr_conflict_deaths.R"),
  list(name = "07b", file = "code/07_ohchr_ukr_conflict_deaths.R"),
  # 07c feeds the UCDP-vs-ACLED reconciliation in table A1; without it 15
  # stops on a missing ukr_rus_acled_source_comparison.rds
  list(name = "07c", file = "code/07_acled_ukr_conflict_deaths.R"),
  list(name = "08", file = "code/08_ualosses.R"),
  list(name = "08b", file = "code/08b_registration_lag.R"),
  list(name = "09", file = "code/09_ualosses_missing_analysis.R"),
  # 09j's table A36 is one of 15's deliverables, and 09k's tests are cited in the
  # documentation; both read 09's cached linkage, so neither needs the register files
  list(name = "09j", file = "code/09j_held_out_window.R"),
  list(name = "09k", file = "code/09k_cohort_overlap_test.R"),
  # 09l reports how the bound on the missing alive moves with the windows behind it
  list(name = "09l", file = "code/09l_bound_stability.R"),
  list(name = "10", file = "code/10_dists_tables.R"),
  list(name = "11", file = "code/11_estimation_conflict_allcause_mortality.R"),
  list(name = "12", file = "code/12_plotting_estimates.R"),
  list(name = "13", file = "code/13_sensitivity_migration_contributions.R"),
  list(name = "13b", file = "code/13b_sensitivity_missing_alive.R"),
  list(name = "13c", file = "code/13c_sensitivity_migration_range.R"),
  list(name = "13d", file = "code/13d_sensitivity_registration_lag.R"),
  list(name = "13e", file = "code/13e_sensitivity_donetsk_luhansk.R"),
  list(name = "13f", file = "code/13f_sensitivity_2022_timing.R"),
  list(name = "13g", file = "code/13g_sensitivity_baseline_pert.R"),
  list(name = "13h", file = "code/13h_sensitivity_civilian_age_profile.R"),
  list(name = "13i", file = "code/13i_sensitivity_donbas_fighters.R"),
  list(name = "13j", file = "code/13j_sensitivity_coherent_forecast.R"),
  # 13l checks itself against 13d, and 13m against 13e, so they follow them
  list(name = "13l", file = "code/13l_sensitivity_lag_tail.R"),
  list(name = "13m", file = "code/13m_sensitivity_base_coherent.R"),
  list(name = "13n", file = "code/13n_sensitivity_covid_carryover.R"),
  list(name = "13o", file = "code/13o_sensitivity_register_coverage.R"),
  list(name = "14", file = "code/14_e0_loss_decomposition_by_cause.R"),
  list(name = "14b", file = "code/14b_yll_45q15.R"),
  list(name = "14c", file = "code/14c_fixed_counterfactual.R"),
  list(name = "14d", file = "code/14d_drawn_population_base.R"),
  # 14e reads the draws (14) and every sensitivity step's output (13b-13n)
  list(name = "14e", file = "code/14e_structural_uncertainty.R"),
  # 15 builds every manuscript figure and table. It was missing from this
  # list, so a full pipeline run produced none of the deliverables.
  list(name = "15", file = "code/15_paper_figures_tables.R"),
  # 16 draws figure A0, the pipeline itself, from the release list and n_sim
  list(name = "16", file = "code/16_pipeline_figure.R")
)

for (step in steps_to_run) {
  if (file.exists(step$file)) {
    cat("\n[", format(Sys.time(), "%H:%M:%S"), "] Step", step$name, ":", step$file, "\n")
    tryCatch({
      source(step$file)
      cat("  ✓ Complete\n")
    }, error = function(e) {
      cat("  ✗ ERROR:", e$message, "\n")
      stop(e)
    })
  } else {
    cat("\nStep", step$name, "SKIPPED (file not found):", step$file, "\n")
  }
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n========================================\n")
cat("✓ Pipeline complete!\n")
cat("Total time:", round(as.numeric(elapsed), 1), "minutes\n")
cat("========================================\n")
