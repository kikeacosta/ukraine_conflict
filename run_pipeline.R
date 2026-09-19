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
  list(name = "09", file = "code/09_ualosses_missing_analysis.R"),
  list(name = "10", file = "code/10_dists_tables.R"),
  list(name = "11", file = "code/11_estimation_conflict_allcause_mortality.R"),
  list(name = "12", file = "code/12_plotting_estimates.R"),
  list(name = "13", file = "code/13_sensitivity_migration_contributions.R"),
  list(name = "13b", file = "code/13b_sensitivity_alpha_missing.R"),
  list(name = "13c", file = "code/13c_sensitivity_migration_range.R"),
  list(name = "14", file = "code/14_e0_loss_decomposition_by_cause.R"),
  # 15 builds every manuscript figure and table. It was missing from this
  # list, so a full pipeline run produced none of the deliverables.
  list(name = "15", file = "code/15_paper_figures_tables.R")
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
