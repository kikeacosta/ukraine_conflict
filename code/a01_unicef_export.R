# ==============================================================================
# STEP A01 - UNICEF export: all 5,000 iterations, by cause of death
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Writes the full simulation output in the flat format UNICEF consumes: one
# row per year x sex x age group x cause x iteration, carrying the death
# count (dx) and the death rate (mx).
#
# CAUSES OF DEATH
# ---------------
# The previous export (260720) carried three causes: all, conflict, expected.
# This version splits the conflict deaths into their three components, so the
# five causes are:
#
#   expected             deaths that would have occurred without the war
#                        (Lee-Carter counterfactual, see 04)
#   civilians            civilian conflict deaths
#                        (UCDP totals, OHCHR age-sex profile)
#   combatants_reported  combatant deaths individually confirmed in the
#                        ualosses register
#   combatants_imputed   additional combatant deaths imputed from the
#                        personnel recorded as missing (see 09)
#   all                  all-cause = expected + the three above
#
# The three conflict components sum exactly to the old "conflict" cause, and
# expected + conflict sums exactly to "all". Both identities are asserted.
#
# AGE GROUPS
# ----------
# The pipeline itself works in SINGLE YEARS of age throughout, and the stored
# estimates stay that way. The abridging below happens ONLY in this export,
# because that is the format UNICEF receives: 0, 1, 5, 10, ... 75, 80+
# (18 groups), matching the previous file.
#
# Deaths and exposure are summed within a group and mx is recomputed as
# dx / exposure — NOT averaged across the single-year rates, which would not
# be consistent with the counts.
#
# INPUT    data_inter/ukr_sim_draws_2022_2025_n<n_sim>.rds   (from 11)
# OUTPUT   data_inter/unicef/<yymmdd>_ukraine_mx_estimates_2022_2025_
#            conflict_and_allcause.csv  (+ .zip)
#          Too large for git; excluded in .gitignore.
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

# named for the simulation size set in 00_setup.R, so this reads the draws that
# match the configured n_sim rather than whatever a previous run left behind
draws_file <- sprintf("data_inter/ukr_sim_draws_2022_2025_n%d.rds", n_sim)
if (!file.exists(draws_file)) {
  stop(
    "Simulation draws not found at n_sim = ", n_sim, ". Run step 11 first:\n  ",
    draws_file,
    call. = FALSE
  )
}

out_dir <- "data_inter/unicef"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.Date(), "%y%m%d")
base <- paste0(
  stamp, "_ukraine_mx_estimates_2022_2025_conflict_and_allcause"
)
csv_path <- file.path(out_dir, paste0(base, ".csv"))
zip_path <- file.path(out_dir, paste0(base, ".zip"))

# 1. LOAD AND GROUP AGES =======================================================
d <- as.data.table(readRDS(draws_file))

# abridged groups: 0, 1, 5, 10, ..., 75, 80+
brks <- c(0, 1, seq(5, 80, 5), Inf)
labs <- c(0, 1, seq(5, 80, 5))
d[, age_grp := as.integer(as.character(
  cut(age, breaks = brks, labels = labs, right = FALSE)
))]

stopifnot(
  !any(is.na(d$age_grp)),
  length(unique(d$age_grp)) == 18
)

# 2. AGGREGATE WITHIN AGE GROUPS ===============================================
# Exposure and each death component are summed; rates are derived afterwards.
agg <- d[
  ,
  .(
    pop = sum(pop),
    expected = sum(expected),
    civilians = sum(civilian),
    combatants_reported = sum(combatant_confirmed),
    combatants_imputed = sum(combatant_imputed)
  ),
  by = .(year, sex, age = age_grp, it = sim_id)
]

agg[, all := expected + civilians + combatants_reported + combatants_imputed]

# 3. RESHAPE TO THE EXPORT LAYOUT ==============================================
CAUSES <- c(
  "expected", "civilians", "combatants_reported", "combatants_imputed", "all"
)

out <- melt(
  agg,
  id.vars = c("year", "sex", "age", "it", "pop"),
  measure.vars = CAUSES,
  variable.name = "cause",
  value.name = "dx"
)
out[, cause := as.character(cause)]
out[, mx := dx / pop]

# Exposure is dropped: it is identical across causes within a cell, so
# carrying it would add ~40 MB to the file for no new information. The column
# set therefore matches the previous export exactly. Where a cause has
# dx > 0, exposure is recoverable as dx / mx.
setcolorder(out, c("year", "sex", "age", "cause", "it", "dx", "mx"))
out[, pop := NULL]
setorder(out, year, sex, age, cause, it)

# 4. CHECKS ====================================================================
n_exp <- length(unique(out$year)) * length(unique(out$sex)) *
  18L * length(CAUSES) * length(unique(out$it))

cat("rows:", nrow(out), " expected:", n_exp, "\n")
stopifnot(nrow(out) == n_exp, !any(is.na(out$dx)), !any(is.na(out$mx)))

# the components must reconstruct the all-cause total exactly
chk <- dcast(out, year + sex + age + it ~ cause, value.var = "dx")
gap <- max(abs(
  chk$all - (chk$expected + chk$civilians +
    chk$combatants_reported + chk$combatants_imputed)
))
cat("max |all - sum(components)| :", gap, "\n")
stopifnot(gap < 1e-9)

cat("\ntotal deaths by cause, summed over all iterations / 5000:\n")
print(out[, .(deaths = sum(dx) / length(unique(out$it))), by = cause])

# 5. WRITE =====================================================================
fwrite(out, csv_path)
cat(
  "\nwritten:", csv_path,
  "(", round(file.size(csv_path) / 1e6, 1), "MB )\n"
)

# The zip package writes archives from R itself, so this does not depend on a
# zip.exe being on PATH (there is none on Windows without Rtools on PATH).
# Falls back to gzip, which fwrite does natively, if that is unavailable too.
# zipr() stores the file at the top level of the archive without needing a
# `root` argument, so every path here stays relative to the project directory.
zipped <- FALSE
if (requireNamespace("zip", quietly = TRUE)) {
  zip::zipr(
    zipfile = zip_path,
    files = csv_path,
    compression_level = 9
  )
  zipped <- file.exists(zip_path)
}

if (isTRUE(zipped)) {
  cat(
    "written:", zip_path,
    "(", round(file.size(zip_path) / 1e6, 1), "MB )\n"
  )
} else {
  gz_path <- paste0(csv_path, ".gz")
  fwrite(out, gz_path, compress = "gzip")
  cat(
    "zip package unavailable; wrote", gz_path,
    "(", round(file.size(gz_path) / 1e6, 1), "MB ) instead\n"
  )
}

message("Done.")
