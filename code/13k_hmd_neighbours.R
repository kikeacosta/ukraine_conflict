# ==============================================================================
# STEP 13k - Belarus (and Russia) from the Human Mortality Database
# ==============================================================================
#
# WHY
# ---
# The coherent counterfactual of 13j pools Ukraine with eight EU neighbours from
# Eurostat. The post-Soviet neighbours closest to Ukraine's mortality history are
# not in that series. The Human Mortality Database holds Belarus to 2018 and
# Russia to 2014, and does not hold Moldova.
#
# Neither reaches 2019. The pooled rate is the sum of deaths over the sum of
# exposures, so a country that covers part of the window alone moves the common
# index by its own weight rather than by mortality - and for Belarus the year it
# is missing is the jump-off year itself, where a kink in the index bends the
# drift the forecast rests on. 13j therefore leaves out any country that does
# not cover the whole of fit_years, and says so.
#
# To use what this step fetches, end fit_years where the country's series ends:
# 2018 for Belarus, 2014 for Russia. That buys a wider base for the trend at the
# price of a jump-off earlier than 04's, so it is a choice about the design, not
# a repair - which is why this step only fetches, and changes nothing.
#
# CREDENTIALS
# -----------
# The HMD needs an account. This step reads it from the environment and never
# stores it:
#
#   HMD_USER=you@example.org
#   HMD_PASS=...
#
# Put those two lines in ~/.Renviron (outside this repository), restart R, and
# run the step. Nothing is written to the repository but the counts themselves.
#
# LICENSING
# ---------
# The HMD user agreement governs what may be redistributed. Until that is
# checked, the cache is written to a gitignored folder, unlike the Eurostat one.
#
# INPUTS   the HMD web files Deaths_1x1 and Exposures_1x1, by country
# OUTPUTS  data_input/hmd_neighbours/hmd_deaths_exposures.csv (gitignored):
#          geo, sex, age (0-100, last open), year, deaths, exposure
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

countries <- c(BLR = "Belarus")      # add RUS = "Russia" once the pool is balanced
fit_years <- 2000:2019
out_dir <- "data_input/hmd_neighbours"
out_csv <- file.path(out_dir, "hmd_deaths_exposures.csv")

user <- Sys.getenv("HMD_USER")
pass <- Sys.getenv("HMD_PASS")
if (!nzchar(user) || !nzchar(pass))
  stop("Set HMD_USER and HMD_PASS in ~/.Renviron (see the header of this file).")

# 1. FETCH =====================================================================
# The HMD serves one text file per country and statistic, behind a form login.
# One session, then one request per file.
login <- request("https://www.mortality.org/Account/Login") |>
  req_body_form(Email = user, Password = pass) |>
  req_error(is_error = \(resp) FALSE) |>
  req_perform()
if (resp_status(login) >= 400)
  stop("The HMD login failed (", resp_status(login), "). Check HMD_USER and HMD_PASS.")

hmd_file <- function(cntry, item) {
  url <- paste0("https://www.mortality.org/File/GetDocument/hmd.v6/", cntry,
                "/STATS/", item, ".txt")
  txt <- request(url) |>
    req_cookie_preserve(login$cookies) |>
    req_perform() |>
    resp_body_string()
  if (grepl("<html", substr(txt, 1, 200), ignore.case = TRUE))
    stop("The HMD returned a page rather than ", item, " for ", cntry,
         ": the session is not authenticated.")
  # the files carry two header lines before the column names
  read_table(I(txt), skip = 2, na = ".", show_col_types = FALSE)
}

# 2. TIDY ======================================================================
# Single ages 0-110+ in the HMD, closed at 100 here as the rest of the pool is.
tidy_one <- function(df, value) {
  df |>
    mutate(Age = as.integer(str_remove(Age, "\\+")), Year = as.integer(Year)) |>
    filter(Year %in% fit_years, !is.na(Age)) |>
    pivot_longer(c(Female, Male), names_to = "sex", values_to = "value") |>
    mutate(sex = if_else(sex == "Female", "f", "m"),
           age = pmin(Age, 100)) |>
    summarise("{value}" := sum(value, na.rm = TRUE), .by = c(sex, age, year = Year))
}

hmd <- imap(countries, function(name, code) {
  message("  ", name, " ...")
  dx <- tidy_one(hmd_file(code, "Deaths_1x1"), "deaths")
  ex <- tidy_one(hmd_file(code, "Exposures_1x1"), "exposure")
  full_join(dx, ex, by = c("sex", "age", "year")) |> mutate(geo = code, .before = 1)
}) |> list_rbind()

stopifnot(nrow(hmd) > 0, all(hmd$deaths >= 0), all(hmd$exposure >= 0),
          n_distinct(hmd$age) == 101)

# 3. WRITE =====================================================================
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(hmd, out_csv)
cat("\nWritten:", out_csv, "\n")
hmd |> summarise(years = paste(min(year), max(year), sep = "-"),
                 deaths = sum(deaths), .by = geo) |> print()
cat("\n13j reads this file when it is present.\n")
