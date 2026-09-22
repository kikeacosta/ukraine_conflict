# ==============================================================================
# STEP 13j - A coherent counterfactual: Ukraine forecast with its neighbours
# ==============================================================================
#
# WHY
# ---
# The counterfactual is a Lee-Carter forecast fitted to Ukraine's own
# 2000-2019 rates (04). Its index for men moved erratically over those years -
# the crisis of 2005-2008 and the recovery after it - and a random walk fitted
# to twenty points carries that into a wide forecast band: from 2023 the
# forecast holds most of the variance of the male loss (table A39). A coherent
# forecast borrows strength from populations with similar mortality histories
# (Li and Lee 2005): a common age pattern of change, fitted to the pooled
# rates of the group, carries the trend, and Ukraine's own departure from it
# is forecast to fade. Its drift rests on a wider base, and its variance on
# more than Ukraine's own twenty years.
#
# THE MODEL (each sex on its own, ages 0-100, 2000-2019)
#   group        log m_G(x,t) = A(x) + B(x) K(t): Lee-Carter on the pooled
#                rates of Ukraine and eight EU neighbours - Poland, Slovakia,
#                Hungary, Romania, Bulgaria, Lithuania, Latvia, Estonia.
#                Belarus, Russia and Moldova are not in Eurostat's series. The
#                Human Mortality Database holds Belarus to 2018 and Russia to
#                2014, and not Moldova; since the pooled rate is the sum of
#                deaths over the sum of exposures, a country covering part of
#                the window alone would move K(t) by its own weight, so either
#                the pool or the fit window would have to be balanced first.
#   Ukraine      log m_U(x,t) = a_U(x) + B(x) K(t) + b(x) k(t), with b and k
#                the first singular component of what the common factor
#                leaves.
#   forecast     K a random walk with drift; k an AR(1) without constant, so
#                Ukraine's departure from the group decays; jump-off at the
#                fitted 2019 rates, as 04. 5,000 simulated paths give the mean
#                rates (the point forecast, as 04 reports the mean) and the
#                spread of the counterfactual e0.
# The point forecast is put through the deterministic projection at the
# central value of every conflict and migration input, as 13g does for its
# windows; the 95% band of the counterfactual e0 in 2025 is carried to the
# loss by projecting the paths at its 2.5th and 97.5th percentiles. Ages
# 0-100: the neighbours' last age is open (100+), as Ukraine's is.
#
# DATA  Eurostat demo_magec (deaths by age, completed years) and demo_pjan
#       (population on 1 January), 2000-2020, fetched once from Eurostat's
#       dissemination API into data_input/eurostat_neighbours/ and read from
#       there afterwards. Exposure is the mean of the populations on 1 January
#       of the year and the next.
#
# INPUTS   data_inter/ukr_life_tables_1989_2021.rds (03), data_input/DataDxEx.csv,
#          data_input/eurostat_neighbours/eurostat_deaths_population_2000_2020.csv,
#          the static inputs at the central values
# OUTPUTS  data_inter/ukr_coherent_forecast_e0.rds
# ==============================================================================

rm(list = ls())
gc()
source("code/00_setup.R")

mi <- mode_projection_inputs()
neighbours <- c(PL = "Poland", SK = "Slovakia", HU = "Hungary", RO = "Romania",
                BG = "Bulgaria", LT = "Lithuania", LV = "Latvia", EE = "Estonia")
fit_years <- 2000:2019

# 1. THE NEIGHBOURS' DEATHS AND POPULATION =====================================
# One request per dataset and country to Eurostat's JSON-stat API; the cells
# are laid out with the last dimension varying fastest.
eurostat_json <- function(dataset, geo, years) {
  url <- paste0("https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/", dataset,
                "?lang=EN&geo=", geo, "&sex=M&sex=F&sinceTimePeriod=", min(years),
                "&untilTimePeriod=", max(years))
  r <- httr::GET(url, httr::timeout(300))
  httr::stop_for_status(r)
  js <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  ids <- unlist(js$id)
  codes <- lapply(ids, \(d) {
    idx <- unlist(js$dimension[[d]]$category$index)
    names(idx)[order(idx)]
  })
  names(codes) <- ids
  grid <- expand.grid(rev(codes), stringsAsFactors = FALSE)[, ids]
  value <- rep(NA_real_, nrow(grid))
  value[as.integer(names(js$value)) + 1] <- as.numeric(unlist(js$value))
  as_tibble(grid) |> mutate(value = value)
}
age_of <- function(code) {
  single <- str_detect(code, "^Y[0-9]+$")
  out <- rep(NA_real_, length(code))
  out[single] <- as.numeric(sub("Y", "", code[single]))
  out[code == "Y_LT1"] <- 0
  out[code == "Y_OPEN"] <- 100
  out
}
nb_file <- "data_input/eurostat_neighbours/eurostat_deaths_population_2000_2020.csv"
if (!file.exists(nb_file)) {
  dir.create(dirname(nb_file), showWarnings = FALSE, recursive = TRUE)
  fetched <- map_dfr(names(neighbours), \(g) {
    bind_rows(
      eurostat_json("demo_magec", g, 2000:2020) |> mutate(measure = "deaths"),
      eurostat_json("demo_pjan", g, 2000:2020) |> mutate(measure = "population")
    )
  })
  nb_wide <-
    fetched |>
    mutate(age = age_of(age), year = as.integer(time)) |>
    filter(!is.na(age), sex %in% c("M", "F")) |>
    select(geo, measure, sex, year, age, value) |>
    arrange(geo, measure, sex, year, age) |>
    pivot_wider(names_from = age, values_from = value, names_prefix = "age")
  write_csv(nb_wide, nb_file)
}
nb <-
  read_csv(nb_file, show_col_types = FALSE) |>
  pivot_longer(starts_with("age"), names_to = "age", values_to = "value") |>
  mutate(age = as.numeric(sub("age", "", age)), sex = if_else(sex == "M", "m", "f"))
stopifnot(setequal(unique(nb$geo), names(neighbours)))

# exposure: the mean of the populations on 1 January of the year and the next.
# Some early years lack cells: Latvia's deaths in 2000-2001, and single ages
# from 85 or 90 in a few countries' early years, where the open age group
# starts earlier and so holds more than the 100+ it is read as. A country-year
# with any cell missing at 85-99 is left out of the pooled rates from 85 up,
# and any other missing cell is left out on its own. Every cell of the pool
# keeps Ukraine, the largest population in it; below 85 it keeps seven or all
# eight neighbours, and from 85 up before 2013 fewer (Hungary's single ages
# start in 2013, and in 2000-2001 only two or three neighbours remain).
nb_dx <-
  nb |>
  pivot_wider(names_from = measure, values_from = value) |>
  arrange(geo, sex, age, year) |>
  mutate(exposure = (population + lead(population)) / 2, .by = c(geo, sex, age)) |>
  filter(year %in% fit_years) |>
  mutate(old_complete = all(!is.na(deaths[age %in% 85:99]) & !is.na(exposure[age %in% 85:99])),
         .by = c(geo, sex, year)) |>
  filter(old_complete | age < 85, !is.na(deaths), !is.na(exposure)) |>
  select(geo, sex, age, year, deaths, exposure)

# Countries from the Human Mortality Database, when 13k has fetched them: the
# same columns, already closed at 100. A country that covers part of the window
# only moves the pooled rate by its own weight, so one that does is taken out
# here rather than allowed to bend the common index.
hmd_file <- "data_input/hmd_neighbours/hmd_deaths_exposures.csv"
if (file.exists(hmd_file)) {
  hmd <- read_csv(hmd_file, show_col_types = FALSE) |>
    filter(year %in% fit_years) |>
    select(geo, sex, age, year, deaths, exposure)
  span <- hmd |> summarise(last = max(year), n = n_distinct(year), .by = geo)
  short <- span$geo[span$n < length(fit_years)]
  if (length(short))
    message("  left out of the pool, covering part of ", min(fit_years), "-",
            max(fit_years), " only: ",
            paste0(span$geo[match(short, span$geo)], " (to ",
                   span$last[match(short, span$geo)], ")", collapse = ", "),
            ". Ending fit_years at that year would let it in, at the cost of a ",
            "jump-off earlier than 04's.")
  nb_dx <- bind_rows(nb_dx, hmd |> filter(!geo %in% short))
  message("  pool: ", paste(sort(unique(nb_dx$geo)), collapse = ", "))
}
stopifnot(all(nb_dx$exposure >= 0), all(nb_dx$deaths >= 0))
pool_size <- nb_dx |> count(sex, age, year)
stopifnot(nrow(pool_size) == 2 * 101 * length(fit_years), min(pool_size$n) >= 2,
          min(pool_size$n[pool_size$age < 85]) >= 7)

# Ukraine, as 04 and 13g read it: the rates of the regions with complete
# registration and their exposures
ukr_exp <-
  read.csv("data_input/DataDxEx.csv", header = TRUE) |>
  rename_with(tolower) |>
  filter(data == "Ex") |>
  pivot_longer(starts_with("age"), names_to = "age", values_to = "exposure") |>
  mutate(age = as.numeric(gsub("age", "", age)), sex = str_sub(sex, 1, 1)) |>
  select(year, sex, age, exposure)
ukr_dx <-
  read_rds("data_inter/ukr_life_tables_1989_2021.rds") |>
  mutate(sex = str_sub(sex, 1, 1)) |>
  select(year, sex, age, mx) |>
  filter(year %in% fit_years) |>
  left_join(ukr_exp, by = c("year", "sex", "age")) |>
  mutate(deaths = mx * exposure, geo = "UA")
stopifnot(!anyNA(ukr_dx$exposure), nrow(ukr_dx) == 2 * 101 * length(fit_years))

# 2. THE COHERENT MODEL ========================================================
# matrices ages x years of log rates
log_rates <- function(d) {
  d |>
    select(age, year, lm) |>
    arrange(year, age) |>
    pivot_wider(names_from = year, values_from = lm) |>
    arrange(age) |>
    select(-age) |>
    as.matrix()
}
first_component <- function(M) {
  s <- svd(M)
  u <- s$u[, 1]
  list(load = u / sum(u), index = s$d[1] * s$v[, 1] * sum(u))
}
fit_coherent <- function(s) {
  group <-
    bind_rows(nb_dx |> filter(sex == s), ukr_dx |> filter(sex == s) |> select(geo, sex, age, year, deaths, exposure)) |>
    summarise(lm = log(sum(deaths) / sum(exposure)), .by = c(age, year))
  stopifnot(all(is.finite(group$lm)))
  G <- log_rates(group)
  A <- rowMeans(G)
  common <- first_component(G - A)
  U <- log_rates(ukr_dx |> filter(sex == s) |> mutate(lm = log(mx)))
  a_u <- rowMeans(U)
  own <- first_component(U - a_u - outer(common$load, common$index))
  K <- common$index
  k <- own$index
  n <- length(K)
  drift <- (K[n] - K[1]) / (n - 1)
  sigma_K <- sqrt(sum((diff(K) - drift)^2) / (n - 2))
  phi <- sum(k[-1] * k[-n]) / sum(k[-n]^2)
  sigma_k <- sqrt(sum((k[-1] - phi * k[-n])^2) / (n - 2))
  list(sex = s, a = a_u, B = common$load, K = K, b = own$load, k = k,
       drift = drift, sigma_K = sigma_K, phi = phi, sigma_k = sigma_k)
}
fits <- map(c(f = "f", m = "m"), fit_coherent)
stopifnot(all(map_dbl(fits, "phi") < 1))

# 5,000 paths of the two indices, 2020-2025, with the uncertainty of the
# drift as a random walk's forecast carries it (h * sigma^2 * (1 + h / (n - 1)))
set.seed(4242)
n_paths <- 5000
horizon <- 6
years_fc <- 2019 + seq_len(horizon)
simulate_rates <- function(ft) {
  n <- length(ft$K)
  drift_draw <- rnorm(n_paths, ft$drift, ft$sigma_K / sqrt(n - 1))
  eps_K <- matrix(rnorm(n_paths * horizon, 0, ft$sigma_K), n_paths)
  eps_k <- matrix(rnorm(n_paths * horizon, 0, ft$sigma_k), n_paths)
  K_path <- ft$K[n] + outer(drift_draw, seq_len(horizon)) + t(apply(eps_K, 1, cumsum))
  k_path <- matrix(NA_real_, n_paths, horizon)
  prev <- rep(ft$k[n], n_paths)
  for (h in seq_len(horizon)) {
    prev <- ft$phi * prev + eps_k[, h]
    k_path[, h] <- prev
  }
  # rates by age (rows) and path (columns), one matrix per forecast year
  map(seq_len(horizon), \(h) exp(ft$a + outer(ft$B, K_path[, h]) + outer(ft$b, k_path[, h])))
}
paths <- map(fits, simulate_rates)

# the point forecast: the mean rate over the paths, as 04 reports the mean
point_fc <-
  imap_dfr(paths, \(pl, s) map_dfr(seq_len(horizon), \(h) {
    tibble(year = years_fc[h], sex = s, age = 0:100, mx = rowMeans(pl[[h]]))
  })) |>
  filter(year %in% 2022:2025)

# the counterfactual e0 of every path in 2025, and the paths at its 2.5th and
# 97.5th percentiles
h25 <- which(years_fc == 2025)
e0_paths <- imap(paths, \(pl, s) lt_cols(pl[[h25]], rep(s, n_paths))$ex[1, ])
path_at <- function(s, p) which.min(abs(e0_paths[[s]] - quantile(e0_paths[[s]], p)))
fc_path <- function(p) {
  map_dfr(c("f", "m"), \(s) {
    j <- path_at(s, p)
    map_dfr(seq_len(horizon), \(h) tibble(year = years_fc[h], sex = s, age = 0:100, mx = paths[[s]][[h]][, j]))
  }) |>
    filter(year %in% 2022:2025)
}

# 3. THE LOSS WITH THE COHERENT COUNTERFACTUAL =================================
with_mx <- function(fc) mi$static_inputs |> select(-mx) |> left_join(fc, by = c("year", "sex", "age"))
label <- "2000-2019, coherent with eight neighbours (Li-Lee)"
project <- function(fc) {
  st <- with_mx(fc)
  stopifnot(all(st$mx > 0), !any(is.na(st$mx)))
  loss_at_mode(mi$draws, st, mi$pop22_ini)
}
baseline <- project(point_fc) |> mutate(window = label, first = 2000, last = 2019, .before = 1)
spread <-
  map_dfr(c(lower = 0.025, upper = 0.975), \(p) project(fc_path(p)), .id = "bound") |>
  mutate(window = label, first = 2000, last = 2019, .before = 1)

# the width of the forecast, against the Lee-Carter forecast of 04 as 11 draws
# it: the standard deviation of the counterfactual e0 in 2025
sd_e0 <- tibble(sex = c("f", "m"), sd_e0_2025 = c(sd(e0_paths$f), sd(e0_paths$m)))
params <-
  map_dfr(fits, \(ft) tibble(sex = ft$sex, drift = ft$drift, sigma_K = ft$sigma_K, phi = ft$phi,
                             sigma_k = ft$sigma_k)) |>
  left_join(sd_e0, by = "sex")
write_rds(list(baseline = baseline, spread = spread, params = params,
               neighbours = neighbours, fit_years = fit_years),
          "data_inter/ukr_coherent_forecast_e0.rds")

cat("\n=== COHERENT COUNTERFACTUAL: PARAMETERS ===\n")
print(as.data.frame(params |> mutate(across(where(is.numeric), \(x) round(x, 4)))))
cat("\n=== COHERENT COUNTERFACTUAL: e0 AND LOSS, 2025 ===\n")
print(as.data.frame(
  bind_rows(baseline |> mutate(bound = "point"), spread) |>
    filter(year == 2025) |>
    mutate(across(c(e0_bsn, loss), \(x) round(x, 2))) |>
    select(bound, sex, e0_bsn, loss) |>
    pivot_wider(names_from = sex, values_from = c(e0_bsn, loss))
))
