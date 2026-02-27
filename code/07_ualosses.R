rm(list = ls())
source("code/00_setup.R")

dt <-
  read_xlsx("data_input/UKR_ualosses_Personnel.xlsx", sheet = "Database")

dt2 <-
  dt |>
  rename_with(tolower) |>
  mutate(
    age = round(age),
    sex = case_when(
      sex == 0 ~ "f",
      sex == 1 ~ "m",
      TRUE ~ NA_character_
    ),
    date_bth = ymd("1900-01-01") + as.numeric(datebirth),
    date_evnt = ymd("1900-01-01") + as.numeric(dateevent)
  )

identify_sex_from_names <- function(patronym, firstname) {
  # normalize inputs
  p <- tolower(stringi::stri_trans_general(
    trimws(as.character(patronym)),
    "Latin-ASCII"
  ))
  f <- tolower(stringi::stri_trans_general(
    trimws(as.character(firstname)),
    "Latin-ASCII"
  ))
  n <- max(length(p), length(f))
  p <- rep_len(p, n)
  f <- rep_len(f, n)

  sex <- rep(NA_character_, n)
  confidence <- rep(NA_character_, n)
  reason <- rep(NA_character_, n)

  # 1. ENHANCED PATRONYMIC PATTERNS
  # Added 'ich' explicitly to catch 'Illich'
  male_pat <- "(ovych|vych|ovich|vich|ich)$"
  female_pat <- "(ivna|yivna|ovna|vna)$"

  is_male_pat <- !is.na(p) & p != "" & p != "na" & str_detect(p, male_pat)
  is_female_pat <- !is.na(p) & p != "" & p != "na" & str_detect(p, female_pat)

  sex[is_male_pat] <- "m"
  confidence[is_male_pat] <- "high"
  reason[is_male_pat] <- "patronymic"

  sex[is_female_pat] <- "f"
  confidence[is_female_pat] <- "high"
  reason[is_female_pat] <- "patronymic"

  # 2. UPDATED FIRST-NAME LISTS
  # Added your specific male exceptions (Illja, Mykyta, Jarema, Zaza, etc.)
  male_names <- c(
    "andriy",
    "andrii",
    "oleksandr",
    "oleksander",
    "vitalii",
    "vitaliy",
    "vladimir",
    "volodymyr",
    "serhiy",
    "serhii",
    "ivan",
    "mykola",
    "dmytro",
    "artem",
    "oleksiy",
    "yuriy",
    "roman",
    "vasyl",
    "anatoliy",
    "oleg",
    "petro",
    "illja",
    "illya",
    "ylja",
    "mykyta",
    "jarema",
    "zaza",
    "mamuka",
    "kakha"
  )

  female_names <- c(
    "olena",
    "kateryna",
    "anastasiya",
    "iryna",
    "irina",
    "nataliya",
    "marina",
    "olga",
    "sophia",
    "sofiya",
    "yana",
    "tetyana",
    "tetiana",
    "lyudmyla",
    "inna",
    "valeriya",
    "oleksandra"
  )

  is_male_first <- is.na(sex) & !is.na(f) & f %in% male_names
  is_female_first <- is.na(sex) & !is.na(f) & f %in% female_names

  sex[is_male_first] <- "m"
  confidence[is_male_first] <- "medium"
  reason[is_male_first] <- "first_name_list"

  sex[is_female_first] <- "f"
  confidence[is_female_first] <- "medium"
  reason[is_female_first] <- "first_name_list"

  # 3. SUFFIX HEURISTICS (with exception handling)
  rem <- which(is.na(sex))
  if (length(rem) > 0) {
    fname_rem <- f[rem]

    # Heuristic: female if ends in 'a' or 'ya', UNLESS it's a known male exception handled above
    is_f_suffix <- !is.na(fname_rem) & str_detect(fname_rem, "(a|ya)$")
    is_m_suffix <- !is.na(fname_rem) & !is_f_suffix & fname_rem != ""

    sex[rem[is_f_suffix]] <- "f"
    confidence[rem[is_f_suffix]] <- "low"
    reason[rem[is_f_suffix]] <- "first_name_suffix"

    sex[rem[is_m_suffix]] <- "m"
    confidence[rem[is_m_suffix]] <- "low"
    reason[rem[is_m_suffix]] <- "first_name_suffix"
  }

  # Final Cleanup
  unk <- which(is.na(sex) | sex == "")
  if (length(unk) > 0) {
    sex[unk] <- "Unknown"
    confidence[unk] <- "unknown"
    reason[unk] <- "no_data"
  }

  tibble::tibble(
    sex_pred = sex,
    sex_confidence = confidence,
    sex_reason = reason
  )
}

# Apply to dt3 (assuming dt3 has columns 'patronym' and 'firstname')
dt3 <-
  dt2 |>
  mutate(
    # create normalized lower-case helper columns only if you want them; not required
    across(c(patronym, firstname), ~.x),
    # bind the prediction columns
    bind_cols(identify_sex_from_names(patronym, firstname)),
    tst = sex_pred == sex
  )

dt3 |>
  filter(!tst) |>
  select(lastname, patronym, firstname) |>
  print(n = 25)

dt4 <-
  dt3 |>
  filter(nationality == "Ukraine") %>%
  mutate(
    year = year(date_evnt),
    age = floor(interval(date_bth, date_evnt) / years(1))
  ) |>
  filter(nationality == "Ukraine", year %in% 2022:2025) %>%
  select(status, year, sex = sex_pred, age)

dts <-
  dt4 %>%
  summarise(dx = n(), .by = c(status, year, sex, age)) |>
  arrange(status, year, sex, age, dx)

# Imputing age and year of the event
dts2 <-
  dts |>
  group_by(status, year, sex) |>
  mutate(dts_tot = sum(dx)) |>
  filter(!is.na(age)) |>
  mutate(dts_sum = sum(dx)) |>
  ungroup() |>
  mutate(dx = dx * dts_tot / dts_sum) |>
  group_by(status, sex, age) |>
  mutate(dts_tot = sum(dx)) |>
  filter(!is.na(year)) |>
  mutate(dts_sum = sum(dx)) |>
  ungroup() |>
  select(-dts_tot, -dts_sum)

unique(dts2$age)
unique(dts2$year)

unique(dts2$age) |> sort()

dts3 <-
  dts2 |>
  complete(status, year, sex, age, fill = list(dx = 0))

dts3 |>
  summarise(dx = sum(dx), .by = c(year))

dt2 |>
  mutate(year = year(date_evnt)) |>
  filter(nationality == "Ukraine", year %in% 2022:2025) %>%
  summarise(dx = n(), .by = c(year)) |>
  arrange(year)

dts3 |>
  filter(year %in% 2022:2024, status != "prisioner") |>
  summarise(dx = sum(dx))

write_rds(dts3, "data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")
