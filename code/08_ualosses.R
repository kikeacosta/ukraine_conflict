rm(list = ls())
source("code/00_setup.R")

# data compiled by Olivier Hubert
# in https://www.kaggle.com/datasets/ol4ubert/confirmed-ukrainian-military-personnel-losses
#
# The register is ~31 MB of individual records and is not tracked in git.
# Everything below the cache boundary is aggregated to status-year-sex-age
# counts before being written out, so the cached extract is a few kB AND
# carries no personal data (no names, no dates of birth).

# The most recent release, v19: published on 21 July 2026, with events to
# 18 June 2026 (ual_releases, 00_setup.R). Every age-sex count comes
# from it; the earlier releases are read only by 09, to follow the missing.
ual_file <- "data_input/ualosses_hubert_datasets/260919_UKR_ualosses_Personnel_v19.xlsx"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sex from Ukrainian patronymics / first names.
# Only used as a FALLBACK: the register records sex directly for ~99.5% of
# entries, so the heuristic is applied to the residual few hundred only.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  male_names <- c(
    "andriy", "andrii", "oleksandr", "oleksander", "vitalii", "vitaliy",
    "vladimir", "volodymyr", "serhiy", "serhii", "ivan", "mykola", "dmytro",
    "artem", "oleksiy", "yuriy", "roman", "vasyl", "anatoliy", "oleg", "oleh",
    "petro", "illja", "illya", "ylja", "mykyta", "jarema", "zaza", "mamuka",
    "kakha"
  )

  female_names <- c(
    "olena", "kateryna", "anastasiya", "iryna", "irina", "nataliya", "marina",
    "olga", "sophia", "sofiya", "yana", "tetyana", "tetiana", "lyudmyla",
    "inna", "valeriya", "oleksandra"
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

    # Heuristic: female if ends in 'a' or 'ya', UNLESS a known male exception
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
    sex[unk] <- NA_character_
    confidence[unk] <- "unknown"
    reason[unk] <- "no_data"
  }

  tibble::tibble(
    sex_pred = sex,
    sex_confidence = confidence,
    sex_reason = reason
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# raw register -> anonymised status/year/sex/age counts (cached)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dts <- cache_rds("data_inter/ualosses_counts_2022_2025.rds", {
  dt2 <-
    read_xlsx(require_raw(ual_file), sheet = "Database") |>
    rename_with(tolower) |>
    mutate(
      sex_rec = case_when(
        sex == 0 ~ "f",
        sex == 1 ~ "m",
        TRUE ~ NA_character_
      ),
      date_bth = excel_date(datebirth),
      date_evnt = excel_date(dateevent)
    )

  dt3 <-
    dt2 |>
    bind_cols(identify_sex_from_names(dt2$patronym, dt2$firstname)) |>
    # the recorded field wins; the name heuristic only fills the gaps
    mutate(sex = coalesce(sex_rec, sex_pred))

  message(
    "  sex recorded: ", sum(!is.na(dt3$sex_rec)),
    " | filled from names: ", sum(is.na(dt3$sex_rec) & !is.na(dt3$sex_pred)),
    " | still unknown: ", sum(is.na(dt3$sex)),
    " | heuristic agrees with record: ",
    round(100 * mean(dt3$sex_pred == dt3$sex_rec, na.rm = TRUE), 1), "%"
  )

  dt3 |>
    # ual_is_ukrainian() (00_setup.R) also keeps the records whose nationality
    # field holds a Ukrainian place of origin instead of a country
    filter(ual_is_ukrainian(nationality), !is.na(sex)) |>
    mutate(
      year = year(date_evnt),
      age = floor(interval(date_bth, date_evnt) / years(1))
    ) |>
    summarise(dx = n(), .by = c(status, year, sex, age)) |>
    arrange(status, year, sex, age)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# raw register -> anonymised counts by place of origin (cached)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The From field gives a place of origin as "settlement, community, district,
# oblast", the city of Kyiv as "Kyiv" alone, or "Unknown". The register does
# not say whether it is the place of birth or of residence. The oblast is the
# last part of the field. The Ukrainian dead and missing of 2022-2025 are
# counted by event year, status and origin group: Donetsk, Luhansk, another
# oblast, Kyiv city, unknown. 13q takes from these counts the share of the
# combatants who come from outside Donetsk and Luhansk. Only the counts are
# written: no names, dates of birth or individual records.
origin <- cache_rds("data_inter/ualosses_origin_2022_2025.rds", {
  ori <-
    read_xlsx(require_raw(ual_file), sheet = "Database", col_types = "text") |>
    transmute(
      status = Status,
      year = as.integer(year(excel_date(DateEvent))),
      ukrainian = ual_is_ukrainian(Nationality),
      place = str_trim(str_extract(From, "[^,]*$"))
    ) |>
    filter(ukrainian, status %in% c("dead", "missing"), year %in% 2022:2025)
  # every value is an oblast, Kyiv city or unknown: a place of another kind
  # (Crimea, Sevastopol) would otherwise be counted as unknown without notice
  stopifnot(all(str_detect(ori$place, " Oblast$") | ori$place %in% c("Kyiv", "Unknown") |
                  is.na(ori$place)))
  ori |>
    mutate(origin = case_when(
      place == "Donetsk Oblast" ~ "Donetsk",
      place == "Luhansk Oblast" ~ "Luhansk",
      place == "Kyiv" ~ "Kyiv city",
      str_detect(place, " Oblast$") ~ "another oblast",
      .default = "unknown"
    )) |>
    count(year, status, origin, name = "n") |>
    complete(year = 2022:2025, status = c("dead", "missing"),
             origin = c("Donetsk", "Luhansk", "another oblast", "Kyiv city", "unknown"),
             fill = list(n = 0L)) |>
    arrange(year, status, origin)
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# imputing unknown age and unknown year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Two stages, in this order, because a record can be missing either field or
# both. The previous version filtered on !is.na(age) BEFORE handling the year,
# which silently deleted every record missing both (~2.7% of the register)
# instead of redistributing it, and then omitted the second rescaling step
# altogether.
known <- dts |> filter(!is.na(year), !is.na(age))
no_age <- dts |> filter(!is.na(year), is.na(age))
no_year <- dts |> filter(is.na(year))

# (1) known year, unknown age -> over the ages of the same status-sex-year
step1 <-
  known |>
  left_join(
    no_age |> summarise(add = sum(dx), .by = c(status, sex, year)),
    by = c("status", "sex", "year")
  ) |>
  replace_na(list(add = 0)) |>
  mutate(dx = dx + add * dx / sum(dx), .by = c(status, sex, year)) |>
  select(-add)

# (2) unknown year (these are also missing age) -> over the whole known
#     year-age distribution of the same status-sex
dts2 <-
  step1 |>
  left_join(
    no_year |> summarise(add = sum(dx), .by = c(status, sex)),
    by = c("status", "sex")
  ) |>
  replace_na(list(add = 0)) |>
  mutate(dx = dx + add * dx / sum(dx), .by = c(status, sex)) |>
  select(-add)

# the imputation must be mass-preserving
cat(
  "records in:", sum(dts$dx),
  "| after imputation:", round(sum(dts2$dx)),
  "| difference:", round(sum(dts$dx) - sum(dts2$dx), 6), "\n"
)
stopifnot(abs(sum(dts2$dx) - sum(dts$dx)) < 1e-6)

dts3 <-
  dts2 |>
  filter(year %in% 2022:2025) %>%
  complete(status, year, sex, age, fill = list(dx = 0)) |>
  # released_prisoner (v19 onwards) is someone returned from captivity: alive,
  # neither a death nor missing, and kept only so the counts stay complete
  mutate(status = factor(status, levels = c("missing", "dead", "prisoner", "released_prisoner")))

dts3 |>
  summarise(dx = sum(dx), .by = status)

write_rds(dts3, "data_inter/ukr_ualosses_conflict_deaths_sex_age_2022_2025.rds")

copy_this(
  dts3 |>
    summarise(dx = sum(dx), .by = c(year, status)) |>
    spread(status, dx) |>
    mutate(rt = missing / (missing + dead))
)

cols <- c("#457b9d", "#e63946")

