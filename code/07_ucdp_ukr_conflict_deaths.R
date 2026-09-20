# ==============================================================================
# STEP 07 (UCDP) - Conflict death TOTALS with uncertainty bounds
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Turns UCDP georeferenced event records into yearly totals of civilian and
# combatant deaths, carrying UCDP's own low/best/high bounds through.
#
#   1. Events are split into inter-state conflicts (both sides are
#      governments) and everything else (intra-state).
#   2. UCDP attributes some deaths to an unknown side. These are redistributed
#      over side A / side B / civilians in proportion to the deaths that ARE
#      attributed, within each dyad and year.
#   3. Combatant deaths in inter-state conflicts are reassigned to the
#      combatant's OWN country, so the Ukrainian total excludes Russian losses.
#   4. The low and high bounds are scaled by the same shares as the point
#      estimate.
#
# These totals become the civilian min/mode/max in the parameter table (10).
#
# NOTE: an earlier version of this step carried no uncertainty bounds and was
# retired. This file was called ..._invals.R only to tell the two apart; with
# the old one gone the suffix has been dropped.
#
# The output keeps the name ukr_ucdp_invals.rds, where "invals" means
# intervals: unlike the other sources, UCDP supplies low and high bounds.
#
# INPUTS   data_input/ucdp/GEDEvent_v26_1.csv   (GED 26.1, 1989-2025, 274 MB;
#          data extracted by UCDP on 30 March 2026, so 2025 is final rather
#          than candidate events). Cached as data_inter/ucdp_ged_events_slim.rds,
#          so the raw download is only needed once. See data_input/README.md.
# OUTPUT   data_inter/ukr_ucdp_invals.rds   <- used by 10
# ==============================================================================

rm(list = ls())
source("code/00_setup.R")

# UCDP Georeferenced Event Dataset (GED)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The raw GED file is 274 MB, over GitHub's hard limit, and is not tracked in
# git. Only the twelve columns below are ever used, so the cached extract is a
# few MB and is what the repository actually ships. See cache_rds() in
# 00_setup.R.
all2 <- cache_rds("data_inter/ucdp_ged_events_slim.rds", {
  keep <- c(
    "year", "type_of_violence", "conflict_name", "side_a", "side_b",
    "country", "adm_1", "deaths_a", "deaths_b", "deaths_civilians",
    "deaths_unknown", "best", "low", "high"
  )
  read_csv(
    require_raw("data_input/ucdp/GEDEvent_v26_1.csv"),
    show_col_types = FALSE
  ) %>%
    select(all_of(keep)) %>%
    select(
      year,
      country,
      adm1 = adm_1,
      conflict_name,
      type_of_violence,
      side_a,
      side_b,
      a = deaths_a,
      b = deaths_b,
      c = deaths_civilians,
      u = deaths_unknown,
      t = best,
      t_l = low,
      t_u = high
    ) %>%
    replace_na(list(a = 0, b = 0, c = 0, t = 0)) %>%
    filter(t_l + t > 0)
})

# Crimea and Sevastopol are outside the population this study measures, so the
# civilians killed there are outside its numerator. Civilians are counted by the
# place of the event, so the events themselves are dropped; combatants are
# counted by the side they fought for, wherever they died, and are unaffected -
# the few Ukrainian personnel killed there stay in the Ukrainian total, as they
# are in the register (Methods 1, Supplementary S1.3).
crimea_adm1 <- regex("crimea|sevastopol", ignore_case = TRUE)
in_crimea <- all2$country == "Ukraine" & str_detect(coalesce(all2$adm1, ""), crimea_adm1)
crimea_civ <- sum(all2$c[in_crimea])
cat(sprintf("\nCrimea and Sevastopol: %s civilian deaths in %d events left out of the count\n",
            format(crimea_civ, big.mark = ","), sum(in_crimea)))
# the civilian deaths come off the event's totals and bounds before anything is
# aggregated, so the redistribution of the unknown-side deaths sees the same book
all2$t[in_crimea]   <- all2$t[in_crimea]   - all2$c[in_crimea]
all2$t_l[in_crimea] <- pmax(0, all2$t_l[in_crimea] - all2$c[in_crimea])
all2$t_u[in_crimea] <- pmax(0, all2$t_u[in_crimea] - all2$c[in_crimea])
all2$c[in_crimea]   <- 0
stopifnot(all(all2$c >= 0), all(all2$t >= 0))

all_sum <-
  all2 |>
  # filter(country == "Ukraine", year >= 2022) |>
  summarise(
    civilians = sum(c),
    combatants = sum(a + b),
    unk = sum(u),
    tot2 = sum(t),
    tot_l = sum(t_l),
    tot_u = sum(t_u),
    .by = c(year, country)
  ) |>
  # Spread the deaths of unknown side over civilians/combatants in proportion
  # to the deaths that ARE attributed, and carry the same shares onto the low
  # and high bounds.
  #
  # NOTE: mutate() evaluates sequentially, so the shares have to be computed
  # into their own columns first. Deriving civilians_l from an already
  # rescaled `civilians` (as this block used to do) applies tot/tot_sum twice
  # and inflates both bounds.
  mutate(
    tot = civilians + combatants + unk,
    tot_sum = civilians + combatants,
    sh_civ = if_else(tot_sum > 0, civilians / tot_sum, 0),
    sh_cmb = if_else(tot_sum > 0, combatants / tot_sum, 0),
    civilians = sh_civ * tot,
    civilians_l = sh_civ * tot_l,
    civilians_u = sh_civ * tot_u,
    combatants = sh_cmb * tot,
    combatants_l = sh_cmb * tot_l,
    combatants_u = sh_cmb * tot_u
  ) |>
  select(-unk, -starts_with("tot"), -sh_civ, -sh_cmb)

unique(all2$country) %>% sort

all_sum2 <-
  all_sum |>
  select(country, year, civilians, civilians_l, civilians_u) |>
  rename(dts = civilians, dts_l = civilians_l, dts_u = civilians_u) |>
  mutate(role = "civilians") |>
  bind_rows(
    all_sum |>
      select(country, year, combatants, combatants_l, combatants_u) |>
      rename(dts = combatants, dts_l = combatants_l, dts_u = combatants_u) |>
      mutate(role = "combatants")
  )


# write_rds(all_sum, "data_inter/dt_ucdp_unadjusted.rds")

# identifying civilians and foreign and local soldiers in inter-state conflicts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inter_state <-
  all2 %>%
  # filtering only inter-state conflicts where the two sides are governments
  filter(type_of_violence == 1 & str_detect(side_b, "Government"))

inter_state2 <-
  inter_state %>%
  summarise(
    a = sum(a),
    b = sum(b),
    c = sum(c),
    u = sum(u),
    t = sum(t),
    t_l = sum(t_l),
    t_u = sum(t_u),
    .by = c(country, side_a, side_b, year)
  ) %>%
  group_by(country, side_a, side_b) %>%
  mutate(
    a2 = case_when(
      u == 0 ~ a,
      u != 0 & a != 0 & b != 0 & c != 0 ~ a + u * (a / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ a +
        u * sum(a) / (sum(a) + sum(b) + sum(c))
    ),
    b2 = case_when(
      u == 0 ~ b,
      u != 0 & a != 0 & b != 0 & c != 0 ~ b + u * (b / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ b +
        u * sum(b) / (sum(a) + sum(b) + sum(c))
    ),
    c2 = case_when(
      u == 0 ~ c,
      u != 0 & a != 0 & b != 0 & c != 0 ~ c + u * (c / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ c +
        u * sum(c) / (sum(a) + sum(b) + sum(c))
    ),
    t2 = a2 + b2 + c2,
    diff = t - t2
  ) %>%
  ungroup() %>%
  select(country, side_a, side_b, year, a = a2, b = b2, c = c2, t, t_l, t_u)

# identifying civilians and foreign and local soldiers in intra-state conflicts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
intra_state <-
  all2 %>%
  anti_join(inter_state) %>%
  group_by(country, year) %>%
  summarise(
    a = sum(a),
    b = sum(b),
    c = sum(c),
    u = sum(u),
    t = sum(t),
    t_l = sum(t_l),
    t_u = sum(t_u)
  ) %>%
  ungroup() %>%
  group_by(country) %>%
  mutate(
    a2 = case_when(
      u == 0 ~ a,
      u != 0 & a != 0 & b != 0 & c != 0 ~ a + u * (a / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ a +
        u * sum(a) / (sum(a) + sum(b) + sum(c))
    ),
    b2 = case_when(
      u == 0 ~ b,
      u != 0 & a != 0 & b != 0 & c != 0 ~ b + u * (b / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ b +
        u * sum(b) / (sum(a) + sum(b) + sum(c))
    ),
    c2 = case_when(
      u == 0 ~ c,
      u != 0 & a != 0 & b != 0 & c != 0 ~ c + u * (c / (a + b + c)),
      u != 0 & (a == 0 | b == 0 | c == 0) ~ c +
        u * sum(c) / (sum(a) + sum(b) + sum(c))
    ),
    t2 = a2 + b2 + c2,
    diff = t - t2
  ) |>
  ungroup()

intra_state2 <-
  intra_state %>%
  select(country, year, a = a2, b = b2, c = c2, t, t_l, t_u)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identifying combatants nationalities and civilians
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combatants killed in inter_state conflict
comb_inter <-
  inter_state2 %>%
  mutate(
    country_a = str_replace(side_a, "Government of ", ""),
    country_b = str_replace(side_b, "Government of ", "")
  ) %>%
  summarise(
    dts_a = sum(a),
    dts_a_l = sum(a * t_l / t),
    dts_a_u = sum(a * t_u / t),
    dts_b = sum(b),
    dts_b_l = sum(b * t_l / t),
    dts_b_u = sum(b * t_u / t),
    .by = c(country_a, country_b, year)
  )

# reasigning killed combatants to their own country
comb_inter2 <-
  bind_rows(
    comb_inter %>%
      select(
        country = country_a,
        year,
        dts = dts_a,
        dts_l = dts_a_l,
        dts_u = dts_a_u
      ),
    comb_inter %>%
      select(
        country = country_b,
        year,
        dts = dts_b,
        dts_l = dts_b_l,
        dts_u = dts_b_u
      )
  ) %>%
  mutate(role = "combatants")

# combatants killed in intra_state conflict
comb_intra <-
  intra_state2 %>%
  summarise(
    dts = sum(a + b),
    dts_l = sum((a + b) * t_l / t),
    dts_u = sum((a + b) * t_u / t),
    .by = c(country, year)
  ) %>%
  mutate(role = "combatants")

civils <-
  bind_rows(inter_state2, intra_state2) %>%
  summarise(
    dts = sum(c),
    dts_l = sum(c * t_l / t),
    dts_u = sum(c * t_u / t),
    .by = c(country, year)
  ) %>%
  mutate(role = "civilians")

# testing consistency before and after imputation ~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dts <-
  bind_rows(comb_inter2, comb_intra, civils) %>%
  summarise(
    dts = sum(dts),
    dts_l = sum(dts_l),
    dts_u = sum(dts_u),
    .by = c(country, year, role)
  ) %>%
  mutate(
    country = ifelse(country == "Yemen (North Yemen)", "Yemen", country),
    code = countrycode(country, origin = "country.name", destination = "iso3c")
  ) %>%
  drop_na(dts)

dts %>%
  summarise(dts = sum(dts))

all2 %>%
  summarise(dts = sum(t))

dts %>%
  summarise(dts_u = sum(dts_u))

all2 %>%
  summarise(dts_u = sum(t_u))

# (the full event-level table used to be written to data_inter/dt_ucdp_invals.rds
# here and then overwritten further down by a completely different object. Both
# writes were unread by the rest of the pipeline; the slim cached extract at
# data_inter/ucdp_ged_events_slim.rds now serves that purpose.)


# in Ukraine ====
# ~~~~~~~~~~~~~~~
ukr <-
  dts %>%
  filter(country == "Ukraine", year >= 2022)

copy_this(
  ukr %>%
    summarise(
      dts = sum(dts),
      dts_l = sum(dts_l),
      dts_u = sum(dts_u),
      .by = c(role, year)
    )
)

copy_this(
  all2 |>
    filter(country == "Ukraine", year >= 2022) |>
    summarise(
      civilians = sum(c),
      combatants = sum(a + b),
      unk = sum(u),
      total = sum(t)
    )
)

write_rds(ukr, "data_inter/ukr_ucdp_invals.rds")

# How UCDP attributes this war's deaths: by the side the dead fought for, and by
# the place the event happened. The two do not coincide, which is why residents
# of the occupied east who died in Russian-controlled forces cannot be recovered
# from it (Methods 1, Supplementary S1.3), and why Ukrainian losses inside Russia
# are in the Ukrainian total although the events are not in Ukraine.
dyad <- all2 |> filter(year >= 2022, str_detect(conflict_name, "Russia - Ukraine"))
crimea <- regex("crimea|sevastopol", ignore_case = TRUE)
attribution <- list(
  # deaths of Ukraine's side (side B of the state dyad), by where the event was
  ukr_side_by_country = dyad |> summarise(deaths = sum(b), .by = c(country, year)),
  # civilians, by where the event was: the estimand counts those in Ukraine
  civilians_by_country = dyad |> summarise(deaths = sum(c), .by = c(country, year)),
  # the fighting in the occupied east is located in Ukraine, under Kyiv's names
  by_adm1 = dyad |>
    filter(country == "Ukraine") |>
    summarise(events = n(), deaths = sum(t), .by = adm1) |>
    arrange(desc(deaths)),
  # Ukraine's own losses in events outside Ukraine, by region: the Kursk
  # incursion of August 2024 onwards is almost all of it
  ukr_side_outside = dyad |>
    filter(country != "Ukraine") |>
    summarise(events = n(), deaths = sum(b), .by = c(country, adm1)) |>
    arrange(desc(deaths)),
  # Crimea and Sevastopol are outside the estimand's territory
  crimea_civilians = dyad |>
    filter(country == "Ukraine", str_detect(coalesce(adm1, ""), crimea)) |>
    summarise(events = n(), civilians = sum(c), deaths = sum(t))
)
write_rds(attribution, "data_inter/ukr_ucdp_attribution.rds")
cat("\n=== HOW UCDP ATTRIBUTES THE DEAD: SIDE AGAINST PLACE ===\n")
print(as.data.frame(attribution$ukr_side_by_country |>
                      summarise(`Ukraine's side` = sum(deaths), .by = country)))
print(as.data.frame(attribution$crimea_civilians))

# Unadjusted vs adjusted totals for Ukraine and Russia, used by the source
# tables assembled in 15.
#
#   unadjusted - every death recorded in events located in that country,
#                in UCDP's own categories. Combatant deaths here include
#                soldiers of BOTH sides, and the deaths of unknown side are
#                still held separately.
#   adjusted   - after combatant deaths are reassigned to the combatant's own
#                nationality and the unknowns are redistributed.
#
# UCDP publishes low/best/high per EVENT, not per category, so the bounds for
# a category are obtained by scaling it by that event's low/best and high/best
# ratio - the same convention used throughout this script.
write_rds(
  bind_rows(
    all2 |>
      filter(country %in% c("Ukraine", "Russia (Soviet Union)"), year >= 2022) |>
      # the ratio is taken event by event here, and a handful of events carry
      # best = 0, so the division needs a guard or the whole sum goes NaN
      summarise(
        civilians = sum(c),
        civilians_l = sum(if_else(t > 0, c * t_l / t, 0)),
        civilians_u = sum(if_else(t > 0, c * t_u / t, 0)),
        combatants = sum(a + b),
        combatants_l = sum(if_else(t > 0, (a + b) * t_l / t, 0)),
        combatants_u = sum(if_else(t > 0, (a + b) * t_u / t, 0)),
        unknown = sum(u),
        .by = country
      ) |>
      mutate(source = "UCDP", type = "unadjusted"),
    dts |>
      filter(country %in% c("Ukraine", "Russia (Soviet Union)"), year >= 2022) |>
      summarise(
        across(c(dts, dts_l, dts_u), sum),
        .by = c(country, role)
      ) |>
      pivot_wider(
        names_from = role,
        values_from = c(dts, dts_l, dts_u),
        names_glue = "{role}_{.value}"
      ) |>
      rename(
        civilians = civilians_dts, civilians_l = civilians_dts_l,
        civilians_u = civilians_dts_u,
        combatants = combatants_dts, combatants_l = combatants_dts_l,
        combatants_u = combatants_dts_u
      ) |>
      mutate(unknown = NA_real_, source = "UCDP", type = "adjusted")
  ) |>
    mutate(country = if_else(str_starts(country, "Russia"), "Russia", country)),
  "data_inter/ukr_rus_ucdp_source_comparison.rds"
)

unique(dts$country)

# in Russia ====
# ~~~~~~~~~~~~~~
dts %>%
  filter(country == "Russia (Soviet Union)", year >= 2022)

dts %>%
  filter(country == "Russia (Soviet Union)", year >= 2022) |>
  summarise(
    dts = sum(dts),
    dts_l = sum(dts_l),
    dts_u = sum(dts_u),
    .by = c(role)
  )

all2 |>
  filter(country == "Russia (Soviet Union)", year >= 2022) |>
  summarise(
    civilians = sum(c),
    combatants = sum(a + b),
    unk = sum(u),
    total = sum(t)
  )

# adjusted vs unadjusted, long form, for inspection
tst <-
  dts |>
  mutate(type = "adjusted") |>
  bind_rows(all_sum2 |> mutate(type = "unadjusted"))

# adjusted totals by country/year/role, used by the diagnostic below
cmp_adj <-
  dts |>
  summarise(dts_adj = sum(dts), .by = c(country, year, role)) |>
  bind_rows(
    dts |>
      summarise(dts_adj = sum(dts), .by = c(country, year)) |>
      mutate(role = "total")
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tst2 <-
  all2 |>
  # filter(year == 2023) |>
  summarise(
    civilians = sum(c),
    combatants = sum(a + b),
    unk = sum(u),
    total = sum(t),
    .by = c(country, year)
  ) |>
  gather(-country, -year, key = role, value = dts) |>
  left_join(cmp_adj, by = c("country", "year", "role"))


tst2 |>
  filter(role == "total") |>
  mutate(diff = round(dts - dts_adj)) |>
  arrange(-diff)

copy_this(
  tst2 |>
    filter(role == "total") |>
    mutate(diff = round(dts - dts_adj)) |>
    arrange(-diff)
)


# Colombia ====
# ~~~~~~~~~~~~~
all2 |>
  filter(country == "Colombia") |>
  summarise(
    total = sum(t),
    .by = c(year)
  ) |>
  print(n = 40)

all2 |>
  filter(country == "Colombia") |>
  summarise(
    total = sum(t)
  )

all2 |>
  filter(country == "Colombia") |>
  summarise(
    total = sum(t),
    .by = c(year)
  ) |>
  ggplot() +
  geom_line(aes(year, total))
