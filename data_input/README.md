# Raw inputs

**Every input file and data source, with provenance, references and URLs, is
catalogued in `documents/data_sources.md`, which is kept with the other project
documents outside the repository.**
Migration inputs, including saved copies of the web sources behind them, are in
[`migration/`](migration/README.md).

Most files in this folder are tracked in git. The large ones are **not** — see
`../.gitignore`.

## You probably do not need to download anything

Every script that reads a large raw file wraps that step in `cache_rds()`
(defined in `code/00_setup.R`). The small cached extract is written to
`data_inter/`, is tracked in git, and is reused on every later run. So a fresh
clone can run the entire pipeline **without** any of the large downloads.

The heavy step only re-runs if you delete the corresponding `.rds`, or
pass `refresh = TRUE`.

| Cached extract in `data_inter/` | Replaces | Raw size |
|---|---|---|
| `ucdp_ged_events_slim.rds` (1.7 MB) | `ucdp/GEDEvent_v26_1.csv` | 274 MB |
| `ukr_asfr_wpp_2022_2025.rds` (4 kB) | `WPP2024_FERT_F01_*.xlsx` | 78 MB |
| `ualosses_counts_2022_2025.rds` (3 kB) | the ualosses register, release v19 (step 08) | 30 MB |
| `ualosses_registration_by_month.rds` (3 kB) | the four ualosses releases, for registration lag (step 08b) | 110 MB |
| `ualosses_window_transitions.rds` (8 kB) | the four ualosses releases, for the linkage of the missing (step 09) | 110 MB |
| `acled_fatalities_slim.rds` (14 kB) | `acled/*.xlsx` | 88 kB |

## Where to get the untracked files

Only needed if you want to rebuild a cache from scratch.

### UCDP Georeferenced Event Dataset → `ucdp/`

- `GEDEvent_v26_1.csv` — UCDP GED **global** release, v26.1 (1989–2025; data
  extracted by UCDP on 30 March 2026), unzipped from
  <https://ucdp.uu.se/downloads/ged/ged261-csv.zip> (39 MB zipped). Codebook:
  <https://ucdp.uu.se/downloads/ged/ged261.pdf>.

### UN World Population Prospects 2024

- `WPP2024_FERT_F01_FERTILITY_RATES_BY_SINGLE_AGE_OF_MOTHER.xlsx`
  <https://population.un.org/wpp/downloads>

  Note: WPP2024 publishes observed fertility only through **2023**. The
  projection carries the 2023 schedule forward to 2024 and 2025.

### Ukrainian military personnel losses (ualosses)

Compiled by Olivier Hubert:
<https://www.kaggle.com/datasets/ol4ubert/confirmed-ukrainian-military-personnel-losses>

Four releases, in `ualosses_hubert_datasets/`, with their paths and dates in
`ual_releases` (`code/00_setup.R`):

- `250916_UKR_ualosses_Personnel_v14.xlsx` — v14, 16 September 2025
- `251204_UKR_ualosses_Personnel_v16.xlsx` — v16, 4 December 2025
- `260423_UKR_ualosses_Personnel_v18.xlsx` — v18, 23 April 2026
- `260919_UKR_ualosses_Personnel_v19.xlsx` — v19, 21 July 2026: every age–sex
  count, and the last release in the linkage

Each release is dated by its version on Kaggle, where the version history gives
the dates; each file's latest events fall a few days to a few weeks before them.
The v19 file is Kaggle's version 19 of 21 July 2026 (30.93 MB, latest events of
18 June 2026); its name carries the date it was saved, not the release date.
Versions 15 (2 November 2025) and 17 (6 February 2026) are on Kaggle too and are
not used.

`260514_UKR_ualosses_Personnel.xlsx` is the same list as v18, saved again;
nothing reads it.

⚠️ **These contain individual-level personal data (names, patronymics, dates
of birth) and must not be committed or redistributed.** The pipeline
aggregates them to anonymous counts before anything is written to
`data_inter/`.

## Tracked (small) inputs

| File | Source |
|---|---|
| `sssu_ukr_data.xlsx` | State Statistics Service of Ukraine — population and deaths. "Continental" means the whole country without Crimea; the Donetsk and Luhansk sheets are subsets of it, not complements |
| `LifeTables.xlsx` | Independent life tables. No longer read by any step — the SSSU-vs-life-table rate comparison in `01` was retired because nothing consumed its output. Kept for reference |
| `DataDxEx.csv` | Deaths and exposures 1989–2021, for the Kannisto fit and Lee-Carter model; supplied by SSSU |
| `a0.csv` | Infant mortality separation factors; supplied by SSSU |
| `ohchr_civilian_deaths.xlsx` | OHCHR annual civilian death totals, and the monthly count of civilians killed that times 2022's civilian deaths (step 10); supplied by OHCHR |
| `UKR_2022-2025_OHCHR.xlsx` | OHCHR civilian deaths by age and sex, supplied by OHCHR |
| `refugees_eurostat/` | Eurostat temporary protection stock by age and sex (read by 06), and the TP decisions series (evaluated, not used) |
| `eurostat_neighbours/` | Eurostat deaths by single year of age (`demo_magec`) and population on 1 January (`demo_pjan`), 2000–2020, for Poland, Slovakia, Hungary, Romania, Bulgaria, Lithuania, Latvia and Estonia: the reference group of the coherent counterfactual (13j). Fetched by 13j from Eurostat's dissemination API when absent, then read from here |
| `migration/` | UNHCR Data Finder extract, UNHCR's 2022 border crossings by month (the timing of 2022's net outflow), CES figures, the Canada and USA programme counts, and reference copies — see `migration/README.md` |
| `official_figures.csv` | Figures used only to check the estimates: official statements of military deaths, US and CSIS estimates, returns from captivity and repatriated bodies (Coordination Headquarters for the Treatment of Prisoners of War), prisoners of war held, registered deaths and births (Ministry of Justice), and deaths in the units of the self-proclaimed Donetsk and Luhansk republics; one row per figure, with its source and URL (step 15) |
| `Pozniak_migration.csv` | SBGSU border crossings 2022–2026. Not read by code; kept for reference |
| `ukr_centre_for_economic_strategy_age_sex_migrants.csv` | CES fifth-wave age-sex shares. Not read by code; kept for reference |
| `acled/` | ACLED fatality counts, comparison only. Not tracked - the cached extract ships instead |
