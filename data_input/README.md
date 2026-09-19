# Raw inputs

**Every input file and data source, with provenance, references and URLs, is
catalogued in [`../documents/data_sources.md`](../documents/data_sources.md).**
Migration inputs, including saved copies of the web sources behind them, are in
[`migration/`](migration/README.md).

> ⚠️ `migration/national_programme_arrivals_PLACEHOLDER.csv` approximates
> Canada and USA arrivals until year-end series are added. Search the code for
> `TODO-PLACEHOLDER`.

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
| `ucdp_ged_events_slim.rds` (1.6 MB) | `ucdp/GEDEvent_*.csv` | 259 MB |
| `ukr_asfr_wpp_2022_2025.rds` (4 kB) | `WPP2024_FERT_F01_*.xlsx` | 78 MB |
| `ualosses_counts_2022_2025.rds` (5 kB) | `260514_UKR_ualosses_Personnel.xlsx` | 27 MB |
| `ualosses_status_transitions.rds` (2 kB) | both ualosses registers | 71 MB |
| `acled_fatalities_slim.rds` (14 kB) | `acled/*.xlsx` | 88 kB |

## Where to get the untracked files

Only needed if you want to rebuild a cache from scratch.

### UCDP Georeferenced Event Dataset → `ucdp/`

- `GEDEvent_v25_1.csv` — UCDP GED **global** release, v25.1 (1989–2024).
  <https://ucdp.uu.se/downloads/>
- `GEDEvent_v25_01_25_12.csv` — UCDP **Candidate Events** dataset,
  January–December 2025. This is a separate product from the annual release:
  <https://ucdp.uu.se/downloads/candidateged/GEDEvent_v25_01_25_12.csv>

### UN World Population Prospects 2024

- `WPP2024_FERT_F01_FERTILITY_RATES_BY_SINGLE_AGE_OF_MOTHER.xlsx`
  <https://population.un.org/wpp/downloads>

  Note: WPP2024 publishes observed fertility only through **2023**. Steps 11
  and 13 carry the 2023 schedule forward to 2024 and 2025.

### Ukrainian military personnel losses (ualosses)

Compiled by Olivier Hubert:
<https://www.kaggle.com/datasets/ol4ubert/confirmed-ukrainian-military-personnel-losses>

- `260514_UKR_ualosses_Personnel.xlsx` — register v18, May 2026
- `ualosses_hubert_datasets/250916_UKR_ualosses_Personnel_v14.xlsx` — v14,
  September 2025, used only to observe how "missing" cases resolve over time

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
| `ohchr_civilian_deaths.xlsx` | OHCHR annual civilian death totals; supplied by OHCHR |
| `UKR_2022-2025_OHCHR.xlsx` | OHCHR civilian deaths by age and sex, supplied by OHCHR |
| `refugees_eurostat/` | Eurostat temporary protection stock by age and sex (read by 06), and the TP decisions series (evaluated, not used) |
| `migration/` | UNHCR Data Finder extract, CES figures, Canada/USA placeholder anchors, and reference copies — see `migration/README.md` |
| `Pozniak_migration.csv` | SBGSU border crossings 2022–2026. Not read by code; kept for reference |
| `ukr_centre_for_economic_strategy_age_sex_migrants.csv` | CES fifth-wave age-sex shares. Not read by code; kept for reference |
| `acled/` | ACLED fatality counts, comparison only. Not tracked - the cached extract ships instead |
