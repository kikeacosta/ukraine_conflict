# Migration inputs

Everything step 06 needs to build net emigration and its PERT bounds, saved
locally so the pipeline does not depend on web pages that can change. The full
reasoning is in `documents/migration_methodology.md`; every file in the project
is catalogued in `documents/data_sources.md`.

## Files read by the code

| File | What it is |
|---|---|
| `unhcr/unhcr_population_coo_UKR_2021_2025.csv` | UNHCR Refugee Data Finder: refugees, asylum seekers and other populations originating from Ukraine, by country of asylum and year. Flattened from the raw API response next to it |
| `unhcr/unhcr_population_coo_UKR_2021_2025_page1.json` | Raw API response, retrieved 2026-09-16 |
| `unhcr/unhcr_border_crossings_2022.csv` | UNHCR Operational Data Portal, Ukraine refugee situation: border crossings out of Ukraine and into it, cumulative since 24 February 2022 (population groups 5460 and 5472), at 1 March and each month end of 2022, read on 2026-09-19. Step 10 takes the timing of 2022's net outflow from it |
| `ces_2026_figures.csv` | Figures transcribed from the CES fifth-wave report, with printed page and quote. `used_by_code = yes` marks the ones step 06 reads |
| `national_programme_counts.csv` | Dated counts for Canada (CUAET) and the USA (Uniting for Ukraine), each with its source: programme arrivals in the first years, when few had yet left, and people present later. Step 06 interpolates between them to year ends |

The Eurostat temporary-protection stock that step 06 also reads lives in
`../refugees_eurostat/`.

## Canada and USA: what would improve the series

Neither country publishes a year-end count of Ukrainians present, so
`national_programme_counts.csv` holds the dated counts that are published, and
step 06 interpolates between them; step 13c moves the two series 30% up and down
(table A8). A better series would need, for **each of 31 Dec 2022, 2023, 2024 and
2025**:

- **Canada:** cumulative CUAET holders who *arrived* in Canada since 17 March
  2022 (not applications or approvals). Source: archived versions of IRCC's
  "Canada-Ukraine authorization for emergency travel: Key figures" page
  (Wayback Machine), IRCC open data, or an access-to-information request.
- **USA:** cumulative Ukrainians *paroled into* the USA under Uniting for
  Ukraine, ideally plus those paroled at ports of entry from 24 Feb to 25 Apr
  2022. Source: DHS Office of Homeland Security Statistics, USCIS, or a FOIA
  request.

Add them as rows to the CSV; no code changes are needed (see the note at the top
of `code/06_net_migration.R`).

## Reference copies (not read by code)

| Folder / file | What it is |
|---|---|
| `references/ces_2026_ukrainian_refugees_fifth_wave.pdf` | CES, *Ukrainian Refugees. Fifth Wave of Research*, 26 Feb 2026 |
| `references/pozniak_2023_ceemr_forced_migrants_ukraine.pdf` | Pozniak (2023), CEEMR 12(1), CC BY 4.0 |
| `references/ueffing_et_al_2023_jrc_ukraine_population_future.pdf` | Ueffing et al. (2023), JRC, CC BY 4.0 |
| `references/sbgsu_border_crossings_2022_2026.docx` | SBGSU border crossings table with coverage notes (same data as `../Pozniak_migration.csv`) |
| `web_snapshots/uscis_foia_2023-02-06_ukrainian_parolees_kassa.pdf` | USCIS letter: 93,928 U4U arrivals by 13 Dec 2022 |
| `web_snapshots/uscis_foia_2023-03-10_ukrainian_parolees_pocan.pdf` | USCIS letter on pre-U4U port-of-entry parolees |
| `web_snapshots/cicnews_2023-02_cuaet_arrivals_20260916.html` | CIC News, 24 Feb 2023: just over 167,000 CUAET arrivals |
| `web_snapshots/ces_fifth_wave_page_20260916.html` | CES web page for the fifth wave |
| `web_snapshots/ircc_qp_note_2024-00010_20260916.html` | IRCC Question Period note, 4 Apr 2024 |

### Still to save by hand

canada.ca refuses automated downloads, so these two sources behind the Canadian
counts are not yet in `web_snapshots/`. Open each in a browser
and save it as PDF:

- <https://www.canada.ca/en/immigration-refugees-citizenship/services/immigrate-canada/ukraine-measures/key-figures.html> (298,128 arrivals to 1 Apr 2024)
- <https://www.canada.ca/en/immigration-refugees-citizenship/corporate/transparency/committees/cimm-dec-05-2023/ukraine.html> (over 202,900 arrivals)

### Licensing

Pozniak (2023) and Ueffing et al. (2023) are CC BY 4.0. The CES report, news
articles and government pages are kept as research copies; check their terms
before redistributing the repository publicly.
