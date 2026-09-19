# Conflict and all-cause mortality in Ukraine, 2022–2025

Estimates the direct mortality impact of the war in Ukraine: how many deaths
were caused by conflict violence, how many years of life expectancy they cost,
and how that loss divides between civilian and combatant deaths — all with
uncertainty intervals from a Monte Carlo simulation over the disputed inputs.

Conflict deaths are added to a modelled counterfactual; no all-cause mortality
is observed for 2022–2025, so indirect effects of the war are not estimated.

## Running it

Open `ukraine_conflict.Rproj` (paths are relative to the project root) and run
the whole pipeline, which executes every step below in order and ends with the
manuscript figures and tables:

```r
source("run_pipeline.R")
```

`00_setup.R` is sourced by each step and installs anything missing. Steps can
also be run one at a time in numerical order.

**You do not need the large raw data files.** Every step that reads one caches
a small extract into `data_inter/`, and those extracts are in the repository.
Steps re-run the heavy processing only if their cache is missing. See
[`data_input/README.md`](data_input/README.md).

The one exception is step 11, which regenerates the simulation draws because
they are too large to track. Their number is set once, as `n_sim` in
`code/00_setup.R`: 1,000 while the pipeline is under construction (about a
minute), 20,000 for reported estimates (about half an hour). The cache is named
for the size, so the two never mix, and `tables/_run_provenance.csv` records the
size every table was built at.

## The pipeline

| Step | Does | Key output |
|---|---|---|
| `00_setup.R` | Packages, caching helper, life table and Arriaga functions, shared assumptions | — |
| `01` | SSSU population and deaths by age, sex, region | `ukr_pop_sssu.rds` |
| `02` | Kannisto extrapolation of old-age mortality | `ukr_mx_1989_2021_adj.rds` |
| `03` | Period life tables 1989–2021 | `ukr_life_tables_1989_2021.rds` |
| `04` | Lee-Carter forecast → counterfactual "no war" mortality | `ukr_mxs_obs_plus_frcst_1989_2025.rds` |
| `05` | WPP2024 age-specific fertility | `ukr_asfr_wpp_2022_2025.rds` |
| `06` | Net migration by age and sex | `ukr_migrants_...rds` |
| `07_ucdp_ukr_conflict_deaths` | Conflict death **totals** with low/high bounds | `ukr_ucdp_invals.rds` |
| `07_ohchr` | Age-sex **profile** of civilian deaths | `ukr_ohchr_civilian_casualties.rds` |
| `07_acled` | ACLED totals, comparison only | `ukr_acled.rds` |
| `08` | Combatant deaths by age and sex from the ualosses register | `ukr_ualosses_..._sex_age_...rds` |
| `09` | Linkage between register releases, Markov imputation of how the "missing" resolve, and the range of the share alive | `ukr_ualosses_..._imputed_...rds`, `ukr_alpha_missing.rds` |
| `10` | min / mode / max parameter table | `ukr_param_table.rds` |
| `11` | Monte Carlo cohort-component projection | `ukr_sim_draws_2022_2025_n<n_sim>.rds` |
| `12` | Figures for the mortality estimates | `figures/mort_rates_*.png` |
| `13` | Sensitivity: the migration denominator effect | `ukr_migration_decomposition.rds` |
| `13b` | Sensitivity: the share of the missing who are alive, and sampling error in the resolution rates | `ukr_alpha_sensitivity_*.rds` |
| `13c` | Sensitivity: the range of each net migration component | `ukr_migration_sensitivity_e0.rds` |
| `14` | **Life expectancy loss decomposed by cause** | `ukr_e0_loss_by_cause_*.rds` |
| `15` | **All manuscript figures and tables** | `figures/fig*.png`, `tables/table*.csv` |
| `a01` | UNICEF export (not part of the paper): every iteration, by cause | `data_inter/unicef/<yymmdd>_ukraine_mx_estimates_*.csv` |

Steps 01→04 build the counterfactual. 05→09 assemble the conflict deaths and
the demographic components. 10→11 run the simulation. 12→15 report. `a01` is a separate deliverable for UNICEF, not part of the paper.

All estimates are produced and stored in **single years of age**. The only
place they are abridged is the UNICEF export in step a01, which aggregates to
the 18 groups (0, 1, 5, 10, … 80+) that file has always used.

## Method in brief

Population is projected from 1 January 2022 one year at a time, with mid-year
timing for net migration, conflict deaths and expected deaths. The counterfactual
is a Lee-Carter forecast fitted to **2000–2019** — 2020 and 2021 are excluded
so COVID excess mortality does not enter the "no war" baseline.

Conflict death totals and migration are uncertain, so they are drawn from PERT
distributions and every draw is projected through the full accounting. Civilian
deaths are drawn independently by year. Military deaths follow from one draw of
the share of the never-resolved missing who are alive, and each migration
component takes one quantile per draw, both held across all four years, because
their uncertainty is one systematic quantity rather than four separate ones. Step
14 then builds two life tables per
draw (baseline and observed), runs an Arriaga decomposition, and splits each
age's contribution between causes in proportion to that age's conflict deaths.

That split is exact rather than approximate: both life tables share the same
exposure denominator, so `mx_all − mx_expected` *is* the conflict death rate
and conflict deaths are additive across causes. Step 14 asserts this on every
run, along with the requirement that the decomposition reproduces the life
expectancy gap.

## Assumptions worth knowing

| Where | Assumption |
|---|---|
| `00_setup`, `09` | The share of the never-resolved missing who are alive cannot be estimated from the register, which records deaths but not returns. It is drawn between 0.009 and 0.203, centred on about 7,000 prisoners of war held by Russia over the 81,759 never-resolved (0.086), on the assumption that a long-term missing soldier who is alive is a prisoner. It carries about 90% of the imputed deaths; step 13b shows the results across [0, 1] (table A6, figure A4). |
| `09` | A person recorded as missing who cannot be found in the later register release is treated as still missing, not as resurfaced alive. |
| `11` | Imputed combatant deaths take the age–sex profile of the missing; registered deaths that of the confirmed dead. |
| `00_setup` | Sex ratio at birth `srb = 1.06`. |
| `11`, `13` | WPP2024 publishes fertility only to 2023; the 2023 schedule is carried forward to 2024–2025. |
| `11` | Starts from continental Ukraine (SSSU region `cnt`): the whole country without Crimea and Sevastopol, **including** Donetsk and Luhansk. SSSU population estimates for those two regions assume complete registration, which has not held since 2015. |
| `13` | Conflict deaths are absolute counts, so net migration shrinks the denominator and raises the rates. Step 13 quantifies how much of the loss this accounts for; step 13c how far the loss moves across the range of each migration component. |

## Layout

```
code/              pipeline, numbered in execution order
  deprecated/      superseded scripts, kept for reference
data_input/        raw sources (large ones untracked — see its README)
data_inter/        intermediate results and cached extracts
data_deprecated/   retired inputs and outputs, untracked
figures/           manuscript figures, written only by step 15
  exploratory/     diagnostic plots from steps 02-14
  deprecated/      figures from previous versions, untracked
tables/            manuscript tables as .csv, written only by step 15
```

## Figures and tables

Everything the paper needs comes from **step 15**, so the whole set can be
regenerated in one run and the styling stays consistent. Steps 02–14 still
draw their own diagnostic plots, but those go to `figures/exploratory/` — so
anything sitting directly in `figures/` is a manuscript deliverable.

| Output | Content |
|---|---|
| `fig1_combatant_deaths_missing.png` | Registered deaths vs unresolved disappearances, counts and shares |
| `fig2_mortality_rates_by_age.png` | All-cause vs counterfactual death rates by age, with 95% bands |
| `fig3_life_expectancy_loss.png` | Distribution of the e0 loss across draws |
| `fig4_decomposition_migration_mortality.png` | Mortality vs the migration denominator effect |
| `fig5_decomposition_by_cause.png` | e0 loss by civilians / registered / imputed combatants |
| `fig6_uncertainty_shares.png` | How much of the spread in the loss each input accounts for |
| `figA1_pert_draw_distributions.png` | Beta-PERT draw densities against their input bounds |
| `figA2_pert_draw_distributions_migration.png` | The same for the two migration components |
| `figA3_cumulative_migration_draws.png` | Cumulative net migration across draws |
| `figA4_alpha_sensitivity_missing.png` | Military deaths and e0 loss across the share of the missing who are alive, with the evidence range and the 95% interval of the draws |
| `figA5_migration_sensitivity.png` | e0 loss across the range of each net migration component |
| `table1_source_totals.csv` | What each source reports, with bounds and age-sex availability |
| `table2_missing_imputation.csv` | Missing combatants imputed to dead or alive |
| `table3_pert_input_bounds.csv` | Year-specific PERT min / mode / max |
| `table4_conflict_deaths_by_cause.csv` | Deaths by year, cause and sex, with intervals |
| `table5_life_expectancy_loss.csv` | e0 expected, observed and lost, by year and sex |
| `table6_totals_by_cause.csv` | Deaths by cause, summed over years |
| `table7_totals_by_year.csv` | Deaths by year, summed over causes |
| `tableA1_source_reconciliation.csv` | UCDP vs ACLED, Ukraine and Russia |
| `tableA2_status_transitions.csv` | Observed transitions of the missing between registers |
| `tableA3_e0_loss_by_cause.csv` | Years of life expectancy lost, by cause |
| `tableA4_uncertainty_shares.csv` | How much of the spread in the loss each input accounts for |
| `tableA5_military_reconciliation.csv` | Ukrainian military deaths: UCDP against the register and this study |
| `tableA6_alpha_sensitivity.csv` | Military deaths and e0 loss across the share of the missing who are alive |
| `tableA7_migration_sensitivity.csv` | 2025 e0 loss across the range of each net migration component, and with no migration |
| `_run_provenance.csv` | Draws, seed, draws file and commit the tables were built from |
