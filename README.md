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
`code/00_setup.R`: 20,000 for the reported estimates (about half an hour), or
1,000 for a quick run while the code changes (about a minute). The cache is named
for the size, so the two never mix, and `tables/_run_provenance.csv` records the
size every table was built at.

## The pipeline

| Step | Does | Key output |
|---|---|---|
| `00_setup.R` | Packages, caching helper, life table and Arriaga functions, shared assumptions | — |
| `01` | SSSU population and deaths by age, sex, region | `ukr_pop_sssu.rds` |
| `02` | Kannisto extrapolation of old-age mortality | `ukr_mx_1989_2021_adj.rds` |
| `03` | Period life tables 1989–2021 | `ukr_life_tables_1989_2021.rds` |
| `04` | Lee-Carter forecast → counterfactual "no war" mortality, and the forecast's error for the simulation to draw | `ukr_mxs_obs_plus_frcst_1989_2025.rds`, `ukr_lc_forecast_error.rds` |
| `04b` | Not part of the estimates: the territory of the 1989-2021 deaths and exposures (one territory in every year: Ukraine without Crimea, Sevastopol, Donetsk and Luhansk), and a test of the Lee-Carter index for a level shift at 2014-15 | `ukr_lc_level_shift.rds` |
| `05` | WPP2024 age-specific fertility | `ukr_asfr_wpp_2022_2025.rds` |
| `06` | Net migration by age and sex: a register stock differenced along cohorts, lightly smoothed across age within each Eurostat band | `ukr_migrants_...rds` |
| `07_ucdp_ukr_conflict_deaths` | Conflict death **totals** with low/high bounds | `ukr_ucdp_invals.rds` |
| `07_ohchr` | Age-sex **profile** of civilian deaths | `ukr_ohchr_civilian_casualties.rds` |
| `07_acled` | ACLED totals, comparison only | `ukr_acled.rds` |
| `08` | Combatant deaths by age and sex from the ualosses register (release v19) | `ukr_ualosses_..._sex_age_...rds` |
| `08b` | Registration lag: how much each event month of v19 will still grow, from the six releases - the dead net of records dropped, the missing by the new persons listed, a record being new only if its key is absent from the whole earlier release, so that an undated death that gains a date is not counted again - with 2,000 resampled sets of the completion factors | `ukr_registration_completion.rds` |
| `09` | Linkage of everyone listed as missing across the six releases (v14 to v19), the composition of the register's first resolutions by event year and the bound on the missing alive it gives, the imputation of the missing completed for registration lag, the evidence on the prisoners of war among them, the resolution hazards by months since disappearance as a description of the register, and the comparison of alternative rules and readings of the evidence | `ukr_ualosses_..._imputed_...rds`, `ukr_alive_missing.rds`, `ukr_military_inputs.rds`, `ukr_ualosses_projection_table.rds` |
| `09f`, `09g` | Investigation, not part of the estimates: per-person histories across the releases and a competing-risks multistate model, with the pass/fail test that keeps it out of production | `ukr_ualosses_multistate_test.rds` |
| `09h` | Not part of the estimates: a random sample of 600 links and drop-outs for clerical review by hand. Its output holds names and dates of birth, is written to the gitignored `documents/clerical_review/` and must never be committed | (local only) |
| `09i` | Superseded by the six-release linkage of 09, and kept for the record: it measured the duration model against v15 and v17 when those two releases were outside every fit | `ukr_ualosses_out_of_sample.rds`, `ualosses_window_transitions_six.rds` |
| `09j` | Not part of the estimates: fits the duration model on four of the five release windows and predicts the fifth with the average of the fitted windows' multipliers, each window held out in turn, and reports the error by outcome and by months since the event. Reads 09's cached linkage, so it needs no register files | `ukr_ualosses_held_out_window.rds` |
| `09k` | Not part of the estimates: whether an event-year cohort adds anything to duration in the register's resolutions, by quasi-F tests on cohort multipliers, and what the 2022 cohort's exits from the register would imply if they were carried to the later cohorts | `ukr_ualosses_cohort_overlap.rds` |
| `09l` | Not part of the estimates: how the bound on the missing alive moves with the windows behind it - on each cumulative span of releases, with each window left out and with the windows resampled - and the drop-outs and deaths by event year and window per 1,000 person-months of the missing at risk, the evidence for giving the events of 2022 the bound of the later cohorts. Reads tracked data only | `ukr_bound_stability.rds` |
| `10` | min / mode / max parameter table, and the timing of 2022's deaths and net outflow within the year | `ukr_param_table.rds`, `ukr_timing_absent.rds` |
| `11` | Monte Carlo cohort-component projection, drawing every uncertain input that has a distribution | `ukr_sim_draws_2022_2025_n<n_sim>.rds` |
| `12` | Figures for the mortality estimates | `figures/mort_rates_*.png` |
| `13` | Sensitivity: the migration denominator effect | `ukr_migration_decomposition.rds` |
| `13b` | Sensitivity: the share of the missing who are alive, beyond the range the evidence allows, and the linkage rules | `ukr_alive_sensitivity_*.rds`, `ukr_linkage_rules_e0.rds` |
| `13c` | Sensitivity: the range of each net migration component, and three specification checks of the migration input | `ukr_migration_sensitivity_e0.rds`, `ukr_migration_specification_e0.rds` |
| `13d` | Sensitivity: the horizon of the registration-lag correction, alone and crossed with the share alive | `ukr_registration_lag.rds` |
| `13e` | Sensitivity: the population base of Donetsk and Luhansk, 0.5 million lower to taken out | `ukr_denominator_donetsk_luhansk.rds` |
| `13f` | Sensitivity: the timing of 2022's deaths and departures against the mid-year convention | `ukr_timing_2022.rds` |
| `13g` | Sensitivity: the Lee-Carter window, the forecast one standard deviation either way, and the shape of the PERT distributions | `ukr_baseline_window_e0.rds`, `ukr_pert_shape_e0.rds` |
| `13h` | Sensitivity: an older age profile for the civilian deaths beyond OHCHR's verified count, mostly in occupied territory and besieged cities | `ukr_civilian_age_profile_e0.rds` |
| `13i` | Sensitivity: residents of occupied Donbas killed in Russian-controlled forces, counted to September 2023 or continued to 2025 | `ukr_donbas_fighters_e0.rds` |
| `13j` | Sensitivity: a coherent (Li-Lee) counterfactual, Ukraine forecast with eight EU neighbours from Eurostat, with its 95% band | `ukr_coherent_forecast_e0.rds` |
| `13l` | Sensitivity: the registration lag beyond four years as a fitted curve, exponential and inverse power, in place of a truncation or a constant rate, with its error from resampled event months | `ukr_registration_lag_tail.rds` |
| `13m` | Sensitivity: a smaller base at ages 20-64 outside Donetsk and Luhansk, with the Lee-Carter model refitted on exposures smaller by the same share, so that base and counterfactual stay coherent | `ukr_base_coherent_e0.rds` |
| `13n` | Sensitivity: a quarter and a half of 2021's excess mortality carried into the 2022 counterfactual | `ukr_covid_carryover_e0.rds` |
| `13o` | Sensitivity: the combatant total and the loss if the register never lists 10%, 25% or 40% of combatant deaths; a stated size for a bias nothing available measures, kept out of the second interval | `ukr_register_coverage_e0.rds` |
| `13k` | Not part of the estimates: fetches Belarus and Russia from the Human Mortality Database for `13j`'s pool, reading an account from `HMD_USER` and `HMD_PASS` in `~/.Renviron`. Neither series reaches 2019, so `13j` leaves out any country that does not cover its whole fit window; using them means ending that window earlier. The cache is gitignored: the HMD agreement governs redistribution | `data_input/hmd_neighbours/` (local) |
| `14` | **Life expectancy loss decomposed by cause** | `ukr_e0_loss_by_cause_*.rds` |
| `14b` | Years of life lost and adult mortality (45q15) | `ukr_yll_45q15_summary.rds` |
| `14c` | The loss with the counterfactual fixed at its point forecast, beside it drawn: the two tiers of uncertainty | `ukr_fixed_counterfactual_e0.rds` |
| `14d` | The loss with the population base of Donetsk and Luhansk drawn, PERT(0, 0.5, 1.0 million fewer), beside the base at the published estimates: the second tier of the structural uncertainty. The range is a judgement of the team, not a measurement, so it is reported beside the estimate and not inside it | `ukr_drawn_base_e0.rds` |
| `14e` | The intervals with the structural choices folded in: one alternative per choice sampled in every draw with equal weights and its shift from the central projection added, reported beside the interval of the drawn inputs. The weights are a convention, not a measurement | `ukr_structural_uncertainty.rds` |
| `15` | **All manuscript figures and tables** | `figures/fig*.png`, `tables/table*.csv` |
| `16` | Figure A0: the pipeline from sources to estimates, each step's box drawn dashed if its output is drawn in the simulation and solid if it is held; the release list and the number of draws come from the run | `figures/figA0_pipeline.png` |
| `a01` | UNICEF export (not part of the paper): every iteration, by cause | `data_inter/unicef/<yymmdd>_ukraine_mx_estimates_*.csv` |

Steps 01→04 build the counterfactual. 05→09 assemble the conflict deaths and
the demographic components. 10→11 run the simulation. 12→15 report. `a01` is a separate deliverable for UNICEF, not part of the paper.

All estimates are produced and stored in **single years of age**. The only
place they are abridged is the UNICEF export in step a01, which aggregates to
the 18 groups (0, 1, 5, 10, … 80+) that file has always used.

## Method in brief

Population is projected from 1 January 2022 one year at a time, with mid-year
timing for net migration, conflict deaths and expected deaths, except in 2022,
when the timing of deaths and net outflow is taken from the months of the
events. Each row of the
projection is one birth cohort: rates and death profiles given by completed age
are carried onto the cohorts as each lives the year, half at each of its two
ages, and the results are returned by completed age, so the life tables of
steps 12 and 14 reproduce the input rates. The counterfactual is a Lee-Carter
forecast fitted to **2000–2019** — 2020 and 2021 are excluded so COVID excess
mortality does not enter the "no war" baseline.

Every uncertain input that can be given a distribution is drawn, and every draw
is projected through the full accounting. Civilian deaths are drawn
independently by year from PERT distributions. Military deaths follow from one
draw of the evidence on how many of the missing are alive - the prisoners of war
among them, and where between none and its bound the share alive outside
captivity sits - with a draw of the weights on the release windows behind that
bound, and of the registration-lag factors.
The western
net outflow follows from one weight blending the two sources' four-year paths,
and Russia and Belarus from one draw. Each of these is held across the four
years, because its uncertainty is one systematic quantity rather than four
separate ones. The counterfactual follows each draw's own path of the
Lee-Carter index. The estimate is the mean of the draws, with their 2.5th and
97.5th percentiles as its interval: means add up across ages, sexes, years and
components, and agree with the central projection the sensitivity steps use.
Step
14 then builds two life tables per draw (baseline and observed), runs an
Arriaga decomposition, and splits each age's contribution between causes in
proportion to that age's conflict deaths: civilians, registered combatants,
their late registrations and the imputed missing.

That split is exact rather than approximate: both life tables share the same
exposure denominator, so `mx_all − mx_expected` *is* the conflict death rate
and conflict deaths are additive across causes. Step 14 asserts this on every
run, along with the requirement that the decomposition reproduces the life
expectancy gap.

## Assumptions worth knowing

| Where | Assumption |
|---|---|
| `00_setup`, `09` | A missing person is counted as dead unless there is evidence of life. The prisoners of war among the missing come from the official figures - about 7,000 held and 7,291 military personnel returned by February 2026, less the 11,325 the register records as prisoners - times the share of the returned whom the register had listed as missing (47.6% for 2024-2025 events, the mode). Of the rest, the share alive outside captivity lies between none and the share of an event year's first resolutions outside captivity that leave the register beyond list maintenance - beyond the rate at which the dead of the same cohort leave it in the same window - and is drawn uniformly between the two, one draw for the four event years. The events of 2022 take the bound of the events of 2023-2025: their own drop-outs are one clean-up of their records (09l), and their own bound (48%) is an alternative in table A20. The bound assumes the missing not yet resolved are no more often alive than those resolved. |
| `00_setup`, `09` | A person recorded as missing who is no longer listed — absent from a later release and from every release after it — is resolved alive. Before that, each release is searched under a corrected name or date of birth. A person whose first resolution is a return from captivity was a prisoner the register had not recorded, and is resolved alive, to captivity. |
| `08b`, `09` | The register's dead and missing are brought to the completeness events reach at four years; the late registrations are reported apart. The missing are completed by the new persons listed only, since those who drop off the register are resolutions. How the missing resolve by months since the event is described, not projected: the estimate counts them by the identity above. |
| `11` | Imputed combatant deaths take the age–sex profile of the missing; registered deaths that of the confirmed dead. |
| `07_ohchr`, `11` | Civilian deaths take the age–sex profile of the deaths OHCHR verified, mostly in government-controlled territory. Step 13h gives the deaths beyond OHCHR's verified count, mostly in occupied territory and besieged cities, an older profile. |
| `00_setup` | Sex ratio at birth `srb = 1.06`. |
| `11`, `13` | WPP2024 publishes fertility only to 2023; the 2023 schedule is carried forward to 2024–2025. |
| `10`, `11` | Deaths and net outflow are spread evenly within each year (the mid-year convention), except in 2022, when their timing is taken from the months of the events: OHCHR's civilians killed by month, the register's dead and missing by month, and UNHCR's net border crossings by month. |
| `11` | Starts from continental Ukraine (SSSU region `cnt`): the whole country without Crimea and Sevastopol, **including** Donetsk and Luhansk. SSSU population estimates for those two regions assume complete registration, which has not held since 2015. |
| `13` | Conflict deaths are absolute counts, so net migration shrinks the denominator and raises the rates. Step 13 quantifies how much of the loss this accounts for; step 13c how far the loss moves across the range of each migration component. |
| `06`, `11` | The western net outflow is bracketed by source: each draw blends the border-crossing reading's four-year path with the register reading's by one weight, so the four-year total stays between the two sources' own totals (4.29–5.27 million). |

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
documents/         methodology reports and the draft paper, kept locally, untracked
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
| `figA0_pipeline.png` | The pipeline from sources to estimates: what each source contributes, in which order, and which steps' outputs are drawn |
| `figA4_pert_draw_distributions.png` | Beta-PERT draw densities against their input bounds |
| `figA2_pert_draw_distributions_migration.png` | Draws of the two migration components: the western blend of the two readings, and Russia and Belarus |
| `figA3_cumulative_migration_draws.png` | Cumulative net migration across draws |
| `figA6_missing_alive_sensitivity.png` | Military deaths and e0 loss across the share of the missing who are alive, with the evidence range and the 95% interval of the draws |
| `figA7_migration_sensitivity.png` | e0 loss across the range of each net migration component |
| `figA1_duration_lexis.png` | The durations since disappearance each event month is observed at, window by window, and the projection to 48 months |
| `figA8_multistate_cif_by_cohort.png` | The multistate model of `09g`, evaluated and not used: cumulative incidence of each resolution of the missing by event-year cohort, drawn from its tracked curves |
| `table1_source_totals.csv` | What each source reports, with bounds and age-sex availability |
| `table2_missing_imputation.csv` | Missing combatants imputed to dead or alive |
| `table3_pert_input_bounds.csv` | Year-specific PERT min / mode / max |
| `tableA12_missing_alive_inputs.csv` | The inputs on the missing alive: the share of the unrecorded prisoners among the missing, the prisoners held, and the share of the unresolved alive for other reasons, with the means the deterministic analyses use |
| `table4_conflict_deaths_by_cause.csv` | Deaths by year, cause and sex, with intervals |
| `table5_life_expectancy_loss.csv` | e0 expected, observed and lost, by year and sex |
| `table6_totals_by_cause.csv` | Deaths by cause, summed over years |
| `table7_totals_by_year.csv` | Deaths by year, summed over causes |
| `table8_structural_sensitivity.csv` | How far the military total and the loss move with each structural choice the intervals do not cover |
| `tableA1_source_reconciliation.csv` | UCDP vs ACLED, Ukraine and Russia |
| `tableA8_status_transitions.csv` | Resolution of the missing between the first and the last release, composed from the six register releases |
| `tableA38_e0_loss_by_cause.csv` | Years of life expectancy lost, by cause |
| `tableA39_uncertainty_shares.csv` | How much of the spread in the loss each input accounts for |
| `tableA2_military_reconciliation.csv` | Ukrainian military deaths: UCDP against the register and this study |
| `tableA18_missing_alive_sensitivity.csv` | Military deaths and e0 loss across the share of the missing who are alive |
| `tableA26_migration_sensitivity.csv` | 2025 e0 loss across the range of each net migration component, and with no migration |
| `tableA27_migration_specification.csv` | e0 loss under alternative specifications of the migration input, and the Canada and USA series +/-30% |
| `tableA20_linkage_and_chain.csv` | Military deaths and e0 loss under alternative rules of the imputation and the linkage, models of resolution projected in its place, and readings of the prisoner-of-war evidence |
| `tableA22_registration_lag.csv` | The same across the horizon of the registration-lag correction |
| `tableA28_population_base_and_timing.csv` | e0 loss with a smaller Donetsk-Luhansk base, and with other timings of 2022's deaths and net outflow |
| `tableA32_counterfactual_window.csv` | e0 loss and counterfactual e0 by Lee-Carter window, with a coherent forecast with eight neighbours, and with the forecast one standard deviation either way; the 95% band of each forecast |
| `tableA34_pert_shape.csv` | e0 loss and its interval by PERT shape |
| `tableA40_years_of_life_lost.csv` | Years of life lost by year and cause |
| `tableA41_adult_mortality_45q15.csv` | 45q15 without and with conflict deaths |
| `tableA23_missing_alive_by_lag_horizon.csv` | Military deaths and the 2025 male loss over the missing alive and the registration-lag horizon together |
| `tableA43_military_triangulation.csv` | Military deaths against official statements and other estimates at the dates they were made |
| `tableA6_returned_prisoners_prior_status.csv` | Where the people the register records as returned from captivity were listed before |
| `tableA5_register_dropout.csv` | How often the missing and the dead leave the register, window by window |
| `tableA7_resolution_hazards.csv` | Monthly hazards of resolution by months since disappearance, and the release windows' multipliers |
| `tableA44_civil_register_check.csv` | The projection's deaths and births against those the Ministry of Justice registered |
| `tableA30_civilian_age_profile.csv` | e0 loss with an older age profile for the civilian deaths beyond OHCHR's verified count |
| `tableA31_donbas_fighters.csv` | e0 loss with the residents of occupied Donbas killed in Russian-controlled forces added |
| `tableA24_registration_lag_tail.csv` | Military deaths and e0 loss with the register's growth beyond four years modelled by a fitted curve |
| `tableA29_base_coherent.csv` | e0 loss with a smaller base at working ages outside Donetsk and Luhansk, the counterfactual as fitted and refitted on exposures smaller by the same share |
| `tableA33_covid_carryover.csv` | The 2022 loss with a quarter and a half of 2021's excess mortality carried into the counterfactual |
| `tableA46_held_out_window.csv` | The duration model against each release window held out of its fit, by outcome |
| `tableA15_structural_intervals.csv` | Military deaths, all conflict deaths and the e0 loss with the interval of the drawn inputs, and with the structural choices folded in |
| `tableA42_population_present.csv` | The population the rates are computed on, by draw: men aged 20-50, men, and both sexes, with 95% intervals (11, 15) |
| `tableA3_register_status_by_year.csv` | The register's records of 2022-2025 events by event year and status in v19, the cross-section behind the linkage (15) |
| `tableA25_register_coverage.csv` | Combatant deaths and the loss if the register misses 10%, 25% or 40% of combatant deaths (13o) |
| `figA5_structural_tornado.png` | The range of each structural choice beside the 95% interval of the drawn inputs |
| `_run_provenance.csv` | Draws, seed, draws file and commit the tables were built from |
