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
for the size, so the two never mix, and `outputs/manifest.csv` records the
size every table and figure of the paper was built at.

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
| `06` | Net migration by age and sex: UNHCR's year-end totals of refugees from Ukraine recorded outside Russia, given Eurostat's age-sex structure, differenced along cohorts and lightly smoothed across age within each Eurostat band, bracketed with net border crossings; Russia as one 2022 component | `ukr_migrants_...rds`, `ukr_migration_bounds.rds` |
| `07_ucdp_ukr_conflict_deaths` | Conflict death **totals** with low/high bounds | `ukr_ucdp_invals.rds` |
| `07_ohchr` | Age-sex **profile** of civilian deaths | `ukr_ohchr_civilian_casualties.rds` |
| `07_acled` | ACLED totals, comparison only | `ukr_acled.rds` |
| `08` | Combatant deaths by age and sex from the ualosses register (release v19) | `ukr_ualosses_..._sex_age_...rds` |
| `08b` | Registration lag: how much each event month of v19 will still grow, from the six releases - the dead net of records dropped, the missing by the new persons listed, a record being new only if its key is absent from the whole earlier release, so that an undated death that gains a date is not counted again, and a missing person only if no earlier record is theirs under a corrected or re-spelt key. Each release is dated by the data it holds, four days after its latest event. The growth is kept by release window, so that the completion factors can be rebuilt under the weights each simulation puts on the windows. The rates are held non-increasing by lag within each status: as measured, the missing's rise again after 24 months, on the 2022-23 cohorts alone, whose lists the register rebuilt from the Interior Ministry's database in 2025. The rates as measured, and the factors under two other measures of the growth (per person listed, and with the releases dated by Kaggle version), are kept for `13d` | `ukr_registration_completion.rds` |
| `09` | Linkage of everyone listed as missing across the six releases (v14 to v19), the composition of the register's first resolutions by event year and the bound on the missing alive it gives, the imputation of the missing completed for registration lag, the evidence on the prisoners of war among them, the resolution hazards by months since disappearance as a description of the register, and the comparison of alternative rules and readings of the evidence | `ukr_ualosses_..._imputed_...rds`, `ukr_alive_missing.rds`, `ukr_military_inputs.rds`, `ukr_ualosses_projection_table.rds` |
| `09f`, `09g` | Investigation, not part of the estimates: per-person histories across the releases and a competing-risks multistate model, with the pass/fail test that keeps it out of production | `ukr_ualosses_multistate_test.rds` |
| `09h` | Not part of the estimates: a random sample of 600 links and drop-outs for clerical review by hand. Its output holds names and dates of birth, is written to the gitignored `documents/team/clerical_review/` and must never be committed | (local only) |
| `09i` | Superseded by the six-release linkage of 09, and kept for the record: it measured the duration model against v15 and v17 when those two releases were outside every fit | `ukr_ualosses_out_of_sample.rds`, `ualosses_window_transitions_six.rds` |
| `09j` | Not part of the estimates: fits the duration model on four of the five release windows and predicts the fifth with the average of the fitted windows' multipliers, each window held out in turn, and reports the error by outcome and by months since the event. Reads 09's cached linkage, so it needs no register files | `ukr_ualosses_held_out_window.rds` |
| `09k` | Not part of the estimates: whether an event-year cohort adds anything to duration in the register's resolutions, by quasi-F tests on cohort multipliers, and what the 2022 cohort's exits from the register would imply if they were carried to the later cohorts | `ukr_ualosses_cohort_overlap.rds` |
| `09l` | Not part of the estimates: how the bound on the missing alive moves with the windows behind it - on each cumulative span of releases, with each window left out and with the windows resampled - and the drop-outs and deaths by event year and window per 1,000 person-months of the missing at risk, the evidence for pooling the bound of the events of 2023-2025 and for not reading the events of 2022 on their own. Reads tracked data only | `ukr_bound_stability.rds` |
| `10` | The parameter table of the simulation, and the timing of 2022's deaths and net outflow within the year. Combatants: the resolution-based scenario, with the lower and the upper bound as minimum and maximum. Civilians: UCDP's low, best and high values of each year, the minimum, mode and maximum of a Beta-PERT; and the draws of each year's deaths from uncertainUCDP (Vesco et al. 2026), event by event, for the alternative. Net migration: the two readings by component | `ukr_param_table.rds`, `ukr_civilians_uncertain_ucdp.rds`, `ukr_timing_absent.rds` |
| `11` | Monte Carlo cohort-component projection, drawing every uncertain input that has a distribution | `ukr_sim_draws_2022_2025_n<n_sim>.rds` |
| `12` | Figures for the mortality estimates | `figures/mort_rates_*.png` |
| `13` | Sensitivity: the migration denominator effect | `ukr_migration_decomposition.rds` |
| `13b` | Sensitivity: the share of the missing who are alive, beyond the range the evidence allows, and the linkage rules | `ukr_alive_sensitivity_*.rds`, `ukr_linkage_rules_e0.rds` |
| `13c` | Sensitivity: the range of each net migration component, and three specification checks of the migration input | `ukr_migration_sensitivity_e0.rds`, `ukr_migration_specification_e0.rds` |
| `13d` | Sensitivity: the horizon of the registration-lag correction, alone and crossed with the share alive, and the two other measures of the register's growth at the horizon used | `ukr_registration_lag.rds` |
| `13e` | Sensitivity: the population base of Donetsk and Luhansk, 0.5 million lower to taken out | `ukr_denominator_donetsk_luhansk.rds` |
| `13f` | Sensitivity: the timing of 2022's deaths and departures against the mid-year convention | `ukr_timing_2022.rds` |
| `13g` | Sensitivity: the Lee-Carter window, the forecast one standard deviation either way, and the shape of the PERT distributions | `ukr_baseline_window_e0.rds`, `ukr_pert_shape_e0.rds` |
| `13h` | Sensitivity: an older age profile for the civilian deaths beyond OHCHR's verified count, mostly in occupied territory and besieged cities | `ukr_civilian_age_profile_e0.rds` |
| `13i` | Sensitivity: residents of occupied Donbas killed in Russian-controlled forces, counted to September 2023 or continued to 2025 | `ukr_donbas_fighters_e0.rds` |
| `13j` | Sensitivity: a coherent (Li-Lee) counterfactual, Ukraine forecast with eight EU neighbours from Eurostat, with its 95% band | `ukr_coherent_forecast_e0.rds` |
| `13l` | Sensitivity: the registration lag beyond four years as a fitted curve, exponential and inverse power, in place of a truncation or a constant rate, with its error from the release windows re-weighted | `ukr_registration_lag_tail.rds` |
| `13m` | Sensitivity: a smaller base at ages 20-64 outside Donetsk and Luhansk, with the Lee-Carter model refitted on exposures smaller by the same share, so that base and counterfactual stay coherent | `ukr_base_coherent_e0.rds` |
| `13n` | Sensitivity: a quarter and a half of 2021's excess mortality carried into the 2022 counterfactual | `ukr_covid_carryover_e0.rds` |
| `13o` | Sensitivity: the combatant total and the loss if the register never lists 10%, 25% or 40% of combatant deaths; a stated size for a bias nothing available measures, kept out of the second interval | `ukr_register_coverage_e0.rds` |
| `13p` | Sensitivity: the 2022 civilian total and the loss by the rule that splits one UCDP event, the siege of Mariupol, whose 15,996 deaths of unknown status are 94% of the year's - none, the dyad-year share (used), Mariupol's own share or all of them civilian; a row of Table 6, not drawn and kept out of the second interval | `ukr_mariupol_split_e0.rds` |
| `13q` | Three figures reported beside the headline, each in the central projection: the loss without the imputed missing, the change per 10 points of the missing alive (from `13s`'s sweep), and a territory-consistent second estimate for Ukraine without Crimea, Donetsk and Luhansk, with the combatants split by the register's place of origin (counted in `08`) and the civilians by UCDP's location | `ukr_headline_companions.rds` |
| `13r` | Sensitivity: the register stock built bottom-up - Eurostat temporary protection with the clean-ups of invalid statuses back-cast, the Data Finder for other countries, Canada and the USA from programme counts - with Russia and Belarus as the separate component, and the loss under each stock; also writes the clean-ups for the docs | `ukr_bottom_up_stock.rds`, `ukr_bottom_up_stock_e0.rds`, `ukr_migration_breaks.rds` |
| `13s` | The estimate at every share of the missing counted dead, from none to all, every other input at its central value: the loss of male life expectancy by year and for 2022-2025 as one period, combatant deaths by year and the four-year total, with the lower bound, the resolution-based scenario and the upper bound marked (main-text figure) and, for the Supplementary Materials, the bodies scenario, a heavy-desertion reading, the crossing at which 2025 overtakes 2022 as the deadliest year for men, the rule scenario, and the register's resolution models projected in the scenario's place; bars from the draws with the share fixed. Runs before `13q`, which reads its sweep | `ukr_missing_dead_share_sweep.rds`, `ukr_missing_dead_share_models.rds`, `ukr_missing_dead_share_bars.rds`, `figures/exploratory/missing_dead_share_main.png`, `missing_dead_share_sm.png` |
| `13t` | Civilian deaths under the two methods: a Beta-PERT on UCDP's low, best and high values, which the estimate uses, and uncertainUCDP's distributions, with the loss under whichever is not used (`civilian_source` in `00_setup`); both beside OHCHR's verified counts, and the civilian deaths by the reported size of the event. Also the scenario with every death from UCDP: Ukraine's combatants drawn from uncertainUCDP as the civilians are, in place of the register, with the register's ages, and the same with the combatants at UCDP's best value | `ukr_uncertain_ucdp.rds` |
| `13k` | Not part of the estimates: fetches Belarus and Russia from the Human Mortality Database for `13j`'s pool, reading an account from `HMD_USER` and `HMD_PASS` in `~/.Renviron`. Neither series reaches 2019, so `13j` leaves out any country that does not cover its whole fit window; using them means ending that window earlier. The cache is gitignored: the HMD agreement governs redistribution | `data_input/hmd_neighbours/` (local) |
| `14` | **Life expectancy loss decomposed by cause** | `ukr_e0_loss_by_cause_*.rds` |
| `14b` | Years of life lost and adult mortality (45q15) | `ukr_yll_45q15_summary.rds` |
| `14c` | The loss with the counterfactual fixed at its point forecast, beside it drawn: the two tiers of uncertainty | `ukr_fixed_counterfactual_e0.rds` |
| `14d` | The loss with the population base of Donetsk and Luhansk drawn, PERT(0, 0.5, 1.0 million fewer), beside the base at the published estimates: the second tier of the structural uncertainty. The range is a judgement of the team, not a measurement, so it is reported beside the estimate and not inside it | `ukr_drawn_base_e0.rds` |
| `14e` | The intervals with the structural choices folded in: one alternative per choice sampled in every draw with equal weights and its shift from the central projection added, reported beside the interval of the drawn inputs. The weights are a convention, not a measurement | `ukr_structural_uncertainty.rds` |
| `17` | **The tables and figures of the paper**, main text and Supplement, and only those (`15`, which wrote the working tables of every analysis, is no longer run) | `outputs/` |
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
is projected through the full accounting. Civilian deaths are drawn by year
from a Beta-PERT on UCDP's low, best and high values; uncertainUCDP's
distributions of each UCDP event's deaths (Vesco et al. 2026) are the
alternative reported beside the estimate. Combatant deaths are the register's dead,
completed for registration lag, plus a share of the missing: the share the
register's own resolved cases recorded dead, which is held fixed, so that the
draws of the combatant total carry the registration-lag factors alone, rebuilt
under a draw of the weights on the release windows. The estimate is reported
between two bounds, none and all of the missing dead. The rule scenario of the
appendix (no evidence of life means dead) has its own inputs drawn beside it:
the prisoners of war among the missing, and where between none and its bound
the share alive outside captivity sits.
The western
net outflow follows from one weight, uniform between them, blending the two sources' four-year paths,
and Russia from one draw. Each of these is held across the four
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
| `00_setup`, `09`, `13s` | No source says how many of the missing are dead. The estimate is the resolution-based scenario: the missing not yet resolved are taken to resolve as the register's resolved cases did (dead, captivity and exits beyond list maintenance over the five windows, pooled), so the share those cases recorded dead is counted dead in every event year (`scenario_military()`; `data_inter/ukr_scenarios.rds`). It is reported between two bounds, none and all of the missing dead, which are the logical limits of that share and not plausible values. The simulation holds the share fixed and draws everything else. `13s` draws the estimate at every share from none to all, with the bounds and the estimate marked (figure 7, figure A6). Two further scenarios are kept for the appendix. The bodies scenario adds to the named dead the bodies awaiting identification (about 17,000 unidentified in September 2026, nine in ten taken as combatants listed as missing, which is a judgement; `official_figures.csv`), placed across event years by months since the event. The rule scenario counts a missing person as dead unless there is evidence of life, and its inputs are drawn beside the estimate. In the rule scenario the prisoners of war among the missing come from the official figures - about 7,000 held and 7,291 military personnel returned by February 2026, less the 11,325 the register records as prisoners - times the share of the returned whom the register had listed as missing (47.6% for 2024-2025 events, the mode). Of the rest, the share alive outside captivity lies between none and the share of the first resolutions outside captivity of the events of 2023-2025, pooled, that leave the register beyond list maintenance - beyond the rate at which the dead of the same cohort leave it in the same window - and is drawn uniformly between the two, one draw for the four event years. One bound serves every event year: the exits behind each later cohort's own bound fall in one window (72-100% of them in v18-v19), and the events of 2022 are not read on their own because their drop-outs are clean-ups of their records (09l); each cohort's own bound (6.6%, 1.9% and 6.2% for 2023-2025, 51% for 2022) is an alternative in table A19. The bound assumes the missing not yet resolved are no more often alive than those resolved, and holds only if the living leave the list at least as fast as the dead are identified. |
| `00_setup`, `09` | A person recorded as missing who is no longer listed — absent from a later release and from every release after it — is resolved alive. Before that, each release is searched under a corrected name or date of birth. A person whose first resolution is a return from captivity was a prisoner the register had not recorded, and is resolved alive, to captivity. |
| `08b`, `09` | The register's dead and missing are brought to the completeness events reach at four years, with the growth rates by lag held non-increasing within each status; the late registrations are reported apart. The missing are completed by the new persons listed only, since those who drop off the register are resolutions. How the missing resolve by months since the event is described, not projected: the estimate counts them by the identity above. |
| `11` | Imputed combatant deaths take the age–sex profile of the missing; registered deaths that of the confirmed dead. |
| `07_ohchr`, `11` | Civilian deaths take the age–sex profile of the deaths OHCHR verified, mostly in government-controlled territory. Step 13h gives the deaths beyond OHCHR's verified count, mostly in occupied territory and besieged cities, an older profile. |
| `00_setup` | Sex ratio at birth `srb = 1.06`. |
| `11`, `13` | WPP2024 publishes fertility only to 2023; the 2023 schedule is carried forward to 2024–2025. |
| `10`, `11` | Deaths and net outflow are spread evenly within each year (the mid-year convention), except in 2022, when their timing is taken from the months of the events: OHCHR's civilians killed by month, the register's dead and missing by month, and UNHCR's net border crossings by month. |
| `11` | Starts from continental Ukraine (SSSU region `cnt`): the whole country without Crimea and Sevastopol, **including** Donetsk and Luhansk. SSSU population estimates for those two regions assume complete registration, which has not held since 2015. |
| `13` | Conflict deaths are absolute counts, so net migration shrinks the denominator and raises the rates. Step 13 quantifies how much of the loss this accounts for; step 13c how far the loss moves across the range of each migration component. |
| `06`, `11` | The western net outflow is bracketed by source: each draw blends the border-crossing reading's four-year path with the register reading's by one weight, uniform between them, so the four-year total stays between the two sources' own totals (4.29–5.67 million). |

## Layout

```
code/              pipeline, numbered in execution order
  deprecated/      superseded scripts, kept for reference
data_input/        raw sources (large ones untracked — see its README)
data_inter/        intermediate results and cached extracts
data_deprecated/   retired inputs and outputs, untracked
outputs/           the tables and figures of the paper, written only by step 17
figures/           working figures of the analyses
  exploratory/     diagnostic plots from steps 02-14
  deprecated/      figures from previous versions, untracked
documents/         methodology reports and the draft paper, kept locally, untracked
```

## Figures and tables

Everything in the paper, main text and Supplement, comes from **step 17**
(`code/17_paper_outputs.R`) and goes to `outputs/`, one version of each item:
the folder is emptied and written again on every run. `outputs/manifest.csv`
lists every item with its number in the paper, its file, the files it is built
from, the number of draws and the commit. The analysis steps (02-14) still draw
their own diagnostic plots in `figures/exploratory/`; those are not the paper's.

| Item | File in `outputs/` | Content |
|---|---|---|
| Table 1 | `tables/table1_conflict_deaths_inputs.csv` | Conflict deaths as they enter the estimate: civilians at UCDP's adjusted low, best and high; combatants at the lower bound, the resolution-based scenario and the upper bound |
| Table 2 | `tables/table2_register_combatants.csv` | The register by year of the event: registered deaths, late registrations, the missing, and the missing counted dead in the resolution-based scenario |
| Table 3 | `tables/table3_deaths_by_type_sex.csv` | Conflict deaths 2022-2025 by type and sex, mean and 95% interval; late registrations are part of the registered combatants |
| Table 4 | `tables/table4_life_expectancy_loss.csv` | Life expectancy at birth expected and with war, the loss with its interval, and the loss at the two bounds |
| Figure 1 | `figures/fig1_register_dead_missing.png` | Combatants registered as dead and as missing, counts and shares |
| Figure 2 | `figures/fig2_death_rates_by_age.png` | Expected and with-war death rates by age, sex and year |
| Figure 3 | `figures/fig3_loss_of_life_expectancy.png` | The loss between the two bounds, with the resolution-based scenario and its interval |
| Figure 4 | `figures/fig4_loss_mortality_net_migration.png` | Shares of the loss due to conflict mortality and to net migration |
| Figure 5 | `figures/fig5_loss_by_type_of_death.png` | Shares of the loss by type of death |
| Table A1 | `tables/tableA1_sources.csv` | The sources as published (ACLED, UCDP, OHCHR, UALosses) and UCDP as adjusted, by year |
| Table A2 | `tables/tableA2_life_expectancy_2000_2025.csv` | Life expectancy at birth: observed 2000-2021, expected and with war 2022-2025 |
| Table A3 | `tables/tableA3_registration_lag.csv` | Combatant deaths and the male loss without the correction for registration lag and with it to 48, 60 and 72 months |
| Table A4 | `tables/tableA4_resolution_of_missing.csv` | What became of the missing whose status a later release settled: the source of the 69.6% |
| Table A5 | `tables/tableA5_net_migration.csv` | Net migration by year: border crossings, UNHCR's totals, their midpoint, and Russia |
| Table A6 | `tables/tableA6_deaths_by_year_type_sex.csv` | Conflict deaths by year, type and sex |
| Table A7 | `tables/tableA7a_ucdp_only_deaths.csv`, `tables/tableA7b_ucdp_only_loss.csv` | Every death from UCDP, in two panels: the deaths under uncertainUCDP beside the estimate, and the loss by sex with civilians, combatants and both replaced |
| Figure A1 | `figures/figA1_loss_by_share_of_missing.png` | The loss and combatant deaths at every share of the missing counted dead (drawn by 13s) |
| Figure A2 | `figures/figA2_loss_ucdp_only.png` | Figure 3 with the scenario that takes every death from UCDP |

Step 15 (`code/15_paper_figures_tables.R`) is no longer run by the pipeline. It
wrote the working tables of every analysis to `tables/`, which has been
removed; it is kept for reference and still runs by hand.
