# Author checks for the data-compilation methods

The draft, now integrated into `main.tex`, describes executable transformations in `analysis/00_raw_data.R` and the helpers in `R/functions.R`. It covers the 33 datasets included by `analysis/01_clean_data.R`, rather than every dataset ever processed. It does not establish that the present script produced the existing analysis objects: the raw files at the configured paths were unavailable here, so the raw workflow was not rerun. No analysis code was changed. These notes are not included in the compiled Supplementary Information.

## Checks affecting density, QMD, or observation years

1. **RAINFOR area units** — `analysis/00_raw_data.R:3173–3191` renames `number_stems` to `density` and `total_BA` to `ba`, without dividing either by `PlotArea`. Confirm from the input documentation whether both are already per hectare. If both are plot totals, their ratio still gives the same QMD, but density needs division by plot area; if their area bases differ, QMD also changes. The draft describes the assignments explicitly and does not claim a conversion that is absent from the code.

2. **Palanam and SERC diameter units** — The basal-area calculations at `analysis/00_raw_data.R:2491` and `analysis/00_raw_data.R:2594` treat the source `DBH` values as centimetres. Unlike Luquillo, SCBI, and Pasoh, there is no preceding multiplication by 0.1. Verify the units in these particular input files. The draft says that the code treats the values as centimetres, rather than claiming that the raw metadata confirm this.

3. **Wytham 2008** — `analysis/00_raw_data.R:2669–2688` reshapes `dbh_2008` but only status columns for 2010, 2016, and 2021. The 2008 rows therefore receive missing status and are discarded by `filter(status != "D")`. The draft lists the reshaped years but does not assert that 2008 survives. Confirm whether dropping that census is intended. @lauramarques

4. **Missing diameters and multiple stems** — `R/functions.R:333–375` counts input rows for density (`n() / plotsize`) but omits missing basal areas from the basal-area sum. Several dataset blocks do not explicitly remove missing diameter after selecting living trees. Nor does the common helper collapse multiple stem records to unique tree tags. Verify whether each source has already dealt with these cases. The draft explicitly distinguishes retained records/stems from unique biological individuals.

5. **Calendar year versus census** — `from_tree_data()` groups by calendar year, not the source census identifier. Pasoh derives year from each tree's exact measurement date (`analysis/00_raw_data.R:2787–2798`), while Palanam and SERC override only selected census years. If a census spans multiple years, partial within-census counts could be divided by the entire plot area. RAINFOR rounds decimal census dates. Confirm the intended time unit and whether within-census years were standardised upstream.

6. **Dataset-specific numerical exclusions** — The draft includes the actual filters: Palanam and SERC retain `logQMD > 1`; Pasoh retains `logDensity > 7.5`; ForestPlots retains `logQMD < 4.5`. The script provides little justification beyond outlier/error comments. Confirm the basis for these cutoffs before presenting them as quality-control decisions. Eight Swiss reserve plots are explicitly excluded: 7007028, 7005006, 7005003, 7005002, 7001001, 7001002, 7001003, and 7002003.

7. **Subplot boundaries and identifiers** — UNITO processes only parent plots `9_1`, `9_2`, and `9_3`, using breaks from `seq(min, max, by = 20)`. If the final edge is below the maximum coordinate, some records can receive a missing grid identifier; the code then groups this identifier. UNITBV groups by grid cell without including the original plot identifier. SERC uses a fixed 0–400 m grid. Confirm domain coverage, whether all cells have the assigned 0.04 ha area, and whether UNITBV contains a single parent plot as assumed in the comments.

8. **ForestPlots diameter fields** — `analysis/00_raw_data.R:2950–2954` takes the arithmetic mean of D0–D4 before converting to centimetres. Confirm the meaning of these fields. If they represent separate stems, averaging diameters would not preserve either the number of stems or summed basal area. The draft reports the calculation without assigning an unsupported interpretation to the fields. No live/dead filter is applied in this block or in IBERBAS; verify any upstream selection. @lauramarques

## Dataset identity and metadata

- **Table S1 versus current object names:** the draft retains the table identifiers and gives aliases where useful. Current script names include `fvabw` (table: `forst`), `unitbv` (`incds`), `npvbw` (`nbw`), `inrae_lessem` (`france`), and `fep_gre_stand` (`greece_stand`). `from_tree_data()` and related helpers derive the dataset field from the input object's name; consequently ForestPlots and RAINFOR currently receive `forestplots` and `rainfor_plots`, whereas downstream code sometimes checks `df_forestplots` and `df_rainfor`. Reconcile these names separately from the prose edit.
- **BCI:** explicitly excluded from the combined list at `analysis/01_clean_data.R:247`, consistent with the raw-script comment about unresolved management history. It is mentioned as an exclusion, not counted as a 34th constituent dataset.
  - beni: ok. Excluded from SI.
- **BNP:** the input is commented as species-level, but the executed code uses `from_stand_data()`, with no species summation. Confirm that the input is already one stand summary per plot-year. The draft follows the executed operation.
- **Projected coordinates:** the draft reports the CRS codes selected by `biomes_coords_utm()` (Spain: 25830; FVA-BW/LWF: 25832; NW-FVA: 31467; UNITO: 32632). Their correctness against the original coordinate metadata was not independently verified.
- **Site coordinates:** several ForestGEO datasets use one geographic location for all quadrats; RAINFOR fills missing plot coordinates by site-code prefix. The draft makes these approximations explicit.
- **Forestplots:** Information about data provenance is insufficient. The datasets are identified as coming from different countries. Are locations of plots known? Are standardised plot names/IDs used? Are YM and HZ ok with being listed as data providers?

## Management and workflow interpretation

- **Swiss NFI:** the block sums strictly increasing `LETZTENU` values, although the variable is described as years since intervention at each census (`analysis/00_raw_data.R:493–496`). This is not generally the same as time since intervention at the first census. Any occurrence of 999 overrides missing values in other censuses; otherwise missing or non-increasing values produce `management = 1`. There is no explicit 70-year threshold matching information in Table S1. Confirm the intended calculation and table wording. The manuscript describes the implemented screening variable without endorsing it as a first-census interval. (@lauramarques)
- **Australia:** all plots are assigned `management = 0`, even when `year_last_management` is present. The date changes `management_cat` only (missing date: 2; recorded date: 3), and the elapsed-time variable remains missing. Thus, recorded thinning does not exclude a plot from the common zero-value management filter, but category 3 excludes it from the stricter history subsets. Reconcile this with the table's within-monitoring description.
- **Source criteria versus executable filters:** the table caption now distinguishes reported source selection from row removal in the raw-data script. Most datasets receive `management = 0` as an assignment. Sweden retains the source flag; FIA uses the first retained condition record per plot and `RESERVCD - 1`; Swiss NFI derives the flag from its interval calculation. Verify upstream selection where the input itself contains no intervention history.
- **French INRAE LESSEM:** comments describe subtracting a last-harvest date from first measurement, but the executable block does not calculate the interval. `from_tree_data()` requires the elapsed-time column to exist in the joined input. Verify that the supplied metadata contain it and document its derivation; otherwise this block cannot run as written.
- **Fixed intervals and special codes:** Norway (64 years), CZU (69), IBERBAS (70), and SCBI (53) receive dataset-wide constants rather than intervals derived from plot-specific dates. Confirm their provenance and relevance to the first census; the CZU comment explicitly requests clarification. Uholka is category 1 with marker 999 and therefore passes the 30- and 100-year thresholds but not the category-2-only filter. The manuscript distinguishes 999 from a measured interval.
- **Comments versus assignments:** UNITO and WULS are category 2 with missing intervals, rather than the possible 50- and 20-year intervals in comments. NPV-BW and BNP use date differences rather than the commented 42- and 6-year constants. NW-FVA also uses date differences; its comments flag disagreement with a 52-year interval in supplementary information. BNP is subsequently reclassified to category 2 in `analysis/02_create_filters.R:77–78`. LWF remains category 3 despite the table's reserve designation. Dataset paragraphs now describe these distinctions explicitly.
- **Two filtering routes:** `data_unm_fc()` in `R/functions.R`, called by `analysis/01_clean_data.R`, retains intervals of at least 30 years **or missing intervals**, alongside `management == 0`. The current `analysis/02_create_filters.R` instead reads `df_all.rds` and applies `management == 0` without a common age threshold. Its downstream sensitivity filters in `R/analyse_biome.R` retain `(category 1 AND interval >= 30 or 100 years) OR category 2`, or category 2 alone. The manuscript describes this latter route; confirm that the reported results use it.
- The FIA download/read branch contains a malformed tree-file `read.csv()` call (`analysis/00_raw_data.R:302`), and its cache existence check uses a different filename from the one subsequently loaded. The combining script assigns `data_all` and then continues from `df_all` (`analysis/01_clean_data.R:221–260`). These existing execution issues were not changed by the writing task.
- The ForestPlots/RAINFOR deduplication comment in `analysis/01_clean_data.R` is not sufficient evidence of cross-source deduplication: `plotID` is prefixed with dataset before `distinct(country, plotID, year)`. The draft makes no claim that overlapping plots across sources were removed.

## Changes to the previous methods text

- Moved the common calculations into Data compilation and removed the separate, duplicative Data processing section.
- Distinguished a count of retained records (`n`) from density per hectare (`N`), included plot area in the basal-area equation, and included the factor of 100 needed for QMD in centimetres.
- Retained the equation labels `eq:ba`, `eq:qmd`, and `eq:qmd2`.
- Described source-specific diameter filters without asserting a common individual-tree diameter threshold. The later QMD ≥ 10 cm requirement is a stand-level filter.
- Used source-location comments before each dataset subsection for traceability. No new external references or unverified census counts were added.

# Remaining revisions

Remaining revisions of `00_raw_data.R`.

## To do

- Do not calculate the logarithm of N and QMD nd do not append this as a new variable to the data at this step.
- For the tree-level aggregation, do not omit missing individual basal areas from basal-area sums. Instead, impute missing basal area values based on the diameter at breast height (DBH) for the respective tree. If this variable is missing as well, do not calculate the basal area and make it `NA` instead.
- Avoid any use of the factor 0.00007853982. Instead, use R's internal `pi` explicitly in all respective calculations. For example $0.00007853982 = (1/(100\cdot\sqrt{4/\pi}))^2$.
- For US FIA data processing, explicitly remove plots that were not classified as forest reserve.
- Add column in table providing information of data provider. Name should be part of the author list.
- Add to constituent dataset tables information as additional columns:
  - Names of data providers or reference if taken from published sources
- Harmonise names, also in table and plot labels: 
  - `fep_gre_stand` instead of `greece_stand`
  - `inrae_lessem` instead of `france`
  - `fvabw` instead of `forst`
  - `unitbv` instead of `incds`
  - `lwf` instead of `lwf_tree`
  - `nbw` instead of `npvbw`
  - `nwfva` instead of `nwfva_tree`
  - `tuzvo` instead of `tuzvo_tree`
  - `ul` instead of `ul_tree`
  - `bci` instead of `scbi`
  - `palanan` instead of `palanam`
  - `forestplots` instead of `df_forestplots`
<!-- - Make the BA-QMD conversion into a function, using the exact value of $\pi$ (`100 * sqrt((4 * ba)/(pi * density))` instead of `QMD = sqrt(ba / (0.0000785 * density))`), and apply it to each constituent dataset. -->
<!-- - Avoid creating new variables for `logQMD` and `logDensity`. These are now calculated with each model fit later in the workflow and specified as a "recipe". -->

## To check Laura

(@lauramarques)

- Confirm processing of Norwegian National Forest Inventory: All plots were assumed unmanaged and were given `management_cat = 1` and a fixed elapsed-time value of 64 years to every plot.
- Confirm: Uholka is changed from `management_cat = 1` and `management_since_census1_yrs = 999` to `management_cat = 2` (primary forest).
- Should EuFoRia plots all categorised as `management_cat = 2`? Now, at least some are `management_cat = 1` and time since last management is specified (inlcuding `bnp`, `iberbas`, `lwf`).
- EuFoRia Unito: comment in script `00_raw_data.R`: "In Table S1, assign to 50 years. Need to confirm with data providers." Now set to `management_cat = 2` (primary/old growth).
- EuFoRia urk: Comment in script: "In Table S1, assign to NA years. It can be considered a primary/old-growth forests.". Therefore, I changed this to `management_cat = 2` (primary/old growth). Ok?
- Treatment of `management_since_census1_yrs` in `aus_plots` was revised. Now interpreting `year_last_management`: if *after* the first available census per plot, then set `management = 1`, otherwise `management = 0` and set `management_since_census1_yrs = min(year, na.rm = TRUE) - year_last_management`. Please confirm. It was unclear why all plots for which `year_last_management` was *not* NA were set as `management_cat = 3` (unknown).
