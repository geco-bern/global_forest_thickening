# Increased mature forests' carbon sink across biomes

## Abstract

Global forest growth has increased in recent decades, but its role in sustaining biomass gains remains uncertain, especially in mature forests. Biomass increases may result from accelerated growth but also from forest recovery, and disentangling these effects is difficult due to confounding factors like stand history and environmental effects. This study explores the growth-biomass relationship by examining shifts in stand carrying capacity. We find increasing stand density for a given mean tree size in unmanaged closed-canopy forests in five major forest biomes over the past decades. We tested if that change in the size-density relationship was influenced by different environmental drivers. A higher shift in the stand carrying capacity is observed in drier and carbon-poor areas and regions with higher content of nitrogen and phosphorus. Based on the observed upward shift we estimated the mature forest carbon sink in aboveground biomass. Our analyses align with previous reports and model predictions of persistent biomass increases and provide a better understanding of whether and how growth enhancements drive increasing land carbon storage.

![](https://github.com/geco-bern/global_forest_thickening/blob/master/manuscript/figures/fig_4.png)

## Workflow

### Prepare data

`analysis/00_raw_data.R`

`analysis/00_collect_global_drivers.R`

### Process data

`analysis/01_clean_data.R`

- Initial data selection, variable definition, and unit conversions.

`analysis/02_create_filters.R`

- Create filters to reduce to self-thinning dominated stands.

### Analyse data: negative size-density relationship (SDR)

`analysis/03_analyse_filter_effects.R`

- Analyses sensitivity of SDR to management-history by constituent dataset.

`analysis/04_stl_year.R`

- Analyses SDR and sensitivity to management-history by biome.

### Model preprocessing

The filtered input data contain predictors in their original units.
`fit_model()` uses a `recipes::step_normalize()` recipe to center and scale
`logQMD` and `year` on the complete training rows of each filter stage, biome,
QMD bin, or bootstrap sample. The response and plot identifiers are not scaled.
The fitted recipe is stored as `attr(fit, "recipe")`; the returned object remains
an `lqmm` or `lmerMod` model, with standardized coefficients named `logQMD_sc`
and `year_sc`.

Use `predict_model(fit, new_data, ...)` with predictors in original units to apply
the fitted preprocessing before prediction. `coef_original_units(fit)` converts
the fixed effects, including the intercept, back to original units.
The local bootstrap helper in `analyse_biome()` calls `fit_model()` and returns
these converted coefficients, so its `year` estimate is per calendar year and annual percentage change is
`100 * expm1(estimate)`. Bootstrap resampling remains restricted to `no_badslope`,
with 500 replicates by default (`analyse_biome(..., n_boot = 500)`).

Re-run `analysis/02_create_filters.R` and `analysis/04_stl_year.R` to regenerate
the input data and fitted results. Previously saved models do not contain the
trained recipes. The model checks can be run with
`Rscript tests/test_model_preprocessing.R`; they require `recipes`, `lqmm`,
`rsample`, `lme4`, `testthat`, and the plotting packages used by the analysis.


## Repository use

This repository includes the R scripts to analyse the links between growth and biomass from global forest inventory data.

The repository is organized as follows:
- analysis - includes the scripts for data analysis.
- data - includes model outputs.
- manuscript - includes figures.
- vignettes - contains dynamic notebooks, i.e. R markdown files.
