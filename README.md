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


## Repository use

This repository includes the R scripts to analyse the links between growth and biomass from global forest inventory data.

The repository is organized as follows:
- analysis - includes the scripts for data analysis.
- data - includes model outputs.
- manuscript - includes figures.
- vignettes - contains dynamic notebooks, i.e. R markdown files.
