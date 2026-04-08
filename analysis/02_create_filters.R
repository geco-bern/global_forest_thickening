# Implements additional data filters to retain only self-thinning affected plots
# Workflow overview:
#  
# By forest plot
# - Use unmanaged plots only (read data file created by first step, done by Laura)
# - Remove disturbance-affected plots (simultaneous QMD and N drop)
# - Remove ingrowth-affected plots (simultaneous QMD drop and N increase)
#
# By dataset (major)
# Small dataset should be grouped to new classification `dataset_major`
# - Remove years (at head of time series) with substantially lower QMD
#
# By biome
# - Filter based on biome-specific slope
#     - overall slope by plot
#     - for long plot-level time series (N≥5), remove full plot if segment (either first half or second half) has outlier slope
# 
# Workflow/code design:
# - Each step returns a flat data frame with the filter criterium added as an additional column, but no rows removed.
#  
## Loading ------------------------------------------------------------------
library(tidyverse)
library(here)
library(segmented)

source(here("R/identify_disturbed_plots.R"))
source(here("R/identify_ingrowth_plots.R"))
source(here("R/identify_badbins.R"))

# df_all contains information of management:
# - management_since_census1_yrs (before called years_last_management)
# - management_cat = 1,2 or 3.
# - management = 0 (indicating no management between census)
# 
# Among management_cat:
# 1. Recorded history
#    The time since the last management before the first census is known and recorded.
#    the number of years since the last management before the first census, or
#    the calendar year of the last management intervention.
# 2. Pristine/primary/old-growth/protected
#    The time since the last management before the first census is unknown. The forest was considered pristine, primary, old-growth, protected, or minimally used at the time of the first census.
# 3. Unrecorded history
#    The time since the last management before the first census is unknown. The forest may have a history of substantial past wood harvesting and/or visible signs of forest use before the first census.

# Read data prepared with 01_clean_data.R, line 290
df_all <- read_rds(here("data/inputs/df_all.rds"))

# Apply filters relevant for all analyses
df_unm <- df_all |>
    # filter for min qmd
    filter(QMD >= 10) |>
  
    # filter for forest type
    filter(type == "Forest") |>
  
    # filter for unmanaged plots (between censi)
    filter(management == 0) |>
  
    # Filter by minimum 3 censi
    group_by(dataset, plotID) |>
    mutate(n_census_unm = n()) |>
    filter(n_census_unm >= 3) |>
  
    # remove plots with no change in ln(N)
    mutate(var_logdensity = diff(range(logDensity))) |>
    filter(var_logdensity > 0.001) |> 
  
    # clean up
    ungroup() |>
    relocate(n_census_unm, .after = n_census)

# optionally subset
do_subset_primary <- FALSE  # <- manually adjust here

if (do_subset_primary){
  suffix_subset <- "_SUBSET"

  table_s1 <- read_csv(here("data/table_s1_exported.csv"))
  vec_datasets_subset <- table_s1 |>
    filter(`Can be considered primary/old growth forest?` == "YES") |>
    pull(Dataset)

  df_unm <- df_unm |>
    filter(dataset %in% vec_datasets_subset)

} else {
  suffix_subset <- ""
}

## By forest plot --------------------------------------------------------------

### Identify disturbed plots ---------------------------------------------------
# adds column 'disturbed'
df_unm <- df_unm |>
  identify_disturbed_plots()

### Identify ingrowth-affected plots -------------------------------------------
# adds column 'ingrowth'
df_unm <- df_unm |>
  identify_ingrowth_plots()

## By (major) dataset ----------------------------------------------------------
### Define major datasets ----------
# grouping into major based on analysis/01_clean_data.R
df_unm <- df_unm |> 
  mutate(
    dataset_major = ifelse(
      dataset %in% c("bnp", "czu", "forst", "iberbas", "incds", "lwf", "nbw", "nfr", "nwfva", "tuzvo", "ul", "unito", "urk", "wuls", "tuzvo_tree", "nwfva_tree", "ul_tree"),
      "euforia",
      ifelse(
        dataset %in% c("luquillo", "bci", "scbi", "palanam", "serc", "pasoh", "mudumalai"),
        "forestgeo",
        ifelse(
          dataset == "df_forestplots", 
          "forestplots",
          ifelse(
            dataset == "df_rainfor",
            "rainfor",
            dataset
          )
        )
      )
    )
  )

tmp <- df_unm |> 
  group_by(dataset_major) |> 
  summarise(n = n())

View(tmp)

### Identify years (in 10-year bins) based on deviating QMD distribution ------------
# adds column 'badqmdbin'
tmp_badbins <- df_unm |> 
  # filter(dataset_major == "aus_plots") |> 
  group_by(dataset_major) |> 
  nest() |> 
  mutate(out = purrr::map2(data, dataset_major, ~identify_badbins(.x, .y))) |> 
  mutate(
    gg = purrr::map(out, "gg"),
    data = purrr::map(out, "df")
  ) |> 
  dplyr::select(-out)

gg_qmdfilter <- cowplot::plot_grid(
  plotlist = tmp_badbins$gg,
  ncol = 4,
  labels = letters[1:18]
)

ggsave(
  here("fig/qmdfilter.pdf"),
  plot = gg_qmdfilter,
  width = 12,
  height = 10
)

df_unm <- tmp_badbins |> 
  dplyr::select(-gg) |> 
  unnest(data) |> 
  ungroup()

### Identify non-self-thinning plots based on outlying slope ----------
# separately for each (major) dataset
# adds column 'badslope'
tmp_slopefilter <- df_unm |> 
  group_by(dataset_major) |> 
  nest() |> 
  mutate(out = purrr::map2(data, dataset_major, ~filter_stl_slope2(.x, .y))) |> 
  mutate(
    gg = purrr::map(out, "gg"),
    data = purrr::map(out, "df")
  ) |> 
  dplyr::select(-out)

gg_slopefilter <- cowplot::plot_grid(
  plotlist = tmp_slopefilter$gg,
  ncol = 4,
  labels = letters[1:18]
)

ggsave(
  here("fig/slopefilter.pdf"),
  plot = gg_slopefilter,
  width = 12,
  height = 10
)

df_unm <- tmp_slopefilter |> 
  dplyr::select(dataset_major, data) |> 
  unnest(data)

write_rds(df_unm, here("data/inputs/df_unm_withfilters.rds"))



