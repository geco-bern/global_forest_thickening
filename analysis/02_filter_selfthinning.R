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

## Loading ------------------------------------------------------------------
library(tidyverse)
library(here)

source(here("R/identify_disturbed_plots.R"))
source(here("R/identify_ingrowth_plots.R"))
source(here("R/identify_badbins.R"))

# Load data (prepared in analysis/01_clean_data.R)
# This is after applying function data_unm_fc(), see R/functions.R:
# - filter for minimum QMD (>=10)
# - filter for forest type
# - filter for unmanaged plots
# - filter for years since last management >= 30 (or no info available)
# - filter by minimum 3 censuses
# - remove plots with no change in ln(N)
data_unm <- read_rds(here("data/inputs/data_unm.rds"))

# optionally subset
do_subset_primary <- FALSE  # <- manually adjust here

if (do_subset_primary){
  suffix_subset <- "_SUBSET"

  table_s1 <- read_csv(here("data/table_s1_exported.csv"))
  vec_datasets_subset <- table_s1 |>
    filter(`Can be considered primary/old growth forest?` == "YES") |>
    pull(Dataset)

  data_unm <- data_unm |>
    filter(dataset %in% vec_datasets_subset)

} else {
  suffix_subset <- ""
}

## By forest plot --------------------------------------------------------------

### Remove disturbed plots ---------------------------------------------------
data_unm_und <- data_unm |>
  identify_disturbed_plots()

### Remove ingrowth-affected plots -------------------------------------------
data_unm_und_noin <- data_unm_und |>
  identify_ingrowth_plots()

## By (major) dataset ----------------------------------------------------------
### Define major datasets ----------
# grouping into major based on analysis/01_clean_data.R
data_unm_und_noin <- data_unm_und_noin |> 
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

df_summ_dataset_major <- data_unm_und_noin |> 
  group_by(dataset_major) |> 
  summarise(n = n())

View(df_summ_dataset_major)

### Remove years based on QMD distribution ------------
# (at head or tail of time series)
# Compare pairwise between N-year bins:
data_unm_und_noin_qmdfilt <- data_unm_und_noin |> 
  # filter(dataset_major == "aus_plots") |> 
  group_by(dataset_major) |> 
  nest() |> 
  mutate(out = purrr::map2(data, dataset_major, ~identify_badbins(.x, .y))) |> 
  mutate(
    gg_badbins = purrr::map(out, "gg"),
    data = purrr::map(out, "df_binned")
)

gg_qmdfilter <- cowplot::plot_grid(
  plotlist = data_unm_und_noin_qmdfilt$gg_badbins,
  ncol = 4,
  labels = letters[1:18]
)

ggsave(
  here("fig/qmdfilter.pdf"),
  plot = gg_qmdfilter,
  width = 12,
  height = 10
)

