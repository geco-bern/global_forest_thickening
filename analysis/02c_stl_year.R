# Additional tests, visualisations, and alternative biome grouping

# Load packages ----------------------------------------------------------------
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(rFIA)
library(lme4)
library(lmerTest)
library(ggeffects)
library(effects)
library(sjPlot)
library(measurements)
library(lqmm)
library(ggforce)
library(MuMIn)
library(DescTools)
library(here)
library(viridis)
library(purrr)
library(rsample)
library(cowplot)

# Load functions ---------------------------------------------------------------
source(here("R/functions.R"))
source(here("R/plot_stl_bybiome.R"))
source(here("R/identify_disturbed_plots.R"))
source(here("R/get_breaks.R"))
source(here("R/plot_lqmm_bybiome.R"))
source(here("R/calc_lqmm_byqmdbin.R"))
source(here("R/create_table_latex.R"))
source(here("R/process_cite_lines.R"))
source(here("R/wrap_fit_lqmm.R"))
source(here("R/calc_percent_change.R"))
source(here("R/plot_disturbed.R"))
source(here("R/filter_stl_slope.R"))

# Load data -----
data_unm <- readRDS(here("data/inputs/data_unm.rds"))

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

# Biome 1 + 2 Tropical & Subtropical Broadleaf Forests ----------------------
# Combining all tropical and subtropical
data_unm_biome <- read_rds(
  file = here(paste0("data/data_unm_undist_biome1", suffix_subset, ".rds"))) |>
  bind_rows(
    read_rds(file = here(paste0("data/data_unm_undist_biome2", suffix_subset, ".rds")))
  )

## Plot differences between biomes ------------------------
ggplot() +
  geom_point(
    data = data_unm_biome,
    aes(x = logQMD, y = logDensity, color = biome),
    alpha = 0.4
  ) +
  labs(
    x = expression(ln(QMD)),
    y = expression(ln(italic(N))),
    title = "Tropical forests",
    color  = NULL
  ) +
  scale_color_manual(values = c(viridis(3)[2], viridis(3)[3])) +
  theme_bw() +
  scale_x_continuous(limits = c(2.2, 4.7), breaks = seq(3,4,1)) +
  scale_y_continuous(limits = c(2.9,9.3), breaks = seq(4,8,2)) +
  theme(
    legend.position = "bottom"
  )

## LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.9, # c(0.75, 0.90),
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(
    LP_max_iter = 500, # inner loop iterations
    LP_tol_ll   = 1e-05, # inner loop tolerance
    startQR     = TRUE
  )
)
# summary(fit_lqmm)

percent_change <- calc_percent_change(data_unm_biome, fit_lqmm)

gg_lqmm_biome1 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("a") ~ ~"Tropical & Subtropical Broadleaf Forests")
)
gg_lqmm_biome1

### Within QMD bins ------------------------------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome1_byqmdbin <- plot_lqmm_byqmdbin(
  df_lqmm_byqmdbin$df,
  df_lqmm_byqmdbin_including_disturbed$df
)

gg_lqmm_biome1_both <- cowplot::plot_grid(
  gg_lqmm_biome1,
  gg_lqmm_biome1_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  align = "v",
  # labels = c("",  "g"),
  label_y = 1.1
)
gg_lqmm_biome1_both


# Biome 4 + 5 Temperate Forests ----------------------
# Combining all broadleaf and coniferous
# Read unmanaged and undisturbed data. 
data_unm_biome <- read_rds(
  file = here(paste0("data/data_unm_undist_biome4", suffix_subset, ".rds"))) |>
  bind_rows(
    read_rds(file = here(paste0("data/data_unm_undist_biome5", suffix_subset, ".rds")))
  ) |> 
  group_by(dataset) |> 
  nest()

## Plot subset ------------------------
analyse_stl_bydataset <- function(dataset_name, df, multiplier_sd){

  ### Fit 90% quantile regression ------------------------------------------------
  # and its temporal shift on all data
  set.seed(123)
  fit_lqmm_dataset <- lqmm(
    logDensity ~ logQMD_sc + year_sc,
    random = ~1,
    group = plotID,
    tau = 0.9,
    data = df,
    type = "normal",
    control = lqmmControl(
      LP_max_iter = 500, # inner loop iterations
      LP_tol_ll   = 1e-05, # inner loop tolerance
      startQR     = TRUE
    )
  )

  slope_lqmm <- coef(fit_lqmm_dataset)["logQMD_sc"]

  ### Plot STL from LQMM ---------------------------------------------------------
  gg_stl_lqmm <- plot_lqmm_bybiome(
    df,
    fit_lqmm_dataset,
    name = dataset_name
  )
  gg_stl_lqmm

  ### Within QMD bins ------------------------------------------------------------
  # Test whether upward shift of 90% quantile is significant within logQMD-bins
  # returns data frame with pval indicating significance level of a positive
  # effect of year.
  df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(df)
  gg_lqmm_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df)

  ## Plot QMD over time
  gg_qmd <- df |>
    mutate(inbin = cut(year, breaks = seq(1950, 2025, by = 5))) |>
    group_by(inbin) |> 
    ggplot(aes(inbin, logQMD)) +
    geom_boxplot(fill = "grey70") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90)
    )

  # Filters data based on whether a plot's regression slope of logN ~ logQMD is
  # within a certain margin of the regression slope obtained from the quantile
  # regression mixed effects model.
  df_byplot <- df |>
    group_by(plotID) |>
    nest() |>
    mutate(nobs = purrr::map_int(data, ~nrow(.))) |>
    filter(nobs >= 3) |>  # minimum satisfied by design of data selection, but here for safety
    mutate(linmod = purrr::map(data, ~lm(logDensity ~ logQMD_sc, data = .))) |>
    mutate(slope = purrr::map_dbl(linmod, ~coef(.)[2]))

  # remove gg_outliers
  limits <- quantile(df_byplot$slope, probs = c(0.01, 0.99))
  df_byplot <- df_byplot |> 
    filter(slope > limits[1] & slope < limits[2])
  sd_slopes <- sd(df_byplot$slope)

  # plot distribution of slopes
  gg_slopes <- df_byplot |>
    ggplot() +
    geom_histogram(
      aes(slope, after_stat(density)),
      fill = "grey70",
      color = "black"
    ) +
    geom_vline(
      xintercept = slope_lqmm,
      color = "red"
    ) +
    geom_vline(
      xintercept = c(slope_lqmm - multiplier_sd * sd_slopes, slope_lqmm + multiplier_sd * sd_slopes),
      color = "red",
      linetype = "dashed"
    ) +
    xlim(-4, 3) +
    theme_classic() +
    labs(
      title = dataset_name,
      x = "Slope (unitless)",
      y = "Density"
    )


  ### Apply slope filter ---------------------------------------------------------
  df_slopefiltered <- filter_stl_slope(df, fit_lqmm_dataset, multiplier_sd = multiplier_sd) |> 
    unnest(data)

  # fit again
  set.seed(123)
  fit_lqmm_dataset_slopefiltered <- lqmm(
    logDensity ~ logQMD_sc + year_sc,
    random = ~1,
    group = plotID,
    tau = 0.9,
    data = df_slopefiltered,
    type = "normal",
    control = lqmmControl(
      LP_max_iter = 500, # inner loop iterations
      LP_tol_ll   = 1e-05, # inner loop tolerance
      startQR     = TRUE
    )
  )

  # and plot
  gg_stl_lqmm <- plot_lqmm_bybiome(
    df_slopefiltered,
    fit_lqmm_dataset_slopefiltered,
    name = dataset_name
  )
  gg_stl_lqmm

  # Test whether upward shift of 90% quantile is significant within logQMD-bins
  # returns data frame with pval indicating significance level of a positive
  # effect of year.
  df_lqmm_byqmdbin_slopefiltered <- calc_lqmm_byqmdbin(df_slopefiltered)
  gg_lqmm_byqmdbin_slopefiltered <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df)

  # plot with points and lines connecting points
  gg_stl_slopefiltered <- ggplot() +
    geom_line(
      data = df_slopefiltered,
      aes(x = logQMD, y = logDensity, group = plotID),
      alpha = 300/nrow(df_slopefiltered)
    ) +
    labs(
      title = dataset_name,
      x = expression(ln(QMD)),
      y = expression(ln(italic(N))),
      color  = NULL
    ) +
    scale_color_viridis_c() +
    theme_classic()

  gg_stl_slopefiltered

  gg_out <- cowplot::plot_grid(
    gg_stl_lqmm, 
    gg_slopes,
    gg_stl_slopefiltered,
    gg_lqmm_byqmdbin,
    gg_qmd,
    gg_lqmm_byqmdbin_slopefiltered,
    ncol = 3
  )

  ggsave(
    here(paste0("manuscript/figures/panel_", dataset_name, ".pdf")),
    plot = gg_out,
    width = 12,
    height = 8
  )

  return(gg_out)

}

tmp <- data_unm_biome |> 
  slice(1:5) |> 
  mutate(gg_combined = purrr::map2(dataset, data, ~analyse_stl_bydataset(.x, .y, multiplier_sd = 1.0)))

multiplier_sd <- 1.0
tmp <- data_unm_biome$data[[2]]
dataset_name <- data_unm_biome$dataset[[2]]


xxxxxx


# just Swiss Forest Reserves
data_unm_biome <- data_unm_biome |>
  filter(dataset == "nfr_swi")



# get fitted quantile regression model object from corresponding biome
fit_lqmm <- read_rds(here("data/gg_outputs/fit_lqmm_biome4.rds"))
slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]
sd_slopes <- sd(df_byplot$slope)


data_unm_biome_slopefilter <- df_byplot |>

  # remove gg_outliers (top and bottom 1% of data)
  filter(slope > quantile(df_byplot$slope, 0.01) & slope < quantile(df_byplot$slope, 0.99)) |>

  # subset data based on slope: must be within a certain margin of slope obtained from quantile regression
  # margin defined as half a standard deviation across all slopes
  mutate(filter_slope = slope < slope_lqmm + 0.5 * sd_slopes & slope > slope_lqmm - 0.5 * sd_slopes) |>
  filter(filter_slope) |>
  unnest(data)

# plot data after applying the slope filter
gg_stl_subset_slopefilter <- ggplot() +
  geom_point(
    data = data_unm_biome_slopefilter,
    aes(x = logQMD, y = logDensity, group = plotID, color = year),
    size = 0.7
  ) +
  geom_line(
    data = data_unm_biome_slopefilter,
    aes(x = logQMD, y = logDensity, group = plotID),
    alpha = 0.4
  ) +
  labs(
    subtitle = "Swiss Forest Reserves, broadleaf subset",
    title = "Slope-filtered",
    x = expression(ln(QMD)),
    y = expression(ln(italic(N))),
    color  = NULL
  ) +
  scale_color_viridis_c() +
  theme_bw()

gg_stl_subset_slopefilter

cowplot::plot_grid(
  gg_stl_subset,
  gg_slopes_subset,
  gg_stl_subset_slopefilter,
  labels = letters[1:3],
  nrow = 1,
  rel_widths = c(1, 0.7, 1)
)

ggsave(here("manuscript/figures/demo_swiss_forest_reserves.pdf"), width = 12, height = 4)

## LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.9, # c(0.75, 0.90),
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(
    LP_max_iter = 5000, # inner loop iterations
    LP_tol_ll   = 1e-05, # inner loop tolerance
    startQR     = TRUE
  )
)
# summary(fit_lqmm)

percent_change <- calc_percent_change(data_unm_biome, fit_lqmm)

gg_lqmm_biome2 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("b") ~ ~"Temperate Forests")
)
gg_lqmm_biome2

### Within QMD bins ------------------------------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome2_byqmdbin <- plot_lqmm_byqmdbin(
  df_lqmm_byqmdbin$df,
  df_lqmm_byqmdbin_including_disturbed$df
)

gg_lqmm_biome2_both <- cowplot::plot_grid(
  gg_lqmm_biome2,
  gg_lqmm_biome2_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  align = "v",
  # labels = c("",  "g"),
  label_y = 1.1
)
gg_lqmm_biome2_both

