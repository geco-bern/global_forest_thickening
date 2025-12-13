# Global C sink ----------------------------------------------------------------
# Calculates the global C sink as a consequence of the STL change, scaled up to
# global forests with environmental covariates. Uncertainty estimation with
# bootstrapping

library(dplyr)
library(readr)
library(lme4)
library(purrr)
library(tidyr)
library(here)
library(rsample)
library(recipes)
library(multidplyr)
library(tictoc)
library(tidyterra)
library(rnaturalearth)
library(ggplot2)
library(scales)

## Load data -------------------------------------------------------------------
### Plot data, filtered upper edge ---------------------------------------------
use_slopefilter <- TRUE

if (use_slopefilter){

  lab_filter <- "slopefilter"

  # filtered by slope
  data_forest_plots <- read_rds(here("data/data_unm_slopefilter.rds")) |>
    mutate(NQMD2 = density * QMD^2)

} else {

  lab_filter <- "quantilefilter"

  # Load and engineer data with environmental factors
  # plot-level data for model fitting, filtered based on quantiles by bin
  data_forest_plots <- read_rds(here::here("data/inputs/data_fil75_biomes.rds")) |>
    mutate(NQMD2 = density * QMD^2)

}

### Maps of environmental covariates -------------------------------------------
# Load data for upscaling: maps of environmental factors
grid_drivers <- read_rds(here("data/df_grid.rds")) |>
  as_tibble() |>
  mutate(lon_i = round(lon * 4), lat_i = round(lat * 4))

## Function definitions --------------
### Fit LMM for STL per bootstrap --------
# this is hard-coded here based on model selection results from 03_env_drivers.R
fit_stl_byboot <- function(df) {
  vars_to_scale <- c("logQMD", "year", "tavg", "ai", "ndep", "ORGC", "CNrt", "PBR")

  rec <- recipe(
    logDensity ~ logQMD + year + tavg + ai + ndep + ORGC + CNrt + PBR + dataset + plotID + species,
    data = df |>
      drop_na(
        all_of(
          c(
            "logDensity",
            "logQMD",
            "year",
            "tavg",
            "ai",
            "ndep",
            "ORGC",
            "CNrt",
            "PBR",
            "dataset",
            "plotID",
            "species"
          )
        )
      )
  ) %>%
    step_center(all_of(vars_to_scale)) %>%
    step_scale(all_of(vars_to_scale))

  # fit the recipe (stores the sds and sds)
  rec_prep <- prep(rec)

  # get the preprocessed training data
  df_scaled <- bake(rec_prep, new_data = NULL)

  # fit your model using these preprocessed columns
  model <- lmer(
    logDensity ~ logQMD +
      year * tavg +
      year * ai +
      year * ndep +
      year * ORGC +
      year * CNrt +
      year * PBR +
      (1 | dataset / plotID) + (1 | species),
    data = df_scaled
  )

  # useful to combine
  model_bundle <- list(recipe = rec_prep, model = model)

  return(model_bundle)
}

### Fit LMM for biomass per bootstrap --------
fit_biomass_byboot <- function(df) {
  # biomass data is in kg-DM/ha
  df <- df |>
    drop_na(
      all_of(
        c(
          "biomass",
          "NQMD2",
          "dataset",
          "plotID"
        )
      )
    )
  if (nrow(df) == 0) {
    return(NA)
  } else {
    model <- lmer(
      biomass ~ NQMD2 + 0 + (1 | dataset / plotID),
      data = df
    )
    return(model)
  }
}

### Predict density change (dn) with fitted model
predict_dn <- function(bundle_stl, vec_qmd, df) {
  # sample QMD from empirical distribution and add dummy random factor levels
  # (ignored for prediction but required to avoid error)
  df <- df |>
    mutate(
      qmd = sample(vec_qmd, nrow(df), replace = TRUE),
      year = 2000,
      dataset = "dummy_dataset",
      plotID = "dummy_plot",
      species = "dummy_species"
    ) |>
    mutate(logQMD = log(qmd))

  # apply same preprocessing as for the STL data
  new_df_scaled <- bake(bundle_stl$recipe, new_data = df)

  # predict N, given QMD, year, and environmental factors for all forest gridcells globally
  df$logDensity_0 <- predict(
    bundle_stl$model,
    newdata = new_df_scaled,
    re.form = NA
  ) # Predict ignoring random effects

  df <- df |>
    mutate(year = 2001)

  # apply same preprocessing as for the STL data
  new_df_scaled <- bake(bundle_stl$recipe, new_data = df)

  # predict N, given QMD, year, and environmental factors for all forest gridcells globally
  df$logDensity_1 <- predict(
    bundle_stl$model,
    newdata = new_df_scaled,
    re.form = NA
  ) # Predict ignoring random effects

  df <- df |>
    mutate(dn = exp(logDensity_1) - exp(logDensity_0)) |>
    mutate(dnqmd2 = dn * qmd^2) |>
    select(lon_i, lat_i, dn, dnqmd2)

  return(df)
}

predict_db <- function(df, model) {
  if (identical(model, NA)) {
    df <- df |>
      mutate(db = NA)
  } else {
    # add dummy random factor levels
    # (ignored for prediction but required to avoid error)
    df <- df |>
      mutate(
        dataset = "dummy_dataset",
        plotID  = "dummy_plot"
      )

    # predict N, given QMD, year, and environmental factors for all forest gridcells globally
    df$db <- predict(
      model,
      newdata = df |>
        rename(NQMD2 = dnqmd2), # predict one with difference rather than twice with absolute
      re.form = NA # Predict ignoring random effects
    )

    df <- df |>
      select(-dataset, -plotID)
  }

  return(df)
}

calc_csink_global <- function(df, df_info){
  # biomass data is in kg-DM/ha
  df |>
    left_join(
      df_info,
      by = join_by(lon_i, lat_i)
    ) |>
    # convert from kg-DM ha-2 on forest land to gC (integral across gridcell)
    mutate(db_gc_per_gridcell = db * 1e3 * 0.5 * fcf * area_ha) |>

    # sum across all gridcells globally and convert to PgC
    summarise(db_pgc = 1e-15 * sum(db_gc_per_gridcell, na.rm = TRUE)) |>
    pull(db_pgc)
}

calc_csink_bylat <- function(df, df_info){
  # biomass data is in kg-DM/ha
  df |>
    left_join(
      df_info,
      by = join_by(lon_i, lat_i)
    ) |>
    # convert from kg-DM ha-2 on forest land to gC (integral across gridcell)
    mutate(db_gc_per_gridcell = db * 1e3 * 0.5 * fcf * area_ha) |>

    # sum across all gridcells by latitudinal band and convert to TgC
    group_by(lat_i) |>
    summarise(db_tgc = 1e-12 * sum(db_gc_per_gridcell, na.rm = TRUE))
}


## Create bootstraps ------------
set.seed(1982)
n_boot <- 3000
boot_resamples <- bootstraps(data_forest_plots, times = n_boot)

### Un-parallel version --------------------------
# tic()
# df_boot <- boot_resamples |>
#   mutate(id = row_number()) |>
#   mutate(
#     bundle_stl = map(splits, ~ fit_stl_byboot(analysis(.x))),
#     model_biomass = map(splits, ~ fit_biomass_byboot(analysis(.x)))
#   ) |>
#   select(-splits) |>
#   mutate(
#     grid_predictions = map(
#       bundle_stl,
#       ~ predict_dn(., vec_qmd = data_forest_plots$QMD, df = grid_drivers)
#     )) |>
#   select(-bundle_stl) |>
#   mutate(
#     grid_predictions = map2(
#       grid_predictions,
#       model_biomass,
#       ~ predict_db(.x, .y)
#     )
#   ) |>
#   select(-model_biomass) |>
#
#   # calculate C sink by latitude (retain to keep sprad across bootstraps)
#   mutate(db_pgc_global = map_dbl(grid_predictions, ~calc_csink_global(., grid_drivers)))
#
# toc()

### Parallel version --------------------------
# ncores <- 4 # parallel::detectCores() - 2
#
# cl <- new_cluster(n = ncores) |>
#   cluster_library(
#     packages = c(
#       "dplyr",
#       "tidyr",
#       "lme4",
#       "purrr",
#       "recipes",
#       "rsample"
#     )
#   ) |>
#   cluster_assign(
#     fit_stl_byboot = fit_stl_byboot,
#     fit_biomass_byboot = fit_biomass_byboot,
#     predict_dn = predict_dn,
#     predict_db = predict_db,
#     data_forest_plots = data_forest_plots,
#     grid_drivers = grid_drivers,
#     calc_csink_global = calc_csink_global
#   )
#
# tic()
# df_boot <- boot_resamples |>
#   mutate(id = row_number()) |>
#   partition(cl) |>
#   mutate(
#     bundle_stl = map(splits, ~ fit_stl_byboot(analysis(.x))),
#     model_biomass = map(splits, ~ fit_biomass_byboot(analysis(.x)))
#   ) |>
#   select(-splits) |>
#   mutate(
#     grid_predictions = map(
#       bundle_stl,
#       ~ predict_dn(., vec_qmd = data_forest_plots$QMD, df = grid_drivers)
#     )) |>
#   select(-bundle_stl) |>
#   mutate(
#     grid_predictions = map2(
#       grid_predictions,
#       model_biomass,
#       ~ predict_db(.x, .y)
#     )
#   ) |>
#   select(-model_biomass) |>
#
#   # calculate C sink by latitude (retain to keep sprad across bootstraps)
#   mutate(db_pgc_global = map_dbl(grid_predictions, ~calc_csink_global(., grid_drivers))) |>
#
#   # drop full information to avoid excessive memory use - uncomment only when n_boot is very large
#   select(-grid_predictions) |>
#
#   collect()
# toc()
#
# write_rds(
#   df_boot,
#   file = here(paste0("data/df_boot_", lab_filter, "_nboot_", as.character(n_boot),".rds"))
# )

df_boot <- read_rds(here("data/df_boot_slopefilter_nboot_2000.rds"))

### Summarise across bootstraps ------------------------------------------------
# # stack predictions from all bootstrap samples into (very) long vector
# df_summ <- df_boot |>
#   select(-db_pgc_global) |>
#   unnest(grid_predictions) |>
#   group_by(lat_i, lon_i) |>
#   summarise(
#     db_mean = mean(db, na.rm = TRUE),
#     db_median = median(db, na.rm = TRUE),
#     db_sd = sd(db, na.rm = TRUE)
#   ) |>
#   drop_na() |>
#   ungroup() |>
#   left_join(
#     grid_drivers,
#     by = join_by(lon_i, lat_i)
#   ) |>
#   mutate(lon = lon_i/4, lat = lat_i/4) |>
#   mutate(
#     across(
#       starts_with("db"),
#       list(
#         gc_per_ha_forest   = ~ .x * 1e3 * 0.5,
#         gc_per_ha_gridcell = ~ (.x * 1e3 * 0.5) * fcf
#       ),
#       .names = "{.col}_{.fn}"
#     )
#   )
#
# write_rds(
#   df_summ,
#   here(paste0("data/df_summ_", lab_filter, "_nboot_", as.character(n_boot),".rds"))
# )

df_summ <- read_rds(here("data/df_summ_slopefilter_nboot_30.rds"))

## Visualisations --------------------------------------------------------------
### Distribution of global C sink estimates ------------------------------------
# sum across all gridcells, distribution across bootstraps
gg_hist_csink_boot <- df_boot |>
  ggplot(aes(db_pgc_global, after_stat(density))) +
  geom_histogram(fill = "grey", color = "black", bins = 50) +
  labs(
    x = expression(paste("PgC ", yr^-1)),
    y = "Density"
  ) +
  theme_bw()

ggsave(
  here("manuscript/figures/gg_hist_csink_boot.pdf"),
  plot = gg_hist_csink_boot,
  width = 5,
  height = 4
)

# numbers for paper
mean(df_boot$db_pgc_global)
quantile(df_boot$db_pgc_global, probs = c(0.025, 0.975))

# # across gridcells (median across bootstraps)
# gg_hist_db_boot <- df_summ |>
#   ggplot(aes(db_median_gc_per_ha_forest * 1e-6, after_stat(density))) +
#   geom_histogram(fill = "grey", color = "black", bins = 50) +
#   labs(
#     x = expression(paste("MgC ha"^-1, "yr"^-1)),
#     y = "Density"
#   ) +
#   theme_bw()
#
# ggsave(
#   here("manuscript/figures/gg_hist_db_boot.pdf"),
#   plot = gg_hist_db_boot,
#   width = 5,
#   height = 4
#   )

### Map ---------------------------------------------------------------
#### Per forest area -----------
##### Mean --------------
# df_summ <- read_rds(file = here("data/df_summ.rds"))

coast <- rnaturalearth::ne_coastline(
  scale = 110,
  returnclass = "sf"
)

layer_ocean <- rnaturalearth::ne_download( # ne_load(
  scale = 110,
  type = "ocean",
  category = "physical",
  returnclass = "sf",
  destdir = here("data/")
)

gg_map_sink_perforestarea <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_median_gc_per_ha_forest * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  # scale_fill_stepsn(
  #   colours = rev(khroma::color("berlin")(50)), # viridisLite::viridis(7),
  #   guide = guide_coloursteps(),
  #   na.value = "grey30", # <- missing data color
  #   limits = c(-1.5, 1.5),
  #   oob = squish, # clamp values outside limits,
  #   name = expression(paste("MgC ha"^-1, "yr"^-1)),
  # ) +
  khroma::scale_fill_berlin(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    limits = c(-1.5, 1.5),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1)),
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey40", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_perforestarea.pdf"),
  plot = gg_map_sink_perforestarea,
  width = 8,
  height = 6
)

##### SD --------------
gg_map_sink_perforestarea_sd <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_sd_gc_per_ha_forest * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  khroma::scale_fill_lajolla(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    limits = c(0, 0.2),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey40", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_perforestarea_sd.pdf"),
  plot = gg_map_sink_perforestarea_sd,
  width = 8,
  height = 4
)

gg_map_sink_perforestarea_combined <- cowplot::plot_grid(
  gg_map_sink_perforestarea,
  gg_map_sink_perforestarea_sd,
  ncol = 1,
  labels = letters[1:2]
)

ggsave(
  here("manuscript/figures/gg_map_sink_perforestarea_combined.pdf"),
  plot = gg_map_sink_perforestarea_combined,
  width = 8.2,
  height = 5.5
)

#### Per grid area -----------
##### Mean --------------
gg_map_sink_pergridarea <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_median_gc_per_ha_gridcell * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  khroma::scale_fill_berlin(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    limits = c(-1.5, 1.5),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey30", color = NA)) # ocean = light grey

gg_map_sink_pergridarea

ggsave(
  here("manuscript/figures/gg_map_sink_pergridarea.pdf"),
  plot = gg_map_sink_pergridarea,
  width = 8,
  height = 6
)

##### SD --------------
gg_map_sink_pergridgridtarea_sd <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_sd_gc_per_ha_gridcell * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  khroma::scale_fill_lajolla(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    limits = c(0, 0.15),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey40", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_pergridgridtarea_sd.pdf"),
  plot = gg_map_sink_pergridgridtarea_sd,
  width = 8,
  height = 4
)

gg_map_sink_pergridarea_combined <- cowplot::plot_grid(
  gg_map_sink_pergridarea,
  gg_map_sink_pergridgridtarea_sd,
  ncol = 1,
  labels = letters[1:2]
)

ggsave(
  here("manuscript/figures/gg_map_sink_pergridarea_combined.pdf"),
  plot = gg_map_sink_pergridarea_combined,
  width = 8.2,
  height = 5.5
)

### Sink across latitude ---------
df_boot_small <- read_rds(here("data/df_boot_slopefilter_nboot_30.rds"))

df_bylat <- df_boot_small |>
  select(-db_pgc_global) |>
  mutate(csink_tgc_bylat = map(grid_predictions, ~calc_csink_bylat(., grid_drivers))) |>
  select(-grid_predictions) |>
  unnest(csink_tgc_bylat) |>
  group_by(lat_i) |>
  summarise(
    db_pgc_bylat_median = median(db_tgc, na.rm = TRUE),
    db_pgc_bylat_q05 = quantile(db_tgc, probs = 0.5, na.rm = TRUE),
    db_pgc_bylat_q95 = quantile(db_tgc, probs = 0.95, na.rm = TRUE),
  ) |>
  mutate(lat = lat_i / 4)

gg_csink_bylat <- df_bylat |>
  ggplot() +   # multiplication by two because 0.5 degree resolution
  geom_ribbon(
    aes(
      x = lat,
      ymin = db_pgc_bylat_q05 * 2,
      ymax = db_pgc_bylat_q95 * 2),
    fill = "grey50") +
  geom_line(aes(lat, db_pgc_bylat_median * 2)) +
  theme_classic() +
  labs(
    x = "Latitude (°)",
    y = expression(paste("C sink (TgC yr"^{-1}, " °" ^{-1},")"))
  ) +
  coord_flip()

gg_csink_bylat

## Publication figure --------------------------
gg_map_sink_all <- cowplot::plot_grid(
  gg_csink_bylat,
  gg_map_sink_pergridarea,
  gg_hist_csink_boot,
  gg_map_sink_pergridgridtarea_sd,
  ncol = 2,
  labels = letters[1:4],
  rel_widths = c(0.3, 1)
)

ggsave(
  here("manuscript/figures/gg_map_sink_all.pdf"),
  plot = gg_map_sink_all,
  width = 11,
  height = 5.5
)

