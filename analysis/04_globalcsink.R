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
# Load and engineer data with environmental factors
# plot-level data for model fitting
data_forest_plots <- read_rds(here::here("data/inputs/data_fil75_biomes.rds")) |>
  # filter(year > 1980) |> # XXX why this filter?
  mutate(NQMD2 = density * QMD^2)

### Maps of environmental covariates -------------------------------------------
# Load data for upscaling: maps of environmental factors
grid_drivers <- read_rds(here::here("data/global_drivers.rds")) |>
  as_tibble() |>
  mutate(lon_i = round(lon * 4), lat_i = round(lat * 4))

### Forest cover fraction ------------------------------------------------------
# load modis fraction forest cover raster
# r_fcf <- terra::rast("/home/laura/data/forest_fraction/MODIS_ForestCoverFraction.nc")
r_fcf <- terra::rast(
  "~/data/archive/forestcovermodis_dimiceli_2015/data/MODIS-C006_MOD44B_ForestCoverFraction/MODIS-TERRA_C6__MOD44B__ForestCoverFraction__LPDAAC__GLOBAL__0.5degree__UHAM-ICDC__20100306__fv0.02.nc",
  lyrs = "forestcoverfraction")

df_fcf <- as.data.frame(r_fcf, xy = TRUE, na.rm = FALSE) |>
  as_tibble() |>
  rename(lon = x, lat = y) |>
  mutate(
    lon_i = round(lon * 4),
    lat_i = round(lat * 4),
    forestcoverfraction = forestcoverfraction * 1e-2
    )

# combine with environmental covariates df
grid_drivers <- grid_drivers |>
  left_join(
    df_fcf |>
      select(-lon, -lat),
    by = join_by(lon_i, lat_i)
  ) |>
  rename(fcf = forestcoverfraction)

## Function definitions --------------
### Fit LMM for STL per bootstrap --------
fit_stl_byboot <- function(df) {
  vars_to_scale <- c("logQMD", "year", "ai", "ndep", "ORGC", "PBR")

  rec <- recipe(
    logDensity ~ logQMD + year + ai + ndep + ORGC + PBR + dataset + plotID + species,
    data = df |>
      drop_na(
        all_of(
          c(
            "logDensity",
            "logQMD",
            "year",
            "ai",
            "ndep",
            "ORGC",
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
      year * ai +
      year * ndep +
      year * ORGC +
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

## Create bootstraps ------------
n_boot <- 300 # Will have to increase this
boot_resamples <- bootstraps(data_forest_plots, times = n_boot)

### Un-parallel version --------------------------
tic()
df_boot <- boot_resamples %>%
  slice(1:30) |>
  mutate(id = row_number()) |>
  mutate(
    bundle_stl = map(splits, ~ fit_stl_byboot(analysis(.x))),
    model_biomass = map(splits, ~ fit_biomass_byboot(analysis(.x)))
  ) |>
  mutate(
    grid_predictions = map(
      bundle_stl,
      ~ predict_dn(., vec_qmd = data_forest_plots$QMD, df = grid_drivers)
    ),
    grid_predictions = map2(
      grid_predictions,
      model_biomass,
      ~ predict_db(.x, .y)
    )
  )
toc()

### Parallel version --------------------------
ncores <- 4 # parallel::detectCores() - 2

cl <- new_cluster(n = ncores) |>
  cluster_library(
    packages = c(
      "dplyr",
      "tidyr",
      "lme4",
      "purrr",
      "recipes",
      "rsample"
    )
  ) |>
  cluster_assign(
    fit_stl_byboot = fit_stl_byboot,
    fit_biomass_byboot = fit_biomass_byboot,
    predict_dn = predict_dn,
    predict_db = predict_db,
    data_forest_plots = data_forest_plots,
    grid_drivers = grid_drivers
  )

tic()
df_boot <- boot_resamples %>%
  mutate(id = row_number()) |>
  partition(cl) |>
  mutate(
    bundle_stl = map(splits, ~ fit_stl_byboot(analysis(.x))),
    model_biomass = map(splits, ~ fit_biomass_byboot(analysis(.x)))
  ) |>
  mutate(
    grid_predictions = map(
      bundle_stl,
      ~ predict_dn(., vec_qmd = data_forest_plots$QMD, df = grid_drivers)
    ),
    grid_predictions = map2(
      grid_predictions,
      model_biomass,
      ~ predict_db(.x, .y)
    )
  ) |>
  # mean across gridcells within bootstraps
  mutate(
    dn_mean = map_dbl(grid_predictions, ~ mean(.$dn, na.rm = TRUE)),
    db_mean = map_dbl(grid_predictions, ~ mean(.$db, na.rm = TRUE))
  ) |>
  collect()
toc()

## Global C sink calculation ---------------------------------------------------
calc_global_csink <- function(df, df_info){
  # biomass data is in kg-DM/ha
  df |>
    left_join(
      df_info,
      by = join_by(lon_i, lat_i)
    ) |>
    # convert from kg-DM ha-2 on forest land to gC (integral across gridcell)
    mutate(db_gc_per_gridcell = db * 1e3 * 0.5 * fcf * area_ha) |>
    summarise(db_gc_global = sum(db_gc_per_gridcell, na.rm = TRUE)) |>
    pull(db_gc_global)
}

# calculate global C sink in PgC on each bootstrap
df_boot <- df_boot |>
  mutate(db_pgc_global = 1e-15 * map_dbl(grid_predictions, ~calc_global_csink(., grid_drivers)))

write_rds(df_boot, file = here("data/df_boot.rds"))

### Summarise across bootstraps ------------------------------------------------
# stack predictions from all bootstrap samples into (very) long vector
df_summ <- df_boot |>
  select(id_boot = id, grid_predictions) |>
  unnest(grid_predictions) |>
  group_by(lon_i, lat_i) |>
  summarise(
    db_mean = mean(db, na.rm = TRUE),
    db_median = median(db, na.rm = TRUE),
    db_sd = sd(db, na.rm = TRUE)
  ) |>
  drop_na() |>
  ungroup() |>
  left_join(
    grid_drivers,
    by = join_by(lon_i, lat_i)
  ) |>
  mutate(lon = lon_i/4, lat = lat_i/4) |>
  mutate(
    across(
      starts_with("db"),
      list(
        gc_per_ha_forest   = ~ .x * 1e3 * 0.5,
        gc_per_ha_gridcell = ~ (.x * 1e3 * 0.5) * fcf
      ),
      .names = "{.col}_{.fn}"
    )
  )

write_rds(df_summ, file = here("data/df_summ.rds"))

## Visualisations --------------------------------------------------------------
### Distribution of global C sink estimates ------------------------------------
# sum across all gridcells, distribution across bootstraps
gg_hist_csink_boot <- df_boot |>
  ggplot(aes(db_pgc_global, ..density..)) +
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

# across gridcells (median across bootstraps)
gg_hist_db_boot <- df_summ |>
  ggplot(aes(db_median_gc_per_ha_forest * 1e-6, ..density..)) +
  geom_histogram(fill = "grey", color = "black", bins = 50) +
  labs(
    x = expression(paste("MgC ha"^-1, "yr"^-1)),
    y = "Density"
  ) +
  theme_bw()

ggsave(
  here("manuscript/figures/gg_hist_db_boot.pdf"),
  plot = gg_hist_db_boot,
  width = 5,
  height = 4
  )

### Map ---------------------------------------------------------------
#### Per forest area -----------
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
  khroma::scale_fill_berlin(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey20", # <- missing data color
    limits = c(-1.5, 1.5),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey20", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_perforestarea.pdf"),
  plot = gg_map_sink_perforestarea,
  width = 8,
  height = 6
  )

#### Per grid area -----------
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
    na.value = "grey20", # <- missing data color
    limits = c(-1, 1),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey20", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_pergridarea.pdf"),
  plot = gg_map_sink_pergridarea,
  width = 8,
  height = 6
)


