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

## Load data -----------
# Load and engineer data with environmental factors
# plot-level data for model fitting
data_forest_plots <- read_rds(here::here("data/inputs/data_fil75_biomes.rds")) |>
  # filter(year > 1980) |> # XXX why this filter?
  mutate(NQMD2 = density * QMD^2)

# Load data for upscaling: maps of environmental factors
grid_drivers <- read_rds(here::here("data/global_drivers.rds")) |>
  as_tibble()

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
    select(lon, lat, area_ha, dn, dnqmd2)

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
n_boot <- 500 # Will have to increase this
boot_resamples <- bootstraps(data_forest_plots, times = n_boot)

### Single-core version --------------------------
# Bootstrap STL and biomass model fitting
tic()
df_boot <- boot_resamples %>%
  slice(1:3) |>
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
  mutate(
    dn_mean = map_dbl(grid_predictions, ~ mean(.$dn, na.rm = TRUE)),
    db_mean = map_dbl(grid_predictions, ~ mean(.$db, na.rm = TRUE))
  )
toc()

#### Plot ---------------------
# distribution of mean changes across bootstrap samples
hist(df_boot$grid_predictions[[1]]$dn)
hist(df_boot$grid_predictions[[1]]$db)

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
df_boot_parallel <- boot_resamples %>%
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
  mutate(
    dn_mean = map_dbl(grid_predictions, ~ mean(.$dn, na.rm = TRUE)),
    db_mean = map_dbl(grid_predictions, ~ mean(.$db, na.rm = TRUE))
  ) |>
  collect()
toc()

# unit conversions
df_boot_parallel <- df_boot_parallel |>
  mutate(grid_predictions = map(
    grid_predictions,
    ~ mutate(
      .,
      dB_Mg_ha = db * 10^-3,
      dC_Mg_ha = dB_Mg_ha * 0.5
    )
  ))

#### Plot ---------------------
# distribution of mean changes across bootstrap samples
hist(df_boot_parallel$grid_predictions[[1]]$dn)
hist(df_boot_parallel$grid_predictions[[1]]$db)

write_rds(df_boot_parallel, file = here("data/df_boot_parallel.rds"))

### Summarise across bootstraps ------------------------------------------------
# stack predictions from all bootstrap samples into (very) long vector
df_boot_unnested <- df_boot_parallel |>
  select(id, grid_predictions) |>
  unnest(grid_predictions)

write_rds(df_boot_unnested, file = here("data/df_boot_unnested.rds"))

df_summ <- df_boot_unnested |>
  mutate(lon_i = round(lon * 4), lat_i = round(lat * 4)) |>
  group_by(lon_i, lat_i) |>
  summarise(
    db_mean = mean(db, na.rm = TRUE),
    db_median = median(db, na.rm = TRUE),
    db_sd = sd(db, na.rm = TRUE)
  ) |>
  ungroup() |>
  mutate(lon = lon_i/4, lat = lat_i/4) |>
  drop_na()

write_rds(df_summ, file = here("data/df_summ.rds"))

## Visualisations --------------------------------------------------------------
### Distribution ---------------------------------------------------------------
df_boot_unnested <- read_rds(file = here("data/df_boot_unnested.rds"))

df_boot_unnested |>
  ggplot(aes(dC_Mg_ha, ..density..)) +
  geom_histogram(fill = "grey", color = "black", bins = 50) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  # theme_classic() +
  labs(x = expression(paste("Mg C ", ha^-1, " ", yr^-1))) +
  theme_bw() +
  xlim(-3, 3)

ggsave(here("manuscript/figures/histogram_db.pdf"))

### Map ---------------------------------------------------------------
df_summ <- read_rds(file = here("data/df_summ.rds"))

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

df_summ |>
  mutate(
    dB_Mg_ha = db_median * 10^-3,
    dC_Mg_ha = dB_Mg_ha * 0.5
  ) |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = dC_Mg_ha),
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
    oob = squish # clamp values outside limits
  ) +
  # scale_fill_viridis_c(
  #   # name =  expression(paste("Mg C ", ha^-1, " ", yr^-1))
  # ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey20", color = NA)) # ocean = light grey
# labs(
#   #subtitle = "Global forest carbon increase per ha and year"
# )

ggsave(here("manuscript/figures/fig4.pdf"))

## Global C sink calculation ---------------------------------------------------
### Forest cover fraction weighing ---------------------------------------------

# load modis fraction forest cover raster
# r_fcf <- terra::rast("/home/laura/data/forest_fraction/MODIS_ForestCoverFraction.nc")
r_fcf <- terra::rast(
  "/data/archive/forestcovermodis_dimiceli_2015/data/MODIS-C006_MOD44B_ForestCoverFraction/MODIS-TERRA_C6__MOD44B__ForestCoverFraction__LPDAAC__GLOBAL__0.5degree__UHAM-ICDC__20100306__fv0.02.nc",
  lyrs = "forestcoverfraction")

df_fcf <- as.data.frame(r_fcf, xy = TRUE, na.rm = FALSE)

plot(r_fcf)

# Convert df to SpatVector
points <- vect(results, geom = c("lon", "lat"), crs = crs(r_fcf))

# Extract raster values at given points
extracted <- extract(r_fcf, points, fun = NULL, na.rm = TRUE, touches = TRUE)

# Combine extracted values with polygon attributes
results$fcf <- extracted$forestcoverfraction  # Assuming 'layer' contains the raster values

db_Pg_yr = dB_Mg_ha*1e-9*area_ha*fcf*1e-2

