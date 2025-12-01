library(dplyr)
library(lme4)
library(MASS)  # for mvrnorm
library(sf)
library(purrr)
library(tidyr)
library(here)
library(lmerTest) 
library(ggplot2)
library(ggeffects)
library(effects)
library(viridis)
library(tidymodels)
library(terra)
library(ncdf4)
library(patchwork)
#library(ggokabeito)
# library(ingestr)
library(multidplyr)
library(tictoc)
library(readr)

# Prepare data for global C sink estimate ---------
# This script evaluates the mature forests C sink

# Load functions ----
source(here::here("R/functions.R"))
source(here::here("R/get_drivers_by_biome.R"))
source(here("R/csink_global.R"))
source(here("R/calc_db.R"))

# Collect global environmental covariates --------------
filn <- here::here("data/global_drivers.rds")

if (!file.exists(filn)){
  
  # From the WWF Ecoregions data
  # to read as data.frame
  # v_biomes <- st_read(file.path(here::here(), "/data/wwf/wwf_terr_ecos.shp")) 
  v_biomes <- st_read("/data/scratch/bstocker/biomes/wwf_ecoregions/official/wwf_terr_ecos.shp")  # XXX deposit on data_archive
  
  # Filter only the forest biomes, BIOME == 1,2,3,4,5,6,11,12 XXX where is this documented which indexes are forest?
  v_biomes_forests <- v_biomes |>
    dplyr::filter(BIOME==1|BIOME==2|BIOME==3|BIOME==4|BIOME==5|BIOME==6|BIOME==12)
  
  global_drivers <- get_drivers_by_biome(v_biomes_forests)
  
  # Save stand-level data
  saveRDS(global_drivers, file = file.path(here::here(), "/data/inputs/global_drivers.rds"))
  
} else {
  
  global_drivers <- readRDS(here::here("data/global_drivers.rds"))
  
}

global_drivers <- global_drivers |> 
  as_tibble()

## Visualisation of drivers --------
coast <- rnaturalearth::ne_coastline(
  scale = 110,
  returnclass = "sf"
)

### AI ---------------
global_drivers |> 
  ggplot() +
  geom_raster(
    aes(x = lon, y = lat, fill = ai),
    show.legend = TRUE
  ) +
  geom_sf(
    data = coast,
    colour = 'black',
    linewidth = 0.3
  )  +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  scale_fill_viridis_c(
    option = "magma"
  ) +
  theme_void()

### PBR -----------------
global_drivers |> 
  ggplot() +
  geom_raster(
    aes(x = lon, y = lat, fill = PBR),
    show.legend = TRUE
  ) +
  geom_sf(
    data = coast,
    colour = 'black',
    linewidth = 0.3
  )  +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  scale_fill_viridis_c(
    option = "cividis"
  ) +
  theme_void()

### ORGC -----------------
global_drivers |> 
  ggplot() +
  geom_raster(
    aes(x = lon, y = lat, fill = ORGC),
    show.legend = TRUE
  ) +
  geom_sf(
    data = coast,
    colour = 'black',
    linewidth = 0.3
  )  +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  scale_fill_viridis_c(
    option = "cividis"
  ) +
  theme_void()

### N-deposition  ---------------------
global_drivers |> 
  ggplot() +
  geom_raster(
    aes(x = lon, y = lat, fill = ndep),
    show.legend = TRUE
  ) +
  geom_sf(
    data = coast,
    colour = 'black',
    linewidth = 0.3
  )  +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  scale_fill_viridis_c(
    option = "viridis"
  ) +
  theme_void()


# Global C sink estimate --------------
# Number of samples of model coefficients (global, fully correlated across gridcells)
n_coef <- 50

# Number of samples of QMD within each gridcell (gridcell-specific, uncorrelated across gridcells)
n_qmd <- 3

## Load data ----------------
# plot-level data for model fitting
data_forest_plots <- readRDS(here::here("data/data_fil_biomes.rds")) |>
  filter(year > 1980) |>    # XXX why this filter?
  mutate(NQMD2 = density * QMD^2)

data_forest_plots_selfthinning <- data_forest_plots |> 
  drop_na(all_of(c("logDensity", "logQMD", "year", "ai", "ndep", "ORGC", "PBR", "dataset", "plotID", "species")))

data_forest_plots_biomass <- data_forest_plots |> 
  drop_na(all_of(c("biomass", "NQMD2", "dataset", "plotID")))

# scale (normalise) newdata with the parameters from the data used for model fitting
# repeat the single row for matiching dimensions of tmp_t0 and tmp_t1
data_forest_plots_selfthinning_means <- data_forest_plots_selfthinning |> 
  summarise(across(all_of(c("logQMD", "year", "ai", "ndep", "ORGC", "PBR")), mean))
# slice(rep(1:n(), times = nrow(tmp_t0)))

data_forest_plots_selfthinning_sds <- data_forest_plots_selfthinning |> 
  summarise(across(all_of(c("logQMD", "year", "ai", "ndep", "ORGC", "PBR")), sd))
# slice(rep(1:n(), times = nrow(tmp_t0)))

# global fields of predictors for upscaling
global_drivers <- readRDS(here::here("data/global_drivers.rds"))

## Load fitted self-thinning model ----------------
# produced with analysis/03_env_drivers.R
fit_selfthinning <- readRDS(file = here::here("data/mod_lmer_env.rds"))

# fit_selfthinning = lmer(
#   logDensity ~ scale(logQMD) + 
#     scale(year) * scale(ai) + 
#     scale(year) * scale(ndep) + 
#     scale(year) * scale(ORGC) + 
#     scale(year) * scale(PBR) + 
#     (1|dataset/plotID) + (1|species),  
#   data = data_forest_plots_selfthinning
#   )

## Fit relationship of biomass ~ N*QMD^2 ---------------------
fit_biomass = lmer(
  biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
  data = data_forest_plots_biomass
)

# not a lot of data for all biomes!
data_forest_plots_biomass |> 
  dplyr::group_by(biomeID, biome, dataset) |> 
  drop_na(biomass) |> 
  summarise(n = n())

# ### Alternative fit with biome as a grouping variable -----------
# shouldn't this be biome-specific? Problem: then we cannot predict 
# for gridcells that belong to biome for which we had no data for model fitting.
# fit_biomass2 = lmer(
#   biomass ~ NQMD2 + 0 + (1|biomeID) + (1|dataset/plotID), 
#   data = data_forest_plots_biomass
# )

### Inspect distribution  ---------------
# overall
data_forest_plots_biomass |> 
  ggplot(aes(QMD, color = biome, fill = biome)) +
  geom_density(alpha = 0.5) +
  khroma::scale_fill_okabeito() +
  khroma::scale_color_okabeito() +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )

### Inspect relationship biomass N*QMD2 ---------------
# overall
data_forest_plots_biomass |> 
  ggplot(aes(NQMD2, biomass)) +
  geom_hex(bins = 50, show.legend = FALSE) +
  geom_smooth(method = "lm", ) +
  khroma::scale_fill_batlowW(trans = "log", reverse = TRUE) +
  theme_classic() +
  geom_abline(intercept = 0, slope = coef(fit_biomass)$dataset$NQMD2[1], linetype = "dotted") +
  coord_fixed()

# by biome
data_forest_plots_biomass |> 
  ggplot(aes(NQMD2, biomass)) +
  geom_hex(bins = 50, show.legend = FALSE) +
  khroma::scale_fill_batlowW(trans = "log", reverse = TRUE) +
  theme_classic() +
  geom_abline(intercept = 0, slope = coef(fit_biomass)$dataset$NQMD2[1], linetype = "dotted") +
  coord_fixed() +
  facet_wrap(~biome, ncol = 2)

# by just moist broadleaved tropical
# XXX correct? 
data_forest_plots_biomass |> 
  filter(biomeID == 1) |> 
  ggplot(aes(NQMD2, biomass, color = dataset)) +
  geom_point() +
  theme_classic() +
  coord_fixed()

## Construct sampling ----------------
### Sample coefficients -----------------
fixef_means_selfthinning <- fixef(fit_selfthinning)
vcov_matrix_selfthinning <- vcov(fit_selfthinning)  # variance-covariance matrix of fixed effects

fixef_means_biomass <- fixef(fit_biomass)
vcov_matrix_biomass <- vcov(fit_biomass)  # variance-covariance matrix of fixed effects

# Sample from multivariate normal distribution
coef_samples_selfthinning <- MASS::mvrnorm(
  n = n_coef, 
  mu = fixef_means_selfthinning, 
  Sigma = vcov_matrix_selfthinning
  )

coef_samples_biomass <- MASS::mvrnorm(
  n = n_coef, 
  mu = fixef_means_biomass, 
  Sigma = vcov_matrix_biomass
)

# Combine to single data frame assuming the two are independent
coef_samples <- coef_samples_selfthinning |> 
  bind_cols(
    coef_samples_biomass
  )

### Sample QMD within each gridcell separately --------------
generate_samples_by_gridcell <- function(idx, global_drivers, data_forest_plots, n_qmd){
  global_drivers |> 
    slice(idx) |> 
    slice(rep(1, n_qmd)) |> 
    mutate(idx_qmd = row_number()) |> 
    mutate(logQMD = sample(data_forest_plots$logQMD, size = n_qmd))
}

# slow!
df_samples_qmd <- purrr::map_dfr(
  as.list(1:nrow(global_drivers)),  # doing it only for three gridcells
  ~generate_samples_by_gridcell(., global_drivers, data_forest_plots, n_qmd) 
)

### Cross the two data frames ---------
# For each sample of the coefficients
add_coefs <- function(idx, df_samples_qmd, coef_samples){
  coef_samples_row <- coef_samples |> 
    slice(idx)
  
  # Repeat B to match rows in A
  coef_samples_row_expanded <- coef_samples_row[rep(1, nrow(df_samples_qmd)), ]
  
  # Combine A and repeated B
  bind_cols(df_samples_qmd, coef_samples_row_expanded)    
}

df_samples <- purrr::map_dfr(
  as.list(1:n_coef),
  ~add_coefs(., df_samples_qmd, coef_samples)
)


# Calculate biomass difference for each sample ---------------------------------
filnam <- here::here(paste0("data/df_db_n_coef_", n_coef, "_n_qmd_", n_qmd, ".rds"))

if (!file.exists(filnam)){
  # ## Unparallel version ------------
  # tic()
  # df_db <- pmap(slice(df_samples, 1:30), function(...) {
  #   row_df <- tibble(...)  # reconstruct one-row data frame
  #   calc_db(
  #     row_df,
  #     data_forest_plots_selfthinning_means,
  #     data_forest_plots_selfthinning_sds,
  #     coef_samples_selfthinning,
  #     coef_samples_biomass
  #     )
  # }) |>
  #   unlist() |>
  #   as_tibble() |>
  #   setNames("dB") |>
  #   mutate(dB_Mg_ha = dB * 10^-3) |>
  #   bind_cols(
  #     df_samples |>
  #       slice(1:30) |>
  #       dplyr::select(plotID, lon, lat, area_ha)
  #   ) |>
  #   arrange(plotID)
  # toc()
  
  ## Parallelised version -------------
  ncores <- 25 # parallel::detectCores() - 2
  
  cl <- new_cluster(n = ncores) |> 
    cluster_library(packages = c("dplyr")) |> 
    cluster_assign(
      calc_db = calc_db,
      data_forest_plots_selfthinning_means = data_forest_plots_selfthinning_means,
      data_forest_plots_selfthinning_sds = data_forest_plots_selfthinning_sds,
      coef_samples_selfthinning = coef_samples_selfthinning,
      coef_samples_biomass = coef_samples_biomass
    )
  
  tic()
  df_db <- df_samples |>
    mutate(id = row_number()) |>
    partition(cl) |> 
    mutate(result = purrr::pmap_dbl(
      across(), 
      function(...) calc_db(
        tibble(...),
        data_forest_plots_selfthinning_means,
        data_forest_plots_selfthinning_sds,
        coef_samples_selfthinning,
        coef_samples_biomass
      )
    )) |> 
    collect() |> 
    dplyr::select(plotID, lon, lat, area_ha, dB = result) |> 
    mutate(dB_Mg_ha = dB * 10^-3) |> 
    arrange(plotID)
  toc()
  
  readr::write_rds(df_db, file = filnam)
  
} else {
  df_db <- read_rds(filnam)
}

## Visualisations ---------
df_db |> 
  ggplot(aes(dB_Mg_ha)) +
  geom_density() +
  theme_classic()
  
### map ---------------
df_db |> 
    ggplot() +
    geom_raster(
      aes(x = lon, y = lat, fill = dB_Mg_ha),
      show.legend = TRUE
    ) +
    geom_sf(
      data = coast,
      colour = 'black',
      linewidth = 0.3
    )  +
    coord_sf(
      ylim = c(-60, 85),
      expand = FALSE
    ) +
    scale_fill_viridis_c(
      option = "magma"
    ) +
    theme_void()

## Weighing with forest cover fraction

# XXX please complement, adopting code from 05_globalsink.R
