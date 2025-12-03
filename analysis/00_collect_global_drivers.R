# This creates the data frame for all grids globally (0.5 deg) containing data
# on all environmental predictors determined in 03_env_drivers.R and used in
# 04_globalcsink.R for global upscaling.

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tidyterra)
library(here)
library(ingestr)

# Forest cover fraction --------------------------------------------------------
r_fcf <- terra::rast(
  "~/data/archive/forestcovermodis_dimiceli_2015/data/MODIS-C006_MOD44B_ForestCoverFraction/MODIS-TERRA_C6__MOD44B__ForestCoverFraction__LPDAAC__GLOBAL__0.5degree__UHAM-ICDC__20100306__fv0.02.nc",
  lyrs = "forestcoverfraction"
)

df_fcf <- as.data.frame(r_fcf, xy = TRUE, na.rm = FALSE) |>
  as_tibble() |>
  rename(lon = x, lat = y, fcf = forestcoverfraction) |>
  mutate(
    lon_i = round(lon * 4),
    lat_i = round(lat * 4),
    fcf = fcf * 1e-2 # to make it a fraction (not percent)
  )

# Forest biomes ----------------------------------------------------------------
# Rasterise biomes map and select forest biomes for determining the maximum
# extent of global upscaling.
# sf::sf_use_s2(FALSE)

# Biomes -------------
# Get data from Olson, D. M, 2020, "Terrestrial ecoregions of the world",
# https://doi.org/10.7910/DVN/WTLNRG, Harvard Dataverse, V1
# WWF Ecoregions data
# Load shapefile
biomes <- sf::read_sf("~/data/biomes/olson_harvard_dataverse/tnc_terr_ecoregions.shp")

# get biome names and code
df_biomes_codes <- biomes |>
  select(WWF_MHTNAM, WWF_MHTNUM) |>
  as.data.frame() |>
  select(-geometry) |>
  distinct() |>
  arrange(WWF_MHTNUM)

# write biome codes to file (for later use)
write_csv(df_biomes_codes, file = here("data/df_biomes_codes.csv"))

# subset only forest biomes
df_forest_biomes <- df_biomes_codes |>
  filter(WWF_MHTNUM %in% c(1, 2, 4, 5, 6, 12))

# Rasterise
# perform shapefile to raster conversion, use forest cover fraction raster
# as template
biomes_raster <- terra::rasterize(biomes, r_fcf, field = "WWF_MHTNAM")

# Mask to keep only forest biomes in raster
mask <- biomes_raster %in% df_forest_biomes$WWF_MHTNAM
biomes_raster_forest <- mask(biomes_raster, mask, maskvalues = 0)

# # write to file as NetCDF
# terra::writeCDF(
#   biomes_raster,
#   filename = here("data/biomes_raster_0.5deg.nc"),
#   overwrite = TRUE,
#   varname = "biome"
# )

df_grid <- as.data.frame(biomes_raster_forest, xy = TRUE, na.rm = FALSE) |>
  as_tibble() |>
  rename(lon = x, lat = y) |>
  left_join(
    df_forest_biomes,
    by = join_by(WWF_MHTNAM)
  ) |>
  rename(biome_name = WWF_MHTNAM, biome_number = WWF_MHTNUM) |>
  drop_na()

# Environmental drivers --------------------------------------------------------
## MAT -------------------------------------------------------------------------
# Read monthly mean temperature into one raster stack
mat_stack <- rast(
  list.files(
    "~/data/archive/worldclim_fick_2017/data_10m/wc2.1_10m_tavg/",
    pattern = "\\.tif$",
    full.names = TRUE
    )
  )

# Aggregate across layers using mean
mat_mean <- mean(mat_stack)

# Regrid to template (forest cover fraction, half degree)
mat_mean_halfdeg <- resample(mat_mean, r_fcf, method = "bilinear")

# extract values into data frame
df_tavg <- as.data.frame(mat_mean_halfdeg, xy = TRUE, na.rm = FALSE) |>
  as_tibble() |>
  rename(lon = x, lat = y, tavg = mean) |>
  drop_na()

## N deposition ----------------------------------------------------------------
rasta_ndep <- rast("~/data/ndep_lamarque/Ndep_Lamarque11cc_historical_halfdeg.nc")

# mask to retain only forest biome pixels (mask defined above)
rasta_ndep_mask <- mask(rasta_ndep, mask, maskvalues = 0)

# take sum of both species
noy <- rasta_ndep_mask[[1:160]]
nhx <- rasta_ndep_mask[[161:320]]
rasta_ndep <- noy + nhx

# subset to years 1970-2009
years <- 1850:2009
target_years <- 1970:2009
idx <- which(years %in% target_years)
rasta_ndep <- rasta_ndep[[idx]]

# Mean N deposition 1970–2009
rasta_ndep_mean <- mean(rasta_ndep)

# Linear trend (slope) per pixel, 1970–2009
# time vector (years)
t <- target_years

trend_fun <- function(y) {
  # y is the vector of values for a single pixel across layers
  if (all(is.na(y))) return(NA)
  coef(lm(y ~ t))[2]   # slope per year
}

ndep_trend <- app(rasta_ndep, trend_fun)

# extract mean N-dep values into data frame
df_ndep <- as.data.frame(rasta_ndep_mean, xy = TRUE, na.rm = FALSE) |>
  as_tibble() |>
  rename(lon = x, lat = y, ndep = mean) |>
  drop_na()

## Moisture index --------------------------------------------------------------
rasta_mi <- rast("~/data/aridityindex_zomer_2022/Global-AI_ET0_v3_annual/aridityindex_p_over_pet_zomeretal2022_v3_yr_halfdeg.nc")

# mask to retain only forest biome pixels (mask defined above)
rasta_mi_mask <- mask(rasta_mi, mask, maskvalues = 0)

# extract mean N-dep values into data frame
df_mi <- as.data.frame(rasta_mi_mask, xy = TRUE, na.rm = FALSE) |>
  as_tibble() |>
  rename(lon = x, lat = y) |>
  drop_na()

## Soil C:N --------------------------------------------------------------------
siteinfo <- df_grid |>
  mutate(sitename = 1:n()) |>
  select(sitename, lon, lat)

settings_wise <- get_settings_wise(varnam = c("CNrt", "ORGC"), layer = 1:3)

df_wise <- ingest(
  siteinfo,
  source    = "wise",
  settings  = settings_wise,
  dir       = "~/data/soil/wise") |>
  unnest(cols = data) |>
  right_join(
    df_grid |>
      mutate(sitename = 1:n()) |>
      select(sitename, lon, lat),
    by = "sitename"
  ) |>
  ungroup() |>
  drop_na() |>
  select(-sitename)

# Clean and combine data frames ------------------------------------------------
clean_latlon_halfdeg <- function(df){
  df |>
    mutate(
      lon = round(lon*4)/4,
      lat = round(lat*4)/4,
      lon_i = as.integer(round(lon*4)),
      lat_i = as.integer(round(lat*4))
      )
}

df_grid <- df_grid |>
  clean_latlon_halfdeg() |>

  # add forest cover fraction
  left_join(
    df_fcf |>
      clean_latlon_halfdeg() |>
      select(-lat, -lon),
    by = join_by(lat_i, lon_i)
  ) |>

  # add MAT
  left_join(
    df_tavg |>
      clean_latlon_halfdeg() |>
      select(-lat, -lon),
    by = join_by(lat_i, lon_i)
  ) |>

  # add Ndep
  left_join(
    df_ndep |>
      clean_latlon_halfdeg() |>
      select(-lat, -lon),
    by = join_by(lat_i, lon_i)
  ) |>

  # add MI
  left_join(
    df_mi |>
      clean_latlon_halfdeg() |>
      select(-lat, -lon),
    by = join_by(lat_i, lon_i)
  ) |>

  # add CN and ORGC
  left_join(
    df_wise |>
      clean_latlon_halfdeg() |>
      select(-lat, -lon),
    by = join_by(lat_i, lon_i)
  )

# write to file
write_rds(
  df_grid,
  file = here("data/df_grid.rds")
)

# df_grid |>
#   ggplot(aes(lon, lat, fill = CNrt)) +
#   geom_raster() +
#   coord_sf() +
#   scale_fill_viridis_c()

# Map sites --------------------------------------------------------------------
# get coastline
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

df_plots <- read_rds(here("data/inputs/data_fil75_biomes.rds")) |>
  distinct(plotID, lat, lon)

# df_plots |>
#   filter(dataset %in% "lwf")
#
# data_unm_euforia <- df_plots |>
#   filter(dataset %in% c(
#     "bnp", "czu", "forst", "iberbas", "incds", "lwf", "nbw", "uholka",
#     "nfr", "nwfva", "tuzvo", "ul", "unito", "urk", "wuls"
#   )) |>
#   distinct(plotID, dataset, lat, lon)

# # Keep only the first occurrence of each name
# df_labels <- data_unm_euforia %>%
#   group_by(dataset) %>%
#   slice(1) %>%
#   ungroup()

## Site density ----------------------------------------------------------------
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")

gg_sitedensity <- ggplot() +
  # world country outlines
  geom_sf(data = world, fill = "gray95", color = "gray70", size = 0.2) +

  # hex bin layer: count of points per hex
  stat_bin_hex(
    data = df_plots,
    mapping = aes(x = lon, y = lat),
    bins = 60,
    color = NA,
    alpha = 0.9
  ) +

  # color (count) scale
  scale_fill_viridis_c(
    name = "Sites\ncount",
    option = "D",
    trans = "sqrt",
    na.value = "transparent"
  ) +

  # coordinate system (preserves lat/lon aspect)
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE # to draw map strictly bounded by the specified extent
  ) +

  # labels and theme
  labs(
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "gray90", size = 0.2),
    legend.position = "right"
  )

ggsave(
  here("manuscript/figures/map_sitedensity.pdf"),
  plot = gg_sitedensity,
  width = 8,
  height = 4
  )

## Site points -----------------------------------------------------------------
gg_sitepoints <- ggplot() +
  tidyterra::geom_spatraster(data = biomes_raster_forest) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.2,
    fill = "grey70"
  ) +
  scale_fill_manual(
    values = c(
      "Boreal Forests/Taiga"                                         = "dodgerblue4",
      # "Tundra"                                                       = "lightcyan3"
      # "Deserts and Xeric Shrublands"                                 = "#FFD3A0",
      # "Flooded Grasslands and Savannas"                              = "indianred3",
      # "Inland Water"                                                 = "azure",
      # "Mangroves"                                                    = "violetred",
      "Mediterranean Forests, Woodlands and Scrub"                   = "orangered3",
      # "Montane Grasslands and Shrublands"                            = "steelblue3",
      # "Rock and Ice"                                                 = "azure4",
      "Temperate Broadleaf and Mixed Forests"                        = "darkgreen",
      "Temperate Conifer Forests"                                    = "lightseagreen",
      # "Temperate Grasslands, Savannas and Shrublands"                = "goldenrod3",
      # "Tropical and Subtropical Coniferous Forests"                  = "#31A278",
      "Tropical and Subtropical Dry Broadleaf Forests"               = "goldenrod4",
      # "Tropical and Subtropical Grasslands, Savannas and Shrublands" = "darkolivegreen",
      "Tropical and Subtropical Moist Broadleaf Forests"             = "springgreen3"
    ),
    na.value = "transparent",
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  geom_point(
    aes(lon, lat),
    data = df_plots,
    color = "red",
    fill = "white",
    alpha = .7,
    shape = 21,
    size = 1.2
    ) +
  # geom_point(
  #   aes(lon, lat),
  #   data = data_unm_euforia,
  #   color = "blue",
  #   fill = "white",
  #   alpha = .7,
  #   shape = 21,
  #   size = 1.2
  #   ) +
  # geom_count(aes(lon, lat), data = data_biomes_fil_unique_plots, color="red",alpha = .6) +
  # scale_size_area() +
  # geom_hex(aes(lon, lat, color = biome), data = data_biomes_fil, bins = 50, linewidth = 1) +
  # set extent in longitude and latitude
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE # to draw map strictly bounded by the specified extent
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid = element_line(color = "gray90", size = 0.2),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    # panel.background = element_rect(fill = "white")
  )

gg_sitepoints

ggsave(
  here("manuscript/figures/gg_sitepoints.pdf"),
  plot = gg_sitepoints,
  width = 10,
  height = 5
)

## Plot plots in biomes --------------------------------------------------------
# remotes::install_github("https://github.com/valentinitnelav/plotbiomes")
# library(plotbiomes)
