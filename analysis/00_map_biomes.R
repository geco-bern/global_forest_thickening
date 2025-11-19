# taken from analysis/plot_map_biomes.R

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyterra)
library(see)

#sf::sf_use_s2(FALSE)

# Biomes shapefile -------------
# Get data from Olson, D. M, 2020, "Terrestrial ecoregions of the world (Copy to use in GapAnalysis R package)", 
# https://doi.org/10.7910/DVN/WTLNRG, Harvard Dataverse, V1
# WWF Ecoregions data
biomes <- sf::read_sf(file.path(here::here(), "data/biomes/olson_harvard_dataverse/tnc_terr_ecoregions.shp"))

# get biome names and code
df_biomes_codes <- biomes |> 
  select(WWF_MHTNAM, WWF_MHTNUM) |> 
  as.data.frame() |> 
  select(-geometry) |> 
  distinct() |> 
  arrange(WWF_MHTNUM)

write_csv(df_biomes_codes, file = here::here("data/biomes/df_biomes_codes.csv"))

# Rasterise ---------------
# Load template raster available in this repository
rasta <- terra::rast(here::here("data/biomes/modis_landcover__LPDAAC__v5.1__0.1deg__2010.nc"))

# perform shapefile to raster conversion
biomes_raster <- terra::rasterize(biomes, rasta, field = "WWF_MHTNAM")

# write to file as NetCDF
terra::writeRaster(biomes_raster, filename = here::here("data/biomes/biomes_raster_0.1deg.tif"), overwrite = TRUE)

# get coastline -------------------
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# Plot map ----------------
data_unm <- readRDS(file.path(here::here(), "/data/data_unm.rds"))
summary(data_unm$plotsize)
data_unm_unique_plots <- data_unm |>
  distinct(plotID, lat, lon)

data_unm |>
  filter(dataset %in% "lwf")

data_unm_euforia <- data_unm |>
  filter(dataset %in% c("bnp", "czu", "forst", "iberbas", "incds", "lwf", "nbw", "uholka", 
                        "nfr", "nwfva", "tuzvo", "ul", "unito", "urk", "wuls")) |>
  distinct(plotID, dataset, lat, lon)

# Keep only the first occurrence of each name
df_labels <- data_unm_euforia %>% group_by(dataset) %>% slice(1) %>% ungroup()

gg <- ggplot() +
  #stat_summary_hex(aes(lon, lat, z = density), data = data_biomes_fil,fun = mean, bins = 50) +
  #coord_sf(crs = pull_crs(biomes_raster, alpha=0.5)) +
  tidyterra::geom_spatraster(data = biomes_raster, alpha=0.5) +
  geom_sf(data = coast,
          colour = 'black',
          linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "Boreal Forests/Taiga"                                         = "dodgerblue4", 
      #"Tundra"                                                       = "lightcyan3"
      #"Deserts and Xeric Shrublands"                                 = "#FFD3A0", 
      #"Flooded Grasslands and Savannas"                              = "indianred3", 
      #"Inland Water"                                                 = "azure", 
      #"Mangroves"                                                    = "violetred", 
      "Mediterranean Forests, Woodlands and Scrub"                   = "orangered3", 
      #"Montane Grasslands and Shrublands"                            = "steelblue3", 
      #"Rock and Ice"                                                 = "azure4", 
      "Temperate Broadleaf and Mixed Forests"                        = "darkgreen", 
      "Temperate Conifer Forests"                                    = "lightseagreen", 
      #"Temperate Grasslands, Savannas and Shrublands"                = "goldenrod3", 
      #"Tropical and Subtropical Coniferous Forests"                  = "#31A278", 
      "Tropical and Subtropical Dry Broadleaf Forests"               = "goldenrod4", 
      #"Tropical and Subtropical Grasslands, Savannas and Shrublands" = "darkolivegreen",
      "Tropical and Subtropical Moist Broadleaf Forests"             = "springgreen3"),
      na.value = "white",
      breaks = ~ .x[!is.na(.x)]) + 
  geom_point(aes(lon, lat), data = data_unm_unique_plots, color="red",fill = "white",alpha = .7, shape=21,size=1.2) + 
  geom_point(aes(lon, lat), data = data_unm_euforia, color="blue",fill = "white",alpha = .7, shape=21,size=1.2) +
  #geom_count(aes(lon, lat), data = data_biomes_fil_unique_plots, color="red",alpha = .6) +
  #scale_size_area() +
  #geom_hex(aes(lon, lat, color = biome), data = data_biomes_fil, bins = 50, linewidth = 1) +
  # set extent in longitude and latitude
  coord_sf(ylim = c(-60, 85),
           expand = FALSE   # to draw map strictly bounded by the specified extent
  ) +
  labs(title = "", x = "", y = "") + 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.5))
gg
ggsave(paste0(here::here(), "/manuscript/figures/fig_S1euforia.png"), width = 13, height = 8, dpi=300)
ggsave(paste0(here::here(), "/manuscript/figures/map.pdf"), width = 5, height = 4.5, dpi=300)

