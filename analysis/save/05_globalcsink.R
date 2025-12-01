# This script evaluates the mature forests C sink

# load packages
library(dplyr)
library(lme4) 
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
library(ingestr)
library(sf)
library(purrr)
library(tidyr)
#library(here)

# Load functions ----
source(here::here("R/functions.R"))
source(here::here("R/get_drivers_by_biome.R"))
source(here("R/csink_global.R"))

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


# Global csink ----
# For each forest grid in the map:
# 1. Sample QMD from its distribution in the total dataset which includes all stands (with and without info about biomass) → QMDj and log-transform
# 2. Estimate mean N given QMDj and two consecutive years (e.g. 2000, 2001) from the LMM relating N and QMD using the total dataset
# 3. Using the subset of data with biomass information, estimate the change in biomass per ha given the  QMDj, N0 and N1 as dB = ak * QMDj^2 * (N1 - N0), 
# and forcing the relationship through the origin (adding +0 to the LMM), where ak is a sampled value from the fitted a, considering its standard error and a normal distribution.
# 4. Repeat steps 1-3 multiple times. This gives the distribution of mature forest biomass change per unit area. 
# the function csink includes the steps 1-3. We run it 1e5 calling the fc using purrr::map_dfr

## Load data ----------------
data_forest_plots <- readRDS(here::here("data/data_fil_biomes.rds")) |>
  filter(year > 1980) |>    # XXX why this filter?
  mutate(NQMD2 = density * QMD^2)

# data_all <- data_fil_biomes
# 
# data_biomass <- data_fil_biomes |>
#   mutate(NQMD2 = density * QMD^2)

## Fit self-thinning relationship ---------------------
# modified by time and environmental covariates
fit_selfthinning = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) * scale(ai) + 
    scale(year) * scale(ndep) + 
    scale(year) * scale(ORGC) + 
    scale(year) * scale(PBR) + 
    (1|dataset/plotID) + (1|species),  
  data = data_fil_biomes
  )


## Fit relationship of biomass ~ N*QMD^2 ---------------------
fit_biomass = lmer(
  biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
  data = data_fil_biomes, 
  na.action = "na.exclude"
  )

a_mean <- summary(fit_biomass)$coefficient[1,1]
a_sd <- summary(fit_biomass)$coefficient[1,2]


# sample logQMD for each gridcell separately  
generate_samples_by_gridcell <- function(idx, global_drivers, data_forest_plots, n_qmd){
  global_drivers |> 
    slice(idx) |> 
    slice(rep(1, n_qmd)) |> 
    mutate(idx_qmd = row_number()) |> 
    mutate(logQMD = sample(data_forest_plots$logQMD, size = n_qmd))
}

n_qmd <- 3
df_samples_qmd <- purrr::map_dfr(
  as.list(1:5),  # doing it only for three gridcells
  ~generate_samples_by_gridcell(., global_drivers, data_forest_plots, n_qmd) 
)

# All get the sampled coefficients of the two fitted models
# Sample coefficients independently
n_coefs <- 2
b_ai <- rnorm(n_coefs, coef_ai_mean, coef_ai_sd)

# Repeat and tag
df_samples_coef <- map_dfr(b_ai, function(b_ai) {
  df_samples_qmd |> 
    mutate(b_ai = b_ai)
})


df_samples <- df_samples |> 
  mutate(
    b_ai = rnorm(1, coef_ai_mean, coef_ai_sd),
    b_ndep = rnorm(1, coef_ndep_mean, coef_ndep_sd),
    b_orgc = rnorm(1, coef_orgc_mean, coef_orgc_sd),
    b_pbr = rnorm(1, coef_pbr_mean, coef_pbr_sd),
    b_aiyear = rnorm(1, coef_aiyear_mean, coef_aiyear_sd),
    b_ndepyear = rnorm(1, coef_ndepyear_mean, coef_ndepyear_sd),
    b_orgcyear = rnorm(1, coef_orgcyear_mean, coef_orgcyear_sd),
    b_pbryear = rnorm(1, coef_pbryear_mean, coef_pbryear_sd)
  )





# global_drivers <- readRDS(file.path(here::here(), "/data/inputs/global_drivers.rds")) |> 
#   drop_na()

# global_drivers <- global_drivers[1:2,] # test

# Function to apply csink 30 times for each row using purrr::map_dfr
apply_csink_global_n_times <- function(global_drivers_row, data_fil_biomes, n_times) {
  map_dfr(
    1:n_times, 
    ~csink_global(
      global_drivers_row, 
      data_fil_biomes,
      coef_ai_mean, coef_ai_sd,
      coef_ndep_mean, coef_ndep_sd,
      coef_orgc_mean, coef_orgc_sd,
      coef_pbr_mean, coef_pbr_sd,
      coef_aiyear_mean, coef_aiyear_sd,
      coef_ndepyear_mean, coef_ndepyear_sd,
      coef_orgcyear_mean, coef_orgcyear_sd,
      coef_pbryear_mean, coef_pbryear_sd), 
    .id = "iteration")
}

# Using pmap_dfr to apply the function to each row of data_to_iterate (global_drivers) and iterate 30 times
system.time(
  results <- global_drivers %>%
    pmap_dfr(function(...) {
      global_drivers_row <- list(...)
      apply_csink_global_n_times(global_drivers_row, data_all, n_times = 100)
    })
)

# View the final results
# print(results)
saveRDS(results, file = file.path(here::here(), "/data/inputs/results_global_csink_v100.rds"))

#results2 <- readRDS(file.path(here::here(), "data/inputs/final_old/results_global_csink.rds"))
results <- readRDS(file.path(here::here(), "data/inputs/results_global_csink.rds"))
results <- readRDS(file.path(here::here(), "data/inputs/results_global_csink_v50.rds"))
results <- readRDS(file.path(here::here(), "data/inputs/results_global_csink_v100.rds"))

# db_Mg_ha ----

agg_results <- results |>
  filter(dB_Mg_ha > 0) |>
  group_by(lon, lat, area_ha) |>
  summarise(dB_Mg_ha = mean(dB_Mg_ha, na.rm=T))

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = dB_Mg_ha), data = agg_results, size = 0.5, alpha=0.5) +
  scale_color_viridis()

# dB per yr ----
# Total mature C sink

# load modis fraction forest cover raster
r_fcf <- terra::rast("/home/laura/data/forest_fraction/MODIS_ForestCoverFraction.nc")

# select only the forestcoverfraction
r_fcf <- r_fcf[[1]]
plot(r_fcf)

# Convert df to SpatVector
points <- vect(agg_results, geom = c("lon", "lat"), crs = crs(r_fcf))

# Extract raster values at given points
extracted <- extract(r_fcf, points, fun = NULL, na.rm = TRUE, touches = TRUE)

# Combine extracted values with polygon attributes
agg_results$fcf <- extracted$forestcoverfraction  # Assuming 'layer' contains the raster values

agg_results <- agg_results |>
  mutate(db_Pg_yr = dB_Mg_ha*1e-9*area_ha*fcf*1e-2)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = db_Pg_yr), data = agg_results, size = 0.5, alpha=0.5) +
  scale_color_viridis()

# Changes in C
agg_results <- agg_results |>
  mutate(dC_Mg_ha = dB_Mg_ha*0.5,
         dC_Pg_yr = db_Pg_yr*0.5)

## Value for poster! 
sum(agg_results$dC_Pg_yr)

fig4 <- rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = dC_Mg_ha), data = agg_results, size = 0.5, alpha=0.5) +
  scale_color_viridis(breaks = seq(0, 2, 1), limits = c(0,2.5)) +
  labs(title = "Forest carbon increase per ha and year",
       color = expression(paste("Mg C ", ha^-1, " ", yr^-1))) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),axis.title = element_text(size = 13),
        axis.text.y = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),legend.title = element_text(size = 10),
        plot.title = element_text(size = 12),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.1, 0.25),
        legend.direction="vertical",
        legend.box = "horizontal",
        legend.margin = margin(2, 2, 2, 2),
        legend.key.size = unit(.4, 'cm'),
        legend.box.margin = margin(1, 1, 1, 1)) 
fig4

#ggsave(paste0(here::here(), "/manuscript/figures/fig_4.png"), width = 13, height = 8, dpi=300)

# Figure 4 ----
coast <- rnaturalearth::ne_coastline(
  scale = 110,
  returnclass = "sf"
)

fig4 <- agg_results |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = dC_Mg_ha),
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
    name =  expression(paste("Mg C ", ha^-1, " ", yr^-1)),
    #option = "cividis",
    limits = c(0, 2.5), breaks = seq(0, 2, 1)
  ) +
  theme_void() +
  labs(
    #subtitle = "Global forest carbon increase per ha and year"
  ) +  
  theme(legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.1, 0.25),
        legend.direction="vertical",
        legend.box = "horizontal",
        legend.margin = margin(2, 2, 2, 2),
        legend.key.size = unit(.6, 'cm'),
        legend.box.margin = margin(1, 1, 1, 1)) 
fig4
ggsave(paste0(here::here(), "/manuscript/figures/fig_4.png"),  plot = fig4, width = 13, height = 8, dpi=300)

### Get distributions ----

agg_results 
results

# load modis fraction forest cover raster
r_fcf <- terra::rast("/home/laura/data/forest_fraction/MODIS_ForestCoverFraction.nc")

# select only the forestcoverfraction
r_fcf <- r_fcf[[1]]
plot(r_fcf)

# Convert df to SpatVector
points <- vect(results, geom = c("lon", "lat"), crs = crs(r_fcf))

# Extract raster values at given points
extracted <- extract(r_fcf, points, fun = NULL, na.rm = TRUE, touches = TRUE)

# Combine extracted values with polygon attributes
results$fcf <- extracted$forestcoverfraction  # Assuming 'layer' contains the raster values

results <- results |>
  mutate(db_Pg_yr = dB_Mg_ha*1e-9*area_ha*fcf*1e-2,
         dC_Mg_ha = dB_Mg_ha*0.5,
         dC_Pg_yr = db_Pg_yr*0.5) 

# Sum across the same row_id from each group
agg_results_total <- results %>%
  group_by(lon, lat) |>
  mutate(row_id = row_number()) |>
  ungroup() |>
  group_by(row_id) %>%
  summarise(dC_Mg_ha_global = mean(dC_Mg_ha, na.rm = T),
            dC_Pg_yr_global = sum(dC_Pg_yr, na.rm = T)) 

dC_Pg_yr_global <- agg_results_total |> 
  ggplot(aes(dC_Pg_yr_global, ..density..)) +
  geom_density(color="#D55E00",fill="#D55E00", alpha =0.7) +
  theme_classic() +
  labs(#title = "Total forest carbon per year",
    x = expression(paste("Pg C ", yr^-1)), y = "Density") +
  scale_x_continuous(limits = c(1.36,1.39)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),axis.title = element_text(size = 12),
        axis.text.y = element_text(hjust = 0.5)) 
dC_Pg_yr_global

mean(agg_results_total$dC_Pg_yr_global)
sd(agg_results_total$dC_Pg_yr_global)


dC_Pg_yr_both <- ggplot() +
  geom_density(data= out_csink_agg, aes(dC_Pg_yr_allbiomes, ..density..), color="#999999",fill="#999999", alpha =0.7) +
  geom_density(data= agg_results_total, aes(dC_Pg_yr_global, ..density..), color="#D55E00",fill="#D55E00", alpha =0.7) +
  theme_classic() +
  labs(#title = "Total forest carbon per year",
    x = expression(paste("Pg C ", yr^-1)), y = "Density") +
  scale_x_continuous(limits = c(0.6,1.6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),axis.title = element_text(size = 12),
        axis.text.y = element_text(hjust = 0.5)) 
dC_Pg_yr_both

# Figure 3 ----
fig3 <- dC_Mg_ha_allbiomes  + dC_Pg_yr + dC_Pg_yr_allbiomes + dC_Pg_yr_global +
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "a",tag_suffix = ")")
fig3
ggsave(paste0(here::here(), "/manuscript/figures/fig3.png"), width = 11, height = 8, dpi=300)

