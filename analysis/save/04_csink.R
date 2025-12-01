# This script evaluates the mature forests C sink for each forest  biome represented in the data

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

# load functions ----
source(here::here("R/csink.R"))

# load data by biome ----

# biome 1
# Tropical & Subtropical Moist Broadleaf Forests
data_fil_biome1 <- readRDS(here::here("data/data_fil_biome1.rds"))

data_biomass_biome1 <- data_fil_biome1 |>
  mutate(NQMD2 = density * QMD^2) # new variable for the csink model

data_biomass_biome1 |>
  distinct(dataset)

data_biomass_biome1 |>
  ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + 
  geom_point() + 
  coord_obs_pred()

# biome 4
# Temperate Broadleaf & Mixed Forests  
data_fil_biome4 <- readRDS(here::here("data/data_fil_biome4.rds"))

data_biomass_biome4 <- data_fil_biome4 |>
  mutate(NQMD2 = density * QMD^2) # new variable for the csink model

data_biomass_biome4 |>
  distinct(dataset) 

data_biomass_biome4 |>
  ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + 
  geom_point() + 
  coord_obs_pred()
  
# biome 5
# Temperate Conifer Forests  
data_fil_biome5 <- readRDS(here::here("data/data_fil_biome5.rds"))

data_biomass_biome5 <- data_fil_biome5 |>
  mutate(NQMD2 = density * QMD^2) # new variable for the csink model

data_biomass_biome5 |>
  distinct(dataset) 

data_biomass_biome5 |>
  ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + 
  geom_point() + 
  coord_obs_pred()

# biome 6
# Boreal Forests/Taiga
data_fil_biome6 <- readRDS(here::here("data/data_fil_biome6.rds"))

data_biomass_biome6 <- data_fil_biome6 |>
  mutate(NQMD2 = density * QMD^2) # new variable for the csink model

data_biomass_biome6 |>
  distinct(dataset) 

data_biomass_biome6 |>
  ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + geom_point() + 
  coord_obs_pred()

# biome 12
# Mediterranean Forests, Woodlands & Scrub
data_fil_biome12 <- readRDS(here::here("data/data_fil_biome12.rds"))

data_biomass_biome12 <- data_fil_biome12 |>
  mutate(NQMD2 = density * QMD^2) # new variable for the csink model

data_biomass_biome12 |>
  distinct(dataset) 

data_biomass_biome12 |>
  ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + geom_point() + coord_obs_pred()

# csink funtion ----
# For each biome:
# 1. Sample QMD from its distribution in the total dataset which includes all stands (with and without info about biomass) → QMDj and log-transform
# 2. Estimate mean N given QMDj and two consecutive years (e.g. 2000, 2001) from the LMM relating N and QMD using the total dataset
# 3. Using the subset of data with biomass information, estimate the change in biomass per ha given the  QMDj, N0 and N1 as dB = ak * QMDj^2 * (N1 - N0), 
# and forcing the relationship through the origin (adding +0 to the LMM), where ak is a sampled value from the fitted a, considering its standard error and a normal distribution.
# 4. Repeat steps 1-3 multiple times. This gives the distribution of mature forest biomass change per unit area. 
# the function csink includes the steps 1-3. We run it 1e5 calling the fc using purrr::map_dfr

## Biome 1 ---------
# Tropical & Subtropical Moist Broadleaf Forests
data_all <- data_fil_biome1
data_biomass <- data_biomass_biome1 

fit1 = lmer(
  log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
  data = data_all, 
  na.action = "na.exclude"
)

fit2 = lmer(
  biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
  data = data_biomass, 
  na.action = "na.exclude"
)

a_mean <- summary(fit2)$coefficient[1,1]
a_sd <- summary(fit2)$coefficient[1,2]

filn <- here::here("data/out_csink_biome1.rds")
if (!file.exists(filn)){
  out <- purrr::map_dfr(
    as.list(seq(1e5)), #1e5
    ~csink(data_all, a_mean, a_sd))
  saveRDS(out, file = filn)
}

## Biome 4 ---------
# Temperate Broadleaf & Mixed Forests 
data_all <- data_fil_biome4
data_biomass <- data_biomass_biome4

fit1 = lmer(
  log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
  data = data_all, 
  na.action = "na.exclude"
)

fit2 = lmer(
  biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
  data = data_biomass, 
  na.action = "na.exclude"
)

a_mean <- summary(fit2)$coefficient[1,1]
a_sd <- summary(fit2)$coefficient[1,2]

filn <- here::here("data/out_csink_biome4.rds")
if (!file.exists(filn)){
  out <- purrr::map_dfr(
    as.list(seq(1e5)), #1e5
    ~csink(data_all, fit1, a_mean, a_sd))
  saveRDS(out, file = filn)
}

## Biome 5 ---------
# Temperate Conifer Forests  
data_all <- data_fil_biome5
data_biomass <- data_biomass_biome5

fit1 = lmer(log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), data = data_all, na.action = "na.exclude")

fit2 = lmer(biomass ~ NQMD2 + 0 + (1|dataset/plotID), data = data_biomass, na.action = "na.exclude")
a_mean <- summary(fit2)$coefficient[1,1]
a_sd <- summary(fit2)$coefficient[1,2]

filn <- here::here("data/out_csink_biome5.rds")
if (!file.exists(filn)){
  out <- purrr::map_dfr(
    as.list(seq(1e5)), #1e5
    ~csink(data_all, fit1, a_mean, a_sd))
  saveRDS(out, file = filn)
}

## Biome 6 ---------
# Boreal Forests/Taiga
data_all <- data_fil_biome6
data_biomass <- data_biomass_biome6

fit1 = lmer(log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), data = data_all, na.action = "na.exclude")

fit2 = lmer(biomass ~ NQMD2 + 0 + (1|plotID), data = data_biomass, na.action = "na.exclude")
a_mean <- summary(fit2)$coefficient[1,1]
a_sd <- summary(fit2)$coefficient[1,2]

filn <- here::here("data/out_csink_biome6.rds")
if (!file.exists(filn)){
  out <- purrr::map_dfr(
    as.list(seq(1e5)), #1e5
    ~csink(data_all, fit1, a_mean, a_sd))
  saveRDS(out, file = filn)
}

## Biome 12 ---------
# Mediterranean Forests, Woodlands & Scrub
data_all <- data_fil_biome12
data_biomass <- data_biomass_biome12

fit1 = lmer(log(density) ~ logQMD + year + (1|plotID) + (1|species), data = data_all, na.action = "na.exclude")

fit2 = lmer(biomass ~ NQMD2 + 0 + (1|plotID), data = data_biomass, na.action = "na.exclude")
a_mean <- summary(fit2)$coefficient[1,1]
a_sd <- summary(fit2)$coefficient[1,2]

filn <- here::here("data/out_csink_biome12.rds")
if (!file.exists(filn)){
  out <- purrr::map_dfr(
    as.list(seq(1e5)), #1e5
    ~csink(data_all, fit1, a_mean, a_sd))
  saveRDS(out, file = filn)
} 

# dB per ha ----
# Biomass change

# Biome 1: Tropical & Subtropical Moist Broadleaf Forests
filn <- here::here("data/out_csink_biome1.rds")

out_csink_biome1 <- out_csink_biome1 |>
  mutate(biome = "Biome 1")

dB_Mg_ha_b1 <- out_csink_biome1 |> 
  ggplot(aes(dB_Mg_ha, ..density..)) +
  geom_density() +
  theme_classic() +
  labs(title = "Biomass change",
       subtitle = "Tropical & Subtropical Moist Broadleaf Forests")

dB_Mg_ha_b1

# Biome 4: Temperate Broadleaf & Mixed Forests 
out_csink_biome4 <- readRDS(here::here("data/out_csink_biome4.rds"))
out_csink_biome4 <- out_csink_biome4 |>
  mutate(biome = "Biome 4")
dB_Mg_ha_b4 <- out_csink_biome4 |> 
  ggplot(aes(dB_Mg_ha, ..density..)) +
  geom_density() +
  theme_classic()  +
  labs(title = "Biomass change",
       subtitle = "Temperate Broadleaf & Mixed Forests")
dB_Mg_ha_b4

# Biome 5: Temperate Conifer Forests  
out_csink_biome5 <- readRDS(here::here("data/out_csink_biome5.rds"))
out_csink_biome5 <- out_csink_biome5 |>
  mutate(biome = "Biome 5")
dB_Mg_ha_b5 <- out_csink_biome5 |> 
  ggplot(aes(dB_Mg_ha, ..density..)) +
  geom_density() +
  theme_classic()  +
  labs(title = "Biomass change",
       subtitle = "Temperate Conifer Forests")
dB_Mg_ha_b5

# Biome 6: Boreal Forests/Taiga
out_csink_biome6 <- readRDS(here::here("data/out_csink_biome6.rds"))
out_csink_biome6 <- out_csink_biome6 |>
  mutate(biome = "Biome 6")
dB_Mg_ha_b6 <- out_csink_biome6 |> 
  ggplot(aes(dB_Mg_ha, ..density..)) +  
  geom_density() +
  theme_classic() +
  labs(title = "Biomass change",
       subtitle = "Boreal Forests/Taiga")
dB_Mg_ha_b6

# Biome 12: Mediterranean forests
out_csink_biome12 <- readRDS(here::here("data/out_csink_biome12.rds"))
out_csink_biome12 <- out_csink_biome12 |>
  mutate(biome = "Biome 12")
dB_Mg_ha_b12 <- out_csink_biome12 |> 
  ggplot(aes(dB_Mg_ha, ..density..)) +  
  geom_density() +
  theme_classic()  +
  labs(title = "Biomass change",
       subtitle = "Mediterranean forests")
dB_Mg_ha_b12

# dB per yr ----
# Total mature C sink

# Get forest fraction cover for each biome

# load biomes vector
# From the WWF Ecoregions data 
v_biomes <- terra::vect(file.path(here::here(), "/data/wwf/wwf_terr_ecos.shp")) # XXX path on GECO WS
values(v_biomes)
# Check the CRS
crs_info <- crs(v_biomes)

# load modis fraction forest cover raster
r_fcf <- terra::rast("/home/laura/data/forest_fraction/MODIS_ForestCoverFraction.nc") # XXX path on GECO WS
names(r_fcf)
# select only the forestcoverfraction
r_fcf <- r_fcf[[1]]
plot(r_fcf)

# Ensure CRS Alignment
# Reproject polygon to raster CRS, if needed
v_biomes <- project(v_biomes, crs(r_fcf))  

# Calculate area of the grid in projected units
v_biomes$area_m2 <- expanse(v_biomes, unit = "m")
values(v_biomes)

# Extract raster values for each polygon
extracted <- extract(r_fcf, v_biomes, fun = NULL, na.rm = TRUE, touches = TRUE)

# Combine extracted values with polygon attributes
v_biomes$fcf <- extracted$forestcoverfraction  # Assuming 'layer' contains the raster values
values(v_biomes)

# View combined polygon attributes
biomes_fcf <- as.data.frame(v_biomes) |>
  mutate(fcf_perc = fcf*1e-2,
        forestcover_area_ha = area_m2*1e-4*fcf_perc) |>
  group_by(BIOME) |>
  summarise(total_forestcover_area_ha = sum(forestcover_area_ha, na.rm = T)) |>
  mutate(forestcover_scinot = (format(total_forestcover_area_ha, scientific = TRUE)))

# 5. Taken the dB estimates per biome, multiple this by the total forest area within this biome 
# to get the distribution of the biome-level total mature forest C sink.
# Mg = 1e-9 Petagrams. 1PG of C = 1e15 grams

# # Biome 1: Tropical & Subtropical Moist Broadleaf Forests
biomes_fcf_b1 <- biomes_fcf |> filter(BIOME==1)

out_csink_biome1 <- readRDS(here::here("data/out_csink_biome1.rds"))
out_csink_biome1 <- out_csink_biome1 |>
  mutate(dB_Mg_total = dB_Mg_ha * (biomes_fcf |> filter(BIOME==1))$total_forestcover_area_ha,
         dB_Pg_yr = dB_Mg_total*1e-9) |>
  mutate(biome = "Tropical Moist Broadleaf Forests")

dB_Pg_yr_b1 <- out_csink_biome1 |> 
  ggplot(aes(dB_Pg_yr, ..density..)) +
  geom_density() +
  theme_classic()  +
  labs(title = "Mature forest C sink",
       subtitle = "Tropical Moist Broadleaf Forests")
dB_Pg_yr_b1

# Biome 4: Temperate Broadleaf & Mixed Forests 
biomes_fcf_b4 <- biomes_fcf |> filter(BIOME==4)

out_csink_biome4 <- readRDS(here::here("data/out_csink_biome4.rds"))
out_csink_biome4 <- out_csink_biome4 |>
  mutate(dB_Mg_total = dB_Mg_ha * (biomes_fcf |> filter(BIOME==4))$total_forestcover_area_ha,
         dB_Pg_yr = dB_Mg_total*1e-9) |>
  mutate(biome = "Temperate Broadleaf & Mixed Forests")
dB_Pg_yr_b4 <- out_csink_biome4 |> 
  ggplot(aes(dB_Pg_yr, ..density..)) +
  geom_density() +
  theme_classic()  +
  labs(title = "Mature forest C sink",
       subtitle = "Temperate Broadleaf & Mixed Forests")
dB_Pg_yr_b4

# Biome 5: Temperate Conifer Forests  
biomes_fcf_b5 <- biomes_fcf |> filter(BIOME==5)

out_csink_biome5 <- readRDS(here::here("data/out_csink_biome5.rds"))
out_csink_biome5 <- out_csink_biome5 |>
  mutate(dB_Mg_total = dB_Mg_ha * (biomes_fcf |> filter(BIOME==5))$total_forestcover_area_ha,
         dB_Pg_yr = dB_Mg_total*1e-9) |>
  mutate(biome = "Temperate Conifer Forests")
dB_Pg_yr_b5 <- out_csink_biome5 |> 
  ggplot(aes(dB_Pg_yr, ..density..)) +
  geom_density() +
  theme_classic() +
  labs(title = "Mature forest C sink",
       subtitle = "Temperate Conifer Forests")
dB_Pg_yr_b5

# Biome 6: Boreal Forests/Taiga
biomes_fcf_b6 <- biomes_fcf |> filter(BIOME==6)

out_csink_biome6 <- readRDS(here::here("data/out_csink_biome6.rds"))
out_csink_biome6 <- out_csink_biome6 |>
  mutate(dB_Mg_total = dB_Mg_ha * (biomes_fcf |> filter(BIOME==6))$total_forestcover_area_ha,
         dB_Pg_yr = dB_Mg_total*1e-9) |>
  mutate(biome = "Boreal Forests/Taiga")
dB_Pg_yr_b6 <- out_csink_biome6 |> 
  ggplot(aes(dB_Pg_yr, ..density..)) +
  geom_density() +
  theme_classic() +
  labs(title = "Mature forest C sink",
       subtitle = "Boreal Forests/Taiga")
dB_Pg_yr_b6

# Biome 12: Mediterranean forests
biomes_fcf_b12 <- biomes_fcf |> filter(BIOME==12)

out_csink_biome12 <- readRDS(here::here("data/out_csink_biome12.rds"))
out_csink_biome12 <- out_csink_biome12 |>
  mutate(dB_Mg_total = dB_Mg_ha * (biomes_fcf |> filter(BIOME==12))$total_forestcover_area_ha,
         dB_Pg_yr = dB_Mg_total*1e-9) |>
  mutate(biome = "Mediterranean forests")
dB_Pg_yr_b12 <- out_csink_biome12 |> 
  ggplot(aes(dB_Pg_yr, ..density..)) +
  geom_density() +
  theme_classic() +
  labs(title = "Mature forest C sink",
       subtitle = "Mediterranean forests")
dB_Pg_yr_b12

# All biomes
out_csink <- out_csink_biome1 |>
  bind_rows(out_csink_biome4) |>
  bind_rows(out_csink_biome5) |>
  bind_rows(out_csink_biome6) |>
  bind_rows(out_csink_biome12) |>
  mutate(dC_Mg_ha = dB_Mg_ha*0.5) |>
  mutate(dC_Pg_yr = dB_Pg_yr*0.5) |>
  mutate(biome = factor(biome, levels =c("Tropical Moist Broadleaf Forests","Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
                                         "Boreal Forests/Taiga","Mediterranean forests")))

# Biomass changes ----
dB_Mg_ha <- out_csink |> 
  ggplot(aes(dB_Mg_ha, ..density.., col=factor(biome), fill=factor(biome))) +
  geom_density() +
  theme_classic() +
  labs(title = "Biomass change per ha and year",
       x = expression(paste("Mg wood ", ha^-1, " ",yr^-1)), y = "Density") +
  scale_fill_okabe_ito(name = "Biome", alpha = .7) +
  scale_color_okabe_ito(name = "Biome", alpha = .7) +
  scale_x_continuous(limits = c(0,4.0)) +
  theme(legend.position = "bottom")
dB_Mg_ha

dB_Pg_yr <- out_csink |> 
  ggplot(aes(dB_Pg_yr, ..density.., col=biome, fill=biome)) +
  geom_density() +
  theme_classic() +
  labs(title = "Total biomass change per yr",
       x = expression(paste("Pg wood per biome ", yr^-1)), y = "Density") +
  scale_fill_okabe_ito(name = "Biome", alpha = .7) +
  scale_color_okabe_ito(name = "Biome", alpha = .7) +
  scale_x_continuous(limits = c(0,2)) +
  theme(legend.position = "bottom")
dB_Pg_yr

# C sink changes ----
dC_Mg_ha <- out_csink |> 
  ggplot(aes(dC_Mg_ha, ..density.., col=factor(biome), fill=factor(biome))) +
  geom_density() +
  theme_classic() +
  labs(#title = "Forest carbon per ha and year",
       x = expression(paste("Mg C ", ha^-1, " ",yr^-1)), y = "Density") +
  scale_fill_okabe_ito(name = "Biome", alpha = .7) +
  scale_color_okabe_ito(name = "Biome", alpha = .7) +
  scale_x_continuous(limits = c(0,2)) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       axis.text = element_text(size = 12),axis.title = element_text(size = 12),
       axis.text.y = element_text(hjust = 0.5),
       legend.text = element_text(size = 10),legend.title = element_text(size = 10),
       plot.title = element_text(size = 12),
       legend.key = element_rect(fill = NA, color = NA),
       legend.position = c(0.65, 0.75),
       legend.direction="vertical",
       legend.box = "horizontal",
       legend.margin = margin(2, 2, 2, 2),
       legend.key.size = unit(.6, 'cm'),
       legend.box.margin = margin(1, 1, 1, 1)) 
dC_Mg_ha

dC_Pg_yr <- out_csink |> 
  ggplot(aes(dC_Pg_yr, ..density.., col=biome, fill=biome)) +
  geom_density() +
  theme_classic() +
  labs(#title = "Forest carbon per biome and year",
       x = expression(paste("Pg C ", yr^-1)), y = "Density") +
  scale_fill_okabe_ito(name = "Biome", alpha = .7) +
  scale_color_okabe_ito(name = "Biome", alpha = .7) +
  scale_x_continuous(limits = c(0,1)) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),axis.title = element_text(size = 12),
        axis.text.y = element_text(hjust = 0.5),
        legend.text = element_text(size = 9),legend.title = element_text(size = 10),
        plot.title = element_text(size = 12),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = "none", #c(0.65, 0.75),
        legend.direction="vertical",
        legend.box = "horizontal",
        legend.margin = margin(2, 2, 2, 2),
        legend.key.size = unit(.5, 'cm'),
        legend.box.margin = margin(1, 1, 1, 1)) 
dC_Pg_yr

# Sum across the same row_id from each group
biomes_fcf_observed <- biomes_fcf |> 
  filter(BIOME==1|BIOME==4|BIOME==5|BIOME==6|BIOME==12)

out_csink_agg <- out_csink |>
  group_by(biome) |>
  mutate(row_id = row_number()) |>
  ungroup() |>
  group_by(row_id) |>
  summarise(dB_Mg_ha_meanbiomes0 = mean(dB_Mg_ha),
            dB_Mg_ha_allbiomes = sum(dB_Mg_total)/sum(biomes_fcf_observed$total_forestcover_area_ha), #weighted mean 
            dB_Pg_yr_allbiomes = sum(dB_Pg_yr)) |>
  mutate(dC_Mg_ha_allbiomes = dB_Mg_ha_allbiomes*0.5,
         dC_Pg_yr_allbiomes = dB_Pg_yr_allbiomes*0.5)
out_csink_agg

# add the agreggated all biomes as new raws for plotting
out_csink_plot <- out_csink |>
  bind_rows(data.frame(dC_Mg_ha = out_csink_agg$dC_Mg_ha_allbiomes, biome = "Observed biomes")) |>
  mutate(biome = factor(biome, levels =c("Observed biomes","Tropical Moist Broadleaf Forests","Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
                                       "Boreal Forests/Taiga","Mediterranean forests")))
unique(out_csink_plot$biome)

# means for paper:
out_csink_plot |>
  group_by(biome) |>
  summarise(dC_Mg_ha_mean = mean(dC_Mg_ha),
            dC_Mg_ha_sd = sd(dC_Mg_ha),
            dC_Pg_yr_mean = mean(dC_Pg_yr),
            dC_Pg_yr_sd = sd(dC_Pg_yr))

dC_Mg_ha_allbiomes <- ggplot() +
  geom_density(data = out_csink_plot, aes(dC_Mg_ha, ..density.., col=factor(biome), fill=factor(biome)), alpha= 0.7) +
  theme_classic() +
  labs(#title = "Forest carbon per ha and year",
    x = expression(paste("Mg C ", ha^-1, " ",yr^-1)), y = "Density") +
  #scale_fill_okabe_ito(name = "Biome", alpha = .7) +
  #scale_color_okabe_ito(name = "Biome", alpha = .7) +
  scale_fill_manual(name = "Biome", values = c("Observed biomes" = "#999999",
                                "Tropical Moist Broadleaf Forests" = "#E69F00", 
                                "Temperate Broadleaf & Mixed Forests" = "#56B4E9",
                                "Temperate Conifer Forests" = "#009E73", 
                                "Boreal Forests/Taiga" = "#F5C710", 
                                "Mediterranean forests" = "#0072B2")) +
  scale_color_manual(name = "Biome", values = c("Observed biomes" = "#999999",
                               "Tropical Moist Broadleaf Forests" = "#E69F00", 
                               "Temperate Broadleaf & Mixed Forests" = "#56B4E9",
                               "Temperate Conifer Forests" = "#009E73", 
                               "Boreal Forests/Taiga" = "#F5C710", 
                               "Mediterranean forests" = "#0072B2")) +
  scale_x_continuous(limits = c(0,2)) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),axis.title = element_text(size = 12),
        axis.text.y = element_text(hjust = 0.5),
        legend.text = element_text(size = 10),legend.title = element_text(size = 10),
        plot.title = element_text(size = 12),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.65, 0.75),
        legend.direction="vertical",
        legend.box = "horizontal",
        legend.margin = margin(2, 2, 2, 2),
        legend.key.size = unit(.5, 'cm'),
        legend.box.margin = margin(1, 1, 1, 1))
dC_Mg_ha_allbiomes
  
dC_Pg_yr_allbiomes <- out_csink_agg |> 
  ggplot(aes(dC_Pg_yr_allbiomes, ..density..)) +
  geom_density(color="#999999",fill="#999999", alpha =0.7) +
  theme_classic() +
  labs(#title = "Total forest carbon per year",
       x = expression(paste("Pg C ", yr^-1)), y = "Density") +
  scale_x_continuous(limits = c(0.6,1.6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),axis.title = element_text(size = 12),
        axis.text.y = element_text(hjust = 0.5)) 
dC_Pg_yr_allbiomes

mean(out_csink_agg$dC_Pg_yr_allbiomes)
sd(out_csink_agg$dC_Pg_yr_allbiomes)

# Figure 3abc ----
fig3abc <- dC_Mg_ha_allbiomes + dC_Pg_yr + dC_Pg_yr_allbiomes +
  plot_layout(ncol = 2)
fig3abc
