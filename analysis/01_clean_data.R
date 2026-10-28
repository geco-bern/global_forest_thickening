# Compile stand records for all datasets, before environmental extraction.

# load packages ----
# library(renv)
library(tidyverse)

# load functions ----
source(here::here("R/functions.R"))

# NFI Spain ----
df_nfi_spain <- readRDS(here::here("data/inputs/df_nfi_spa.rds"))

# NFI Sweeden ----
df_nfi_sweeden <- readRDS(here::here("data/inputs/df_nfi_swe.rds"))

# FIA US ----
df_fia_us <- readRDS(here::here("data/inputs/df_fia_us.rds"))

# NFI Switzerland ----
df_nfi_switzerland <- readRDS(here::here("data/inputs/df_nfi_swi.rds"))

# NFI Norway ----
df_nfi_norway <- readRDS(here::here("data/inputs/df_nfi_nor.rds"))

# EFM Switzerland ----
df_efm_swi <- readRDS(here::here("data/inputs/df_efm_swi.rds"))

# Uholka plot ----
df_uholka <- readRDS(here::here("data/inputs/df_uholka.rds"))

# Greece sites ----
df_fep_gre <- readRDS(here::here("data/inputs/df_fep_gre.rds"))

# France sites ----
df_inrae_lessem <- readRDS(here::here("data/inputs/df_inrae_lessem.rds"))

# EuFoRia plots ----

## bnp ----
df_bnp <- readRDS(here::here("data/inputs/df_bnp.rds"))

## czu ----
df_czu <- readRDS(here::here("data/inputs/df_czu.rds"))

## fvabw ----
df_fvabw <- readRDS(here::here("data/inputs/df_fvabw.rds"))

## iberbas ----
df_iberbas <- readRDS(here::here("data/inputs/df_iberbas.rds"))

## unitbv ----
df_unitbv <- readRDS(here::here("data/inputs/df_unitbv.rds"))

## lwf ----
df_lwf <- readRDS(here::here("data/inputs/df_lwf.rds"))

## npvbw ----
df_npvbw <- readRDS(here::here("data/inputs/df_npvbw.rds"))

## nfr ----
df_nfr_swi <- readRDS(here::here("data/inputs/df_nfr_swi.rds"))

## nwfva ----
df_nwfva <- readRDS(here::here("data/inputs/df_nwfva.rds"))

## tuzvo ----
df_tuzvo <- readRDS(here::here("data/inputs/df_tuzvo.rds"))

## ul ----
df_ul <- readRDS(here::here("data/inputs/df_ul.rds"))

## unito ----
df_unito <- readRDS(here::here("data/inputs/df_unito.rds"))

## urk ----
df_urk <- readRDS(here::here("data/inputs/df_urk.rds"))

## wuls ----
df_wuls <- readRDS(here::here("data/inputs/df_wuls.rds"))

# ForestGEO ----

## Luquillo ----
df_luquillo <- readRDS(here::here("data/inputs/df_luquillo.rds"))

## BCI ----
df_bci <- readRDS(here::here("data/inputs/df_bci.rds"))

## SCBI ----
df_scbi <- readRDS(here::here("data/inputs/df_scbi.rds"))

## Palanam ----
df_palanam <- readRDS(here::here("data/inputs/df_palanam.rds"))

## SERC ----
df_serc <- readRDS(here::here("data/inputs/df_serc.rds"))

## Wytham woods ----
df_wytham <- readRDS(here::here("data/inputs/df_wytham.rds"))

## Pasoh woods ----
df_pasoh <- readRDS(here::here("data/inputs/df_pasoh.rds"))

## Mudumalai ----
df_mudumalai <- readRDS(here::here("data/inputs/df_mudumalai.rds"))

# Forest plots ----
df_forestplots <- readRDS(here::here("data/inputs/df_forestplots.rds"))

# Australian plots ----
df_aus <- readRDS(here::here("data/inputs/df_aus.rds"))

# RAINFOR plots ----
df_rainfor <- readRDS(here::here("data/inputs/df_rainfor.rds"))

# Write data from all plots ----

# join all datasets
data_all <- bind_rows(
  df_nfi_spain,
  df_nfi_sweeden,
  df_fia_us,
  df_nfi_switzerland,
  df_nfi_norway,
  df_efm_swi,
  df_uholka,
  df_fep_gre,
  df_inrae_lessem,
  df_bnp,
  df_czu,
  df_fvabw,
  df_iberbas,
  df_unitbv,
  df_lwf,
  df_npvbw,
  df_nfr_swi,
  df_nwfva,
  df_tuzvo,
  df_ul,
  df_unito,
  df_urk,
  df_wuls,
  df_luquillo,
  # df_bci,  # not included. Comment in `00_raw_data.R`: finally not added to the global dataset due to difficulty to exclude management (communication with  Salomón Aguilar, Research Technician)
  df_scbi,
  df_palanam,
  df_serc,
  df_wytham,
  df_pasoh,
  df_mudumalai,
  df_forestplots,
  df_aus,
  df_rainfor
)

df_all <- data_all |>
  mutate(dataset = harmonise_dataset_names(dataset)) |>
  # unite plotID and dataset to ensure there are not matches form different sources
  rename(plotIDD = plotID) |>
  unite(plotID, dataset, plotIDD, sep = "_", remove = FALSE) |>
  select(-plotIDD) |>
  # Unify Forestplots and RAINFOR datasets
  # Since Forestplots are deta provided by Y. Mahli and RAINFOR are taken from Esquivel et al., some plots can 
  # be duplicated. We need to check the plots and years that are present in both datasets and keep the most updated ones.
  distinct(country, plotID, year, .keep_all = TRUE)
#drop_na(density) |>
#drop_na(QMD) 

saveRDS(df_all, file = here::here("data/inputs/df_all.rds"))

# Next run the plot-extraction section in 00_collect_global_drivers.R,
# then apply the common analysis filters in 02_create_filters.R.

# Historical analyses retained for reference; superseded by 02_create_filters.R.
# They require environmental and log columns that this compilation no longer adds.
run_legacy_analyses <- FALSE
if (run_legacy_analyses) {
# Write unmanaged data ----
# NEW APPROACH:
# Management categories 1,2,3 depending on the knowledge of the forest use history.

# THIS WAS THE INITIAL APPROACH:
# Select only unmanaged forests. These filters are now applied in 02_create_filters.R
df_unm <- data_unm_fc(df_all)
df_unm
df_unm |>
  nrow()

plot_stl(df_unm)
plot_map(df_unm)

saveRDS(df_unm, file = here::here("data/inputs/df_unm.rds"))

# Figures for SI -------------------------


# Generate filtered by biomes ----
# The unmanaged data is divided by biomes, and then for each biome we apply the filter of the upper quantile

# read data written just above
df_unm <- read_rds(here::here("data/inputs/df_unm.rds"))

## Biome 1 ----
# Tropical & Subtropical Moist Broadleaf Forests
df_unm_biome1 <- df_unm |>
  filter(biomeID == 1)
df_unm_biome1 |>
  nrow()

plot_stl(df_unm_biome1)

# filter data with upper quantile
df_fil75_biome1 <- df_filter75_fc(df_unm_biome1)
plot_stl(df_fil75_biome1)
saveRDS(df_fil75_biome1, file = here::here("data/inputs/df_fil75_biome1.rds"))

df_fil90_biome1 <- df_filter90_fc(df_unm_biome1)
plot_stl(df_fil90_biome1)
saveRDS(df_fil90_biome1, file = here::here("data/inputs/df_fil90_biome1.rds"))

df_fil55_biome1 <- df_filter55_fc(df_unm_biome1)
plot_stl(df_fil55_biome1)
saveRDS(df_fil55_biome1, file = here::here("data/inputs/df_fil55_biome1.rds"))

## Biome 2 ----
# Tropical & Subtropical Dry Broadleaf Forests Forest
df_unm_biome2 <- df_unm |> 
  filter(biomeID==2)
df_unm_biome2 |> 
  nrow()
plot_stl(df_unm_biome2)

# filter data with upper quantile
df_fil75_biome2 <- df_filter75_fc(df_unm_biome2)
plot_stl(df_fil75_biome2)
saveRDS(df_fil75_biome2, file = here::here("data/inputs/df_fil75_biome2.rds"))

df_fil90_biome2 <- df_filter90_fc(df_unm_biome2)
plot_stl(df_fil90_biome2)
saveRDS(df_fil90_biome2, file = here::here("data/inputs/df_fil90_biome2.rds"))

df_fil55_biome2 <- df_filter55_fc(df_unm_biome2)
plot_stl(df_fil55_biome2)
saveRDS(df_fil55_biome2, file = here::here("data/inputs/df_fil55_biome2.rds"))

## Biome 3 ----
# Tropical & Subtropical Coniferous Forests Forest
# NOTE: There are no enough data

## Biome 4 ----
# Temperate Broadleaf & Mixed Forests Forest
df_unm_biome4 <- df_unm |>
  filter(biomeID == 4)
plot_stl(df_unm_biome4)

# filter data with upper quantile
df_fil75_biome4 <- df_filter75_fc(df_unm_biome4)
plot_stl(df_fil75_biome4)
saveRDS(df_fil75_biome4, file = here::here("data/inputs/df_fil75_biome4.rds"))

df_fil90_biome4 <- df_filter90_fc(df_unm_biome4)
plot_stl(df_fil90_biome4)
saveRDS(df_fil90_biome4, file = here::here("data/inputs/df_fil90_biome4.rds"))

df_fil55_biome4 <- df_filter55_fc(df_unm_biome4)
plot_stl(df_fil55_biome4)
saveRDS(df_fil55_biome4, file = here::here("data/inputs/df_fil55_biome4.rds"))

## Biome 5 ----
# Temperate Conifer Forests Forest
df_unm_biome5 <- df_unm |>
  filter(biomeID == 5)
plot_stl(df_unm_biome5)

# filter data with upper quantile
df_fil75_biome5 <- df_filter75_fc(df_unm_biome5)
plot_stl(df_fil75_biome5)
saveRDS(df_fil75_biome5, file = here::here("data/inputs/df_fil75_biome5.rds"))

df_fil90_biome5 <- df_filter90_fc(df_unm_biome5)
plot_stl(df_fil90_biome5)
saveRDS(df_fil90_biome5, file = here::here("data/inputs/df_fil90_biome5.rds"))

df_fil55_biome5 <- df_filter55_fc(df_unm_biome5)
plot_stl(df_fil55_biome5)
saveRDS(df_fil55_biome5, file = here::here("data/inputs/df_fil55_biome5.rds"))

## Biome 6 ----
# Boreal Forests/Taiga Forest
df_unm_biome6 <- df_unm |>
  filter(biomeID == 6)
plot_stl(df_unm_biome6)

# filter data with upper quantile
df_fil75_biome6 <- df_filter75_fc(df_unm_biome6)
plot_stl(df_fil75_biome6)
saveRDS(df_fil75_biome6, file = here::here("data/inputs/df_fil75_biome6.rds"))

df_fil90_biome6 <- df_filter90_fc(df_unm_biome6)
plot_stl(df_fil90_biome6)
saveRDS(df_fil90_biome6, file = here::here("data/inputs/df_fil90_biome6.rds"))

df_fil55_biome6 <- df_filter55_fc(df_unm_biome6)
plot_stl(df_fil55_biome6)
saveRDS(df_fil55_biome6, file = here::here("data/inputs/df_fil55_biome6.rds"))

## Biome 12 ----
# Mediterranean Forests
df_unm_biome12 <- df_unm |>
  filter(biomeID == 12)
plot_stl(df_unm_biome12)

# filter data with upper quantile
df_fil75_biome12 <- df_filter75_fc(df_unm_biome12)
plot_stl(df_fil75_biome12)
saveRDS(df_fil75_biome12, file = here::here("data/inputs/df_fil75_biome12.rds"))

df_fil90_biome12 <- df_filter90_fc(df_unm_biome12)
plot_stl(df_fil90_biome12)
saveRDS(df_fil90_biome12, file = here::here("data/inputs/df_fil90_biome126.rds"))

df_fil55_biome12 <- df_filter55_fc(df_unm_biome12)
plot_stl(df_fil55_biome12)
saveRDS(df_fil55_biome12, file = here::here("data/inputs/df_fil55_biome12.rds"))

# Write filtered data ----
# join filtered datasets
df_fil75_biomes <- bind_rows(
  df_fil75_biome1,
  df_fil75_biome2,
  df_fil75_biome4,
  df_fil75_biome5,
  df_fil75_biome6,
  df_fil75_biome12
)
saveRDS(df_fil75_biomes, file = here::here("data/inputs/df_fil75_biomes.rds"))

df_fil90_biomes <- bind_rows(
  df_fil90_biome1,
  df_fil90_biome2,
  df_fil90_biome4,
  df_fil90_biome5,
  df_fil90_biome6,
  df_fil90_biome12
)
saveRDS(df_fil90_biomes, file = here::here("data/inputs/df_fil90_biomes.rds"))

df_fil55_biomes <- bind_rows(
  df_fil55_biome1,
  df_fil55_biome2,
  df_fil55_biome4,
  df_fil55_biome5,
  df_fil55_biome6,
  df_fil55_biome12
)
saveRDS(df_fil55_biomes, file = here::here("data/inputs/df_fil55_biomes.rds"))

# Plot filtered data --------
df_fil75_biomes <- read_rds(here::here("data/inputs/df_fil75_biomes.rds"))

df_fil75_biomes |>
  drop_na(years_since_management) |>
  distinct(dataset)

# Check DBH distributions
ggplot() +
  geom_density(data = df_fil75_biomes, aes(QMD, ..density.., col = factor(biome), fill = factor(biome)), alpha = 0.7) +
  theme_classic() +
  # color palette from okabeito_colors()
  scale_fill_manual(name = "Biome", values = c(
    "Tropical & Subtropical Moist Broadleaf Forests" = "#E69F00", 
    "Tropical & Subtropical Dry Broadleaf Forests" = "#56B4E9",
    "Temperate Broadleaf & Mixed Forests" = "#009E73",
    "Temperate Conifer Forests" = "#F5C710",
    "Boreal Forests/Taiga" = "#0072B2",
    "Mediterranean Forests, Woodlands & Scrub" = "#D55E00"
  )) +
  scale_color_manual(name = "Biome", values = c(
    "Tropical & Subtropical Moist Broadleaf Forests" = "#E69F00",
    "Tropical & Subtropical Dry Broadleaf Forests" = "#56B4E9",
    "Temperate Broadleaf & Mixed Forests" = "#009E73",
    "Temperate Conifer Forests" = "#F5C710",
    "Boreal Forests/Taiga" = "#0072B2",
    "Mediterranean Forests, Woodlands & Scrub" = "#D55E00"
  )) +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 80)) 

# Table S1 ----
data <- readRDS(here::here("/data/inputs/df_unm.rds"))
data <- readRDS(here::here("/data/inputs/df_fil75_biomes.rds"))

# Unmanaged years
data |>
  filter(years_since_management != 999) |>
  summarise(mean=mean(years_since_management, na.rm=T),
            min=min(years_since_management, na.rm=T),
            max=max(years_since_management, na.rm=T))
# plot size
data |>
  summarise(mean=mean(plotsize, na.rm=T),
            min=min(plotsize, na.rm=T),
            max=max(plotsize, na.rm=T))

# Mean year for first and last census per plot
data %>%
  group_by(plotID) %>%
  filter(census == 1) %>%
  ungroup() |>
  summarise(mean=mean(year, na.rm=T))

data %>%
  group_by(plotID) %>%
  filter(census == max(census, na.rm = TRUE)) %>%
  ungroup() |>
  summarise(mean=mean(year, na.rm=T))

species <- as.data.frame(sort(unique(df_unm$species)))
length(unique(df_unm$plotID))
df_unm
summary(data)

df_unm |>
  filter(plotID == "nwfva_tree_03-025_2")

# Table at the plot level
table_s1 <- data %>% group_by(country,dataset) %>% 
  summarise(n_plots=n_distinct(plotID),
            min_lon = min(lon),
            max_lon = max(lon),
            min_lat = min(lat),
            max_lat = max(lat),
            min_year = min(year),
            max_year = max(year),
            timespan=max_year-min_year,
            n_census_mean = as.integer(mean(n_census,na.rm=T)),
            interval=as.integer(mean(period,na.rm=T)),
            interval_sd =sd(period,na.rm=T),
            qmd_mean=mean(QMD,na.rm=T),
            QMD_sd=sd(QMD,na.rm=T),
            density_mean=mean(density,na.rm=T),
            density_sd=sd(density,na.rm=T),
            biomass_mean=mean(biomass,na.rm=T),
            biomass_sd=sd(biomass,na.rm=T))
table_s1
write.csv(table_s1, file = here::here("/data/outputs/tableS1.csv"))

# Table S2 ----
# Table at the biome level
table_s2 <- data %>% 
  group_by(biome) %>% 
  summarise(n_plots=n_distinct(plotID),
            min_year = min(year),
            max_year = max(year),
            timespan=mean(period,na.rm=T),
            timespan_sd =sd(period,na.rm=T),
            qmd_mean=mean(QMD,na.rm=T),
            QMD_sd=sd(QMD,na.rm=T),
            density_mean=mean(density,na.rm=T),
            density_sd=sd(density,na.rm=T))
table_s2
write.csv(table_s2, file = here::here("/data/outputs/tableS2.csv"))

}
