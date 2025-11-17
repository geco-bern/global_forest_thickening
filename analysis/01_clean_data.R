# This script clean data and prepares them for the analyses.

# load packages ----
# library(renv)
library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(rFIA)
library(patchwork)
library(terra)
library(sf)
library(lme4)
library(lmerTest)
library(ggeffects)
library(effects)
library(sjPlot)
library(measurements)
library(sp)
library(lqmm)
library(ggforce)
library(MuMIn)
library(ingestr)
# library(DescTools)

# load functions ----
source(file.path(here::here(), "/R/functions.R"))

# NFI Spain ----
data_nfi_spain <- readRDS(here::here("data/inputs/data_nfi_spa.rds"))
data_nfi_spain
plot_stl(data_nfi_spain)
plot_map(data_nfi_spain)
ggplot(data_nfi_spain) +
  geom_point(aes(year, biomass))

# NFI Sweeden ----
data_nfi_sweeden <- readRDS(here::here("data/inputs/data_nfi_swe.rds"))
data_nfi_sweeden
plot_stl(data_nfi_sweeden)
plot_map(data_nfi_sweeden)

# FIA US ----
data_fia_us <- readRDS(here::here("data/inputs/data_fia_us.rds"))
data_fia_us
plot_stl(data_fia_us)
plot_map(data_fia_us)

# NFI Switzerland ----
data_nfi_switzerland <- readRDS(here::here("data/inputs/data_nfi_swi.rds"))
data_nfi_switzerland
plot_stl(data_nfi_switzerland)
plot_map(data_nfi_switzerland)

# NFI Norway ----
data_nfi_norway <- readRDS(here::here("data/inputs/data_nfi_nor.rds"))
data_nfi_norway
plot_stl(data_nfi_norway)
plot_map(data_nfi_norway)

# EFM Switzerland ----
data_efm_swi <- readRDS(here::here("data/inputs/data_efm_swi.rds"))
data_efm_swi
plot_stl(data_efm_swi)
plot_map(data_efm_swi)

# Uholka plot ----
data_uholka <- readRDS(here::here("data/inputs/data_uholka.rds"))
data_uholka
plot_stl(data_uholka)
plot_map(data_uholka)

# Greece sites ----
data_greece <- readRDS(here::here("data/inputs/data_fp_gre.rds"))
data_greece
plot_stl(data_greece)
plot_map(data_greece)

# France sites ----
data_france <- readRDS(here::here("data/inputs/data_fp_fra.rds"))
data_france
plot_stl(data_france)
plot_map(data_france)

# EuFoRia plots ----

## bnp ----
data_bnp <- readRDS(here::here("data/inputs/data_euf_bnp.rds"))
data_bnp
plot_stl(data_bnp)
plot_map(data_bnp)

## czu ----
data_czu <- readRDS(here::here("data/inputs/data_euf_czu.rds"))
data_czu
plot_stl(data_czu)
plot_map(data_czu)

## forst ----
data_forst <- readRDS(here::here("data/inputs/data_euf_forst.rds"))
data_forst
plot_stl(data_forst)
plot_map(data_forst)

## iberbas ----
data_iberbas <- readRDS(here::here("data/inputs/data_euf_iberbas.rds"))
data_iberbas
plot_stl(data_iberbas)
plot_map(data_iberbas)

## incds ----
data_incds <- readRDS(here::here("data/inputs/data_euf_incds.rds"))
data_incds
plot_stl(data_incds)
plot_map(data_incds)

## lwf ----
data_lwf <- readRDS(here::here("data/inputs/data_euf_lwf.rds"))
data_lwf
plot_stl(data_lwf)
plot_map(data_lwf)

## nbw ----
data_nbw <- readRDS(here::here("data/inputs/data_euf_nbw.rds"))
data_nbw
plot_stl(data_nbw)
plot_map(data_nbw)

## nfr ----
data_nfr_swi <- readRDS(here::here("data/inputs/data_nfr_swi.rds"))
data_nfr_swi
plot_stl(data_nfr_swi)
plot_map(data_nfr_swi)

## nwfva ----
data_nwfva <- readRDS(here::here("data/inputs/data_euf_nwfva.rds"))
data_nwfva
plot_stl(data_nwfva)
plot_map(data_nwfva)

## tuzvo ----
data_tuzvo <- readRDS(here::here("data/inputs/data_euf_tuzvo.rds"))
data_tuzvo
plot_stl(data_tuzvo)
plot_map(data_tuzvo)

## ul ----
data_ul <- readRDS(here::here("data/inputs/data_euf_ul.rds"))
data_ul
plot_stl(data_ul)
plot_map(data_ul)

## unito ----
data_unito <- readRDS(here::here("data/inputs/data_euf_unito.rds"))
data_unito
plot_stl(data_unito)
plot_map(data_unito)

## urk ----
data_urk <- readRDS(here::here("data/inputs/data_euf_urk.rds"))
data_urk
plot_stl(data_urk)
plot_map(data_urk)

## wuls ----
data_wuls <- readRDS(here::here("data/inputs/data_euf_wuls.rds"))
data_wuls
plot_stl(data_wuls)
plot_map(data_wuls)

# ForestGEO ----

## Luquillo ----
data_luquillo <- readRDS(here::here("data/inputs/data_luquillo.rds"))
data_luquillo
plot_stl(data_luquillo)
plot_map(data_luquillo)

## BCI ----
data_bci <- readRDS(here::here("data/inputs/data_bci.rds"))
data_bci
plot_stl(data_bci)
plot_map(data_bci)

## SCBI ----
data_scbi <- readRDS(here::here("data/inputs/data_scbi.rds"))
data_scbi
plot_stl(data_scbi)
plot_map(data_scbi)

## Palanam ----
data_palanam <- readRDS(here::here("data/inputs/data_palanam.rds"))
data_palanam
plot_stl(data_palanam)
plot_map(data_palanam)

## SERC ----
data_serc <- readRDS(here::here("data/inputs/data_serc.rds"))
data_serc
plot_stl(data_serc)
plot_map(data_serc)

## Wytham woods ----
data_wytham <- readRDS(here::here("data/inputs/data_wytham.rds"))
data_wytham
plot_stl(data_wytham)
plot_map(data_wytham)

## Pasoh woods ----
data_pasoh <- readRDS(here::here("data/inputs/data_pasoh.rds"))
data_pasoh
plot_stl(data_pasoh)
plot_map(data_pasoh)
sort(unique(data_pasoh$species))

## Mudumalai ----
data_mudumalai <- readRDS(here::here("data/inputs/data_mudumalai.rds"))
data_mudumalai
plot_stl(data_mudumalai)
plot_map(data_mudumalai)

# Forest plots ----
data_forestplots <- readRDS(here::here("data/inputs/data_mz_forestplots.rds"))
data_forestplots
plot_stl(data_forestplots)
plot_map(data_forestplots)

# Australian plots ----
data_aus <- readRDS(here::here("data/inputs/data_fp_aus.rds"))
data_aus
plot_stl(data_aus)
plot_map(data_aus)

# RAINFOR plots ----
data_rainfor <- readRDS(file.path(here::here(), "/data/inputs/data_df_rainfor.rds"))
data_rainfor
plot_stl(data_rainfor)
plot_map(data_rainfor)

# All plots ----

# join all datasets
data_all <- bind_rows(
  data_nfi_spain,
  data_nfi_sweeden,
  data_fia_us,
  data_nfi_switzerland,
  data_nfi_norway,
  data_efm_swi,
  data_uholka,
  data_greece,
  data_france,
  data_bnp,
  data_czu,
  data_forst,
  data_iberbas,
  data_incds,
  data_lwf,
  data_nbw,
  data_nfr_swi,
  data_nwfva,
  # data_silava,
  data_tuzvo,
  data_ul,
  data_unito,
  data_urk,
  data_wuls,
  data_luquillo,
  # data_bci,
  data_scbi,
  data_palanam,
  data_serc,
  data_wytham,
  data_pasoh,
  data_mudumalai,
  data_forestplots,
  data_aus,
  data_rainfor
)

data_all <- data_all |>
  mutate(ndep = noy + nhx) |>
  # unite plotID and dataset to ensure there are not matches form different sources
  rename(plotIDD = plotID) |>
  unite(plotID, dataset, plotIDD, sep = "_", remove = FALSE) |>
  select(-plotIDD)
# drop_na(density) |>
# drop_na(QMD) |>
# drop_na(logQMD) |>
# drop_na(CNrt) |>
# drop_na(PBR) |>
# drop_na(ORGC) |>
# drop_na(ndep)

saveRDS(data_all, file = here::here("data/inputs/agg/data_all.rds"))

# Unmanaged ----
# Select only unmanaged forests
data_unm <- data_unm_fc(data_all)
data_unm
data_unm |>
  nrow()
plot_stl(data_unm)
plot_map(data_unm)

saveRDS(data_unm, file = here::here("data/inputs/agg/data_unm.rds"))

# Filtered by biomes ----
# The unmanaged data is divided by biomes, and then for each biome we apply the filter of the upper quantile.

## Biome 1 ----
# Tropical & Subtropical Moist Broadleaf Forests
data_unm_biome1 <- data_unm |>
  filter(biomeID == 1)
data_unm_biome1 |>
  nrow()
plot_stl(data_unm_biome1)

# filter data with upper quantile
data_fil_biome1 <- data_filter_fc(data_unm_biome1)
plot_stl(data_fil_biome1)
saveRDS(data_fil_biome1, file = here::here("data/inputs/data_fil_biome1.rds"))

## Biome 2 ----
# Tropical & Subtropical Dry Broadleaf Forests Forest
# NOTE: There are no enough data (but see Mudumalai)

## Biome 3 ----
# Tropical & Subtropical Coniferous Forests Forest
# NOTE: There are no enough data

## Biome 4 ----
# Temperate Broadleaf & Mixed Forests Forest
data_unm_biome4 <- data_unm |>
  filter(biomeID == 4)
plot_stl(data_unm_biome4)

# filter data with upper quantile
data_fil_biome4 <- data_filter_fc(data_unm_biome4)
plot_stl(data_fil_biome4)
saveRDS(data_fil_biome4, file = here::here("data/inputs/data_fil_biome4.rds"))

## Biome 5 ----
# Temperate Conifer Forests Forest
data_unm_biome5 <- data_unm |>
  filter(biomeID == 5)
plot_stl(data_unm_biome5)

# filter data with upper quantile
data_fil_biome5 <- data_filter_fc(data_unm_biome5)
plot_stl(data_fil_biome5)
saveRDS(data_fil_biome5, file = here::here("data/inputs/data_fil_biome5.rds"))

## Biome 6 ----
# Boreal Forests/Taiga Forest
data_unm_biome6 <- data_unm |>
  filter(biomeID == 6)
plot_stl(data_unm_biome6)

# filter data with upper quantile
data_fil_biome6 <- data_filter_fc(data_unm_biome6)
plot_stl(data_fil_biome6)
saveRDS(data_fil_biome6, file = here::here("data/inputs/data_fil_biome6.rds"))

## Biome 12 ----
# Mediterranean Forests
data_unm_biome12 <- data_unm |>
  filter(biomeID == 12)
plot_stl(data_unm_biome12)

# filter data with upper quantile
data_fil_biome12 <- data_filter_fc(data_unm_biome12)
plot_stl(data_fil_biome12)
saveRDS(data_fil_biome12, file = here::here("data/inputs/data_fil_biome12.rds"))

# load data
data_fil_biome1 <- readRDS(here::here("data/inputs/data_fil_biome1.rds"))
data_fil_biome2 <- readRDS(here::here("data/inputs/data_fil_biome2.rds"))
data_fil_biome4 <- readRDS(here::here("data/inputs/data_fil_biome4.rds"))
data_fil_biome5 <- readRDS(here::here("data/inputs/data_fil_biome5.rds"))
data_fil_biome6 <- readRDS(here::here("data/inputs/data_fil_biome6.rds"))
data_fil_biome12 <- readRDS(here::here("data/inputs/data_fil_biome12.rds"))

# Filtered ----
# join filtered datasets
data_fil_biomes <- bind_rows(
  data_fil_biome1,
  data_fil_biome2,
  data_fil_biome4,
  data_fil_biome5,
  data_fil_biome6,
  data_fil_biome12
)

saveRDS(data_fil_biomes, file = here::here("data/inputs/data_fil_biomes.rds"))

data_fil_biomes |>
  drop_na(years_since_management) |>
  distinct(dataset)

mean(data_fil_biomes$years_since_management, na.rm = T)

# Check DBH distributions
ggplot() +
  geom_density(data = data_fil_biomes, aes(QMD, ..density.., col = factor(biome), fill = factor(biome)), alpha = 0.7) +
  theme_classic() +
  scale_fill_manual(name = "Biome", values = c(
    "Tropical & Subtropical Moist Broadleaf Forests" = "#E69F00",
    "Temperate Broadleaf & Mixed Forests" = "#56B4E9",
    "Temperate Conifer Forests" = "#009E73",
    "Boreal Forests/Taiga" = "#F5C710",
    "Mediterranean Forests, Woodlands & Scrub" = "#0072B2"
  )) +
  scale_color_manual(name = "Biome", values = c(
    "Tropical & Subtropical Moist Broadleaf Forests" = "#E69F00",
    "Temperate Broadleaf & Mixed Forests" = "#56B4E9",
    "Temperate Conifer Forests" = "#009E73",
    "Boreal Forests/Taiga" = "#F5C710",
    "Mediterranean Forests, Woodlands & Scrub" = "#0072B2"
  )) +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 80))

# Table S1 ----
data <- readRDS(here::here("data/inputs/data_unm.rds"))
data <- readRDS(here::here("data/inputs/data_fil_biomes.rds"))
summary(data)

data |>
  filter(years_since_management != 999) |>
  summarise(mean = mean(years_since_management))
mean(data$plotsize, na.rm = T)
species <- as.data.frame(sort(unique(data_unm$species)))
length(unique(data_unm$plotID))
length(unique(data_unm$species))
summary(data_unm)
table_s1 <- data %>%
  group_by(country, dataset) %>%
  summarise(
    n_plots = n_distinct(plotID),
    min_lon = min(lon),
    max_lon = max(lon),
    min_lat = min(lat),
    max_lat = max(lat),
    min_year = min(year),
    max_year = max(year),
    timespan = max_year - min_year,
    n_census_mean = as.integer(mean(n_census, na.rm = T)),
    interval = as.integer(mean(period, na.rm = T)),
    interval_sd = sd(period, na.rm = T),
    qmd_mean = mean(QMD, na.rm = T),
    QMD_sd = sd(QMD, na.rm = T),
    density_mean = mean(density, na.rm = T),
    density_sd = sd(density, na.rm = T),
    biomass_mean = mean(biomass, na.rm = T),
    biomass_sd = sd(biomass, na.rm = T)
  )
summary(table_s1)
write.csv(table_s1, file = here::here("data/outputs/tableS1.csv"))

# Table S2 ----
table_s2 <- data %>%
  group_by(biome) %>%
  summarise(
    n_plots = n_distinct(plotID),
    min_year = min(year),
    max_year = max(year),
    timespan = mean(period, na.rm = T),
    timespan_sd = sd(period, na.rm = T),
    qmd_mean = mean(QMD, na.rm = T),
    QMD_sd = sd(QMD, na.rm = T),
    density_mean = mean(density, na.rm = T),
    density_sd = sd(density, na.rm = T)
  )
write.csv(table_s2, file = here::here("data/outputs/tableS2.csv"))
