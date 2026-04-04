# This script reads raw data from the original forest data files and process them to create the input data.

# To avoid subjective decisions regarding management and forest use history, 
# we propose classifying each dataset into the following categories based on the information available:

# Variable management_cat:

# 1. Recorded history
# The time since the last management before the first census is known and recorded.
# the number of years since the last management before the first census, or
# the calendar year of the last management intervention is indicated.

# 2. Pristine/primary/old-growth/protected
# The time since the last management before the first census is unknown. The forest was considered pristine, 
# primary, old-growth, protected, or minimally used at the time of the first census.

# 3. Unrecorded history
# The time since the last management before the first census is unknown. The forest may have a history 
# of substantial past wood harvesting and/or visible signs of forest use before the first census.

# load packages ----
# library(renv)
library(readr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)
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
library(BIOMASS)
library(rbeni)
library(ingestr)
library(lubridate)

# load functions ----
source(file.path(here::here(), "/R/functions.R"))

# NFI Spain ----
# Data providers: Paloma Ruiz-Benito and Veronica Cruz-Alonso

# Management: 
# Original info: Only plots were management was not registered during the inventory were provided.
# Management category assigned to 3 for all plots.

# stand-level for species
nfi_spain <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_spa/nfi_spa.csv")
#nfi_spain <- read.csv("~/data/nfi_spa/nfi_spa.csv")

# rename variables
nfi_spain <- nfi_spain |>
  rename(
    plotID = IDPC234,
    lon = CX_ETRS89,
    lat = CY_ETRS89,
    census = CensusID,
    year = Year,
    density = dens,
    species = Nombre234,
    biomass = bioa
  ) |>
  mutate(plotsize = 0.20)

# aggregate data from species to stand data
df_nfi_spain <- from_species_data(nfi_spain) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year) |>
  mutate(
    management = 0,
    management_since_census1_yrs = NA,
    management_cat = 3,
    country = "Spain"
  )

# add coords and biomes
df_nfi_spain <- biomes_coords_utm(df_nfi_spain) |>
  as_tibble()

# add coords and aridity index
df_nfi_spain <- ai_coords_latlon(df_nfi_spain)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_nfi_spain <- lai_coords_latlon(df_nfi_spain)

# add coords and N deposition (Lamarque 2011)
df_nfi_spain <- ndep_coords_latlon(df_nfi_spain)

# add coords and C:N ratio (ISRIC WISE)
df_nfi_spain <- cn_coords_latlon(df_nfi_spain)

# add coords and Phosphorus P - Bray (PBR)
df_nfi_spain <- phos_coords_latlon(df_nfi_spain)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_nfi_spain <- orgc_coords_latlon(df_nfi_spain)

ggplot() +
  geom_point(data = df_nfi_spain, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_nfi_spain, size = 0.5, alpha = 0.5)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = ai), data = df_nfi_spain, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_nfi_spain, file = file.path(here::here(), "/data/inputs/df_nfi_spa.rds"))

# NFI Sweeden ----
# Data providers: Julian Tijerin-Triviño

# Management: 
# Original info: Plots include the management variable (0-unmanaged, 1-managed).
# Management category assigned to 3 for all plots.
# Other variable that could be used: stand age. Plots with stand age > 80 years could indeed be considered as 
# potentially unmanaged or at least without recent major interventions (Management category 2). 
# This approach has limitations—some old forests may have experienced selective cuts or other human disturbances 
# in the past that are not reflected in the current stand age.

# Tree-level for species
species_code <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_swe/swe_sp_code.csv", sep = ",")
#species_code <- read.csv("~/data/nfi_swe/swe_sp_code.csv")

species_code <- species_code |>
  filter(country == "ES") |>
  rename(species = acceptedname) |>
  select(code, species)

nfi_sweeden_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_swe/nfi_swe_tree.csv", sep = ",")
#nfi_sweeden_tree <- read.csv("~/data/nfi_swe/nfi_swe_tree.csv", sep = ",")

nfi_sweeden_tree <- nfi_sweeden_tree |>
  select(-c(speciescode2, speciescode3)) |>
  rename(
    code = speciescode1,
    plotID = plotcode
  ) |>
  pivot_longer(
    cols = c("ba_ha1", "ba_ha2", "ba_ha3"),
    names_to = "census",
    values_to = "ba_tree"
  ) |>
  mutate(census = parse_number(census)) |>
  left_join(species_code)

# Stand-level
nfi_sweeden <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_swe/nfi_swe_stand.csv", sep = ",")
#nfi_sweeden <- read.csv("~/data/nfi_swe/nfi_swe_stand.csv", sep = ",")

# rename variable
nfi_sweeden <- nfi_sweeden |>
  rename(
    plotID = plotcode,
    ba = basal_area
  ) |>
  mutate(
    dbh = NA,
    plotID = parse_number(plotID)
  ) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year) |>
  mutate(
    plotsize = NA,
    biomass = NA,
    management_since_census1_yrs = NA,
    management_cat = 3,
    country = "Sweeden"
  )

# aggregate data from stand to stand data
df_nfi_sweeden <- from_stand_tree_data(nfi_sweeden_tree, nfi_sweeden)

# add coords and biomes
df_nfi_sweeden <- biomes_coords_latlon(df_nfi_sweeden)

# add coords and aridity index
df_nfi_sweeden <- ai_coords_latlon(df_nfi_sweeden)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_nfi_sweeden <- lai_coords_latlon(df_nfi_sweeden)

# add coords and N deposition (Lamarque 2011)
df_nfi_sweeden <- ndep_coords_latlon(df_nfi_sweeden)

# add coords and C:N ratio (ISRIC WISE)
df_nfi_sweeden <- cn_coords_latlon(df_nfi_sweeden)

# add coords and Phosphorus P - Bray (PBR)
df_nfi_sweeden <- phos_coords_latlon(df_nfi_sweeden)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_nfi_sweeden <- orgc_coords_latlon(df_nfi_sweeden)

ggplot() +
  geom_point(data = df_nfi_sweeden, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_nfi_sweeden, size = 0.5, alpha = 0.5)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = ai), data = df_nfi_sweeden, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_nfi_sweeden, file = file.path(here::here(), "/data/inputs/df_nfi_swe.rds"))

# FIA US ----
# Data providers: rFIA

# Tree-level data
filn <- file.path(here::here(), "data/inputs/data_fia_us.rds")

if (!file.exists(filn)) {
  # Download FIA data ---
  # for all states
  # the dataset needed: COND, PLOT, TREE
  states <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs/states.csv")
  st <- states$State.abbreviation

  # Data unavailable for:  DC, MH
  st <- st[-which(st %in% c("DC", "MH"))]

  for (i in st) {
    getFIA(states = i, dir = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs", tables = "COND", load = FALSE)
    getFIA(states = i, dir = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs", tables = "PLOT", load = FALSE)
    options(timeout = 3600)
    getFIA(states = i, dir = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs", tables = "TREE", load = FALSE)
  }

  # Read data ---
  # UNITCD Survey unit code
  # STATECD State code
  # COUNTYCD County code
  # PLOT Plot number

  ## PLOT table ---
  # meta info for forest plots
  setwd("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs")

  data_plot <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs", pattern = "*_PLOT.csv") |>
    purrr::map(read.csv) |>
    lapply(\(x) mutate(x, across(ECO_UNIT_PNW, as.character))) |>
    bind_rows() |>
    # Make a unique ID for each plot, irrespective of time
    mutate(plotID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = "_"))

  ## COND table ---
  # used for filtering unmanaged forest plots (reserves)
  data_cond <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs", pattern = "*_COND.csv") |>
    purrr::map(read.csv) |>
    # lapply(read_csv) |>
    lapply(\(x) mutate(x, across(HABTYPCD1, as.character))) |>
    lapply(\(x) mutate(x, across(HABTYPCD2, as.character))) |>
    lapply(\(x) mutate(x, across(HABTYPCD1_DESCR_PUB_CD, as.character))) |>
    lapply(\(x) mutate(x, across(HABTYPCD2_DESCR_PUB_CD, as.character))) |>
    lapply(\(x) mutate(x, across(HABTYPCD1_PUB_CD, as.character))) |>
    lapply(\(x) mutate(x, across(HABTYPCD2_PUB_CD, as.character))) |>
    bind_rows() |>
    # Make a unique ID for each plot, irrespective of time
    mutate(plotID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = "_"))

  ## TREE table ---
  # Aggregate from tree to forest plot level
  # Note: R session crushes due to the big size of the files. So, we read the files and save only the summaries at plot level.

  # plot size = 1 acre = 0.4 ha
  # DRYBIO_AG - aboveground biomass (DRYBIO_AG) contained in each tree (in pounds, libs)
  # TPA_UNADJ - trees per acre each tree represents
  # DIA - DBH (inches)
  # 1 lb = 0.453 kg
  # 1 acre = 0.405 ha
  # 1 inch = 2.54 cm
  # 1 sq inch = 0.00064516 sq meter
  # abg_biomass: from pounds per acre to kg per ha = *0.453/0.405
  # density: from indiv per acre to indiv per ha = */0.405
  # dbh: from inches per acre to cm per ha = *2.54/0.405
  # BA: from sq. inches per acre to sq. m per ha = *0.00064516/0.405
  species <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs/species_code.csv")
  states <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs/states.csv")
  st <- states$State.abbreviation
  # Data unavailable for:  DC, MH
  st <- st[-which(st %in% c("DC", "MH"))]

  data_stand <- data.frame()
  for (i in st) {
    currentDF <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/obs/", i, "_TREE.csv")
    currentDF <- currentDF |>
      left_join(species) |>
      relocate(species, .after = SPCD) |>
      mutate(plotID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = "_")) |>
      # create census
      group_by(plotID) |>
      mutate(census = match(INVYR, unique(INVYR))) |>
      ungroup() |>
      relocate(census, .after = INVYR)

    dom_species <- currentDF |>
      mutate(plotID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = "_")) |>
      group_by(plotID, INVYR, species) |>
      summarize(ba = sum((pi * DIA * DIA / 4) * TPA_UNADJ * 0.00064516 / 0.405, na.rm = TRUE)) |>
      slice_max(ba, with_ties = FALSE) |>
      ungroup() |>
      select(plotID, INVYR, species)

    currentDF <- currentDF |>
      # Filter trees alive
      filter(STATUSCD == 1) |>
      # aggregate data
      group_by(plotID, INVYR, census) |>
      summarize(
        abg_biomass_kg_ha = sum(DRYBIO_AG * TPA_UNADJ * 0.453 / 0.405, na.rm = TRUE),
        density = sum(TPA_UNADJ / 0.405, na.rm = TRUE),
        dbh = mean(DIA * 2.54, na.rm = TRUE),
        ba = sum((pi * DIA * DIA / 4) * TPA_UNADJ * 0.00064516 / 0.405, na.rm = TRUE)
      ) |>
      ungroup() |>
      mutate(
        QMD = sqrt(ba / (0.0000785 * density)),
        logDensity = log(density),
        logQMD = log(QMD),
        dataset = "fia_us"
      ) |>
      # calculate period lengths
      arrange(plotID, INVYR) |>
      group_by(plotID) |>
      mutate(period = INVYR - lag(INVYR)) |>
      relocate(period, .after = INVYR) |>
      # calculate basal area increment
      mutate(ba_inc = (ba - lag(ba)) / period) |>
      relocate(ba_inc, .after = ba) |>
      # calculate number of censuses
      mutate(n_census = n()) |>
      ungroup() |>
      filter(INVYR != 9999) |>
      # join dominant species to data
      left_join(dom_species) |>
      ungroup()

    data_stand <- rbind(data_stand, currentDF)
  }
  saveRDS(data_stand, file = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs/data_stand_us.rds")

  data_stand <- readRDS("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fia_us/obs/data_stand_us.rds")

  # We want to filter the unmanaged plots. For that we select those plots classified as Reserves.
  # COND table RESERVCD==1 represents the reserves, where no interventions have been carried out.
  data_cond_sel <- data_cond |>
    select(plotID, UNITCD, STATECD, COUNTYCD, PLOT, RESERVCD) |>
    distinct(plotID, .keep_all = T)

  data_plot_sel <- data_plot |>
    select(plotID, LAT, LON, ELEV) |>
    distinct(plotID, .keep_all = T)

  # Join tables
  data_fia_us <- data_stand |>
    left_join(data_cond_sel) |>
    left_join(data_plot_sel) |>
    # rename variable
    rename(
      year = INVYR,
      lat = LAT,
      lon = LON,
      biomass = abg_biomass_kg_ha
    ) |>
    # convert RESERVCD=1 to management=0 = "no managed"
    mutate(
      management = RESERVCD - 1,
      plotsize = 0.4,
      country = "USA",
      years_since_management = NA
    ) |>
    # Remove entries with density = 0
    filter(density > 0) |>
    # Remove entries with QMD = 0
    filter(QMD > 0)

  saveRDS(data_fia_us, file = file.path(here::here(), "/data/inputs/data_fia_us.rds"))
} else {
  data_fia_us <- readRDS(file.path(here::here(), "/data/inputs/data_fia_us.rds"))
}

# add coords and biomes
data_fia_us <- biomes_coords_latlon(data_fia_us)

# add coords and aridity index
data_fia_us <- ai_coords_latlon(data_fia_us)

# add coords and LAI Modis or NDVI (0.5 degrees)
data_fia_us <- lai_coords_latlon(data_fia_us)

# add coords and N deposition (Lamarque 2011)
data_fia_us <- ndep_coords_latlon(data_fia_us)

# add coords and C:N ratio (ISRIC WISE)
data_fia_us <- cn_coords_latlon(data_fia_us)

# add coords and Phosphorus P - Bray (PBR)
data_fia_us <- phos_coords_latlon(data_fia_us)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
data_fia_us <- orgc_coords_latlon(data_fia_us)

ggplot() +
  geom_point(data = data_fia_us, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = data_fia_us, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(data_fia_us, file = file.path(here::here(), "/data/inputs/data_fia_us.rds"))

# NFI Switzerland ----
# Data providers: Brigitte Rohner
# Management: 
# Original info: Each plots provides the variable LETZTENU, which provides information about
# the number of years since last management before each census (LFI1 to LFI5).
# LETZTENU: Number of years since the last silvicultural intervention according to the oral survey at the Forest Service (derivation LETZTENU in the Swiss NFI database). 
# If it is known that no intervention has been carried out on this sample area for a very long time or
# that it is almost certain that no intervention has ever been carried out, LETZTENU=999. 
# Calculated variables:
# - management_since_census1_yrs, sum of LETZTENU years, only when value increase over censuses. If values are not strictly increasing → NA.
# - management = 1 if management_since_census1_yrs=NA, and management = 0 otherwise.
# - management_cat = 1 for all plots. (Recorded history)

# Plot characteristics
lfi_plot_constant <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_swi/lfi_plot_constant.csv", sep = ";")
#lfi_plot_constant <- read.csv("~/data/nfi_swi/lfi_plot_constant.csv", sep = ";")

# Species names
lfi_species_names <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_swi/lfi_species_names.csv", sep = ",")
#lfi_species_names <- read.csv("~/data/nfi_swi/lfi_species_names.csv", sep = ",")

## calculate dom_species from tree-level data
lfi_tree_census <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_swi/lfi_tree_census.csv", sep = ",")
#lfi_tree_census <- read.csv("~/data/nfi_swi/lfi_tree_census.csv", sep = ",")

dom_species <- lfi_tree_census |>
  mutate(year = str_sub(CENSUS_DATE, 7, 10)) |>
  mutate(year = as.numeric(year)) |>
  relocate(year, .after = CENSUS_DATE) |>
  left_join(lfi_species_names) |>
  rename(species = SPECIES_NAME) |>
  group_by(PLOTID, year, species) |>
  add_tally() |>
  mutate(density = REPRESENTATION * n) |>
  mutate(ba_tree = pi * (DBH * 0.01)^2 / 4 * density) |>
  summarise(ba = sum(ba_tree, na.rm = T)) |>
  slice_max(ba, with_ties = FALSE) |>
  ungroup() |>
  select(PLOTID, year, species)

# Stand-level data
lfi_plot_census <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_swi/lfi_plot_census.csv", sep = ",")
#lfi_plot_census <- read.csv("~/data/nfi_swi/lfi_plot_census.csv", sep = ",")

# prepare data
nfi_switzerland <- lfi_plot_census |>
  left_join(lfi_plot_constant[, c(1:3)]) |>
  # create Year variable from CENSUS_DATE
  mutate(year = year(as.Date(CENSUS_DATE, format="%d.%m.%Y"))) |>
  left_join(dom_species) |>
  # rename variable
  rename(
    plotID = PLOTID,
    lon = LONGITUDE,
    lat = LATITUDE,
    census = CENSUSID,
    ba = BASAL_AREA_HA,
    density = NPH,
    dbh = MEAN_DBH_HA,
    biomass = BIOMASS_VPPS_ABOVEGROUND
  ) |>
  group_by(plotID) |>
  mutate(
    management_since_census1_yrs =
      if (any(LETZTENU == 999, na.rm = TRUE)) 999        # unmanaged
    else if (any(is.na(LETZTENU))) NA_real_              # missing data
    else if (all(diff(LETZTENU) > 0)) sum(LETZTENU)      # strictly increasing
    else NA_real_                                        # management occurred
  ) |>
  ungroup() |>
  mutate(
    management = ifelse(is.na(management_since_census1_yrs), 1, 0),
    management_cat = 1, 
    plotsize = PLOT_AREA_LARGE * 10^-4,
    country = "Switzerland")

# aggregate data from stand to stand data
df_nfi_switzerland <- from_stand_data(nfi_switzerland)

# add coords and biomes
df_nfi_switzerland <- biomes_coords_latlon(df_nfi_switzerland)

# add coords and aridity index
df_nfi_switzerland <- ai_coords_latlon(df_nfi_switzerland)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_nfi_switzerland <- lai_coords_latlon(df_nfi_switzerland)

# add coords and N deposition (Lamarque 2011)
df_nfi_switzerland <- ndep_coords_latlon(df_nfi_switzerland)

# add coords and C:N ratio (ISRIC WISE)
df_nfi_switzerland <- cn_coords_latlon(df_nfi_switzerland)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_nfi_switzerland <- phos_coords_latlon(df_nfi_switzerland)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_nfi_switzerland <- orgc_coords_latlon(df_nfi_switzerland)

ggplot() +
  geom_point(data = df_nfi_switzerland, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_nfi_switzerland, size = 0.5, alpha = 0.5) + theme(legend.position = "bottom")

# Save stand-level data
saveRDS(df_nfi_switzerland, file = file.path(here::here(), "/data/inputs/df_nfi_swi.rds"))

# NFI Norway ----
# Data providers: Oliver Moen Snoksrud and Johannes Breidenbach
# Management:
# Original info: 
# - management_since_census1_yrs = 64,
# - management_cat = 1 for all plots. (Recorded history)

# Stand-level data
nfi_norway <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfi_nor/nfi_nor.csv", sep = ",")
#nfi_norway <- read.csv("~/data/nfi_nor/nfi_nor.csv", sep = ",")

# rename variable
nfi_norway <- nfi_norway |>
  filter(class == "F_Forest") |>
  mutate( # convert to cm
    dbh = dbh / 10,
    qmd = qmd / 10
  ) |> 
  rename(
    ba = basal_area,
    plotsize = area
  ) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year) |>
  mutate(
    management = 0,
    management_cat = 1, 
    management_since_census1_yrs = 64,
    country = "Norway"
  )

# aggregate data from stand
df_nfi_norway <- from_stand_data(nfi_norway) 

# add coords and biomes
df_nfi_norway <- biomes_coords_latlon(df_nfi_norway)

# add coords and aridity index
df_nfi_norway <- ai_coords_latlon(df_nfi_norway)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_nfi_norway <- lai_coords_latlon(df_nfi_norway)

# add coords and N deposition (Lamarque 2011)
df_nfi_norway <- ndep_coords_latlon(df_nfi_norway)

# add coords and C:N ratio (ISRIC WISE)
df_nfi_norway <- cn_coords_latlon(df_nfi_norway)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_nfi_norway <- phos_coords_latlon(df_nfi_norway)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_nfi_norway <- orgc_coords_latlon(df_nfi_norway)

ggplot() +
  geom_point(data = df_nfi_norway, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_nfi_norway, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_nfi_norway, file = file.path(here::here(), "/data/inputs/df_nfi_nor.rds"))

# EFM Switzerland ----
# David Forrester and Jonas Glatthorn
# Management:
# Original: The dataset efm_management includes the year_last_intervention variable
# Calculate: management_since_census1_yrs from year_last_intervention and first_measurement

# Plot characteristics
efm_metadata <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/efm/raw/VFL_LISTmanual.csv")
#efm_metadata <- read.csv("~/data/efm/raw/VFL_LISTmanual.csv")
efm_metadata <- efm_metadata |>
  mutate(FNUM = as.character(FNUM)) |>
  select(-BA)

# Plot area
efm_area <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/efm/raw/EFM_plot_area.csv")
#efm_area <- read.csv("~/data/efm/raw/EFM_plot_area.csv")
efm_area <- efm_area |>
  mutate(FNUM = as.character(FNUM))

efm_locations <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/efm/raw/efm_plot_locations.csv")
#efm_locations <- read.csv("~/data/efm/raw/efm_plot_locations.csv")
efm_locations <- efm_locations |>
  mutate(FNUM = as.character(FNUM))

# Last management intervention
efm_management <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/efm/raw/EFM_last_intervention.csv")
#efm_management <- read.csv("~/data/efm/raw/EFM_last_intervention.csv")
efm_management <- efm_management |>
  mutate(FNUM = as.character(FNUM),
         management_since_census1_yrs = first_measurement - year_last_intervention) |>
  select(FNUM, management_since_census1_yrs)
  
# Stand level data (per plot, year and species)
efm_stand <- readRDS("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/efm/raw/efm_stand.RDS") 
#efm_stand <- readRDS("~/data/efm/raw/efm_stand.RDS")
efm_stand <- efm_stand |>
  mutate(FNUM = as.character(FNUM)) |>
  select(-BA)

# Select stand level data for all species: "All species combined"
dom_species_efm <- efm_stand |>
  filter(Latin != "All species combined") |>
  group_by(FNUM, AJ) |>
  top_n(1, BasalAreaAHC1_2_m2perha) |>
  rename(species = Latin)

efm_swi <- efm_stand |>
  filter(Latin == "All species combined") |>
  left_join(efm_metadata) |>
  left_join(efm_area) |>
  left_join(efm_locations) |>
  left_join(efm_management) |>
  left_join(dom_species_efm[, c(1:3)]) |>
  rename(
    plotID = FNUM,
    year = AJ,
    ba = BasalAreaAHC1_2_m2perha,
    dbh = DBHqAHC1_2_cm,
    plotsize = PlotArea_ha,
    density = TreesPerHectareAHC1_2
  ) |>
  mutate(
    qmd = sqrt(ba / (0.0000785 * density)),
    biomass = (AbovegroundAHC1_2_Mgperha + RootmassAHC1_2_Mgperha) * 1000
  ) |>
  filter(density != 0) |>
  filter(is.na(density) == FALSE) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year) |>
  # create management variable
  mutate(
    management = 0,
    management_cat = 1, 
    country = "Switzerland"
  )

# aggregate data from stand
df_efm_swi <- from_stand_data(efm_swi)

# add coords and biomes
df_efm_swi <- biomes_coords_latlon(df_efm_swi)

# add coords and aridity index
df_efm_swi <- ai_coords_latlon(df_efm_swi)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_efm_swi <- lai_coords_latlon(df_efm_swi)

# add coords and N deposition (Lamarque 2011)
df_efm_swi <- ndep_coords_latlon(df_efm_swi)

# add coords and C:N ratio (ISRIC WISE)
df_efm_swi <- cn_coords_latlon(df_efm_swi)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_efm_swi <- phos_coords_latlon(df_efm_swi)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_efm_swi <- orgc_coords_latlon(df_efm_swi)

ggplot() +
  geom_point(data = df_efm_swi, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_efm_swi, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_efm_swi, file = file.path(here::here(), "/data/inputs/df_efm_swi.rds"))

# uholka ----
# Data providers: Jonas Stillhard
# status 1 Tree found (1) or not found (0) during inventory
# status 2 Tree alive (1) or dead (0)
# status 3 Tree standing (1) or lying (0) Boolean
# status 4 Only for dead trees: entire tree (1) or broken tree (snag), (0)

# Management: 
# Original info: Plots are considered to be primary in the sense that hey have not been managed in the past. 
# Derived info: Since history is known, we assign management_cat = 1 and management_since_census1_yrs = 999.

# Tree-level data
uholka <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/uholka/uholka_tree.csv", sep = ",")
uholka <- read.csv("~/data/uholka/uholka_tree.csv", sep = ",")

# prepare data
uholka <- uholka |>
  filter(status_2 == 1) |> # filter alive trees
  rename(
    x = x_local,
    y = y_local,
    treeID = tree_nr
  ) |>
  # calculate dbh and basal area
  mutate(
    dbh = (dbh_1 + dbh_2) / 2 / 10,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotID = as.numeric(ifelse(nchar(treeID) == 4, substr(treeID, 1, 1), substr(treeID, 1, 2))),
    plotsize = 0.25, # ha
    biomass = NA,
    country = "Ukraine",
    management_since_census1_yrs = 999
  ) |>
  drop_na(dbh) 

# aggregate data from stand to stand data
df_uholka <- from_tree_data(uholka) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year) |>
  mutate(management = 0,
         management_cat = 1)

# add coords and biomes
df_uholka <- biomes_coords_latlon(df_uholka)

# add coords and aridity index
df_uholka <- ai_coords_latlon(df_uholka)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_uholka <- lai_coords_latlon(df_uholka)

# add coords and N deposition (Lamarque 2011)
df_uholka <- ndep_coords_latlon(df_uholka)

# add coords and C:N ratio (ISRIC WISE)
df_uholka <- cn_coords_latlon(df_uholka)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_uholka <- phos_coords_latlon(df_uholka)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_uholka <- orgc_coords_latlon(df_uholka)

ggplot() +
  geom_point(data = df_uholka, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_uholka, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_uholka, file = file.path(here::here(), "/data/inputs/df_uholka.rds"))

# fep gre ----
# Greece
# Data providers: Gavriil Spyroglou and Nikolaos Fyllas
# Management: 
# The dataset provides year_last_intervention and year of first measurement (plot establishment)
# Calculated: management_since_census1_yrs = first_measurement - year_last_intervention.
# Note: Data also includes years_since_management. We did not use this variable.

# tree-level data to estimate dominant species
fep_gre_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fp_gre/fp_gre_tree.csv", sep = ",")
fep_gre_tree <- read.csv("~/data/fp_gre/fp_gre_tree.csv", sep = ",")

dom_species <- fep_gre_tree |>
  filter(status == "Alive") |> # filter alive trees
  rename(
    year = census_year,
    plotID = plot_id
  ) |>
  group_by(plotID, year, species) |>
  add_tally() |>
  mutate(density = n,
         ba_tree = pi * (dbh * 0.01)^2 / 4 * density) |>
  summarise(ba = sum(ba_tree, na.rm = T)) |>
  slice_max(ba, with_ties = FALSE) |>
  ungroup() |>
  select(plotID, year, species)

# stand-level data
fep_gre_stand <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fp_gre/fp_gre_stand.csv", sep = ",")
fep_gre_stand <- read.csv("~/data/fp_gre/fp_gre_stand.csv", sep = ",")

# rename variable
fep_gre_stand <- fep_gre_stand |>
  rename(
    plotID = plot_id,
    lon = Longtidute,
    lat = Latitude,
    year = census_year,
    ba = Basal_area_m2_ha.1,
    density = N_ha,
    dbh = dbh_m,
    plotsize = plot_size
  ) |>
  mutate(biomass = biomass.t.ha.1 * 10^3) |>
  left_join(dom_species) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year)

# aggregate data from stand
df_fep_gre <- from_stand_data(fep_gre_stand) |>
  mutate(
    management = 0,
    management_cat = 1,
    management_since_census1_yrs = first_measurement - year_last_intervention,
    country = "Greece"
  )

# add coords and biomes
df_fep_gre <- biomes_coords_latlon(df_fep_gre)

# add coords and aridity index
df_fep_gre <- ai_coords_latlon(df_fep_gre)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_fep_gre <- lai_coords_latlon(df_fep_gre)

# add coords and N deposition (Lamarque 2011)
df_fep_gre <- ndep_coords_latlon(df_fep_gre)

# add coords and C:N ratio (ISRIC WISE)
df_fep_gre <- cn_coords_latlon(df_fep_gre)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_fep_gre <- phos_coords_latlon(df_fep_gre)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_fep_gre <- orgc_coords_latlon(df_fep_gre)

ggplot() +
  geom_point(data = df_fep_gre, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_fep_gre, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_fep_gre, file = file.path(here::here(), "/data/inputs/df_fep_gre.rds"))

# INRAE LESSEM ----
# France
# Data providers: Georges Kunstler
# Management:
# Original: no info on management in the table.
# New information about management included in fp_fra_stand: year_last_intervention and year of first measurement (plot establishment)
# Calculated: management_since_census1_yrs = first_measurement - year_last_intervention.
# Note: Data also includes years_since_management. We did not use this variable.

# Tree- and -stand level data
inrae_lessem_plot <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fp_fra/fp_fra_stand.csv", sep = ",") |> as_tibble()
inrae_lessem_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fp_fra/fp_fra_tree.csv", sep = ",") |> as_tibble()
inrae_lessem_species <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fp_fra/fra_sp_code.csv", sep = ",") |> as_tibble()
inrae_lessem_status <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/fp_fra/fra_sta_code.csv", sep = ",") |> as_tibble()

inrae_lessem_plot <- read.csv("~/data/fp_fra/fp_fra_stand.csv", sep = ",")
inrae_lessem_tree <- read.csv("~/data/fp_fra/fp_fra_tree.csv", sep = ",")
inrae_lessem_species <- read.csv("~/data/fp_fra/fra_sp_code.csv", sep = ",")
inrae_lessem_status <- read.csv("~/data/fp_fra/fra_sta_code.csv", sep = ",")

# prepare data
inrae_lessem <- inrae_lessem_tree |>
  left_join(inrae_lessem_plot, by = "plot_id") |>
  left_join(inrae_lessem_species) |>
  left_join(inrae_lessem_status) |>
  # filter alive trees
  filter(status == "alive") |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    biomass = NA,
    country = "France"
  ) |>
  rename(
    plotsize = area_ha,
    lon = long,
    species = Latin.name,
    plotID = plot_id
  )

# aggregate data from stand to stand data
df_inrae_lessem <- from_tree_data(inrae_lessem) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
  management = 0,
  management_cat = 1
)

# add coords and biomes
df_inrae_lessem <- biomes_coords_latlon(df_inrae_lessem)

# add coords and aridity index
df_inrae_lessem <- ai_coords_latlon(df_inrae_lessem)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_inrae_lessem <- lai_coords_latlon(df_inrae_lessem)

# add coords and N deposition (Lamarque 2011)
df_inrae_lessem <- ndep_coords_latlon(df_inrae_lessem)

# add coords and C:N ratio (ISRIC WISE)
df_inrae_lessem <- cn_coords_latlon(df_inrae_lessem)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_inrae_lessem <- phos_coords_latlon(df_inrae_lessem)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_inrae_lessem <- orgc_coords_latlon(df_inrae_lessem)

ggplot() +
  geom_point(data = df_inrae_lessem, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_inrae_lessem, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_inrae_lessem, file = file.path(here::here(), "/data/inputs/df_inrae_lessem.rds"))

# EuFoRia plots ----

## bnp ----
# Data providers: Michael Maroschek and Rupert Seidl
# Management:
# Original: In dataset, year_last_management. In Table S1, last management intervention prior census 1 = 6

# Stand-level data by species
bnp <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/bnp/euf_bnp.csv", sep = ",")
bnp <- read.csv("~/data/euforia/bnp/euf_bnp.csv", sep = ",")

# prepare data
bnp <- bnp |>
  rename(
    lon = longitude_wgs84,
    lat = latitude_wgs84,
    dbh = mean_dbh,
    ba = basal_area,
    species = dominant_species
  ) 

# aggregate data from stand
df_bnp <- from_stand_data(bnp) |>
  mutate(
    biomass = NA,
    plotsize = 0.05,
    management = 0,
    management_cat = 1,
    country = "Germany"
  ) |>
group_by(plotID) %>%
  mutate(
    management_since_census1_yrs = min(year, na.rm = TRUE) - year_last_management
  ) |>
  ungroup()

# add coords and biomes
df_bnp <- biomes_coords_latlon(df_bnp)

# add coords and aridity index
df_bnp <- ai_coords_latlon(df_bnp)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_bnp <- lai_coords_latlon(df_bnp)

# add coords and N deposition (Lamarque 2011)
df_bnp <- ndep_coords_latlon(df_bnp)

# add coords and C:N ratio (ISRIC WISE)
df_bnp <- cn_coords_latlon(df_bnp)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_bnp <- phos_coords_latlon(df_bnp)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_bnp <- orgc_coords_latlon(df_bnp)

ggplot() +
  geom_point(data = df_bnp, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_bnp, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_bnp, file = file.path(here::here(), "/data/inputs/df_bnp.rds"))

## czu ----
# Data providers: Miroslav Svoboda
# Data at tree level not used
# Management:
# Original: "Nature reserve established 1955. Low intensity salvage logging of single trees or group of trees till 1990"
# In Table S1. indicated management_since_census1_yrs = 69. Need to clarify with author.

# Tree- and stand-level data by species
czu_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/czu/euf_czu_tree.csv", sep = ",")
czu_stand <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/czu/euf_czu_stand.csv", sep = ",")

czu_tree <- read.csv("~/data/euforia/czu/euf_czu_tree.csv", sep = ",") 
czu_stand <- read.csv("~/data/euforia/czu/euf_czu_stand.csv", sep = ",")

# prepare data
czu <- czu_stand |>
  rename(
    lon = longitude,
    lat = latitude,
    dbh = mean_dbh_mm,
    ba = basal_area_m2_ha,
    species = dominant_species
  ) |>
  mutate(
    dbh = dbh / 10,
    plotsize = plot_size_m2 / 10000,
    biomass = NA
  )

# aggregate data from stand
df_czu <- from_stand_data(czu) |>
  mutate(
    management = 0,
    management_cat = 1,
    management_since_census1_yrs = 69,
    country = "Czech Republic"
  )

# add coords and biomes
df_czu <- biomes_coords_latlon(df_czu)

# add coords and aridity index
df_czu <- ai_coords_latlon(df_czu)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_czu <- lai_coords_latlon(df_czu)

# add coords and N deposition (Lamarque 2011)
df_czu <- ndep_coords_latlon(df_czu)

# add coords and C:N ratio (ISRIC WISE)
df_czu <- cn_coords_latlon(df_czu)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_czu <- phos_coords_latlon(df_czu)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_czu <- orgc_coords_latlon(df_czu)

ggplot() +
  geom_point(data = df_czu, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_czu, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_czu, file = file.path(here::here(), "/data/inputs/df_czu.rds"))

## fvabw ----
# Data providers: Yannek Käber and Lucia Seebach
# Management:
# Original: In Table S1 indicated that last management was before the time of designation
# Confirm what is tsc variable? 
# Assigned management_cat = 2

# Stand-level data by species
fvabw <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/forst/euf_forst.csv", sep = ",")
fvabw <- read.csv("~/data/euforia/forst/euf_forst.csv", sep = ",")

# prepare data
fvabw <- fvabw |>
  rename(
    plotID = plot_id,
    lon = utm_x,
    lat = utm_y,
    density = stem_count_ha,
    dbh = meanDBH,
    ba = basal_area_ha,
    plotsize = plot_size_ha
  ) |>
  mutate(biomass = NA)

# aggregate data from stand to stand data
df_fvabw <- from_species_data(fvabw) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management = 0,
    management_cat = 2,
    management_since_census1_yrs = NA,
    country = "Germany",
  )

# add coords and biomes
df_fvabw <- biomes_coords_utm(df_fvabw)

# add coords and aridity index
df_fvabw <- ai_coords_latlon(df_fvabw)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_fvabw <- lai_coords_latlon(df_fvabw)

# add coords and N deposition (Lamarque 2011)
df_fvabw <- ndep_coords_latlon(df_fvabw)

# add coords and C:N ratio (ISRIC WISE)
df_fvabw <- cn_coords_latlon(df_fvabw)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_fvabw <- phos_coords_latlon(df_fvabw)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_fvabw <- orgc_coords_latlon(df_fvabw)

ggplot() +
  geom_point(data = df_fvabw, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_fvabw, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_fvabw, file = file.path(here::here(), "/data/inputs/df_fvabw.rds"))

## iberbas ----
# Data providers: Tzvetan Zlatanov
# Management:
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = 70

# Tree-level data
iberbas <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/iberbas/euf_iberbas.csv", sep = ",")
iberbas <- read.csv("~/data/euforia/iberbas/euf_iberbas.csv", sep = ",")

# prepare data
iberbas <- iberbas |>
  filter(dbh != 0) |>
  rename(plotsize = plot_size) |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    biomass = NA,
    management_since_census1_yrs = 70,
    country = "Bulgaria")

# aggregate data from stand to stand data
df_iberbas <- from_tree_data(iberbas) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup()  |>
  mutate(
    management = 0,
    management_cat = 1
  )

# add coords and biomes
df_iberbas <- biomes_coords_latlon(df_iberbas)

# add coords and aridity index
df_iberbas <- ai_coords_latlon(df_iberbas)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_iberbas <- lai_coords_latlon(df_iberbas)

# add coords and N deposition (Lamarque 2011)
df_iberbas <- ndep_coords_latlon(df_iberbas)

# add coords and C:N ratio (ISRIC WISE)
df_iberbas <- cn_coords_latlon(df_iberbas)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_iberbas <- phos_coords_latlon(df_iberbas)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_iberbas <- orgc_coords_latlon(df_iberbas)

ggplot() +
  geom_point(data = df_iberbas, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_iberbas, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_iberbas, file = file.path(here::here(), "/data/inputs/df_iberbas.rds"))

## unitbv ----
# Data providers: Any Mary Petritan, Cătălin Petritan
# Management:
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = NA

# Tree-level data
unitbv <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/incds/euf_incds.csv", sep = ",")
unitbv <- read.csv("~/data/euforia/incds/euf_incds.csv", sep = ",")

unitbv <- unitbv |>
  select(-c(Height2003, Height2013, Height2023)) |>
  pivot_longer(
    cols = c(dbh_2003, dbh_2013, dbh_2023, status_2003, status_2013, status_2023, biomass_2003, biomass_2013, biomass_2023),
    names_to = c(".value", "year"),
    names_sep = "_",
    values_drop_na = FALSE
  )

# Divide plot into grids of different size given the coordinates
# plotsize = 100 * 100 = 1 ha
# New plots = 20x20 = 0.04

unitbv$grid20 <- interaction(
  cut(unitbv$X,
    breaks = seq(0, 100, by = 20),
    include.lowest = TRUE
  ),
  cut(unitbv$Y,
    breaks = seq(0, 100, by = 20),
    include.lowest = TRUE
  ),
  sep = "X"
)
# prepare data
unitbv <- unitbv |>
  rename(plotIDD = plotID) |>
  group_by(grid20) |>
  mutate(plotID = cur_group_id()) |>
  ungroup()
ggplot(unitbv) +
  geom_point(aes(X, Y, col = plotID))
length(unique(unitbv$plotID))

# prepare data
unitbv <- unitbv |>
  filter(status == "alive") |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    year = as.numeric(year),
    plotsize = 0.04,
    country = "Rumania",
    management_since_census1_yrs = NA
  )

# aggregate data from stand to stand data
df_unitbv <- from_tree_data(unitbv) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 2,
    management = 0
  )

# add coords and biomes
df_unitbv <- biomes_coords_latlon(df_unitbv)

# add coords and aridity index
df_unitbv <- ai_coords_latlon(df_unitbv)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_unitbv <- lai_coords_latlon(df_unitbv)

# add coords and N deposition (Lamarque 2011)
df_unitbv <- ndep_coords_latlon(df_unitbv)

# add coords and C:N ratio (ISRIC WISE)
df_unitbv <- cn_coords_latlon(df_unitbv)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_unitbv <- phos_coords_latlon(df_unitbv)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_unitbv <- orgc_coords_latlon(df_unitbv)

ggplot() +
  geom_point(data = df_unitbv, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_unitbv, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_unitbv, file = file.path(here::here(), "/data/inputs/df_unitbv.rds"))

## lwf ----
# Data providers: Markus Blaschke
# Management:
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = NA

# Stand-level data
lwf_stand <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/lwf/euf_lwf_stand.csv", sep = ";")
lwf_stand <- read.csv("~/data/euforia/lwf/euf_lwf_stand.csv", sep = ";")

lwf_stand <- lwf_stand |>
  rename(
    plotID = Plot_ID,
    lon = x_Koord_etrs32,
    lat = y_Koord_etrs32,
    plotsize = Plot_size
  ) |>
  select(plotID, lon, lat, plotsize) |>
  distinct()

# Tree-level data
lwf_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/lwf/euf_lwf_tree.csv", sep = ",")
lwf_tree <- read.csv("~/data/euforia/lwf/euf_lwf_tree.csv", sep = ",")

# prepare data
lwf_tree <- lwf_tree |>
  filter(Status == "alive") |>
  mutate(
    ba_tree = pi * (DBH * 0.01 / 2)^2,
    biomass = NA,
    country = "Germany",
    management_since_census1_yrs = NA
  ) |>
  rename(
    plotID = Plot_ID,
    year = Census,
    species = Species,
    status = Status,
    dbh = DBH
  ) |>
  left_join(lwf_stand) 

# aggregate data from tree to stand data
df_lwf <- from_tree_data(lwf_tree) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 3,
    management = 0
  )

# add coords and biomes
df_lwf <- biomes_coords_utm(df_lwf)

# add coords and aridity index
df_lwf <- ai_coords_latlon(df_lwf)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_lwf <- lai_coords_latlon(df_lwf)

# add coords and N deposition (Lamarque 2011)
df_lwf <- ndep_coords_latlon(df_lwf)

# add coords and C:N ratio (ISRIC WISE)
df_lwf <- cn_coords_latlon(df_lwf)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_lwf <- phos_coords_latlon(df_lwf)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_lwf <- orgc_coords_latlon(df_lwf)

ggplot() +
  geom_point(data = df_lwf, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_lwf, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_lwf, file = file.path(here::here(), "/data/inputs/df_lwf.rds"))

## npvbw ----
# Data from Marco Heurich and Isabelle Klein
# Management:
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = 42

# tree-level data
npvbw <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/nbw/euf_nbw_tree.csv", sep = ",")
npvbw <- read.csv("~/data/euforia/nbw/euf_nbw_tree.csv", sep = ",")

# prepare data
npvbw <- npvbw |>
  filter(status == "alive") |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    biomass = NA,
    country = "Germany",
    #management_since_census1_yrs = 42
  ) |>
  rename(
    plotsize = plot_size,
    year = census
  ) |>
  group_by(plotID) |>
  mutate(
    lat = mean(latitude),
    lon = mean(longitude),
    management_since_census1_yrs = min(year, na.rm = TRUE) - year_last_management
  ) |>
  ungroup()

# aggregate data from stand to stand data
df_npvbw <- from_tree_data(npvbw) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 1,
    management = 0
  )

# add coords and biomes
df_npvbw <- biomes_coords_latlon(df_npvbw)

# add coords and aridity index
df_npvbw <- ai_coords_latlon(df_npvbw)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_npvbw <- lai_coords_latlon(df_npvbw)

# add coords and N deposition (Lamarque 2011)
df_npvbw <- ndep_coords_latlon(df_npvbw)

# add coords and C:N ratio (ISRIC WISE)
df_npvbw <- cn_coords_latlon(df_npvbw)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_npvbw <- phos_coords_latlon(df_npvbw)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_npvbw <- orgc_coords_latlon(df_npvbw)

ggplot() +
  geom_point(data = df_npvbw, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_npvbw, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_npvbw, file = file.path(here::here(), "/data/inputs/df_npvbw.rds"))

## nfr ----
# Martina Hobi and Harald Bugmann
# Management:
# Original dataset includes info on management.

# Metadata
nfr_metadata <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfr/raw/NFR_metadata.csv")
nfr_metadata <- read.csv("~/data/euforia/nfr/raw/NFR_metadata.csv", sep = ",")

nfr_metadata <- nfr_metadata |>
  mutate(fg = as.character(fg)) |>
  select(fg, lat, long, ele, temp, precip) |>
  distinct(fg, .keep_all = TRUE)

# Plot area
nfr_plot_area <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfr/raw/NFR_plot_area.csv")
nfr_plot_area <- read.csv("~/data/euforia/nfr/raw/NFR_plot_area.csv", sep = ",")

nfr_plot_area <- nfr_plot_area |>
  mutate(
    fg = as.character(fg),
    FNUM = as.character(FNUM)
  ) |>
  left_join(nfr_metadata)

# Last management intervention
nfr_last_intervention <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfr/raw/NFR_last_intervention.csv")
nfr_last_intervention <- read.csv("~/data/euforia/nfr/raw/NFR_last_intervention.csv", sep = ",")

nfr_last_intervention <- nfr_last_intervention |>
  mutate(fg = as.character(fg))

# Stand level data (per plot, year and species)
nfr_stand <- readRDS("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/nfr/raw/NFR_stand_data.RDS") # 291 plots From David Forrester data
nfr_stand <- readRDS("~/data/euforia/nfr/raw/NFR_stand_data.RDS")

# Select Stand level data for all species: "All species combined"
dom_species_nfr <- nfr_stand |>
  mutate(FNUM = as.character(FNUM)) |>
  filter(Latin != "All species combined") |>
  group_by(FNUM, AJ) |>
  top_n(1, BasalAreaAHC1_2_m2perha) |>
  rename(species = Latin)

nfr_swi <- nfr_stand |>
  filter(Latin == "All species combined") |>
  mutate(FNUM = as.character(FNUM)) |>
  left_join(nfr_plot_area) |>
  left_join(nfr_last_intervention[, c(1, 5, 6)]) |>
  left_join(dom_species_nfr[, c(1, 4, 5)]) |>
  rename(
    plotID = FNUM,
    year = AJ,
    ba = BasalAreaAHC1_2_m2perha,
    dbh = DBHqAHC1_2_cm,
    lon = long,
    plotsize = PlotArea_ha,
    density = TreesPerHectareAHC1_2
  ) |>
  mutate(
    qmd = sqrt(ba / (0.0000785 * density)),
    biomass = (AbovegroundAHC1_2_Mgperha + RootmassAHC1_2_Mgperha) * 1000
  ) |>
  filter(density != 0) |>
  filter(is.na(density) == FALSE) |>
  # Remove plotID because it has a disproportionately high BA and Biomass increment, as suggested by David Forrester
  filter(plotID != 7007028, plotID != 7005006, plotID != 7005003, plotID != 7005002, plotID != 7001001, plotID != 7001002, plotID != 7001003, plotID != 7002003) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year) |>
  # create management variable
  mutate(
    country = "Switzerland",
    management_since_census1_yrs = first_measurement - year_last_intervention
  )

# aggregate data from stand
df_nfr_swi <- from_stand_data(nfr_swi) |>
  mutate(
    management_cat = 1,
    management = 0
  )

# add coords and biomes
df_nfr_swi <- biomes_coords_latlon(df_nfr_swi)

# add coords and aridity index
df_nfr_swi <- ai_coords_latlon(df_nfr_swi)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_nfr_swi <- lai_coords_latlon(df_nfr_swi)

# add coords and N deposition (Lamarque 2011)
df_nfr_swi <- ndep_coords_latlon(df_nfr_swi)

# add coords and C:N ratio (ISRIC WISE)
df_nfr_swi <- cn_coords_latlon(df_nfr_swi)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_nfr_swi <- phos_coords_latlon(df_nfr_swi)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_nfr_swi <- orgc_coords_latlon(df_nfr_swi)

ggplot() +
  geom_point(data = df_nfr_swi, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_nfr_swi, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_nfr_swi, file = file.path(here::here(), "/data/inputs/df_nfr_swi.rds"))

## nwfva ----
# Data providers: Peter Meyer
# Original dataset includes info on management. Max management_since_census1_yrs = 10
# In Table S1, management_since_census1_yrs = 52. Need to contact data providers.

# Stand-level data
nwfva_stand <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/nwfva/euf_nwfva_stand.csv", sep = ",")
nwfva_stand <- read.csv("~/data/euforia/nwfva/euf_nwfva_stand.csv", sep = ",")

# Tree-level data
nwfva_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/nwfva/euf_nwfva_tree.csv", sep = ",")
nwfva_tree <- read.csv("~/data/euforia/nwfva/euf_nwfva_tree.csv", sep = ",")

# prepare data
nwfva_stand <- nwfva_stand |>
  as_tibble() |>
  unite(plotID, c("reserve_id", "plot_id")) |>
  select(plotID, year, census, year_last_management) |>
  group_by(plotID) |>
  mutate(management_since_census1_yrs = min(year, na.rm = TRUE) - year_last_management)

nwfva_tree <- nwfva_tree |>
  filter(status == "A") |> # alive trees
  as_tibble() |>
  unite(plotID, c("reserve_id", "plot_id")) |>
  left_join(nwfva_stand) |>
  rename(
    lon = long_plot,
    lat = lati_plot,
    plotsize = plot_size,
    treeID = treeID1
  ) |>
  mutate(
    dbh = dbhmm / 10,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    year = as.numeric(year),
    country = "Germany",
    azimuth_rad = (90 - azimuth) * pi / 180,
    x = distance * cos(azimuth),
    y = distance * sin(azimuth),
    biomass = NA
  ) 

# aggregate data from tree to stand data
df_nwfva <- from_tree_data(nwfva_tree) |>
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 1,
    management = 0
  )

# add coords and biomes
df_nwfva <- biomes_coords_utm(df_nwfva)

# add coords and aridity index
df_nwfva <- ai_coords_latlon(df_nwfva)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_nwfva <- lai_coords_latlon(df_nwfva)

# add coords and N deposition (Lamarque 2011)
df_nwfva <- ndep_coords_latlon(df_nwfva)

# add coords and C:N ratio (ISRIC WISE)
df_nwfva <- cn_coords_latlon(df_nwfva)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_nwfva <- phos_coords_latlon(df_nwfva)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_nwfva <- orgc_coords_latlon(df_nwfva)

ggplot() +
  geom_point(data = df_nwfva, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_nwfva, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_nwfva, file = file.path(here::here(), "/data/inputs/df_nwfva.rds"))

## tuzvo ----
# Data providers: Stanislav Kucbel and Peter Jalovia
# Original dataset does not include info on management. 
# In Table S1, management_since_census1_yrs = NA, management_cat = 2

# Stand-level data
tuzvo_stand <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/tuzvo/euf_tuzvo_stand.csv", sep = ",")
tuzvo_stand <- read.csv("~/data/euforia/tuzvo/euf_tuzvo_stand.csv", sep = ",")

tuzvo_stand <- tuzvo_stand |>
  rename(
    plotID = Plot_ID,
    lon = Longitude,
    lat = Latitude,
    plotsize = Plot_size
  ) |>
  select(plotID, lon, lat, plotsize) |>
  distinct()

# Tree-level data
tuzvo_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/tuzvo/euf_tuzvo_tree.csv", sep = ",")
tuzvo_tree <- read.csv("~/data/euforia/tuzvo/euf_tuzvo_tree.csv", sep = ",")

# prepare data
# we use wood density from the BIOMASS pkg to estimate biomass given volume (m3)
WD <- getWoodDensity(
  genus = word(tuzvo_tree$Species, 1),
  species = word(tuzvo_tree$Species, 2)
)

tuzvo_tree <- tuzvo_tree |>
  rename(
    plotID = Plot_ID,
    year = Inventory_year,
    species = Species,
    status = Status,
    dbh = Tree_DBH
  ) |>
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    management_since_census1_yrs = NA,
    country = "Slovakia",
    biomass = volume * WD$meanWD * 10^3
  ) |>
  filter(status == "alive") |> # alive trees
  left_join(tuzvo_stand)

# aggregate data from stand to stand data
df_tuzvo <- from_tree_data(tuzvo_tree) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 2,
    management = 0
  )

# add coords and biomes
df_tuzvo <- biomes_coords_latlon(df_tuzvo)

# add coords and aridity index
df_tuzvo <- ai_coords_latlon(df_tuzvo)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_tuzvo <- lai_coords_latlon(df_tuzvo)

# add coords and N deposition (Lamarque 2011)
df_tuzvo <- ndep_coords_latlon(df_tuzvo)

# add coords and C:N ratio (ISRIC WISE)
df_tuzvo <- cn_coords_latlon(df_tuzvo)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_tuzvo <- phos_coords_latlon(df_tuzvo)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_tuzvo <- orgc_coords_latlon(df_tuzvo)

ggplot() +
  geom_point(data = df_tuzvo, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_tuzvo, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_tuzvo, file = file.path(here::here(), "/data/inputs/df_tuzvo.rds"))

## ul ----
# Thomas Nagel
# Management: 
# Original info: Plots are considered to be old forest with no management.
# Derived info: We assign management_cat = 2 and management_since_census1_yrs = NA

# tree-level data
ul_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/ul/euf_ul_tree.csv", sep = ",")
ul_tree <- read.csv("~/data/euforia/ul/euf_ul_tree.csv", sep = ",")

# metadata
ul_metadata <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/ul/euf_ul_meta.csv", sep = ",")
ul_metadata <- read.csv("~/data/euforia/ul/euf_ul_meta.csv", sep = ",")

ul_metadata <- ul_metadata |>
  select(plotid, plot_size) |>
  distinct() |>
  rename(plotsize = plot_size)

ul_coords <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/ul/euf_ul_coor.csv", sep = ",")
ul_coords <- read.csv("~/data/euforia/ul/euf_ul_coor.csv", sep = ",")

ul_coords <- ul_coords |>
  select(plotid, lat, lon)

# prepare data
ul_tree <- ul_tree |>
  left_join(ul_metadata) |>
  left_join(ul_coords) |>
  rename(plotID = plotid) |>
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    management_since_census1_yrs = NA,
    country = "Slovenia",
    biomass = NA
  ) |>
  filter(status == "alive") # alive trees

# aggregate data from stand to stand data
df_ul <- from_tree_data(ul_tree) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 2,
    management = 0
  )

# add coords and biomes
df_ul <- biomes_coords_latlon(df_ul)

# add coords and aridity index
df_ul <- ai_coords_latlon(df_ul)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_ul <- lai_coords_latlon(df_ul)

# add coords and N deposition (Lamarque 2011)
df_ul <- ndep_coords_latlon(df_ul)

# add coords and C:N ratio (ISRIC WISE)
df_ul <- cn_coords_latlon(df_ul)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_ul <- phos_coords_latlon(df_ul)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_ul <- orgc_coords_latlon(df_ul)

ggplot() +
  geom_point(data = df_ul, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_ul, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_ul, file = file.path(here::here(), "/data/inputs/df_ul.rds"))

## unito ----
# Renzo Motta

# Management: 
# Original info: No information
# Derived info: In Table S1, assign to 50 years. Need to confirm with data providers.

# Tree-level data
unito <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/unito/euf_unito_tree.csv", sep = ",")
unito <- read.csv("~/data/euforia/unito/euf_unito_tree.csv", sep = ",")

# Divide plot into grids of different size given the coordinates
# plotsize = 100 * 100 = 1 ha
# New plots = 20x20 = 0.04
unito <- unito |>
  filter(lat != 0) |>
  unite(plotID, plotID, LPI_ID, sep = "_")

ggplot(unito) +
  geom_point(aes(lon, lat, col = plotID))

unito_1 <- unito |>
  filter(plotID == "9_1")
ggplot(unito_1) +
  geom_point(aes(lon, lat, col = plotID))

unito_1$grid20 <- interaction(
  cut(unito_1$lon,
    breaks = seq(min(unito_1$lon), max(unito_1$lon), by = 20),
    include.lowest = TRUE
  ),
  cut(unito_1$lat,
    breaks = seq(min(unito_1$lat), max(unito_1$lat), by = 20),
    include.lowest = TRUE
  ),
  sep = "X"
)
unito_1 <- unito_1 |>
  rename(plotIDD = plotID) |>
  group_by(grid20) |>
  mutate(plotID = cur_group_id()) |>
  ungroup() |>
  unite(plotID, plotID, plotIDD, sep = "_")
ggplot(unito_1) +
  geom_point(aes(lon, lat, col = plotID))

# rewrite lat and lon to have info at plot level not tree level
unito_1 <- unito_1 |>
  mutate(
    lon = lon_plot,
    lat = lat_plot
  )
ggplot(unito_1) +
  geom_point(aes(lon, lat))

unito_2 <- unito |>
  filter(plotID == "9_2")
ggplot(unito_2) +
  geom_point(aes(lon, lat, col = plotID))

unito_2$grid20 <- interaction(
  cut(unito_2$lon,
    breaks = seq(min(unito_2$lon), max(unito_2$lon), by = 20),
    include.lowest = TRUE
  ),
  cut(unito_2$lat,
    breaks = seq(min(unito_2$lat), max(unito_2$lat), by = 20),
    include.lowest = TRUE
  ),
  sep = "X"
)
unito_2 <- unito_2 |>
  rename(plotIDD = plotID) |>
  group_by(grid20) |>
  mutate(plotID = cur_group_id()) |>
  ungroup() |>
  unite(plotID, plotID, plotIDD, sep = "_")
ggplot(unito_2) +
  geom_point(aes(lon, lat, col = plotID))

# rewrite lat and lon to have info at plot level not tree level
unito_2 <- unito_2 |>
  mutate(
    lon = lon_plot,
    lat = lat_plot
  )
ggplot(unito_2) +
  geom_point(aes(lon, lat))

unito_3 <- unito |>
  filter(plotID == "9_3")
ggplot(unito_3) +
  geom_point(aes(lon, lat, col = plotID))

unito_3$grid20 <- interaction(
  cut(unito_3$lon,
    breaks = seq(min(unito_3$lon), max(unito_3$lon), by = 20),
    include.lowest = TRUE
  ),
  cut(unito_3$lat,
    breaks = seq(min(unito_3$lat), max(unito_3$lat), by = 20),
    include.lowest = TRUE
  ),
  sep = "X"
)
unito_3 <- unito_3 |>
  rename(plotIDD = plotID) |>
  group_by(grid20) |>
  mutate(plotID = cur_group_id()) |>
  ungroup() |>
  unite(plotID, plotID, plotIDD, sep = "_")
ggplot(unito_3) +
  geom_point(aes(lon, lat, col = plotID))

# rewrite lat and lon to have info at plot level not tree level
unito_3 <- unito_3 |>
  mutate(
    lon = lon_plot,
    lat = lat_plot
  )
ggplot(unito_3) +
  geom_point(aes(lon, lat))

unito <- unito_1 |>
  bind_rows(unito_2) |>
  bind_rows(unito_3)

# check # plots
unito |> distinct(plotID)

# prepare data
# we use wood density from the BIOMASS pck to estimate biomass given volume (m3)
WD <- getWoodDensity(
  genus = word(unito$species, 1),
  species = word(unito$species, 2)
)

unito <- unito |>
  rename(dbh = DBH) |>
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    management_since_census1_yrs = NA,
    country = "Italy",
    management = 0,
    biomass = volume * WD$meanWD * 10^3,
    plotsize = 0.04
  ) |>
  filter(status == "A") # alive trees

# aggregate data from stand to stand data
df_unito <- from_tree_data(unito) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 2,
    management = 0
  )

# add coords and biomes
df_unito <- biomes_coords_utm(df_unito)

# add coords and aridity index
df_unito <- ai_coords_latlon(df_unito)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_unito <- lai_coords_latlon(df_unito)

# add coords and N deposition (Lamarque 2011)
df_unito <- ndep_coords_latlon(df_unito)

# add coords and C:N ratio (ISRIC WISE)
df_unito <- cn_coords_latlon(df_unito)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_unito <- phos_coords_latlon(df_unito)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_unito <- orgc_coords_latlon(df_unito)

ggplot() +
  geom_point(data = df_unito, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_unito, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_unito, file = file.path(here::here(), "/data/inputs/df_unito.rds"))

## urk ----
# Data providers: Srdjan Keren and Zbigniew Maciejewski
# Management:
# Original info: year_last_management
# Derived info: In Table S1, assign to NA years. It can be considered a primary/old-growth forests.

# Stand-level data
urk <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/urk/euf_urk_stand.csv", sep = ",")
urk <- read.csv("~/data/euforia/urk/euf_urk_stand.csv", sep = ",")

# aggregate data from stand
df_urk <- from_stand_data(urk) |>
  mutate(
    biomass = NA,
    management = 0,
    management_cat = 1,
    country = "Poland"
  ) |>
  rename(plotsize = area_ha) |>
  group_by(plotID) %>%
  mutate(
    management_since_census1_yrs = min(year, na.rm = TRUE) - year_last_management
  ) |>
  ungroup()

# add coords and biomes
df_urk <- biomes_coords_latlon(df_urk)

# add coords and aridity index
df_urk <- ai_coords_latlon(df_urk)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_urk <- lai_coords_latlon(df_urk)

# add coords and N deposition (Lamarque 2011)
df_urk <- ndep_coords_latlon(df_urk)

# add coords and C:N ratio (ISRIC WISE)
df_urk <- cn_coords_latlon(df_urk)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_urk <- phos_coords_latlon(df_urk)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_urk <- orgc_coords_latlon(df_urk)

ggplot() +
  geom_point(data = df_urk, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_urk, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_urk, file = file.path(here::here(), "/data/inputs/df_urk.rds"))

## wuls ----
# Data providers: Bogdan Brzeziecki

# Management:
# Original info: not management info.
# Derived info: In Table S1, assign management_since_census1_yrs = ca. 20 years. It can be considered a primary/old-growth forests.

# Tree-level data
wuls_tree <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/wuls/euf_wuls_tree.csv", sep = ",")
wuls_census <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/wuls/euf_wuls_census.csv", sep = ",")
wuls_areas <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/wuls/euf_wuls_area.csv", sep = ",")
wuls_coords <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/wuls/euf_wuls_coor.csv", sep = ",")
wuls_species <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/euforia/wuls/euf_wuls_sp.csv", sep = ",")

wuls_tree <- read.csv("~/data/euforia/wuls/euf_wuls_tree.csv", sep = ",")
wuls_census <- read.csv("~/data/euforia/wuls/euf_wuls_census.csv", sep = ",")
wuls_areas <- read.csv("~/data/euforia/wuls/euf_wuls_area.csv", sep = ",")
wuls_coords <- read.csv("~/data/euforia/wuls/euf_wuls_coor.csv", sep = ",")
wuls_species <- read.csv("~/data/euforia/wuls/euf_wuls_sp.csv", sep = ",")

wuls <- wuls_tree |>
  select(-c(H1, H2, H3, H4, H5, H6, H7, H8, V1, V2, V3, V4, V5, V6, V7, V8)) |>
  pivot_longer(
    cols = c(
      dbh_1, dbh_2, dbh_3, dbh_4, dbh_5, dbh_6, dbh_7, dbh_8,
      status_1, status_2, status_3, status_4, status_5, status_6, status_7, status_8
    ),
    names_to = c(".value", "census"),
    names_sep = "_",
    values_drop_na = FALSE
  ) |>
  filter(status != 0 & status != 3) |> # 0=absent or 3=dead
  rename(code = species) |>
  mutate(
    dbh = dbh / 10, # from mm to cm
    census = as.numeric(census),
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    biomass = NA,
    management_since_census1_yrs = NA,
    country = "Poland"
  ) |>
  unite(plotIDD, transect, plotID, sep = "_", remove = FALSE) |>
  left_join(wuls_census) |>
  left_join(wuls_areas) |>
  left_join(wuls_coords) |>
  left_join(wuls_species) |>
  rename(
    plotID = plotIDD,
    plot_ID = plotID,
    plotsize = ha
  )

# aggregate data from stand to stand data
df_wuls <- from_tree_data(wuls) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 2,
    management = 0
  )

# add coords and biomes
df_wuls <- biomes_coords_latlon(df_wuls)

# add coords and aridity index
df_wuls <- ai_coords_latlon(df_wuls)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_wuls <- lai_coords_latlon(df_wuls)

# add coords and N deposition (Lamarque 2011)
df_wuls <- ndep_coords_latlon(df_wuls)

# add coords and C:N ratio (ISRIC WISE)
df_wuls <- cn_coords_latlon(df_wuls)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_wuls <- phos_coords_latlon(df_wuls)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_wuls <- orgc_coords_latlon(df_wuls)

ggplot() +
  geom_point(data = df_wuls, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_wuls, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_wuls, file = file.path(here::here(), "/data/inputs/df_wuls.rds"))

## Luquillo ----
# Contact: Jess Zimmerman
# Management:
# Not info about management
# See https://forestgeo.si.edu/sites/north-america/luquillo

# Tree-level data
luquillo <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/luquillo", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

luquillo <- list.files(path = "~/data/forestgeo/luquillo", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

## Prepare data
luquillo <- luquillo |>
  rename(
    plotID = Quadrat,
    census = Census,
    dbh = DBH,
    species = Latin,
    status = Status
  ) |>
  mutate(
    lon = -65.8160,
    lat = 18.3262,
    dbh = dbh * 0.1,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotsize = 0.04,
    date = as.Date(Date, "%d.%m.%y"),
    year = year(date),
    year = ifelse(census == 1, 1990, year),
    year = ifelse(census == 2, 1994, year),
    year = ifelse(census == 3, 2000, year),
    year = ifelse(census == 4, 2005, year),
    year = ifelse(census == 5, 2011, year),
    year = ifelse(census == 6, 2016, year),
    biomass = NA,
    management_since_census1_yrs = NA,
    country = "Puerto Rico",
    species = ifelse(species == "Ficus spp", "Ficus spp.", species)
  ) |>
  filter(status == "alive") |>
  filter(dbh > 0) |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

# aggregate data from stand to stand data
df_luquillo <- from_tree_data(luquillo) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 2,
    management = 0
  )

# add coords and biomes
df_luquillo <- biomes_coords_latlon(df_luquillo)

# add coords and aridity index
df_luquillo <- ai_coords_latlon(df_luquillo)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_luquillo <- lai_coords_latlon(df_luquillo)

# add coords and N deposition (Lamarque 2011)
df_luquillo <- ndep_coords_latlon(df_luquillo)

# add coords and C:N ratio (ISRIC WISE)
df_luquillo <- cn_coords_latlon(df_luquillo)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_luquillo <- phos_coords_latlon(df_luquillo)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_luquillo <- orgc_coords_latlon(df_luquillo)

ggplot() +
  geom_point(data = df_luquillo, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_luquillo, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_luquillo, file = file.path(here::here(), "/data/inputs/df_luquillo.rds"))

## BCI ----
# contact: Salomon Aguilar, technician in BCI
# finally not added to the global dataset due to difficulty to exclude management (communication with  Salomón Aguilar, Research Technician)
# https://forestgeo.si.edu/sites/neotropics/barro-colorado-island

bci <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/bci", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

bci <- list.files(path = "~/data/forestgeo/bci", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows() 

## Prepare data
bci <- bci |>
  select(-species) |>
  drop_na(DBH) |>
  rename(
    plotID = Quadrat,
    census = Census,
    dbh = DBH,
    species = Latin,
    status = Status
  ) |>
  mutate(
    lon = -79.8461,
    lat = 9.1543,
    dbh = dbh * 0.1,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotsize = 0.04,
    date = as.Date(Date, "%d.%m.%y"),
    year = year(date),
    year = ifelse(census == 1, 1980, year),
    year = ifelse(census == 2, 1985, year),
    year = ifelse(census == 3, 1990, year),
    year = ifelse(census == 4, 1995, year),
    year = ifelse(census == 5, 2000, year),
    year = ifelse(census == 6, 2005, year),
    year = ifelse(census == 7, 2010, year),
    year = ifelse(census == 8, 2015, year),
    biomass = NA,
    management_since_census1_yrs = NA,
    country = "Panama"
  ) |>
  filter(dbh > 0) |>
  filter(status == "alive") |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year) |>
  # Filter those plots for mature stands (see map orange area, adviced by Salomón Aguilar)
  # 500*(10/25) = 320
  filter(PY < 320)

# aggregate data from tree to stand data
df_bci <- from_tree_data(bci) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 3,
    management = 0
  )

# check that lianas are excluded
sp_lianas <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/bci/sp_lianas.csv")
df_bci <- df_bci %>%
  filter(!df_bci$species %in% sp_lianas$species)

# add coords and biomes
df_bci <- biomes_coords_latlon(df_bci)

# add coords and aridity index
df_bci <- ai_coords_latlon(df_bci)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_bci <- lai_coords_latlon(df_bci)

# add coords and N deposition (Lamarque 2011)
df_bci <- ndep_coords_latlon(df_bci)

# add coords and C:N ratio (ISRIC WISE)
df_bci <- cn_coords_latlon(df_bci)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_bci <- phos_coords_latlon(df_bci)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_bci <- orgc_coords_latlon(df_bci)

ggplot() +
  geom_point(data = df_bci, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_bci, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_bci, file = file.path(here::here(), "/data/inputs/df_bci.rds"))

## SCBI ----
# Contact: Kristina J. Anderson-Teixeira, William J. McShea, Norman A. Bourg
# https://forestgeo.si.edu/sites/north-america/smithsonian-conservation-biology-institute
# Management:
# No info in dataset
# In Table S1, management_since_census1_yrs = 53, category = 1

scbi <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/scbi", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  lapply(\(x) mutate(x, across(Tag, as.double))) |>
  lapply(\(x) mutate(x, across(StemTag, as.double))) |>
  bind_rows()

scbi <- list.files(path = "~/data/forestgeo/scbi", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  lapply(\(x) mutate(x, across(Tag, as.double))) |>
  lapply(\(x) mutate(x, across(StemTag, as.double))) |>
  bind_rows() 

## Prepare data
scbi <- scbi |>
  rename(
    plotID = Quadrat,
    census = Census,
    dbh = DBH,
    species = Latin,
    status = Status
  ) |>
  mutate(
    lon = -78.1454,
    lat = 38.8935,
    dbh = dbh * 0.1,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotsize = 0.04,
    date = as.Date(Date, "%d.%m.%y"),
    year = year(date),
    year = ifelse(census == 1, 2008, year),
    year = ifelse(census == 2, 2013, year),
    year = ifelse(census == 3, 2018, year),
    biomass = NA,
    management_since_census1_yrs = 53,
    country = "USA",
    species = ifelse(species == "Carya sp", "Carya spp.", species)
  ) |>
  filter(dbh > 0) |>
  filter(status == "alive") |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

# aggregate data from stand to stand data
df_scbi <- from_tree_data(scbi) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 1,
    management = 0
  )

# add coords and biomes
df_scbi <- biomes_coords_latlon(df_scbi)

# add coords and aridity index
df_scbi <- ai_coords_latlon(df_scbi)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_scbi <- lai_coords_latlon(df_scbi)

# add coords and N deposition (Lamarque 2011)
df_scbi <- ndep_coords_latlon(df_scbi)

# add coords and C:N ratio (ISRIC WISE)
df_scbi <- cn_coords_latlon(df_scbi)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_scbi <- phos_coords_latlon(df_scbi)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_scbi <- orgc_coords_latlon(df_scbi)

ggplot() +
  geom_point(data = df_scbi, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_scbi, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_scbi, file = file.path(here::here(), "/data/inputs/df_scbi.rds"))

## Palanam ----
# Contact: Perry S. Ong

palanam <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/palanam", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  lapply(\(x) mutate(x, across(Tag, as.double))) |>
  lapply(\(x) mutate(x, across(StemTag, as.double))) |>
  bind_rows()

palanam <- list.files(path = "~/data/forestgeo/palanam", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  lapply(\(x) mutate(x, across(Tag, as.double))) |>
  lapply(\(x) mutate(x, across(StemTag, as.double))) |>
  bind_rows()  

## Prepare data
palanam <- palanam |>
  rename(
    plotID = Quadrat,
    census = Census,
    dbh = DBH,
    species = Latin,
    status = Status
  ) |>
  mutate(
    lon = 122.3880,
    lat = 17.0402,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotsize = 0.04,
    date = as.Date(Date, "%d.%m.%y"),
    year = year(date),
    year = ifelse(census == 1, 1994, year),
    biomass = NA,
    management_since_census1_yrs = NA,
    country = "Philippines"
  ) |>
  filter(dbh > 0) |>
  filter(status == "alive") |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

# aggregate data from stand to stand data
df_palanam <- from_tree_data(palanam) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  filter(logQMD > 1)  |>
  mutate(
    management_cat = 2,
    management = 0
  )

# add coords and biomes
df_palanam <- biomes_coords_latlon(df_palanam)

# add coords and aridity index
df_palanam <- ai_coords_latlon(df_palanam)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_palanam <- lai_coords_latlon(df_palanam)

# add coords and N deposition (Lamarque 2011)
df_palanam <- ndep_coords_latlon(df_palanam)

# add coords and C:N ratio (ISRIC WISE)
df_palanam <- cn_coords_latlon(df_palanam)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_palanam <- phos_coords_latlon(df_palanam)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_palanam <- orgc_coords_latlon(df_palanam)

ggplot() +
  geom_point(data = df_palanam, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_palanam, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_palanam, file = file.path(here::here(), "/data/inputs/df_palanam.rds"))

## SERC ----
# Contact: Sean McMahon
# https://forestgeo.si.edu/sites/north-america/smithsonian-environmental-research-center

serc <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/serc", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()
unique(serc$Census)
length(unique(serc$Quadrat))

serc <- list.files(path = "~/data/forestgeo/serc", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

# Divide plot into grids of different size given the coordinates
# 20x20m = 0.04 ha
serc$grid20 <- interaction(
  cut(serc$PX,
    breaks = seq(0, 400, by = 20),
    include.lowest = TRUE
  ),
  cut(serc$PY,
    breaks = seq(0, 400, by = 20),
    include.lowest = TRUE
  ),
  sep = "X"
)
serc <- serc |>
  group_by(grid20) |>
  mutate(plotID = cur_group_id()) |>
  ungroup()
length(unique(serc$plotID))
ggplot(serc) +
  geom_point(aes(PX, PY, col = plotID))

## Prepare data
serc <- serc |>
  rename(
    census = Census,
    dbh = DBH,
    species = Latin,
    status = Status
  ) |>
  mutate(
    lon = -76.5594,
    lat = 38.8891,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotsize = 0.04,
    date = as.Date(Date, "%d.%m.%y"),
    year = year(date),
    year = ifelse(census == 1, 2010, year),
    year = ifelse(census == 3, 2019, year),
    biomass = NA,
    years_since_management = NA,
    management = 0,
    country = "USA"
  ) |>
  filter(dbh > 0) |>
  filter(status == "alive") |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)
unique(serc$year)
unique(serc$census)

# aggregate data from stand to stand data
data_serc <- from_tree_data(serc) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  filter(logQMD > 1)

# add coords and biomes
data_serc <- biomes_coords_latlon(data_serc)

# add coords and aridity index
data_serc <- ai_coords_latlon(data_serc)

# add coords and LAI Modis or NDVI (0.5 degrees)
data_serc <- lai_coords_latlon(data_serc)

# add coords and N deposition (Lamarque 2011)
data_serc <- ndep_coords_latlon(data_serc)

# add coords and C:N ratio (ISRIC WISE)
data_serc <- cn_coords_latlon(data_serc)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
data_serc <- phos_coords_latlon(data_serc)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
data_serc <- orgc_coords_latlon(data_serc)

ggplot() +
  geom_point(data = data_serc, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = data_serc, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(data_serc, file = file.path(here::here(), "/data/inputs/data_serc.rds"))

## Wytham ----
# Data from Yadvinder Malhi

# tree-level data
wytham <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/wytham/wytham.csv")
wytham_sp <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/wytham/wytham_sp.csv")

# prepare data
wytham <- wytham |>
  pivot_longer(
    cols = c(dbh_2008, dbh_2010, dbh_2016, dbh_2021, codes_2010, codes_2016, codes_2021),
    names_to = c(".value", "year"),
    names_sep = "_",
    values_drop_na = FALSE
  ) |>
  left_join(wytham_sp) |>
  rename(status = codes) |>
  mutate(
    lon = -1.3379,
    lat = 51.7743,
    dbh = dbh * 0.1,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotsize = 0.04,
    year = as.integer(year),
    biomass = NA,
    years_since_management = NA,
    management = 0,
    country = "UK"
  ) |>
  filter(dbh > 0) |>
  filter(status != "D") |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

# aggregate data from stand to stand data
data_wytham <- from_tree_data(wytham) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup()

# add coords and biomes
data_wytham <- biomes_coords_latlon(data_wytham)

# add coords and aridity index
data_wytham <- ai_coords_latlon(data_wytham)

# add coords and LAI Modis or NDVI (0.5 degrees)
data_wytham <- lai_coords_latlon(data_wytham)

# add coords and N deposition (Lamarque 2011)
data_wytham <- ndep_coords_latlon(data_wytham)

# add coords and C:N ratio (ISRIC WISE)
data_wytham <- cn_coords_latlon(data_wytham)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
data_wytham <- phos_coords_latlon(data_wytham)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
data_wytham <- orgc_coords_latlon(data_wytham)

ggplot() +
  geom_point(data = data_wytham, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = data_wytham, size = 0.5, alpha = 0.5)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = ai), data = data_wytham, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(data_wytham, file = file.path(here::here(), "/data/inputs/data_wytham.rds"))

## Pasho ----
# Contact: Musalmah Nasardin

# tree-level data
pasoh_sp <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh_sp.rds")
pasoh_sp <- pasoh_sp |>
  separate(Latin, c("first", "second"), " ") |>
  mutate(second = tolower(second)) |>
  unite(Latin, first, second, sep = " ")
pasoh1 <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh1.rds")
pasoh2 <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh2.rds")
pasoh3 <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh3.rds")
pasoh4 <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh4.rds")
pasoh5 <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh5.rds")
pasoh6 <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh6.rds")
pasoh7 <- readRDS("/home/laura/data/forestgeo/pasoh/pasoh7.rds")

pasoh <- pasoh1 |>
  bind_rows(pasoh2) |>
  bind_rows(pasoh3) |>
  bind_rows(pasoh4) |>
  bind_rows(pasoh5) |>
  bind_rows(pasoh6) |>
  bind_rows(pasoh7) |>
  left_join(pasoh_sp)

pasoh <- pasoh |>
  rename(
    plotID = quadrat,
    census = CensusID,
    species = Latin,
    PX = gx,
    PY = gy
  ) |>
  mutate(
    lon = 102.30816760,
    lat = 2.97956222,
    dbh = dbh * 0.1, # in cm
    ba_tree = pi * (dbh * 0.01 / 2)^2, # in m2
    plotsize = 0.04,
    date = as.Date(ExactDate, "%Y-%m-%d"),
    year = year(date),
    biomass = agb * 10^3,
    # year = ifelse(census == 1, 1990, year),
    years_since_management = NA,
    management = 0,
    country = "Malaysia"
  ) |>
  filter(DFstatus == "alive") |>
  filter(dbh > 0) |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

# aggregate data from stand to stand data
data_pasoh <- from_tree_data(pasoh) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(
    census = match(year, unique(year)),
    biomass = ifelse(biomass == 0, NA, biomass)
  ) |>
  ungroup()

# add coords and biomes
data_pasoh <- biomes_coords_latlon(data_pasoh)

# add coords and aridity index
data_pasoh <- ai_coords_latlon(data_pasoh)

# add coords and LAI Modis or NDVI (0.5 degrees)
data_pasoh <- lai_coords_latlon(data_pasoh)

# add coords and N deposition (Lamarque 2011)
data_pasoh <- ndep_coords_latlon(data_pasoh)

# add coords and C:N ratio (ISRIC WISE)
data_pasoh <- cn_coords_latlon(data_pasoh)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
data_pasoh <- phos_coords_latlon(data_pasoh)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
data_pasoh <- orgc_coords_latlon(data_pasoh)

# Remove the low density plots
data_pasoh <- data_pasoh |>
  filter(density > 2000)

# Check number of plots before filtering
length(unique(data_pasoh$plotID))

# Check min dbh
min(data_pasoh$QMD)

ggplot() +
  geom_point(data = data_pasoh, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = data_pasoh, size = 0.5, alpha = 0.5)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = ai), data = data_pasoh, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(data_pasoh, file = file.path(here::here(), "/data/inputs/data_pasoh.rds"))

## Mudumalai ----
# Contact: Prof. Sukumar

# Stand-level data
mudumalai <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestgeo/mudumalai/mudumalai_stand.csv", sep = ",")
mudumalai <- mudumalai |>
  mutate(density = nindiv / plotsize) |>
  mutate(ba = ba / plotsize)

# aggregate data from stand
data_mudumalai <- from_stand_data(mudumalai) |>
  # create census
  group_by(plotID) |>
  # mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management = 0,
    country = "India",
    years_since_management = NA,
    biomass = NA
  )

# add coords and biomes
data_mudumalai <- biomes_coords_latlon(data_mudumalai)
data_mudumalai |>
  distinct(biomeID)

# Clasify biome manually as Tropical & Subtropical Dry Broadleaf Forests
data_mudumalai <- data_mudumalai |>
  mutate(
    biomeID = 2,
    biome = "Tropical & Subtropical Dry Broadleaf Forests"
  )

# add coords and aridity index
data_mudumalai <- ai_coords_latlon(data_mudumalai)

# add coords and LAI Modis or NDVI (0.5 degrees)
data_mudumalai <- lai_coords_latlon(data_mudumalai)

# add coords and N deposition (Lamarque 2011)
data_mudumalai <- ndep_coords_latlon(data_mudumalai)

# add coords and C:N ratio (ISRIC WISE)
data_mudumalai <- cn_coords_latlon(data_mudumalai)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
data_mudumalai <- phos_coords_latlon(data_mudumalai)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
data_mudumalai <- orgc_coords_latlon(data_mudumalai)

# Check number of plots before filtering
length(unique(data_mudumalai$plotID))

# Check min dbh
min(data_mudumalai$QMD)

ggplot() +
  geom_point(data = data_mudumalai, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = data_mudumalai, size = 0.5, alpha = 0.5)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = ai), data = data_mudumalai, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(data_mudumalai, file = file.path(here::here(), "/data/inputs/data_mudumalai.rds"))

# Forests plots MZ ----
# Data provided by Yadvinder Mahli and Huanyuan Zhang

# 0) Metadata
meta_fp <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestplots/mahli_zhang/fp_mz_metadata.csv")
meta_fp <- meta_fp |>
  select(Plot_code, Country, Longitude, Latitude, Elevation..m., Plot.Size..ha.) |>
  rename(
    plotID = Plot_code,
    country = Country,
    lon = Longitude,
    lat = Latitude,
    elevation = Elevation..m.,
    plotsize = Plot.Size..ha.
  )

# 1) A selection of countries and plots together
multiplots <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestplots/mahli_zhang/fp_mz_multiplots.csv")

# select and rename variables
multiplots <- multiplots |>
  select(Country, Plot.Code, TreeID, Species, Census.Date, D0, D1, D2, D3, D4) |>
  rename(
    plotID = Plot.Code,
    treeID = TreeID,
    year = Census.Date,
    country = Country,
    species = Species
  ) |>
  mutate(year = as.integer(year)) |>
  mutate(
    dbh = rowMeans(across(D0:D4), na.rm = TRUE),
    dbh = dbh / 10
  ) |>
  filter(dbh != 0) |>
  select(-c(D0, D1, D2, D3, D4)) |>
  left_join(meta_fp |> select(-country), by = "plotID") |>
  as_tibble()
unique(multiplots$country)

# 2) Gabon
gabon <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestplots/mahli_zhang/gabon", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

# select and rename variables
gabon <- gabon |>
  select(country, plot_code, treeID, species, year, dbh) |>
  rename(plotID = plot_code) |>
  mutate(year = as.integer(year)) |>
  mutate(dbh = dbh / 10) |>
  filter(dbh != 0) |>
  left_join(meta_fp |> select(-country), by = "plotID")

# 3) Ghana
ghana <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestplots/mahli_zhang/ghana", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

# select and rename variables
ghana <- ghana |>
  select(country, plotID, treeID, species, year, dbh) |>
  mutate(year = as.integer(year)) |>
  mutate(dbh = dbh / 10) |>
  filter(dbh != 0) |>
  left_join(meta_fp |> select(-country), by = "plotID")

# 4) Malaysia
malaysia <- list.files(path = "/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/forestplots/mahli_zhang/malaysia", full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

# select and rename variables
malaysia <- malaysia |>
  select(country, plotID, treeID, species, year, dbh) |>
  mutate(year = as.integer(year)) |>
  mutate(dbh = dbh / 10) |>
  filter(dbh != 0) |>
  left_join(meta_fp |> select(-country), by = "plotID")

# Join all plots
df_forestplots <- multiplots |>
  bind_rows(gabon) |>
  bind_rows(ghana) |>
  bind_rows(malaysia) |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    biomass = NA,
    management = 0,
    years_since_management = NA
  )

# aggregate data from stand to stand data
data_forestplots <- from_tree_data(df_forestplots) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup()

# add coords and biomes
data_forestplots <- biomes_coords_latlon(data_forestplots)

# add coords and aridity index
data_forestplots <- ai_coords_latlon(data_forestplots)

# add coords and LAI Modis or NDVI (0.5 degrees)
data_forestplots <- lai_coords_latlon(data_forestplots)

# add coords and N deposition (Lamarque 2011)
data_forestplots <- ndep_coords_latlon(data_forestplots)

# add coords and C:N ratio (ISRIC WISE)
data_forestplots <- cn_coords_latlon(data_forestplots)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
data_forestplots <- phos_coords_latlon(data_forestplots)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
data_forestplots <- orgc_coords_latlon(data_forestplots)

# check outliers
data_forestplots <- data_forestplots |>
  filter(logQMD < 4.5)

ggplot() +
  geom_point(data = data_forestplots, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = data_forestplots, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(data_forestplots, file = file.path(here::here(), "/data/inputs/data_mz_forestplots.rds"))

# Australian plots ####
# Data providers: David Forrester and co.

# Stand-level data
aus_plots <- readRDS("/home/laura/data/fp_aus/fp_aus_stand.RDS")

# prepare data
aus_plots <- aus_plots |>
  rename(
    lon = longitude,
    lat = latitude,
    density = TreesPerHectareAHC1_2,
    dbh = DBHqAHC1_2_cm,
    ba = BasalAreaAHC1_2_m2perha,
    plotsize = plot_size,
    species = dominant_species
  ) |>
  mutate(
    biomass = AbovegroundAHC1_2_Mgperha * 10^3,
    country = "Australia",
    years_since_management = year(Sys.Date()) - year_last_management,
    year = as.integer(year)
  )

# aggregate data from stand to stand data
data_aus <- from_stand_data(aus_plots) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(management = 0)

# add coords and biomes
data_aus <- biomes_coords_latlon(data_aus)

# add coords and aridity index
data_aus <- ai_coords_latlon(data_aus)

# add coords and LAI Modis or NDVI (0.5 degrees)
data_aus <- lai_coords_latlon(data_aus)

# add coords and N deposition (Lamarque 2011)
data_aus <- ndep_coords_latlon(data_aus)

# add coords and C:N ratio (ISRIC WISE)
data_aus <- cn_coords_latlon(data_aus)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
data_aus <- phos_coords_latlon(data_aus)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
data_aus <- orgc_coords_latlon(data_aus)

ggplot() +
  geom_point(data = data_aus, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = data_aus, size = 0.5, alpha = 0.5) +
  theme(legend.position = "bottom")

# save stand-level data
saveRDS(data_aus, file = file.path(here::here(), "/data/inputs/data_fp_aus.rds"))

# RAINFOR Amazon Forest Inventory Network ----
# Data from Esquivel-Muelbert, et al. 2025 data published
# Coordinates of the sites from Bennett et al. 2023

# Stand-level data
rainfor_plots <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/rainfor/esquivel_etal/full_final_dataset.csv", sep = ",")
rainfor_bennet <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/rainfor/bennett_etal/table_s3.csv", sep = ",")
rainfor_brienen <- read.csv("/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw/rainfor/brienen_etal/table_s1.csv", sep = ",")

# select variables from coordinates data
rainfor_bennet <- rainfor_bennet |>
  select(plot_code, country, lon, lat)
rainfor_brienen <- rainfor_brienen |>
  select(plot_code, country, lon, lat)

rainfor_coords <- rainfor_bennet |>
  bind_rows(rainfor_brienen) |>
  arrange(plot_code) |>
  distinct(plot_code, .keep_all = TRUE)

# rename variables
df_rainfor <- rainfor_plots |>
  rename(
    plot_code = Plot.Code,
    ba = total_BA,
    mean_ba = mean_BA,
    plotsize = PlotArea,
    density = number_stems,
    plot_identi = PlotID
  ) |>
  mutate(
    year = round(census_date),
    years_since_management = NA,
    species = NA,
    biomass = NA,
    dbh = 100 * sqrt((4 * mean_ba) / pi)
  ) |> # DBH in centimetres
  left_join(rainfor_coords, by = "plot_code") |>
  rename(plotID = plot_code)

# check how many plots are missing coordinates
df_rainfor %>%
  filter(is.na(lat) | is.na(lon)) %>%
  distinct(plotID)

# keep only the alphabetic part of the plotID to join the rest of coordinates
# since not all plot numbers are included in the coordinates info, we exclude the numeric part of the plotID
rainfor_prefix <- rainfor_coords |>
  mutate(plot_prefix = str_extract(plot_code, "^[^-]+")) |>
  distinct(plot_prefix, .keep_all = TRUE)

df_rainfor <- df_rainfor |>
  mutate(plot_prefix = str_extract(plotID, "^[^-]+")) |>
  left_join(rainfor_prefix |> select(plot_prefix, lon, lat), by = "plot_prefix") |>
  mutate(
    lon = if_else(is.na(lon.x), lon.y, lon.x),
    lat = if_else(is.na(lat.x), lat.y, lat.x)
  ) |>
  select(-lon.x, -lon.y, -lat.x, -lat.y)

# check again how many plots are missing coordinates
df_rainfor %>%
  filter(is.na(lat) | is.na(lon)) %>%
  distinct(plotID)

# aggregate data from stand
df_rainfor <- from_stand_data(df_rainfor) |>
  mutate(management = 0)

# add coords and biomes
df_rainfor <- biomes_coords_latlon(df_rainfor)

# add coords and aridity index
df_rainfor <- ai_coords_latlon(df_rainfor)
# df_rainfor <- ai_coords_latlon_opt(df_rainfor)

# add coords and LAI Modis or NDVI (0.5 degrees)
df_rainfor <- lai_coords_latlon(df_rainfor)

# add coords and N deposition (Lamarque 2011)
df_rainfor <- ndep_coords_latlon(df_rainfor)

# add coords and C:N ratio (ISRIC WISE)
df_rainfor <- cn_coords_latlon(df_rainfor)

# add coords and Phosphorus P - Bray (PBR) or Olsen (POL)
df_rainfor <- phos_coords_latlon(df_rainfor)

# add coords and ORGC - Organic carbon content (g kg-1) (ISRIC WISE)
df_rainfor <- orgc_coords_latlon(df_rainfor)

ggplot() +
  geom_point(data = df_rainfor, aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rbeni::plot_map_simpl() +
  geom_point(aes(lon, lat, color = biome), data = df_rainfor, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_rainfor, file = file.path(here::here(), "/data/inputs/data_df_rainfor.rds"))
