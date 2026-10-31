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

# Execution location ----
# Set to "local" on your computer or "server" on the remote server.
execution_environment <- "local"

data_archive_dir <- switch(
  match.arg(execution_environment, c("local", "server")),
  local = path.expand("~/data/archive_restricted/GFDYglobe_marques_2025"),
  server = "/data/archive_restricted/GFDYglobe_marques_2025"
)
raw_data_dir <- file.path(data_archive_dir, "GFDYglobe_data_raw")

# load packages ----
# library(renv)
library(readr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(lubridate)
library(rFIA)
library(sf)
library(rgeco)
library(here)

# load functions ----
# Aggregation helpers harmonise dataset identifiers for all downstream labels.
# Retain missing tree DBH values through the eligibility filters: if individual
# basal area cannot be recovered, the stand basal-area sum must remain NA.
source(here::here("R/functions.R"))

# NFI Spain ----
# Data providers: Paloma Ruiz-Benito and Veronica Cruz-Alonso

# Management:
# Original info: Only plots were management was not registered during the inventory were provided.
# Management category assigned to 3 for all plots.

# stand-level for species
nfi_spain <- read.csv(file.path(raw_data_dir, "nfi_spa/nfi_spa.csv"))

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
  mutate(plotsize = 0.20) # beni: all have the same value, in ha?

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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_nfi_spain <- coords_to_lonlat(df_nfi_spain, crs = 25830) |>
  standardise_stand_data() |>
  as_tibble()

ggplot() +
  geom_point(
    data = df_nfi_spain,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.5,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_nfi_spain)

# Save stand-level data
saveRDS(
  df_nfi_spain,
  file = here::here("data/inputs/df_nfi_spa.rds")
)

# NFI Sweeden ----
# Data providers: Julian Tijerin-Triviño

# Management:
# Original info: Plots include the management variable (0-unmanaged, 1-managed).
# Retain whole plots only when management is recorded as 0 at every census.
# Assign management_cat = 1 (recorded history). Use stand_age at the first
# available census as a proxy for years since management before that census,
# keeping this baseline interval constant across the plot's later censuses.

# Tree-level for species
species_code <- read.csv(file.path(raw_data_dir, "nfi_swe/swe_sp_code.csv"), sep = ",")

species_code <- species_code |>
  filter(country == "ES") |>
  rename(species = acceptedname) |>
  select(code, species)

nfi_sweeden_tree <- read.csv(file.path(raw_data_dir, "nfi_swe/nfi_swe_tree.csv"), sep = ",")

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
  left_join(
    species_code,
    by = join_by(code)
    )

# Stand-level
nfi_sweeden <- read.csv(file.path(raw_data_dir, "nfi_swe/nfi_swe_stand.csv"), sep = ",")

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

  # Select whole plots with no management at any available census. Missing
  # management indicators cannot establish that the plot remained unmanaged.
  arrange(plotID, year) |>
  group_by(plotID) |>
  filter(all(!is.na(management) & management == 0)) |>
  mutate(
    census = match(year, unique(year)),
    management_since_census1_yrs = first(stand_age)
  ) |>
  ungroup() |>
  relocate(census, .after = year) |>
  mutate(
    plotsize = NA,
    biomass = NA,
    management_cat = 1,
    country = "Sweeden"
  )

# aggregate data from tree to stand data
df_nfi_sweeden <- from_stand_tree_data(nfi_sweeden_tree, nfi_sweeden)

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_nfi_sweeden <- standardise_stand_data(df_nfi_sweeden)

ggplot() +
  geom_line(
    data = df_nfi_sweeden,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.1,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_nfi_sweeden)

# Save stand-level data
saveRDS(
  df_nfi_sweeden,
  file = here::here("data/inputs/df_nfi_swe.rds")
)

# FIA US ----
# Data providers: rFIA R package (doi:10.32614/CRAN.package.rFIA)

# Cache unfiltered stand totals, so changing management criteria can restore plots.
filn <- file.path(raw_data_dir, "fia_us/obs/data_stand_us_v2026.rds")

if (!file.exists(filn)) {
  # Download FIA data ---
  # for all states
  # the dataset needed: COND, PLOT, TREE
  states <- read.csv(file.path(raw_data_dir, "fia_us/obs/states.csv"))
  st <- states$State.abbreviation

  # Data unavailable for:  DC, MH
  st <- st[-which(st %in% c("DC", "MH"))]

  # Don't re-download. They are available locally.
  # for (i in st) {
  #   getFIA(states = i, dir = file.path(raw_data_dir, "fia_us/obs"), tables = "COND", load = FALSE)
  #   getFIA(states = i, dir = file.path(raw_data_dir, "fia_us/obs"), tables = "PLOT", load = FALSE)
  #   options(timeout = 3600)
  #   getFIA(states = i, dir = file.path(raw_data_dir, "fia_us/obs"), tables = "TREE", load = FALSE)
  # }

  # Read data ---
  # UNITCD Survey unit code
  # STATECD State code
  # COUNTYCD County code
  # PLOT Plot number

  ## TREE table ---
  # Aggregate from tree to forest plot level
  # Note: R session crushes due to the big size of the files. So, we read the files and save only the summaries at plot level.

  # Sampling area is already encoded in each tree's TPA_UNADJ. FIA uses nested
  # microplots/subplots/macroplots and legacy designs, so there is no single
  # sampled area applicable to all trees. Do not divide these totals by plotsize.
  # FIA database guide, section 1.2.5 (tree expansion factors):
  # https://research.fs.usda.gov/sites/default/files/2025-08/wo-v9-4_Aug2025_UG_FIADB_database_description_NFI.pdf
  # DRYBIO_AG - aboveground biomass (DRYBIO_AG) contained in each tree (in pounds, libs)
  # TPA_UNADJ - trees per acre each tree represents
  # DIA - DBH (inches)
  # 1 lb = 0.453 kg
  # 1 acre = 0.405 ha
  # 1 inch = 2.54 cm
  # 1 sq inch = 0.00064516 sq meter
  # abg_biomass: from pounds per acre to kg per ha = *0.453/0.405
  # density: from indiv per acre to indiv per ha = */0.405
  # dbh: from inches to cm = *2.54 (no area conversion)
  # BA: from sq. inches per acre to sq. m per ha = *0.00064516/0.405
  species <- read.csv(file.path(raw_data_dir, "fia_us/obs/species_code.csv"))
  states <- read.csv(file.path(raw_data_dir, "fia_us/obs/states.csv"))
  st <- states$State.abbreviation

  # Data unavailable for:  DC, MH
  st <- st[-which(st %in% c("DC", "MH"))]

  data_stand <- data.frame()
  for (i in st) {
    currentDF <- read.csv(file.path(raw_data_dir, "fia_us/obs", paste0(i, "_TREE.csv")))
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
      summarize(ba = sum((pi * DIA * DIA / 4) * TPA_UNADJ * 0.00064516 / 0.405, na.rm = FALSE)) |>
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
        density = sum(TPA_UNADJ / 0.405, na.rm = FALSE),
        dbh = mean(DIA * 2.54, na.rm = TRUE),
        ba = sum((pi * DIA * DIA / 4) * TPA_UNADJ * 0.00064516 / 0.405, na.rm = FALSE)
      ) |>
      ungroup() |>
      mutate(
        QMD = ba_to_qmd(ba, density),
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

  saveRDS(data_stand, file = filn)
} else {
  # Rebuild this cache from TREE tables when changing aggregation/missing-BA rules.
  data_stand <- readRDS(filn)
}

# Refresh derived diameters in older caches, retaining every unfiltered plot.
data_stand <- data_stand |>
  select(-any_of(c("logQMD", "logDensity"))) |>
  mutate(QMD = ba_to_qmd(ba, density))

## PLOT table ---
# meta info for forest plots
data_plot <- list.files(path = file.path(raw_data_dir, "fia_us/obs"), pattern = "*_PLOT.csv", full.names = TRUE) |>
  purrr::map(read.csv) |>
  lapply(\(x) mutate(x, across(ECO_UNIT_PNW, as.character))) |>
  bind_rows() |>

  # Make a unique ID for each plot, irrespective of time
  mutate(plotID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = "_"))

data_plot_sel <- data_plot |>
  select(plotID, LAT, LON, ELEV) |>
  distinct(plotID, .keep_all = T)

# Join tables
df_fia_us <- data_stand |>
  left_join(data_plot_sel) |>

  # rename variable
  rename(
    year = INVYR,
    lat = LAT,
    lon = LON,
    biomass = abg_biomass_kg_ha
  ) |>
  
  # Management classification is assigned from COND history below.
  mutate(
    # PLOT has design codes and COND has area proportions, not a single area
    # in hectares for all sampled diameter classes. Keep this metadata missing;
    # density and ba already use the supplied tree-specific expansion factors.
    plotsize = NA_real_,
    country = "USA"
  ) |>

  # Remove entries with density = 0
  filter(is.na(density) | density > 0) |>

  # Remove entries with QMD = 0
  filter(is.na(QMD) | QMD > 0)

# Retain forest reserves with no evidence of management between censuses.
# Category 2: old growth at the first census; category 1: dated prior treatment
# or a valid first-census stand age (used as a proxy); category 3: unknown history.
# Original COND records are required even when using cached stand totals.
source(here::here("R/prepare_fia_management.R"))
fia_condition_history <- read_fia_condition_history(file.path(raw_data_dir, "fia_us/obs"))
df_fia_us <- prepare_fia_management(df_fia_us, fia_condition_history)

message(
  "FIA management selection: ",
  n_distinct(df_fia_us$plotID),
  " plots, ",
  nrow(df_fia_us),
  " plot-census records retained."
)

# check final data availability
df_fia_us |>
  filter(management == 0) |>
  group_by(management_cat) |>
  summarise(
    n_plots = n_distinct(plotID),
    n_plot_censuses = n(),
    .groups = "drop"
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_fia_us <- standardise_stand_data(df_fia_us)

ggplot() +
  geom_line(
    data = df_fia_us |> 
      filter(management_cat == 2 | (management_cat == 1 & management_since_census1_yrs > 30)),
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.1,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_fia_us)

# Save stand-level data
saveRDS(df_fia_us, file = here::here("data/inputs/df_fia_us_v2026.rds"))

# NFI Switzerland ----
# Data providers: Brigitte Rohner
# Management:
# Original info: Each plots provides the variable LETZTENU, which provides information about
# the number of years since last management before each census (LFI1 to LFI5).
# LETZTENU: Number of years since the last silvicultural intervention according to the oral survey at the Forest Service (derivation LETZTENU in the Swiss NFI database).
# If it is known that no intervention has been carried out on this sample area for a very long time or
# that it is almost certain that no intervention has ever been carried out, LETZTENU=999.
# MANAGEMENT_INDICATOR describes treatment since the preceding inventory; 1 = no treatment.
# Check every census after the first available census, excluding the baseline
# indicator because it concerns history before the observation period (LFI1 uses -1).
# Calculated variables:
# - management = 0 only if all follow-up indicators are known and equal to 1;
#   1 if any follow-up reports treatment; NA if there is no follow-up or status is unknown.
# - management_since_census1_yrs = LETZTENU at the earliest dated census, held
#   constant across censuses. Missing baseline values remain missing; 999 remains
#   a special marker, not a measured interval. Later reports do not fill the baseline.
# - management_cat = 1 for all plots. (Recorded history)
# LETZTENU is relative to the field-measurement year (LFI4 results report, p. 154):
# https://www.lfi.ch/media/documents/publikationen/ergebnisberichte/LFI4_Ergebnisbericht.pdf

# Plot characteristics
lfi_plot_constant <- read.csv(file.path(raw_data_dir, "nfi_swi/lfi_plot_constant.csv"), sep = ";")

# Species names
lfi_species_names <- read.csv(file.path(raw_data_dir, "nfi_swi/lfi_species_names.csv"), sep = ",")

## calculate dom_species from tree-level data
lfi_tree_census <- read.csv(file.path(raw_data_dir, "nfi_swi/lfi_tree_census.csv"), sep = ",")

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
  summarise(ba = sum(ba_tree, na.rm = FALSE)) |>
  slice_max(ba, with_ties = FALSE) |>
  ungroup() |>
  select(PLOTID, year, species)

# Stand-level data
lfi_plot_census <- read.csv(file.path(raw_data_dir, "nfi_swi/lfi_plot_census.csv"), sep = ",")

# prepare data
nfi_switzerland <- lfi_plot_census |>
  left_join(lfi_plot_constant[, c(1:3)]) |>

  # create Year variable from CENSUS_DATE
  mutate(
    census_date = as.Date(CENSUS_DATE, format = "%d.%m.%Y"),
    year = year(census_date)
  ) |>
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
  arrange(plotID, census_date) |>
  group_by(plotID) |>
  mutate(
    management = {
      followup <- MANAGEMENT_INDICATOR[row_number() > 1L]
      if (any(is.na(census_date)) || length(followup) == 0L) NA_integer_
      else if (any(followup > 1L, na.rm = TRUE)) 1L
      else if (all(!is.na(followup) & followup == 1L)) 0L
      else NA_integer_
    },
    management_since_census1_yrs =
      if (any(is.na(census_date))) NA_real_ else as.numeric(first(LETZTENU))
  ) |>
  ungroup() |>
  select(-census_date) |>
  filter(management == 0) |> 
  mutate(
    management_cat = 1,
    plotsize = PLOT_AREA_LARGE * 10^-4,
    country = "Switzerland"
  )

# aggregate data from stand to stand data
df_nfi_switzerland <- from_stand_data(nfi_switzerland)

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_nfi_switzerland <- standardise_stand_data(df_nfi_switzerland)

df_nfi_switzerland |>
  filter(management == 0 & management_since_census1_yrs == 999) |>
  ggplot() +
  geom_hex(
    aes(x = log(QMD), y = log(density)),
    bins = 30
  )

plot_point_map(df_nfi_switzerland)

# Save stand-level data
saveRDS(
  df_nfi_switzerland,
  file = here::here("data/inputs/df_nfi_swi_v2026.rds")
)

# NFI Norway ----
# Data providers: Oliver Moen Snoksrud and Johannes Breidenbach
# Management:
# Original info:
# - management_since_census1_yrs = 64,
# - management_cat = 1 for all plots. (Recorded history)

# Stand-level data
nfi_norway <- read.csv(file.path(raw_data_dir, "nfi_nor/nfi_nor.csv"), sep = ",")

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
    management_cat = 3,
    management_since_census1_yrs = NA, # Beni: revised this. Was previously set to 64 for all plots, which seems implausible.
    country = "Norway"
  )

# Preserve the provider's QMD, converted from mm to cm above.
df_nfi_norway <- from_stand_data(nfi_norway, qmd_col = "qmd") # beni: retaining dataset's own QMD values, not overwriting from own calculation

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_nfi_norway <- standardise_stand_data(df_nfi_norway)

ggplot() +
  geom_line(
    data = df_nfi_norway,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.2,
    # size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_nfi_norway)

# Save stand-level data
saveRDS(
  df_nfi_norway,
  file = here::here("data/inputs/df_nfi_nor_v2026.rds")
)

# EFM Switzerland ----
# David Forrester and Jonas Glatthorn
# Management:
# Original: The dataset efm_management includes the year_last_intervention variable
# Calculate: management_since_census1_yrs from year_last_intervention and first_measurement

# Plot characteristics
efm_metadata <- read.csv(file.path(raw_data_dir, "efm/raw/VFL_LISTmanual.csv"))
efm_metadata <- efm_metadata |>
  mutate(FNUM = as.character(FNUM)) |>
  select(-BA)

# Plot area
efm_area <- read.csv(file.path(raw_data_dir, "efm/raw/EFM_plot_area.csv"))
# efm_area <- read.csv(file.path(raw_data_dir, "efm/raw/PlotAreas.csv"))
efm_area <- efm_area |>
  mutate(FNUM = as.character(FNUM))

efm_locations <- read.csv(file.path(raw_data_dir, "efm/raw/efm_plot_locations.csv"))
efm_locations <- efm_locations |>
  mutate(FNUM = as.character(FNUM))

# Last management intervention
efm_management <- read.csv(file.path(raw_data_dir, "efm/raw/EFM_last_intervention.csv"))
efm_management <- efm_management |>
  mutate(FNUM = as.character(FNUM),
         management_since_census1_yrs = first_measurement - year_last_intervention) |>
  select(FNUM, management_since_census1_yrs)

# Stand level data (per plot, year and species)
efm_stand <- readRDS(file.path(raw_data_dir, "efm/raw/efm_stand_v2026.rds"))
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
    qmd = ba_to_qmd(ba, density),
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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_efm_swi <- standardise_stand_data(df_efm_swi)

ggplot() +
  geom_line(
    data = df_efm_swi,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    # size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_efm_swi)

# Save stand-level data
saveRDS(df_efm_swi, file = here::here("data/inputs/df_efm_swi_v2026.rds"))

# uholka ----
# Data providers: Jonas Stillhard
# status 1 Tree found (1) or not found (0) during inventory
# status 2 Tree alive (1) or dead (0)
# status 3 Tree standing (1) or lying (0) Boolean
# status 4 Only for dead trees: entire tree (1) or broken tree (snag), (0)

# Management:
# Original info: Plots are considered to be primary in the sense that hey have not been managed in the past.
# Derived info: Since history is known, we assign management_cat = 2 and management_since_census1_yrs = 999.

# Tree-level data
uholka <- read.csv(file.path(raw_data_dir, "uholka/uholka_tree.csv"), sep = ",")

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
    plotID = as.numeric(ifelse(
      nchar(treeID) == 4,
      substr(treeID, 1, 1),
      substr(treeID, 1, 2)
    )),
    plotsize = 0.25, # ha
    biomass = NA,
    country = "Ukraine"
  )

# aggregate data from stand to stand data
df_uholka <- from_tree_data(uholka) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  relocate(census, .after = year) |>
  mutate(
    management = 0,
    management_cat = 2,
    management_since_census1_yrs = 999
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_uholka <- standardise_stand_data(df_uholka)

ggplot() +
  geom_line(
    data = df_uholka,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    # size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_uholka)

# Save stand-level data
saveRDS(df_uholka, file = here::here("data/inputs/df_uholka_v2026.rds"))

# fep gre ----
# Greece
# Data providers: Gavriil Spyroglou and Nikolaos Fyllas
# Management:
# Control-plot labels (including the source typo "Contol") and numeric
# info_management values indicate no management; thinning labels indicate treatment.
# A plot qualifies as unmanaged only if all its records indicate no management.
# management_since_census1_yrs uses years_since_management at the earliest census,
# held constant across the plot's records, including a missing baseline value.

# tree-level data to estimate dominant species
fep_gre_tree <- read.csv(file.path(raw_data_dir, "fp_gre/fp_gre_tree.csv"), sep = ",")

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
  summarise(ba = sum(ba_tree, na.rm = FALSE)) |>
  slice_max(ba, with_ties = FALSE) |>
  ungroup() |>
  select(plotID, year, species)

# stand-level data
fep_gre_stand <- read.csv(file.path(raw_data_dir, "fp_gre/fp_gre_stand.csv"), sep = ",")

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
  mutate(
    biomass = biomass.t.ha.1 * 10^3,
    management = {
      info <- trimws(as.character(info_management))
      case_when(
        info %in%
          c("Control plot (no thinning)", "Contol plot (no thinning)") ~ 0L,
        is.finite(suppressWarnings(as.numeric(info))) ~ 0L,
        grepl("thinning", info, ignore.case = TRUE) ~ 1L,
        TRUE ~ NA_integer_
      )
    }
  ) |>
  left_join(dom_species) |>

  # Use the first available census as the reference for prior management history.
  arrange(plotID, year) |>
  group_by(plotID) |>
  mutate(
    census = match(year, unique(year)),
    management = if (any(management == 1L, na.rm = TRUE)) {
      1L
    } else if (any(is.na(management))) {
      NA_integer_
    } else {
      0L
    },
    management_since_census1_yrs = if (any(is.na(year))) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(as.character(first(years_since_management))))
    }
  ) |>
  ungroup() |>
  relocate(census, .after = year) |> 

  # retain only unmanaged
  filter(management == 0)

# aggregate data from stand
df_fep_gre <- from_stand_data(fep_gre_stand) |>
  mutate(
    management_cat = 1,
    country = "Greece"
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_fep_gre <- standardise_stand_data(df_fep_gre)

ggplot() +
  geom_line(
    data = df_fep_gre,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    # size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_fep_gre)

# Save stand-level data
saveRDS(df_fep_gre, file = here::here("data/inputs/df_fep_gre_v2026.rds"))

# INRAE LESSEM ----
# France
# Data providers: Georges Kunstler
# Management:
# Original: no info on management in the table.
# New information about management included in fp_fra_stand: year_last_intervention and year of first measurement (plot establishment)
# Includes management info: management_since_census1_yrs = first_measurement - date_last_harvest

# Tree- and -stand level data
inrae_lessem_plot <- read.csv(file.path(raw_data_dir, "fp_fra/fp_fra_stand.csv"), sep = ",") |> as_tibble()
inrae_lessem_tree <- read.csv(file.path(raw_data_dir, "fp_fra/fp_fra_tree.csv"), sep = ",") |> as_tibble()
inrae_lessem_species <- read.csv(file.path(raw_data_dir, "fp_fra/fra_sp_code.csv"), sep = ",") |> as_tibble()
inrae_lessem_status <- read.csv(file.path(raw_data_dir, "fp_fra/fra_sta_code.csv"), sep = ",") |> as_tibble()

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
    country = "France",
    management_since_census1_yrs = NA
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
    management_since_census1_yrs = NA,
    management = 0,
    management_cat = 3
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_inrae_lessem <- standardise_stand_data(df_inrae_lessem)

ggplot() +
  geom_line(
    data = df_inrae_lessem,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    # size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_inrae_lessem)

# Save stand-level data
saveRDS(
  df_inrae_lessem,
  file = here::here("data/inputs/df_inrae_lessem_v2026.rds")
)

# EuFoRia plots ----

## bnp ----
# Data providers: Michael Maroschek and Rupert Seidl
# Management:
# Original: In dataset, year_last_management. In Table S1, last management intervention prior census 1 = 6

# Stand-level data by species
bnp <- read.csv(file.path(raw_data_dir, "bnp/euf_bnp.csv"), sep = ",")

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
    management_cat = 2, # Define all bnp plots as management_cat == 2 (primary), see email Rupert Seidl, 16.04.2026
    country = "Germany"
  ) |>
  group_by(plotID) %>%
  mutate(
    management_since_census1_yrs = min(year, na.rm = TRUE) -
      year_last_management
  ) |>
  ungroup()

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_bnp <- standardise_stand_data(df_bnp)

ggplot() +
  geom_line(
    data = df_bnp,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.1,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_bnp)

# Save stand-level data
saveRDS(df_bnp, file = here::here("data/inputs/df_bnp_v2026.rds"))

## czu ----
# Data providers: Miroslav Svoboda
# Data at tree level not used
# Management:
# Original: "Nature reserve established 1955. Low intensity salvage logging of single trees or group of trees till 1990"
# In Table S1. indicated management_since_census1_yrs = 69. Need to clarify with author.

# Tree- and stand-level data by species
czu_tree <- read.csv(file.path(raw_data_dir, "czu/euf_czu_tree.csv"), sep = ",")
# czu_stand <- read.csv(file.path(raw_data_dir, "czu/euf_czu_stand.csv"), sep = ",")

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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_czu <- standardise_stand_data(df_czu)

ggplot() +
  geom_line(
    data = df_czu,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_czu)

# Save stand-level data
saveRDS(df_czu, file = here::here("data/inputs/df_czu_v2026.rds"))

## fvabw ----
# Data providers: Yannek Käber and Lucia Seebach
# Before called forst
# Management:
# Original: In Table S1 indicated that last management was before the time of designation
# Confirm what is tsc variable?
# Assigned management_cat = 3; beni: changed from 2 to 3 in view of very mixed pattern in logN-logQMD plot. It's not convincing that these are all mature stands with density-driven mortality.

# Stand-level data by species
fvabw <- read.csv(file.path(raw_data_dir, "forst/euf_forst.csv"), sep = ",")

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
df_fvabw_utm <- from_species_data(fvabw) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management = 0,
    management_cat = 3,
    management_since_census1_yrs = NA,
    country = "Germany",
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_fvabw <- coords_to_lonlat(df_fvabw_utm, crs = 25832) |>
  standardise_stand_data()

ggplot() +
  geom_point(
    data = df_fvabw,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.5,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_fvabw)

# Save stand-level data
saveRDS(df_fvabw, file = here::here("data/inputs/df_fvabw_v2026.rds"))

## iberbas ----
# Data providers: Tzvetan Zlatanov
# Management:
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = 70

# Tree-level data
iberbas <- read.csv(file.path(raw_data_dir, "iberbas/euf_iberbas.csv"), sep = ",")

# prepare data
iberbas <- iberbas |>
  filter(is.na(dbh) | dbh != 0) |>
  rename(plotsize = plot_size) |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    biomass = NA,
    management_since_census1_yrs = 70,
    country = "Bulgaria"
  )

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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_iberbas <- standardise_stand_data(df_iberbas)

ggplot() +
  geom_line(
    data = df_iberbas,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_iberbas)

# Save stand-level data
saveRDS(df_iberbas, file = here::here("data/inputs/df_iberbas_v2026.rds"))

## unitbv ----
# Data providers: Any Mary Petritan, Cătălin Petritan
# Before called incds
# Management:
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = NA

# Tree-level data
unitbv <- read.csv(file.path(raw_data_dir, "incds/euf_incds.csv"), sep = ",")

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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_unitbv <- standardise_stand_data(df_unitbv)

ggplot() +
  geom_line(
    data = df_unitbv,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_unitbv)

# Save stand-level data
saveRDS(df_unitbv, file = here::here("data/inputs/df_unitbv_v2026.rds"))

## lwf ----
# Data providers: Markus Blaschke
# Management:
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = NA

# Stand-level data
lwf_stand <- read.csv(file.path(raw_data_dir, "lwf/euf_lwf_stand.csv"), sep = ";")

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
lwf_tree <- read.csv(file.path(raw_data_dir, "lwf/euf_lwf_tree.csv"), sep = ",")

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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_lwf <- coords_to_lonlat(df_lwf, crs = 25832) |>
  standardise_stand_data()

ggplot() +
  geom_line(
    data = df_lwf,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_lwf)

# Save stand-level data
saveRDS(df_lwf, file = here::here("data/inputs/df_lwf_v2026.rds"))

## npvbw ----
# Data from Marco Heurich and Isabelle Klein
# Management:
# Use the recorded last-management year for each plot. An intervention within
# the first-to-last census years (inclusive) flags every census as managed.
# The baseline interval is only known when last management predates census 1.

# tree-level data
npvbw <- read.csv(file.path(raw_data_dir, "nbw/euf_nbw_tree.csv"), sep = ",")

# Derive plot history before removing dead trees or aggregating tree records.
npvbw_management <- npvbw |>
  group_by(plotID) |>
  summarise(
    first_year = min(census, na.rm = TRUE),
    last_year = max(census, na.rm = TRUE),
    last_management = if (all(is.na(year_last_management))) NA_real_ else
      max(year_last_management, na.rm = TRUE),
    management = case_when(
      any(year_last_management >= first_year &
            year_last_management <= last_year, na.rm = TRUE) ~ 1L,
      is.na(last_management) ~ NA_integer_,
      TRUE ~ 0L
    ),
    management_since_census1_yrs = if_else(
      management == 0L & last_management < first_year,
      first_year - last_management, NA_real_
    ),
    .groups = "drop"
  ) |>
  select(plotID, management, management_since_census1_yrs)

# prepare data
npvbw <- npvbw |>
  filter(status == "alive") |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    biomass = NA,
    country = "Germany"
  ) |>
  rename(
    plotsize = plot_size,
    year = census
  ) |>
  group_by(plotID) |>
  mutate(
    lat = mean(latitude),
    lon = mean(longitude)
  ) |>
  ungroup() |>
  left_join(npvbw_management, by = "plotID")

# aggregate data from stand to stand data
df_npvbw <- from_tree_data(npvbw) |>
  # The generic tree aggregator drops management and averages the interval;
  # restore the plot-level fields, including explicit NA rather than NaN.
  select(-management_since_census1_yrs) |>
  left_join(npvbw_management, by = "plotID") |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(management_cat = 1)

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_npvbw <- standardise_stand_data(df_npvbw)

ggplot() +
  geom_line(
    data = df_npvbw,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_npvbw)

# Save stand-level data
saveRDS(df_npvbw, file = here::here("data/inputs/df_npvbw_v2026.rds"))

## nfr ----
# Martina Hobi and Harald Bugmann
# Management:
# Original dataset includes info on management.

# Metadata
nfr_metadata <- read.csv(file.path(raw_data_dir, "nfr/raw/NFR_metadata.csv"))

nfr_metadata <- nfr_metadata |>
  mutate(fg = as.character(fg)) |>
  select(fg, lat, long, ele, temp, precip) |>
  distinct(fg, .keep_all = TRUE)

# Plot area
nfr_plot_area <- read.csv(file.path(raw_data_dir, "nfr/raw/NFR_plot_area.csv"))

nfr_plot_area <- nfr_plot_area |>
  mutate(
    fg = as.character(fg),
    FNUM = as.character(FNUM)
  ) |>
  left_join(nfr_metadata)

# Last management intervention
nfr_last_intervention <- read.csv(file.path(raw_data_dir, "nfr/raw/NFR_last_intervention.csv"))

nfr_last_intervention <- nfr_last_intervention |>
  mutate(fg = as.character(fg))

# Stand level data (per plot, year and species)
nfr_stand <- readRDS(file.path(raw_data_dir, "nfr/raw/NFR_stand_data_v2026.rds")) # 291 plots From David Forrester data

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
    qmd = ba_to_qmd(ba, density),
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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_nfr_swi <- standardise_stand_data(df_nfr_swi)

ggplot() +
  geom_point(data = df_nfr_swi, aes(x = log(QMD), y = log(density)), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rgeco::plot_map_simpl() +
  geom_point(aes(lon, lat), data = df_nfr_swi, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_nfr_swi, file = here::here("data/inputs/df_nfr_swi_v2026.rds"))

## nwfva ----
# Data providers: Peter Meyer
# Original dataset includes info on management. Max management_since_census1_yrs = 10
# In Table S1, management_since_census1_yrs = 52. Need to contact data providers.

# Stand-level data
nwfva_stand <- read.csv(file.path(raw_data_dir, "nwfva/euf_nwfva_stand.csv"), sep = ",")

# Tree-level data
nwfva_tree <- read.csv(file.path(raw_data_dir, "nwfva/euf_nwfva_tree.csv"), sep = ",")

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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_nwfva <- coords_to_lonlat(df_nwfva, crs = 31467) |>
  standardise_stand_data()

ggplot() +
  geom_line(
    data = df_nwfva,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_nwfva)

# Save stand-level data
saveRDS(df_nwfva, file = here::here("data/inputs/df_nwfva_v2026.rds"))

## tuzvo ----
# Data providers: Stanislav Kucbel and Peter Jalovia
# Original dataset does not include info on management.
# In Table S1, management_since_census1_yrs = NA, management_cat = 2

# Stand-level data
tuzvo_stand <- read.csv(file.path(raw_data_dir, "tuzvo/euf_tuzvo_stand.csv"), sep = ",")

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
tuzvo_tree <- read.csv(file.path(raw_data_dir, "tuzvo/euf_tuzvo_tree.csv"), sep = ",")

# prepare data
# we use wood density from the BIOMASS pkg to estimate biomass given volume (m3)
WD <- BIOMASS::getWoodDensity(
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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_tuzvo <- standardise_stand_data(df_tuzvo)

ggplot() +
  geom_line(
    data = df_tuzvo,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_tuzvo)

# Save stand-level data
saveRDS(df_tuzvo, file = here::here("data/inputs/df_tuzvo_v2026.rds"))

## ul ----
# Thomas Nagel
# Management:
# Original info: Plots are considered to be old forest with no management.
# Derived info: We assign management_cat = 2 and management_since_census1_yrs = NA

# tree-level data
ul_tree <- read.csv(file.path(raw_data_dir, "ul/euf_ul_tree.csv"), sep = ",")

# metadata
ul_metadata <- read.csv(file.path(raw_data_dir, "ul/euf_ul_meta.csv"), sep = ",")

ul_metadata <- ul_metadata |>
  select(plotid, plot_size) |>
  distinct() |>
  rename(plotsize = plot_size)

ul_coords <- read.csv(file.path(raw_data_dir, "ul/euf_ul_coor.csv"), sep = ",")

ul_coords <- ul_coords |>
  select(plotid, lat, lon)

# prepare data
ul_tree <- ul_tree |>
  left_join(ul_metadata, by = "plotid") |>
  left_join(ul_coords, by = "plotid") |>
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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_ul <- standardise_stand_data(df_ul)

ggplot() +
  geom_line(
    data = df_ul,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_ul)

# Save stand-level data
saveRDS(df_ul, file = here::here("data/inputs/df_ul_v2026.rds"))

## unito ----
# Renzo Motta

# Management:
# Original info: No information
# Derived info: In Table S1, assign to 50 years. Need to confirm with data providers.

# Tree-level data
unito <- read.csv(file.path(raw_data_dir, "unito/euf_unito_tree.csv"), sep = ",")

# Divide plot into grids of different size given the coordinates
# plotsize = 100 * 100 = 1 ha
# New plots = 50 x 50 = 0.25 ha
unito <- unito |>
  filter(lat != 0) |>
  unite(plotID, plotID, LPI_ID, sep = "_")

ggplot(unito) +
  geom_point(aes(lon, lat, col = plotID))

unito_1 <- unito |>
  filter(plotID == "9_1")

ggplot(unito_1) +
  geom_point(aes(lon, lat, col = plotID))

unito_1$grid50 <- interaction(
  cut(unito_1$lon,
    breaks = seq(min(unito_1$lon), max(unito_1$lon), by = 50),
    include.lowest = TRUE
  ),
  cut(unito_1$lat,
    breaks = seq(min(unito_1$lat), max(unito_1$lat), by = 50),
    include.lowest = TRUE
  ),
  sep = "X"
)

unito_1 <- unito_1 |>
  rename(plotIDD = plotID) |>
  group_by(grid50) |>
  mutate(plotID = cur_group_id()) |>
  ungroup() |>
  unite(plotID, plotID, plotIDD, sep = "_")

ggplot(unito_1) +
  geom_point(aes(lon, lat, col = plotID))

# drop ID 5_9_1 since it's at the margins
unito_1 <- unito_1 |>
  filter(plotID != "5_9_1")

# rewrite lat and lon to have info at plot level not tree level
unito_1 <- unito_1 |>
  mutate(
    lon = lon_plot,
    lat = lat_plot
  )

# ggplot(unito_1) +
#   geom_point(aes(lon, lat))

unito_2 <- unito |>
  filter(plotID == "9_2")

ggplot(unito_2) +
  geom_point(aes(lon, lat, col = plotID))

unito_2$grid50 <- interaction(
  cut(unito_2$lon,
    breaks = seq(min(unito_2$lon), max(unito_2$lon), by = 50),
    include.lowest = TRUE
  ),
  cut(unito_2$lat,
    breaks = seq(min(unito_2$lat), max(unito_2$lat), by = 50),
    include.lowest = TRUE
  ),
  sep = "X"
)

unito_2 <- unito_2 |>
  rename(plotIDD = plotID) |>
  group_by(grid50) |>
  mutate(plotID = cur_group_id()) |>
  ungroup() |>
  unite(plotID, plotID, plotIDD, sep = "_")

ggplot(unito_2) +
  geom_point(aes(lon, lat, col = plotID))

# drop ID 5_9_2 since it's at the margins
unito_2 <- unito_2 |>
  filter(plotID != "5_9_2")

# rewrite lat and lon to have info at plot level not tree level
unito_2 <- unito_2 |>
  mutate(
    lon = lon_plot,
    lat = lat_plot
  )

unito_3 <- unito |>
  filter(plotID == "9_3")

ggplot(unito_3) +
  geom_point(aes(lon, lat, col = plotID))

unito_3$grid50 <- interaction(
  cut(unito_3$lon,
    breaks = seq(min(unito_3$lon), max(unito_3$lon), by = 50),
    include.lowest = TRUE
  ),
  cut(unito_3$lat,
    breaks = seq(min(unito_3$lat), max(unito_3$lat), by = 50),
    include.lowest = TRUE
  ),
  sep = "X"
)

unito_3 <- unito_3 |>
  rename(plotIDD = plotID) |>
  group_by(grid50) |>
  mutate(plotID = cur_group_id()) |>
  ungroup() |>
  unite(plotID, plotID, plotIDD, sep = "_")

ggplot(unito_3) +
  geom_point(aes(lon, lat, col = plotID))

# drop ID 5_9_2 since it's at the margins
unito_3 <- unito_3 |>
  filter(plotID != "5_9_3")

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
WD <- BIOMASS::getWoodDensity(
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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_unito <- coords_to_lonlat(df_unito, crs = 32632) |>
  standardise_stand_data()

ggplot() +
  geom_line(
    data = df_unito,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_unito)

# Save stand-level data
saveRDS(df_unito, file = here::here("data/inputs/df_unito_v2026.rds"))

## urk ----
# Data providers: Srdjan Keren and Zbigniew Maciejewski
# Management:
# Original info: year_last_management
# Derived info: In Table S1, assign to NA years. It can be considered a primary/old-growth forests.

# Stand-level data
urk <- read.csv(file.path(raw_data_dir, "urk/euf_urk_stand.csv"), sep = ",")

# aggregate data from stand
df_urk <- from_stand_data(urk) |>
  mutate(
    biomass = NA,
    management = 0,
    management_cat = 2,
    country = "Poland"
  ) |>
  rename(plotsize = area_ha)

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_urk <- standardise_stand_data(df_urk)

ggplot() +
  geom_line(
    data = df_urk,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_urk)

# Save stand-level data
saveRDS(df_urk, file = here::here("data/inputs/df_urk_v2026.rds"))

## wuls ----
# Data providers: Bogdan Brzeziecki

# Management:
# Original info: not management info.
# Derived info: In Table S1, assign management_since_census1_yrs = ca. 20 years. It can be considered a primary/old-growth forests.

# Tree-level data
wuls_tree <- read.csv(file.path(raw_data_dir, "wuls/euf_wuls_tree.csv"), sep = ",")
wuls_census <- read.csv(file.path(raw_data_dir, "wuls/euf_wuls_census.csv"), sep = ",")
wuls_areas <- read.csv(file.path(raw_data_dir, "wuls/euf_wuls_area.csv"), sep = ",")
wuls_coords <- read.csv(file.path(raw_data_dir, "wuls/euf_wuls_coor.csv"), sep = ",")
wuls_species <- read.csv(file.path(raw_data_dir, "wuls/euf_wuls_sp.csv"), sep = ",")

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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_wuls <- standardise_stand_data(df_wuls)

ggplot() +
  geom_line(
    data = df_wuls,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.2,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_wuls)

# Save stand-level data
saveRDS(df_wuls, file = here::here("data/inputs/df_wuls_v2026.rds"))

# ForestGeo plots -----
## Luquillo ----
# Contact: Jess Zimmerman
# Management:
# Not info about management
# See https://forestgeo.si.edu/sites/north-america/luquillo

# Tree-level data
luquillo <- list.files(path = file.path(raw_data_dir, "luquillo"), full.names = TRUE, pattern = "\\.csv$") |>
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
  filter(is.na(dbh) | dbh > 0) |>
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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_luquillo <- standardise_stand_data(df_luquillo)

ggplot() +
  geom_line(
    data = df_luquillo,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.1,
    col = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_luquillo,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.1,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_luquillo)

# Save stand-level data
saveRDS(
  df_luquillo,
  file = here::here("data/inputs/df_luquillo_v2026.rds")
)

## BCI ----
# contact: Salomon Aguilar, technician in BCI
# finally not added to the global dataset due to difficulty to exclude management (communication with  Salomón Aguilar, Research Technician)
# https://forestgeo.si.edu/sites/neotropics/barro-colorado-island

bci <- list.files(path = file.path(raw_data_dir, "bci"), full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

## Prepare data
bci <- bci |>
  select(-species) |>
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
  filter(is.na(dbh) | dbh > 0) |>
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
sp_lianas <- read.csv(file.path(raw_data_dir, "bci/sp_lianas.csv"))

df_bci <- df_bci %>%
  filter(!df_bci$species %in% sp_lianas$species)

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_bci <- standardise_stand_data(df_bci)

ggplot() +
  geom_point(
    data = df_bci,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.5,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_bci)

# Save stand-level data
saveRDS(df_bci, file = here::here("data/inputs/df_bci_v2026.rds"))

## SCBI ----
# Contact: Kristina J. Anderson-Teixeira, William J. McShea, Norman A. Bourg
# https://forestgeo.si.edu/sites/north-america/smithsonian-conservation-biology-institute
# Management:
# No info in dataset
# In Table S1, management_since_census1_yrs = 53, category = 1

scbi <- list.files(path = file.path(raw_data_dir, "scbi"), full.names = TRUE, pattern = "\\.csv$") |>
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
  filter(is.na(dbh) | dbh > 0) |>
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

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_scbi <- standardise_stand_data(df_scbi)

ggplot() +
  geom_point(
    data = df_scbi,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.1,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  ) +
  geom_line(
    data = df_scbi,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.1,
    col = "black",
    inherit.aes = FALSE
  ) +
  theme_classic()

plot_point_map(df_scbi)

# Save stand-level data
saveRDS(df_scbi, file = here::here("data/inputs/df_scbi_v2026.rds"))

## Palanam ----
# Contact: Perry S. Ong

palanam <- list.files(path = file.path(raw_data_dir, "palanam"), full.names = TRUE, pattern = "\\.csv$") |>
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
  filter(is.na(dbh) | dbh > 0) |>
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
  filter(is.na(QMD) | log(QMD) > 1) |>
  mutate(
    management_cat = 2,
    management = 0
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_palanam <- standardise_stand_data(df_palanam)

ggplot() +
  geom_point(
    data = df_palanam,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.5,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

rgeco::plot_map_simpl() +
  geom_point(aes(lon, lat), data = df_palanam, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_palanam, file = here::here("data/inputs/df_palanam_v2026.rds"))

## SERC ----
# Contact: Sean McMahon
# https://forestgeo.si.edu/sites/north-america/smithsonian-environmental-research-center

serc <- list.files(path = file.path(raw_data_dir, "serc"), full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()

unique(serc$Census)
length(unique(serc$Quadrat))

# Divide plot into grids of different size given the coordinates
# 20 x 20 m = 0.04 ha
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
    management_since_census1_yrs = NA,
    country = "USA"
  ) |>
  filter(is.na(dbh) | dbh > 0) |>
  filter(status == "alive") |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

unique(serc$year)
unique(serc$census)

# aggregate data from stand to stand data
df_serc <- from_tree_data(serc) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  filter(is.na(QMD) | log(QMD) > 1) |>
  mutate(
    management_cat = 2,
    management = 0
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_serc <- standardise_stand_data(df_serc)

ggplot() +
  geom_line(
    data = df_serc,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.2,
    col = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_serc,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.2,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  )

plot_point_map(df_serc)

# Save stand-level data
saveRDS(df_serc, file = here::here("data/inputs/df_serc_v2026.rds"))

## Wytham ----
# Data from Yadvinder Malhi
# https://forestgeo.si.edu/sites/europe/wytham-woods
# Management: no information

# tree-level data
wytham <- read.csv(file.path(raw_data_dir, "wytham/wytham.csv"))
wytham_sp <- read.csv(file.path(raw_data_dir, "wytham/wytham_sp.csv"))

# prepare data
wytham <- wytham |>
  pivot_longer(
    cols = c(
      dbh_2008,
      dbh_2010,
      dbh_2016,
      dbh_2021,
      codes_2010,
      codes_2016,
      codes_2021
    ),
    names_to = c(".value", "year"),
    names_sep = "_",
    values_drop_na = FALSE
  ) |>
  left_join(wytham_sp, by = "spID") |>
  rename(status = codes) |>
  mutate(
    lon = -1.3379,
    lat = 51.7743,
    dbh = dbh * 0.1,
    ba_tree = pi * (dbh * 0.01 / 2)^2,
    plotsize = 0.04,
    year = as.integer(year),
    biomass = NA,
    management_since_census1_yrs = NA,
    country = "UK"
  ) |>
  filter(is.na(dbh) | dbh > 0) |>
  filter(status != "D") |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

# aggregate data from stand to stand data
df_wytham <- from_tree_data(wytham) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management_cat = 2,
    management = 0
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_wytham <- standardise_stand_data(df_wytham)

ggplot() +
  geom_line(
    data = df_wytham,
    aes(x = log(QMD), y = log(density), group = plotID),
    alpha = 0.2,
    col = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_wytham,
    aes(x = log(QMD), y = log(density)),
    alpha = 0.2,
    size = 1.5,
    col = "black",
    inherit.aes = FALSE
  ) +
  theme_classic()

plot_point_map(df_wytham)

# Save stand-level data
saveRDS(df_wytham, file = here::here("data/inputs/df_wytham_v2026.rds"))

## Pasho ----
# Contact: Musalmah Nasardin, Yao Tze Leong
# https://forestgeo.si.edu/sites/asia/pasoh

# tree-level data
pasoh_sp <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh_sp_v2026.rds"))

pasoh_sp <- pasoh_sp |>
  separate(Latin, c("first", "second"), " ") |>
  mutate(second = tolower(second)) |>
  unite(Latin, first, second, sep = " ")

pasoh1 <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh1.rds"))
pasoh2 <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh2.rds"))
pasoh3 <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh3.rds"))
pasoh4 <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh4.rds"))
pasoh5 <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh5.rds"))
pasoh6 <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh6.rds"))
pasoh7 <- readRDS(file.path(raw_data_dir, "forestgeo/pasoh/pasoh7.rds"))


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
    management_since_census1_yrs = NA,
    country = "Malaysia"
  ) |>
  filter(DFstatus == "alive") |>
  filter(is.na(dbh) | dbh > 0) |>
  drop_na(PX) |>
  drop_na(PY) |>
  drop_na(year)

# aggregate data from stand to stand data
df_pasoh <- from_tree_data(pasoh) |>
  # create census at the stand level
  group_by(plotID) |>
  mutate(
    census = match(year, unique(year)),
    biomass = ifelse(biomass == 0, NA, biomass)
  ) |>
  ungroup()  |>
  mutate(
    management_cat = 2,
    management = 0
  ) |>
  filter(log(density) > 7.5) # outliers, possible errors

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_pasoh <- standardise_stand_data(df_pasoh)

ggplot() +
  geom_point(data = df_pasoh, aes(x = log(QMD), y = log(density)), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rgeco::plot_map_simpl() +
  geom_point(aes(lon, lat), data = df_pasoh, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_pasoh, file = here::here("data/inputs/df_pasoh.rds"))

## Mudumalai ----
# Contact: Raman Sukumar
# https://forestgeo.si.edu/sites/asia/mudumalai

# Stand-level data
mudumalai <- read.csv(file.path(raw_data_dir, "forestgeo/mudumalai/mudumalai_stand.csv"), sep = ",")


mudumalai <- mudumalai |>
  mutate(density = nindiv / plotsize) |>
  mutate(ba = ba / plotsize)

# aggregate data from stand
df_mudumalai <- from_stand_data(mudumalai) |>
  # create census
  group_by(plotID) |>
  # mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    biomass = 0,
    management_since_census1_yrs = NA,
    country = "India",
    management_cat = 2,
    management = 0
  )

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_mudumalai <- standardise_stand_data(df_mudumalai)

ggplot() +
  geom_point(data = df_mudumalai, aes(x = log(QMD), y = log(density)), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rgeco::plot_map_simpl() +
  geom_point(aes(lon, lat), data = df_mudumalai, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(
  df_mudumalai,
  file = here::here("data/inputs/df_mudumalai.rds")
)

# Forests plots MZ ----
# Data provided by Yadvinder Mahli and Huanyuan Zhang

# 0) Metadata
meta_fp <- read.csv(file.path(raw_data_dir, "forestplots/mahli_zhang/fp_mz_metadata.csv"))


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
multiplots <- read.csv(file.path(raw_data_dir, "forestplots/mahli_zhang/fp_mz_multiplots.csv"))


# select and rename variables
multiplots <- multiplots |>
  select(
    Country,
    Plot.Code,
    TreeID,
    Species,
    Census.Date,
    D0,
    D1,
    D2,
    D3,
    D4
  ) |>
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
    dbh = dbh / 10 # must be in cm
  ) |>
  filter(is.na(dbh) | dbh != 0) |>
  select(-c(D0, D1, D2, D3, D4)) |>
  left_join(meta_fp |> select(-country), by = "plotID") |>
  as_tibble()
unique(multiplots$country)

# 2) Gabon
gabon <- list.files(path = file.path(raw_data_dir, "forestplots/mahli_zhang/gabon"), full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()


# select and rename variables
gabon <- gabon |>
  select(country, plot_code, treeID, species, year, dbh) |>
  rename(plotID = plot_code) |>
  mutate(year = as.integer(year)) |>
  mutate(dbh = dbh / 10) |>  # must be in cm
  filter(is.na(dbh) | dbh != 0) |>
  left_join(
    meta_fp |>
      select(-country),
    by = "plotID"
  )

# 3) Ghana
ghana <- list.files(path = file.path(raw_data_dir, "forestplots/mahli_zhang/ghana"), full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()


# select and rename variables
ghana <- ghana |>
  select(country, plotID, treeID, species, year, dbh) |>
  mutate(year = as.integer(year)) |>
  mutate(dbh = dbh / 10) |> # must be in cm
  filter(is.na(dbh) | dbh != 0) |>
  left_join(
    meta_fp |>
      select(-country),
    by = "plotID"
  )

# 4) Malaysia
malaysia <- list.files(path = file.path(raw_data_dir, "forestplots/mahli_zhang/malaysia"), full.names = TRUE, pattern = "\\.csv$") |>
  lapply(read_csv) |>
  bind_rows()


# select and rename variables
malaysia <- malaysia |>
  select(country, plotID, treeID, species, year, dbh) |>
  mutate(year = as.integer(year)) |>
  mutate(dbh = dbh / 10) |>  # must be in cm
  filter(is.na(dbh) | dbh != 0) |>
  left_join(meta_fp |> 
  select(-country), by = "plotID")

# Join all plots
forestplots <- multiplots |>
  bind_rows(gabon) |>
  bind_rows(ghana) |>
  bind_rows(malaysia) |>
  # calculate basal area
  mutate(
    ba_tree = pi * (dbh * 0.01 / 2)^2,  # ba in m2, dbh in cm
    biomass = NA,
    management_since_census1_yrs = NA
  )

# aggregate data from stand to stand data
df_forestplots <- from_tree_data(forestplots) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup() |>
  mutate(
    management = 0,
    management_cat = 2
    ) |>
  filter(is.na(QMD) | log(QMD) < 4.5) # remove outliers

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_forestplots <- standardise_stand_data(df_forestplots)

ggplot() +
  geom_point(data = df_forestplots, aes(x = log(QMD), y = log(density)), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rgeco::plot_map_simpl() +
  geom_point(aes(lon, lat), data = df_forestplots, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(
  df_forestplots,
  file = here::here("data/inputs/df_forestplots.rds")
)

# Australian plots ----
# Data providers: David Forrester and collaborators:
# Sustainable Timber Tasmania: David Bartlett or Shaun Suitor or Crispen Marunda
# Forestry Corporation of NSW: Mike Sutton
# Queensland data: Michael Ngug or John Neldner
# Victoria: Andrew Clark
# TERN: Open
# Includes management info.
# Data has last management intervention. This shows the last year when some of the trees were labelled as thinned.
# Many plots did not have any trees labelled as being thinned.
# If year_last_management = NA, no management among censuses, management_cat = 2
# If year_last_management is numeric, thinning among censuses, management_cat = 3

# Stand-level data
aus_plots <- readRDS(file.path(raw_data_dir, "fp_aus/fp_aus_stand.RDS"))


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
    year = as.integer(year),
    management = ifelse(year_last_management < min(year, na.rm = TRUE), 0, 1), # beni: added this to be sure
    management_since_census1_yrs = ifelse(
      management == 0,
      min(year, na.rm = TRUE) - year_last_management,
      NA
    ), # beni: added this instead of considering it NA
    # management_cat = ifelse(is.na(year_last_management), 2, 3)  # beni: changed this to cat = 1: known history
    management_cat = ifelse(is.na(year_last_management), 2, 1)
  )

# aggregate data from stand to stand data
df_aus <- from_stand_data(aus_plots) |>
  # create census
  group_by(plotID) |>
  mutate(census = match(year, unique(year))) |>
  ungroup()

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_aus <- standardise_stand_data(df_aus)

ggplot() +
  geom_point(data = df_aus, aes(x = log(QMD), y = log(density)), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rgeco::plot_map_simpl() +
  geom_point(aes(lon, lat), data = df_aus, size = 0.5, alpha = 0.5) +
  theme(legend.position = "bottom")

# save stand-level data
saveRDS(df_aus, file = here::here("data/inputs/df_aus.rds"))

# RAINFOR Amazon Forest Inventory Network ----
# Data from Esquivel-Muelbert, et al. 2025 data published
# Coordinates of the sites from Bennett et al. 2023
# no management info

# Stand-level data
rainfor_plots <- read.csv(file.path(raw_data_dir, "rainfor/esquivel_etal/full_final_dataset.csv"), sep = ",")
rainfor_bennet <- read.csv(file.path(raw_data_dir, "rainfor/bennett_etal/table_s3.csv"), sep = ",")
rainfor_brienen <- read.csv(file.path(raw_data_dir, "rainfor/brienen_etal/table_s1.csv"), sep = ",")


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
rainfor_plots <- rainfor_plots |>
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
    management_since_census1_yrs = NA,
    management = 0,
    management_cat = 2,  # beni: changed from 3 to 2
    species = NA,
    biomass = NA,
    dbh = ba_to_qmd(mean_ba, density = 1) # mean basal area per stem
  ) |> # DBH in centimetres
  left_join(rainfor_coords, by = "plot_code") |>
  rename(plotID = plot_code)

# check how many plots are missing coordinates
rainfor_plots %>%
  filter(is.na(lat) | is.na(lon)) %>%
  distinct(plotID)

# keep only the alphabetic part of the plotID to join the rest of coordinates
# since not all plot numbers are included in the coordinates info, we exclude the numeric part of the plotID
rainfor_prefix <- rainfor_coords |>
  mutate(plot_prefix = str_extract(plot_code, "^[^-]+")) |>
  distinct(plot_prefix, .keep_all = TRUE)

rainfor_plots <- rainfor_plots |>
  mutate(plot_prefix = str_extract(plotID, "^[^-]+")) |>
  left_join(rainfor_prefix |> select(plot_prefix, lon, lat), by = "plot_prefix") |>
  mutate(
    lon = if_else(is.na(lon.x), lon.y, lon.x),
    lat = if_else(is.na(lat.x), lat.y, lat.x)
  ) |>
  select(-lon.x, -lon.y, -lat.x, -lat.y)

# check again how many plots are missing coordinates
rainfor_plots %>%
  filter(is.na(lat) | is.na(lon)) %>%
  distinct(plotID)

# aggregate data from stand
df_rainfor <- from_stand_data(rainfor_plots)

# Standardise stand columns; environmental fields are added in 00_collect_global_drivers.R
df_rainfor <- standardise_stand_data(df_rainfor)

ggplot() +
  geom_point(data = df_rainfor, aes(x = log(QMD), y = log(density)), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE)

rgeco::plot_map_simpl() +
  geom_point(aes(lon, lat), data = df_rainfor, size = 0.5, alpha = 0.5)

# Save stand-level data
saveRDS(df_rainfor, file = here::here("data/inputs/df_rainfor.rds"))
