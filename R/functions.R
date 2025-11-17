# This script runs the functions to process and plot the data used in the analyses

# plot self-thinning line
plot_stl <- function(data) {
  plot <- ggplot(data) +
    geom_point(aes(x = logQMD, y = logDensity), alpha = 0.5, size = 1.5, col = "black", inherit.aes = FALSE) +
    geom_smooth(aes(x = logQMD, y = logDensity), method = "lm", se = F, fullrange = TRUE)
  return(plot)
}

plot_stl_biome <- function(data) {
  plot <- ggplot(data) +
    geom_point(aes(x = logQMD, y = logDensity, col = biome), alpha = 0.5, size = 1.5, inherit.aes = FALSE) +
    geom_smooth(aes(x = logQMD, y = logDensity), method = "lm", se = F, fullrange = TRUE)
  return(plot)
}

# filter unmanaged data
data_unm_fc <- function(data) {
  data_unm <- data |>
    # filter for min qmd
    filter(QMD >= 10) |>
    # filter for forest type
    filter(type == "Forest") |>
    # filter for unmanaged plots
    filter(management == 0) |>
    filter(years_since_management >= 30 | is.na(years_since_management)) |>
    # Filter by min 3 censuses
    group_by(dataset, plotID) |>
    mutate(n_census_unm = n()) |>
    ungroup() |>
    relocate(n_census_unm, .after = n_census) |>
    filter(n_census_unm >= 3)

  return(data_unm)
}

# filter the upper quantile data points and remove the outliers
data_filter_fc <- function(data_unm) {
  data_fil <- data_unm |>
    # filter by upper 70th percentile
    # create QMD bins
    # mutate(QMD_bins = cut(QMD, breaks = 600, include.lowest = TRUE))  |>
    mutate(QMD_bins = cut(QMD, breaks = seq(0, max(QMD), by = 0.25), include.lowest = TRUE)) |>
    # calculate the 75th percentile of density for each qmd
    group_by(QMD_bins) |>
    mutate(upper70 = quantile(density, c(0.70))) |>
    ungroup() |>
    # select upper quantile by QMD bins, i.e., those plots with higher density
    filter(density >= upper70) |>
    # filter removing outliers
    mutate(
      q25_density = quantile(logDensity, c(0.25), na.rm = FALSE),
      q75_density = quantile(logDensity, c(0.75), na.rm = FALSE),
      IQR_density = IQR(logDensity),
      low_density = q25_density - 1.5 * IQR_density,
      upp_density = q75_density + 1.5 * IQR_density,
      q25_QMD = quantile(logQMD, c(0.25), na.rm = FALSE),
      q75_QMD = quantile(logQMD, c(0.75), na.rm = FALSE),
      IQR_QMD = IQR(logQMD),
      low_QMD = q25_QMD - 1.5 * IQR_QMD,
      upp_QMD = q75_QMD + 1.5 * IQR_QMD
    ) |>
    filter(logDensity > low_density & logDensity < upp_density) |>
    filter(logQMD > low_QMD & logQMD < upp_QMD)

  return(data_fil)
}

# plot map biome
plot_map <- function(data) {
  library(rnaturalearth)
  library(rnaturalearthdata)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  map <- ggplot(data = world) +
    geom_sf(color = "grey") +
    geom_point(aes(lon, lat, color = biome), data = data, size = 0.1, alpha = 0.5) +
    theme(legend.position = "bottom") +
    # some layout modifications
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(
      axis.ticks.y.right = element_line(),
      axis.ticks.x.top = element_line(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.position = "bottom"
    )
  return(map)
}

# plot map ai
library(viridis)
library(wesanderson)

plot_map_ai <- function(data) {
  library(rnaturalearth)
  library(rnaturalearthdata)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  map <- ggplot(data = world) +
    geom_sf(color = "grey") +
    geom_point(aes(lon, lat, color = ai), data = data, size = 0.1, alpha = 0.5) +
    scale_color_viridis(direction = -1, option = "magma") +
    # scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous", direction = -1))+
    theme(legend.position = "bottom") +
    # some layout modifications
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(
      axis.ticks.y.right = element_line(),
      axis.ticks.x.top = element_line(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.position = "bottom"
    )
  return(map)
}

# prepare_stand_data
from_species_data <- function(data) {
  # get name of dataset
  datanm <- deparse(substitute(data))

  # identify dominant species
  dom_species <- data |>
    group_by(plotID, year) |>
    slice_max(ba, with_ties = FALSE) |>
    ungroup() |>
    select(plotID, year, species)

  # Aggregate data at stand level
  data_agg <- data |>
    group_by(plotID, plotsize, lon, lat, year) |>
    summarise(
      density = sum(density, na.rm = T),
      ba = sum(ba, na.rm = T),
      dbh = mean(dbh, na.rm = T),
      biomass = sum(biomass, na.rm = T)
    ) |>
    ungroup() |>
    # calculate QMD and logs
    mutate(
      QMD = sqrt(ba / (0.0000785 * density)),
      logDensity = log(density),
      logQMD = log(QMD),
      dataset = datanm
    ) |>
    # calculate period lengths
    arrange(plotID, year) |>
    group_by(plotID) |>
    mutate(period = year - lag(year)) |>
    relocate(period, .after = year) |>
    # calculate basal area increment
    mutate(ba_inc = (ba - lag(ba)) / period) |>
    relocate(ba_inc, .after = ba) |>
    # calculate number of censuses
    mutate(n_census = n()) |>
    ungroup() |>
    filter(density > 0) |>
    # join dominant species to data
    left_join(dom_species) |>
    ungroup()

  return(data_agg)
}

# prepare_stand_data
from_stand_tree_data <- function(data_tree, data_stand) {
  # get name of dataset
  datanm <- deparse(substitute(data_stand))

  # identify dominant species
  dom_species <- data_tree |>
    group_by(plotID, census, species) |>
    summarise(ba = sum(ba_tree)) |>
    distinct() |>
    ungroup() |>
    group_by(plotID, census) |>
    slice_max(ba, with_ties = FALSE) |>
    ungroup() |>
    select(plotID, census, species)

  # Aggregate data at stand level
  data_agg <- data_stand |>
    # calculate QMD and logs
    mutate(
      QMD = sqrt(ba / (0.0000785 * density)),
      logDensity = log(density),
      logQMD = log(QMD),
      dataset = datanm
    ) |>
    # calculate period lengths
    arrange(plotID, year) |>
    group_by(plotID) |>
    mutate(period = year - lag(year)) |>
    relocate(period, .after = year) |>
    # calculate basal area increment
    mutate(ba_inc = (ba - lag(ba)) / period) |>
    relocate(ba_inc, .after = ba) |>
    # calculate number of censuses
    mutate(n_census = n()) |>
    ungroup() |>
    filter(density > 0) |>
    # join dominant species to data
    left_join(dom_species)

  return(data_agg)
}

# prepare_stand_data
from_stand_data <- function(data) {
  # get name of dataset
  datanm <- deparse(substitute(data))

  if (deparse(substitute(data)) == "nfi_norway") {
    # Aggregate data at stand level
    data_agg <- data |>
      # calculate QMD and logs
      rename(QMD = qmd) |>
      mutate(
        logDensity = log(density),
        logQMD = log(QMD),
        dataset = datanm
      ) |>
      # calculate period lengths
      arrange(plotID, year) |>
      group_by(plotID) |>
      mutate(period = year - lag(year)) |>
      relocate(period, .after = year) |>
      # calculate basal area increment
      mutate(ba_inc = (ba - lag(ba)) / period) |>
      relocate(ba_inc, .after = ba) |>
      # calculate number of censuses
      mutate(n_census = n()) |>
      ungroup() |>
      filter(density > 0)
  } else {
    # Aggregate data at stand level
    data_agg <- data |>
      # calculate QMD and logs
      mutate(
        QMD = sqrt(ba / (0.0000785 * density)),
        logDensity = log(density),
        logQMD = log(QMD),
        dataset = datanm
      ) |>
      # calculate period lengths
      arrange(plotID, year) |>
      group_by(plotID) |>
      mutate(period = year - lag(year)) |>
      relocate(period, .after = year) |>
      # calculate basal area increment
      mutate(ba_inc = (ba - lag(ba)) / period) |>
      relocate(ba_inc, .after = ba) |>
      # calculate number of censuses
      mutate(n_census = n()) |>
      ungroup() |>
      filter(density > 0)
  }
  return(data_agg)
}

# prepare_stand_data
from_tree_data <- function(data) {
  # get name of dataset
  datanm <- deparse(substitute(data))

  # identify dominant species
  dom_species <- data |>
    group_by(plotID, year, species) |>
    reframe(ba = sum(ba_tree) / plotsize) |>
    distinct() |>
    ungroup() |>
    group_by(plotID, year) |>
    slice_max(ba, with_ties = FALSE) |> # top_n(1, ba)
    ungroup() |>
    select(plotID, year, species)

  # Aggregate data at stand level
  data_agg <- data |>
    group_by(country, plotID, lon, lat, year, plotsize) |>
    reframe(
      density = n() / plotsize,
      ba = sum(ba_tree) / plotsize,
      dbh = mean(dbh, na.rm = T),
      biomass = sum(biomass, na.rm = F) / plotsize,
      years_since_management = mean(years_since_management, na.rm = F),
      management = mean(management, na.rm = F)
    ) |>
    ungroup() |>
    distinct() |>
    # calculate QMD and logs
    mutate(
      QMD = sqrt(ba / (0.0000785 * density)),
      logDensity = log(density),
      logQMD = log(QMD),
      dataset = datanm
    ) |>
    # calculate period lengths
    arrange(plotID, year) |>
    group_by(plotID) |>
    mutate(period = year - lag(year)) |>
    relocate(period, .after = year) |>
    # calculate basal area increment
    mutate(ba_inc = (ba - lag(ba)) / period) |>
    relocate(ba_inc, .after = ba) |>
    # calculate number of censuses
    mutate(n_census = n()) |>
    ungroup() |>
    filter(density > 0) |>
    # join dominant species to data
    left_join(dom_species) |>
    ungroup()

  return(data_agg)
}

# add coords and biomes
biomes_coords_utm <- function(data) {
  # reprojecting data to lon/lat (epsg:4326)
  coords <- data |>
    rename(
      lonUTM = lon,
      latUTM = lat
    ) |>
    select(lonUTM, latUTM) |>
    as.data.frame()

  # sf transformation of df sets coordinate columns and projection epsg
  if (deparse(substitute(data)) == "data_nfi_spain") {
    coords_sf <- sf::st_as_sf(
      coords,
      coords = c("lonUTM", "latUTM"),
      crs = "epsg:25830"
    ) # 25832 #32630 #EPSG:25830 epsg:32630
  }

  if (deparse(substitute(data)) == "data_nwfva") {
    coords_sf <- sf::st_as_sf(
      coords,
      coords = c("lonUTM", "latUTM"),
      crs = "epsg:31467"
    ) # 25832 #32630
  }

  if (deparse(substitute(data)) == "data_unito") {
    coords_sf <- sf::st_as_sf(
      coords,
      coords = c("lonUTM", "latUTM"),
      crs = "epsg:32632"
    )
  }

  if (deparse(substitute(data)) == "data_forst" |
    deparse(substitute(data)) == "data_lwf") {
    coords_sf <- sf::st_as_sf(
      coords,
      coords = c("lonUTM", "latUTM"),
      crs = "epsg:25832"
    ) # 25832 #32630 #31467
  }

  # transformation to lon/lat coordinate system
  coords_sf_lonlat <- coords_sf |>
    st_transform(crs = "epsg:4326")

  # return to simple data frame
  coords_lonlat <- coords_sf_lonlat |>
    dplyr::as_tibble() |>
    dplyr::mutate(
      lon = st_coordinates(geometry)[, "X"],
      lat = st_coordinates(geometry)[, "Y"]
    ) |>
    dplyr::select(-geometry)

  # join latlon to stand data
  data_agg <- data |>
    rename(
      lonUTM = lon,
      latUTM = lat
    ) |>
    cbind(coords_lonlat)

  # Identify biomes for each plot using coordinates
  # WWF Ecoregions data
  wwf_eco <- terra::vect(
    file.path(here::here(), "/data/wwf/wwf_terr_ecos.shp")
  )

  # extract coordinates from data
  coordinates <- data_agg |>
    select(lon, lat)

  # extract map values for coordinates
  coordinates$wwf <- terra::extract(wwf_eco, coordinates)

  # Unnest the data.frames
  coordinates <- coordinates |>
    unnest(cols = c(wwf)) |>
    # select variables
    select(lon, lat, BIOME) |>
    rename(biomeID = BIOME) |>
    # take unique rows of coordinates a
    distinct()

  # read legend for wwf biomes
  legend_wwf <- read.csv(file.path(here::here(), "/data/wwf/legend.csv"))

  # Join biomes to data
  data_agg <- data_agg |>
    left_join(coordinates, by = c("lon", "lat")) |>
    left_join(legend_wwf, by = "biomeID") |>
    drop_na(biomeID) |>
    drop_na(biome) |>
    mutate(
      biome = ifelse(biome == "Tundra", "Boreal Forests/Taiga", biome),
      biomeID = ifelse(biomeID == 11, 6, biomeID)
    ) |>
    distinct() |>
    filter(type == "Forest")

  # Select variables to join datasets
  data_agg <- data_agg |>
    select(dataset, country, plotID, plotsize, lon, lat, census, year, management, years_since_management, ba, dbh, QMD, density, logQMD, logDensity, period, ba_inc, biomass, n_census, species, biomeID, biome, type) %>%
    relocate(any_of(c("dataset", "country", "plotID", "plotsize", "lon", "lat", "census", "year", "management", "years_since_management", "ba", "dbh", "QMD", "density", "logQMD", "logDensity", "period", "ba_inc", "biomass", "n_census", "species", "biomeID", "biome", "type"))) |>
    mutate(
      plotID = as.character(plotID),
      census = as.character(census)
    )

  return(data_agg)
}

# add coords and biomes
biomes_coords_latlon <- function(data) {
  # Identify biomes for each plot using coordinates
  # WWF Ecoregions data
  wwf_eco <- terra::vect(
    file.path(here::here(), "/data/wwf/wwf_terr_ecos.shp")
  )

  # extract coordinates from data
  coordinates <- data |>
    select(lon, lat)

  # extract map values for "coordinates
  coordinates$wwf <- terra::extract(wwf_eco, coordinates)

  # Unnest the data.frames
  coordinates <- coordinates |>
    as_tibble() |>
    unnest(cols = c(wwf)) |>
    # select variables
    select(lon, lat, BIOME) |>
    rename(biomeID = BIOME) |>
    # take unique rows of coordinates a
    distinct()

  # read legend for wwf biomes
  legend_wwf <- read.csv(file.path(here::here(), "/data/wwf/legend.csv"))

  # Join biomes to data
  data_agg <- data |>
    left_join(coordinates, by = c("lon", "lat")) |>
    left_join(legend_wwf) |>
    drop_na(biomeID) |>
    drop_na(biome) |>
    mutate(
      biome = ifelse(biome == "Tundra", "Boreal Forests/Taiga", biome),
      biomeID = ifelse(biomeID == 11, 6, biomeID)
    ) |>
    distinct() |>
    filter(type == "Forest")

  # Select variables to join datasets
  data_agg <- data_agg |>
    select(dataset, country, plotID, plotsize, lon, lat, census, year, management, years_since_management, ba, dbh, QMD, density, logQMD, logDensity, period, ba_inc, biomass, n_census, species, biomeID, biome, type) %>%
    relocate(any_of(c("dataset", "country", "plotID", "plotsize", "lon", "lat", "census", "year", "management", "years_since_management", "ba", "dbh", "QMD", "density", "logQMD", "logDensity", "period", "ba_inc", "biomass", "n_census", "species", "biomeID", "biome", "type"))) |>
    mutate(
      plotID = as.character(plotID),
      census = as.character(census)
    )

  return(data_agg)
}

# add coords and aridity index (later called Moisture Index for being more intuitive)
# global map from (Zomer et al. 2022)

ai_coords_latlon <- function(data) {
  # read aridity index data
  rasta <- rast("~/data/aridityindex_zomer_2022/Global-AI_ET0_v3_annual/ai_v3_yr.tif")

  # values are provided as integers, multiplied by 1e-4
  values(rasta) <- values(rasta) * 1e-4

  # aggregate to 0.1 deg
  rasta_agg <- aggregate(
    rasta,
    fact = res(rasta)[1]^-1 * 0.1,
    fun = "mean"
  )

  # extract coordinates from data
  coordinates <- data |>
    select(lon, lat)

  # extract map values for coordinates
  coordinates$ai <- terra::extract(rasta_agg, coordinates, xy = FALSE, ID = FALSE)

  # Unnest the data.frames
  coordinates <- coordinates |>
    as_tibble() |>
    unnest(cols = c(ai)) |>
    rename(ai = awi_pm_sr_yr)

  # Join aridity index to data
  data_agg <- data |>
    bind_cols(coordinates[, 3])

  return(data_agg)
}

# add coords and LAI

lai_coords_latlon <- function(data) {
  coordinates0 <- data |>
    select(lon, lat) |>
    as.data.frame()

  for (y in 2000:2023) {
    nc_all <- list.files(path = paste0("~/data/lai/", y), full.names = TRUE, pattern = "\\.nc$") |>
      lapply(rast)

    # extract coordinates from data
    data_coord <- data |>
      select(lon, lat) |>
      as.data.frame()

    for (i in 1:length(nc_all)) {
      # extract map values for coordinates
      i_lai <- terra::extract(nc_all[[i]]$lai, data_coord[, 1:2], xy = FALSE, ID = FALSE)
      colnames(i_lai) <- paste0("lai", i)
      data_coord <- cbind(data_coord, i_lai)
    }

    lai_max <- data_coord[, 3:length(data_coord)] |>
      summarise(across(everything(), sum, na.rm = T)) |>
      which.max() |>
      as.numeric()

    data_coord <- data_coord[, lai_max + 2] |>
      as.data.frame()
    colnames(data_coord) <- paste0("lai_", y)

    coordinates0 <- cbind(coordinates0, data_coord)
  }
  coordinates <- coordinates0 |>
    mutate(lai = rowMeans(select(coordinates0, lai_2000:lai_2023))) |>
    select(lai)

  # Join aridity index to data
  data_agg <- data |>
    bind_cols(coordinates)

  return(data_agg)
}

# add coords and ndep
# This reads nitrogen deposition from global annual maps by Lamarque et al. (2011).
# This provides annual data separately for NHx and NOy (reduced and oxidised) in gN m−2 yr−1
# from a global map provided at half-degree resolution and covering years 1860-2009.

ndep_coords_latlon <- function(data) {
  siteinfo <- data |>
    select(plotID, lon, lat) |>
    rename(sitename = plotID) |>
    mutate(
      year_start = 1970,
      year_end = 2009
    ) |>
    distinct()

  df_ndep <- ingest(
    siteinfo,
    source    = "ndep",
    timescale = "y",
    dir       = "~/data/ndep_lamarque/",
    verbose   = FALSE
  )

  df_ndep <- df_ndep |>
    unnest(cols = data) |>
    group_by(sitename) |>
    # Calculate the mean for the last 40 years available (1970 to 2009)
    summarise(noy = mean(noy, na.rm = T), nhx = mean(nhx, na.rm = T)) |>
    rename(plotID = sitename)

  # Join CN ratio to data
  data_agg <- data |>
    left_join(df_ndep) |>
    mutate(ndep = noy + nhx)

  return(data_agg)
}

# add coords and C:N ratio
cn_coords_latlon <- function(data) {
  siteinfo <- data |>
    select(plotID, lon, lat) |>
    rename(sitename = plotID) |>
    distinct()

  settings_wise <- get_settings_wise(varnam = c("CNrt"), layer = 1:3)

  df_wise <- ingest(
    siteinfo,
    source    = "wise",
    settings  = settings_wise,
    dir       = "~/data/soil/wise"
  )

  df_wise <- df_wise |>
    unnest(cols = data) |>
    rename(plotID = sitename)

  # Join CN ratio to data
  data_agg <- data |>
    left_join(df_wise)

  return(data_agg)
}

# add coords and phosphorus
# The amount of P using the Bray1 method ppm of weight
phos_coords_latlon <- function(data) {
  siteinfo <- data |>
    select(plotID, lon, lat) |>
    rename(sitename = plotID) |>
    distinct()

  settings_gsde <- list(varnam = c("PBR"), layer = 1:3)

  df_gsde <- ingest(
    siteinfo,
    source    = "gsde",
    settings  = settings_gsde,
    dir       = "~/data/soil/shangguan"
  )

  df_gsde <- df_gsde |>
    unnest(cols = data) |>
    rename(plotID = sitename)

  # Join PBR ratio to data
  data_agg <- data |>
    left_join(df_gsde)

  return(data_agg)
}

# add coords and ORGC - Organic carbon content (g kg-1)
orgc_coords_latlon <- function(data) {
  siteinfo <- data |>
    select(plotID, lon, lat) |>
    rename(sitename = plotID) |>
    distinct()

  settings_wise <- get_settings_wise(varnam = c("ORGC"), layer = 1:3)

  df_wise <- ingest(
    siteinfo,
    source    = "wise",
    settings  = settings_wise,
    dir       = "~/data/soil/wise"
  )

  df_wise <- df_wise |>
    unnest(cols = data) |>
    rename(plotID = sitename)

  # Join CN ratio to data
  data_agg <- data |>
    left_join(df_wise)

  return(data_agg)
}

# global C sink

csink_global_v1 <- function(global_drivers_row, data_all) {
  # 1*. Sample the row and het QMD and the random variables
  random_row_index <- sample(1:nrow(data_all), 1)
  random_values <- data_all[random_row_index, c("QMD", "species", "dataset", "plotID")]
  QMDj <- random_values$QMD

  # 2. Estimate mean N given QMDj using the LMM  (lnN ~ lnQMD + year + 1|...) and two years (for example year = 2000, 2001) => N0, N1.
  # Single new data point t0
  point_t0 <- random_values |>
    mutate(
      logQMD = log(QMD),
      year = 2000,
      ai = global_drivers_row$ai,
      PBR = global_drivers_row$PBR,
      ORGC = global_drivers_row$ORGC,
      ndep = global_drivers_row$ndep
    )

  # Predict for the new point, including random effects
  logN0_mean <- predict(fit1, newdata = point_t0, re.form = NULL, allow.new.levels = TRUE)
  N0_mean <- exp(logN0_mean)

  # Single new data point t1
  point_t1 <- point_t0 |> mutate(year = 2001)
  # Predict for the new point, including random effects
  logN1_mean <- predict(fit1, newdata = point_t1, re.form = NULL, allow.new.levels = TRUE)
  N1_mean <- exp(logN1_mean)

  # 3. Estimate the biomass change, given the QMDj, N0 and N1 as: dB = ak * QMDj^2 * (N1 - N0),
  # where ak is a sampled value from the fitted a, considering its standard error and a normal distribution.
  a_mean <- summary(fit2)$coefficient[1, 1]
  a_sd <- summary(fit2)$coefficient[1, 2]

  ak <- rnorm(1, a_mean, a_sd)
  dB <- ak * QMDj^2 * (N1_mean - N0_mean)

  out <- tibble(
    QMDj = QMDj,
    N0 = N0_mean,
    N1 = N1_mean,
    a = ak,
    dB_Mg_ha = dB * 10^-3,
    lon = global_drivers_row$lon,
    lat = global_drivers_row$lat,
    area_ha = global_drivers_row$area_ha
  )

  return(out)
}
