# This script runs the functions to process and plot the data used in the analyses

source(here::here("R/harmonise_dataset_names.R"))

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
    filter(management_since_census1_yrs >= 30 | is.na(management_since_census1_yrs)) |>
    # Filter by min 3 censuses
    group_by(dataset, plotID) |>
    mutate(n_census_unm = n()) |>
    # remove plots with no change in ln(N)
    mutate(var_logdensity = diff(range(logDensity))) |>
    ungroup() |>
    relocate(n_census_unm, .after = n_census) |>
    filter(n_census_unm >= 3) |>
    filter(var_logdensity > 0.001)

  return(data_unm)
}

# filter the upper quantile data points and remove the outliers
data_filter75_fc <- function(data_unm) {
  data_fil <- data_unm |>
    # create QMD bins
    # mutate(QMD_bins = cut(QMD, breaks = 600, include.lowest = TRUE))  |>
    mutate(QMD_bins = cut(QMD, breaks = seq(0, max(QMD), by = 0.25), include.lowest = TRUE)) |>
    # calculate the 75th percentile of density for each qmd
    group_by(QMD_bins) |>
    mutate(upper75 = quantile(density, c(0.75))) |>
    ungroup() |>
    # select upper quantile by QMD bins, i.e., those plots with higher density
    filter(density >= upper75) |>
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

data_filter90_fc <- function(data_unm) {
  data_fil <- data_unm |>
    # create QMD bins
    # mutate(QMD_bins = cut(QMD, breaks = 600, include.lowest = TRUE))  |>
    mutate(QMD_bins = cut(QMD, breaks = seq(0, max(QMD), by = 0.25), include.lowest = TRUE)) |>
    # calculate the 90th percentile of density for each qmd
    group_by(QMD_bins) |>
    mutate(upper90 = quantile(density, c(0.90))) |>
    ungroup() |>
    # select upper quantile by QMD bins, i.e., those plots with higher density
    filter(density >= upper90) |>
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

data_filter55_fc <- function(data_unm) {
  data_fil <- data_unm |>
    # create QMD bins
    # mutate(QMD_bins = cut(QMD, breaks = 600, include.lowest = TRUE))  |>
    mutate(QMD_bins = cut(QMD, breaks = seq(0, max(QMD), by = 0.25), include.lowest = TRUE)) |>
    # calculate the 90th percentile of density for each qmd
    group_by(QMD_bins) |>
    mutate(upper55 = quantile(density, c(0.55))) |>
    ungroup() |>
    # select upper quantile by QMD bins, i.e., those plots with higher density
    filter(density >= upper55) |>
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

# Convert basal area (m2/ha) and stem density (stems/ha) to QMD (cm).
# Both inputs must describe the same census and tree population. They may also
# be totals for the same plot area; use density = 1 for mean basal area per stem.
ba_to_qmd <- function(ba, density) {
  100 * sqrt((4 * ba) / (pi * density))
}

# prepare_stand_data
from_species_data <- function(data) {
  # get name of dataset
  datanm <- harmonise_dataset_names(deparse(substitute(data)))

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
      # Species-level density and basal area are already per-hectare totals.
      # Mean species DBH cannot recover a missing basal-area total.
      density = sum(density, na.rm = FALSE),
      ba = sum(ba, na.rm = FALSE),
      dbh = mean(dbh, na.rm = T),
      biomass = sum(biomass, na.rm = FALSE)
    ) |>
    ungroup() |>

    # calculate QMD
    mutate(
      QMD = ba_to_qmd(ba, density),
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
    filter(is.na(density) | density > 0) |>

    # join dominant species to data
    left_join(dom_species) |>
    ungroup()

  return(data_agg)
}

# prepare_stand_data
from_stand_tree_data <- function(data_tree, data_stand) {

  # get name of dataset
  datanm <- harmonise_dataset_names(deparse(substitute(data_stand)))

  # identify dominant species
  dom_species <- data_tree |>
    group_by(plotID, census, species) |>
    summarise(ba = sum(ba_tree, na.rm = FALSE)) |>
    distinct() |>
    ungroup() |>
    group_by(plotID, census) |>
    slice_max(ba, with_ties = FALSE) |>
    ungroup() |>
    select(plotID, census, species)

  # Aggregate data at stand level
  data_agg <- data_stand |>
    select(-any_of(c("logQMD", "logDensity"))) |>

    # calculate QMD
    mutate(
      QMD = ba_to_qmd(ba, density),
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
  datanm <- harmonise_dataset_names(deparse(substitute(data)))

  # Aggregate data at stand level
  data_agg <- data |>
    select(-any_of(c("logQMD", "logDensity"))) |>

    # calculate QMD
    mutate(
      QMD = ba_to_qmd(ba, density),
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

  return(data_agg)
}

# prepare_stand_data
from_tree_data <- function(data) {
  # get name of dataset
  datanm <- harmonise_dataset_names(deparse(substitute(data)))

  # Preserve supplied individual basal areas (m2). Fill only missing values
  # from that tree's DBH (cm); if both are missing, retain NA in the total.
  data <- data |>
    mutate(ba_tree = coalesce(ba_tree, pi * (dbh / 200)^2, NA_real_))

  # identify dominant species
  dom_species <- data |>
    group_by(plotID, year, species) |>
    reframe(ba = sum(ba_tree, na.rm = FALSE) / plotsize) |>
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
      ba = sum(ba_tree, na.rm = FALSE) / plotsize,
      dbh = mean(dbh, na.rm = T),
      management_since_census1_yrs = mean(management_since_census1_yrs, na.rm = T),
      biomass = sum(biomass, na.rm = FALSE) / plotsize) |>
    ungroup() |>
    distinct() |>
    # calculate QMD
    mutate(
      QMD = ba_to_qmd(ba, density),
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

# Coordinate preparation and combined environmental extraction.
source(here::here("R/extract_plot_environment.R"))

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
