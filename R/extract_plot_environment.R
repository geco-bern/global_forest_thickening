# Prepare common stand columns without extracting or filtering global fields.
standardise_stand_data <- function(data) {
  data |>
    dplyr::select(dplyr::any_of(c(
      "dataset", "country", "plotID", "plotsize", "lon", "lat", "census", "year",
      "management", "management_cat", "management_since_census1_yrs", "ba", "dbh",
      "QMD", "density", "period", "ba_inc", "biomass", "n_census", "species"
    ))) |>
    dplyr::mutate(plotID = as.character(plotID), census = as.character(census))
}

# Raw sources with projected coordinates must be transformed before compilation.
# Missing coordinates remain missing; no stand records are removed here.
coords_to_lonlat <- function(data, crs) {
  valid <- is.finite(data$lon) & is.finite(data$lat)
  coordinates <- matrix(NA_real_, nrow(data), 2)
  if (any(valid)) {
    points <- sf::st_as_sf(data[valid, c("lon", "lat")],
                           coords = c("lon", "lat"), crs = crs)
    coordinates[valid, ] <- sf::st_coordinates(sf::st_transform(points, 4326))
  }
  data$lon <- coordinates[, 1]
  data$lat <- coordinates[, 2]
  data
}

plot_environment_paths <- function() {
  list(
    biomes = here::here("data/wwf/wwf_terr_ecos.shp"),
    ai = path.expand("~/data/aridityindex_zomer_2022/Global-AI_ET0_v3_annual/ai_v3_yr.tif"),
    ndep = path.expand("~/data/ndep_lamarque"),
    wise = path.expand("~/data/soil/wise"),
    gsde = path.expand("~/data/soil/shangguan"),
    worldclim = path.expand("~/data/archive/worldclim_fick_2017/data")
  )
}

# WWF codes/names and the existing analysis classification. Keeping the lookup
# here avoids depending on the untracked data/wwf/legend.csv file.
wwf_biome_legend <- function() {
  tibble::tibble(
    biomeID = 1:14,
    biome = c(
      "Tropical & Subtropical Moist Broadleaf Forests",
      "Tropical & Subtropical Dry Broadleaf Forests",
      "Tropical & Subtropical Coniferous Forests",
      "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests",
      "Boreal Forests/Taiga",
      "Tropical & Subtropical Grasslands, Savannas & Shrublands",
      "Temperate Grasslands, Savannas & Shrublands",
      "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands",
      "Tundra", "Mediterranean Forests, Woodlands & Scrub",
      "Deserts & Xeric Shrublands", "Mangroves"
    ),
    type = c(rep("Forest", 6), rep("Grassland", 4),
             "Forest", "Forest", "Desert", "Mangroves")
  )
}

extract_biome_locations <- function(locations, path) {
  ecoregions <- terra::vect(path)
  points <- terra::vect(as.data.frame(locations[, c("lon", "lat")]),
                        geom = c("lon", "lat"), crs = "EPSG:4326")
  matches <- terra::extract(ecoregions[, "BIOME"],
                            terra::project(points, terra::crs(ecoregions))) |>
    dplyr::distinct(id.y, BIOME)
  # Duplicate polygons with the same biome are harmless. Conflicting matches
  # at polygon boundaries remain unclassified rather than duplicating records.
  matches <- matches |>
    dplyr::group_by(id.y) |>
    dplyr::summarise(biomeID = {
      ids <- unique(BIOME[!is.na(BIOME)])
      if (length(ids) == 1L) as.integer(ids) else NA_integer_
    }, .groups = "drop")
  tibble::tibble(biomeID = matches$biomeID[match(seq_len(nrow(locations)), matches$id.y)]) |>
    dplyr::mutate(biomeID = dplyr::if_else(biomeID == 11L, 6L, biomeID)) |>
    dplyr::left_join(wwf_biome_legend(), by = "biomeID")
}

extract_raster_locations <- function(raster, locations) {
  points <- terra::vect(as.data.frame(locations[, c("lon", "lat")]),
                        geom = c("lon", "lat"), crs = "EPSG:4326")
  terra::extract(raster, terra::project(points, terra::crs(raster)), ID = FALSE)
}

mean_available <- function(values) {
  result <- rowMeans(as.matrix(values), na.rm = TRUE)
  result[is.nan(result)] <- NA_real_
  result
}

extract_aridity_locations <- function(locations, path) {
  raster <- terra::rast(path)
  # Aggregate the source integers first, without loading the full global
  # raster into R. Scaling after extraction is equivalent to scaling first.
  factor <- as.integer(round(0.1 / terra::res(raster)))
  if (any(factor < 1L)) stop("Aridity input must have resolution <= 0.1 degree.")
  raster <- terra::aggregate(raster, fact = factor, fun = "mean", na.rm = FALSE)
  extract_raster_locations(raster, locations)[[1]] * 1e-4
}

extract_ndep_locations <- function(locations, dir) {
  result <- lapply(c("noy", "nhx"), function(variable) {
    path <- file.path(dir, paste0("ndep_", variable, "_lamarque11cc_historical_halfdeg.nc"))
    raster <- terra::rast(path)
    nc <- ncdf4::nc_open(path)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    time_dimension <- nc$dim[[which(tolower(names(nc$dim)) == "time")]]
    dates <- terra::time(raster)
    # These files encode TIME as decimal calendar years ("year as %Y.%f"),
    # which terra can otherwise interpret as seconds after 1970-01-01.
    years <- if (grepl("^year as", time_dimension$units)) {
      as.integer(floor(time_dimension$vals))
    } else if (inherits(dates, c("Date", "POSIXt"))) {
      as.integer(format(dates, "%Y"))
    } else {
      as.integer(dates)
    }
    if (!all(1970:2009 %in% years)) stop("Nitrogen deposition must cover 1970-2009.")
    mean_available(extract_raster_locations(raster[[which(years %in% 1970:2009)]], locations))
  })
  tibble::tibble(noy = result[[1]], nhx = result[[2]], ndep = noy + nhx)
}

extract_temperature_locations <- function(locations, dir) {
  files <- list.files(dir, pattern = "tavg.*\\.tif$", full.names = TRUE)
  months <- as.integer(sub(".*_([0-9]+)\\.tif$", "\\1", basename(files)))
  if (length(files) != 12L || !setequal(months, 1:12)) {
    stop("WorldClim directory must contain exactly 12 monthly tavg GeoTIFFs: ", dir)
  }
  mean_available(extract_raster_locations(terra::rast(files[order(months)]), locations))
}

# One ingest call for both WISE variables, using the existing top-three-layer
# depth weighting. GSDE retains its existing variable, layers and directory.
extract_soil_locations <- function(locations, paths, include_phosphorus = TRUE) {
  siteinfo <- locations |>
    dplyr::transmute(sitename = as.character(dplyr::row_number()), lon, lat)
  read_soil <- function(source, settings, dir, variables) {
    values <- ingestr::ingest(siteinfo, source = source, settings = settings, dir = dir) |>
      tidyr::unnest(cols = data) |>
      dplyr::ungroup() |>
      dplyr::select(sitename, dplyr::all_of(variables)) |>
      dplyr::mutate(sitename = as.character(sitename))
    if (anyDuplicated(values$sitename)) stop("Multiple soil summaries returned for one location.")
    dplyr::left_join(siteinfo[, "sitename", drop = FALSE], values, by = "sitename") |>
      dplyr::select(-sitename)
  }
  values <- read_soil("wise", ingestr::get_settings_wise(c("CNrt", "ORGC"), layer = 1:3),
                      paths$wise, c("CNrt", "ORGC"))
  if (include_phosphorus) {
    values$PBR <- read_soil("gsde", list(varnam = "PBR", layer = 1:3), paths$gsde, "PBR")$PBR
  } else {
    values$PBR <- NA_real_
  }
  values
}

extract_environment_at_locations <- function(locations, paths, include_phosphorus = TRUE) {
  message("Extracting biomes for ", nrow(locations), " unique plot locations")
  result <- dplyr::bind_cols(locations, extract_biome_locations(locations, paths$biomes))
  message("Extracting aridity index")
  result$ai <- extract_aridity_locations(locations, paths$ai)
  message("Extracting mean nitrogen deposition (1970-2009)")
  result <- dplyr::bind_cols(result, extract_ndep_locations(locations, paths$ndep))
  message("Extracting mean annual temperature (WorldClim, 1970-2000)")
  result$tavg <- extract_temperature_locations(locations, paths$worldclim)
  message("Extracting soil C:N, organic C", if (include_phosphorus) " and phosphorus")
  dplyr::bind_cols(result, extract_soil_locations(locations, paths, include_phosphorus))
}

# Extract once per distinct coordinate pair, then join to every plot-census.
# Coordinates, not plot IDs, are the key: IDs can repeat between datasets and
# some plots have corrected coordinates between censuses. Keep all input rows;
# the common analysis filters in 02_create_filters.R select forest biomes.
extract_plot_environment <- function(data, paths = plot_environment_paths(),
                                     include_phosphorus = TRUE,
                                     extractor = extract_environment_at_locations) {
  stopifnot(all(c("lon", "lat") %in% names(data)),
            is.numeric(data$lon), is.numeric(data$lat))
  fields <- c("biomeID", "biome", "type", "ai", "noy", "nhx", "ndep", "tavg", "CNrt", "ORGC", "PBR")
  locations <- data |>
    dplyr::select(lon, lat) |>
    dplyr::distinct() |>
    dplyr::filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90)
  if (nrow(locations)) {
    values <- extractor(locations, paths, include_phosphorus)
    stopifnot(all(c("lon", "lat", fields) %in% names(values)),
              nrow(values) == nrow(locations),
              !anyDuplicated(values[c("lon", "lat")]))
    if (nrow(dplyr::anti_join(locations, values, by = c("lon", "lat")))) {
      stop("Environmental extraction returned different coordinates.")
    }
    values <- dplyr::select(values, lon, lat, dplyr::all_of(fields))
  } else {
    values <- tibble::tibble(lon = numeric(), lat = numeric(), biomeID = integer(),
      biome = character(), type = character(), ai = numeric(), noy = numeric(),
      nhx = numeric(), ndep = numeric(), tavg = numeric(), CNrt = numeric(),
      ORGC = numeric(), PBR = numeric())
  }
  # Remove previous extracts (including retired LAI) to prevent stale values or
  # .x/.y columns when the stage is rerun on an already enriched df_all.rds.
  result <- data |>
    dplyr::ungroup() |>
    dplyr::select(-dplyr::any_of(c(fields, "lai"))) |>
    dplyr::left_join(values, by = c("lon", "lat"), na_matches = "never")
  
  if ("dataset" %in% names(result)) {
    mudumalai <- !is.na(result$dataset) & result$dataset == "mudumalai"
    result$biomeID[mudumalai] <- 2L
    result$biome[mudumalai] <- "Tropical & Subtropical Dry Broadleaf Forests"
    result$type[mudumalai] <- "Forest"
  }
  stopifnot(nrow(result) == nrow(data))
  result
}
