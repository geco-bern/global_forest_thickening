# Run from the project root: Rscript tests/test_plot_environment.R
source("R/extract_plot_environment.R")

# Repeated censuses and IDs across datasets must not duplicate or misjoin rows.
records <- tibble::tibble(
  row_id = 1:7, dataset = c("a", "a", "b", "b", "c", "c", "mudumalai"),
  plotID = "1", lon = c(8, 8, 9, 10, NA, 181, 8),
  lat = c(47, 47, 48, 49, 47, 47, 47), ai = -999, lai = 99
)
calls <- 0L
mock_extract <- function(locations, paths, include_phosphorus) {
  calls <<- calls + 1L
  stopifnot(nrow(locations) == 3L)
  # Deliberately return another order; location keys must control the join.
  locations[3:1, ] |>
    dplyr::mutate(biomeID = 4L, biome = "test", type = "Forest",
                  ai = lon / 10, noy = 1, nhx = 2, ndep = 3,
                  tavg = lat, CNrt = 10, ORGC = 20, PBR = NA_real_)
}
result <- extract_plot_environment(records, extractor = mock_extract)
stopifnot(
  calls == 1L, identical(result$row_id, records$row_id),
  identical(result$ai, c(.8, .8, .9, 1, NA_real_, NA_real_, .8)),
  !"lai" %in% names(result), !any(grepl("\\.[xy]$", names(result))),
  result$biomeID[7] == 2L,
  result$biome[7] == "Tropical & Subtropical Dry Broadleaf Forests",
  identical(extract_plot_environment(result, extractor = mock_extract), result)
)
empty <- extract_plot_environment(records[FALSE, ], extractor = function(...) stop("Unexpected extraction"))
missing <- extract_plot_environment(records[5:6, ], extractor = function(...) stop("Unexpected extraction"))
stopifnot(nrow(empty) == 0L, nrow(missing) == 2L, all(is.na(missing$ai)))

# Verify every projected CRS used by the constituent sources, with missing rows.
for (epsg in c(25830, 25832, 31467, 32632)) {
  xy <- sf::st_as_sf(data.frame(lon = 9, lat = 48), coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(epsg) |>
    sf::st_coordinates()
  projected <- data.frame(lon = c(xy[1, 1], NA), lat = c(xy[1, 2], NA))
  restored <- coords_to_lonlat(projected, epsg)
  stopifnot(abs(restored$lon[1] - 9) < 1e-6, abs(restored$lat[1] - 48) < 1e-6,
            is.na(restored$lon[2]), is.na(restored$lat[2]))
}

# Small independent rasters check aridity scaling, NA propagation, monthly
# temperature averaging, and extracting coordinates outside the raster extent.
scratch <- tempfile("plot_environment_")
dir.create(scratch)
r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = .2,
                  ymin = 0, ymax = .2, crs = "EPSG:4326")
terra::values(r) <- 10000
r[1] <- NA
ai_path <- file.path(scratch, "ai.tif")
terra::writeRaster(r, ai_path)
locations <- tibble::tibble(lon = c(.025, .175, 1), lat = c(.175, .025, 1))
stopifnot(isTRUE(all.equal(extract_aridity_locations(locations, ai_path), c(NA, 1, NA))))
for (month in 1:12) {
  terra::values(r) <- month
  r[1] <- NA
  terra::writeRaster(r, file.path(scratch, sprintf("wc_tavg_%02d.tif", month)))
}
stopifnot(isTRUE(all.equal(extract_temperature_locations(locations, scratch), c(NA, 6.5, NA))))

# Overlapping polygons must not duplicate rows; Tundra retains its established
# mapping to Boreal, and ocean coordinates remain unclassified.
polygons <- terra::vect(c(
  "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
  "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
  "POLYGON ((2 0, 3 0, 3 1, 2 1, 2 0))",
  "POLYGON ((2 0, 3 0, 3 1, 2 1, 2 0))"
), crs = "EPSG:4326")
polygons$BIOME <- c(11L, 11L, 4L, 5L)
polygon_path <- file.path(scratch, "biomes.gpkg")
terra::writeVector(polygons, polygon_path)
biomes <- extract_biome_locations(tibble::tibble(lon = c(.5, 2.5, 4), lat = .5), polygon_path)
stopifnot(identical(biomes$biomeID, c(6L, NA_integer_, NA_integer_)),
          biomes$biome[1] == "Boreal Forests/Taiga")

# Reproduce the Lamarque decimal-calendar-year time axis. A mistaken selection
# of all years (or terra's parsed 1970 timestamps) changes the expected mean.
lon <- ncdf4::ncdim_def("lon", "degrees_east", c(.25, .75))
lat <- ncdf4::ncdim_def("lat", "degrees_north", c(.25, .75))
time <- ncdf4::ncdim_def("TIME", "year as %Y.%f", 1850:2009)
for (variable in c("noy", "nhx")) {
  definition <- ncdf4::ncvar_def(variable, "gN m-2 yr-1", list(lon, lat, time), -9999)
  path <- file.path(scratch, paste0("ndep_", variable, "_lamarque11cc_historical_halfdeg.nc"))
  nc <- ncdf4::nc_create(path, definition)
  ncdf4::ncvar_put(nc, variable, array(rep(1850:2009, each = 4), dim = c(2, 2, 160)))
  ncdf4::nc_close(nc)
}
ndep <- extract_ndep_locations(tibble::tibble(lon = .25, lat = .25), scratch)
stopifnot(ndep$noy == mean(1970:2009), ndep$nhx == mean(1970:2009),
          ndep$ndep == 2 * mean(1970:2009))
cat("Plot environmental extraction checks passed.\n")
