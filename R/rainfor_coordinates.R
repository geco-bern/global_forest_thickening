# Read released coordinates without changing the archived source files.
# Resolution is a property of each release, not the number of printed decimal
# places in an individual value (e.g. -61.4 may be a rounded 0.01-degree value).
read_rainfor_coordinate_sources <- function(raw_data_dir) {
  read_source <- function(path, id, lon, lat, source, doi, resolution) {
    x <- read.csv(file.path(raw_data_dir, "rainfor", path),
                  stringsAsFactors = FALSE)
    dplyr::tibble(
      plotID = trimws(as.character(x[[id]])),
      lon = as.numeric(x[[lon]]),
      lat = as.numeric(x[[lat]]),
      coord_source = source,
      coord_source_doi = doi,
      coord_resolution_deg = resolution
    )
  }

  dplyr::bind_rows(
    read_source(
      "bennett_etal/plot_latlong.csv", "PlotCode", "long", "lat",
      "Bennett 2023", "10.5521/forestplots.net/2023_2", 0.01
    ),
    read_source(
      "brienen_etal/metadata.csv", "PlotCode", "LongitudeDecimal", "LatitudeDecimal",
      "Brienen 2015", "10.5521/forestplots.net/2014_4", 0.1
    ),
    read_source(
      "sullivan_2020/Data package/AllPlots.csv", "PlotCode", "Long", "Lat",
      "Sullivan 2020", "10.5521/forestplots.net/2020_2", 0.01
    ),
    read_source(
      paste0("sullivan_2025/Variation in wood density across tropical South America/",
             "Wood density data.csv"),
      "PlotCode", "Long", "Lat",
      "Sullivan 2025", "10.5521/forestplots.net/2024_4", 0.01
    ),
    read_source(
      "esquivel_muelbert_2020/AmazonTreeMortality_data/d01_plot_rates.csv",
      "Plot.Code", "Longitude.Decimal", "Latitude.Decimal",
      "Esquivel-Muelbert 2020", "10.5521/forestplots.net/2020_6", 0.1
    )
  )
}

combine_rainfor_coordinates <- function(candidates) {
  required <- c("plotID", "lon", "lat", "coord_source",
                "coord_source_doi", "coord_resolution_deg")
  stopifnot(all(required %in% names(candidates)))

  candidates |>
    dplyr::mutate(plotID = trimws(as.character(plotID))) |>
    # Keep coordinate pairs together; do not assemble a location from two
    # incomplete rows. Invalid finer records must not mask valid coarse records.
    dplyr::filter(
      !is.na(plotID), nzchar(plotID),
      is.finite(lon), is.finite(lat),
      abs(lon) <= 180, abs(lat) <= 90,
      is.finite(coord_resolution_deg), coord_resolution_deg > 0
    ) |>
    dplyr::group_by(plotID) |>
    dplyr::filter(coord_resolution_deg == min(coord_resolution_deg)) |>
    dplyr::group_modify(function(.x, .y) {
      pairs <- dplyr::distinct(.x, lon, lat)
      dplyr::tibble(
        lon = mean(pairs$lon),
        lat = mean(pairs$lat),
        coord_resolution_deg = min(.x$coord_resolution_deg),
        coord_source = paste(sort(unique(.x$coord_source)), collapse = "; "),
        coord_source_doi = paste(sort(unique(.x$coord_source_doi)), collapse = "; "),
        coord_n_pairs = nrow(pairs),
        coord_lon_range_deg = diff(range(pairs$lon)),
        coord_lat_range_deg = diff(range(pairs$lat))
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::arrange(plotID)
}
