# Plot WGS84 coordinates (decimal degrees) on a north-up country map.
# Returns a ggplot, so titles, themes and additional layers can be added with +.
# `lon` and `lat` are column names supplied as strings. `padding` is the
# fraction added on each side; `min_span` is the minimum extent in degrees,
# keeping single locations and small clusters in their geographical context.
# Natural Earth boundaries are bundled with rnaturalearthdata (no download).
plot_point_map <- function(data, lon = "lon", lat = "lat",
                           padding = 0.15, min_span = 2,
                           size = 0.5, alpha = 0.5, colour = "black",
                           country_labels = TRUE) {
  packages <- c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.data.frame(data) || !all(c(lon, lat) %in% names(data))) {
    stop("`data` must contain the longitude and latitude columns.", call. = FALSE)
  }
  if (!is.numeric(data[[lon]]) || !is.numeric(data[[lat]])) {
    stop("Longitude and latitude must be numeric decimal degrees.", call. = FALSE)
  }
  if (!is.numeric(padding) || length(padding) != 1L ||
      !is.finite(padding) || padding < 0 ||
      !is.numeric(min_span) || length(min_span) != 1L ||
      !is.finite(min_span) || min_span <= 0) {
    stop("`padding` must be non-negative and `min_span` must be positive.", call. = FALSE)
  }

  points <- data.frame(lon = data[[lon]], lat = data[[lat]])
  valid <- with(points, is.finite(lon) & is.finite(lat) &
                  lon >= -180 & lon <= 180 & lat >= -90 & lat <= 90)
  if (!any(valid)) {
    stop("No valid WGS84 coordinates to plot.", call. = FALSE)
  }
  if (any(!valid)) {
    warning(sum(!valid), " row(s) with missing or invalid coordinates omitted.", call. = FALSE)
  }
  points <- points[valid, , drop = FALSE]

  padded_range <- function(x, limits) {
    extent <- range(x)
    half_width <- max(diff(extent), min_span) * (0.5 + padding)
    pmax(limits[1], pmin(limits[2], mean(extent) + c(-half_width, half_width)))
  }
  xlim <- padded_range(points$lon, c(-180, 180))
  ylim <- padded_range(points$lat, c(-90, 90))
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  window <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
    crs = sf::st_crs(4326)
  ))
  graticule <- sf::st_graticule(window, lon = pretty(xlim), lat = pretty(ylim))

  map <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = world, fill = "grey96", colour = "grey45", linewidth = 0.3) +
    ggplot2::geom_sf(data = graticule, colour = "grey80", linewidth = 0.2)
  if (country_labels) {
    # Place names inside the visible part of each country. A planar
    # equidistant cylindrical CRS avoids spherical centroid ambiguities.
    visible <- sf::st_intersection(
      sf::st_transform(sf::st_geometry(world), 4087),
      sf::st_transform(window, 4087)
    )
    labels <- sf::st_sf(
      name = world$name[attr(visible, "idx")[, 1]],
      geometry = sf::st_transform(sf::st_point_on_surface(visible), 4326)
    )
    map <- map + ggplot2::geom_sf_text(
      data = labels, ggplot2::aes(label = name),
      fun.geometry = identity,
      colour = "grey40", size = 3, check_overlap = TRUE
    )
  }
  map +
    ggplot2::geom_point(data = points, ggplot2::aes(x = lon, y = lat),
                       size = size, alpha = alpha, colour = colour) +
    ggplot2::coord_sf(crs = sf::st_crs(4326), default_crs = sf::st_crs(4326),
                      xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#e6f2f8"),
      panel.grid.major = ggplot2::element_line(colour = "grey80", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank()
    )
}
