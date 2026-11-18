# Query published polygons only. Raw plot coordinates stay local.
read_aus_arcgis_layer <- function(layer, fields, path, where = "1=1", page_size = 2500L) {
  if (nzchar(path)) return(sf::st_read(path, quiet = TRUE))

  query <- function(...) {
    params <- c(where = where, ...)
    encoded <- vapply(params, utils::URLencode, character(1), reserved = TRUE)
    paste0(layer, "/query?", paste(names(params), encoded, sep = "=", collapse = "&"))
  }
  count <- jsonlite::fromJSON(query(returnCountOnly = "true", f = "json"))$count
  if (length(count) != 1L || count < 1L) stop("Could not count features in ", layer)
  pages <- lapply(seq.int(0L, count - 1L, by = page_size), function(offset) {
    sf::st_read(query(
      outFields = paste(c("OBJECTID", fields), collapse = ","),
      orderByFields = "OBJECTID", resultOffset = offset, resultRecordCount = page_size,
      outSR = "4326", geometryPrecision = "6", f = "geojson"
    ), quiet = TRUE)
  })
  polygons <- do.call(rbind, pages)
  if (nrow(polygons) != count || !"OBJECTID" %in% names(polygons) || anyDuplicated(polygons$OBJECTID)) {
    stop("Incomplete or duplicated polygon download from ", layer)
  }
  polygons
}

read_aus_management_layers <- function() {
  # For a versioned or offline run, set the three *_FILE variables to complete
  # local exports of the corresponding public layers.
  # CAPAD terrestrial areas. Park codes and gazettal fields are defined in
  # https://www.dcceew.gov.au/sites/default/files/documents/capad-2024-tech-spec.pdf
  list(
    protected_areas = read_aus_arcgis_layer(
      "https://gis.environment.gov.au/gispubmap/rest/services/ogc_services/CAPAD/FeatureServer/0",
      c("STATE", "NAME", "TYPE_ABBR", "ENVIRON", "NRS_PA", "IUCN", "GAZ_DATE", "LATEST_GAZ"),
      Sys.getenv("AUS_CAPAD_PARKS_FILE"),
      where = paste(
        "ENVIRON IN ('T','B') AND (",
        "(STATE IN ('NT','WA') AND TYPE_ABBR IN ('CP','NAP','NP','NPC')) OR",
        "(STATE = 'NSW' AND TYPE_ABBR IN ('NP','CCAZ1','NR','SCA','CCAZ3','REP')) OR",
        "(STATE IN ('VIC','TAS') AND NRS_PA = 'Y' AND IUCN = 'II'))"
      )
    ),
    nsw_old_growth = read_aus_arcgis_layer(
      "https://mapprod2.environment.nsw.gov.au/arcgis/rest/services/IFOA/Coastal_IFOA_ESAonly/MapServer/11",
      character(), Sys.getenv("AUS_NSW_OLD_GROWTH_FILE")
    ),
    vic_logging = read_aus_arcgis_layer(
      "https://biod-gis.mapshare.vic.gov.au/arcgis/rest/services/NatureKit/ReferenceMaps/MapServer/28",
      "SEASON", Sys.getenv("AUS_VIC_LOGGING_FILE"), page_size = 10000L
    )
  )
}

aus_plot_points <- function(plots) {
  plots |>
    filter(is.finite(longitude), is.finite(latitude)) |>
    distinct(plotID, longitude, latitude) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_transform(3577)
}

aus_gazettal_year <- function(date) {
  # The REST export stores dates in milliseconds since the Unix epoch.
  if (is.numeric(date)) date <- as.POSIXct(date / 1000, origin = "1970-01-01", tz = "UTC")
  as.integer(substr(as.character(date), 1, 4))
}

match_aus_protected_areas <- function(plots, layers) {
  points <- aus_plot_points(plots)
  parks <- layers$protected_areas
  required <- c("STATE", "NAME", "TYPE_ABBR", "ENVIRON", "NRS_PA", "IUCN", "GAZ_DATE", "LATEST_GAZ")
  if (!all(required %in% names(parks)) || !all(c("NSW", "NT", "WA", "VIC", "TAS") %in% parks$STATE)) {
    stop("Supply a CAPAD layer covering NSW, NT, WA, VIC and TAS with the gazettal and park attributes.")
  }
  parks <- parks |>
    filter(ENVIRON %in% c("T", "B")) |>
    mutate(
      basis = case_when(
        STATE %in% c("NT", "WA") & TYPE_ABBR %in% c("CP", "NAP", "NP", "NPC") ~ "CAPAD park",
        STATE == "NSW" & TYPE_ABBR %in% c("NP", "CCAZ1", "NR", "SCA", "CCAZ3", "REP") ~ "CAPAD park",
        STATE %in% c("VIC", "TAS") & NRS_PA == "Y" & IUCN == "II" ~ "CAPAD protected at first census"
      ),
      original_gazettal = aus_gazettal_year(GAZ_DATE),
      latest_gazettal = aus_gazettal_year(LATEST_GAZ)
    ) |>
    filter(!is.na(basis)) |>
    select(state = STATE, area = NAME, basis, original_gazettal, latest_gazettal) |>
    sf::st_transform(3577) |>
    sf::st_make_valid()
  park_matches <- sf::st_join(points, parks, left = FALSE)

  # Victoria and Tasmania use protection by the first census. The gazettal
  # rule guards against assigning earlier protection to later park additions.
  clean_censuses <- plots |>
    filter(is.finite(year)) |>
    summarise(first_year = floor(min(year)), last_year = floor(max(year)), .by = plotID) |>
    anti_join(
      plots |> filter(!is.na(year_last_management) | TreesPerHectareThinned > 0),
      by = "plotID"
    )
  earlier_protection <- park_matches |>
    filter(basis == "CAPAD protected at first census") |>
    inner_join(clean_censuses, by = "plotID") |>
    filter(original_gazettal < first_year, latest_gazettal < first_year)

  # Victorian Forest Monitoring Program (https://doi.org/10.1071/WF24135).
  # Exclude known logging seasons that start by the last census.
  logging <- layers$vic_logging |>
    mutate(harvest_year = as.integer(substr(SEASON, 1, 4))) |>
    select(harvest_year) |>
    sf::st_transform(3111) |>
    sf::st_make_valid() |>
    sf::st_transform(3577)
  logged <- earlier_protection |>
    filter(state == "VIC") |>
    sf::st_join(logging, left = FALSE) |>
    filter(harvest_year <= last_year) |>
    sf::st_drop_geometry() |>
    distinct(plotID)

  # The NSW EPA old-growth layer provides evidence independently of park tenure.
  # The FCNSW harvest-history polygons do not give plot-specific dates.
  # https://services2.arcgis.com/iCBB4zKDwkw2iwDD/arcgis/rest/services/HFDHarvestHistory/FeatureServer
  nsw_old_growth <- layers$nsw_old_growth |>
    sf::st_transform(3577) |>
    sf::st_make_valid()
  old_growth_matches <- points[lengths(sf::st_intersects(points, nsw_old_growth)) > 0L, ] |>
    mutate(state = "NSW", area = NA_character_, basis = "NSW EPA High Conservation Value Old Growth")

  dplyr::bind_rows(
    sf::st_drop_geometry(filter(park_matches, basis == "CAPAD park")),
    sf::st_drop_geometry(anti_join(earlier_protection, logged, by = "plotID")),
    sf::st_drop_geometry(old_growth_matches)
  ) |>
    tibble::as_tibble() |>
    distinct(plotID, state, area, basis, original_gazettal, latest_gazettal)
}

match_aus_queensland_history <- function(plots) {
  # Graham (2006), Tables 1 and 11, give site names, MGA Zone 55 coordinates
  # and first census years. Match the selected ~0.5 ha plots within 250 m.
  # https://rainforest-crc.jcu.edu.au/publications/permanent_plots_full.pdf
  # Associated inventory data: https://doi.org/10.4225/08/59475C67BE7A4
  sites <- tibble::tribble(
    ~site, ~easting, ~northing, ~first_year, ~category,
    "EP2",   349690,   8103676,       1971,        3L,
    "EP3",   350365,   8110226,       1971,        2L,
    "EP4",   375564,   8120676,       1972,        3L,
    "EP9",   354565,   8107176,       1972,        1L,
    "EP18",  314965,   8172426,       1973,        2L,
    "EP19",  368914,   7954277,       1975,        2L,
    "EP29",  345965,   8059577,       1975,        2L,
    "EP30",  293915,   8199325,       1976,        2L,
    "EP31",  375714,   8061736,       1976,        2L,
    "EP33",  348115,   8088376,       1976,        3L,
    "EP34",  369514,   8075026,       1976,        2L,
    "EP35",  322265,   8191175,       1977,        2L,
    "EP37",  660212,   7649280,       1977,        3L,
    "EP38",  338615,   8073676,       1977,        1L,
    "EP40",  297265,   8199125,       1978,        2L,
    "EP41",  333315,   8215375,       1977,        2L,
    "EP44",  295265,   8205975,       1980,        2L
  ) |>
    sf::st_as_sf(coords = c("easting", "northing"), crs = 28355) |>
    sf::st_transform(3577)
  points <- aus_plot_points(plots)
  first_censuses <- plots |>
    filter(is.finite(year)) |>
    slice_min(year, n = 1L, with_ties = FALSE, by = plotID) |>
    filter(between(plot_size, 0.45, 0.60)) |>
    transmute(plotID, plot_first_year = floor(year))
  csiro <- points |>
    inner_join(first_censuses, by = "plotID") |>
    sf::st_join(sites, join = sf::st_is_within_distance, dist = 250, left = FALSE) |>
    filter(plot_first_year == first_year) |>
    sf::st_drop_geometry() |>
    distinct(plotID, site, category)

  # Laidlaw et al. report the crane location and 680 stems/ha in 2005.
  # https://doi.org/10.1111/j.1442-9993.2007.01739.x
  crane <- sf::st_sfc(sf::st_point(c(145 + 26/60 + 40/3600, -(16 + 6/60 + 20/3600))), crs = 4326) |>
    sf::st_transform(3577)
  crane_census <- plots |>
    filter(floor(year) == 2005, between(plot_size, 0.95, 1.05), round(TreesPerHectareAHC1_2) == 680)
  nearby <- lengths(sf::st_is_within_distance(points, crane, dist = 100)) > 0L
  matches <- dplyr::bind_rows(
    csiro,
    tibble::tibble(plotID = intersect(points$plotID[nearby], crane_census$plotID),
                   site = "Cape Tribulation crane", category = 2L)
  )
  if (anyDuplicated(matches$site) || anyDuplicated(matches$plotID)) {
    stop("Ambiguous match to published Queensland plot histories.")
  }
  # Earlier selective logging at EP9 and EP38 gives a 20-year lower bound.
  # Bauman et al. (2022), https://doi.org/10.1038/s41586-022-04737-7
  matches |>
    tibble::as_tibble() |>
    mutate(published_interval = if_else(category == 1L, 20, NA_real_))
}

classify_aus_management <- function(plots, protected_plot_ids, queensland_history) {
  # Keep the full census date for thinning checks. Management years have
  # annual precision, so only a strictly earlier year establishes prior use.
  plots |>
    group_by(plotID) |>
    summarise(
      first_census = if (all(is.na(year))) NA_real_ else min(year, na.rm = TRUE),
      thinned_after_first = any(year > first_census & TreesPerHectareThinned > 0, na.rm = TRUE),
      dated_management = any(year_last_management >= floor(first_census), na.rm = TRUE),
      has_management_year = any(!is.na(year_last_management)),
      prior_interval = {
        prior_years <- year_last_management[which(year_last_management < floor(first_census))]
        if (length(prior_years)) floor(first_census) - max(prior_years) else NA_real_
      },
      .groups = "drop"
    ) |>
    left_join(queensland_history, by = "plotID", relationship = "many-to-one") |>
    mutate(
      management = case_when(
        thinned_after_first | dated_management ~ 1L,
        is.na(first_census) & has_management_year ~ NA_integer_,
        TRUE ~ 0L
      ),
      management_cat = coalesce(category, case_when(
        !is.na(prior_interval) ~ 1L,
        plotID %in% protected_plot_ids ~ 2L,
        TRUE ~ 3L
      )),
      management_since_census1_yrs = case_when(
        category == 1L ~ published_interval,
        management_cat == 1L ~ prior_interval,
        TRUE ~ NA_real_
      )
    ) |>
    select(plotID, management, management_since_census1_yrs, management_cat)
}
