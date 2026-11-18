# Run from the project root: Rscript tests/test_australian_management.R
suppressPackageStartupMessages(library(dplyr))
source("R/australian_management.R")
expressions <- as.list(parse("analysis/01_process_data_orig.R"))
assignment <- function(name) {
  Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, expressions)
}
no_site_matches <- tibble(plotID = character(), category = integer(), published_interval = numeric())
prepare <- function(data, protected_plot_ids = character(), site_matches = no_site_matches) {
  env <- new.env()
  env$aus_plots <- data
  env$aus_management <- classify_aus_management(data, protected_plot_ids, site_matches)
  eval(assignment("aus_plots")[[2L]], env)
  env$aus_plots
}
stand <- tibble(
  plotID = rep(c("early", "late", "same_year", "during", "after", "mixed"), each = 2),
  year = c(1950, 1960, rep(c(2000, 2010), 5)),
  year_last_management = c(1940, 1940, 1980, 1980, 2000, 2000,
                           2005, 2005, 2020, 2020, 1980, 2005),
  longitude = 150, latitude = -30, TreesPerHectareAHC1_2 = 100,
  TreesPerHectareThinned = 0,
  DBHqAHC1_2_cm = 20, BasalAreaAHC1_2_m2perha = 10,
  plot_size = 1, dominant_species = "test", AbovegroundAHC1_2_Mgperha = 50
)
result <- prepare(stand)
stopifnot(
  identical(result$management, rep(c(0L, 0L, 1L, 1L, 1L, 1L), each = 2)),
  identical(result$management_since_census1_yrs, c(10, 10, 20, 20, rep(NA_real_, 6), 20, 20)),
  identical(result$management_cat, rep(c(1L, 1L, 3L, 3L, 3L, 1L), each = 2)),
  !is_grouped_df(result),
  nrow(result) == nrow(stand),
  identical(result, prepare(stand[nrow(stand):1, ])[nrow(stand):1, ])
)
missing_dates <- bind_rows(
  stand[1:2, ] |> mutate(plotID = "missing", year_last_management = NA_real_),
  stand[1:2, ] |> mutate(plotID = "partly_missing", year_last_management = c(NA, 1940)),
  stand[1:2, ] |> mutate(plotID = "managed_missing", year_last_management = c(NA, 1955))
)
missing_result <- prepare(missing_dates)
stopifnot(
  identical(missing_result$management, c(0L, 0L, 0L, 0L, 1L, 1L)),
  identical(missing_result$management_since_census1_yrs,
            c(NA_real_, NA_real_, 10, 10, NA_real_, NA_real_)),
  identical(missing_result$management_cat, rep(c(3L, 1L, 3L), each = 2))
)
decimal <- stand[1:2, ] |>
  mutate(year = c(2000.25, 2005.5), year_last_management = 2000)
decimal_result <- prepare(decimal)
stopifnot(all(decimal_result$management == 1L),
          all(is.na(decimal_result$management_since_census1_yrs)),
          all(decimal_result$management_cat == 3L))
earlier_decimal <- prepare(decimal |> mutate(year_last_management = 1998))
stopifnot(all(earlier_decimal$management == 0L),
          all(earlier_decimal$management_since_census1_yrs == 2))
# Conflicting earlier dates use the latest event for every census, including
# rows where the date is absent. All three management fields must be constant.
conflicting <- stand[1:2, ] |>
  mutate(year_last_management = c(1930, 1940))
conflicting_result <- prepare(conflicting)
stopifnot(all(conflicting_result$management_since_census1_yrs == 10))

thinning <- bind_rows(
  decimal |> mutate(plotID = "first_only", year_last_management = NA_real_,
                     TreesPerHectareThinned = c(4, NA_real_)),
  decimal |> mutate(plotID = "within_year", year = c(2000.25, 2000.5),
                     year_last_management = NA_real_, TreesPerHectareThinned = c(NA_real_, 4)),
  decimal |> mutate(plotID = "dated_then_thinned", year_last_management = 1980,
                     TreesPerHectareThinned = c(0, 4))
)
thinning_result <- prepare(thinning)
stopifnot(
  identical(thinning_result$management, rep(c(0L, 1L, 1L), each = 2)),
  identical(thinning_result$management_cat, rep(c(3L, 3L, 1L), each = 2)),
  identical(thinning_result$management_since_census1_yrs, c(rep(NA_real_, 4), 20, 20))
)

# Match published Queensland site names to arbitrary compilation IDs.
qld_locations <- tibble(
  plotID = c("field_a", "field_b", "field_c", "field_d", "wrong_census"),
  easting = c(349690, 354565, 314965, 338615, 350365),
  northing = c(8103676, 8107176, 8172426, 8073676, 8110226),
  year = c(1971, 1972, 1973, 1977, 1980)
) |>
  sf::st_as_sf(coords = c("easting", "northing"), crs = 28355) |>
  sf::st_transform(4326)
qld_coords <- sf::st_coordinates(qld_locations)
qld <- bind_rows(
  qld_locations |>
    sf::st_drop_geometry() |>
    mutate(longitude = qld_coords[, 1L], latitude = qld_coords[, 2L],
           plot_size = 0.5, TreesPerHectareAHC1_2 = 100),
  tibble(plotID = "crane", longitude = 145 + 26/60 + 40/3600,
         latitude = -(16 + 6/60 + 20/3600), year = 2005.25,
         plot_size = 1, TreesPerHectareAHC1_2 = 680)
) |>
  mutate(year_last_management = NA_real_, TreesPerHectareThinned = 0)
qld_history <- match_aus_queensland_history(qld)
qld_result <- classify_aus_management(qld, "field_a", qld_history) |>
  arrange(match(plotID, qld$plotID))
stopifnot(
  setequal(qld_history$site, c("EP2", "EP9", "EP18", "EP38", "Cape Tribulation crane")),
  identical(qld_result$management_cat, c(3L, 1L, 2L, 1L, 3L, 2L)),
  identical(qld_result$management_since_census1_yrs, c(NA_real_, 20, NA_real_, 20, NA_real_, NA_real_))
)

box <- function(lon) sf::st_polygon(list(rbind(
  c(lon - 0.1, -30.1), c(lon + 0.1, -30.1), c(lon + 0.1, -29.9),
  c(lon - 0.1, -29.9), c(lon - 0.1, -30.1)
)))
park_table <- tibble::tribble(
  ~lon, ~STATE, ~TYPE_ABBR, ~NAME,              ~LATEST_GAZ,
   130, "NT",   "NAP",     "NT park",          "1980-01-01",
   131, "WA",   "CP",      "WA park",          "1980-01-01",
   132, "NSW",  "SCA",     "NSW conservation", "2015-01-01",
   133, "NSW",  "FLR",     "Flora reserve",    "1980-01-01",
   135, "VIC",  "NP",      "Early park",       "1990-01-01",
   136, "VIC",  "NP",      "Park addition",    "2004-01-01",
   137, "VIC",  "NP",      "Logged park",      "1990-01-01",
   138, "VIC",  "NP",      "Before logging",   "1990-01-01",
   139, "VIC",  "NP",      "Thinned park",     "1990-01-01",
   140, "TAS",  "NP",      "Early park",       "1990-01-01",
   141, "TAS",  "NP",      "Later protection", "2004-01-01",
   142, "NSW",  "NP",      "Dated treatment",  "1980-01-01"
) |>
  mutate(ENVIRON = "T", NRS_PA = "Y", IUCN = "II", GAZ_DATE = as.Date("1980-01-01"))
parks <- sf::st_sf(park_table, geometry = sf::st_sfc(lapply(park_table$lon, box), crs = 4326))
layers <- list(
  protected_areas = parks,
  nsw_old_growth = sf::st_sfc(box(134), crs = 4326),
  vic_logging = sf::st_sf(SEASON = c("2005/06", "2006/07"),
                         geometry = sf::st_sfc(box(137), box(138), crs = 4326))
)
spatial_cases <- tibble::tribble(
  ~plotID,        ~longitude, ~management, ~category,
  "nt_park",             130,          0L,        2L,
  "wa_park",             131,          0L,        2L,
  "nsw_park",            132,          0L,        2L,
  "outside_nsw_park",    133,          0L,        3L,
  "nsw_old_growth",      134,          1L,        2L,
  "early_vic",           135,          0L,        2L,
  "late_vic",            136,          0L,        3L,
  "vic_logged",          137,          0L,        3L,
  "vic_logged_later",    138,          0L,        2L,
  "vic_thinned",         139,          0L,        3L,
  "tas_early",           140,          0L,        2L,
  "tas_late",            141,          0L,        3L,
  "nsw_dated",           142,          0L,        1L,
  "missing_location",     NA,          0L,        3L
)
spatial_plots <- decimal[rep(1:2, nrow(spatial_cases)), ] |>
  mutate(
    plotID = rep(spatial_cases$plotID, each = 2),
    longitude = rep(spatial_cases$longitude, each = 2),
    year_last_management = if_else(plotID == "nsw_dated", 1980, NA_real_),
    TreesPerHectareThinned = if_else(
      (plotID == "nsw_old_growth" & year == 2005.5) | (plotID == "vic_thinned" & year == 2000.25),
      4, 0
    )
  )
protected_matches <- match_aus_protected_areas(spatial_plots, layers)
spatial_result <- prepare(spatial_plots, protected_matches$plotID)
stopifnot(
  identical(spatial_result$management_cat, rep(spatial_cases$category, each = 2)),
  identical(spatial_result$management, rep(spatial_cases$management, each = 2)),
  identical(spatial_result$management_since_census1_yrs[spatial_result$plotID == "nsw_dated"], c(20, 20)),
  n_distinct(protected_matches$basis) == 3L
)

for (data in list(result, missing_result, decimal_result, earlier_decimal,
                  conflicting_result, thinning_result, spatial_result)) {
  labels <- data |> distinct(plotID, management, management_since_census1_yrs, management_cat)
  stopifnot(nrow(labels) == n_distinct(data$plotID))
}
cat("Australian published histories, protected-area and logging overlays, and thinning checks passed.\n")
