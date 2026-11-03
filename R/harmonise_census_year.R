# Most frequent non-missing year. Ties resolve to the earliest year.
# An entirely undated plot-census remains NA.
modal_census_year <- function(year) {
  year <- year[!is.na(year)]
  if (!length(year)) return(NA_real_)
  years <- sort(unique(year))
  years[which.max(tabulate(match(year, years)))]
}

# Defaults match the renamed columns in 00_raw_data.R. To keep an audit trail,
# add year_original = year before calling. Does not merge different censuses.
harmonise_census_year <- function(data, plot = "plotID", census = "census",
                                 year = "year") {
  stopifnot(all(c(plot, census, year) %in% names(data)),
            is.numeric(data[[year]]))
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(plot, census)))) |>
    dplyr::mutate("{year}" := modal_census_year(.data[[year]])) |>
    dplyr::ungroup()
}
