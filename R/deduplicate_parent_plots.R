# Remove duplicate source representations of a parent plot, retaining all
# subplots/censuses from the source with the longest first-to-last census span.
# Different parent codes are never merged merely because they are nearby.
deduplicate_parent_plots <- function(
    data, review_distance_km = 10, measurement_tolerance_pct = 2,
    geographic_datasets = c("forestplots", "multiplots", "gabon", "ghana",
                            "malaysia", "rainfor")) {
  stopifnot(all(c("dataset", "plotID", "year", "lon", "lat", "ba", "density", "plotsize") %in% names(data)),
            length(review_distance_km) == 1L, is.finite(review_distance_km),
            review_distance_km > 0, length(measurement_tolerance_pct) == 1L,
            is.finite(measurement_tolerance_pct), measurement_tolerance_pct >= 0,
            !anyNA(data$dataset), !anyNA(data$plotID))
  finite_min <- function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
  finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  finite_median <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
  distance_km <- function(lon1, lat1, lon2, lat2) {
    rad <- pi / 180
    a <- sin((lat2 - lat1) * rad / 2)^2 +
      cos(lat1 * rad) * cos(lat2 * rad) * sin((lon2 - lon1) * rad / 2)^2
    2 * 6371.0088 * asin(sqrt(pmin(1, pmax(0, a))))
  }
  keyed <- data |>
    dplyr::mutate(
      .parent_code = stringr::str_extract(toupper(trimws(plotID)),
                                         "^[A-Z]{3}-[0-9]{2}(?![0-9])"),
      .parent_id = dplyr::coalesce(.parent_code, as.character(plotID))
    )
  # Non-coded identifiers are retained separately; they can still be included
  # in the geographic review for the selected source datasets.
  parents <- keyed |>
    dplyr::filter(!is.na(.parent_code) | dataset %in% geographic_datasets) |>
    dplyr::group_by(dataset, parent_plotID = .parent_id, parent_code = .parent_code) |>
    dplyr::summarise(
      first_year = finite_min(year), last_year = finite_max(year),
      n_years = dplyr::n_distinct(year[is.finite(year)]),
      n_subplots = dplyr::n_distinct(plotID), n_rows = dplyr::n(),
      .groups = "drop"
    )
  locations <- keyed |>
    dplyr::filter(is.finite(lon), is.finite(lat), abs(lon) <= 180, abs(lat) <= 90) |>
    dplyr::distinct(dataset, .parent_id, lon, lat) |>
    dplyr::group_by(dataset, parent_plotID = .parent_id) |>
    dplyr::summarise(lon = finite_median(lon), lat = finite_median(lat), .groups = "drop")
  parents <- parents |>
    dplyr::left_join(locations, by = c("dataset", "parent_plotID")) |>
    dplyr::mutate(record_years = last_year - first_year)

  # Tie-breaks: more distinct census years, then later last census, then dataset
  # name. Never rank by row count, which would favour subplot-divided sources.
  ranked <- parents |>
    dplyr::filter(!is.na(parent_code)) |>
    dplyr::group_by(parent_code) |>
    dplyr::filter(dplyr::n() > 1L) |>
    dplyr::arrange(dplyr::desc(record_years), dplyr::desc(n_years),
                   dplyr::desc(last_year), dataset, .by_group = TRUE) |>
    dplyr::mutate(retained = dplyr::row_number() == 1L,
                  retained_dataset = dplyr::first(dataset),
                  retained_record_years = dplyr::first(record_years),
                  distance_to_retained_km = distance_km(lon, lat,
                    dplyr::first(lon), dplyr::first(lat)),
                  coordinate_conflict = distance_to_retained_km > review_distance_km) |>
    dplyr::ungroup()
  removed <- ranked |>
    dplyr::filter(!retained) |>
    dplyr::select(dataset, .parent_id = parent_plotID)
  # Review different codes across ForestPlots/RAINFOR sources, including those
  # lacking a recognised code. Compute before removal so discrepancies remain
  # visible; mark which source representations survive code deduplication.
  geo <- parents |>
    dplyr::filter(dataset %in% geographic_datasets, is.finite(lon), is.finite(lat)) |>
    dplyr::mutate(.id = dplyr::row_number())
  pairs <- tidyr::expand_grid(i = seq_len(nrow(geo)), j = seq_len(nrow(geo))) |>
    dplyr::filter(i < j, geo$dataset[i] != geo$dataset[j],
                  geo$parent_plotID[i] != geo$parent_plotID[j]) |>
    dplyr::mutate(distance_km = distance_km(geo$lon[i], geo$lat[i],
                                          geo$lon[j], geo$lat[j])) |>
    dplyr::filter(distance_km <= review_distance_km)
  parents <- parents |>
    dplyr::left_join(ranked |> dplyr::select(dataset, parent_plotID, retained),
                     by = c("dataset", "parent_plotID")) |>
    dplyr::mutate(retained = dplyr::coalesce(retained, TRUE))
  geo <- geo |>
    dplyr::left_join(parents |> dplyr::select(dataset, parent_plotID, retained),
                     by = c("dataset", "parent_plotID"))
  candidates <- dplyr::bind_cols(
    pairs |> dplyr::select(distance_km),
    geo[pairs$i, ] |> dplyr::select(-.id) |> dplyr::rename_with(~paste0(.x, "_1")),
    geo[pairs$j, ] |> dplyr::select(-.id) |> dplyr::rename_with(~paste0(.x, "_2"))
  ) |>
    dplyr::arrange(distance_km, dataset_1, parent_plotID_1, dataset_2, parent_plotID_2) |>
    dplyr::rename(retained_after_code_1 = retained_1, retained_after_code_2 = retained_2) |>
    dplyr::mutate(pair_id = dplyr::row_number())

  # Compare per-hectare parent values, weighting subplots by sampled area.
  # Missing values/areas or multiple records for one subplot-year do not count
  # as agreement. Compare exact shared census years, without interpolation.
  parent_mean <- function(x, area, id) {
    if (any(!is.finite(x) | x < 0 | !is.finite(area) | area <= 0) ||
        anyDuplicated(id)) return(NA_real_)
    sum(x * area) / sum(area)
  }
  relative_difference <- function(x, y) {
    # Symmetric percentage difference: neither source is treated as ground truth.
    dplyr::if_else(x == 0 & y == 0, 0, 200 * abs(x - y) / (abs(x) + abs(y)))
  }
  census <- keyed |>
    dplyr::filter(dataset %in% geographic_datasets, is.finite(year)) |>
    dplyr::group_by(dataset, parent_plotID = .parent_id, year) |>
    dplyr::summarise(
      ba = parent_mean(ba, plotsize, plotID),
      density = parent_mean(density, plotsize, plotID), .groups = "drop")
  comparisons <- candidates |>
    dplyr::select(pair_id, dataset_1, parent_plotID_1, dataset_2, parent_plotID_2) |>
    dplyr::inner_join(census, by = c("dataset_1" = "dataset", "parent_plotID_1" = "parent_plotID"),
                      relationship = "many-to-many") |>
    dplyr::rename(ba_1 = ba, density_1 = density) |>
    dplyr::inner_join(census, by = c("dataset_2" = "dataset", "parent_plotID_2" = "parent_plotID", "year"),
                      relationship = "many-to-one") |>
    dplyr::rename(ba_2 = ba, density_2 = density) |>
    dplyr::mutate(ba_difference_pct = relative_difference(ba_1, ba_2),
                  density_difference_pct = relative_difference(density_1, density_2),
                  comparable = is.finite(ba_difference_pct) & is.finite(density_difference_pct))
  agreement <- comparisons |>
    dplyr::group_by(pair_id) |>
    dplyr::summarise(n_shared_years = dplyr::n(), n_comparable_years = sum(comparable),
      max_ba_difference_pct = finite_max(ba_difference_pct),
      max_density_difference_pct = finite_max(density_difference_pct),
      confirmed_duplicate = all(comparable) &&
        all(ba_difference_pct <= measurement_tolerance_pct &
            density_difference_pct <= measurement_tolerance_pct), .groups = "drop")
  candidates <- candidates |>
    dplyr::left_join(agreement, by = "pair_id") |>
    dplyr::mutate(n_shared_years = dplyr::coalesce(n_shared_years, 0L),
                  n_comparable_years = dplyr::coalesce(n_comparable_years, 0L),
                  confirmed_duplicate = dplyr::coalesce(confirmed_duplicate, FALSE))

  # Of the source-parent representations surviving code matching, retain the
  # longest record and remove only its directly confirmed geographic matches.
  # Do not infer duplicate identity through a chain of nearby plots.
  priority <- parents |>
    dplyr::filter(retained) |>
    dplyr::arrange(dplyr::desc(record_years), dplyr::desc(n_years),
                   dplyr::desc(last_year), dataset, parent_plotID)
  candidates$removed_dataset <- NA_character_
  candidates$removed_parent_plotID <- NA_character_
  for (k in seq_len(nrow(priority))) {
    ds <- priority$dataset[k]
    pid <- priority$parent_plotID[k]
    if (any(removed$dataset == ds & removed$.parent_id == pid)) next
    matches <- which(candidates$confirmed_duplicate &
      ((candidates$dataset_1 == ds & candidates$parent_plotID_1 == pid) |
       (candidates$dataset_2 == ds & candidates$parent_plotID_2 == pid)))
    for (m in matches) {
      other <- if (candidates$dataset_1[m] == ds && candidates$parent_plotID_1[m] == pid) 2 else 1
      other_ds <- candidates[[paste0("dataset_", other)]][m]
      other_pid <- candidates[[paste0("parent_plotID_", other)]][m]
      if (any(removed$dataset == other_ds & removed$.parent_id == other_pid)) next
      removed <- dplyr::bind_rows(removed,
        tibble::tibble(dataset = other_ds, .parent_id = other_pid))
      candidates$removed_dataset[m] <- other_ds
      candidates$removed_parent_plotID[m] <- other_pid
    }
  }
  kept <- keyed |>
    dplyr::anti_join(removed, by = c("dataset", ".parent_id")) |>
    dplyr::select(-.parent_code, -.parent_id)
  final_flags <- removed |>
    dplyr::transmute(dataset, parent_plotID = .parent_id, removed = TRUE)
  parents <- parents |>
    dplyr::left_join(final_flags, by = c("dataset", "parent_plotID")) |>
    dplyr::mutate(retained = !dplyr::coalesce(removed, FALSE)) |>
    dplyr::select(-removed)
  ranked <- ranked |>
    dplyr::left_join(parents |> dplyr::select(dataset, parent_plotID, final_retained = retained),
                     by = c("dataset", "parent_plotID"))
  message("Parent-plot deduplication: removed ", nrow(removed),
          " source-parent representations (", nrow(data) - nrow(kept),
          " rows); ", sum(candidates$confirmed_duplicate), " of ", nrow(candidates),
          " geographic candidate pairs meet the BA/density criterion.")
  list(data = kept, parents = parents, code_matches = ranked,
       geographic_candidates = candidates, geographic_year_comparisons = comparisons)
}
