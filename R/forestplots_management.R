# Explicit interpretations of fp_mz_metadata.csv. Relative ages are assumed to
# describe history at census 1; the metadata do not state their reference date.
# Values marked lower_bound or approximate must not be read as exact dates.
forestplots_history_lookup <- function() {
  dplyr::bind_rows(
    tibble::tibble(
      history_text = c('intact', 'undisturbed', 'old growth', 'old-growth',
                       'intact vegetation, mild hunting', 'old-growth forest', 'native forest'),
      management_cat = 2L, years = NA_real_, basis = 'provider_old_growth',
      event_first = NA_real_, event_last = NA_real_
    ),
    tibble::tibble(
      history_text = c(
        'reduced impact logging: 50 year post-logging',
        'reduced impact logging: 10+ year post-logging',
        'lightly logged 60 yrs ago',
        'reduced impact logging: 20+ year post-logging',
        'old abandoned agroforest (around 50 years)',
        'late successional (40 years since disturbance)',
        'early successional (40 years since disturbance)',
        'selective logged (40 years since disturbance)',
        'lightly logged > 30 years ago',
        '3-400 year old fire/logging',
        '60 years since disturbance', '20 years since disturbance',
        '12 years since disturbance', '5 years since disturbance',
        '5 years since distubance'),
      management_cat = 1L,
      years = c(50, 10, 60, 20, 50, 40, 40, 40, 30, 300, 60, 20, 12, 5, 5),
      basis = c('reported_duration', 'lower_bound', 'reported_duration',
                'lower_bound', 'approximate', rep('reported_duration', 3),
                'lower_bound', 'approximate_lower_bound_300_to_400',
                rep('reported_duration', 5)),
      event_first = NA_real_, event_last = NA_real_
    ),
    tibble::tibble(
      history_text = c('logged- logging experiment 1994',
        'reduced impact logging: logged around 2016', 'logged 2015', 'logged 2017',
        'burned in 2008', 'burned in 2016', 'burned in 2018', 'burned in 2019',
        'burned in 2018 & 2019'),
      management_cat = 1L, years = NA_real_, basis = 'dated_event',
      event_first = c(1994, 2016, 2015, 2017, 2008, 2016, 2018, 2019, 2018),
      event_last = c(1994, 2016, 2015, 2017, 2008, 2016, 2018, 2019, 2019)
    )
  )
}

classify_forestplots_management <- function(
  meta, 
  first_census_year = NULL,
  keep_forest = TRUE
) {
  normalise <- function(x) tolower(trimws(gsub('[[:space:]]+', ' ', as.character(x))))
  stopifnot(all(c('forest_type', 'anthropogenic_effects') %in% names(meta)))
  if (is.null(first_census_year)) {
    # Only accept a standalone year; one metadata cell contains URLs.
    start <- if ('Monitoring.Start.Date' %in% names(meta))
      trimws(as.character(meta$Monitoring.Start.Date)) else rep(NA_character_, nrow(meta))
    first_census_year <- rep(NA_real_, nrow(meta))
    valid <- !is.na(start) & grepl('^(18|19|20)[0-9]{2}$', start)
    first_census_year[valid] <- as.numeric(start[valid])
  }
  stopifnot(length(first_census_year) == nrow(meta))
  lookup <- forestplots_history_lookup()
  i <- match(normalise(meta$anthropogenic_effects), lookup$history_text)
  meta$management_cat <- dplyr::coalesce(lookup$management_cat[i], 3L)
  meta$management_since_census1_yrs <- lookup$years[i]
  meta$management_history_basis <- dplyr::coalesce(lookup$basis[i], 'unknown_history')
  meta$management_event_first_year <- lookup$event_first[i]
  meta$management_event_last_year <- lookup$event_last[i]
  dated <- which(!is.na(meta$management_event_last_year))
  elapsed <- first_census_year[dated] - meta$management_event_last_year[dated]
  usable <- is.finite(elapsed) & elapsed >= 0
  meta$management_since_census1_yrs[dated[usable]] <- elapsed[usable]
  # A dated event after census 1 does not establish time unmanaged at census 1.
  meta$management_cat[dated[!usable]] <- 3L
  meta$management_history_basis[dated[!usable]] <- 'dated_event_no_pre_census_interval'
  meta$management_history_basis[!is.na(i) & lookup$history_text[i] ==
    'reduced impact logging: logged around 2016'] <- 'approximate_event_year_2016'
  if (keep_forest) {
    forest <- c('lowland tropical forest', 'tropical lowland forest',
                'montain forest', 'mountain forest', 'temperate rainforest')
    meta <- meta[normalise(meta$forest_type) %in% forest, , drop = FALSE]
  }
  meta
}

# Join only forest plots with matching metadata. Re-anchor calendar-dated
# events to actual census 1, rather than metadata monitoring-start dates.
join_forestplots_metadata <- function(trees, meta) {
  years <- trees |>
    dplyr::group_by(plotID) |>
    dplyr::summarise(first = if (any(is.finite(year))) min(year[is.finite(year)]) else NA_real_,
                     last = if (any(is.finite(year))) max(year[is.finite(year)]) else NA_real_,
                     .groups = 'drop')
  stopifnot(!anyDuplicated(meta$plotID))
  i <- match(meta$plotID, years$plotID)
  meta <- classify_forestplots_management(meta, years$first[i], keep_forest = TRUE)
  i <- match(meta$plotID, years$plotID)
  during <- (!is.na(meta$management_event_first_year) &
    meta$management_event_first_year >= years$first[i] &
    meta$management_event_first_year <= years$last[i]) |
    (!is.na(meta$management_event_last_year) &
    meta$management_event_last_year >= years$first[i] &
    meta$management_event_last_year <= years$last[i]
  )

  # Preserve the existing source-selection assumption for other histories;
  # this metadata classification does not independently prove no intervention.
  meta$management <- as.integer(dplyr::coalesce(during, FALSE))
  dplyr::inner_join(trees, dplyr::select(meta, -dplyr::any_of('country')),
                   by = 'plotID', relationship = 'many-to-one')
}

# Generic tree aggregation drops management_cat and management. Restore them
# from the classified tree data, after any parent-to-subplot ID conversion.
forestplots_stand_management <- function(stands, trees) {
  history <- trees |>
    dplyr::distinct(plotID, management, management_cat, management_since_census1_yrs)
  stopifnot(!anyDuplicated(history$plotID))
  stands |>
    dplyr::select(-dplyr::any_of(c('management', 'management_cat', 'management_since_census1_yrs'))) |>
    dplyr::left_join(history, by = 'plotID', relationship = 'many-to-one')
}
