# FIA COND definitions: STDAGE, STDORGCD, TRTCD1-3, TRTYR1-3,
# STND_COND_CD_PNWRS and STUMP_CD_PNWRS:
# https://research.fs.usda.gov/sites/default/files/2024-05/wo-v9-2_apr2024_ug_fiadb_database_description_nfi.pdf
read_fia_condition_history <- function(dir) {
  files <- list.files(dir, pattern = "_COND\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No FIA COND files found in ", dir)
  columns <- c("UNITCD", "STATECD", "COUNTYCD", "PLOT", "INVYR", "CONDID",
               "COND_STATUS_CD", "RESERVCD", "STDORGCD", "STDAGE",
               "STND_COND_CD_PNWRS", "STUMP_CD_PNWRS", "TRTCD1", "TRTCD2",
               "TRTCD3", "TRTYR1", "TRTYR2", "TRTYR3")
  lapply(files, function(path) {
    readr::read_csv(path, col_select = dplyr::all_of(columns),
      col_types = readr::cols(.default = readr::col_double(),
                             STUMP_CD_PNWRS = readr::col_character()),
      progress = FALSE) |>
      dplyr::mutate(plotID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = "_"))
  }) |>
    dplyr::bind_rows()
}

prepare_fia_management <- function(data, conditions) {
  required <- c("plotID", "INVYR", "CONDID", "COND_STATUS_CD", "RESERVCD",
                "STDORGCD", "STDAGE", "STND_COND_CD_PNWRS", "STUMP_CD_PNWRS",
                "TRTCD1", "TRTCD2", "TRTCD3", "TRTYR1", "TRTYR2", "TRTYR3")
  stopifnot(all(c("plotID", "year") %in% names(data)),
            all(required %in% names(conditions)))
  if (any(!is.finite(data$year) | data$year == 9999)) {
    stop("FIA stand records must have valid census years before management classification.")
  }
  if (!nrow(data)) {
    return(data |>
      dplyr::mutate(management = numeric(), management_cat = integer(),
                     management_since_census1_yrs = numeric()))
  }
  span <- data |>
    dplyr::group_by(plotID) |>
    dplyr::summarise(first_year = min(year), last_year = max(year), .groups = "drop")
  history <- conditions |>
    dplyr::inner_join(span, by = "plotID") |>
    dplyr::filter(is.finite(INVYR), INVYR != 9999, INVYR <= last_year)
  history$no_interval_treatment <- TRUE
  history$no_baseline_treatment <- TRUE
  history$last_prior_treatment <- NA_real_
  for (i in 1:3) {
    code <- history[[paste0("TRTCD", i)]]
    date <- history[[paste0("TRTYR", i)]]
    # Zero, missing and sentinel years cannot date an intervention.
    dated <- is.finite(date) & date >= 1000 & date < 9998 & date <= history$INVYR
    treatment <- !is.na(code) & code > 0
    prior <- treatment & dated & date <= history$first_year
    history$last_prior_treatment <- pmax(history$last_prior_treatment,
                                         ifelse(prior, date, NA_real_), na.rm = TRUE)
    # A baseline treatment predates that census. At a later census, a positive
    # code is disqualifying unless explicitly dated before the first census.
    # Primary NA is unknown; blank secondary slots mean no additional treatment.
    none <- code %in% 0 | (i > 1 & is.na(code))
    history$no_baseline_treatment <- history$no_baseline_treatment & none
    history$no_interval_treatment <- history$no_interval_treatment &
      (none | (treatment & dated & date < history$first_year))
  }

  observed <- history |> dplyr::filter(INVYR >= first_year)
  baseline <- observed |> dplyr::filter(INVYR == first_year)
  if (!nrow(baseline)) {
    return(data[FALSE, ] |>
      dplyr::mutate(management = numeric(), management_cat = integer(),
                     management_since_census1_yrs = numeric()))
  }
  # Date known pre-census management by condition. Older COND records can inform
  # history, but their age or old-growth codes do not replace baseline values.
  prior_dates <- history |>
    dplyr::group_by(plotID, CONDID) |>
    dplyr::summarise(last_prior_treatment = if (all(is.na(last_prior_treatment))) {
      NA_real_
    } else max(last_prior_treatment, na.rm = TRUE), .groups = "drop")
  baseline_conditions <- baseline |>
    dplyr::left_join(prior_dates, by = c("plotID", "CONDID"), suffix = c("_row", "")) |>
    dplyr::mutate(
      # STDAGE is the stand-age estimate, not the field placeholder FLDAGE.
      # Zero denotes nonstocked; negative values and 998/999/9999 are unknown.
      age = dplyr::if_else(is.finite(STDAGE) & STDAGE > 0 &
                            !STDAGE %in% c(998, 999, 9999), STDAGE, NA_real_),
      # Recorded treatment dates take precedence. Stand age is a fallback proxy
      # for elapsed time, not evidence of the date of the last harvest.
      elapsed = dplyr::coalesce(first_year - last_prior_treatment, age)
    ) |>
    dplyr::group_by(plotID, CONDID) |>
    dplyr::summarise(
      oldgrowth = all(STND_COND_CD_PNWRS %in% 7 & STUMP_CD_PNWRS %in% "N" &
                        STDORGCD %in% 0 & no_baseline_treatment) &
        all(is.na(last_prior_treatment)),
      baseline_stumps = any(STUMP_CD_PNWRS %in% "Y"),
      baseline_planted = any(STDORGCD %in% 1),
      elapsed = if (anyNA(elapsed)) NA_real_ else min(elapsed),
      .groups = "drop"
    )

  # Missing condition metadata must not conceal a managed census. All conditions
  # are checked because the stand aggregation pools trees across the whole plot.
  missing_conditions <- data |>
    dplyr::distinct(plotID, year) |>
    dplyr::anti_join(dplyr::distinct(observed, plotID, INVYR),
                     by = c("plotID", "year" = "INVYR"))
  eligible <- observed |>
    dplyr::left_join(baseline_conditions, by = c("plotID", "CONDID")) |>
    dplyr::group_by(plotID) |>
    dplyr::summarise(
      forest = all(COND_STATUS_CD %in% 1 & !is.na(CONDID)),
      reserved = all(RESERVCD %in% 1),
      no_management = all(INVYR == first_year |
        (no_interval_treatment &
           (!STUMP_CD_PNWRS %in% "Y" | baseline_stumps %in% TRUE) &
           (!STDORGCD %in% 1 | baseline_planted %in% TRUE))),
      .groups = "drop"
    ) |>
    dplyr::filter(forest, reserved, no_management) |>
    dplyr::anti_join(missing_conditions, by = "plotID")

  classification <- baseline_conditions |>
    dplyr::semi_join(eligible, by = "plotID") |>
    dplyr::group_by(plotID) |>
    dplyr::summarise(
      management_cat = if (all(oldgrowth)) 2L else if (!anyNA(elapsed)) 1L else 3L,
      # For mixed conditions use the shortest known interval conservatively;
      # if any condition has unknown history, the plot remains category 3.
      management_since_census1_yrs = if (management_cat == 1L) min(elapsed) else NA_real_,
      .groups = "drop"
    )
  data |>
    dplyr::select(-dplyr::any_of(c("management", "management_cat", "management_since_census1_yrs"))) |>
    dplyr::inner_join(classification, by = "plotID") |>
    dplyr::mutate(management = 0)
}
