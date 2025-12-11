calc_lqmm_byqmdbin <- function(df, binwidth = 0.25, breaks = NA) {
  # -- 1. make sure at least 300 points are on average in each bin
  # nbins <- min(round(length(df$logQMD[!is.na(df$logQMD)]) / 300), 10)
  # if (identical(breaks, NA)){
  #  breaks <- pretty(df$logQMD, n = nbins)
  # }
  # bin_labels <- breaks[1:length(breaks)-1] + (breaks[2] - breaks[1])/2

  # Create QMD breaks at fixed intervals (default 0.25)
  if (identical(breaks, NA)) {
    qmin <- floor(min(df$logQMD, na.rm = TRUE) / binwidth) * binwidth
    qmax <- ceiling(max(df$logQMD, na.rm = TRUE) / binwidth) * binwidth
    breaks <- seq(qmin, qmax, by = binwidth)
  }

  # Bin midpoints
  bin_labels <- head(breaks, -1) + diff(breaks) / 2

  # Assign bins
  df <- df |>
    mutate(bin_lqmm = cut(
      logQMD,
      breaks = breaks,
      include.lowest = TRUE,
      labels = bin_labels
    ))

  # Nest and fit models, including logQMD_sc
  df <- df |>
    group_by(bin_lqmm) |>
    nest() |>
    mutate(nvals = purrr::map_int(data, nrow)) |>
    filter(nvals >= 30) |> # require enough data in a bin

    # scale QMD *within each bin*
    mutate(data = purrr::map(data, ~ .x |>
      mutate(logQMD_sc = scale(logQMD)[, 1]))) |>

    # LQMM fit per bin
    mutate(mod_lqmm = purrr::map(data, ~ lqmm(
      logDensity ~ year_sc + logQMD_sc,
      random = ~1,
      group = plotID,
      tau = 0.90,
      type = "normal",
      data = .,
      control = lqmmControl(
        LP_max_iter = 6000, # inner loop iterations
        UP_max_iter = 6000,
        LP_tol_ll   = 1e-05, # inner loop tolerance
        startQR     = TRUE
      )
    ))) |>
    # extract summaries safely
    mutate(sum_lqmm = purrr::map(mod_lqmm, ~ try(summary(.)))) |>
    mutate(failed = purrr::map_lgl(sum_lqmm, determine_failed)) |>
    filter(!failed) |>
    # extract only YEAR effect
    mutate(
      pval = purrr::map_dbl(sum_lqmm, ~ get_pval(.)),
      coef_year = purrr::map_dbl(sum_lqmm, ~ get_coef_year(.)),
      coef_year_upper = purrr::map_dbl(sum_lqmm, ~ get_coef_year_upper(.)),
      coef_year_lower = purrr::map_dbl(sum_lqmm, ~ get_coef_year_lower(.)),
      upwardshift = coef_year > 0 & pval < 0.01
    )

  return(list(df = df, breaks = breaks))
}

# Robust extractor: finds row in tTable that matches a pattern (e.g., "year")
extract_row <- function(sum, pattern = "year") {
  rn <- rownames(sum$tTable)
  idx <- grep(pattern, rn, ignore.case = TRUE)
  if (length(idx) == 0) {
    return(NA_integer_)
  } # no match
  idx[1]
}

get_pval <- function(sum) {
  idx <- extract_row(sum, "year")
  if (is.na(idx)) {
    return(NA_real_)
  }
  sum$tTable[idx, "Pr(>|t|)"]
}

get_coef_year <- function(sum) {
  idx <- extract_row(sum, "year")
  if (is.na(idx)) {
    return(NA_real_)
  }
  sum$tTable[idx, "Value"]
}

get_coef_year_upper <- function(sum) {
  idx <- extract_row(sum, "year")
  if (is.na(idx)) {
    return(NA_real_)
  }
  sum$tTable[idx, "upper bound"]
}

get_coef_year_lower <- function(sum) {
  idx <- extract_row(sum, "year")
  if (is.na(idx)) {
    return(NA_real_)
  }
  sum$tTable[idx, "lower bound"]
}

determine_failed <- function(out) {
  # summary() returned an error
  if (inherits(out, "try-error")) {
    return(TRUE)
  }

  # summary is not a list
  if (!is.list(out)) {
    return(TRUE)
  }

  # tTable missing
  if (is.null(out$tTable)) {
    return(TRUE)
  }

  # no rows (rare but possible)
  if (nrow(out$tTable) == 0) {
    return(TRUE)
  }

  return(FALSE)
}
