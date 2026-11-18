# Functions for analysing the tree-size distribution within each plot census.
# Weights represent stems per hectare: 1 / plot area for fixed-area censuses,
# or the supplied expansion factors for inventories with nested sampling areas.

# Weighted diameter statistics -------------------------------------------------
weighted_diameter_quantile <- function(d, w, probs) {

  # Invert the weighted empirical distribution without interpolating diameters.
  order_dbh <- order(d)
  d <- d[order_dbh]
  w <- w[order_dbh]
  cumulative_weight <- cumsum(w) / sum(w)

  purrr::map_dbl(probs, function(prob) {
    d[which(cumulative_weight >= prob)[1]]
  })
}

weighted_diameter_gini <- function(d, w) {

  # The sorted-data formula equals the weighted pairwise Gini, while avoiding
  # a large matrix of all possible tree pairs.
  order_dbh <- order(d)
  d <- d[order_dbh]
  w <- w[order_dbh]
  total_weight <- sum(w)

  sum(w * d * (2 * cumsum(w) - w - total_weight)) /
    (total_weight * sum(w * d))
}

# Summarise one plot census ----------------------------------------------------
size_distribution_metrics <- function(dbh, weight, cutoff = 0) {

  # Assess missingness before applying the diameter cutoff. Small but valid
  # stems are excluded from a truncated distribution, not counted as missing.
  valid_record <- is.finite(dbh) & dbh > 0 & is.finite(weight) & weight > 0
  selected <- valid_record & dbh >= cutoff
  d <- dbh[selected]
  w <- weight[selected]

  if (length(d) == 0L) {
    return(tibble::tibble(
      cutoff_cm = cutoff,
      n_stems = 0L,
      n_effective = 0,
      missing_fraction = mean(!valid_record),
      N = 0,
      BA = 0,
      mean_dbh = NA_real_,
      QMD_tree = NA_real_,
      d10 = NA_real_,
      d50 = NA_real_,
      d90 = NA_real_,
      cv_dbh = NA_real_,
      gini_dbh = NA_real_,
      d10_qmd = NA_real_,
      d50_qmd = NA_real_,
      d90_qmd = NA_real_,
      ba_top10 = NA_real_,
      n_small = 0,
      n_medium = 0,
      n_large = 0,
      n_verylarge = 0,
      ba_large_share = NA_real_
    ))
  }

  # Density and diameter moments use the same expansion weights. The effective
  # stem count identifies censuses dominated by a few heavily weighted records.
  N <- sum(w)
  mean_d <- sum(w * d) / N
  qmd <- sqrt(sum(w * d^2) / N)
  quantiles <- weighted_diameter_quantile(d, w, c(0.1, 0.5, 0.9))

  # Allocate exactly 10% of represented stems to the upper tail, including a
  # fractional contribution from the record crossing the percentile boundary.
  order_desc <- order(d, decreasing = TRUE)
  top_weight <- pmax(
    0,
    pmin(w[order_desc], 0.1 * N - (cumsum(w[order_desc]) - w[order_desc]))
  )

  tibble::tibble(
    cutoff_cm = cutoff,
    n_stems = length(d),
    n_effective = N^2 / sum(w^2),
    missing_fraction = mean(!valid_record),
    N = N,
    BA = sum(w * pi * (d / 200)^2),
    mean_dbh = mean_d,
    QMD_tree = qmd,
    d10 = quantiles[1],
    d50 = quantiles[2],
    d90 = quantiles[3],
    cv_dbh = sqrt(sum(w * (d - mean_d)^2) / N) / mean_d,
    gini_dbh = weighted_diameter_gini(d, w),
    d10_qmd = quantiles[1] / qmd,
    d50_qmd = quantiles[2] / qmd,
    d90_qmd = quantiles[3] / qmd,
    ba_top10 = sum(top_weight * d[order_desc]^2) / sum(w * d^2),
    n_small = sum(w[d < 20]),
    n_medium = sum(w[d >= 20 & d < 40]),
    n_large = sum(w[d >= 40]),
    n_verylarge = sum(w[d >= 60]),
    ba_large_share = sum(w[d >= 40] * d[d >= 40]^2) / sum(w * d^2)
  )
}

# Read the existing source-specific tree preparation ---------------------------
read_prepared_tree_section <- function(
  marker,
  variable,
  raw_data_dir,
  forestplots_metadata = FALSE
) {

  # Reuse the reviewed preparation in 01, stopping before stand aggregation.
  # Evaluate only assignment expressions, so its writes and plots are skipped.
  lines <- readLines("analysis/01_process_data_orig.R", warn = FALSE)
  start <- which(startsWith(lines, marker))

  if (length(start) != 1L) {
    stop("Ambiguous source section: ", marker)
  }

  after <- seq.int(start, length(lines))
  end <- after[grepl("^df_.*<- from_tree_data\\(", lines[after])][1] - 1L

  if (!is.finite(end)) {
    stop("No tree aggregation after section ", marker)
  }

  env <- new.env(parent = globalenv())
  env$raw_data_dir <- raw_data_dir

  if (forestplots_metadata) {
    env$forestplots_t1 <- list()

    # Keep read.csv() here to reproduce the column-name conversion used by 01.
    env$meta_fp <- read.csv(
      file.path(raw_data_dir, "forestplots/mahli_zhang/fp_mz_metadata.csv")
    ) |>
      dplyr::rename(
        plotID = Plot_code,
        country = Country,
        lon = Longitude,
        lat = Latitude,
        elevation = Elevation..m.,
        plotsize = Plot.Size..ha.,
        forest_type = Forest.Type,
        anthropogenic_effects = Anthropogenic.Effects,
        key_personnel = Key.Personnel,
        key_reference = Key.Reference
      ) |>
      classify_forestplots_management()
  }

  for (expr in parse(text = lines[start:end])) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-"))) {

      # Wood density is irrelevant to diameter distributions; avoid the lookup.
      if (identical(expr[[2]], as.name("WD"))) {
        env$WD <- tibble::tibble(meanWD = NA_real_)
      } else {
        eval(expr, envir = env)
      }
    }
  }

  tibble::as_tibble(get(variable, envir = env))
}

# Match tree records to the selected stand sample ------------------------------
summarise_selected_trees <- function(
  trees,
  dataset_name,
  stand,
  cutoffs = c(0, 12.7, 20)
) {
  keys <- stand |>
    dplyr::filter(dataset == dataset_name) |>
    dplyr::select(plotID, year)

  trees <- trees |>
    tibble::as_tibble() |>
    dplyr::mutate(
      plotID = paste0(dataset_name, "_", as.character(plotID)),
      year = as.numeric(year)
    ) |>
    dplyr::inner_join(keys, by = c("plotID", "year"))

  if (nrow(trees) == 0L) {
    return(tibble::tibble())
  }

  # Preserve supplied inventory expansion factors. Fixed-area censuses instead
  # give every recorded stem the same per-hectare weight.
  if (!"weight" %in% names(trees)) {
    trees <- trees |>
      dplyr::mutate(weight = 1 / plotsize)
  }

  trees |>
    dplyr::group_by(plotID, year) |>
    dplyr::group_modify(function(df, key) {
      purrr::map(cutoffs, function(threshold) {
        size_distribution_metrics(df$dbh, df$weight, threshold)
      }) |>
        dplyr::bind_rows()
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(dataset = dataset_name)
}

# Pair the first and last eligible census -------------------------------------
make_size_pairs <- function(censuses, beta_qmd, min_stems = 10L, min_years = 5) {

  # One endpoint change per plot avoids counting successive censuses as
  # independent observations. Original stand N and QMD define the STL position;
  # diameter-distribution metrics may refer to a harmonised diameter subset.
  censuses |>
    dplyr::filter(
      n_stems >= min_stems,
      n_effective + 1e-8 >= min_stems,
      missing_fraction <= 0.01,
      native_agreement,
      is.finite(QMD),
      QMD > 0,
      density > 0
    ) |>
    dplyr::group_by(dataset, biome_major, plotID, cutoff_cm) |>
    dplyr::filter(dplyr::n() >= 2L) |>
    dplyr::arrange(year, .by_group = TRUE) |>
    dplyr::group_modify(function(df, key) {
      if (nrow(df) < 2L) {
        return(tibble::tibble())
      }

      first <- dplyr::slice_head(df, n = 1)
      last <- dplyr::slice_tail(df, n = 1)
      years <- last$year - first$year

      if (years < min_years) {
        return(tibble::tibble())
      }

      beta <- unname(beta_qmd[as.character(key$biome_major)])
      metrics <- c(
        "mean_dbh", "d10", "d50", "d90", "cv_dbh", "gini_dbh",
        "d10_qmd", "d50_qmd", "d90_qmd", "ba_top10", "ba_large_share",
        "N", "BA", "n_small", "n_medium", "n_large", "n_verylarge"
      )

      # A positive shift means movement above the fixed biome-specific slope.
      # It is a plot-level index, not a new estimate of the 90th-quantile line.
      out <- tibble::tibble(
        year0 = first$year,
        year1 = last$year,
        years = years,
        qmd_change_pct_decade = 100 * expm1(
          10 * log(last$QMD / first$QMD) / years
        ),
        density_change_pct_decade = 100 * expm1(
          10 * log(last$density / first$density) / years
        ),
        stl_shift_log_decade = 10 * (
          log(last$density / first$density) - beta * log(last$QMD / first$QMD)
        ) / years,
        oldgrowth = all(
          df$management_cat == 2 |
            (df$management_cat == 1 & df$management_since_census1_yrs >= 300),
          na.rm = FALSE
        )
      )

      # Retain both endpoints alongside the per-decade change for later audits
      # and for the additive small/medium/large density decomposition.
      for (metric in metrics) {
        out[[paste0(metric, "_first")]] <- first[[metric]]
        out[[paste0(metric, "_last")]] <- last[[metric]]
        out[[paste0(metric, "_change_decade")]] <-
          10 * (last[[metric]] - first[[metric]]) / years
      }

      out
    }) |>
    dplyr::ungroup()
}
