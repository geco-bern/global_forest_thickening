# Decompose paired density changes into contributions from diameter classes.
# These are changes in stem numbers, not a demographic decomposition into
# recruitment, mortality and growth, or additive parts of a fitted quantile.

# Validate the additive density accounting -------------------------------------
density_contribution_data <- function(pairs, tolerance = 1e-8) {
  components <- c("N", "n_small", "n_medium", "n_large", "n_verylarge")

  # Small (<20), medium (20–40) and large (>=40 cm) form an exhaustive partition.
  # Very large (>=60 cm) is a subset of large, used only in sensitivity checks.
  for (suffix in c("first", "last", "change_decade")) {
    values <- pairs |>
      dplyr::select(dplyr::all_of(paste0(components, "_", suffix))) |>
      as.matrix()

    if (any(!is.finite(values))) {
      stop("Non-finite density components.")
    }

    discrepancy <- values[, 1] - rowSums(values[, 2:4, drop = FALSE])

    if (any(abs(discrepancy) > tolerance * pmax(1, abs(values[, 1])))) {
      stop("Small, medium and large classes do not sum to total stem density.")
    }

    if (suffix != "change_decade" &&
        (any(values < -tolerance) || any(values[, 5] > values[, 4] + tolerance))) {
      stop("Invalid or overlapping size-class counts.")
    }
  }

  if (any(!is.finite(pairs$years) | pairs$years <= 0) || any(pairs$N_first <= 0)) {
    stop("Invalid census interval or baseline density.")
  }

  if (anyDuplicated(dplyr::select(pairs, plotID, cutoff_cm))) {
    stop("Duplicate paired plots.")
  }

  # Check that every class uses the same census interval and per-decade units.
  for (component in components) {
    expected <- 10 * (
      pairs[[paste0(component, "_last")]] - pairs[[paste0(component, "_first")]]
    ) / pairs$years

    if (any(abs(expected - pairs[[paste0(component, "_change_decade")]]) >
            tolerance * pmax(1, abs(expected)))) {
      stop("Inconsistent annualisation.")
    }
  }

  pairs |>
    dplyr::mutate(
      increasing = (N_last - N_first) > tolerance * pmax(1, N_first),
      dlogqmd_decade = log1p(qmd_change_pct_decade / 100)
    )
}

# Estimate the contribution of each size class within a source -----------------
density_source_estimate <- function(df, mode = "all_plots", min_plots = 10L) {
  n_available <- nrow(df)

  if (mode == "increasing_plots") {
    df <- df |>
      dplyr::filter(increasing)
  }

  # This sensitivity keeps the same native-increasing plots at both cutoffs.
  # Otherwise, selecting increases separately can change the comparison cohort.
  if (mode == "native_increasing_matched") {
    df <- df |>
      dplyr::filter(native_increasing %in% TRUE)
  }

  n_plots <- nrow(df)

  if (n_plots == 0L) {
    return(tibble::tibble())
  }

  changes <- df |>
    dplyr::select(
      small = n_small_change_decade,
      medium = n_medium_change_decade,
      large = n_large_change_decade,
      very_large = n_verylarge_change_decade,
      total = N_change_decade
    ) |>
    as.matrix()

  qmd_supported <- NA

  if (mode == "qmd_adjusted") {
    x <- df$dlogqmd_decade

    if (any(!is.finite(x)) || n_plots < 3L || stats::sd(x) < 1e-10) {
      return(tibble::tibble())
    }

    # Use the same design matrix for every response so the fitted intercepts
    # remain additive. The intercept is mean change at zero log-QMD change.
    design <- cbind(1, x)
    fit <- stats::lm.fit(design, changes)

    if (fit$rank < 2L) {
      return(tibble::tibble())
    }

    estimates <- fit$coefficients[1, ]

    # Record extrapolated intercepts but exclude them from across-source tests.
    qmd_supported <- min(x) <= 0 && max(x) >= 0
  } else {

    # Means preserve the additive accounting; medians of classes need not sum
    # to the median total. Each plot receives equal weight within its source.
    estimates <- colMeans(changes)
  }

  stopifnot(
    abs(sum(estimates[c("small", "medium", "large")]) - estimates["total"]) <
      1e-8 * max(1, abs(estimates["total"]))
  )

  tibble::tibble(
    mode = mode,
    n_available = n_available,
    n_plots = n_plots,
    n_increasing = sum(df$increasing),
    qmd_zero_supported = qmd_supported,
    included_in_test = n_plots >= min_plots &&
      (mode != "qmd_adjusted" || qmd_supported),
    small = unname(estimates["small"]),
    medium = unname(estimates["medium"]),
    large = unname(estimates["large"]),
    very_large = unname(estimates["very_large"]),
    total = unname(estimates["total"])
  ) |>
    dplyr::mutate(

      # Ratios of means avoid unstable ratios for plots with almost zero gain.
      # A class may exceed 100% if another class loses stems over the interval.
      small_share = dplyr::if_else(total > 1e-8, small / total, NA_real_),
      medium_share = dplyr::if_else(total > 1e-8, medium / total, NA_real_),
      large_share = dplyr::if_else(total > 1e-8, large / total, NA_real_),
      small_minus_large = small - large,
      small_minus_rest = small - (medium + large),
      below40_minus_above40 = small + medium - large,
      below60_minus_above60 = small + medium + large - 2 * very_large
    )
}

# Test paired small-versus-large contrasts across sources -----------------------
density_contrast_test <- function(contrasts, nboot = 9999L) {
  contrasts <- contrasts[is.finite(contrasts)]
  n_sources <- length(contrasts)
  nonzero <- abs(contrasts) > 1e-8
  n_nonzero <- sum(nonzero)

  # Each source contributes one mean paired contrast. The exact two-sided sign
  # test asks whether positive and negative source contrasts are equally likely.
  # Ties are omitted; trees and site quadrats are not independent replicates.
  sign_test <- NULL
  mean_test <- NULL
  median_boot <- NA_real_

  if (n_nonzero > 0L) {
    sign_test <- stats::binom.test(
      x = sum(contrasts[nonzero] > 0),
      n = n_nonzero,
      p = 0.5
    )
  }

  # The mean contrast and its t interval are a secondary magnitude summary.
  # Resample sources, rather than plots, for the median's bootstrap interval.
  if (n_sources >= 2L && stats::sd(contrasts) > 1e-10) {
    mean_test <- stats::t.test(contrasts, mu = 0)
  }

  if (n_sources >= 2L) {
    median_boot <- purrr::map_dbl(seq_len(nboot), function(i) {
      stats::median(sample(contrasts, n_sources, replace = TRUE))
    })
  }

  tibble::tibble(
    n_sources = n_sources,
    n_small_dominant = sum(contrasts > 1e-8),
    n_large_dominant = sum(contrasts < -1e-8),
    n_ties = n_sources - n_nonzero,
    mean_contrast = if (n_sources > 0L) mean(contrasts) else NA_real_,
    median_contrast = if (n_sources > 0L) stats::median(contrasts) else NA_real_,
    median_q025 = unname(stats::quantile(median_boot, 0.025, na.rm = TRUE)),
    median_q975 = unname(stats::quantile(median_boot, 0.975, na.rm = TRUE)),
    sign_p_two_sided = if (is.null(sign_test)) NA_real_ else sign_test$p.value,
    preference_probability_q025 = if (is.null(sign_test)) NA_real_ else sign_test$conf.int[1],
    preference_probability_q975 = if (is.null(sign_test)) NA_real_ else sign_test$conf.int[2],
    mean_p_two_sided = if (is.null(mean_test)) NA_real_ else mean_test$p.value,
    mean_q025 = if (is.null(mean_test)) NA_real_ else mean_test$conf.int[1],
    mean_q975 = if (is.null(mean_test)) NA_real_ else mean_test$conf.int[2]
  )
}

# Run the primary test and sensitivity analyses --------------------------------
run_density_contribution_analysis <- function(
  pairs,
  out_dir,
  min_plots = 10L,
  nboot = 9999L,
  seed = 731L
) {

  # A >=20 cm population has no small stems by definition. Compare the full
  # recorded population with >=12.7 cm, which excludes FIA's microplot trees.
  df <- density_contribution_data(pairs) |>
    dplyr::filter(cutoff_cm %in% c(0, 12.7))

  ## Match plot histories between the two diameter cutoffs ---------------------
  common <- df |>
    dplyr::count(dataset, plotID, year0, year1) |>
    dplyr::filter(n == 2L) |>
    dplyr::select(-n)

  matched <- df |>
    dplyr::inner_join(common, by = c("dataset", "plotID", "year0", "year1"))

  native <- matched |>
    dplyr::filter(cutoff_cm == 0) |>
    dplyr::select(plotID, year0, year1, native_increasing = increasing)

  matched <- matched |>
    dplyr::left_join(native, by = c("plotID", "year0", "year1"))

  ## Estimate source contributions ---------------------------------------------
  samples <- list(available = df, matched_endpoints = matched)

  sources <- purrr::imap(samples, function(sample_data, sample_name) {
    modes <- c("all_plots", "increasing_plots", "qmd_adjusted")

    if (sample_name == "matched_endpoints") {
      modes <- c(modes, "native_increasing_matched")
    }

    purrr::map(modes, function(mode_name) {
      sample_data |>
        dplyr::group_by(dataset, cutoff_cm) |>
        dplyr::group_modify(function(df, key) {
          density_source_estimate(df, mode_name, min_plots)
        }) |>
        dplyr::ungroup() |>
        dplyr::mutate(sample = sample_name, .before = 1)
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()

  ## Test whether the small-stem contribution is larger ------------------------
  comparisons <- c(
    "small_minus_large", "small_minus_rest",
    "below40_minus_above40", "below60_minus_above60"
  )

  # Isolate bootstrap randomness from other analyses sourcing this workflow.
  tests <- withr::with_seed(seed, {
    sources |>
      dplyr::filter(included_in_test) |>
      dplyr::group_by(sample, cutoff_cm, mode) |>
      dplyr::group_modify(function(df, key) {
        purrr::map(comparisons, function(comparison) {
          density_contrast_test(df[[comparison]], nboot) |>
            dplyr::mutate(
              comparison = comparison,
              mean_total = mean(df$total),
              mean_small = mean(df$small),
              mean_medium = mean(df$medium),
              mean_large = mean(df$large),
              small_share_of_mean_gain = if (mean(df$total) > 1e-8) {
                mean(df$small) / mean(df$total)
              } else {
                NA_real_
              }
            )
        }) |>
          dplyr::bind_rows()
      }) |>
      dplyr::ungroup()
  })

  tests <- tests |>
    dplyr::mutate(
      primary = sample == "available" & cutoff_cm == 0 &
        mode == "increasing_plots" & comparison == "small_minus_large"
    )

  # All additional cutoff/cohort/model tests form one exploratory family.
  # Keep the primary exact p-value separate from their BH-adjusted p-values.
  tests$secondary_sign_p_BH <- NA_real_
  tests$secondary_sign_p_BH[!tests$primary] <- stats::p.adjust(
    tests$sign_p_two_sided[!tests$primary],
    method = "BH"
  )

  ## Save the estimates, tests and reproducibility settings ---------------------
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(sources, file.path(out_dir, "density_contributions_by_source.csv"))
  readr::write_csv(tests, file.path(out_dir, "density_contribution_tests.csv"))

  df |>
    dplyr::select(
      dataset, biome_major, plotID, year0, year1, cutoff_cm, increasing,
      N_change_decade, n_small_change_decade, n_medium_change_decade,
      n_large_change_decade, n_verylarge_change_decade
    ) |>
    readr::write_csv(file.path(out_dir, "density_contributions_by_plot.csv"))

  settings <- list(
    min_plots_per_source = min_plots,
    nboot_source_median = nboot,
    seed = seed,
    primary = "native population; increasing plots; small <20 vs large >=40 cm",
    test = "two-sided exact sign test across source mean paired contrasts",
    run_time = Sys.time()
  )

  readr::write_rds(settings, file.path(out_dir, "density_contribution_settings.rds"))

  ## Plot the additive contributions in density-increasing plots ---------------
  plot_data <- sources |>
    dplyr::filter(
      sample == "available",
      mode == "increasing_plots",
      included_in_test
    ) |>
    dplyr::mutate(
      population = dplyr::if_else(
        cutoff_cm == 0, "All recorded stems", "Stems ≥12.7 cm"
      ),
      label = paste0(dataset, " (n=", n_plots, ")")
    ) |>
    tidyr::pivot_longer(
      c(small_share, medium_share, large_share),
      names_to = "size_class",
      values_to = "share"
    ) |>
    dplyr::mutate(
      size_class = factor(
        size_class,
        levels = c("small_share", "medium_share", "large_share"),
        labels = c("Small: <20 cm", "Medium: 20–40 cm", "Large: ≥40 cm")
      )
    )

  gg_contributions <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(share, label, fill = size_class)
  ) +
    ggplot2::geom_col(orientation = "y", width = 0.7) +
    ggplot2::geom_vline(xintercept = c(0, 1), colour = "grey50", linetype = 2) +
    ggplot2::facet_wrap(~population, scales = "free_y", ncol = 2) +
    ggplot2::scale_x_continuous(labels = scales::label_percent()) +
    ggplot2::scale_fill_manual(values = c("#2b7b69", "#c29c54", "#527ba0")) +
    ggplot2::labs(
      x = "Contribution to mean net density gain",
      y = NULL,
      fill = NULL,
      title = "Which diameter classes account for observed density increases?",
      subtitle = paste(
        "Sources with ≥10 increasing plot series;",
        "net contributions can be negative or exceed 100%"
      ),
      caption = paste0(
        "Shares = mean class change / mean total change. ",
        "Each panel selects increases in its own diameter population.\n",
        "Small means 12.7–20 cm in the right panel. ",
        "Formal replication is by source, not individual tree or quadrat."
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  ggplot2::ggsave(
    filename = file.path(out_dir, "density_gain_contributions.png"),
    plot = gg_contributions,
    width = 12,
    height = 7,
    dpi = 160
  )

  print(dplyr::filter(tests, primary), width = Inf)
  invisible(list(sources = sources, tests = tests))
}
