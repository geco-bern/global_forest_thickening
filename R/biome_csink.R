# Functions for upscaling the conditional 90th density quantile by biome.
# Source fit_model.R first; prediction coefficients are returned to their
# original units before calculating changes in density and biomass.

# Extract the fitted models and matching data ----------------------------------
prepare_biome_sink_inputs <- function(analysis) {
  if (anyDuplicated(analysis$biome_major)) {
    stop("Duplicate biome analyses.")
  }

  inputs <- purrr::map(seq_len(nrow(analysis)), function(i) {
    fit <- analysis$out[[i]]$fit_lqmm[[1]]

    if (!inherits(fit, "lqmm") || length(fit$tau) != 1L || fit$tau != 0.9) {
      stop(
        "A fitted 0.9 LQMM from 05_stl_year.R is required for ",
        analysis$biome_major[i]
      )
    }

    # This also rejects older saved fits without their training recipe.
    beta <- coef_original_units(fit)

    if (!all(is.finite(beta))) {
      stop("Non-finite saved LQMM coefficients.")
    }

    # Match fit_model() eligibility, retaining the sample used in the saved fit.
    df <- analysis$data[[i]] |>
      dplyr::ungroup() |>
      dplyr::select(plotID, year, QMD, density, biomass) |>
      dplyr::filter(
        is.finite(density),
        density > 0,
        is.finite(QMD),
        QMD > 0,
        is.finite(year),
        !is.na(plotID)
      )

    list(data = df, fit = fit)
  })

  purrr::set_names(inputs, analysis$biome_major)
}

# Calculate forest-covered area within the mapped forest biomes -----------------
prepare_biome_sink_grid <- function(grid) {

  # Use the same forest-biome delineation as 04_compile_global_drivers.R.
  codes <- c(1L, 2L, 4L, 5L, 6L, 12L)
  labels <- c(
    "Tropical & Subtropical Broadleaf Forests",
    "Tropical & Subtropical Broadleaf Forests",
    "Temperate Forests", "Temperate Forests",
    "Boreal Forests/Taiga", "Mediterranean Forests"
  )

  grid <- grid |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(biome_number), biome_number %in% codes) |>
    dplyr::select(lon, lat, biome_number, biome_major, fcf, area_ha)

  if (nrow(grid) == 0L) {
    stop("No mapped forest-biome cells.")
  }

  if (anyDuplicated(dplyr::select(grid, lon, lat))) {
    stop("Duplicate global grid cells.")
  }

  if (anyNA(grid$biome_major) ||
      any(grid$biome_major != labels[match(grid$biome_number, codes)])) {
    stop("Biome names do not match the forest-biome map codes.")
  }

  if (any(!is.finite(grid$area_ha) | grid$area_ha <= 0) ||
      any(!is.na(grid$fcf) & (!is.finite(grid$fcf) | grid$fcf < 0 | grid$fcf > 1))) {
    stop("Grid areas must be positive and forest cover must be a fraction in [0, 1].")
  }

  grid <- grid |>
    dplyr::mutate(forest_area_ha = area_ha * fcf)

  # Exclude unknown forest cover from the upscaled area, but retain its extent
  # in the coverage table so missing cells cannot be mistaken for non-forest.
  area <- grid |>
    dplyr::group_by(biome_major) |>
    dplyr::summarise(
      n_grid = dplyr::n(),
      n_missing_cover = sum(is.na(fcf)),
      grid_area_missing_cover_ha = sum(area_ha[is.na(fcf)]),
      forest_area_ha = sum(forest_area_ha, na.rm = TRUE),
      .groups = "drop"
    )

  if (any(area$forest_area_ha <= 0)) {
    stop("A biome has no known forest-covered area.")
  }

  list(grid = grid, area = area)
}

# Fit the biomass conversion within one biome ----------------------------------
fit_biome_biomass <- function(data) {

  # Use the through-origin relationship from 08. Missing or negative biomass
  # cannot contribute to the fit; repeated plots supply the random intercepts.
  df <- data |>
    dplyr::filter(
      is.finite(biomass),
      biomass >= 0,
      is.finite(density),
      density > 0,
      is.finite(QMD),
      QMD > 0,
      !is.na(plotID)
    )

  n_plots <- dplyr::n_distinct(df$plotID)

  if (nrow(df) < 3L || n_plots < 2L || nrow(df) <= n_plots) {
    stop("Insufficient biomass observations for a biome-specific mixed model.")
  }

  # Rescale the design column for the optimiser, then convert the fitted slope
  # back to kg dry matter per hectare per unit of N * QMD^2.
  df <- df |>
    dplyr::mutate(NQMD2_scaled = density * QMD^2 / 1e6)

  model <- lme4::lmer(
    biomass ~ 0 + NQMD2_scaled + (1 | plotID),
    data = df
  )

  biomass_coef <- unname(lme4::fixef(model)["NQMD2_scaled"]) / 1e6

  if (!is.finite(biomass_coef) || biomass_coef <= 0) {
    stop("Nonpositive or non-finite biomass conversion.")
  }

  messages <- model@optinfo$conv$lme4$messages

  if (length(messages) > 0L) {
    warning(paste(messages, collapse = "; "))
  }

  biomass_coef
}

# Predict biomass change associated with the quantile shift ---------------------
predict_biome_quantile_sink <- function(
  beta,
  biomass_coef,
  qmd,
  forest_area_ha,
  year0 = 2000,
  year1 = 2001,
  carbon_fraction = 0.5
) {
  beta <- beta[c("(Intercept)", "logQMD", "year")]

  if (any(!is.finite(beta)) || !is.finite(biomass_coef) || biomass_coef <= 0 ||
      !length(qmd) || any(!is.finite(qmd) | qmd <= 0) ||
      length(forest_area_ha) != 1L || !is.finite(forest_area_ha) || forest_area_ha < 0 ||
      !is.finite(year0) || !is.finite(year1) || year1 <= year0 ||
      !is.finite(carbon_fraction) || carbon_fraction <= 0 || carbon_fraction > 1) {
    stop("Invalid biome sink prediction inputs.")
  }

  # Predict the fixed-effect quantile, with plot random intercept set to zero.
  # A lognormal bias correction would change the quantile estimand.
  logN0 <- beta["(Intercept)"] + beta["logQMD"] * log(qmd) + beta["year"] * year0
  dn <- exp(logN0) * expm1(beta["year"] * (year1 - year0))

  # Hold the same QMD values fixed across years. Average the nonlinear biomass
  # changes over the empirical QMD distribution rather than squaring mean QMD.
  db <- mean(biomass_coef * qmd^2 * dn) / (year1 - year0)

  if (!is.finite(db)) {
    stop("Non-finite quantile-based biomass change.")
  }

  tibble::tibble(
    beta_intercept = unname(beta["(Intercept)"]),
    beta_logQMD = unname(beta["logQMD"]),
    beta_year = unname(beta["year"]),
    year_change_percent = 100 * expm1(unname(beta["year"])),
    biomass_coef = biomass_coef,
    db_kgdm_ha_yr = db,
    dc_mgc_ha_yr = db * carbon_fraction * 1e-3,
    sink_pgc_yr = db * carbon_fraction * forest_area_ha * 1e-12
  )
}

# Refit a biome LQMM, allowing one optimisation retry ----------------------------
fit_biome_quantile <- function(data) {

  # Keep the 05 settings. Retry a reported convergence failure with a larger
  # iteration budget and retain warnings for the bootstrap diagnostics.
  for (max_iter in c(500L, 2000L)) {
    messages <- character()

    fit <- withCallingHandlers(
      withr::with_seed(123L, {
        fit_model(
          data,
          lqmm = TRUE,
          control = list(LP_max_iter = max_iter)
        )
      }),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    nonconverged <- any(grepl(
      "did not converge|failed to converge",
      messages,
      ignore.case = TRUE
    ))

    if (!nonconverged) {
      for (message in unique(messages)) {
        warning(message, call. = FALSE)
      }

      if (max_iter > 500L) {
        warning("LQMM converged after retry with LP_max_iter = 2000.", call. = FALSE)
      }

      return(fit)
    }
  }

  stop("LQMM did not converge after retry: ", paste(unique(messages), collapse = "; "))
}

# Combine the density and biomass relationships --------------------------------
estimate_biome_sink <- function(
  data,
  forest_area_ha,
  fit = NULL,
  year0 = 2000,
  year1 = 2001,
  carbon_fraction = 0.5
) {
  if (is.null(fit)) {

    # fit_model() sets a seed internally; the wrapper isolates it from resampling.
    fit <- fit_biome_quantile(data)
  }

  if (is.null(fit)) {
    stop("LQMM fit unavailable for this sample.")
  }

  predict_biome_quantile_sink(
    beta = coef_original_units(fit),
    biomass_coef = fit_biome_biomass(data),
    qmd = data$QMD,
    forest_area_ha = forest_area_ha,
    year0 = year0,
    year1 = year1,
    carbon_fraction = carbon_fraction
  )
}

# Run one bootstrap replicate --------------------------------------------------
biome_sink_bootstrap_job <- function(job, inputs, area, year0, year1, carbon_fraction) {
  warnings <- character()

  # Retain failed fits as NA, with the error and warnings attached to that draw.
  result <- tryCatch(
    withCallingHandlers(
      {
        data <- inputs[[job$biome_major]]$data |>
          dplyr::slice(job$indices)

        estimate_biome_sink(
          data = data,
          forest_area_ha = area$forest_area_ha[match(job$biome_major, area$biome_major)],
          year0 = year0,
          year1 = year1,
          carbon_fraction = carbon_fraction
        )
      },
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      tibble::tibble(sink_pgc_yr = NA_real_, error = conditionMessage(e))
    }
  )

  if (!"error" %in% names(result)) {
    result <- result |>
      dplyr::mutate(error = NA_character_)
  }

  result |>
    dplyr::mutate(
      biome_major = job$biome_major,
      replicate = job$replicate,
      warnings = paste(unique(warnings), collapse = "; "),
      .before = 1
    )
}

# Bootstrap independently within each biome ------------------------------------
bootstrap_biome_sink <- function(
  inputs,
  area,
  nboot = 30L,
  workers = 1L,
  seed = 1982L,
  year0 = 2000,
  year1 = 2001,
  carbon_fraction = 0.5
) {
  stopifnot(
    length(nboot) == 1L,
    is.finite(nboot),
    nboot >= 1,
    nboot == as.integer(nboot),
    length(workers) == 1L,
    is.finite(workers),
    workers >= 1,
    workers == as.integer(workers)
  )

  if (!setequal(names(inputs), area$biome_major)) {
    stop("Model and mapped biomes differ.")
  }

  # Match the row-resampling scheme in 05. Resample all census variables jointly
  # to preserve covariance between the STL and biomass fits. Draw indices before
  # fitting so a seed reset inside fit_model() cannot duplicate bootstrap draws.
  jobs <- withr::with_seed(seed, {
    purrr::map(sort(names(inputs)), function(biome) {
      n <- nrow(inputs[[biome]]$data)

      purrr::map(seq_len(nboot), function(i) {
        list(
          biome_major = biome,
          replicate = i,
          indices = sample.int(n, n, replace = TRUE)
        )
      })
    }) |>
      unlist(recursive = FALSE)
  })

  # Workers need the compact training data, not the saved point-estimate fits.
  inputs <- purrr::map(inputs, function(input) list(data = input$data))
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (workers == 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = min(workers, length(jobs)))
  }

  result <- furrr::future_map_dfr(
    jobs,
    biome_sink_bootstrap_job,
    inputs = inputs,
    area = area,
    year0 = year0,
    year1 = year1,
    carbon_fraction = carbon_fraction,
    .options = furrr::furrr_options(seed = seed),
    .progress = interactive()
  )

  if (anyNA(result$sink_pgc_yr)) {
    warning("Failed biome bootstrap fits retained as NA; see error column.")
  }

  result
}

# Summarise the biome and global bootstrap draws --------------------------------
summarise_biome_sink <- function(draws) {
  draws |>
    dplyr::group_by(biome_major) |>
    dplyr::summarise(
      nboot = dplyr::n(),
      n_success = sum(is.finite(sink_pgc_yr)),
      mean = if (n_success > 0) mean(sink_pgc_yr, na.rm = TRUE) else NA_real_,
      median = if (n_success > 0) median(sink_pgc_yr, na.rm = TRUE) else NA_real_,
      sd = stats::sd(sink_pgc_yr, na.rm = TRUE),
      q025 = unname(stats::quantile(sink_pgc_yr, 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(sink_pgc_yr, 0.975, na.rm = TRUE)),
      .groups = "drop"
    )
}

sum_biome_sink_draws <- function(draws, biomes) {

  # A failed or missing biome must never become a zero contribution to the sum.
  draws |>
    dplyr::group_by(replicate) |>
    dplyr::summarise(
      sink_pgc_yr = if (
        dplyr::n() == length(biomes) && setequal(biome_major, biomes) &&
          all(is.finite(sink_pgc_yr))
      ) {
        sum(sink_pgc_yr)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(biome_major = "Global")
}
