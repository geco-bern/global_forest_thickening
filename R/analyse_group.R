analyse_group <- function(
  df,
  group_name,
  do_boot = FALSE,
  n_boot = 500,
  boot_workers = parallelly::availableCores(omit = 1L),
  boot_seed = 123L
) {
  # Analyse disturbance ---------
  # trend (before applying disturbed-filter)
  breaks <- get_breaks(df$year)
  gg_fdisturbed <- plot_disturbed(df, group_name, breaks)

  # Plot-level slopes and pooled 90th-quantile slope, using all group data.
  slope_analysis <- analyse_plot_slopes(df, group_name)

  # Define filter stages --------------
  list_df_filtered <- list(
    # all data
    "orig" = df,

    # # no disturbance-affecte plots
    # "no_disturbance" = df |>
    #   filter(ndisturbed == 0),

    # # no data from years with shifted QMD distribution
    # "no_badqmd" = df |>
    #   filter(
    #     # ndisturbed == 0,
    #     !badqmdbin
    #   ),

    # # no data from plots with outlying self-thinning slope
    # "no_badslope" = df |>
    #   filter(
    #     # ndisturbed == 0,
    #     !badqmdbin,
    #     !badslope
    #   ),

    # no data from plots without management history and management <30 years prior to first census
    "no_mgmt_30" = df |>
      filter(
        # ndisturbed == 0,
        # !badqmdbin,
        # !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 30) |
          management_cat == 2
      ),

    # no data from plots without management history and management <100 years prior to first census
    "no_mgmt_100" = df |>
      filter(
        # ndisturbed == 0,
        # !badqmdbin,
        # !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 100) |
          management_cat == 2
      ),

    # only old-growth
    "oldgrowth" = df |>
      filter(
        # ndisturbed == 0,
        # !badqmdbin,
        # !badslope,

        # Pristine/primary/old-growth/protected
        management_cat == 2 |
          (management_cat == 1 & management_since_census1_yrs >= 300)
      )
  )

  # Analyse sensitivity to filters ---------------
  # (year coefficient) subject filter levels
  df_filtereffects <- imap_dfr(list_df_filtered, \(df, step_name) {
    tibble::tibble(data = list(dplyr::ungroup(df))) |>
      mutate(
        n = map_int(data, nrow),

        # store both model + coefficients
        res = map(data, \(d) {
          model <- fit_model(d, lqmm = TRUE)
          coefs <- extract_coef(model, lqmm = TRUE)

          list(
            model = model,
            coefs = coefs
          )
        })
      ) |>
      mutate(
        fit_lqmm = purrr::map(res, "model"),
        coef_lqmm = purrr::map(res, "coefs")
      ) |> # creates columns: model, coefs
      tidyr::unnest(coef_lqmm) |> # expand coefficient table
      dplyr::select(-data) |>
      mutate(step = step_name)
  })

  ## Create filter effects plot --------------
  gg_coef_filters <- plot_filtereffects(
    df_filtereffects,
    group_name = group_name
  )

  # Create "classic" LQMM plot -------------
  ## Filter level 'orig' ----------
  fit_lqmm_orig <- df_filtereffects |>
    filter(step == "orig") |>
    pull(fit_lqmm)

  if (!is.null(fit_lqmm_orig[[1]])) {
    gg_lqmm_orig <- plot_lqmm_bygroup(
      list_df_filtered$orig,
      fit_lqmm_orig[[1]]
    ) +
      labs(title = group_name)
  } else {
    gg_lqmm_orig <- NULL
  }

  ## Filter level no-management-100-years -------------
  fit_lqmm_no_mgmt_100 <- df_filtereffects |>
    filter(step == "no_mgmt_100") |>
    pull(fit_lqmm)

  if (!is.null(fit_lqmm_no_mgmt_100[[1]])) {
    gg_lqmm_no_mgmt_100 <- plot_lqmm_bygroup(
      list_df_filtered$no_mgmt_100,
      fit_lqmm_no_mgmt_100[[1]]
    ) +
      labs(title = group_name)
  } else {
    gg_lqmm_no_mgmt_100 <- NULL
  }

  ## Filter level oldgrowth -----------
  fit_lqmm_oldgrowth <- df_filtereffects |>
    filter(step == "oldgrowth") |>
    pull(fit_lqmm)

  if (!is.null(fit_lqmm_oldgrowth[[1]])) {
    gg_lqmm_oldgrowth <- plot_lqmm_bygroup(
      list_df_filtered$oldgrowth,
      fit_lqmm_oldgrowth[[1]]
    ) +
      labs(title = group_name)
  } else {
    gg_lqmm_oldgrowth <- NULL
  }

  # LQMM by QMD bin ----------------
  # for filter orig only
  df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(list_df_filtered$orig)
  gg_lqmm_byqmdbin <- plot_lqmm_byqmdbin(
    df_lqmm_byqmdbin$df,
    name = group_name
  )

  # Bootstrap fit -----------
  # for filter orig only
  if (do_boot && !is.null(fit_lqmm_orig[[1]])) {
    # Bootstrap LQMM fit for the selected filter level ('orig') only.
    boot_results <- bootstrap_lqmm(
      list_df_filtered$orig,
      n_boot = n_boot,
      workers = boot_workers,
      seed = boot_seed
    )

    df_boot <- boot_results |>
      filter(term == "year") |>
      mutate(group = group_name)

    gg_boot <- plot_boot(df_boot)
  } else {
    gg_boot <- NULL
    df_boot <- NULL
    boot_results <- NULL
  }

  # Function output -----------------
  return(
    list(
      df_filtereffects = df_filtereffects,
      df_slopes = slope_analysis$df_slopes,
      fit_qr = slope_analysis$fit_qr,
      slope_qr = slope_analysis$slope_qr,
      gg_slopes = slope_analysis$gg_slopes,
      fit_lqmm = fit_lqmm_orig,
      boot_results = boot_results,
      df_boot = df_boot,
      gg_fdisturbed = gg_fdisturbed,
      gg_lqmm_orig = gg_lqmm_orig,
      gg_lqmm_no_mgmt_100 = gg_lqmm_no_mgmt_100,
      gg_lqmm_oldgrowth = gg_lqmm_oldgrowth,
      gg_lqmm_byqmdbin = gg_lqmm_byqmdbin,
      gg_coef_filters = gg_coef_filters
    )
  )
}

# Diagnostic only: these slopes do not filter records or change the LQMM fits.
analyse_plot_slopes <- function(df, group_name) {
  # Derive unscaled natural logs locally, ignoring any stale log columns.
  valid <- df |>
    dplyr::ungroup() |>
    dplyr::filter(is.finite(QMD), QMD > 0, is.finite(density), density > 0) |>
    dplyr::mutate(logN = log(density), logQMD = log(QMD))

  per_plot <- valid |>
    dplyr::group_by(plotID) |>
    tidyr::nest() |>
    dplyr::mutate(
      nobs = purrr::map_int(data, nrow),
      slope_lin = purrr::map_dbl(data, function(d) {
        # Same minimum of three observations as filter_stl_slope2().
        if (nrow(d) < 3L || dplyr::n_distinct(d$logQMD) < 2L) return(NA_real_)
        fit <- stats::lm(logN ~ logQMD, data = d)
        unname(stats::coef(fit)["logQMD"])
      })
    ) |>
    dplyr::ungroup() |>
    dplyr::select(plotID, nobs, slope_lin)

  # Preserve one row per input plot, including those without a usable slope.
  df_slopes <- df |>
    dplyr::ungroup() |>
    dplyr::distinct(plotID) |>
    dplyr::left_join(per_plot, by = "plotID") |>
    dplyr::mutate(nobs = dplyr::coalesce(nobs, 0L))

  fit_qr <- NULL
  if (nrow(valid) >= 3L && dplyr::n_distinct(valid$logQMD) >= 2L) {
    fit_qr <- tryCatch(
      quantreg::rq(logN ~ logQMD, tau = 0.9, data = valid),
      error = function(e) {
        warning("Pooled quantile regression failed for ", group_name, ": ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  }

  slope_qr <- if (is.null(fit_qr)) NA_real_ else unname(stats::coef(fit_qr)["logQMD"])

  gg_slopes <- df_slopes |>
    dplyr::filter(is.finite(slope_lin)) |>
    ggplot2::ggplot(ggplot2::aes(x = slope_lin)) +
    ggplot2::geom_histogram(bins = 30, colour = "black", fill = "grey70") +
    ggplot2::labs(x = "Slope", y = "Count", title = group_name) +
    ggplot2::theme_classic() +
    ggplot2::xlim(
      quantile(df_slopes$slope_lin, c(0.10, 0.90), na.rm = TRUE) + c(-10, 10)
    )

    if (is.finite(slope_qr)) {
      gg_slopes <- gg_slopes +
        ggplot2::geom_vline(
          xintercept = slope_qr,
          colour = "red",
          linewidth = 0.7
        )
    }

  list(
    df_slopes = df_slopes,
    fit_qr = fit_qr,
    slope_qr = slope_qr,
    gg_slopes = gg_slopes
  )
}

# Separate helpers avoid exporting analyse_group()'s fitted models and plots
# along with the bootstrap inputs. Requires future, furrr and parallelly.
bootstrap_lqmm <- function(df, n_boot, workers = 1L, seed = 123L) {

  stopifnot(
    length(n_boot) == 1L, is.finite(n_boot), n_boot >= 1,
    n_boot == as.integer(n_boot),
    length(workers) == 1L, is.finite(workers), workers >= 1,
    workers == as.integer(workers)
  )

  # Keep the existing resampling scheme, but transfer only model inputs.
  # Seed resampling explicitly so scheduling/worker count cannot change it.
  model_data <- df |>
    dplyr::ungroup() |>
    dplyr::select(plotID, year, QMD, density) |>
    dplyr::group_by(plotID)

  boot_data <- withr::with_seed(seed, rsample::bootstraps(
    model_data,
    times = n_boot,
    apparent = FALSE
  ))

  # Send row indices per replicate and one shared data argument per chunk,
  # rather than serializing the full data inside every resampling split.
  boot_indices <- purrr::map(boot_data$splits, as.integer)
  boot_data <- dplyr::select(boot_data, -splits)

  # # Positron can inject .ark_* debugger calls when sourcing functions. Those
  # # hooks do not exist in plain R workers. If present, read clean definitions
  # # from disk using parse/eval, bypassing the IDE's source() instrumentation.
  # worker_fit <- fit_bootstrap_sample
  # worker_functions <- list(
  #   fit_bootstrap_sample, fit_model, coef_original_units, model_recipe
  # )
  # has_ide_hooks <- vapply(worker_functions, function(fn) {
  #   any(grepl(".ark_", deparse(body(fn)), fixed = TRUE))
  # }, logical(1))
  # if (any(has_ide_hooks)) {
  #   worker_env <- new.env(parent = asNamespace("dplyr"))
  #   for (file in c("fit_model.R", "analyse_group.R")) {
  #     eval(parse(
  #       text = readLines(here::here("R", file), warn = FALSE),
  #       keep.source = FALSE
  #     ), envir = worker_env)
  #   }
  #   worker_fit <- worker_env$fit_bootstrap_sample
  # }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  workers <- min(workers, n_boot)
  if (workers == 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = workers)
  }

  boot_data$coefs <- furrr::future_map(
    boot_indices,
    fit_bootstrap_sample,
    data = model_data,
    .options = furrr::furrr_options(seed = seed, packages = "dplyr")
  )

  boot_data |>
    dplyr::filter(!purrr::map_lgl(coefs, is.null)) |>
    tidyr::unnest(coefs, ptype = list(
      coefs = tibble::tibble(term = character(), estimate = double())
    ))
}

fit_bootstrap_sample <- function(indices, data) {
  # Drop failed fits; retain convergence warnings for diagnosis.
  tryCatch(
    {
      # Each resample learns its own preprocessing recipe in fit_model().
      model <- fit_model(vctrs::vec_slice(data, indices), lqmm = TRUE)
      if (is.null(model)) return(NULL)
      beta <- coef_original_units(model)

      # The year estimate is per calendar year.
      tibble::tibble(term = names(beta), estimate = unname(beta))
    },
    error = function(e) NULL
  )
}
