## Define workflow -------------
analyse_group <- function(df, group_name, n_boot = 500) {
  # Analyse disturbance trend (before applying disturbed-filter)
  breaks <- get_breaks(df$year)
  gg_fdisturbed <- plot_disturbed(df, group_name, breaks)

  # Define filter stages --------------
  list_df_filtered <- list(
    # all data
    "orig" = df,

    # no disturbance-affecte plots
    "no_disturbance" = df |>
      filter(ndisturbed == 0),

    # no data from years with shifted QMD distribution
    "no_badqmd" = df |>
      filter(ndisturbed == 0, !badqmdbin),

    # no data from plots with outlying self-thinning slope
    "no_badslope" = df |>
      filter(ndisturbed == 0, !badqmdbin, !badslope),

    # no data from plots without management history and management <30 years prior to first census
    "no_mgmt_30" = df |>
      filter(
        ndisturbed == 0,
        !badqmdbin,
        !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 30) |
          management_cat == 2
      ),

    # no data from plots without management history and management <100 years prior to first census
    "no_mgmt_100" = df |>
      filter(
        ndisturbed == 0,
        !badqmdbin,
        !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 100) |
          management_cat == 2
      ),

    # only old-growth
    "oldgrowth" = df |>
      filter(
        ndisturbed == 0,
        !badqmdbin,
        !badslope,

        # Pristine/primary/old-growth/protected
        management_cat == 2 |
          (management_cat == 1 & management_since_census1_yrs >= 300)
      )
  )

  # Analyse sensitivity of fit ---------------
  # (year coefficient) subject filter levels
  df_filtereffects <- imap_dfr(list_df_filtered, \(df, step_name) {
    df |>
      ungroup() |>
      nest() |>
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

  # Create plot --------------  
  gg_coef_filters <- plot_filtereffects(
    df_filtereffects,
    group_name = group_name
  )

  # Create classic LQMM plot -------------
  # for filter level orig
  fit_lqmm <- df_filtereffects |>
    filter(step == "orig") |>
    pull(fit_lqmm)

  # Bootstrap LQMM fit for the selected filter level ('orig') only.
  fit_bootstrap_sample <- function(split) {
    # Drop failed fits; retain convergence warnings for diagnosis.
    tryCatch(
      {
        model <- fit_model(rsample::analysis(split), lqmm = TRUE)
        beta <- coef_original_units(model)

        # The year estimate is per calendar year.
        tibble::tibble(term = names(beta), estimate = unname(beta))
      },
      error = function(e) NULL
    )
  }

  boot_data <- rsample::bootstraps(
    list_df_filtered$orig %>%
      group_by(plotID),
    times = n_boot,
    apparent = FALSE
  )

  boot_results <- boot_data %>%
    mutate(coefs = map(splits, fit_bootstrap_sample)) %>%
    filter(!map_lgl(coefs, is.null)) %>%
    unnest(coefs) |>
    dplyr::select(-splits)

  df_boot <- boot_results |>
    filter(term == "year") |>
    mutate(group = group_name)

  gg_boot <- plot_boot(df_boot)

  # get year effect within QMD bin before applying any filters
  gg_lqmm <- plot_lqmm_bygroup(
    list_df_filtered$orig,
    fit_lqmm[[1]]
  ) +
    labs(title = group_name)

  # LQMM by QMD bin ----------------
  # for filter orig
  df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(list_df_filtered$orig)
  gg_lqmm_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df)

  gg_lqmm_qmdbin_boot <- combine_plot_qmdbin_boot(
    gg_lqmm_byqmdbin,
    gg_boot,
    group_name = group_name
  )

  return(
    list(
      df_filtereffects = df_filtereffects,
      fit_lqmm = fit_lqmm,
      boot_results = boot_results,
      df_boot = df_boot,
      gg_fdisturbed = gg_fdisturbed,
      gg_lqmm = gg_lqmm,
      gg_coef_filters = gg_coef_filters,
      gg_lqmm_qmdbin_boot = gg_lqmm_qmdbin_boot
    )
  )
}
