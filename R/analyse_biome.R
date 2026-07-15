## Define workflow -------------
analyse_biome <- function(df, biome_name){

  # Analyse disturbance trend (before applying disturbed-filter)
  breaks <- get_breaks(df$year)
  gg_fdisturbed <- plot_disturbed(df, biome_name, breaks)

  # Define filter stages --------------
  list_df_filtered <- list(

    # all data
    "raw" = df,

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
        ndisturbed == 0, !badqmdbin, !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 30) | management_cat == 2
      ),

    # no data from plots without management history and management <100 years prior to first census
    "no_mgmt_100" = df |>
      filter(
        ndisturbed == 0, !badqmdbin, !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 100) | management_cat == 2
      ),

    # only old-growth
    "primary" = df |> 
      filter(
        ndisturbed == 0, !badqmdbin, !badslope,

        # Pristine/primary/old-growth/protected
        management_cat == 2
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
      ) |>     # creates columns: model, coefs
      tidyr::unnest(coef_lqmm) |>         # expand coefficient table
      dplyr::select(-data) |>
      mutate(step = step_name)
  })

  # Create plot --------------
  # of coefficient 'year' for all filter levels
  gg_coef_filters <- df_filtereffects |>
    mutate(
      step = fct_relevel(
        step,
        "primary",
        "no_mgmt_100",
        "no_mgmt_30",
        "no_badslope",
        "no_badqmd",
        "no_disturbance",
        "raw",
      )
    ) |> 
    ggplot(aes(x = step, y = estimate)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    # geom_text(aes(label = paste0("n=", n)), hjust = -0.9, size = 3) +
    # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(
      x = "Filtering step",
      y = "Coefficient (year)"
    ) +
    # Add observation counts to the right of points
    geom_text(
      aes(label = n), 
      y = Inf,
      hjust = 1.1,
      size = 3
    ) + 
    theme_classic() +
    ylim(-0.5, 0.5) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # Create classic LQMM plot -------------
  # for filter level no_badslope
  fit_lqmm <- df_filtereffects |> 
      filter(step == "no_badslope") |> 
      pull(fit_lqmm)

  # Bootstrap LQMM fit for the selected filter level only.
  boot_data <- rsample::bootstraps(
    list_df_filtered$no_badslope %>%
      group_by(plotID),
    times = 500,
    apparent = FALSE
  )

  boot_results <- boot_data %>%
    mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
    filter(!map_lgl(coefs, is.null)) %>%
    unnest(coefs) |>
    dplyr::select(-splits)

  df_boot <- boot_results |>
    filter(term == "year") |>
    mutate(biome = biome_name)
  
  gg_lqmm <- plot_lqmm_bybiome(
    list_df_filtered$no_badslope,
    fit_lqmm[[1]]
  )

  # LQMM by QMD bin ----------------
  # for filter no_badslope
  df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(list_df_filtered$no_badslope)
  gg_lqmm_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df)

  return(
    list(
      gg_fdisturbed = gg_fdisturbed,
      df_filtereffects = df_filtereffects,
      gg_lqmm = gg_lqmm,
      fit_lqmm = fit_lqmm,
      boot_results = boot_results,
      df_boot = df_boot,
      gg_coef_filters = gg_coef_filters,
      gg_lqmm_byqmdbin = gg_lqmm_byqmdbin
    )
  )

}
