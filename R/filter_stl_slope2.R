filter_stl_slope2 <- function(df, dataset_name){

  df_byplot <- df |>
    group_by(plotID) |>
    nest() |>
    mutate(nobs = purrr::map_int(data, ~nrow(.))) |>
    filter(nobs >= 3) |>  # minimum satisfied by design of data selection, but here for safety

    # linear model
    mutate(
      linmod = purrr::map(data, ~lm(logDensity ~ logQMD, data = .))
    ) |>
    
    # segmented model (1 breakpoint 2 segments)
    mutate(
      segmod = map(linmod, ~ {
        tryCatch(
          segmented(.x, seg.Z = ~logQMD, npsi = 1),
          error = function(e) NULL
        )
      })
    ) |> 
    
    mutate(
      slope_lin = map_dbl(linmod, ~ coef(.x)[2]),
      slopes_seg = map(segmod, ~ if (is.null(.x)) NA else slope(.x)$logQMD_sc[, "Est."])
    )
  
  # Retain only inner 90% of the data by slope from simple linear regression
  range_accepted <- quantile(df_byplot$slope_lin, c(0.05, 0.95), na.rm = TRUE)

  gg <- df_byplot |> 
    ggplot(aes(slope_lin, after_stat(count))) +
    geom_histogram(color = "black", fill = "grey70") +
    xlim(range_accepted + c(-10, 10)) +
    annotate(
      "rect",
      xmin = range_accepted[1], xmax = range_accepted[2],
      ymin = 0, ymax = Inf,
      fill = "red",
      alpha = 0.3  # transparency
    ) +
    labs(
      x = "slope", 
      y = "Count", 
      title = dataset_name
    )

  df <- df_byplot |> 
    mutate(
      bad_lin = slope_lin < range_accepted[1] | slope_lin > range_accepted[2],
      bad_seg = map_lgl(slopes_seg, ~ any(.x < range_accepted[1] | .x > range_accepted[2], na.rm = TRUE))
    ) |> 
    mutate(
      badslope = bad_lin | bad_seg
    ) %>%
    dplyr::select(plotID, data, badslope) |> 
    unnest(data) |> 
    ungroup()

  return(
    list(
      df = df,
      gg = gg
    ))
}