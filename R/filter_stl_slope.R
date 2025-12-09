filter_stl_slope <- function(df, fit_lqmm){
  # Filters data based on whether a plot's regression slope of logN ~ logQMD is
  # within a certain margin of the regression slope obtained from the quantile
  # regression mixed effects model.
  df_byplot <- df |>
    group_by(plotID, dataset) |>
    nest() |>
    mutate(nobs = purrr::map_int(data, ~nrow(.))) |>
    filter(nobs >= 3) |>  # minimum satisfied by design of data selection, but here for safety
    mutate(linmod = purrr::map(data, ~lm(logDensity ~ logQMD_sc, data = .))) |>
    mutate(slope = purrr::map_dbl(linmod, ~coef(.)[2])) |>
    select(-linmod)

  # remove outliers (top and bottom 1% of data)
  df_byplot <- df_byplot |>
    filter(slope > quantile(df_byplot$slope, 0.01) & slope < quantile(df_byplot$slope, 0.99))

  slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]
  sd_slopes <- sd(df_byplot$slope)

  # subset data based on slope: must be within a certain margin of slope obtained from quantile regression
  # margin defined as half a standard deviation across all slopes
  df_byplot <- df_byplot |>
    mutate(filter_slope = slope < slope_lqmm + 0.5 * sd_slopes & slope > slope_lqmm - 0.5 * sd_slopes)

  return(df_byplot)
}
