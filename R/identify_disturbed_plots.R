# idenify disturbed plots as those where both logQMD and logDensity decline from one census to another
# the returned data frame contains 'ndisturbed' which counts, by plot, the between census intervals identified as disturbance-affected
identify_disturbed_plots <- function(df){
  
  # identify disturbed plots based on whether plot of at least one difference 
  # between inventories looks like disturbance
  plots_disturbed <- df |> 
    group_by(plotID) |> 
    nest() |> 
    mutate(data = purrr::map(data, ~identify_disturbed_byinventory(.))) |> 
    mutate(ndisturbed = purrr::map_int(data, ~get_sum_disturbed(.))) |> 
    unnest(data)
}

get_sum_disturbed <- function(df){
  sum(df$disturbed, na.rm = TRUE)
}

# counts intervals between inventories that look like disturbance was at play
identify_disturbed_byinventory <- function(df){
  df |>
    mutate(
      logQMD = log(QMD),
      logDensity = log(density)
    ) |>
    mutate(
      dlogQMD = logQMD - lag(logQMD),
      dlogDensity = logDensity - lag(logDensity)
    ) |>
    mutate(
      disturbed = ifelse(dlogQMD < 0 & dlogDensity < 0, TRUE, FALSE)
    ) |>
    dplyr::select(-logQMD, -logDensity)    
    
}
