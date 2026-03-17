identify_ingrowth_plots <- function(df){
  
  # identify ingrowth-affected plots based on whether plot of at least one difference 
  # between inventories looks like disturbance
  plots_disturbed <- df |> 
    group_by(plotID) |> 
    nest() |> 
    mutate(data = purrr::map(data, ~identify_ingrowth_byinventory(.))) |> 
    mutate(ningrowth = purrr::map_int(data, ~get_sum_ingrowth(.))) |> 
    unnest(data)
}

get_sum_ingrowth <- function(df){
  sum(df$ingrowth, na.rm = TRUE)
}

# counts intervals between inventories that look like disturbance was at play
identify_ingrowth_byinventory <- function(df){
  df |> 
    mutate(
      dlogQMD = logQMD - lag(logQMD),
      dlogDensity = logDensity - lag(logDensity)
    ) |> 
    mutate(
      ingrowth = ifelse(dlogQMD < 0 & dlogDensity > 0, TRUE, FALSE)
      )
    
}
