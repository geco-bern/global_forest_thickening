plot_disturbed <- function(df, biome_name){

  df_disturbed <- df |>
    mutate(year_bin = cut(
      year,
      breaks = breaks,
      labels = breaks[1:length(breaks) - 1] + 2.5,
      include.lowest = TRUE
    )) |>
    group_by(year_bin) |>
    summarise(
      nplots = length(unique(plotID)),
      ndisturbed = sum(disturbed, na.rm = TRUE)
    ) |>
    mutate(fdisturbed = ndisturbed / nplots) |>
    filter(fdisturbed > 0) |>   # unfortunately needed
    mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed))) |>
    mutate(year_dec = as.numeric(as.character(year_bin)))

  if (nrow(df_disturbed) > 2){
    linmod_fdistrubed <- lm(
      fdisturbed_logit ~ year_dec,
      data = df_disturbed
    )

    tab_lm <- broom::tidy(linmod_fdistrubed, conf.int = TRUE) |>
      filter(term == "year_dec") |>
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  }

  gg <- df_disturbed |>
    ggplot(aes(year_dec, fdisturbed_logit)) +
    geom_point() +
    theme_classic() +
    labs(
      x = "Year",
      title = biome_name,
      subtitle = bquote(
        "Slope: " ~ .(tab_lm$estimate) ~ "[" ~
          .(tab_lm$conf.low) * "," ~
          .(tab_lm$conf.high) * "]"
      )
    ) +
    scale_y_continuous(
      name = expression(logit(Fraction ~ disturbed)),
      sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
    )

  if (nrow(df_disturbed) > 2){
    gg <- gg +
      geom_smooth(method = "lm", color = "red")
  }

  return(gg)
}
