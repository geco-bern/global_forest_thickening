plot_stl_longplots <- function(df, biome_name = "") {

  set.seed(1982)

  df <- df |>
    filter(
      biome_major == biome_name,
      ndisturbed == 0,
      !badqmdbin,
      !badslope,

      # Pristine/primary/old-growth/protected
      # management_cat == 2
    )

    # for demo keep only best
  tmp <- df |>
    group_by(plotID, dataset) |>
    summarise(year_min = min(year), year_max = max(year)) |>
    mutate(len = year_max - year_min) |>
    arrange(desc(len)) |>
    ungroup()

  thresh_len <- floor(quantile(tmp$len, probs = 0.95, na.rm = TRUE))

  df_longplots <- tmp |>
    slice(1:30) |> 
    # filter(len >= thresh_len) |>
    # slice_sample(n = 10) |> # used 3 before
    select(plotID, dataset) |>
    distinct()

  df <- df |>
    mutate(highlight = ifelse(plotID %in% df_longplots$plotID, TRUE, FALSE))

  gg <- ggplot() +
    geom_point(
      aes(x = logQMD, y = logDensity),
      data = df,
      color = "grey80",
      # alpha = 0.1,
      size = 1,
      shape = 16,
      inherit.aes = FALSE
    ) +
    geom_line(
      aes(x = logQMD, y = logDensity, group = plotID), #, color = dataset
      data = df |>
        filter(highlight),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    labs(
      title = biome_name,
      x = expression(ln(QMD)),
      y = expression(ln(italic(N)))
    ) +
    scale_color_viridis_d(name = "Dataset") +
    theme_classic() +
    scale_x_continuous(limits = c(2.2, 4.7), breaks = seq(3, 4, 1)) +
    scale_y_continuous(limits = c(2.9, 9.3), breaks = seq(4, 8, 2))

  return(gg)
}
