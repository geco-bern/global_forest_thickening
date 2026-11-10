plot_logn_logqmd_dataset <- function(df) {

  alpha <- min(1, max(0, 50 / length(unique(df$plotID))))

  # # Rank usable plot records by elapsed years, then net QMD increase.
  # trajectories <- df |>
  #   filter(
  #     !is.na(plotID), is.finite(year),
  #     is.finite(QMD), QMD > 0, is.finite(density), density > 0
  #   ) |>
  #   arrange(plotID, year)
  # if ("census" %in% names(trajectories)) {
  #   trajectories <- trajectories |> arrange(plotID, year, census)
  # }

  # longest_records <- trajectories |>
  #   group_by(plotID) |>
  #   summarise(
  #     record_length = last(year) - first(year),
  #     qmd_increase = last(QMD) - first(QMD),
  #     .groups = "drop"
  #   ) |>
  #   filter(record_length > 0) |>
  #   arrange(desc(record_length), desc(qmd_increase), plotID) |>
  #   slice_head(n = 3)

  # trajectories <- trajectories |>
  #   semi_join(longest_records, by = "plotID")

  ggplot() +
    geom_point(
      data = df,
      aes(x = QMD, y = density),
      alpha = alpha,
      size = 1.5,
      col = "black",
    ) +
    scale_x_log10() +
    scale_y_log10() +
    # # Follow census order, including intervals when QMD decreases.
    # geom_path(
    #   data = trajectories,
    #   aes(x = log(QMD), y = log(density), group = plotID),
    #   alpha = 1,
    #   linewidth = 0.6,
    #   col = "tomato"
    # ) +
    labs(title = as.character(unique(df$dataset))) +
    theme_classic()
}
