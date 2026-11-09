# Year coefficient across filter levels, as a separate plot for one group.
plot_filtereffects <- function(df_filtereffects, group_name) {
  df_filtereffects |>
    dplyr::mutate(
      step = forcats::fct_relevel(
        step,
        "oldgrowth",
        "no_mgmt_100",
        "no_mgmt_30",
        "no_badslope",
        "no_badqmd",
        "no_disturbance",
        "orig",
      )
    ) |>
    ggplot(aes(x = step, y = estimate)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(
      title = group_name,
      x = NULL,
      y = "Coefficient (year)"
    ) +
    # Add observation counts to the right of points
    geom_text(
      aes(label = n),
      y = Inf,
      hjust = 1.1,
      vjust = 1.2,
      size = 3
    ) +
    theme_classic() +
    # ylim(-0.5, 0.5) +
    coord_flip()
}
