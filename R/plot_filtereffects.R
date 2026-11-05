# of coefficient 'year' for all filter levels
plot_filtereffects <- function(df_filtereffects, group_name) {
  df_filtereffects |>
    mutate(
      step = fct_relevel(
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
    # geom_text(aes(label = paste0("n=", n)), hjust = -0.9, size = 3) +
    # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(
      title = group_name,
      x = "Filtering step",
      y = "Coefficient (year)"
    ) +
    # Add observation counts to the right of points
    geom_text(
      aes(label = paste0("italic(N) == ", n)),
      parse = TRUE,
      y = 0.5,
      hjust = 1.1,
      size = 3
    ) +
    theme_classic() +
    ylim(-0.5, 0.5) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}