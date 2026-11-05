plot_boot <- function(df_boot) {
  ggplot(
    aes(
      x = estimate,
      y = after_stat(density)
    ),
    data = df_boot
  ) +
    geom_histogram(bins = 30, fill = "grey70", color = "black") +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_x_continuous(limits = symmetric_limits) +
    coord_flip()
}

symmetric_limits <- function(x) {
  bound <- max(abs(x), na.rm = TRUE)
  c(-bound, bound)
}
