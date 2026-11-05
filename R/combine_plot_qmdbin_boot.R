combine_plot_qmdbin_boot <- function(gg1, gg2, group_name){

  gg_lqmm_combined <- cowplot::plot_grid(
    gg1,
    gg2,
    nrow = 1,
    align = "h",
    axis = "tb"
  )

  title <- ggdraw() +
    draw_label(
      group_name,
      fontface = 'bold',
      size = 14,
      hjust = 0,
      # x = 0.5
    )

  cowplot::plot_grid(
    title,
    gg_lqmm_combined,
    ncol = 1,
    rel_heights = c(0.1, 1) # Adjust height ratio so the title isn't squished
  )
}