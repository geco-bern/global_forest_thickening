select_stl_examples <- function(censuses) {
  if (nrow(censuses) == 0L) return(dplyr::select(censuses, dataset, plotID))

  eligible <- censuses |>
    dplyr::mutate(log_qmd = log10(QMD), log_density = log10(density)) |>
    dplyr::group_by(dataset, plotID) |>
    dplyr::filter(
      dplyr::n_distinct(year) >= 3L,
      diff(range(year)) >= 10,
      diff(range(log_qmd)) >= 0.12,
      stats::sd(log_density) > 0
    )
  if (nrow(eligible) == 0L) return(dplyr::select(dplyr::ungroup(eligible), dataset, plotID))

  candidates <- eligible |>
    dplyr::summarise(
      qmd_min = min(log_qmd),
      qmd_max = max(log_qmd),
      qmd_mean = mean(log_qmd),
      density_mean = mean(log_density),
      slope = stats::cov(log_qmd, log_density) / stats::var(log_qmd),
      density_trend = stats::cor(log_qmd, log_density),
      qmd_trend = stats::cor(year, log_qmd),
      .groups = "drop"
    ) |>
    dplyr::filter(
      dplyr::between(slope, -2.5, -0.65), density_trend <= -0.90, qmd_trend >= 0.85
    )

  # Compare fitted densities where trajectories overlap, spanning the 5th to 95th percentile.
  sets <- lapply(seq(1.05, 1.95, by = 0.025), function(ref) {
    overlapping <- candidates |>
      dplyr::filter(qmd_min + 0.05 <= ref, qmd_max - 0.05 >= ref) |>
      dplyr::mutate(density_at_ref = density_mean + slope * (ref - qmd_mean)) |>
      dplyr::arrange(density_at_ref, dataset, plotID)
    n <- nrow(overlapping)
    ranks <- round(seq(1 + 0.05 * (n - 1), 1 + 0.95 * (n - 1), length.out = min(5L, n)))
    overlapping[ranks, ]
  })

  # Prefer five examples, then maximize the smallest gap between adjacent fits.
  sizes <- vapply(sets, nrow, integer(1))
  gaps <- vapply(sets, function(examples) {
    if (nrow(examples) < 2L) return(0)
    min(diff(examples$density_at_ref))
  }, numeric(1))
  sets[[order(sizes, gaps, decreasing = TRUE)[1L]]] |>
    dplyr::select(dataset, plotID)
}

plot_stl_longplots <- function(df, biome_name = "") {

  df <- df |>
    dplyr::ungroup() |>
    dplyr::filter(is.finite(QMD), QMD > 0, is.finite(density), density > 0)
    # dplyr::filter(
    #   management_cat == 2 |
    #     (management_cat == 1 & management_since_census1_yrs >= 30)
    # )

  if (biome_name != ""){
    df <- df |>
      dplyr::filter(biome_major == biome_name)
  }

  visible_censuses <- df |>
    dplyr::filter(
      is.finite(year), !is.na(plotID), !is.na(dataset),
      dplyr::between(QMD, 10, 100), dplyr::between(density, 20, 11000)
    )
  examples <- select_stl_examples(visible_censuses)
  df_longplots <- visible_censuses |>
    dplyr::semi_join(examples, by = c("dataset", "plotID")) |>
    dplyr::arrange(year)

  gg <- ggplot() +
    geom_hex(
      aes(x = QMD, y = density, fill = after_stat(count)),
      data = df,
      binwidth = c(0.02, 0.055)
    ) +
    # geom_point(
    #   aes(x = QMD, y = density),
    #   data = df,
    #   color = "grey80",
    #   # alpha = 0.1,
    #   size = 1,
    #   shape = 16,
    #   inherit.aes = FALSE
    # ) +
    geom_path(
      aes(x = QMD, y = density, group = interaction(dataset, plotID)),
      data = df_longplots,
      inherit.aes = FALSE,
      color = "red"
    ) +
    labs(
      title = biome_name,
      x = "QMD (cm)",
      y = expression(italic(N))
    ) +
    scale_color_viridis_d(name = "Dataset") +
    theme_classic() +
    scale_x_log10(limits = c(10, 100)) +
    scale_y_log10(limits = c(20, 11000)) +
    scale_fill_gradient(
      low = "#EFEFEF",
      high = "#505050",
      transform = "log10",
      guide = "none"
    )

  return(gg)
}
