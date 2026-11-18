# Run from the project root: Rscript tests/test_plot_lqmm_bygroup.R
suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})
source("R/fit_model.R")
source("R/plot_lqmm_bygroup.R")

set.seed(15)
df <- data.frame(
  plotID = rep(seq_len(40), each = 6),
  QMD = exp(rnorm(240, 3.5, 0.3)),
  year = rep(seq(1980, 2020, by = 8), 40)
)
df$density <- exp(12 - 1.4 * log(df$QMD) + 0.008 * (df$year - 2000) +
  rep(rnorm(40, 0, 0.3), each = 6) + rnorm(240, 0, 0.15))
original <- df
fit <- fit_model(df)
p <- plot_lqmm_bygroup(df, fit)

# Expected predictions from unscaled coefficients, independent of plot code.
beta <- coef_original_units(fit)
lines <- p$layers[[2]]$data
expected <- beta[["(Intercept)"]] + beta[["logQMD"]] * lines$logQMD +
  beta[["year"]] * as.numeric(as.character(lines$year))
stopifnot(identical(df, original), !"logQMD" %in% names(df))
expect_equal(unname(lines$pred_90), unname(expected), tolerance = 1e-10)
expect_equal(p$layers[[1]]$data$logQMD, log(df$QMD))
expect_equal(p$layers[[1]]$data$logDensity, log(df$density))
expect_equal(nrow(lines), 300L)
expect_equal(range(lines$logQMD), range(log(df$QMD)))
# The same point has the same prediction alone: no scaling learned from the grid.
one <- as.data.frame(lines[1, c("logQMD", "year", "plotID")])
one$year <- as.numeric(as.character(one$year))
expect_equal(unname(predict_model(fit, one, level = 0)), unname(lines$pred_90[1]),
             tolerance = 1e-10)

# Raw columns take precedence over stale derived columns; unplottable rows
# cannot create infinite limits or change the prediction grid.
bad <- df[1:3, ]
bad$QMD <- c(NA, 0, -1)
dirty <- bind_rows(df, bad) |>
  mutate(logQMD = 999, logDensity = -999)
p_dirty <- plot_lqmm_bygroup(dirty, fit)
expect_equal(p_dirty$layers[[2]]$data, lines)
expect_equal(nrow(p_dirty$layers[[1]]$data), nrow(df))
expect_error(plot_lqmm_bygroup(bad, fit), "No positive, finite")
# Build both layers to check that aesthetics resolve against the local columns.
built <- suppressWarnings(ggplot_build(p))
expect_equal(nrow(built$data[[2]]), 300L)
cat("Raw-input plot, recipe-scaled predictions, log axes, ID types and invalid-value checks passed.\n")
