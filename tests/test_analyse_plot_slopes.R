# Run from the repository root: Rscript tests/test_analyse_plot_slopes.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})
source("R/fit_model.R")
source("R/analyse_group.R")
set.seed(34)
plot_slopes <- seq(-1.5, -2, length.out = 8)
df <- tibble(
  plotID = rep(letters[1:8], each = 8),
  QMD = exp(rep(seq(2, 4, length.out = 8), 8)),
  density = exp(12 + rep(rnorm(8, sd = 0.35), each = 8) +
                  rep(plot_slopes, each = 8) * log(QMD) + rnorm(64, sd = 0.01)),
  year = rep(2000:2007, 8)
)
r <- analyse_plot_slopes(df, "Example")
scaled <- mutate(df, logDensity = log(density), logQMD_sc = as.numeric(scale(log(QMD))))
set.seed(123)
expected <- lqmm::lqmm(
  logDensity ~ logQMD_sc, random = ~1, group = plotID, tau = 0.9,
  data = as.data.frame(scaled), type = "normal",
  control = lqmm::lqmmControl(LP_max_iter = 500, LP_tol_ll = 1e-05, startQR = TRUE)
)
stopifnot(nrow(r$df_slopes) == 8L, all(r$df_slopes$nobs == 8L),
          max(abs(r$df_slopes$slope_lin - plot_slopes)) < 0.02,
          inherits(r$fit_qr, "lqmm"),
          identical(all.vars(r$fit_qr$mtf), c("logDensity", "logQMD_sc")),
          identical(colnames(r$fit_qr$mmr), "(Intercept)"),
          r$fit_qr$ngroups == length(unique(df$plotID)),
          r$fit_qr$tau == 0.9,
          isTRUE(all.equal(coef(r$fit_qr), coef(expected))),
          isTRUE(all.equal(r$slope_qr, coef(expected)[["logQMD_sc"]] / sd(log(df$QMD)))))
built <- ggplot_build(r$gg_slopes)
median_slope <- median(r$df_slopes$slope_lin)
subtitle <- get_labs(r$gg_slopes)$subtitle
stopifnot(sum(built$data[[1]]$count) == 8,
          identical(built$data[[2]]$xintercept, c(r$slope_qr, median_slope)),
          identical(built$data[[2]]$linetype, c("solid", "dashed")),
          all(built$data[[2]]$colour == "red"),
          grepl(sprintf("Overall quantile regression = %.3f", r$slope_qr), subtitle, fixed = TRUE),
          grepl(sprintf("Median of plot-specific slopes = %.3f", median_slope), subtitle, fixed = TRUE),
          inherits(r$legend_slopes, "gtable"),
          identical(r$gg_slopes$theme$legend.position, "none"))

# Invalid raw values and stale log columns cannot alter either fit.
bad <- df[1:4, ] |> mutate(QMD = c(NA, 0, -1, Inf))
dirty <- bind_rows(df, bad) |> mutate(logN = 99, logQMD = -99)
dirty_before <- dirty
dirty_result <- analyse_plot_slopes(dirty, "Example")
stopifnot(identical(dirty, dirty_before),
          identical(dirty_result$df_slopes, r$df_slopes),
          isTRUE(all.equal(dirty_result$slope_qr, r$slope_qr)))

# Calendar year is irrelevant to the no-year fit, including model preprocessing.
without_year <- analyse_plot_slopes(select(df, -year), "No year column")
missing_year <- analyse_plot_slopes(mutate(df, year = NA_real_), "Missing years")
stopifnot(isTRUE(all.equal(coef(without_year$fit_qr), coef(r$fit_qr))),
          isTRUE(all.equal(coef(missing_year$fit_qr), coef(r$fit_qr))))

# Every input plot remains in the slope table, even when a slope is undefined.
edge <- bind_rows(
  df,
  mutate(df[1:2, ], plotID = "short"),
  mutate(df[1:4, ], plotID = "constant", QMD = 20),
  mutate(df[1:4, ], plotID = "invalid", density = NA_real_)
)
e <- analyse_plot_slopes(edge, "Example")
stopifnot(nrow(e$df_slopes) == 11L,
          all(is.na(filter(e$df_slopes, plotID %in% c("short", "constant", "invalid"))$slope_lin)),
          filter(e$df_slopes, plotID == "invalid")$nobs == 0L)
empty <- suppressWarnings(analyse_plot_slopes(df[FALSE, ], "Empty"))
stopifnot(nrow(empty$df_slopes) == 0L, is.null(empty$fit_qr),
          is.na(empty$slope_qr), length(ggplot_build(empty$gg_slopes)$data) == 1L,
          is.null(empty$legend_slopes))
one_plot <- suppressWarnings(analyse_plot_slopes(mutate(df, plotID = "a"), "One plot"))
stopifnot(is.null(one_plot$fit_qr), is.na(one_plot$slope_qr),
          identical(ggplot_build(one_plot$gg_slopes)$data[[2]]$linetype, "dashed"))

# Exercise analyse_group's wiring without rerunning the unrelated LQMM fits.
get_breaks <- function(...) NULL
plot_disturbed <- function(...) ggplot()
original_fit_model <- fit_model
fit_model <- function(..., include_year = TRUE) {
  if (include_year) return(NULL)
  original_fit_model(..., include_year = FALSE)
}
extract_coef <- function(...) tibble(term = "year", estimate = NA_real_)
plot_filtereffects <- function(...) ggplot()
calc_lqmm_byqmdbin <- function(d) list(df = d)
plot_lqmm_byqmdbin <- function(...) ggplot()
out <- analyse_group(mutate(df, badqmdbin = FALSE, management_cat = 2L,
                            management_since_census1_yrs = NA_real_), "Example")
stopifnot(identical(out$df_slopes, r$df_slopes),
          isTRUE(all.equal(out$slope_qr, r$slope_qr)),
          inherits(out$fit_qr, "lqmm"), inherits(out$gg_slopes, "ggplot"),
          inherits(out$legend_slopes, "gtable"))
cat("Plot slopes, reference lines, subtitle, legend, sparse data and analyse_group outputs passed.\n")
