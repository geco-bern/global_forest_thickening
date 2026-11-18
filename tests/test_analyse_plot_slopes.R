# Run from the repository root: Rscript tests/test_analyse_plot_slopes.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})
source("R/analyse_group.R")
set.seed(34)
df <- tibble(
  plotID = rep(c("a", "b"), each = 8),
  QMD = exp(rep(seq(2, 4, length.out = 8), 2)),
  density = exp(12 + rep(c(-1.5, -2), each = 8) * log(QMD) + rnorm(16, sd = 0.01)),
  year = rep(2000:2007, 2)
)
r <- analyse_plot_slopes(df, "Example")
stopifnot(nrow(r$df_slopes) == 2L, all(r$df_slopes$nobs == 8L),
          max(abs(r$df_slopes$slope_lin - c(-1.5, -2))) < 0.02,
          identical(all.vars(formula(r$fit_qr)), c("logN", "logQMD")),
          r$fit_qr$tau == 0.9,
          isTRUE(all.equal(r$slope_qr, unname(coef(quantreg::rq(
            log(density) ~ log(QMD), tau = 0.9, data = df))[2]))))
built <- ggplot_build(r$gg_slopes)
stopifnot(sum(built$data[[1]]$count) == 2,
          all(built$data[[2]]$xintercept == r$slope_qr))

# Invalid raw values and stale log columns cannot alter either fit.
bad <- df[1:4, ] |> mutate(QMD = c(NA, 0, -1, Inf))
dirty <- bind_rows(df, bad) |> mutate(logN = 99, logQMD = -99)
dirty_before <- dirty
dirty_result <- analyse_plot_slopes(dirty, "Example")
stopifnot(identical(dirty, dirty_before),
          identical(dirty_result$df_slopes, r$df_slopes),
          isTRUE(all.equal(dirty_result$slope_qr, r$slope_qr)))

# Every input plot remains in the slope table, even when a slope is undefined.
edge <- bind_rows(
  df,
  mutate(df[1:2, ], plotID = "short"),
  mutate(df[1:4, ], plotID = "constant", QMD = 20),
  mutate(df[1:4, ], plotID = "invalid", density = NA_real_)
)
e <- analyse_plot_slopes(edge, "Example")
stopifnot(nrow(e$df_slopes) == 5L,
          all(is.na(filter(e$df_slopes, plotID %in% c("short", "constant", "invalid"))$slope_lin)),
          filter(e$df_slopes, plotID == "invalid")$nobs == 0L)
empty <- analyse_plot_slopes(df[FALSE, ], "Empty")
stopifnot(nrow(empty$df_slopes) == 0L, is.null(empty$fit_qr),
          is.na(empty$slope_qr), length(ggplot_build(empty$gg_slopes)$data) == 1L)

# Exercise analyse_group's wiring without rerunning the unrelated LQMM fits.
get_breaks <- function(...) NULL
plot_disturbed <- function(...) ggplot()
fit_model <- function(...) NULL
extract_coef <- function(...) tibble(term = "year", estimate = NA_real_)
plot_filtereffects <- function(...) ggplot()
calc_lqmm_byqmdbin <- function(d) list(df = d)
plot_lqmm_byqmdbin <- function(...) ggplot()
out <- analyse_group(mutate(df, management_cat = 2L,
                            management_since_census1_yrs = NA_real_), "Example")
stopifnot(identical(out$df_slopes, r$df_slopes),
          isTRUE(all.equal(out$slope_qr, r$slope_qr)),
          inherits(out$fit_qr, "rq"), inherits(out$gg_slopes, "ggplot"))
cat("Plot slopes, pooled quantile fit, histogram line, invalid/sparse data and analyse_group outputs passed.\n")
