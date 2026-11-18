# Run from the project root: Rscript tests/test_fit_model_sparse.R
suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(viridis)
})
for (file in c("fit_model", "extract_coef", "analyse_group", "get_breaks",
               "plot_disturbed", "plot_filtereffects", "plot_lqmm_bygroup",
               "calc_lqmm_byqmdbin")) {
  source(paste0("R/", file, ".R"))
}

set.seed(15)
df <- data.frame(
  plotID = rep(seq_len(40), each = 6),
  QMD = exp(rnorm(240, 3.5, 0.03)),
  year = rep(seq(1980, 2020, by = 8), 40)
)
df$density <- exp(12 - 1.4 * log(df$QMD) + 0.008 * (df$year - 2000) +
  rep(rnorm(40, 0, 0.3), each = 6) + rnorm(240, 0, 0.15))

expect_warning(expect_null(fit_model(df[FALSE, ])), "0 usable rows")
expect_warning(expect_null(fit_model(df[1:29, ])), "29 usable rows")
bad <- df[1:34, ]
bad$QMD[1:5] <- c(0, -1, NA, Inf, NaN)
expect_warning(expect_null(fit_model(bad)), "29 usable rows")
bad <- df[1:34, ]
bad$density[1:3] <- c(0, -1, Inf)
bad$year[4] <- Inf
bad$plotID[5] <- NA
expect_warning(expect_null(fit_model(bad)), "29 usable rows")
expect_warning(expect_null(fit_model(transform(df, year = 2000))), "must both vary")
expect_warning(expect_null(fit_model(transform(df, QMD = 30))), "must both vary")
expect_warning(expect_null(fit_model(transform(df, plotID = 1))), "two plots")
expect_s3_class(fit_model(df[1:30, ]), "lqmm")
expect_true(is.na(extract_coef(NULL)$estimate))

# The strict filters retain just three or zero rows. The valid stages must
# still finish, while the skipped stages retain their counts and NA estimates.
df <- transform(df, disturbed = FALSE, ndisturbed = 0L, badqmdbin = FALSE,
                badslope = FALSE, management_cat = 1L,
                management_since_census1_yrs = c(rep(150, 3), rep(50, 237)))
out <- suppressWarnings(analyse_group(df, "Sparse filters", do_boot = FALSE))
expect_equal(nrow(out$df_filtereffects), 7L)
strict <- filter(out$df_filtereffects, step %in% c("no_mgmt_100", "oldgrowth"))
expect_equal(strict$n, c(3L, 0L))
expect_true(all(is.na(strict$estimate)))
expect_true(all(map_lgl(strict$fit_lqmm, is.null)))
expect_true(is.finite(filter(out$df_filtereffects, step == "orig")$estimate))

# A group that is too small even before filtering skips dependent model plots.
small <- suppressWarnings(analyse_group(df[1:20, ], "Small group", do_boot = TRUE))
expect_null(small$gg_lqmm)
expect_null(small$gg_lqmm_byqmdbin)
expect_null(small$boot_results)
expect_equal(nrow(small$df_filtereffects), 7L)
empty_boot <- suppressWarnings(bootstrap_lqmm(df[1:20, ], 2, workers = 1))
expect_equal(nrow(empty_boot), 0L)
expect_true(all(c("id", "term", "estimate") %in% names(empty_boot)))
cat("Sparse/invalid input guards, 30-row boundary, skipped filter stages and downstream handling passed.\n")
