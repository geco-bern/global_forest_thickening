# Run from the project root: Rscript tests/test_bootstrap_lqmm.R
suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
})
source("R/fit_model.R")
source("R/analyse_group.R")

set.seed(15)
df <- data.frame(
  plotID = rep(seq_len(40), each = 6),
  QMD = exp(rnorm(240, 3.5, 0.3)),
  year = rep(seq(1980, 2020, by = 8), 40)
)
df$density <- exp(12 - 1.4 * log(df$QMD) + 0.008 * (df$year - 2000) +
  rep(rnorm(40, 0, 0.3), each = 6) + rnorm(240, 0, 0.15))
df$unused_column <- "not needed by workers"
original <- df
old_plan <- future::plan()

# Real model fits exercise recipe learning and function/package export to
# independent R processes, as well as reproducibility across worker counts.
serial <- bootstrap_lqmm(df, n_boot = 4, workers = 1, seed = 42)
expect_equal(future::plan(), old_plan)
parallel_result <- bootstrap_lqmm(df, n_boot = 4, workers = 2, seed = 42)
expect_equal(future::plan(), old_plan)
expect_equal(parallel_result, serial, tolerance = 1e-10)
expect_equal(nrow(serial), 12L)
expect_equal(sum(serial$term == "year"), 4L)
expect_true(all(is.finite(serial$estimate)))
expect_identical(df, original)

# Compare with the existing serial fitting approach on the same draws.
splits <- withr::with_seed(42, rsample::bootstraps(
  dplyr::group_by(df, plotID), times = 4, apparent = FALSE
))
reference <- purrr::map_dfr(splits$splits, function(split) {
  model <- fit_model(rsample::analysis(split), lqmm = TRUE)
  beta <- coef_original_units(model)
  tibble::tibble(term = names(beta), estimate = unname(beta))
})
expect_equal(serial$estimate, reference$estimate, tolerance = 1e-8)

# Simulate Positron injecting a base::.ark_auto_step hook into exported
# functions. Ordinary R (including workers) has no such debugger helper.
worker_names <- c(
  "fit_bootstrap_sample", "fit_model", "coef_original_units", "model_recipe"
)
original_functions <- mget(worker_names)
for (name in worker_names) {
  fn <- original_functions[[name]]
  body(fn) <- substitute({ base::.ark_auto_step(NULL); ORIGINAL },
                         list(ORIGINAL = body(fn)))
  assign(name, fn)
}
from_ide <- bootstrap_lqmm(df, n_boot = 4, workers = 2, seed = 42)
expect_equal(from_ide, serial, tolerance = 1e-10)
list2env(original_functions, envir = .GlobalEnv)
expect_equal(future::plan(), old_plan)

# The previous future plan must also be restored if mapping itself fails.
saved_fit <- fit_bootstrap_sample
fit_bootstrap_sample <- function(indices, data) stop("test worker failure")
expect_error(bootstrap_lqmm(df, 2, workers = 2), "test worker failure")
fit_bootstrap_sample <- saved_fit
expect_equal(future::plan(), old_plan)

cat("Serial/parallel estimates, Positron-hook regression, input preservation and plan cleanup passed.\n")
