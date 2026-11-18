# Run from the project root: Rscript tests/test_tree_density_contributions.R

# Load packages and functions --------------------------------------------------
library(testthat)
library(dplyr)
library(purrr)

source("R/tree_density_contributions.R")

# Construct paired censuses with exactly known class contributions -------------
make_test_pairs <- function(n = 10L) {
  x <- seq(-0.2, 0.2, length.out = n)

  df <- tibble(
    dataset = "test",
    biome_major = "biome",
    plotID = paste0("test_", seq_len(n)),
    cutoff_cm = 0,
    year0 = 2000,
    year1 = 2010,
    years = 10,
    qmd_change_pct_decade = 100 * expm1(x),
    n_small_first = 100,
    n_medium_first = 60,
    n_large_first = 40,
    n_verylarge_first = 20,
    n_small_last = 100 + 5 + 2 * x,
    n_medium_last = 60 + 2 + 3 * x,
    n_large_last = 40 + 1 + x,
    n_verylarge_last = 20 + 0.5 + 0.5 * x
  ) |>
    mutate(
      N_first = n_small_first + n_medium_first + n_large_first,
      N_last = n_small_last + n_medium_last + n_large_last
    )

  for (component in c("N", "n_small", "n_medium", "n_large", "n_verylarge")) {
    df[[paste0(component, "_change_decade")]] <-
      10 * (df[[paste0(component, "_last")]] - df[[paste0(component, "_first")]]) /
      df$years
  }

  df
}

test_that("class changes reproduce the total and use the correct interval", {
  pairs <- make_test_pairs()
  df <- density_contribution_data(pairs)

  expect_true(all(df$increasing))
  expect_equal(df$dlogqmd_decade, seq(-0.2, 0.2, length.out = 10))
  expect_equal(
    df$N_change_decade,
    df$n_small_change_decade + df$n_medium_change_decade + df$n_large_change_decade
  )

  expect_error(
    density_contribution_data(mutate(pairs, N_last = N_last + 1)),
    "do not sum"
  )
  expect_error(
    density_contribution_data(mutate(pairs, n_verylarge_first = n_large_first + 1)),
    "overlapping"
  )
  expect_error(
    density_contribution_data(mutate(pairs, years = 20)),
    "annualisation"
  )
  expect_error(
    density_contribution_data(mutate(pairs, n_small_last = NA_real_)),
    "Non-finite"
  )
  expect_error(density_contribution_data(bind_rows(pairs, pairs)), "Duplicate")
})

test_that("movement between classes alone is not a total-density increase", {
  pairs <- make_test_pairs(1L)
  changes <- c(n_small = -10, n_medium = 10, n_large = 0, n_verylarge = 0, N = 0)

  for (component in names(changes)) {
    pairs[[paste0(component, "_last")]] <-
      pairs[[paste0(component, "_first")]] + changes[[component]]
    pairs[[paste0(component, "_change_decade")]] <- changes[[component]]
  }

  df <- density_contribution_data(pairs)
  expect_false(df$increasing)
  expect_equal(nrow(density_source_estimate(df, "increasing_plots")), 0L)
})

test_that("mean contributions remain additive and can exceed 100 percent", {
  df <- density_contribution_data(make_test_pairs())
  estimate <- density_source_estimate(df, "increasing_plots")

  expect_equal(estimate$total, estimate$small + estimate$medium + estimate$large)
  expect_equal(estimate$small_minus_large, 4)
  expect_equal(estimate$small_share, 5 / 8)
  expect_true(estimate$included_in_test)
  expect_false(density_source_estimate(slice(df, 1:9))$included_in_test)

  # Replicating plots does not change a source's contrast or give it more votes.
  replicated <- density_source_estimate(bind_rows(df, df))
  expect_equal(replicated$small_minus_large, estimate$small_minus_large)

  offsetting <- df |>
    mutate(n_large_change_decade = -1, N_change_decade = 6)
  shares <- density_source_estimate(offsetting)
  expect_lt(shares$large_share, 0)
  expect_gt(shares$small_share + shares$medium_share, 1)
})

test_that("fixed-QMD contributions recover known intercepts without extrapolation", {
  df <- density_contribution_data(make_test_pairs())
  adjusted <- density_source_estimate(df, "qmd_adjusted")

  expect_equal(adjusted$small, 5)
  expect_equal(adjusted$medium, 2)
  expect_equal(adjusted$large, 1)
  expect_equal(adjusted$total, 8)
  expect_true(adjusted$qmd_zero_supported)
  expect_true(adjusted$included_in_test)

  unsupported <- df |>
    mutate(dlogqmd_decade = dlogqmd_decade + 1) |>
    density_source_estimate(mode = "qmd_adjusted")

  expect_false(unsupported$qmd_zero_supported)
  expect_false(unsupported$included_in_test)
})

test_that("the exact test uses one sign per source and excludes ties", {
  contrasts <- c(1:6, -1)
  result <- density_contrast_test(contrasts, nboot = 99L)
  with_ties <- density_contrast_test(c(contrasts, 0, NA), nboot = 99L)

  expect_equal(result$n_sources, 7L)
  expect_equal(result$n_small_dominant, 6L)
  expect_equal(result$sign_p_two_sided, 0.125)
  expect_equal(with_ties$n_ties, 1L)
  expect_equal(with_ties$sign_p_two_sided, result$sign_p_two_sided)
  expect_true(is.na(density_contrast_test(c(0, 0), nboot = 99L)$sign_p_two_sided))

  first <- withr::with_seed(42, density_contrast_test(contrasts, nboot = 99L))
  second <- withr::with_seed(42, density_contrast_test(contrasts, nboot = 99L))
  expect_equal(first, second)
})

cat("Tree-density contribution tests passed.\n")
