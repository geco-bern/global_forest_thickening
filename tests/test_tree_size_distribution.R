# Run from the project root: Rscript tests/test_tree_size_distribution.R

# Load packages and functions --------------------------------------------------
library(testthat)
library(dplyr)
source("R/tree_size_distribution.R")

# Weighted diameter statistics -------------------------------------------------
# Integer expansion weights allow comparison with an explicit expanded sample.
d <- c(10, 20, 40)
w <- c(1, 2, 1)
m <- size_distribution_metrics(d, w)
expanded <- rep(d, w)

expect_equal(m$N, 4)
expect_equal(m$QMD_tree, sqrt(mean(expanded^2)))
expect_equal(m$BA, sum(pi * (expanded / 200)^2))
expect_equal(m$QMD_tree^2, m$mean_dbh^2 * (1 + m$cv_dbh^2))
expect_equal(
  m$gini_dbh,
  mean(abs(outer(expanded, expanded, "-"))) / (2 * mean(expanded))
)
expect_equal(
  unname(weighted_diameter_quantile(d, w, c(0.1, 0.5, 0.9))),
  c(10, 20, 40)
)
expect_equal(m$ba_top10, 0.4 * 40^2 / sum(w * d^2))
expect_equal(m$N, m$n_small + m$n_medium + m$n_large)

# Shape invariance and missing values ------------------------------------------
expect_equal(size_distribution_metrics(d, w, 20)$N, 3)
expect_equal(size_distribution_metrics(d, w * 10)$gini_dbh, m$gini_dbh)
expect_equal(size_distribution_metrics(d * 2, w)$cv_dbh, m$cv_dbh)
expect_equal(size_distribution_metrics(d * 2, w)$d90_qmd, m$d90_qmd)
expect_equal(size_distribution_metrics(rep(20, 10), rep(1, 10))$gini_dbh, 0)
expect_equal(size_distribution_metrics(rep(20, 10), rep(1, 10))$ba_top10, 0.1)
expect_equal(names(size_distribution_metrics(c(NA, -1), w[1:2])), names(m))
expect_equal(size_distribution_metrics(c(NA, 10), c(1, 1))$missing_fraction, 0.5)

# Match selected plots and expand to per-hectare density ------------------------
trees <- tibble(
  plotID = c("A", "A", "B"),
  year = 2000,
  dbh = c(10, 20, 80),
  plotsize = 0.1
)

stand <- tibble(dataset = "test", plotID = "test_A", year = 2000)
s <- summarise_selected_trees(trees, "test", stand, 0)

expect_equal(nrow(s), 1L)
expect_equal(s$N, 20)
expect_equal(s$QMD_tree, sqrt(250))

# Recover an exactly known STL displacement ------------------------------------
# Movement along N proportional to QMD^b gives zero upward shift. Adding a
# change in the intercept should recover that shift on the original log scale.
x <- bind_rows(m, m) |>
  mutate(
    dataset = "test",
    biome_major = "biome",
    plotID = "A",
    cutoff_cm = 0,
    year = c(2000, 2010),
    QMD = c(20, 25),
    density = 100 * (QMD / 20)^(-1.5),
    n_stems = 20,
    n_effective = 20,
    native_agreement = TRUE,
    management_cat = 2,
    management_since_census1_yrs = NA_real_
  )

p <- make_size_pairs(x, c(biome = -1.5))
expect_equal(p$stl_shift_log_decade, 0, tolerance = 1e-12)

x$density[2] <- x$density[2] * exp(0.1)
expect_equal(make_size_pairs(x, c(biome = -1.5))$stl_shift_log_decade, 0.1)

# A failed reconstruction audit must exclude that census from paired inference.
x$native_agreement[2] <- FALSE
expect_equal(nrow(make_size_pairs(x, c(biome = -1.5))), 0L)

cat("Tree-size distribution tests passed.\n")
