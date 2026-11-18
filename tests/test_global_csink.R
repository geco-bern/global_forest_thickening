# Run from the project root: Rscript tests/test_global_csink.R
suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(tidyr)
  library(recipes)
  library(lme4)
})
source("R/global_csink.R")

set.seed(17)
n <- 600L
d <- data.frame(
  plotID = rep(seq_len(60), each = 10),
  species = rep(letters[1:6], length.out = n), dataset = "test",
  year = rep(seq(1980, 2025, length.out = 10), 60),
  logQMD = rnorm(n, log(30), 0.3), tavg = runif(n, 0, 25),
  ai = runif(n, 0.5, 2), ndep = runif(n, 0.1, 2),
  ORGC = runif(n, 1, 100), CNrt = runif(n, 8, 20), PBR = runif(n, 1, 20)
)
d$logDensity <- 12 - 1.7 * d$logQMD +
  (d$year - 2000) * (0.01 + 0.0002 * d$tavg) +
  rep(rnorm(60, 0, 0.3), each = 10) +
  rep(rnorm(6, 0, 0.1), length.out = n) + rnorm(n, 0, 0.01)
d$QMD <- exp(d$logQMD)
d$NQMD2 <- exp(d$logDensity) * d$QMD^2
d$biomass <- 0.5 * d$NQMD2 + rep(rnorm(60, 0, 100), each = 10) + rnorm(n)
reference <- fit_stl_byboot(d)
g <- d[1:30, ] |> mutate(lon_i = row_number(), lat_i = 1L, fcf = 0.5, area_ha = 100)
g0 <- mutate(g, year = 2000)
g1 <- mutate(g, year = 2001)
predict_bundle <- function(bundle, newdata) {
  predict(bundle$model, newdata = bake(bundle$recipe, new_data = newdata), re.form = NA)
}

# Use identical bootstrap rows but different centering/scaling constants.
# Interactions and both prediction years are included in this invariance test.
for (seed in 1:3) {
  set.seed(seed)
  boot <- d[sample.int(n, n, replace = TRUE), ]
  by_boot <- fit_stl_byboot(boot)
  fixed <- list(recipe = reference$recipe, model = lmer(
    formula(reference$model), data = bake(reference$recipe, new_data = boot)
  ))
  expect_equal(predict_bundle(by_boot, g0), predict_bundle(fixed, g0), tolerance = 1e-5)
  expect_equal(predict_bundle(by_boot, g1) - predict_bundle(by_boot, g0),
               predict_bundle(fixed, g1) - predict_bundle(fixed, g0), tolerance = 1e-5)
}

# A positive annual log-density change and positive biomass conversion cannot
# produce a negative integrated sink when forest areas are nonnegative.
delta <- predict_bundle(reference, g1) - predict_bundle(reference, g0)
expect_true(all(delta > 0))
p <- predict_dn(reference, rep(30, 2), g)
g0$logQMD <- log(30)
g1$logQMD <- log(30)
expect_equal(unname(p$dn), unname(exp(predict_bundle(reference, g1)) -
               exp(predict_bundle(reference, g0))), tolerance = 1e-10)
bm <- fit_biomass_byboot(d)
expect_gt(unname(fixef(bm)["NQMD2"]), 0)
p <- predict_db(p, bm)
expect_true(all(p$db > 0))
expect_gt(calc_csink_global(p, g), 0)
expect_equal(calc_csink_global(p, g),
             sum(p$db * g$fcf * g$area_ha) * 0.5e-12)
expect_equal(sum(calc_csink_bylat(p, g)$db_tgc) / 1000, calc_csink_global(p, g))
expect_true(is.na(calc_csink_global(mutate(p, db = NA_real_), g)))
expect_error(calc_csink_global(mutate(p, db = Inf), g), "Non-finite")
cat("Global carbon sink tests passed.\n")
