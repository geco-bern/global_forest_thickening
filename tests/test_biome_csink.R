# Run from the project root: Rscript tests/test_biome_csink.R
suppressPackageStartupMessages({library(testthat); library(dplyr)})
source("R/fit_model.R")
source("R/biome_csink.R")

# Heterogeneous QMD: test the full nonlinear expectation, not a mean-QMD shortcut.
beta <- c("(Intercept)" = -12, logQMD = -1.5, year = 0.01)
qmd <- c(10, 20, 60)
a <- 0.6
area <- 2e8
prediction <- predict_biome_quantile_sink(beta, a, qmd, area)
n0 <- exp(beta[1] + beta[2] * log(qmd) + beta[3] * 2000)
n1 <- exp(beta[1] + beta[2] * log(qmd) + beta[3] * 2001)
db <- mean(a * qmd^2 * (n1 - n0))
expect_equal(prediction$db_kgdm_ha_yr, unname(db), tolerance = 1e-10)
expect_equal(prediction$sink_pgc_yr, unname(db * 0.5 * area * 1e-12))
expect_equal(prediction$dc_mgc_ha_yr, unname(db * 0.5e-3))
expect_gt(prediction$sink_pgc_yr, 0)
expect_equal(predict_biome_quantile_sink(replace(beta, 3, 0), a, qmd, area)$sink_pgc_yr, 0)
expect_lt(predict_biome_quantile_sink(replace(beta, 3, -0.01), a, qmd, area)$sink_pgc_yr, 0)
expect_equal(predict_biome_quantile_sink(beta, a, qmd, area * 2)$sink_pgc_yr,
             2 * prediction$sink_pgc_yr)
expect_error(predict_biome_quantile_sink(beta, a, qmd, -1), "Invalid")

grid <- data.frame(lon = 1:8, lat = 0, biome_number = c(1, 2, 4, 5, 6, 12, 3, 6),
                   biome_major = c(rep("Tropical & Subtropical Broadleaf Forests", 2),
                                   rep("Temperate Forests", 2), "Boreal Forests/Taiga",
                                   "Mediterranean Forests", "Excluded", "Boreal Forests/Taiga"),
                   area_ha = 100, fcf = c(0.5, 0.2, 0.3, 0, 0.4, 0.1, 1, NA))
mapped <- prepare_biome_sink_grid(grid)
expect_equal(nrow(mapped$grid), 7L)
expect_equal(sum(mapped$area$forest_area_ha), 150)
expect_equal(sum(mapped$area$n_missing_cover), 1L)
expect_error(prepare_biome_sink_grid(bind_rows(grid, grid[1, ])), "Duplicate")
expect_error(prepare_biome_sink_grid(mutate(grid, fcf = 50)), "fraction")

# Exercise actual LQMM recipes, original-unit prediction and joint resampling.
set.seed(8)
d <- data.frame(plotID = rep(1:40, each = 6), year = rep(seq(1980, 2020, by = 8), 40),
                QMD = exp(rnorm(240, 3.3, 0.35)))
d$density <- exp(12 - 1.6 * log(d$QMD) + 0.01 * (d$year - 2000) +
                  rep(rnorm(40, 0, 0.2), each = 6) + rnorm(240, 0, 0.05))
d$biomass <- 0.55 * d$density * d$QMD^2 + rep(rnorm(40, 0, 500), each = 6) + rnorm(240, 0, 20)
fit <- fit_model(d)
saved <- tibble(biome_major = "Temperate Forests", data = list(d),
                out = list(list(fit_lqmm = list(fit))))
inputs <- prepare_biome_sink_inputs(saved)
new0 <- mutate(d, year = 2000, logQMD = log(QMD))
new1 <- mutate(new0, year = 2001)
ref0 <- predict_model(fit, new0, level = 0)
ref1 <- predict_model(fit, new1, level = 0)
a <- fit_biome_biomass(d)
expect_equal(a, 0.55, tolerance = 0.01)
point <- estimate_biome_sink(d, area, fit)
expect_equal(point$db_kgdm_ha_yr, mean(a * d$QMD^2 * (exp(ref1) - exp(ref0))),
             tolerance = 1e-8)

areas <- tibble(biome_major = names(inputs), forest_area_ha = area)
old_plan <- future::plan()
serial <- bootstrap_biome_sink(inputs, areas, nboot = 3, workers = 1)
parallel <- bootstrap_biome_sink(inputs, areas, nboot = 3, workers = 2)
expect_equal(future::plan(), old_plan)
expect_equal(serial, parallel, tolerance = 1e-10)
expect_equal(nrow(serial), 3L)
expect_true(all(is.finite(serial$sink_pgc_yr)))
expect_gt(length(unique(serial$beta_year)), 1L)
expect_equal(sum_biome_sink_draws(serial, names(inputs))$sink_pgc_yr, serial$sink_pgc_yr)

missing_biome <- bind_rows(serial, mutate(serial, biome_major = "Boreal Forests/Taiga",
                                         sink_pgc_yr = NA_real_))
global <- sum_biome_sink_draws(missing_biome, unique(missing_biome$biome_major))
expect_true(all(is.na(global$sink_pgc_yr)))
expect_equal(summarise_biome_sink(global)$n_success, 0L)
expect_true(is.na(summarise_biome_sink(global)$mean))
inputs[[1]]$data$biomass <- NA_real_
expect_warning(failed <- bootstrap_biome_sink(inputs, areas, nboot = 1), "Failed biome")
expect_true(is.na(failed$sink_pgc_yr))
expect_match(failed$error, "Insufficient biomass")
cat("Biome quantile sink tests passed.\n")
