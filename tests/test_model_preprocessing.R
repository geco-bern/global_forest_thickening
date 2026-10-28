library(testthat)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(ggplot2)
library(viridis)

source("R/fit_model.R")
source("R/calc_percent_change.R")
source("R/extract_coef.R")
source("R/calc_lqmm_byqmdbin.R")
source("R/plot_lqmm_bybiome.R")
source("R/plot_disturbed.R")
source("R/get_breaks.R")
source("R/analyse_biome.R")

set.seed(15)
df <- data.frame(
  plotID = rep(seq_len(40), each = 6),
  logQMD = rnorm(240, 3.5, 0.3),
  year = rep(seq(1980, 2020, by = 8), 40)
)
df$logDensity <- 12 - 1.4 * df$logQMD + 0.008 * (df$year - 2000) +
  rep(rnorm(40, 0, 0.3), each = 6) + rnorm(240, 0, 0.15)
original_df <- df
fit <- fit_model(df)
biome_df <- transform(
  df,
  disturbed = plotID %% 10 == 0 & (year <= 2000 | plotID %% 20 == 0),
  ndisturbed = as.integer(plotID %% 10 == 0),
  badqmdbin = plotID %% 7 == 0,
  badslope = plotID %% 13 == 0,
  management_cat = ifelse(plotID %% 3 == 0, 2, 1),
  management_since_census1_yrs = 120
)

test_that("scaling is local to the model and ignores stale scaled input columns", {
  expect_identical(df, original_df)
  baked <- recipes::bake(model_recipe(fit), new_data = df)
  for (variable in c("logQMD", "year")) {
    scaled <- baked[[paste0(variable, "_sc")]]
    expect_equal(mean(scaled), 0, tolerance = 1e-12)
    expect_equal(sd(scaled), 1, tolerance = 1e-12)
    expect_equal(scaled, as.numeric(scale(df[[variable]])))
  }
  expect_identical(baked$logDensity, df$logDensity)
  expect_identical(baked$plotID, df$plotID)
  expect_false(model_recipe(fit)$retained)

  stale <- transform(df, logQMD_sc = -999, year_sc = 999)
  expect_equal(coef(fit_model(stale)), coef(fit))
})

test_that("a bootstrap of the full sample reproduces the ordinary fit", {
  local_mocked_bindings(
    bootstraps = function(data, ...) {
      split <- rsample::make_splits(
        list(analysis = seq_len(nrow(data)), assessment = integer()), data = data
      )
      rsample::manual_rset(list(split), ids = "Apparent")
    },
    .package = "rsample"
  )
  analysis <- analyse_biome(biome_df, "Test biome", n_boot = 1)
  result <- analysis$boot_results
  model <- analysis$fit_lqmm[[1]]
  selected <- biome_df |>
    filter(ndisturbed == 0, !badqmdbin, !badslope)
  beta <- coef_original_units(model)
  expect_identical(result$term, names(beta))
  expect_equal(result$estimate, unname(beta))
  expect_equal(beta[["year"]], coef(model)[["year_sc"]] / sd(selected$year))
  expect_equal(beta[["logQMD"]], coef(model)[["logQMD_sc"]] / sd(selected$logQMD))
})

test_that("scaling is learned only from complete model rows", {
  training <- df
  training$logDensity[1] <- NA_real_
  training$logQMD[2] <- NA_real_
  model <- fit_model(training)
  complete <- training[complete.cases(training), ]
  expected <- fit_model(complete)
  expect_equal(coef(model), coef(expected))
  expect_equal(coef_original_units(model), coef_original_units(expected))
})

test_that("predictions reuse training scaling and survive serialization", {
  new_data <- data.frame(
    plotID = c(100L, 101L), logQMD = c(3.2, 4.0), year = c(1990, 2050)
  )
  pred <- predict_model(fit, new_data, level = 0)
  beta <- coef_original_units(fit)
  expected <- as.matrix(cbind(1, new_data[c("logQMD", "year")])) %*% beta
  expect_equal(unname(pred), as.numeric(expected), tolerance = 1e-10)
  expect_equal(
    unname(predict_model(fit, new_data[1, ], level = 0)), unname(pred[1])
  )

  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  unlink(path)
  expect_equal(predict_model(restored, new_data, level = 0), pred)
  expect_equal(coef_original_units(restored), beta)
})

test_that("annual percentage change is expressed per calendar year", {
  expected <- 100 * expm1(coef(fit)[["year_sc"]] / sd(df$year))
  expect_equal(unname(calc_percent_change(df, fit)), expected, tolerance = 1e-10)
})

test_that("the linear mixed model shares preprocessing", {
  model <- fit_model(df, lqmm = FALSE)
  beta <- coef_original_units(model)
  expected <- as.matrix(cbind(1, df[c("logQMD", "year")])) %*% beta
  expect_equal(
    unname(predict_model(model, df, re.form = NA)),
    as.numeric(expected), tolerance = 1e-10
  )
})

test_that("QMD-bin fits use the shared model and retain raw data", {
  result <- calc_lqmm_byqmdbin(df, breaks = c(0, 10))$df
  expect_equal(nrow(result), 1)
  model <- result$mod_lqmm[[1]]
  expect_s3_class(model_recipe(model), "recipe")
  expect_equal(model$control$LP_max_iter, 6000)
  expected <- fit_model(df, control = list(LP_max_iter = 6000, UP_max_iter = 6000))
  expect_equal(coef(model), coef(expected))
  expect_false(any(c("logQMD_sc", "year_sc") %in% names(result$data[[1]])))
})

test_that("biome bootstrap samples use no_badslope rows and their own scaling", {
  bootstrap <- rsample::bootstraps
  resamples <- NULL
  local_mocked_bindings(
    bootstraps = function(data, ...) {
      resamples <<- bootstrap(data, ...)
      resamples
    },
    .package = "rsample"
  )
  result <- analyse_biome(biome_df, "Test biome", n_boot = 3)
  selected <- biome_df |>
    filter(ndisturbed == 0, !badqmdbin, !badslope)
  expect_equal(coef(result$fit_lqmm[[1]]), coef(fit_model(selected)))
  expect_equal(nrow(result$df_boot), 3)
  expect_true(all(result$df_boot$term == "year"))
  expect_true(all(result$df_boot$biome == "Test biome"))

  for (i in seq_len(nrow(resamples))) {
    sample <- rsample::analysis(resamples$splits[[i]])
    expect_true(all(sample$ndisturbed == 0 & !sample$badqmdbin & !sample$badslope))
    model <- fit_model(sample)
    beta <- coef_original_units(model)
    actual <- result$boot_results |>
      filter(id == resamples$id[i])
    expect_equal(actual$estimate, unname(beta))
    expect_equal(beta[["year"]], coef(model)[["year_sc"]] / sd(sample$year))
  }

  lines <- result$gg_lqmm$layers[[2]]$data
  expect_equal(
    lines$pred_90,
    predict_model(result$fit_lqmm[[1]], transform(lines, year = as.numeric(as.character(year))), level = 0)
  )
})
