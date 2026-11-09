# Fit the same model after learning predictor scaling from its training sample.
fit_model <- function(df, lqmm = TRUE, control = list(), min_rows = 30L) {

  # Count usable model inputs, rather than rows that will be lost to missing
  # values or produce non-finite logarithms. A skipped fit returns NULL.
  df <- as.data.frame(df)[, c("density", "QMD", "year", "plotID")]

  usable <- is.finite(df$density) & df$density > 0 &
    is.finite(df$QMD) & df$QMD > 0 & is.finite(df$year) & !is.na(df$plotID)

  df <- df[usable, , drop = FALSE]

  if (nrow(df) < min_rows) {
    warning(sprintf("Skipping model: only %d usable rows (minimum %d).",
                    nrow(df), min_rows), call. = FALSE)
    return(NULL)
  }

  if (length(unique(df$QMD)) < 2L || length(unique(df$year)) < 2L) {
    warning("Skipping model: QMD and year must both vary.", call. = FALSE)
    return(NULL)
  }
  
  if (length(unique(df$plotID)) < 2L) {
    warning("Skipping mixed model: at least two plots are required.", call. = FALSE)
    return(NULL)
  }

  set.seed(123)
  df$logDensity <- log(df$density)
  df$logQMD <- log(df$QMD)
  df <- df[, c("logDensity", "logQMD", "year", "plotID")]

  rec <- recipes::recipe(logDensity ~ logQMD + year + plotID, data = df) |>
    recipes::update_role(plotID, new_role = "id") |>
    recipes::step_normalize(logQMD, year, id = "normalize_predictors") |>
    recipes::step_rename(logQMD_sc = logQMD, year_sc = year, id = "rename_predictors") |>
    recipes::prep(training = df, retain = FALSE)

  df_scaled <- recipes::bake(rec, new_data = df)
  
  fixed_formula <- logDensity ~ logQMD_sc + year_sc

  if (lqmm){

    control <- utils::modifyList(
      list(LP_max_iter = 500, LP_tol_ll = 1e-05, startQR = TRUE),
      control
    )
    fit <- lqmm::lqmm(
      fixed_formula,
      random = ~1,
      group = plotID,
      tau = 0.9,
      data = df_scaled,
      type = "normal",
      control = do.call(lqmm::lqmmControl, control)
    )

  } else {
    
    fit <- lme4::lmer(
      stats::update(fixed_formula, . ~ . + (1 | plotID)),
      data = df_scaled,
      na.action = "na.exclude"
    )
  }

  attr(fit, "recipe") <- rec
  return(fit)
}

model_recipe <- function(fit) {
  rec <- attr(fit, "recipe")
  if (is.null(rec)) {
    stop("Model has no preprocessing recipe; refit it with fit_model().")
  }
  rec
}

# New data must be in original units; reuse the training means and SDs.
predict_model <- function(fit, new_data, ...) {
  df_scaled <- recipes::bake(model_recipe(fit), new_data = new_data)
  stats::predict(fit, newdata = as.data.frame(df_scaled), ...)
}

# Back-transform the fixed effects, including the centered intercept.
coef_original_units <- function(fit) {
  scaling <- recipes::tidy(model_recipe(fit), number = 1)
  means <- scaling[scaling$statistic == "mean", ]
  sds <- scaling[scaling$statistic == "sd", ]
  centers <- stats::setNames(means$value, means$terms)
  scales <- stats::setNames(sds$value, sds$terms)
  beta <- if (inherits(fit, "lqmm")) stats::coef(fit) else lme4::fixef(fit)
  slopes <- beta[paste0(names(scales), "_sc")] / scales

  c(
    "(Intercept)" = unname(beta["(Intercept)"] - sum(slopes * centers[names(scales)])),
    stats::setNames(unname(slopes), names(scales))
  )
}
