# Fit model for each filter stage
fit_model <- function(df, lqmm = TRUE) {
  set.seed(123)

  if (lqmm){

    fit <- lqmm(
      logDensity ~ logQMD_sc + year_sc,
      random = ~1,
      group = plotID,
      tau = 0.9,
      data = df,
      type = "normal",
      control = lqmmControl(
        LP_max_iter = 500, # xxx TEST: revert back to 500,
        LP_tol_ll   = 1e-05,
        startQR     = TRUE
      )
    )

  } else {
    
    # not yet tested
    fit <- lmer(
      logDensity ~ logQMD_sc + year_sc +
        # (1 | dataset / plotID) + (1 | species),
        (1 | plotID), # + (1 | species),
      data = df,
      na.action = "na.exclude"
    )
  }

  return(fit)
}