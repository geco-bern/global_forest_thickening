calc_percent_change <- function(df, fit){

  year_mean <- mean(df$year, na.rm = TRUE)
  new_data <- data.frame(
    logQMD = mean(df$logQMD, na.rm = TRUE),
    year = c(year_mean, year_mean + 1),
    plotID = df$plotID[1]
  )
  pred <- predict_model(fit, new_data, level = 0)
  100 * expm1(pred[2] - pred[1])
}
