calc_percent_change <- function(df, fit){

  # Pick a fixed value for logQMD_sc (e.g., mean)
  mean_logQMD_sc <- mean(df$logQMD_sc, na.rm = TRUE)
  sd_year <- sd(df$year, na.rm = TRUE)

  # Predict at year_sc = 0 and year_sc = 1 / sd_year (since one calendar year corresponds to 1/sd_year in scaled units)
  newdata1 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 0, plotID = df$plotID[1])
  newdata2 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 1 / sd_year, plotID = df$plotID[1])

  pred1 <- predict(fit, newdata = newdata1)
  pred2 <- predict(fit, newdata = newdata2)

  delta_logDensity <- pred2 - pred1
  percent_change <- (exp(delta_logDensity) - 1) * 100
  return(percent_change)
}
