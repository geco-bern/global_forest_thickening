# Extract coefficient + CI
extract_coef <- function(fit, lqmm = TRUE) {
  
  if (lqmm) {
    cf <- summary(fit)$tTable

    tibble(
      term = rownames(cf),
      estimate = cf[, "Value"],
      se = cf[, "Std. Error"],
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se
    ) |>
      filter(term == "year_sc")
  } else {
    beta <- fixef(fit)
    vc <- vcov(fit)

    idx <- which(names(beta) == "year_sc")

    est <- beta[idx]
    se <- sqrt(vc[idx, idx])

    tibble(
      estimate = est,
      lower = est - 1.96 * se,
      upper = est + 1.96 * se
    ) |>
      mutate(term = "year_sc")
  }
}
