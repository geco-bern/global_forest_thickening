# Run from the project root: Rscript tests/test_australian_management.R
suppressPackageStartupMessages(library(dplyr))
expressions <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("aus_plots"))
}, as.list(parse("analysis/00_raw_data.R")))
management_expression <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("aus_management"))
}, as.list(parse("analysis/00_raw_data.R")))[[1]]
prepare <- function(data) {
  env <- new.env(parent = globalenv())
  env$aus_plots <- data
  eval(management_expression, env)
  eval(expressions[[2]], env)
  env$aus_plots
}
stand <- tibble(
  plotID = rep(c("early", "late", "same_year", "during", "after", "mixed"), each = 2),
  year = c(1950, 1960, rep(c(2000, 2010), 5)),
  year_last_management = c(1940, 1940, 1980, 1980, 2000, 2000,
                           2005, 2005, 2020, 2020, 1980, 2005),
  longitude = 150, latitude = -30, TreesPerHectareAHC1_2 = 100,
  DBHqAHC1_2_cm = 20, BasalAreaAHC1_2_m2perha = 10,
  plot_size = 1, dominant_species = "test", AbovegroundAHC1_2_Mgperha = 50
)
result <- prepare(stand)
stopifnot(
  identical(result$management, rep(c(0L, 0L, 1L, 1L, 1L, 1L), each = 2)),
  identical(result$management_since_census1_yrs, c(10, 10, 20, 20, rep(NA_real_, 8))),
  all(result$management_cat == 1L),
  !is_grouped_df(result),
  nrow(result) == nrow(stand),
  identical(result, prepare(stand[nrow(stand):1, ])[nrow(stand):1, ])
)
missing_dates <- bind_rows(
  stand[1:2, ] |> mutate(plotID = "missing", year_last_management = NA_real_),
  stand[1:2, ] |> mutate(plotID = "partly_missing", year_last_management = c(NA, 1940)),
  stand[1:2, ] |> mutate(plotID = "managed_missing", year_last_management = c(NA, 1955))
)
missing_result <- prepare(missing_dates)
stopifnot(
  identical(missing_result$management, c(0L, 0L, 0L, 0L, 1L, 1L)),
  identical(missing_result$management_since_census1_yrs,
            c(NA_real_, NA_real_, 10, 10, NA_real_, NA_real_)),
  all(missing_result$management_cat == 1L)
)
decimal <- stand[1:2, ] |>
  mutate(year = c(2000.25, 2005.5), year_last_management = 2000)
decimal_result <- prepare(decimal)
stopifnot(all(decimal_result$management == 1L),
          all(is.na(decimal_result$management_since_census1_yrs)),
          all(decimal_result$management_cat == 1L))
earlier_decimal <- prepare(decimal |> mutate(year_last_management = 1998))
stopifnot(all(earlier_decimal$management == 0L),
          all(earlier_decimal$management_since_census1_yrs == 2))
# Conflicting earlier dates use the latest event for every census, including
# rows where the date is absent. All three management fields must be constant.
conflicting <- stand[1:2, ] |>
  mutate(year_last_management = c(1930, 1940))
conflicting_result <- prepare(conflicting)
stopifnot(all(conflicting_result$management_since_census1_yrs == 10))
for (data in list(result, missing_result, decimal_result, earlier_decimal, conflicting_result)) {
  labels <- data |> distinct(plotID, management, management_since_census1_yrs, management_cat)
  stopifnot(nrow(labels) == n_distinct(data$plotID))
}
cat("Australian per-plot first census, whole-plot management, and elapsed-year checks passed.\n")
