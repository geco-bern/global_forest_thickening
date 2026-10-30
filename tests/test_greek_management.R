# Run from the project root: Rscript tests/test_greek_management.R
suppressPackageStartupMessages(library(dplyr))
source("R/functions.R")
source("R/extract_plot_environment.R")

# Exercise the production expressions without file I/O or other datasets.
raw <- as.list(parse("analysis/00_raw_data.R"))
assignments <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("fep_gre_stand"))
}, raw)
finalise <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("df_fep_gre"))
}, raw)
stopifnot(length(assignments) == 2L, length(finalise) == 2L)
prepare <- function(stand) {
  env <- new.env(parent = globalenv())
  env$fep_gre_stand <- stand
  env$dom_species <- stand |>
    transmute(plotID = plot_id, year = census_year, species = "test") |> distinct()
  suppressMessages(eval(assignments[[2]], env))
  for (expr in finalise) eval(expr, env)
  env$df_fep_gre
}

ids <- c("control", "source_typo", "numeric", "thinning", "mixed",
         "unknown", "missing_baseline", "numeric_decimal")
stand <- tibble(
  plot_id = rep(ids, each = 2),
  census_year = rep(c(2000, 2005), length(ids)),
  info_management = c(
    rep("Control plot (no thinning)", 2), rep(" Contol plot (no thinning) ", 2),
    "34", "39", rep("Light thinning (4 years)", 2),
    "Control plot (no thinning)", "Heavy thinning (4 years)",
    "34", NA, "33", "33", " 7.5 ", "12.5"
  ),
  years_since_management = c(100, 105, 100, 105, 34, 39, 50, 55,
                            20, 25, 34, 39, NA, 38, 7.5, 12.5),
  Longtidute = 22, Latitude = 40, Basal_area_m2_ha.1 = 20,
  N_ha = 100, dbh_m = 0.3, plot_size = 0.25, biomass.t.ha.1 = 100
)
result <- prepare(stand)
expected <- tibble(
  plotID = ids, management = c(0L, 0L, 0L, 1L, 1L, NA, 0L, 0L),
  management_since_census1_yrs = c(100, 100, 34, 50, 20, 34, NA, 7.5)
) |> arrange(plotID)
observed <- result |> distinct(plotID, management, management_since_census1_yrs)
stopifnot(
  isTRUE(all.equal(observed, expected)),
  all(result$management_cat == 1),
  all(result$dataset == "fep_gre_stand"),
  identical(result, prepare(stand[nrow(stand):1, ]))
)
cat("Greek management labels, whole-plot status, and baseline intervals passed.\n")
