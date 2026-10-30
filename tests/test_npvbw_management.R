# Run from the project root: Rscript tests/test_npvbw_management.R
suppressPackageStartupMessages(library(dplyr))
source("R/functions.R")

# Exercise the production expressions without reading or writing datasets.
expressions <- as.list(parse("analysis/00_raw_data.R"))
assignments <- function(name) Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name(name))
}, expressions)
prepare <- function(trees) {
  env <- new.env(parent = globalenv())
  env$npvbw <- trees
  eval(assignments("npvbw_management")[[1]], env)
  eval(assignments("npvbw")[[2]], env)
  for (expr in assignments("df_npvbw")) suppressMessages(eval(expr, env))
  env$df_npvbw
}

ids <- c("prior", "during", "first", "last", "after", "unknown", "varying")
trees <- tibble(
  plotID = rep(ids, each = 2), census = rep(c(2000, 2010), length(ids)),
  year_last_management = c(1980, 1980, 2005, 2005, 2000, 2000,
                           2010, 2010, 2015, 2015, NA, NA, 1980, 2005),
  status = "alive", longitude = 13, latitude = 49, plot_size = 0.1,
  species = "test", dbh = 20
)
# Management evidence on a dead tree must survive the alive-tree filter.
trees <- bind_rows(trees, trees[1, ] |>
                    mutate(plotID = "dead_evidence"),
                  trees[2, ] |> mutate(plotID = "dead_evidence"),
                  trees[1, ] |> mutate(plotID = "dead_evidence",
                                       status = "dead", year_last_management = 2005))
result <- prepare(trees)
expected <- tibble(
  plotID = c(ids, "dead_evidence"),
  management = c(0L, 1L, 1L, 1L, 0L, NA_integer_, 1L, 1L),
  management_since_census1_yrs = c(20, rep(NA_real_, 7))
) |> arrange(plotID)
observed <- result |>
  distinct(plotID, management, management_since_census1_yrs) |> arrange(plotID)
stopifnot(
  isTRUE(all.equal(observed, expected)),
  all(result$management_cat == 1),
  all(result$dataset == "nbw"),
  nrow(result) == 2L * nrow(expected),
  !any(is.nan(result$management_since_census1_yrs)),
  identical(result, prepare(trees[nrow(trees):1, ]))
)
cat("NPVBW whole-plot management and baseline interval checks passed.\n")
