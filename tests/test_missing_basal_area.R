# Run from the project root: Rscript tests/test_missing_basal_area.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
})
source("R/functions.R")

# Imputation applies to individual trees before either species or plot sums.
# Supplied basal areas take precedence, even when DBH differs or is missing.
trees <- tibble(
  country = "Test", plotID = rep(c("A", "B"), each = 4),
  lon = 0, lat = 0, plotsize = 0.5,
  year = rep(c(2000, 2000, 2005, 2005), 2),
  species = rep(c("small", "large"), 4),
  dbh = c(20, 40, 20, NA, NA, NA, 10, NA),
  ba_tree = c(NA, NA, 0.02, NA, 0.03, 0.04, 0.2, 0.1),
  management_since_census1_yrs = NA_real_, biomass = 1
)
original <- trees
result <- from_tree_data(trees)
stopifnot(
  identical(trees, original),
  nrow(result) == 4L,
  all(result$density == 4), # Trees with unknown DBH still count.
  isTRUE(all.equal(result$ba, c(pi / 10, NA, 0.14, 0.6))),
  isTRUE(all.equal(result$QMD[1], sqrt(mean(c(20, 40)^2)))),
  is.na(result$QMD[2]),
  is.na(result$ba_inc[2]),
  result$species[1] == "large"
)

# All-missing and mixed NA/NaN measurements must not become zero basal area.
unknown <- trees[1:2, ] |>
  mutate(dbh = c(NA_real_, NaN), ba_tree = c(NaN, NA_real_))
unknown_result <- from_tree_data(unknown)
stopifnot(
  nrow(unknown_result) == 1L,
  identical(unknown_result$ba, NA_real_),
  is.na(unknown_result$QMD),
  unknown_result$density == 4
)
zero <- trees[1:2, ] |> mutate(ba_tree = 0)
stopifnot(from_tree_data(zero)$ba == 0)

# Species-level arithmetic mean DBH must not be used to invent missing totals.
species_data <- trees |>
  mutate(ba = c(1, NA, 1, 2, 3, 4, NA, NA),
         density = c(2, 2, 2, NA, 2, 2, 2, 2))
species_result <- from_species_data(species_data)
stopifnot(
  nrow(species_result) == 4L,
  isTRUE(all.equal(species_result$ba, c(NA, 3, 7, NA))),
  isTRUE(all.equal(species_result$density, c(4, NA, 4, 4))),
  all(is.na(species_result$QMD[c(1, 2, 4)]))
)

# Exercise the actual raw-data preparation pipelines on synthetic records.
# No raw files, rasters, or saved data products are read or overwritten.
raw <- parse("analysis/00_raw_data.R")
prepare <- function(name, data) {
  assignments <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, as.list(raw))
  env <- new.env(parent = globalenv())
  assign(name, data, envir = env)
  eval(tail(assignments, 1)[[1]], envir = env)
  get(name, envir = env)
}
forestgeo <- tibble(
  Quadrat = "A", Census = 1, Latin = "Tree",
  DBH = c(200, NA, 0, NA), Status = c("alive", "alive", "alive", "dead"),
  Date = "01.01.00", PX = 10, PY = 10
)
for (name in c("luquillo", "bci", "scbi", "palanam")) {
  input <- if (name == "bci") mutate(forestgeo, species = "Tree") else forestgeo
  prepared <- prepare(name, input)
  aggregated <- from_tree_data(prepared)
  stopifnot(
    nrow(prepared) == 2L, # Keep the living tree with missing DBH, not dead/zero trees.
    sum(is.na(prepared$dbh)) == 1L,
    nrow(aggregated) == 1L,
    is.na(aggregated$ba),
    is.na(aggregated$QMD),
    aggregated$density == 50
  )
}
# FIA uses expansion factors rather than counting rows. Test the actual
# aggregation expression nested inside the state-processing loop.
find_fia_aggregation <- function(expr) {
  if (!is.call(expr) && !is.expression(expr)) return(list())
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("currentDF")) &&
      any(grepl("STATUSCD", deparse(expr)))) return(list(expr))
  unlist(lapply(as.list(expr), find_fia_aggregation), recursive = FALSE)
}
fia_expr <- find_fia_aggregation(raw)
stopifnot(length(fia_expr) == 1L)
fia <- new.env(parent = globalenv())
fia$currentDF <- tibble(
  plotID = rep(c("complete", "missing_dbh", "missing_weight"), each = 3),
  INVYR = 2000, census = 1, STATUSCD = rep(c(1, 1, 2), 3),
  DIA = c(10, 20, 999, 10, NA, 999, 10, 20, 999),
  TPA_UNADJ = c(5, 10, 99, 5, 10, 99, 5, NA, 99), DRYBIO_AG = 100
)
fia$dom_species <- tibble(
  plotID = c("complete", "missing_dbh", "missing_weight"),
  INVYR = 2000, species = "Tree"
)
eval(fia_expr[[1]], envir = fia)
stopifnot(
  nrow(fia$currentDF) == 3L,
  isTRUE(all.equal(fia$currentDF$ba[1],
                  (pi * (10 * 0.0254 / 2)^2 * 5 +
                   pi * (20 * 0.0254 / 2)^2 * 10) / 0.405)),
  all(is.na(fia$currentDF$ba[2:3])),
  isTRUE(all.equal(fia$currentDF$density, c(15 / 0.405, 15 / 0.405, NA))),
  all(is.na(fia$currentDF$QMD[2:3]))
)
cat("Missing basal-area imputation, propagation, and raw-data retention checks passed.\n")
