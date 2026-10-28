# Run from the project root: Rscript tests/test_swedish_management.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Exercise the actual stand-preparation expression without running file I/O,
# spatial extraction, or other constituent-dataset processing.
raw <- parse("analysis/00_raw_data.R")
assignments <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("nfi_sweeden"))
}, as.list(raw))
stopifnot(length(assignments) == 2L)
prepare_swedish_nfi <- function(stand) {
  env <- new.env(parent = globalenv())
  env$nfi_sweeden <- stand
  eval(assignments[[2]], envir = env)
  env$nfi_sweeden
}

stand <- tibble(
  plotcode = c(rep("SW1", 3), rep(paste0("SW", 2:7), each = 2)),
  year = c(2010, 2000, 2005, 2000, 2005, 2000, 2005,
           2005, 2010, 2000, 2005, 2010, 2000, 2000, 2005),
  stand_age = c(50, 40, 45, 60, 65, 70, 75, 30, 35, NA, 70, 10, 0, 80, 85),
  management = c(0, 0, 0, 0, 1, 0, NA, 0, 0, 0, 0, 0, 0, 1, 0),
  basal_area = 20, density = 100, lon = 15, lat = 60
)
prepared <- prepare_swedish_nfi(stand)
stopifnot(
  setequal(prepared$plotID, c(1, 4, 5, 6)),
  nrow(prepared) == 9L,
  all(prepared$management == 0),
  all(prepared$management_cat == 1),
  identical(prepared$management_since_census1_yrs,
            c(40, 40, 40, 30, 30, NA, NA, 0, 0)),
  identical(prepared$census, c(1L, 2L, 3L, 1L, 2L, 1L, 2L, 1L, 2L))
)
# First census is determined by year, independently of input row order.
stopifnot(identical(prepared, prepare_swedish_nfi(stand[nrow(stand):1, ])))
cat("Swedish NFI whole-plot management and first-census age checks passed.\n")
