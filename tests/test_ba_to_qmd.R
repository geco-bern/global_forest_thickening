library(dplyr)
source("R/functions.R")

# Two censuses of a half-hectare plot with known individual diameters (cm).
trees <- tibble(
  country = "Test", plotID = "A", lon = 0, lat = 0, plotsize = 0.5,
  year = rep(c(2000, 2005), each = 2), census = rep(1:2, each = 2),
  species = rep(c("a", "b"), 2), dbh = c(10, 30, 20, 40),
  management_since_census1_yrs = NA_real_, biomass = 1
) |>
  mutate(ba_tree = pi * (dbh / 200)^2)
expected <- c(sqrt(mean(c(10, 30)^2)), sqrt(mean(c(20, 40)^2)))

stand <- trees |>
  group_by(plotID, year, census) |>
  summarise(ba = sum(ba_tree) / 0.5, density = n() / 0.5, .groups = "drop")

# Check against individual-tree RMS diameters, and invariance to area units.
stopifnot(
  isTRUE(all.equal(ba_to_qmd(stand$ba, stand$density), expected)),
  isTRUE(all.equal(ba_to_qmd(stand$ba * 0.5, stand$density * 0.5), expected)),
  isTRUE(all.equal(ba_to_qmd(stand$ba / stand$density, 1), expected)),
  is.na(ba_to_qmd(NA_real_, 100)),
  is.na(ba_to_qmd(10, NA_real_))
)

species <- trees |>
  mutate(ba = ba_tree / plotsize, density = 1 / plotsize)

results <- list(
  from_tree_data(trees),
  from_species_data(species),
  from_stand_tree_data(trees, stand),
  from_stand_data(stand)
)
for (result in results) {
  stopifnot(
    isTRUE(all.equal(result$QMD, expected)),
    !any(c("logQMD", "logDensity") %in% names(result)),
    identical(result$year, c(2000, 2005)),
    identical(result$period, c(NA_real_, 5))
  )
}

# Norway preserves source QMD in cm, including missing values, after sorting
# and standardisation; other datasets continue to use the BA/density conversion.
source("R/extract_plot_environment.R")
nfi_norway <- stand |>
  mutate(qmd = c(18, NA_real_), QMD = -999, logQMD = -999, logDensity = -999) |>
  arrange(desc(year))
norway <- from_stand_data(nfi_norway, qmd_col = "qmd") |>
  standardise_stand_data()
stopifnot(
  identical(norway$QMD, c(18, NA_real_)),
  identical(norway$dataset, rep("nfi_norway", 2)),
  identical(norway$year, c(2000, 2005)),
  identical(norway$period, c(NA_real_, 5)),
  !any(c("logQMD", "logDensity") %in% names(norway)),
  isTRUE(all.equal(from_stand_data(nfi_norway)$QMD, expected)),
  inherits(try(from_stand_data(stand, qmd_col = "absent"), silent = TRUE), "try-error")
)

invisible(parse(file = "analysis/00_raw_data.R"))
cat("BA-to-QMD checks passed for all aggregation routes.\n")
