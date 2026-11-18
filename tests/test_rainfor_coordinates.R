# Run from the project root: Rscript tests/test_rainfor_coordinates.R
suppressPackageStartupMessages(library(dplyr))
source("R/rainfor_coordinates.R")

candidates <- tibble(
  plotID = c("A", "A", "A", "A", "B", "B", "C", "C", "D", "ZERO", "NOU-PP"),
  lon = c(-60.1, -60.12, -60.14, -60.12, NA, -70.1, NA, -80, 200, 0, -52.68),
  lat = c(-3.1, -3.12, -3.14, -3.12, -2.12, -2.1, 10, NA, 0, 0, 4.08),
  coord_resolution_deg = c(0.1, 0.01, 0.01, 0.01, 0.01, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01),
  coord_source = paste0("source", seq_len(11)),
  coord_source_doi = paste0("doi", seq_len(11))
)
result <- combine_rainfor_coordinates(candidates)
stopifnot(
  setequal(result$plotID, c("A", "B", "ZERO", "NOU-PP")),
  !anyDuplicated(result$plotID),
  isTRUE(all.equal(result$lon[result$plotID == "A"], -60.13)),
  isTRUE(all.equal(result$lat[result$plotID == "A"], -3.13)),
  result$coord_n_pairs[result$plotID == "A"] == 2L,
  result$coord_resolution_deg[result$plotID == "A"] == 0.01,
  result$lon[result$plotID == "B"] == -70.1,
  result$coord_resolution_deg[result$plotID == "B"] == 0.1,
  result$lon[result$plotID == "ZERO"] == 0,
  is.na(left_join(tibble(plotID = "NOU-59s"), result, by = "plotID")$lon),
  identical(result, combine_rainfor_coordinates(candidates[nrow(candidates):1, ]))
)
# Repeated identical coordinate pairs must not change their weighting.
stopifnot(isTRUE(all.equal(
  result |> select(-coord_source, -coord_source_doi),
  combine_rainfor_coordinates(bind_rows(candidates, candidates)) |>
    select(-coord_source, -coord_source_doi)
)))

# Optional integration check of the actual script and local archive.
raw_data_dir <- path.expand("~/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw")
if (file.exists(file.path(raw_data_dir, "rainfor/sullivan_2025/README.md"))) {
  source("R/functions.R")
  lines <- readLines("analysis/00_raw_data.R")
  start <- grep("^# Stand-level data of plotID, N, BA, and year", lines)
  end <- grep("^# aggregate data from stand$", lines)
  end <- min(end[end > start])
  env <- new.env(parent = globalenv())
  eval(parse(text = lines[seq(start, end - 1L)]), env)
  original <- read.csv(file.path(raw_data_dir, "rainfor/esquivel_etal/full_final_dataset.csv"))
  stopifnot(
    nrow(env$rainfor_missingcoords) == 0L,
    !("NOU-59s" %in% env$rainfor_plots$plotID),
    nrow(env$rainfor_plots) == sum(original$Plot.Code != "NOU-59s"),
    n_distinct(env$rainfor_plots$plotID) == n_distinct(original$Plot.Code) - 1L,
    !anyDuplicated(env$rainfor_coords$plotID),
    all(is.finite(env$rainfor_plots$lon)),
    all(is.finite(env$rainfor_plots$lat))
  )
  # Every chosen resolution is the finest complete valid source for that plot.
  best <- env$rainfor_coordinate_candidates |>
    filter(is.finite(lon), is.finite(lat), abs(lon)<=180, abs(lat)<=90) |>
    group_by(plotID) |> summarise(best = min(coord_resolution_deg), .groups="drop")
  check <- left_join(env$rainfor_coords, best, by="plotID")
  stopifnot(all(check$coord_resolution_deg == check$best))
  cat("Integration passed:", n_distinct(env$rainfor_plots$plotID),
      "retained plots with coordinates; no joined-row duplication.\n")
}
cat("RAINFOR coordinate resolution, duplicates, incomplete pairs and exact joins passed.\n")
