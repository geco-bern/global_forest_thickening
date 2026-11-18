# Run from the project root: Rscript tests/test_multiplots_subplots.R
suppressPackageStartupMessages(library(dplyr))
source("R/functions.R")
source("R/forestplots_subplots.R")
row <- function(parent, label, tree, x, y = 1, year = 2000, area = 1) {
  tibble(plotID = parent, source_subplot = label, treeID = tree, X = x, Y = y,
         year = year, country = "Test", lon = 10, lat = 45, plotsize = area,
         species = "test", dbh = 20, ba_tree = pi * 0.1^2, biomass = NA_real_,
         management_since_census1_yrs = NA_real_)
}
trees <- bind_rows(
  # Missing coordinates exclude this cell in ALL censuses, not its neighbour.
  row("A", "one", 1, 1), row("A", "one", 2, NA, year = 2005),
  row("A", "two", 3, 25), row("A", "two", 3, 25, year = 2005),
  # Out-of-bounds tree: identify cell via its source label.
  row("B", "one", 1, 1), row("B", "one", 2, 101), row("B", "two", 3, 25),
  # An ambiguous source subplot overlaps two cells: both excluded; third kept.
  row("C", "wide", 1, 1), row("C", "wide", 2, 25), row("C", "wide", 3, NA),
  row("C", "other", 4, 45),
  # Same tree's other census gives a narrower assignment than the source label.
  row("D", "wide", 1, 1), row("D", "wide", 1, NA, year = 2005),
  row("D", "wide", 2, 25),
  # No tree/subplot mapping: cannot assert that any grid cell is unaffected.
  row("E", "one", 1, 1), row("E", NA_character_, 2, NA),
  # Area error affects the located cell, not every record in the parent.
  row("F", "one", 1, 1, area = NA_real_), row("F", "two", 2, 25),
  # Unknown geometry: report each source subplot, do not invent cell areas.
  row("G", "one", 1, 1), row("G", "two", 2, 25),
  # Boundaries and parent namespaces.
  row("H", "one", 1, 0, 0), row("H", "one", 2, 20, 20),
  row("H", "two", 3, 20.01, 20.01), row("H", "last", 4, 100, 100),
  row("I", "one", 1, 0, 0)
)
allowed <- setdiff(LETTERS[1:9], "G")
result <- forestplots_subplots(trees, allowed)
stopifnot(
  nrow(result$data) == 11L,
  setequal(result$data$parent_plotID, c("A", "B", "C", "D", "F", "H", "I")),
  all(result$data$plotsize == 0.04),
  all(filter(result$data, parent_plotID == "A")$grid_x == 2),
  nrow(filter(result$data, parent_plotID == "A")) == 2L,
  all(filter(result$data, parent_plotID == "C")$grid_x == 3),
  all(filter(result$data, parent_plotID == "D")$grid_x == 2),
  all(!filter(result$audit, parent_plotID == "E")$retained),
  all(filter(result$audit, parent_plotID == "E")$n_unresolved_membership == 1L),
  nrow(filter(result$audit, parent_plotID == "G")) == 2L,
  identical(filter(result$data, parent_plotID == "H")$subplotID,
    c("20m_x01_y01", "20m_x01_y01", "20m_x02_y02", "20m_x05_y05"))
)
reversed <- forestplots_subplots(trees[nrow(trees):1, ], allowed)
stopifnot(identical(arrange(result$data, plotID, treeID, year),
                    arrange(reversed$data, plotID, treeID, year)))
stand <- from_tree_data(result$data)
h <- filter(stand, plotID == "H_20m_x01_y01")
stopifnot(h$density == 50, abs(h$ba - 2 * pi * 0.1^2 / 0.04) < 1e-10,
          abs(h$QMD - 20) < 1e-10)
empty <- forestplots_subplots(filter(trees, plotID == "E"), allowed)
stopifnot(nrow(empty$data) == 0L, nrow(from_tree_data(empty$data)) == 0L)
cat("ForestPlots: cell-level exclusions across censuses, membership mapping, boundaries, stable IDs and aggregation passed.\n")
