suppressPackageStartupMessages(library(dplyr))
source("R/forestplots_subplots.R")
d <- tibble(
  plotID = c(rep("SPD-01", 7), rep("SPD-02", 4), "ESP-01", "MLA-02", "ANK-03"),
  year = c(2009, 2009, 2010, 2010, 2011, 2011, NA, rep(2011, 7)),
  subplot = c("25", "26", "26", "38.r", "39", "0", "01", "23.r", "23", "", "26", "03", "25", "03"),
  plotsize = 9
)
r <- forestplots_t1_subplots(d)
stopifnot(
  identical(r$data$plotID, c("SPD-01_25", "SPD-01_26", "SPD-01_38", "SPD-02_23",
                            "SPD-02_23", "ESP-01_03", "MLA-02_25", "ANK-03_03")),
  identical(r$data$plotsize_major, c(1, 1.52, 1.52, 1, 1, 1, 1, 9)),
  all(r$data$plotsize[1:7] == .04), r$data$plotsize[8] == .0384,
  sum(r$audit$n_tree_census[!r$audit$retained]) == 6L,
  sum(r$audit$n_tree_census) == nrow(d)
)
stopifnot(nrow(forestplots_t1_subplots(d[FALSE, ])$data) == 0L)
source("R/extract_plot_environment.R")
standard <- standardise_stand_data(mutate(r$data, census = 1L))
stopifnot(identical(standard$plotsize_major, r$data$plotsize_major))
cat("T1 suffix normalization, census-dependent area/limits, invalid labels and ANK-03 checks passed.\n")
