# Run from the repository root: Rscript tests/test_deduplicate_parent_plots.R
suppressPackageStartupMessages(library(dplyr))
source("R/deduplicate_parent_plots.R")
row <- function(dataset, id, years, density = 250, ba = 25, area = 1,
                lon = 0, lat = 0) {
  tibble(dataset = dataset, plotID = id, year = years, density = density,
         ba = ba, plotsize = area, lon = lon, lat = lat)
}
# Parent comparison must weight areas, not compare individual subplots or sum
# their per-hectare values. Keep the longer source with all of its census rows.
d <- bind_rows(row("multiplots", "AAA-01_01", 2000:2002, 100, 10, .01),
               row("multiplots", "AAA-01_02", 2000:2002, 300, 30, .03),
               row("rainfor", "BBB-02", 2000:2003))
r <- deduplicate_parent_plots(d)
stopifnot(all(r$data$dataset == "rainfor"), nrow(r$data) == 4,
          r$geographic_candidates$confirmed_duplicate,
          r$geographic_candidates$n_shared_years == 3,
          all(r$geographic_year_comparisons$ba_difference_pct < 1e-12))
# One disagreeing shared census prevents removal, as does a missing value.
d_bad <- d
d_bad$ba[d_bad$dataset == "rainfor" & d_bad$year == 2002] <- 28
stopifnot(nrow(deduplicate_parent_plots(d_bad)$data) == nrow(d_bad))
d_bad$ba[d_bad$dataset == "rainfor" & d_bad$year == 2002] <- NA_real_
stopifnot(nrow(deduplicate_parent_plots(d_bad)$data) == nrow(d_bad))
# No overlapping census years, or distance above the threshold: no match.
d_bad <- d
d_bad$year[d_bad$dataset == "rainfor"] <- 2010:2013
stopifnot(!deduplicate_parent_plots(d_bad)$geographic_candidates$confirmed_duplicate)
d_bad <- d
d_bad$lon[d_bad$dataset == "rainfor"] <- 1
stopifnot(nrow(deduplicate_parent_plots(d_bad)$geographic_candidates) == 0)
# The tolerance is explicit and applies to BOTH metrics.
d_bad <- d
d_bad$density[d_bad$dataset == "rainfor"] <- 254
stopifnot(deduplicate_parent_plots(d_bad)$geographic_candidates$confirmed_duplicate,
          !deduplicate_parent_plots(d_bad, measurement_tolerance_pct = 1)$geographic_candidates$confirmed_duplicate)
# Code matches apply regardless of coordinates; ties favour later last census.
# An unrelated, non-coded identifier is never removed.
d_code <- bind_rows(row("multiplots", "AAA-01_01", 2000:2002),
                    row("multiplots", "AAA-01_02", 2000:2002),
                    row("rainfor", "AAA-01", 2001:2003, lon = 5),
                    row("other", "123", 2000:2001))
r <- deduplicate_parent_plots(d_code)
stopifnot(nrow(r$data) == 5, !any(r$data$dataset == "multiplots"),
          all(r$code_matches$coordinate_conflict[!r$code_matches$retained]))
# Input order must not change decisions.
stopifnot(identical(arrange(r$data, dataset, plotID, year),
                    arrange(deduplicate_parent_plots(d_code[nrow(d_code):1, ])$data,
                            dataset, plotID, year)))
# Empty data should preserve its schema and yield empty audits.
r <- deduplicate_parent_plots(d[FALSE, ])
stopifnot(nrow(r$data) == 0, identical(names(r$data), names(d)),
          nrow(r$geographic_candidates) == 0)
cat("Parent-code selection, weighted comparisons, tolerance, missing data, distance and empty-input checks passed.\n")
