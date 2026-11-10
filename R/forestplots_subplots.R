# Provider-defined T1 subplots, independent of X/Y availability or orientation.
# Only unresolved tree records are excluded, not their entire parent plot.
forestplots_t1_subplots <- function(data) {
  stopifnot(all(c("plotID", "subplot", "plotsize", "year") %in% names(data)))
  # Provider clarification: .r denotes root measurements within the same T1
  # subplot, not a separate spatial unit. Blank and zero labels remain invalid.
  label <- sub("\\.r$", "", trimws(as.character(data$subplot)))
  number <- suppressWarnings(as.integer(label))
  # SPD-01 expanded in 2010 from 25 to 38 standard 20 x 20 m subplots.
  limit <- ifelse(data$plotID == "SPD-01",
                  ifelse(is.finite(data$year), ifelse(data$year >= 2010, 38L, 25L), NA_integer_),
                  25L)
  resolved <- !is.na(label) & grepl("^[0-9]{1,2}$", label) &
    !is.na(number) & !is.na(limit) & number >= 1L & number <= limit
  parent_area <- dplyr::case_when(
    data$plotID == "SPD-01" ~ limit * 0.04,
    data$plotID %in% c("SPD-02", "ESP-01", "MLA-02") ~ 1,
    TRUE ~ as.numeric(data$plotsize)
  )
  audit <- data |>
    dplyr::mutate(subplot = as.character(subplot), subplot_clean = label,
                  max_subplot = limit, plotsize_major = parent_area, retained = resolved) |>
    dplyr::count(plotID, year, subplot, subplot_clean, max_subplot, plotsize_major,
                 retained, name = "n_tree_census")
  trees <- data[resolved, ] |>
    dplyr::mutate(
      parent_plotID = plotID,
      plotsize_major = parent_area[resolved],
      plotsize = dplyr::if_else(plotID == "ANK-03", 40 * 240 / 25 / 10000, 0.04),
      subplot = sprintf("%02d", number[resolved]),
      plotID = paste(plotID, subplot, sep = "_")
    )
  message("ForestPlots T1: excluded ", sum(!resolved), " tree-census records; retained ",
          dplyr::n_distinct(trees$plotID), " subplots in ",
          dplyr::n_distinct(trees$parent_plotID), " parent plots.")
  list(data = trees, audit = audit)
}

# Candidate 100 x 100 m plots: metadata report 1 ha and observed X/Y span
# the 0-100 m parent-plot frame. Other coordinate conventions need review.
forestplots_grid_plots <- c(
  "IVI-02", "ANK-02", paste0("BOB-0", 1:6), "LPG-01", "LPG-02", "TAM-09"
)

# Apply after tree eligibility filters and parent-level metadata/imputation.
# Missing/out-of-bounds coordinates exclude affected CELLS across all censuses.
# Locate affected cells using the same tree's valid records, then the recorded
# source subplot's overlap with the grid. Never guess a source-label ordering.
forestplots_subplots <- function(data, grid_plots = forestplots_grid_plots) {
  stopifnot(all(c("plotID", "treeID", "source_subplot", "X", "Y", "plotsize") %in%
                  names(data)), !anyNA(data$plotID),
            is.numeric(data$X) || all(is.na(data$X)),
            is.numeric(data$Y) || all(is.na(data$Y)))
  data <- dplyr::mutate(data,
    source_subplot = trimws(as.character(source_subplot)),
    source_subplot = dplyr::na_if(source_subplot, ""),
    grid_x = cut(as.numeric(X), seq(0, 100, 20), include.lowest = TRUE, labels = FALSE),
    grid_y = cut(as.numeric(Y), seq(0, 100, 20), include.lowest = TRUE, labels = FALSE)
  )
  # Stable cell indices, including both 0 and 100 m outer boundaries.
  cell <- ifelse(is.na(data$grid_x) | is.na(data$grid_y), NA_integer_,
                 data$grid_x + 5L * (data$grid_y - 1L))
  keep <- rep(FALSE, nrow(data))
  audits <- list()
  for (parent in unique(data$plotID)) {
    rows <- which(data$plotID == parent)
    d <- data[rows, ]
    cells <- cell[rows]
    missing <- !is.finite(d$X) | !is.finite(d$Y)
    outside <- (is.finite(d$X) & (d$X < 0 | d$X > 100)) |
      (is.finite(d$Y) & (d$Y < 0 | d$Y > 100))
    bad_area <- !is.finite(d$plotsize) | abs(d$plotsize - 1) > 1e-8
    if (!parent %in% grid_plots) {
      # No defensible grid/area: report each recorded source subplot separately.
      audits[[parent]] <- d |>
        dplyr::mutate(n_missing_xy = as.integer(missing),
                      n_outside_grid = as.integer(outside)) |>
        dplyr::group_by(source_subplot) |>
        dplyr::summarise(n_tree_census = dplyr::n(),
          n_missing_xy = sum(n_missing_xy), n_outside_grid = sum(n_outside_grid),
          .groups = "drop") |>
        dplyr::mutate(parent_plotID = parent,
          subplotID = paste0("source_", dplyr::coalesce(source_subplot, "unknown")),
          subplot_basis = "source label; grid unresolved", retained = FALSE,
          reason = "Unresolved coordinate system or subplot geometry")
      next
    }
    # An invalid record can affect several grid cells if its source subplot
    # straddles grid boundaries. Exclude each candidate, not unrelated cells.
    n_missing <- n_outside <- n_area <- n_unresolved <- integer(25)
    for (i in which(missing | outside | bad_area)) {
      candidates <- cells[i]
      if (is.na(candidates)) {
        candidates <- if (!is.na(d$treeID[i])) unique(cells[
          !is.na(d$treeID) & d$treeID == d$treeID[i] & !is.na(cells)
        ]) else integer()
        if (!length(candidates) && !is.na(d$source_subplot[i])) {
          candidates <- unique(cells[!is.na(d$source_subplot) &
            d$source_subplot == d$source_subplot[i] & !is.na(cells)])
        }
      }
      unresolved <- !length(candidates) || all(is.na(candidates))
      if (unresolved) candidates <- seq_len(25)
      n_missing[candidates] <- n_missing[candidates] + as.integer(missing[i])
      n_outside[candidates] <- n_outside[candidates] + as.integer(outside[i])
      n_area[candidates] <- n_area[candidates] + as.integer(bad_area[i])
      n_unresolved[candidates] <- n_unresolved[candidates] + as.integer(unresolved)
    }
    audit <- tibble::tibble(
      parent_plotID = parent, grid_cell = seq_len(25),
      subplotID = sprintf("20m_x%02d_y%02d", (seq_len(25) - 1) %% 5 + 1,
                          (seq_len(25) - 1) %/% 5 + 1),
      subplot_basis = "20 m coordinate grid", source_subplot = NA_character_,
      n_tree_census = tabulate(cells, nbins = 25),
      n_missing_xy = n_missing, n_outside_grid = n_outside,
      n_bad_area = n_area, n_unresolved_membership = n_unresolved
    ) |>
      dplyr::mutate(
        retained = n_missing_xy == 0 & n_outside_grid == 0 & n_bad_area == 0 &
          n_unresolved_membership == 0 & n_tree_census > 0,
        reason = dplyr::case_when(
          n_unresolved_membership > 0 ~ "Unlocated tree cannot be attributed to a subset of cells",
          n_bad_area > 0 ~ "Area missing or incompatible with specified grid",
          n_missing_xy > 0 ~ "Missing coordinates in this cell or overlapping source subplot",
          n_outside_grid > 0 ~ "Out-of-bounds coordinates in this cell or overlapping source subplot",
          n_tree_census == 0 ~ "No located tree records",
          TRUE ~ "Retained"
        )
      )
    keep[rows] <- cells %in% audit$grid_cell[audit$retained]
    audits[[parent]] <- audit
  }
  audit <- dplyr::bind_rows(audits)
  trees <- data[keep, ] |>
    dplyr::mutate(parent_plotID = plotID,
      subplotID = sprintf("20m_x%02d_y%02d", grid_x, grid_y),
      plotID = paste(parent_plotID, subplotID, sep = "_"), plotsize = 0.04)
  message("ForestPlots grid: retained ", dplyr::n_distinct(trees$plotID),
          " cells in ", dplyr::n_distinct(trees$parent_plotID), " parent plots.")
  list(data = trees, audit = audit)
}
