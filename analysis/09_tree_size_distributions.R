# Analyse tree-size distribution changes in the unmanaged plots retained for 05.
# Compare each source's native population with DBH >=12.7 and >=20 cm subsets.
# Run from the project root; derived summaries are saved in analysis/diagnostics.

# Load packages ----------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(ggplot2)
library(here)

# Load functions ---------------------------------------------------------------
source(here("R/tree_size_distribution.R"))
source(here("R/tree_density_contributions.R"))
source(here("R/harmonise_census_year.R"))
source(here("R/forestplots_subplots.R"))
source(here("R/forestplots_management.R"))
source(here("R/fit_model.R"))

# Settings ---------------------------------------------------------------------
raw_data_dir <- path.expand(
  "~/data/archive_restricted/GFDYglobe_marques_2025/GFDYglobe_data_raw"
)
out_dir <- here("analysis/diagnostics/tree_size_distributions")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse cached summaries unless the raw inputs or their preparation have changed.
rebuild_tree_metrics <- FALSE
cutoffs <- c(0, 12.7, 20)  # 12.7 cm is FIA's 5-inch subplot threshold.
min_stems <- 10L
min_years <- 5
nboot <- 499L             # Plot resampling for descriptive within-source intervals.
set.seed(1982)

# Read the selected stand sample -----------------------------------------------
stand <- read_rds(here("data/inputs/df_unm_withfilters_v2026.rds")) |>
  ungroup() |>
  select(
    dataset, biome_major, plotID, year, density, QMD, plotsize,
    management_cat, management_since_census1_yrs
  )

stopifnot(!anyDuplicated(select(stand, plotID, year)))

# Read and summarise the tree records ------------------------------------------
## Reuse the source-specific preparation from 01 -------------------------------
# Section markers identify the reviewed tree preparation before aggregation.
# ForestPlots contributes only the retained country extracts listed here.
sections <- tribble(
  ~dataset, ~marker, ~variable,
  "uholka", "# uholka ----", "uholka",
  "inrae_lessem", "# INRAE LESSEM ----", "inrae_lessem",
  "iberbas", "## iberbas ----", "iberbas",
  "unitbv", "## unitbv ----", "unitbv",
  "lwf", "## lwf ----", "lwf_tree",
  "nbw", "## npvbw ----", "npvbw",
  "nwfva", "## nwfva ----", "nwfva_tree",
  "tuzvo", "## tuzvo ----", "tuzvo_tree",
  "ul", "## ul ----", "ul_tree",
  "unito", "## unito ----", "unito",
  "wuls", "## wuls ----", "wuls",
  "luquillo", "## Luquillo ----", "luquillo",
  "scbi", "## SCBI ----", "scbi",
  "serc", "## SERC ----", "serc",
  "wytham", "## Wytham", "wytham",
  "pasoh", "## Pasho ----", "pasoh",
  "forestplots", "## Gabon ----", "gabon",
  "forestplots", "## Ghana ----", "ghana",
  "forestplots", "## Malaysia ----", "malaysia"
)

# Cache only derived metrics, keeping the large tree tables out of later runs.
cache_summary <- function(name, reader, dataset) {
  path <- file.path(cache_dir, paste0(name, ".rds"))

  if (!file.exists(path) || rebuild_tree_metrics) {
    message("Reading trees: ", name)
    trees <- reader()
    result <- summarise_selected_trees(trees, dataset, stand, cutoffs)
    write_rds(result, path)
  } else {
    message("Using metric cache: ", name)
  }

  read_rds(path)
}

results <- list()
failures <- list()

for (i in seq_len(nrow(sections))) {
  row <- slice(sections, i)

  results[[row$variable]] <- tryCatch(
    cache_summary(
      name = row$variable,
      reader = function() {
        trees <- read_prepared_tree_section(
          marker = row$marker,
          variable = row$variable,
          raw_data_dir = raw_data_dir,
          forestplots_metadata = row$dataset == "forestplots"
        )

        # Earlier Wytham preparation excluded explicitly dead-coded stems.
        # Retain ordinary blank codes, including the uncoded 2008 census.
        if (row$dataset == "wytham") {
          trees <- trees |>
            filter(is.na(status) | status != "D")
        }

        trees
      },
      dataset = row$dataset
    ),
    error = function(e) {
      failures[[row$variable]] <<- tibble(
        dataset = row$dataset,
        source = row$variable,
        error = conditionMessage(e)
      )

      message("FAILED ", row$variable, ": ", conditionMessage(e))
      tibble()
    }
  )
}

## Sources supplied separately from their stand summaries ----------------------
results$czu <- cache_summary(
  name = "czu",
  reader = function() {
    read_csv(
      file.path(raw_data_dir, "czu/euf_czu_tree.csv"),
      show_col_types = FALSE
    ) |>
      filter(status == "alive") |>
      mutate(
        dbh = dbh_mm / 10,
        plotsize = plot_size_m2 / 10000
      )
  },
  dataset = "czu"
)

results$fep_gre_stand <- cache_summary(
  name = "fep_gre_stand",
  reader = function() {
    read_csv(
      file.path(raw_data_dir, "fp_gre/fp_gre_tree.csv"),
      show_col_types = FALSE
    ) |>
      filter(status == "Alive") |>
      rename(plotID = plot_id, year = census_year) |>
      left_join(
        stand |>
          filter(dataset == "fep_gre_stand") |>
          mutate(plotID = str_remove(plotID, "^fep_gre_stand_")) |>
          select(plotID, year, plotsize),
        by = c("plotID", "year")
      )
  },
  dataset = "fep_gre_stand"
)

results$nfi_switzerland <- cache_summary(
  name = "nfi_switzerland",
  reader = function() {

    # fread limits the cost of loading large inventory files. Subsequent data
    # selection and transformation use the same tidyverse workflow as above.
    data.table::fread(file.path(raw_data_dir, "nfi_swi/lfi_tree_census.csv")) |>
      as_tibble() |>
      filter(TREE_STATUS == 1) |>
      transmute(
        plotID = PLOTID,
        year = as.integer(str_sub(CENSUS_DATE, 7, 10)),
        dbh = DBH,
        weight = REPRESENTATION
      )
  },
  dataset = "nfi_switzerland"
)

## FIA inventory: process one state at a time ----------------------------------
# Supplied expansion factors account for the different sampling areas used for
# small and large trees. Equal record weights would distort the distribution.
fia_states <- read_csv(
  file.path(raw_data_dir, "fia_us/obs/states.csv"),
  show_col_types = FALSE
)
fia_states <- fia_states |>
  distinct(across(all_of(names(fia_states)[2])), .keep_all = TRUE)

needed_states <- stand |>
  filter(dataset == "fia_us") |>
  pull(plotID) |>
  str_remove("^fia_us_") |>
  str_split("_") |>
  map_int(function(parts) as.integer(parts[2])) |>
  unique()

for (i in which(fia_states[[2]] %in% needed_states)) {
  abbreviation <- fia_states[[3]][i]
  name <- paste0("fia_", abbreviation)

  results[[name]] <- cache_summary(
    name = name,
    reader = function() {
      data.table::fread(
        file.path(raw_data_dir, "fia_us/obs", paste0(abbreviation, "_TREE.csv")),
        select = c(
          "UNITCD", "STATECD", "COUNTYCD", "PLOT", "INVYR", "DIA",
          "STATUSCD", "TPA_UNADJ"
        ),
        showProgress = FALSE
      ) |>
        as_tibble() |>
        filter(STATUSCD == 1) |>
        transmute(
          plotID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = "_"),
          year = INVYR,
          dbh = DIA * 2.54,
          weight = TPA_UNADJ / 0.405
        )
    },
    dataset = "fia_us"
  )
}

write_csv(bind_rows(failures), file.path(out_dir, "read_failures.csv"))

if (length(failures) > 0L) {
  stop("Resolve the recorded source-reading failures before interpreting results.")
}

metrics <- bind_rows(results)

if (anyDuplicated(select(metrics, plotID, year, cutoff_cm))) {
  stop("Duplicate source plot-census metrics.")
}

# Audit tree summaries against the original stand totals -----------------------
# Agreement checks protect against mismatched joins, units, status filters and
# expansion factors. Both native density and QMD must agree within 2%.
native <- metrics |>
  filter(cutoff_cm == 0) |>
  select(
    plotID, year,
    native_N = N,
    native_QMD = QMD_tree,
    native_missing = missing_fraction
  )

censuses <- metrics |>
  left_join(stand, by = c("dataset", "plotID", "year")) |>
  left_join(native, by = c("plotID", "year")) |>
  mutate(
    native_agreement = is.finite(native_N) & is.finite(native_QMD) &
      abs(native_N / density - 1) <= 0.02 & abs(native_QMD / QMD - 1) <= 0.02,
    native_agreement = coalesce(native_agreement, FALSE)
  )

write_csv(censuses, file.path(out_dir, "plot_census_metrics.csv"))

coverage <- stand |>
  count(dataset, name = "selected_censuses") |>
  left_join(
    censuses |>
      filter(cutoff_cm == 0) |>
      group_by(dataset) |>
      summarise(
        tree_censuses = n(),
        plots = n_distinct(plotID),
        agreement_censuses = sum(native_agreement),
        censuses_missing_over_1pct = sum(missing_fraction > 0.01),
        .groups = "drop"
      ),
    by = "dataset"
  ) |>
  mutate(
    availability = case_when(
      dataset == "nfi_sweeden" ~ paste(
        "Tree extract has basal area per ha but no DBH or per-tree expansion factors"
      ),
      is.na(tree_censuses) ~ "No usable tree diameters in supplied inputs",
      .default = "Tree-level analysis"
    )
  )

write_csv(coverage, file.path(out_dir, "dataset_coverage.csv"))

# Load biome STL slopes and pair the censuses ----------------------------------
# Undo predictor scaling before using the saved 05 coefficients. Cache only the
# small coefficient vectors, rather than the large collection of plots and fits.
coef_path <- file.path(cache_dir, "biome_coefs.rds")
analysis_path <- here("data/df_biome_analysis.rds")

if (!file.exists(coef_path) || file.mtime(coef_path) < file.mtime(analysis_path)) {
  message("Loading biome STL slopes from 05")
  analyses <- read_rds(analysis_path)

  biome_coefs <- analyses$out |>
    map(function(out) coef_original_units(out$fit_lqmm[[1]])) |>
    set_names(analyses$biome_major)

  write_rds(biome_coefs, coef_path)
  rm(analyses)
  invisible(gc())
} else {
  biome_coefs <- read_rds(coef_path)
}

beta_qmd <- map_dbl(biome_coefs, function(beta) unname(beta["logQMD"]))
pairs <- make_size_pairs(censuses, beta_qmd, min_stems, min_years)
write_csv(pairs, file.path(out_dir, "paired_plot_changes.csv"))

# Test small-versus-large contributions to density increases --------------------
# The helper also saves the per-source accounting, sensitivity tests and figure.
# This step can be rerun independently with 09b_tree_density_contributions.R.
density_results <- run_density_contribution_analysis(pairs, out_dir)

# Summarise paired changes in distribution metrics -----------------------------
# One value per plot: endpoint change divided by elapsed years, per decade.
# Intervals resample plots within a source, not individual tree records.
change_metrics <- c(
  "d10", "d50", "d90", "cv_dbh", "gini_dbh", "ba_top10",
  "ba_large_share", "n_small", "n_large", "d10_qmd", "d50_qmd", "d90_qmd"
)

summary <- pairs |>
  group_by(dataset, cutoff_cm) |>
  group_modify(function(df, key) {
    metric_names <- c(
      "stl_shift_log_decade", "qmd_change_pct_decade",
      paste0(change_metrics, "_change_decade")
    )

    # Each resample contains the full vector of changes for every sampled plot.
    boot <- replicate(nboot, {
      indices <- sample.int(nrow(df), nrow(df), replace = TRUE)

      df |>
        slice(indices) |>
        select(all_of(metric_names)) |>
        map_dbl(median, na.rm = TRUE)
    })

    tibble(
      metric = metric_names,
      n_plots = nrow(df),
      median = map_dbl(select(df, all_of(metric_names)), median, na.rm = TRUE),
      q025 = apply(boot, 1, quantile, 0.025, na.rm = TRUE),
      q975 = apply(boot, 1, quantile, 0.975, na.rm = TRUE)
    )
  }) |>
  ungroup()

write_csv(summary, file.path(out_dir, "paired_change_summary.csv"))

# Relate distribution changes to the density–QMD shift --------------------------
# These within-source associations are descriptive: counts and diameters enter
# both the shift index and the distribution metrics, creating shared terms.
association_metrics <- c(
  "cv_dbh", "gini_dbh", "ba_top10", "n_small", "n_large", "d10_qmd", "d90_qmd"
)

associations <- pairs |>
  group_by(dataset, cutoff_cm) |>
  group_modify(function(df, key) {
    map(association_metrics, function(metric) {
      complete <- df |>
        transmute(
          y = .data[[paste0(metric, "_change_decade")]],
          x = stl_shift_log_decade
        ) |>
        filter(is.finite(y), is.finite(x))

      rho <- NA_real_

      if (nrow(complete) >= 10L && sd(complete$x) > 0 && sd(complete$y) > 0) {
        rho <- cor(complete$x, complete$y, method = "spearman")
      }

      tibble(metric = metric, n_plots = nrow(complete), rho = rho)
    }) |>
      bind_rows()
  }) |>
  ungroup()

write_csv(associations, file.path(out_dir, "shift_shape_associations.csv"))

## Adjust for changes in QMD ---------------------------------------------------
# Regress per-decade metric change on per-decade log-QMD change. The intercept
# estimates mean change at fixed QMD, with HC3 heteroskedasticity-robust SEs.
# This mean relationship is not a replacement for the 90th-quantile model.
adjusted <- pairs |>
  group_by(dataset, cutoff_cm) |>
  group_modify(function(df, key) {
    if (nrow(df) < 10L) {
      return(tibble())
    }

    df <- df |>
      mutate(dlogqmd = log1p(qmd_change_pct_decade / 100))

    map(c("cv_dbh", "gini_dbh", "ba_top10", "d10_qmd", "d90_qmd"), function(metric) {
      complete <- df |>
        transmute(y = .data[[paste0(metric, "_change_decade")]], x = dlogqmd) |>
        drop_na()

      if (nrow(complete) < 10L || sd(complete$x) < 1e-10) {
        return(tibble())
      }

      fit <- lm(y ~ x, data = complete)
      design <- model.matrix(fit)

      # HC3 inflates each residual by its leverage before the sandwich product.
      bread <- solve(crossprod(design))
      score <- design * as.vector(residuals(fit) / (1 - hatvalues(fit)))
      se <- sqrt(diag(bread %*% crossprod(score) %*% bread))[1]
      intercept <- unname(coef(fit)[1])
      critical <- qt(0.975, df.residual(fit))

      tibble(
        metric = metric,
        n_plots = nrow(complete),
        change_at_fixed_qmd = intercept,
        q025 = intercept - critical * se,
        q975 = intercept + critical * se
      )
    }) |>
      bind_rows()
  }) |>
  ungroup()

write_csv(adjusted, file.path(out_dir, "qmd_adjusted_shape_changes.csv"))

# Sensitivity analyses ---------------------------------------------------------
## Lower stem-count minimum and old-growth selection ---------------------------
sensitivity <- bind_rows(
  make_size_pairs(censuses, beta_qmd, 5L, min_years) |>
    mutate(selection = "min_5_stems"),
  pairs |>
    filter(oldgrowth %in% TRUE) |>
    mutate(selection = "oldgrowth_min_10_stems")
) |>
  group_by(selection, dataset, cutoff_cm) |>
  summarise(
    n_plots = n(),
    cv_change = median(cv_dbh_change_decade),
    gini_change = median(gini_dbh_change_decade),
    stl_shift = median(stl_shift_log_decade),
    .groups = "drop"
  )

write_csv(sensitivity, file.path(out_dir, "sensitivity.csv"))

## Diameter cutoffs compared on identical plot histories -----------------------
# Match both the plot and the census endpoints. This separates truncation of
# the diameter distribution from changes in which plot histories qualify.
common_endpoints <- pairs |>
  count(dataset, plotID, year0, year1) |>
  filter(n == length(cutoffs)) |>
  select(-n)

cutoff_sensitivity <- pairs |>
  inner_join(common_endpoints, by = c("dataset", "plotID", "year0", "year1")) |>
  group_by(dataset, cutoff_cm) |>
  summarise(
    n_plots = n(),
    cv_change = median(cv_dbh_change_decade),
    gini_change = median(gini_dbh_change_decade),
    d90_change = median(d90_change_decade),
    .groups = "drop"
  )

write_csv(
  cutoff_sensitivity,
  file.path(out_dir, "cutoff_same_endpoint_sensitivity.csv")
)

# Save input provenance and settings -------------------------------------------
input_files <- c(
  here("data/inputs/df_unm_withfilters_v2026.rds"),
  here("analysis/01_process_data_orig.R"),
  here("R/tree_size_distribution.R")
)

provenance <- tibble(
  path = input_files,
  modified = as.character(file.info(input_files)$mtime),
  md5 = unname(tools::md5sum(input_files))
)

write_csv(provenance, file.path(out_dir, "input_provenance.csv"))

settings <- list(
  cutoffs_cm = cutoffs,
  min_effective_stems = min_stems,
  min_years = min_years,
  nboot = nboot,
  seed = 1982,
  biome_qmd_slopes = beta_qmd,
  date = Sys.time()
)

write_rds(settings, file.path(out_dir, "settings.rds"))

# Figures ----------------------------------------------------------------------
## Changes above the common 12.7 cm cutoff -------------------------------------
labels <- c(
  d10_change_decade = "D10 (cm/decade)",
  d50_change_decade = "Median DBH (cm/decade)",
  d90_change_decade = "D90 (cm/decade)",
  cv_dbh_change_decade = "DBH CV (/decade)",
  gini_dbh_change_decade = "DBH Gini (/decade)",
  ba_top10_change_decade = "Top 10% BA share (/decade)"
)

gg_changes <- summary |>
  filter(cutoff_cm == 12.7, metric %in% names(labels), n_plots >= 10) |>
  mutate(metric = factor(metric, levels = names(labels), labels = labels)) |>
  ggplot(aes(median, reorder(dataset, median), xmin = q025, xmax = q975)) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = 2) +
  geom_errorbar(orientation = "y", width = 0.25, colour = "#527b70") +
  geom_point(size = 1.8, colour = "#194e40") +
  facet_wrap(~metric, scales = "free_x", ncol = 3) +
  labs(
    x = "Median within-plot change per decade (plot-bootstrap 95% interval)",
    y = NULL,
    title = "How tree diameter distributions change in repeatedly measured plots",
    subtitle = "DBH ≥12.7 cm; ≥10 measured stems at both endpoints; interval ≥5 years",
    caption = paste(
      "Equal plot weight within each dataset. ForestGEO quadrats are subsamples",
      "of single sites; intervals are not among-site uncertainty."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  filename = file.path(out_dir, "diameter_distribution_changes.png"),
  plot = gg_changes,
  width = 13,
  height = 9,
  dpi = 160
)

## Changes across each source's full recorded population -----------------------
gg_native <- summary |>
  filter(cutoff_cm == 0, metric %in% names(labels), n_plots >= 10) |>
  mutate(metric = factor(metric, levels = names(labels), labels = labels)) |>
  ggplot(aes(median, reorder(dataset, median), xmin = q025, xmax = q975)) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = 2) +
  geom_errorbar(orientation = "y", width = 0.25, colour = "#527b70") +
  geom_point(size = 1.8, colour = "#194e40") +
  facet_wrap(~metric, scales = "free_x", ncol = 3) +
  labs(
    x = "Median within-plot change per decade (plot-bootstrap 95% interval)",
    y = NULL,
    title = "Changes across each dataset's full recorded diameter distribution",
    subtitle = "Native inventory thresholds; ≥10 effective stems per endpoint; interval ≥5 years",
    caption = paste(
      "Diameter inclusion thresholds differ among sources.",
      "Site quadrats do not represent independent forest sites."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  filename = file.path(out_dir, "native_diameter_distribution_changes.png"),
  plot = gg_native,
  width = 13,
  height = 9,
  dpi = 160
)

## Sensitivity to the lower diameter cutoff ------------------------------------
gg_cutoff <- cutoff_sensitivity |>
  filter(n_plots >= 10) |>
  pivot_longer(c(cv_change, gini_change), names_to = "metric", values_to = "change") |>
  mutate(
    cut = factor(
      cutoff_cm,
      levels = c(0, 12.7, 20),
      labels = c("All recorded stems", "DBH ≥12.7 cm", "DBH ≥20 cm")
    ),
    metric = recode(metric, cv_change = "Diameter CV", gini_change = "Diameter Gini")
  ) |>
  ggplot(aes(change, dataset, colour = cut)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(position = position_dodge(width = 0.55), size = 2.4) +
  facet_wrap(~metric, scales = "free_x") +
  scale_colour_manual(values = c("#a64429", "#267464", "#3c64a0")) +
  labs(
    x = "Median within-plot change per decade",
    y = NULL,
    colour = NULL,
    title = "The small-tree component can reverse the apparent change in size inequality",
    subtitle = "Identical plots and census endpoints across all three diameter cutoffs",
    caption = paste(
      "All recorded stems means each source's native inventory population.",
      "Only sources with ≥10 shared plot series are shown."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  filename = file.path(out_dir, "diameter_cutoff_sensitivity.png"),
  plot = gg_cutoff,
  width = 11,
  height = 7,
  dpi = 160
)

## Association between STL displacement and distribution change ----------------
gg_associations <- associations |>
  filter(cutoff_cm == 12.7, n_plots >= 10) |>
  mutate(
    metric = factor(
      metric,
      levels = association_metrics,
      labels = c(
        "DBH CV", "DBH Gini", "Top 10% BA share", "Small-stem density",
        "Large-stem density", "D10 / QMD", "D90 / QMD"
      )
    )
  ) |>
  ggplot(aes(metric, dataset, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(
    low = "#b34e32",
    mid = "white",
    high = "#216859",
    limits = c(-1, 1),
    name = "Spearman ρ"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Which size-distribution changes accompany upward density–QMD shifts?",
    subtitle = "Correlation of paired per-decade changes across plots within each source; DBH ≥12.7 cm",
    caption = paste(
      "Shift = Δlog(N) − biome STL slope × Δlog(QMD). Shared diameter/count",
      "terms cause mathematical coupling; this is not causal attribution."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  filename = file.path(out_dir, "shift_shape_associations.png"),
  plot = gg_associations,
  width = 11,
  height = 8,
  dpi = 160
)

message("Analysis saved in ", out_dir)
