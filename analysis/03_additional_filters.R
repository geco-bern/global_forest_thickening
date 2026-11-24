# Implements additional data filters to retain only self-thinning affected plots
# Workflow overview:
#  
# By forest plot
# - Use unmanaged plots only (read data file created by first step, done by Laura)
# - Remove disturbance-affected plots (simultaneous QMD and N drop)
# - Remove ingrowth-affected plots (simultaneous QMD drop and N increase)
#
# By dataset (major)
# Small dataset should be grouped to new classification `dataset_major`
# - Remove years (at head of time series) with substantially lower QMD
#
# By biome
# - Filter based on biome-specific slope
#     - overall slope by plot
#     - for long plot-level time series (N≥5), remove full plot if segment (either first half or second half) has outlier slope
# 
# Workflow/code design:
# - Each step returns a flat data frame with the filter criterium added as an additional column, but no rows removed.
#  
## Loading ------------------------------------------------------------------
library(tidyverse)
library(here)
library(segmented)
library(lqmm)

source(here("R/harmonise_dataset_names.R"))
source(here("R/identify_disturbed_plots.R"))
source(here("R/identify_ingrowth_plots.R"))
source(here("R/identify_badbins.R"))
source(here("R/filter_stl_slope2.R"))
source(here("R/fit_model.R"))
source(here("R/plot_stl_longplots.R"))
source(here("R/analyse_group.R"))

# df_all contains information of management:
# - management_since_census1_yrs (before called years_last_management)
# - management_cat = 1,2 or 3.
# - management = 0 (indicating no management between census)
#
# Among management_cat:
# 1. Recorded history
#    The time since the last management before the first census is known and recorded.
#    the number of years since the last management before the first census, or
#    the calendar year of the last management intervention.
# 2. Pristine/primary/old-growth/protected
#    The time since the last management before the first census is unknown. The forest was considered pristine, primary, old-growth, protected, or minimally used at the time of the first census.
# 3. Unrecorded history
#    The time since the last management before the first census is unknown. The forest may have a history of substantial past wood harvesting and/or visible signs of forest use before the first census.

# read data created by 02_combine_data_unm.R
df_unm <- readRDS(here::here("data/inputs/df_unm_v2026.rds"))
 
## By forest plot --------------------------------------------------------------

### Identify disturbed plots ---------------------------------------------------
# adds column 'disturbed'
df_unm <- df_unm |>
  identify_disturbed_plots()

# ### Identify ingrowth-affected plots -------------------------------------------
# # adds column 'ingrowth'
# df_unm <- df_unm |>
#   identify_ingrowth_plots()

## By (major) dataset ----------------------------------------------------------

### Identify years (in 10-year bins) based on deviating QMD distribution ------------
# adds column 'badqmdbin'
tmp_badbins <- df_unm |> 
  group_by(dataset_major) |> 
  nest() |> 
  mutate(out = purrr::map2(data, dataset_major, ~identify_badbins(.x, .y))) |> 
  mutate(
    gg = purrr::map(out, "gg"),
    data = purrr::map(out, "df")
  ) |> 
  dplyr::select(-out)

gg_qmdfilter <- cowplot::plot_grid(
  plotlist = tmp_badbins$gg,
  ncol = 4,
  labels = letters[1:18]
)

ggsave(
  here("fig/qmdfilter.pdf"),
  plot = gg_qmdfilter,
  width = 12,
  height = 10
)

df_unm <- tmp_badbins |> 
  dplyr::select(-gg) |> 
  unnest(data) |> 
  ungroup()

## By (major) biome ----------------------------------------------------------
# ### Identify non-self-thinning plots based on outlying slope ----------
# # separately for each (major) biome
# # adds column 'badslope'
# tmp_slopefilter <- df_unm |>
#   group_by(biome_major) |>
#   nest() |>
#   filter(!is.na(biome_major)) |>
#   mutate(
#     mod = purrr::map(data, ~ fit_model(.)),
#     out = purrr::pmap(
#       list(data, mod, biome_major),
#       ~ filter_stl_slope2(..1, ..2, ..3, multiplier_sd = 1.0)
#     )
#   ) |>
#   mutate(
#     gg = purrr::map(out, "gg"),
#     data = purrr::map(out, "df")
#   ) |>
#   dplyr::select(-out)

# gg_slopefilter <- cowplot::plot_grid(
#   plotlist = tmp_slopefilter$gg,
#   ncol = 2,
#   labels = letters[1:18]
# )

# ggsave(
#   here("manuscript/figures/slopefilter.pdf"),
#   plot = gg_slopefilter,
#   width = 12,
#   height = 8
# )

# # For manuscript text:
# # Median slope by biome:
# tmp_slopefilter |>
#   mutate(slope_med = purrr::map_dbl(data, ~median(.$slope_lin, na.rm = TRUE)))

# # Make df flat again
# df_unm_withfilters <- tmp_slopefilter |>
#   dplyr::select(biome_major, data) |>
#   unnest(data) |>
#   unnest(data) |>
#   dplyr::select(-logDensity, -logQMD)

# because code above is commented out
df_unm_withfilters <- df_unm 

write_rds(df_unm_withfilters, here("data/inputs/df_unm_withfilters_v2026.rds"))

# Data coverage tables for SI ---------------------------------------
df_unm_withfilters <- read_rds(here("data/inputs/df_unm_withfilters_v2026.rds")) |> 
  tidyr::drop_na(QMD, density)

coverage_filter_subsets <- coefficient_filter_subsets(df_unm_withfilters)[
  c("no_mgmt_30", "no_mgmt_100", "oldgrowth")
]

summarise_coverage <- function(
  group_column,
  data = df_unm_withfilters,
  subsets = coverage_filter_subsets
) {
  coverage <- data |>
    group_by(.data[[group_column]], plotID) |>
    summarise(
      n_observations = n_distinct(year),
      first_year = min(year),
      last_year = max(year),
      .groups = "drop_last"
    ) |>
    summarise(
      n_observations = sum(n_observations),
      n_plots = n(),
      mean_period = mean(last_year - first_year),
      first_year = min(first_year),
      last_year = max(last_year),
      .groups = "drop"
    ) |>
    arrange(desc(n_observations))

  plot_counts <- imap(subsets, \(subset, stage) {
    subset |>
      group_by(.data[[group_column]]) |>
      summarise("{stage}" := n_distinct(plotID), .groups = "drop")
  })

  reduce(
    plot_counts,
    \(table, counts) left_join(table, counts, by = group_column),
    .init = coverage
  ) |>
    mutate(across(all_of(names(subsets)), \(n) replace_na(n, 0L))) |>
    relocate(all_of(names(subsets)), .after = n_plots) |>
    relocate(mean_period, .after = last_year)
}

mean_plot_areas <- df_unm_withfilters |>
  filter(is.finite(plotsize), plotsize > 0) |>
  group_by(dataset, plotID) |>
  summarise(plot_area_ha = mean(plotsize), .groups = "drop_last") |>
  summarise(mean_plot_area_ha = mean(plot_area_ha), .groups = "drop")

df_coverage_dataset <- summarise_coverage("dataset") |>
  left_join(mean_plot_areas, by = "dataset")

df_coverage_biome <- bind_rows(
  summarise_coverage("biome_major"),
  summarise_coverage(
    "biome_major",
    data = mutate(df_unm_withfilters, biome_major = "Total"),
    subsets = map(coverage_filter_subsets, \(subset) mutate(subset, biome_major = "Total"))
  )
)

write_coverage_table <- function(
  coverage,
  group_label,
  caption,
  label,
  filename,
  bold_last_row = FALSE
) {
  display <- coverage |>
    mutate(
      across(
        all_of(c("n_observations", "n_plots", names(coverage_filter_subsets))),
        \(n) formatC(n, format = "d", big.mark = ",")
      ),
      mean_period = sprintf("%.1f", mean_period),
      across(any_of("mean_plot_area_ha"), \(area) {
        if_else(is.na(area), "--", sprintf("%.3f", area))
      })
    )

  headings <- c(
    group_label,
    "Obs.",
    "Plots",
    "30 yr",
    "100 yr",
    "Old-growth",
    "First",
    "Last",
    "Span (yr)"
  )
  if ("mean_plot_area_ha" %in% names(display)) {
    headings <- c(headings, "Area (ha)")
  }

  tex <- knitr::kable(
    display,
    format = "latex",
    longtable = TRUE,
    booktabs = TRUE,
    align = c("@{}p{3.2cm}", rep("r", ncol(display) - 2L), "r@{}"),
    col.names = headings,
    caption = caption,
    label = label,
    linesep = ""
  )
  if (bold_last_row) {
    tex <- kableExtra::row_spec(tex, nrow(display), bold = TRUE)
  }

  writeLines(
    c(
      "% Generated by analysis/03_additional_filters.R.",
      "\\begingroup",
      "\\sffamily\\small",
      "\\setlength{\\tabcolsep}{3pt}",
      tex,
      "\\endgroup"
    ),
    here("manuscript", "tables", filename)
  )
}

write_coverage_table(
  df_coverage_dataset, "Dataset", "Inventory coverage by constituent dataset.",
  "coverage-dataset", "coverage_by_dataset.tex"
)
write_coverage_table(
  df_coverage_biome, "Major biome", "Inventory coverage by major forest biome.",
  "coverage-biome", "coverage_by_biome.tex", bold_last_row = TRUE
)

# Figures for SI -----------------
## Length of repeated observations ----------------------------------
df_len <- df_unm_withfilters |>
  group_by(plotID, biome_major) |>
  summarise(start = min(year), end = max(year)) |>
  mutate(len = end - start)

# numbers for paper
df_len |>
  ungroup() |>
  group_by(biome_major) |>
  summarise(
    len_median = median(len),
    len_mean = mean(len)
  )

df_len |>
  ungroup() |>
  summarise(
    len_median = median(len),
    len_mean = mean(len)
  )

df_len |>
  ggplot(aes(x = len, color = biome_major, fill = biome_major)) +
  geom_density(adjust = 3, alpha = 0.5) +
  scale_fill_manual(
    values = c(
      "Boreal Forests/Taiga" = "dodgerblue4",
      "Mediterranean Forests" = "orangered3",
      "Temperate Forests" = "darkgreen",
      "Tropical & Subtropical Broadleaf Forests" = "springgreen3"
    ),
    na.value = NA,
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  scale_color_manual(
    values = c(
      "Boreal Forests/Taiga" = "dodgerblue4",
      "Mediterranean Forests" = "orangered3",
      "Temperate Forests" = "darkgreen",
      "Tropical & Subtropical Broadleaf Forests" = "springgreen3"
    ),
    na.value = NA,
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = "Length (years)",
    y = "Density"
  )

ggsave(
  filename = here(
    "manuscript/figures/distribution_length.pdf"
  ),
  width = 10,
  height = 5
)

## Barplot: observations by sources over time ----------------
plot_hist_year_biome <- function(df, biome_name, xmin) {
  df |>
    mutate(dataset = harmonise_dataset_names(dataset)) |>
    filter(biome_major == biome_name) |>
    ggplot(aes(x = year, fill = dataset)) +
    geom_histogram(
      color = "black",
      binwidth = 1,
      position = "stack",
      linewidth = 0.3
    ) +
    theme_bw() +
    coord_cartesian(xlim = c(xmin, 2024)) +
    labs(
      title = biome_name,
      x = "Year",
      y = "Count",
      fill = ""
    ) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    theme(
      legend.position = "right",
      legend.key.height = grid::unit(0.35, "cm"),
      legend.text = element_text(size = 9)
    )
}

gg_hist_year_tropical <- plot_hist_year_biome(
  df_unm_withfilters,
  "Tropical & Subtropical Broadleaf Forests",
  xmin = 1970
)

gg_hist_year_temperate <- plot_hist_year_biome(
  df_unm_withfilters,
  "Temperate Forests",
  xmin = 1960
)

gg_hist_year_boreal <- plot_hist_year_biome(
  df_unm_withfilters,
  "Boreal Forests/Taiga",
  xmin = 1980
)

gg_hist_year_mediterranean <- plot_hist_year_biome(
  df_unm_withfilters,
  "Mediterranean Forests",
  xmin = 1980
)

# Combine
row1 <- cowplot::plot_grid(
  gg_hist_year_tropical,
  ncol = 1,
  labels = letters[1]
)

row2 <- cowplot::plot_grid(
  gg_hist_year_temperate,
  ncol = 1,
  labels = letters[2]
)

row3 <- cowplot::plot_grid(
  gg_hist_year_boreal,
  gg_hist_year_mediterranean,
  ncol = 2,
  labels = letters[3:4]
)

fig_hist_year <- cowplot::plot_grid(
  row1,
  row2,
  row3,
  ncol = 1
)
fig_hist_year

ggsave(
  filename = here("manuscript/figures/fig_hist_year.pdf"),
  plot = fig_hist_year,
  width = 9,
  height = 10
)

## Longest-observation --------------------------------------------------------

# This figure requires the optional slope filter.
if ("badslope" %in% names(df_unm_withfilters)) {

  gg_longplots_tropical <- plot_stl_longplots(
    df_unm_withfilters,
    biome_name = "Tropical & Subtropical Broadleaf Forests"
  )

  gg_longplots_temperate <- plot_stl_longplots(
    df_unm_withfilters,
    biome_name = "Temperate Forests"
  )

  gg_longplots_boreal <- plot_stl_longplots(
    df_unm_withfilters,
    biome_name = "Boreal Forests/Taiga"
  )

  gg_longplots_mediterranean <- plot_stl_longplots(
    df_unm_withfilters,
    biome_name = "Mediterranean Forests"
  )

  fig_stl_longplots <- cowplot::plot_grid(
    gg_longplots_tropical,
    gg_longplots_temperate,
    gg_longplots_boreal,
    gg_longplots_mediterranean,
    ncol = 2,
    labels = letters[1:4]
  )

  fig_stl_longplots

  ggsave(
    filename = here("manuscript/figures/stl_longplots.pdf"),
    plot = fig_stl_longplots,
    width = 9,
    height = 7
  )
}

# Numbers for main text ------
n_observations <- nrow(unique(df_unm_withfilters[c("plotID", "year")]))
n_plots <- length(unique(df_unm_withfilters$plotID))

# with text from the article:
# data size
cat(sprintf(
  paste0(
    "We compiled inventory data from unmanaged forest stands of six major ",
    "forest biomes (Fig. S1, Tab. S1) into a harmonized dataset, ",
    "comprising %s observations from %s plots.\n"
  ),
  formatC(n_observations, format = "d", big.mark = ","),
  formatC(n_plots, format = "d", big.mark = ",")
))

n_obs_availableinfo <- df_unm_withfilters |>
  distinct(plotID, management_cat) |>
  filter(management_cat %in% c(1,2)) |>
  nrow()

n_sites_oldgrowth <- df_unm_withfilters |>
  distinct(plotID, management_cat) |>
  filter(management_cat == 2) |>
  nrow()

n_obs_oldgrowth <- df_unm_withfilters |> 
  distinct(plotID, year, management_cat) |> 
  filter(management_cat == 2) |> 
  nrow()

cat(sprintf(
  paste0(
    "For %s plots, forest management prior to the first ",
    "available census was specified, of which %s were ",
    "classified as old-growth (%s observations in total) with absent signs of wood extraction ..., \n"
  ),
  formatC(n_obs_availableinfo, format = "d", big.mark = ","),
  formatC(n_sites_oldgrowth, format = "d", big.mark = ","),
  formatC(n_obs_oldgrowth, format = "d", big.mark = ",")
))

# years covered
firstyear <- min(df_unm_withfilters$year)
lastyear <- max(df_unm_withfilters$year)

cat(sprintf(
  paste0(
    "We investigated the constancy of the NDSR over the past decades (covering %s-%s, Figs. S2, S3),\n"
  ),
  formatC(firstyear, format = "d"),
  formatC(lastyear, format = "d")
))


