# Estimate the biome-by-biome C sink associated with a shift of the 90th
# density–QMD quantile, using the fitted LQMMs from 05_stl_year.R.
# Run 05_stl_year.R and 04_compile_global_drivers.R first.
#
# Environmental covariates do not enter this analysis. The result describes
# biomass change associated with the upper quantile, rather than a direct
# estimate of the change in average observed forest biomass.

# Load packages ----------------------------------------------------------------
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(here)

# Load functions ---------------------------------------------------------------
source(here("R/fit_model.R"))
source(here("R/biome_csink.R"))

# Settings ---------------------------------------------------------------------
nboot <- 2000L  # Increase for the final uncertainty analysis.
workers <- min(4L, parallelly::availableCores(omit = 1L))
seed <- 1982L
year0 <- 2000
year1 <- 2001
carbon_fraction <- 0.5

# Keep outputs separate from 08. The option allows diagnostic runs to save
# locally when the project's data directory is a symlink to external storage.
output_dir <- getOption(
  "forest_thickening.biome_sink_output",
  here("data", paste0("csink_biome_lqmm_nboot_", nboot))
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the fitted models and their matching samples ----------------------------
analysis_file <- here("data/df_biome_analysis.rds")
grid_file <- here("data/df_grid.rds")
message("Loading saved biome fits from 05_stl_year.R")

df_biome_analysis <- read_rds(analysis_file)
biome_inputs <- prepare_biome_sink_inputs(df_biome_analysis)

# Release the large plots and auxiliary fits once the required inputs are saved.
rm(df_biome_analysis)
invisible(gc())

# fit_lqmm is the 'orig' fit in analyse_group(). Retain its matching sample;
# the old-growth selection in the environmental upscaling is a different fit.
mapped <- prepare_biome_sink_grid(read_rds(grid_file))
forest_area <- mapped$area

if (!setequal(names(biome_inputs), forest_area$biome_major)) {
  stop("The saved biome fits and forest map must cover exactly the same biomes.")
}

print(forest_area)

if (any(forest_area$n_missing_cover > 0)) {
  message(
    "Cells with missing forest cover are excluded; ",
    "coverage is recorded in forest_area.csv."
  )
}

# Estimate the sink from the saved 0.9-quantile fits -----------------------------
df_point <- imap(biome_inputs, function(input, biome) {
  estimate_biome_sink(
    data = input$data,
    forest_area_ha = forest_area$forest_area_ha[match(biome, forest_area$biome_major)],
    fit = input$fit,
    year0 = year0,
    year1 = year1,
    carbon_fraction = carbon_fraction
  ) |>
    mutate(
      biome_major = biome,
      n_observations = nrow(input$data),
      .before = 1
    )
}) |>
  bind_rows()

# Bootstrap the STL and biomass fits within each biome -------------------------
# Resampling both relationships together retains their covariance. Global draws
# sum one replicate across all biomes, retaining failed fits as missing values.
message("Fitting ", nboot, " bootstrap replicates per biome on ", workers, " workers")

df_boot <- bootstrap_biome_sink(
  inputs = biome_inputs,
  area = forest_area,
  nboot = nboot,
  workers = workers,
  seed = seed,
  year0 = year0,
  year1 = year1,
  carbon_fraction = carbon_fraction
)

df_boot_global <- sum_biome_sink_draws(df_boot, forest_area$biome_major)

df_summary <- bind_rows(df_boot, df_boot_global) |>
  summarise_biome_sink() |>
  left_join(
    bind_rows(
      df_point |>
        transmute(biome_major, point = sink_pgc_yr),
      tibble(biome_major = "Global", point = sum(df_point$sink_pgc_yr))
    ),
    by = "biome_major"
  )

print(df_summary)

# Allocate the biome sink to mapped forest area --------------------------------
# The expected rate per forest hectare is constant within a biome. Grid-cell
# totals vary with forest cover fraction multiplied by cell area.
df_grid_sink <- mapped$grid |>
  left_join(
    df_summary |>
      filter(biome_major != "Global"),
    by = "biome_major"
  ) |>
  left_join(
    forest_area |>
      select(biome_major, biome_forest_area_ha = forest_area_ha),
    by = "biome_major"
  ) |>
  mutate(
    rate_mean_mgc_ha_yr = mean / biome_forest_area_ha * 1e9,
    sink_mean_pgc_yr = mean * forest_area_ha / biome_forest_area_ha,
    sink_sd_pgc_yr = sd * forest_area_ha / biome_forest_area_ha
  )

# Save estimates, failures, coverage and provenance -----------------------------
settings <- list(
  nboot = nboot,
  workers = workers,
  seed = seed,
  year0 = year0,
  year1 = year1,
  carbon_fraction = carbon_fraction,
  tau = 0.9,
  sample = "05_stl_year.R: orig",
  resampling_unit = "row within biome (as in 05)",
  qmd_integration = "empirical distribution of each biome bootstrap sample",
  forest_biome_codes = c(1, 2, 4, 5, 6, 12),
  inputs = tibble(
    path = c(analysis_file, grid_file),
    modified = file.info(c(analysis_file, grid_file))$mtime
  ),
  run_time = Sys.time()
)

write_rds(
  list(
    settings = settings,
    point = df_point,
    bootstrap = df_boot,
    global = df_boot_global,
    summary = df_summary,
    forest_area = forest_area,
    grid = df_grid_sink
  ),
  file.path(output_dir, "results.rds")
)

write_csv(df_point, file.path(output_dir, "point_estimates.csv"))
write_csv(df_boot, file.path(output_dir, "bootstrap_by_biome.csv"))
write_csv(df_boot_global, file.path(output_dir, "bootstrap_global.csv"))
write_csv(df_summary, file.path(output_dir, "summary.csv"))
write_csv(forest_area, file.path(output_dir, "forest_area.csv"))

# Figures ----------------------------------------------------------------------
## Distribution of bootstrap sink estimates ------------------------------------
gg_distribution <- bind_rows(df_boot, df_boot_global) |>
  ggplot(aes(sink_pgc_yr)) +
  geom_histogram(bins = 15, fill = "grey75", colour = "white", na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~biome_major, scales = "free", ncol = 2) +
  labs(
    x = expression(paste("C sink associated with quantile shift (PgC ", yr^-1, ")")),
    y = "Bootstrap replicates"
  ) +
  theme_bw()

ggsave(
  filename = file.path(output_dir, "sink_distributions.png"),
  plot = gg_distribution,
  width = 10,
  height = 8,
  dpi = 180
)

## Mean sink per forest hectare ------------------------------------------------
gg_map <- df_grid_sink |>
  ggplot(aes(lon, lat, fill = rate_mean_mgc_ha_yr)) +
  geom_tile(width = 0.5, height = 0.5) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85)) +
  scale_fill_gradient2(
    low = "#b34e32",
    mid = "white",
    high = "#216859",
    midpoint = 0,
    na.value = "transparent"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = expression(paste("MgC ", ha^-1, " ", yr^-1)),
    title = "Biomass change associated with the 90th-quantile shift",
    subtitle = "Bootstrap mean per forest hectare; mapped forest biomes only"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(output_dir, "sink_map.png"),
  plot = gg_map,
  width = 10,
  height = 5,
  dpi = 180
)

message("Results saved in ", output_dir)
