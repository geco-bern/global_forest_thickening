# Global C sink ----------------------------------------------------------------
# Calculates the global C sink as a consequence of the STL change, scaled up to
# global forests with environmental covariates. Uncertainty estimation with
# bootstrapping

library(dplyr)
library(readr)
library(lme4)
library(purrr)
library(tidyr)
library(here)
library(rsample)
library(recipes)
library(multidplyr)
library(tictoc)
library(tidyterra)
library(rnaturalearth)
library(ggplot2)
library(scales)

source(here("R/harmonise_dataset_names.R"))

## Load data -------------------------------------------------------------------
df <- read_rds(here("data/inputs/df_unm_withfilters_v2026.rds"))|>

  dplyr::ungroup() |>
  dplyr::filter(is.finite(QMD), QMD > 0, is.finite(density), density > 0) |>
  dplyr::mutate(logDensity = log(density), logQMD = log(QMD)) |>
  # mutate(across(any_of(c("dataset", "dataset_major")), harmonise_dataset_names)) |>

  # variable substitute, needed here
  mutate(NQMD2 = density * QMD^2) |>

    # use only "primary" for global model
    # # use only "old-growth" for global model
    filter(
      management_cat == 2 |
        (management_cat == 1 & management_since_census1_yrs >= 300)
    )

    ### Plot data, filtered upper edge ---------------------------------------------

### Maps of environmental covariates -------------------------------------------
# Load data for upscaling: maps of environmental factors
grid_drivers <- read_rds(here("data/df_grid.rds")) |>
  as_tibble() |>
  mutate(
    lon_i = round(lon * 4),
    lat_i = round(lat * 4)
  )

# Reuse each fit's recipe for both calendar years, including all interactions.
source(here("R/global_csink.R"))

# Negative soil values occur in the cached inputs (including averaged missing
# codes). Report them explicitly; do not silently change the scientific sample.
for (input_name in c("df", "grid_drivers")) {
  soil <- get(input_name)[, c("ORGC", "CNrt", "PBR")]
  if (any(as.matrix(soil) < 0, na.rm = TRUE)) {
    warning(input_name, ": negative soil predictors present. Rebuild soil inputs with ",
            "missing codes masked before aggregation before interpreting the soil model.")
  }
}

## Create bootstraps ------------
set.seed(1982)
n_boot <- 9
boot_resamples <- bootstraps(df, times = n_boot)

## Run once, serially (ncores = 1) or in parallel --------------------------
ncores <- 4 # parallel::detectCores() - 2

if (ncores > 1L) {
  cl <- new_cluster(n = ncores) |>
    cluster_library(
      packages = c(
        "dplyr",
        "tidyr",
        "lme4",
        "purrr",
        "recipes",
        "rsample"
      )
    ) |>
    cluster_assign(
      fit_stl_byboot = fit_stl_byboot,
      fit_biomass_byboot = fit_biomass_byboot,
      predict_dn = predict_dn,
      predict_db = predict_db,
      df = df,
      grid_drivers = grid_drivers,
      calc_csink_global = calc_csink_global,
      sum_csink_available = sum_csink_available
    )
}

tic()
df_boot <- boot_resamples |> mutate(id = row_number())

if (ncores > 1L) df_boot <- partition(df_boot, cl)

df_boot <- df_boot |>
  mutate(
    bundle_stl = map(splits, ~ fit_stl_byboot(analysis(.x))),
    model_biomass = map(splits, ~ fit_biomass_byboot(analysis(.x)))
  ) |>
  select(-splits) |>
  mutate(
    grid_predictions = map2(
      bundle_stl, id,
      ~ withr::with_seed(1982L + .y,
                        predict_dn(.x, vec_qmd = df$QMD, df = grid_drivers))
    )) |>
  select(-bundle_stl) |>
  mutate(
    grid_predictions = map2(
      grid_predictions,
      model_biomass,
      ~ predict_db(.x, .y)
    )
  ) |>
  select(-model_biomass) |>

  # calculate C sink by latitude (retain to keep spread across bootstraps)
  mutate(db_pgc_global = map_dbl(grid_predictions, ~calc_csink_global(., grid_drivers)))

if (ncores > 1L) df_boot <- collect(df_boot)

df_boot <- arrange(df_boot, id)

toc()

write_rds(
  df_boot,
  file = here(paste0("data/df_boot_nboot_", as.character(n_boot),".rds"))
)

# df_boot <- read_rds(here("data/df_boot_nboot_2000.rds"))

### Summarise across bootstraps ------------------------------------------------
# stack predictions from all bootstrap samples into (very) long vector
df_summ <- df_boot |>
  select(-db_pgc_global) |>
  unnest(grid_predictions) |>
  group_by(lat_i, lon_i) |>
  summarise(
    db_mean = mean(db, na.rm = TRUE),
    db_median = median(db, na.rm = TRUE),
    db_sd = sd(db, na.rm = TRUE)
  ) |>
  drop_na() |>
  ungroup() |>
  left_join(
    grid_drivers,
    by = join_by(lon_i, lat_i)
  ) |>
  mutate(lon = lon_i/4, lat = lat_i/4) |>
  mutate(
    across(
      starts_with("db"),
      list(
        gc_per_ha_forest   = ~ .x * 1e3 * 0.5,
        gc_per_ha_gridcell = ~ (.x * 1e3 * 0.5) * fcf
      ),
      .names = "{.col}_{.fn}"
    )
  )

write_rds(
  df_summ,
  here(paste0("data/df_summ_nboot_", as.character(n_boot),".rds"))
)

# df_summ <- read_rds(here("data/df_summ_nboot_300.rds"))

## Visualisations --------------------------------------------------------------
### Distribution of global C sink estimates ------------------------------------
# sum across all gridcells, distribution across bootstraps
gg_hist_csink_boot <- df_boot |>
  ggplot(aes(db_pgc_global, after_stat(density))) +
  geom_histogram(fill = "grey", color = "black", bins = 30) +
  labs(
    x = expression(paste("PgC ", yr^-1)),
    y = "Density"
  ) +
  theme_bw()

ggsave(
  here("manuscript/figures/gg_hist_csink_boot.pdf"),
  plot = gg_hist_csink_boot,
  width = 5,
  height = 4
)

# numbers for paper
mean(df_boot$db_pgc_global)
quantile(df_boot$db_pgc_global, probs = c(0.025, 0.975))

# # across gridcells (median across bootstraps)
# gg_hist_db_boot <- df_summ |>
#   ggplot(aes(db_median_gc_per_ha_forest * 1e-6, after_stat(density))) +
#   geom_histogram(fill = "grey", color = "black", bins = 50) +
#   labs(
#     x = expression(paste("MgC ha"^-1, "yr"^-1)),
#     y = "Density"
#   ) +
#   theme_bw()
#
# ggsave(
#   here("manuscript/figures/gg_hist_db_boot.pdf"),
#   plot = gg_hist_db_boot,
#   width = 5,
#   height = 4
#   )

### Map ---------------------------------------------------------------
#### Per forest area -----------
##### Mean --------------
# df_summ <- read_rds(file = here("data/df_summ.rds"))

coast <- rnaturalearth::ne_coastline(
  scale = 110,
  returnclass = "sf"
)

layer_ocean <- rnaturalearth::ne_download( # ne_load(
  scale = 110,
  type = "ocean",
  category = "physical",
  returnclass = "sf",
  destdir = here("data/")
)

gg_map_sink_perforestarea <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_median_gc_per_ha_forest * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  # scale_fill_stepsn(
  #   colours = rev(khroma::color("berlin")(50)), # viridisLite::viridis(7),
  #   guide = guide_coloursteps(),
  #   na.value = "grey30", # <- missing data color
  #   limits = c(-1.5, 1.5),
  #   oob = squish, # clamp values outside limits,
  #   name = expression(paste("MgC ha"^-1, "yr"^-1)),
  # ) +
  khroma::scale_fill_berlin(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    # limits = c(-1.5, 1.5),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1)),
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey40", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_perforestarea.pdf"),
  plot = gg_map_sink_perforestarea,
  width = 8,
  height = 6
)

##### SD --------------
gg_map_sink_perforestarea_sd <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_sd_gc_per_ha_forest * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  khroma::scale_fill_lajolla(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    limits = c(0, 0.2),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey40", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_perforestarea_sd.pdf"),
  plot = gg_map_sink_perforestarea_sd,
  width = 8,
  height = 4
)

gg_map_sink_perforestarea_combined <- cowplot::plot_grid(
  gg_map_sink_perforestarea,
  gg_map_sink_perforestarea_sd,
  ncol = 1,
  labels = letters[1:2]
)

ggsave(
  here("manuscript/figures/gg_map_sink_perforestarea_combined.pdf"),
  plot = gg_map_sink_perforestarea_combined,
  width = 8.2,
  height = 5.5
)

#### Per grid area -----------
##### Mean --------------
gg_map_sink_pergridarea <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_median_gc_per_ha_gridcell * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  khroma::scale_fill_berlin(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    # limits = c(-1.5, 1.5),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey30", color = NA)) # ocean = light grey

gg_map_sink_pergridarea

ggsave(
  here("manuscript/figures/gg_map_sink_pergridarea.pdf"),
  plot = gg_map_sink_pergridarea,
  width = 8,
  height = 6
)

##### SD --------------
gg_map_sink_pergridgridtarea_sd <- df_summ |>
  ggplot() +
  geom_raster(
    aes(lon, lat, fill = db_sd_gc_per_ha_gridcell * 1e-6),
    show.legend = TRUE
  ) +
  geom_sf(
    data = layer_ocean,
    color = NA,
    fill = "grey60"
  ) +
  geom_sf(
    data = coast,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  khroma::scale_fill_lajolla(
    reverse = TRUE,
    midpoint = 0,
    na.value = "grey30", # <- missing data color
    limits = c(0, 0.15),
    oob = squish, # clamp values outside limits,
    name = expression(paste("MgC ha"^-1, "yr"^-1))
  ) +
  theme_void() +
  theme(panel.background = element_rect(fill = "grey40", color = NA)) # ocean = light grey

ggsave(
  here("manuscript/figures/gg_map_sink_pergridgridtarea_sd.pdf"),
  plot = gg_map_sink_pergridgridtarea_sd,
  width = 8,
  height = 4
)

gg_map_sink_pergridarea_combined <- cowplot::plot_grid(
  gg_map_sink_pergridarea,
  gg_map_sink_pergridgridtarea_sd,
  ncol = 1,
  labels = letters[1:2]
)

ggsave(
  here("manuscript/figures/gg_map_sink_pergridarea_combined.pdf"),
  plot = gg_map_sink_pergridarea_combined,
  width = 8.2,
  height = 5.5
)

### Sink across latitude ---------
# Use the same run as the map and global histogram.
df_bylat <- df_boot |>
  select(-db_pgc_global) |>
  mutate(csink_tgc_bylat = map(grid_predictions, ~calc_csink_bylat(., grid_drivers))) |>
  select(-grid_predictions) |>
  unnest(csink_tgc_bylat) |>
  group_by(lat_i) |>
  summarise(
    db_tgc_bylat_mean = mean(db_tgc, na.rm = TRUE),
    db_tgc_bylat_median = median(db_tgc, na.rm = TRUE),
    db_tgc_bylat_q05 = quantile(db_tgc, probs = 0.05, na.rm = TRUE),
    db_tgc_bylat_q95 = quantile(db_tgc, probs = 0.95, na.rm = TRUE),
  ) |>
  mutate(lat = lat_i / 4)

# number for paper
csink_tropical <- df_bylat |>
  mutate(is_topical = ifelse(lat > -30 & lat < 30, TRUE, FALSE)) |>
  group_by(is_topical) |>
  summarise(db_tgc = sum(db_tgc_bylat_mean, na.rm = TRUE))

csink_tropical

gg_csink_bylat <- df_bylat |>
  ggplot() +   # multiplication by two because 0.5 degree resolution
  geom_ribbon(
    aes(
      x = lat,
      ymin = db_tgc_bylat_q05 * 2,
      ymax = db_tgc_bylat_q95 * 2),
    fill = "grey50",
    color = NA) +
  # geom_line(aes(lat, db_tgc_bylat_median * 2), color = "red") +
  theme_classic() +
  labs(
    x = "Latitude (°)",
    y = expression(paste("C sink (TgC yr"^{-1}, " °" ^{-1},")"))
  ) +
  coord_flip()

gg_csink_bylat

ggsave(
  here("manuscript/figures/gg_csink_bylat.pdf"),
  plot = gg_csink_bylat,
  width = 5,
  height = 4
)

## Publication figure --------------------------
gg_map_sink_all <- cowplot::plot_grid(
  gg_csink_bylat,
  gg_map_sink_pergridarea,
  gg_hist_csink_boot,
  gg_map_sink_pergridgridtarea_sd,
  ncol = 2,
  labels = letters[1:4],
  rel_widths = c(0.3, 1)
)

ggsave(
  here("manuscript/figures/gg_map_sink_all.pdf"),
  plot = gg_map_sink_all,
  width = 11,
  height = 5.5
)

