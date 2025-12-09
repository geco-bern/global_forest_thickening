# This script identifies and analyses plots subject to self-thinning
# Requires 02_stl_year.R to be run first.
source(here("R/filter_stl_slope.R"))
source(here("R/plot_stl_longplots.R"))

# Biome 4 Temperate Broadleaf & Mixed Forests ---------------------------------
## Read and filter data------------------------------
# data filtered by disturbance and file written in 02_stl_year.R
data_unm_biome <- read_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_biome4.rds")
)

fit_lqmm <- read_rds(here("data/outputs/fit_lqmm_biome4.rds"))

data_unm_biome <- filter_stl_slope(data_unm_biome, fit_lqmm)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_slopefilter_biome4.rds")
)

## Visualisation ---------------------------------------------
### Distribution of slopes --------
# distribution of slopes across plots
slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]

gg_slope_distribution_biome4 <- ggplot() +
  geom_histogram(
    aes(slope, after_stat(density)),
    data = data_unm_biome |>
      select(plotID, slope) |>
      distinct(),
    fill = "grey70",
    color = "black"
  ) +
  geom_vline(
    xintercept = slope_lqmm,
    color = "red"
  ) +
  theme_bw() +
  labs(
    title = "Temperate Broadleaf & Mixed Forests",
    x = "Slope (unitless)",
    y = "Density"
  )

### Longest-running forest plots --------------
# Select plots with longest time series
gg_longplots_biome4 <- plot_stl_longplots(data_unm_biome, biome_name = "Temperate Broadleaf & Mixed Forests")
