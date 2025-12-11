# This script identifies and analyses plots subject to self-thinning
# Requires 02_stl_year.R to be run first.
library(tidyverse)
library(here)

source(here("R/filter_stl_slope.R"))
source(here("R/plot_stl_longplots.R"))

# Biome 1 Tropical & Subtropical Moist Broadleaf Forests ---------------------------------
## Read and filter data------------------------------
# data filtered by disturbance and file written in 02_stl_year.R
data_unm_biome <- read_rds(
  file = here("data/data_unm_undist_biome1.rds")
)

fit_lqmm <- read_rds(here("data/outputs/fit_lqmm_biome1.rds"))

data_unm_biome <- filter_stl_slope(data_unm_biome, fit_lqmm)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_slopefilter_biome1.rds")
)

data_unm_slopefilter <- data_unm_biome |>
  filter(filter_slope) |>
  unnest(data)

## Visualisation ---------------------------------------------
### Distribution of slopes --------
# distribution of slopes across plots
slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]

gg_slope_distribution_biome1 <- ggplot() +
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
    title = "Tropical & Subtropical Moist Broadleaf Forests",
    x = "Slope (unitless)",
    y = "Density"
  )

### Longest-running forest plots --------------
# Select plots with longest time series
gg_longplots_biome1 <- plot_stl_longplots(data_unm_biome, biome_name = "Tropical & Subtropical Moist Broadleaf Forests")

# Biome 2 Tropical & Subtropical Moist Broadleaf Forests ---------------------------------
## Read and filter data------------------------------
# data filtered by disturbance and file written in 02_stl_year.R
data_unm_biome <- read_rds(
  file = here("data/data_unm_undist_biome2.rds")
)

fit_lqmm <- read_rds(here("data/outputs/fit_lqmm_biome2.rds"))

data_unm_biome <- filter_stl_slope(data_unm_biome, fit_lqmm)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_slopefilter_biome2.rds")
)

data_unm_slopefilter <- data_unm_slopefilter |>
  bind_rows(
    data_unm_biome |>
      filter(filter_slope) |>
      unnest(data)
  )

## Visualisation ---------------------------------------------
### Distribution of slopes --------
# distribution of slopes across plots
slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]

gg_slope_distribution_biome2 <- ggplot() +
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
    title = "Tropical & Subtropical Dry Broadleaf Forests",
    x = "Slope (unitless)",
    y = "Density"
  )

### Longest-running forest plots --------------
# Select plots with longest time series
gg_longplots_biome2 <- plot_stl_longplots(
  data_unm_biome,
  biome_name = "Tropical & Subtropical Dry Broadleaf Forests"
  )

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

data_unm_slopefilter <- data_unm_slopefilter |>
  bind_rows(
    data_unm_biome |>
      filter(filter_slope) |>
      unnest(data)
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

# Biome 5 Temperate Conifer Forests Forest ---------------------------------
## Read and filter data------------------------------
# data filtered by disturbance and file written in 02_stl_year.R
data_unm_biome <- read_rds(
  file = here("data/data_unm_undist_biome5.rds")
)

fit_lqmm <- read_rds(here("data/outputs/fit_lqmm_biome5.rds"))

data_unm_biome <- filter_stl_slope(data_unm_biome, fit_lqmm)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_slopefilter_biome5.rds")
)

data_unm_slopefilter <- data_unm_slopefilter |>
  bind_rows(
    data_unm_biome |>
      filter(filter_slope) |>
      unnest(data)
  )

## Visualisation ---------------------------------------------
### Distribution of slopes --------
# distribution of slopes across plots
slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]

gg_slope_distribution_biome5 <- ggplot() +
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
    title = "Temperate Conifer Forests Forest",
    x = "Slope (unitless)",
    y = "Density"
  )

### Longest-running forest plots --------------
# Select plots with longest time series
gg_longplots_biome5 <- plot_stl_longplots(
  data_unm_biome,
  biome_name = "Temperate Conifer Forests Forest"
)


# Biome 6 Boreal Forests/Taiga ---------------------------------
## Read and filter data------------------------------
# data filtered by disturbance and file written in 02_stl_year.R
data_unm_biome <- read_rds(
  file = here("data/data_unm_undist_biome6.rds")
)

fit_lqmm <- read_rds(here("data/outputs/fit_lqmm_biome6.rds"))

data_unm_biome <- filter_stl_slope(data_unm_biome, fit_lqmm)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_slopefilter_biome6.rds")
)

data_unm_slopefilter <- data_unm_slopefilter |>
  bind_rows(
    data_unm_biome |>
      filter(filter_slope) |>
      unnest(data)
  )

## Visualisation ---------------------------------------------
### Distribution of slopes --------
# distribution of slopes across plots
slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]

gg_slope_distribution_biome6 <- ggplot() +
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
    title = "Boreal Forests/Taiga",
    x = "Slope (unitless)",
    y = "Density"
  )

### Longest-running forest plots --------------
# Select plots with longest time series
gg_longplots_biome6 <- plot_stl_longplots(
  data_unm_biome,
  biome_name = "Boreal Forests/Taiga"
)


# Biome 6 Mediterranean Forests ---------------------------------
## Read and filter data------------------------------
# data filtered by disturbance and file written in 02_stl_year.R
data_unm_biome <- read_rds(
  file = here("data/data_unm_undist_biome12.rds")
)

fit_lqmm <- read_rds(here("data/outputs/fit_lqmm_biome12.rds"))

data_unm_biome <- filter_stl_slope(data_unm_biome, fit_lqmm)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_slopefilter_biome12.rds")
)

data_unm_slopefilter <- data_unm_slopefilter |>
  bind_rows(
    data_unm_biome |>
      filter(filter_slope) |>
      unnest(data)
  )

## Visualisation ---------------------------------------------
### Distribution of slopes --------
# distribution of slopes across plots
slope_lqmm <- coef(fit_lqmm)["logQMD_sc"]

gg_slope_distribution_biome12 <- ggplot() +
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
    title = "Mediterranean Forests",
    x = "Slope (unitless)",
    y = "Density"
  )

### Longest-running forest plots --------------
# Select plots with longest time series
gg_longplots_biome12 <- plot_stl_longplots(
  data_unm_biome,
  biome_name = "Mediterranean Forests"
)

## SI Figures --------------
### Distribution of slopes plots -----------------------------------------------
cowplot::plot_grid(
  gg_slope_distribution_biome1,
  gg_slope_distribution_biome2,
  gg_slope_distribution_biome4,
  gg_slope_distribution_biome5,
  gg_slope_distribution_biome6,
  gg_slope_distribution_biome12,
  ncol = 2,
  labels = letters
)

ggsave(
  filename = here("manuscript/figures/slope_distribution.pdf"),
  width = 9,
  height = 9
)

### Longplots STL  --------------------------------------------------
cowplot::plot_grid(
  gg_longplots_biome1,
  gg_longplots_biome2,
  gg_longplots_biome4,
  gg_longplots_biome5,
  gg_longplots_biome6,
  gg_longplots_biome12,
  ncol = 2,
  labels = letters
)

ggsave(
  filename = here("manuscript/figures/stl_longplots.pdf"),
  width = 9,
  height = 9
)

## Write data filtered by slope ------------------------------------------------
write_rds(
  data_unm_slopefilter,
  file = here("data/data_unm_slopefilter.rds")
)
