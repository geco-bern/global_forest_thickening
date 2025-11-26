# This script analyses the changes in the STLs over time (calendar year) by biome.

# Load packages ----------------------------------------------------------------
# library(renv)
library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(rFIA)
library(patchwork)
library(terra)
library(sf)
library(lme4)
library(lmerTest)
library(ggeffects)
library(effects)
library(sjPlot)
library(measurements)
library(sp)
library(lqmm)
library(ggforce)
library(MuMIn)
# library(ingestr)
library(DescTools)
library(here)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(purrr)
library(rsample)
library(cowplot)

# Load functions ---------------------------------------------------------------
source(here("R/functions.R"))
source(here("R/plot_stl_bybiome.R"))
source(here("R/identify_disturbed_plots.R"))
source(here("R/get_breaks.R"))
source(here("R/plot_lqmm_bybiome.R"))
source(here("R/calc_lqmm_byqmdbin.R"))

# Quantile regression ----------------------------------------------------------

# load data
data_unm <- readRDS(here::here("data/inputs/data_unm.rds"))

## Biome 1 Tropical & Subtropical Moist Broadleaf Forests ----------------------
data_unm_biome <- data_unm |>
  filter(biomeID == 1) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  )

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |>
  identify_disturbed_plots()

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome1 <- ggplot(data = data_unm_biome, aes(x = year)) +
  geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
  theme_bw() +
  labs(title = "Tropical Moist Broadleaf Forests") +
  xlim(1960, 2024) +
  labs(x = "Year", y = "Count")

gg_hist_year_biome1

### Plot disturbed plots -------------------------------------------------------
breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |>
  mutate(year_bin = cut(
    year,
    breaks = breaks,
    labels = breaks[1:length(breaks) - 1] + 2.5,
    include.lowest = TRUE
  )) |>
  group_by(year_bin) |>
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

gg_fdisturbed_biome1 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",
    title = bquote(bold("a") ~ ~"Tropical & Subtropical Moist Broadleaf Forests")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.9, # c(0.75, 0.90),
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(
    LP_max_iter = 500, # inner loop iterations
    LP_tol_ll   = 1e-05, # inner loop tolerance
    startQR     = TRUE
  )
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome1.rds"))

#### STL shift -----------------------------------------------------------------
# Opt 2: predicted logDensity at two years differing by one year (in scaled units)
percent_change <- calc_percent_change(data_unm_biome, fit_lqmm)

#### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>%
    group_by(plotID),
  times = 500,
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome1.rds"))

df_boot <- boot_results |>
  filter(term == "year") |>
  mutate(biome = "Tropical & Subtropical Moist Broadleaf Forests")

# boot_results |>
#   filter(term == "year") |>
#   mutate(percent_change = 100*(exp(estimate) - 1)) |>
#   ggplot(aes(x = percent_change)) +
#   geom_density()

# # summarise across bootstraps
# summary_stats <- boot_results |>
#   filter(term == "year") %>%
#   group_by(term) %>%
#   summarise(
#     estimate = mean(estimate),
#     std.error = sd(estimate),
#     ci_low = quantile(estimate, 0.025),
#     ci_high = quantile(estimate, 0.975),
#     .groups = "drop"
#   ) |>
#   mutate(biome = "Tropical & Subtropical Moist Broadleaf Forests")

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome1 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("a") ~ ~"Tropical Moist Broadleaf Forests")
)
gg_lqmm_biome1

### Within QMD bins ------------------------------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome1_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

gg_lqmm_biome1_both <- cowplot::plot_grid(
  gg_lqmm_biome1,
  gg_lqmm_biome1_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  align = "v",
  # labels = c("",  "g"),
  label_y = 1.1
)
gg_lqmm_biome1_both

## Biome 2 Tropical & Subtropical Dry Broadleaf Forests ------------------------
data_unm_biome <- data_unm |>
  filter(biomeID == 2) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  )

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |>
  identify_disturbed_plots()

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome2 <- ggplot(data = data_unm_biome, aes(x = year)) +
  geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
  theme_bw() +
  labs(title = "Tropical Dry Broadleaf Forests") +
  xlim(1980, 2024) +
  labs(x = "Year", y = "Count")

gg_hist_year_biome2

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(
    LP_max_iter = 500, # inner loop iterations
    LP_tol_ll   = 1e-05, # inner loop tolerance
    startQR     = TRUE
  )
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome2.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>%
    group_by(plotID),
  times = 500,
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome2.rds"))

df_boot <- bind_rows(
  df_boot,
  boot_results |>
    filter(term == "year") |>
    mutate(biome = "Tropical Dry Broadleaf Forests")
)

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome2 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("b") ~ ~"Tropical Dry Broadleaf Forests")
)
gg_lqmm_biome2

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome2_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

gg_lqmm_biome2_both <- cowplot::plot_grid(
  gg_lqmm_biome2,
  gg_lqmm_biome2_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  align = "v",
  # labels = c("",  "h"),
  label_y = 1.1
)
gg_lqmm_biome2_both

## Biome 4 Temperate Broadleaf & Mixed Forests ---------------------------------
data_unm_biome <- data_unm |>
  filter(biomeID == 4) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  )

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |>
  identify_disturbed_plots()

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome4 <- ggplot(data = data_unm_biome, aes(x = year)) +
  geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
  theme_bw() +
  labs(title = "Temperate Broadleaf & Mixed Forests") +
  xlim(1930, 2024) +
  labs(x = "Year", y = "Count")

gg_hist_year_biome4

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(startQR = TRUE)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome4.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>%
    group_by(plotID),
  times = 500,
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome4.rds"))

df_boot <- bind_rows(
  df_boot,
  boot_results |>
    filter(term == "year") |>
    mutate(biome = "Temperate Broadleaf & Mixed Forests")
)
### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome4 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("c") ~ ~"Temperate Broadleaf & Mixed Forests")
)
gg_lqmm_biome4

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome4_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

gg_lqmm_biome4_both <- cowplot::plot_grid(
  gg_lqmm_biome4,
  gg_lqmm_biome4_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  # labels = c("",  "i"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome4_both

## Biome 5  Temperate Conifer Forests Forest -----------------------------------
data_unm_biome <- data_unm |>
  filter(biomeID == 5) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  )

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |>
  identify_disturbed_plots()

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome5 <- ggplot(data = data_unm_biome, aes(x = year)) +
  geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
  theme_bw() +
  labs(title = "Temperate Conifer Forests") +
  xlim(1960, 2024) +
  labs(x = "Year", y = "Count")

gg_hist_year_biome5

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(startQR = TRUE)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome5.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>%
    group_by(plotID),
  times = 500,
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome5.rds"))

df_boot <- bind_rows(
  df_boot,
  boot_results |>
    filter(term == "year") |>
    mutate(biome = "Temperate Conifer Forests")
)

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome5 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("d") ~ ~"Temperate Conifer Forest")
)
gg_lqmm_biome5

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome5_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

gg_lqmm_biome5_both <- cowplot::plot_grid(
  gg_lqmm_biome5,
  gg_lqmm_biome5_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  # labels = c("",  "j"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome5_both

## Biome 6 Boreal Forests/Taiga ------------------------------------------------
data_unm_biome <- data_unm |>
  filter(biomeID == 6) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  )

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |>
  identify_disturbed_plots()

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome6 <- ggplot(data = data_unm_biome, aes(x = year)) +
  geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
  theme_bw() +
  labs(title = "Boreal Forests/Taiga") +
  xlim(1980, 2024) +
  labs(x = "Year", y = "Count")

gg_hist_year_biome6

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(
    LP_max_iter = 500, # increase max iterations
    LP_tol_ll = 1e-4, # relax tolerance slightly (default is 1e-5)
    startQR = TRUE # good to keep this TRUE
  )
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome6.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>%
    group_by(plotID),
  times = 500,
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome6.rds"))

df_boot <- bind_rows(
  df_boot,
  boot_results |>
    filter(term == "year") |>
    mutate(biome = "Boreal Forests/Taiga")
)

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome6 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("e") ~ ~"Boreal Forests/Taiga")
)
gg_lqmm_biome6

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome6_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

gg_lqmm_biome6_both <- cowplot::plot_grid(
  gg_lqmm_biome6,
  gg_lqmm_biome6_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  # labels = c("",  "k"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome6_both

## Biome 12 Mediterranean Forests ----------------------
data_unm_biome <- data_unm |>
  filter(biomeID == 12) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  )

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |>
  identify_disturbed_plots()

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome12 <- ggplot(data = data_unm_biome, aes(x = year)) +
  geom_histogram(color = "black", fill = "grey80", binwidth = 1) +
  theme_bw() +
  labs(title = "Mediterranean Forests") +
  xlim(1960, 2024) +
  labs(x = "Year", y = "Count")

gg_hist_year_biome12

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(startQR = TRUE)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome12.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>%
    group_by(plotID),
  times = 500,
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome12.rds"))

df_boot <- bind_rows(
  df_boot,
  boot_results |>
    filter(term == "year") |>
    mutate(biome = "Mediterranean Forests")
)

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome12 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm,
  name = bquote(bold("f") ~ ~"Mediterranean Forests")
)
gg_lqmm_biome12

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome12_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

gg_lqmm_biome12_both <- cowplot::plot_grid(
  gg_lqmm_biome12,
  gg_lqmm_biome12_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.6),
  # labels = c("",  "l"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome12_both

# Publication figures ----------------------------------------------------------

## Figure 1 --------------------------------------------------------------------

legend <- get_legend(
  gg_lqmm_biome1 +
    theme(legend.position = "right")
)

# Arrange the 9 plots in a 3x3 grid
fig1_lqmm <- cowplot::plot_grid(
  gg_lqmm_biome1_both,
  gg_lqmm_biome2_both,
  gg_lqmm_biome4_both,
  gg_lqmm_biome5_both,
  gg_lqmm_biome6_both,
  gg_lqmm_biome12_both,
  ncol = 3
)
fig1_lqmm

# Combine grid and legend
fig1_lqmm <- cowplot::plot_grid(
  fig1_lqmm,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

ggsave(
  filename = here::here("manuscript/figures/fig1_lqmm.pdf"),
  plot = fig1_lqmm,
  width = 11,
  height = 12
)

ggsave(
  filename = here::here("manuscript/figures/fig1_lqmm.png"),
  plot = fig1_lqmm,
  width = 11,
  height = 12
)

## SI Figure: histogram over years ---------------------------------------------
fig_hist_year <- cowplot::plot_grid(
  gg_hist_year_biome1,
  gg_hist_year_biome2,
  gg_hist_year_biome4,
  gg_hist_year_biome5,
  gg_hist_year_biome6,
  gg_hist_year_biome12,
  ncol = 3,
  labels = letters
)
fig_hist_year

ggsave(
  filename = here::here("manuscript/figures/fig_hist_year.pdf"),
  plot = fig_hist_year,
  width = 11,
  height = 6
)

ggsave(
  filename = here::here("manuscript/figures/fig_hist_year.png"),
  plot = fig_hist_year,
  width = 11,
  height = 6
)

## SI Figure: Bootstrapped percent change of N per year ------------------------
boot_results |>
  filter(term == "year") |>
  mutate(percent_change = 100*(exp(estimate) - 1)) |>
  ggplot(aes(x = percent_change, group = biome, color = biome)) +
  geom_density()

