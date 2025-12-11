# This script analyses the changes in the STLs over time (calendar year) by biome.

# Load packages ----------------------------------------------------------------
# library(renv)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(rFIA)
library(lme4)
library(lmerTest)
library(ggeffects)
library(effects)
library(sjPlot)
library(measurements)
library(lqmm)
library(ggforce)
library(MuMIn)
library(DescTools)
library(here)
library(viridis)
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
source(here("R/create_table_latex.R"))
source(here("R/process_cite_lines.R"))
source(here("R/wrap_fit_lqmm.R"))
source(here("R/calc_percent_change.R"))

# Quantile regression ----------------------------------------------------------

# load data
data_unm <- readRDS(here("data/inputs/data_unm.rds"))

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

### Histogram of data over years -----------------------------------------------
# gg_hist_year_biome1 <- ggplot(data = data_unm_biome, aes(x = year)) +
#   geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
#   theme_bw() +
#   labs(title = "Tropical Moist Broadleaf Forests") +
#   xlim(1970, 2024) +
#   labs(x = "Year", y = "Count")

# distinguishing datasets
gg_hist_year_biome1 <- ggplot(data = data_unm_biome, aes(x = year, fill = dataset)) +
  geom_histogram(color = "black", binwidth = 1, position = "stack", linewidth = 0.3) +
  theme_bw() +
  labs(title = "Tropical Moist Broadleaf Forests") +
  xlim(1970, 2024) +
  labs(x = "Year", y = "Count", fill = "") +
  khroma::scale_fill_okabeito() +
  theme(
    legend.position = "right"
  )

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
  filter(nplots > 30) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

gg_fdisturbed_biome1 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",
    title = "Tropical & Subtropical Moist Broadleaf Forests"
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

gg_fdisturbed_biome1

### Remove disturbed plots -----------------------------------------------------
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_biome1.rds")
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

write_rds(fit_lqmm, file = here("data/outputs/fit_lqmm_biome1.rds"))

#### STL shift -----------------------------------------------------------------
# Opt 2: predicted logDensity at two years differing by one year (in scaled units)
percent_change <- calc_percent_change(data_unm_biome, fit_lqmm)

#### Bootstrapping LQMM fit ----------------------------------------------------
# boot_data <- rsample::bootstraps(
#   data_unm_biome %>%
#     group_by(plotID),
#   times = 500,
#   apparent = FALSE
# )
#
# # Apply model to each bootstrap sample
# boot_results <- boot_data %>%
#   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
#   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
#   unnest(coefs) |>
#   dplyr::select(-splits)
#
# write_rds(boot_results, file = here("data/boot_results_biome1.rds"))
#
# df_boot <- boot_results |>
#   filter(term == "year") |>
#   mutate(biome = "Tropical & Subtropical Moist Broadleaf Forests")
#
# boot_results <- read_rds(file = here("data/boot_results_biome1.rds"))

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

### Histogram of data over years -----------------------------------------------
# gg_hist_year_biome2 <- ggplot(data = data_unm_biome, aes(x = year)) +
#   geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
#   theme_bw() +
#   labs(title = "Tropical Dry Broadleaf Forests") +
#   xlim(1990, 2024) +
#   labs(x = "Year", y = "Count")

# distinguishing datasets
gg_hist_year_biome2 <- ggplot(data = data_unm_biome, aes(x = year, fill = dataset)) +
  geom_histogram(color = "black", binwidth = 1, position = "stack", linewidth = 0.3) +
  theme_bw() +
  labs(title = "Tropical Dry Broadleaf Forests") +
  xlim(1990, 2024) +
  labs(x = "Year", y = "Count", fill = "") +
  khroma::scale_fill_okabeito() +
  theme(
    legend.position = "right"
  )

gg_hist_year_biome2

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
  filter(nplots > 30) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

gg_fdisturbed_biome2 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  # geom_smooth(method = "lm", color = "red") + # only two points
  theme_classic() +
  labs(
    x = "Year",
    title = "Tropical Dry Broadleaf Forests"
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

gg_fdisturbed_biome2

### Remove disturbed plots -----------------------------------------------------
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_biome2.rds")
)

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

write_rds(fit_lqmm, file = here("data/outputs/fit_lqmm_biome2.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
# boot_data <- rsample::bootstraps(
#   data_unm_biome %>%
#     group_by(plotID),
#   times = 500,
#   apparent = FALSE
# )
#
# # Apply model to each bootstrap sample
# boot_results <- boot_data %>%
#   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
#   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
#   unnest(coefs) |>
#   dplyr::select(-splits)
#
# write_rds(boot_results, file = here("data/boot_results_biome2.rds"))
#
# df_boot <- bind_rows(
#   df_boot,
#   boot_results |>
#     filter(term == "year") |>
#     mutate(biome = "Tropical Dry Broadleaf Forests")
# )

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

### Histogram of data over years -----------------------------------------------
# gg_hist_year_biome4 <- ggplot(data = data_unm_biome, aes(x = year)) +
#   geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
#   theme_bw() +
#   labs(title = "Temperate Broadleaf & Mixed Forests") +
#   xlim(1930, 2024) +
#   labs(x = "Year", y = "Count")

# distinguishing datasets
gg_hist_year_biome4 <- ggplot(data = data_unm_biome, aes(x = year, fill = dataset)) +
  geom_histogram(color = "black", binwidth = 1, position = "stack", linewidth = 0.3) +
  theme_bw() +
  labs(title = "Temperate Broadleaf & Mixed Forests") +
  xlim(1960, 2024) +
  labs(x = "Year", y = "Count", fill = "") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme(
    legend.position = "right"
  )

gg_hist_year_biome4

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
  filter(nplots > 30) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

gg_fdisturbed_biome4 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",
    title = "Temperate Broadleaf & Mixed Forests"
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

gg_fdisturbed_biome4

### Remove disturbed plots -----------------------------------------------------
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_biome4.rds")
  )

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

write_rds(fit_lqmm, file = here("data/outputs/fit_lqmm_biome4.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
# boot_data <- rsample::bootstraps(
#   data_unm_biome %>%
#     group_by(plotID),
#   times = 500,
#   apparent = FALSE
# )
#
# # Apply model to each bootstrap sample
# boot_results <- boot_data %>%
#   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
#   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
#   unnest(coefs) |>
#   dplyr::select(-splits)
#
# write_rds(boot_results, file = here("data/boot_results_biome4.rds"))
#
# df_boot <- bind_rows(
#   df_boot,
#   boot_results |>
#     filter(term == "year") |>
#     mutate(biome = "Temperate Broadleaf & Mixed Forests")
# )

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

### Histogram of data over years -----------------------------------------------
# gg_hist_year_biome5 <- ggplot(data = data_unm_biome, aes(x = year)) +
#   geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
#   theme_bw() +
#   labs(title = "Temperate Conifer Forests") +
#   xlim(1960, 2024) +
#   labs(x = "Year", y = "Count")

# distinguishing datasets
gg_hist_year_biome5 <- ggplot(data = data_unm_biome, aes(x = year, fill = dataset)) +
  geom_histogram(color = "black", binwidth = 1, position = "stack", linewidth = 0.3) +
  theme_bw() +
  labs(title = "Temperate Conifer Forests") +
  xlim(1970, 2024) +
  labs(x = "Year", y = "Count", fill = "") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme(
    legend.position = "right"
  )

gg_hist_year_biome5

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
  filter(nplots > 30) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

gg_fdisturbed_biome5 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",
    title = "Temperate Conifer Forests"
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

gg_fdisturbed_biome5

### Remove disturbed plots -----------------------------------------------------
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_biome5.rds")
)

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

write_rds(fit_lqmm, file = here("data/outputs/fit_lqmm_biome5.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
# boot_data <- rsample::bootstraps(
#   data_unm_biome %>%
#     group_by(plotID),
#   times = 500,
#   apparent = FALSE
# )
#
# # Apply model to each bootstrap sample
# boot_results <- boot_data %>%
#   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
#   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
#   unnest(coefs) |>
#   dplyr::select(-splits)
#
# write_rds(boot_results, file = here("data/boot_results_biome5.rds"))
#
# df_boot <- bind_rows(
#   df_boot,
#   boot_results |>
#     filter(term == "year") |>
#     mutate(biome = "Temperate Conifer Forests")
# )

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

### Histogram of data over years -----------------------------------------------
# gg_hist_year_biome6 <- ggplot(data = data_unm_biome, aes(x = year)) +
#   geom_histogram(fill = "grey80", color = "black", binwidth = 1) +
#   theme_bw() +
#   labs(title = "Boreal Forests/Taiga") +
#   xlim(1980, 2024) +
#   labs(x = "Year", y = "Count")

# distinguishing datasets
gg_hist_year_biome6 <- ggplot(data = data_unm_biome, aes(x = year, fill = dataset)) +
  geom_histogram(color = "black", binwidth = 1, position = "stack", linewidth = 0.3) +
  theme_bw() +
  labs(title = "Boreal Forests/Taiga") +
  xlim(1980, 2024) +
  labs(x = "Year", y = "Count", fill = "") +
  khroma::scale_fill_okabeito() +
  theme(
    legend.position = "right"
  )

gg_hist_year_biome6

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
  filter(nplots > 30) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

gg_fdisturbed_biome6 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",
    title = "Boreal Forests/Taiga"
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

gg_fdisturbed_biome6

### Remove disturbed plots -----------------------------------------------------
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_biome6.rds")
)

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

write_rds(fit_lqmm, file = here("data/outputs/fit_lqmm_biome6.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
# boot_data <- rsample::bootstraps(
#   data_unm_biome %>%
#     group_by(plotID),
#   times = 500,
#   apparent = FALSE
# )
#
# # Apply model to each bootstrap sample
# boot_results <- boot_data %>%
#   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
#   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
#   unnest(coefs) |>
#   dplyr::select(-splits)
#
# write_rds(boot_results, file = here("data/boot_results_biome6.rds"))
#
# df_boot <- bind_rows(
#   df_boot,
#   boot_results |>
#     filter(term == "year") |>
#     mutate(biome = "Boreal Forests/Taiga")
# )

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

### Histogram of data over years -----------------------------------------------
# gg_hist_year_biome12 <- ggplot(data = data_unm_biome, aes(x = year)) +
#   geom_histogram(color = "black", fill = "grey80", binwidth = 1) +
#   theme_bw() +
#   labs(title = "Mediterranean Forests") +
#   xlim(1980, 2024) +
#   labs(x = "Year", y = "Count")

# distinguishing datasets
gg_hist_year_biome12 <- ggplot(data = data_unm_biome, aes(x = year, fill = dataset)) +
  geom_histogram(color = "black", binwidth = 1, position = "stack", linewidth = 0.3) +
  theme_bw() +
  labs(title = "Mediterranean Forests") +
  xlim(1980, 2024) +
  labs(x = "Year", y = "Count", fill = "") +
  khroma::scale_fill_okabeito() +
  theme(
    legend.position = "right"
  )

gg_hist_year_biome12

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
  filter(nplots > 30) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

gg_fdisturbed_biome12 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",
    title = "Mediterranean Forests"
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

gg_fdisturbed_biome12

### Remove disturbed plots -----------------------------------------------------
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |>
  filter(ndisturbed == 0)

write_rds(
  data_unm_biome,
  file = here("data/data_unm_undist_biome12.rds")
)

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

write_rds(fit_lqmm, file = here("data/outputs/fit_lqmm_biome12.rds"))

#### Bootstrapping LQMM fit -----------------------------------------------------
# boot_data <- rsample::bootstraps(
#   data_unm_biome %>%
#     group_by(plotID),
#   times = 500,
#   apparent = FALSE
# )
#
# # Apply model to each bootstrap sample
# boot_results <- boot_data %>%
#   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
#   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
#   unnest(coefs) |>
#   dplyr::select(-splits)
#
# write_rds(boot_results, file = here("data/boot_results_biome12.rds"))
#
# df_boot <- bind_rows(
#   df_boot,
#   boot_results |>
#     filter(term == "year") |>
#     mutate(biome = "Mediterranean Forests")
# )

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

# Publication figures and tables  ----------------------------------------------

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
  filename = here("manuscript/figures/fig1_lqmm.pdf"),
  plot = fig1_lqmm,
  width = 11,
  height = 12
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm.png"),
  plot = fig1_lqmm,
  width = 11,
  height = 12
)

## SI Figure: Histogram over years ---------------------------------------------
row1 <- cowplot::plot_grid(
  gg_hist_year_biome1,
  gg_hist_year_biome2,
  ncol = 2,
  rel_widths = c(1, 0.7),
  labels = letters[1:2]
)

row2 <- cowplot::plot_grid(
  gg_hist_year_biome4,
  ncol = 1,
  labels = letters[3]
)

row3 <- cowplot::plot_grid(
  gg_hist_year_biome5,
  ncol = 1,
  labels = letters[4]
)

row4 <- cowplot::plot_grid(
  gg_hist_year_biome6,
  gg_hist_year_biome12,
  ncol = 2,
  labels = letters[5:6]
)

fig_hist_year <- cowplot::plot_grid(
  row1,
  row2,
  row3,
  row4,
  ncol = 1
)
fig_hist_year

ggsave(
  filename = here("manuscript/figures/fig_hist_year.pdf"),
  plot = fig_hist_year,
  width = 9,
  height = 15
)

## SI Figure: Length of repeated observations ----------------------------------
df_len <- data_unm |>
  group_by(plotID, biome, biomeID) |>
  summarise(start = min(year), end = max(year)) |>
  mutate(len = end - start)

# numbers for paper
df_len |>
  ungroup() |>
  group_by(biomeID, biome) |>
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
  ggplot(aes(x = len, color = biome, fill = biome)) +
  geom_density(adjust = 3, alpha = 0.5) +
  scale_fill_manual(
    values = c(
      "Boreal Forests/Taiga"                                       = "dodgerblue4",
      "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
      "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
      "Temperate Conifer Forests"                                  = "lightseagreen",
      "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
      "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
    ),
    na.value = NA,
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  scale_color_manual(
    values = c(
      "Boreal Forests/Taiga"                                       = "dodgerblue4",
      "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
      "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
      "Temperate Conifer Forests"                                  = "lightseagreen",
      "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
      "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
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
  filename = here("manuscript/figures/distribution_length.pdf"),
  width = 10,
  height = 5
)

# ## SI Figure: Bootstrapped percent change of N per year ------------------------
# write_rds(df_boot, file = here("data/df_boot.rds"))
#
# df_boot |>
#   mutate(percent_change = 100*(exp(estimate) - 1)) |>
#   ggplot(aes(x = percent_change, group = biome, color = biome, fill = biome)) +
#   geom_density(alpha = 0.5) +
#   scale_fill_manual(
#     values = c(
#       "Boreal Forests/Taiga"                                       = "dodgerblue4",
#       "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
#       "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
#       "Temperate Conifer Forests"                                  = "lightseagreen",
#       "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
#       "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
#     ),
#     na.value = NA,
#     breaks = ~ .x[!is.na(.x)],
#     name = ""
#   ) +
#   scale_color_manual(
#     values = c(
#       "Boreal Forests/Taiga"                                       = "dodgerblue4",
#       "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
#       "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
#       "Temperate Conifer Forests"                                  = "lightseagreen",
#       "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
#       "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
#     ),
#     na.value = NA,
#     breaks = ~ .x[!is.na(.x)],
#     name = ""
#   ) +
#   theme_classic() +
#   labs(
#     x = expression(paste("Change in density (%"^{-yr}, ")")),
#     y = "Density",
#     color = "",
#     fill = ""
#     ) +
#   theme(
#     legend.position = "bottom"
#   )
#
# ggsave(
#   filename = here("manuscript/figures/distribution_percent_change.pdf"),
#   width = 8,
#   height = 4
# )
#
# # table of bootstrapped estimates for coefficient and percent change
# summary_stats <- df_boot |>
#   mutate(percent_change = 100*(exp(estimate) - 1)) |>
#   ungroup() |>
#   group_by(biome) %>%
#   summarise(
#     estimate_mean = mean(estimate),
#     estimate_sd = sd(estimate),
#     estimate_ci_low = quantile(estimate, 0.025),
#     estimate_ci_high = quantile(estimate, 0.975),
#
#     percent_change_mean = mean(percent_change),
#     percent_change_sd = sd(percent_change),
#     percent_change_ci_low = quantile(percent_change, 0.025),
#     percent_change_ci_high = quantile(percent_change, 0.975),
#
#     .groups = "drop"
#   )
#
# write_rds(summary_stats, file = here("data/summary_stats.csv"))
#
# create_table_latex(
#   summary_stats |>
#     select(
#       Biome = biome,
#       Mean = percent_change_mean,
#       SD = percent_change_sd
#       ),
#     caption = "Percentage change of forest stand density.",
#     filn = here("manuscript/tables/table_percentage_change.tex")
#     # align = c("p{0.1cm}", "p{5cm}", "p{7cm}")
#     )

## SI Figure: Disturbed plots --------------------------------------------------
cowplot::plot_grid(
  gg_fdisturbed_biome1,
  gg_fdisturbed_biome2,
  gg_fdisturbed_biome4,
  gg_fdisturbed_biome5,
  gg_fdisturbed_biome6,
  gg_fdisturbed_biome12,
  ncol = 2,
  labels = letters
)

ggsave(
  filename = here("manuscript/figures/fdisturbed.pdf"),
  width = 9,
  height = 9
)

## Table: datasets
df_datasets <- data_unm |>
  group_by(dataset) |>
  summarise(count = n()) |>
  arrange(-count) |>
  rename(Dataset = dataset, N = count) |>
  mutate(Description = "", Reference = "")

create_table_latex(
  df_datasets,
  caption = "Constituent forest dataset sizes and descriptions.",
  filn = here("manuscript/tables/datasets.tex")
  # align = c("p{0.1cm}", "p{5cm}", "p{7cm}")
)
