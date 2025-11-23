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
    LP_max_iter = 5000, # inner loop iterations
    LP_tol_ll   = 1e-05, # inner loop tolerance
    startQR     = TRUE
  )
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome1.rds"))

#### STL shift -----------------------------------------------------------------
# Estimated change in N per unit increase in year

out <- summary(fit_lqmm)

# Extract model coefficient
beta_year_sc <- out$tTable[c("year_sc"), "Value"]

# SD of the original (unscaled) year variable
sd_year <- sd(data_unm_biome$year, na.rm = TRUE)

# Change in logDensity per calendar year
real_coef_year <- beta_year_sc * sd_year

# Convert to % change in tree density per year
percent_change_per_year <- (exp(real_coef_year) - 1) * 100

# Opt 2: predicted logDensity at two years differing by one year (in scaled units)

# Pick a fixed value for logQMD_sc (e.g., mean)
mean_logQMD_sc <- mean(data_unm_biome$logQMD_sc, na.rm = TRUE)
sd_year <- sd(data_unm_biome$year, na.rm = TRUE)

# Predict at year_sc = 0 and year_sc = 1 / sd_year (since one calendar year corresponds to 1/sd_year in scaled units)
newdata1 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 0, plotID = data_unm_biome$plotID[1])
newdata2 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 1 / sd_year, plotID = data_unm_biome$plotID[1])

pred1 <- predict(fit_lqmm, newdata = newdata1)
pred2 <- predict(fit_lqmm, newdata = newdata2)

delta_logDensity <- pred2 - pred1
percent_change <- (exp(delta_logDensity) - 1) * 100

#### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>%
    group_by(plotID),
  times = 5000,
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

#write_rds(boot_results, file = here::here("data/boot_results_biome1.rds"))

boot_results |>
  ggplot(aes(estimate)) +
  geom_density() +
  facet_wrap(~term, scales = "free", ncol = 1)

# summarise across bootstraps
summary_stats <- boot_results %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    std.error = sd(estimate),
    ci_low = quantile(estimate, 0.025),
    ci_high = quantile(estimate, 0.975),
    .groups = "drop"
  )

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
  rel_heights = c(1, 0.4),
  align = "v",
  # labels = c("",  "g"),
  label_y = 1.1
)
gg_lqmm_biome1_both

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome1 <- ggplot(data_unm_biome, aes(x=year)) + 
  geom_histogram(color="#FFDB6D", fill="#FFDB6D") + theme_bw() + 
  labs(title = "Tropical Moist Broadleaf Forests") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10)) + 
  scale_x_continuous("Year", limits = c(1960,2024), breaks = seq(1940,2020,20)) +
  scale_y_continuous("Frequency", limits = c(0,650), breaks = seq(0,600,200))
gg_hist_year_biome1

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
    LP_max_iter = 5000, # inner loop iterations
    LP_tol_ll   = 1e-05, # inner loop tolerance
    startQR     = TRUE
  )
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome2.rds"))

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
  rel_heights = c(1, 0.4),
  align = "v",
  # labels = c("",  "h"),
  label_y = 1.1
)
gg_lqmm_biome2_both

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome2 <- ggplot(data_unm_biome, aes(x=year)) + 
  geom_histogram(color="#FFDB6D", fill="#FFDB6D") + theme_bw() + 
  labs(title = "Tropical Dry Broadleaf Forests") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10)) + 
  scale_x_continuous("Year", limits = c(1980,2024), breaks = seq(1940,2020,20)) +
  scale_y_continuous("Frequency", limits = c(0,60), breaks = seq(0,100,20))
gg_hist_year_biome2

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
  rel_heights = c(1, 0.4),
  # labels = c("",  "i"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome4_both

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome4 <- ggplot(data_unm_biome, aes(x=year)) + 
  geom_histogram(color="#FFDB6D", fill="#FFDB6D") + theme_bw() + 
  labs(title = "Temperate Broadleaf & Mixed Forests") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10)) + 
  scale_x_continuous("Year", limits = c(1930,2024), breaks = seq(1940,2020,20)) +
  scale_y_continuous("Frequency", limits = c(0,4100), breaks = seq(0,4000,1000))
gg_hist_year_biome4

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
  rel_heights = c(1, 0.4),
  # labels = c("",  "j"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome5_both

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome5 <- ggplot(data_unm_biome, aes(x=year)) + 
  geom_histogram(color="#FFDB6D", fill="#FFDB6D") + theme_bw() + 
  labs(title = "Temperate Conifer Forests") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10)) + 
  scale_x_continuous("Year", limits = c(1960,2024), breaks = seq(1940,2020,20)) +
  scale_y_continuous("Frequency", limits = c(0,4100), breaks = seq(0,4000,1000))
gg_hist_year_biome5

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
    LP_max_iter = 5000, # increase max iterations
    LP_tol_ll = 1e-4, # relax tolerance slightly (default is 1e-5)
    startQR = TRUE # good to keep this TRUE
  )
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/outputs/fit_lqmm_biome6.rds"))

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
  rel_heights = c(1, 0.4),
  # labels = c("",  "k"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome6_both

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome6 <- ggplot(data_unm_biome, aes(x=year)) + 
  geom_histogram(color="#FFDB6D", fill="#FFDB6D") + theme_bw() + 
  labs(title = "Boreal Forests/Taiga") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10)) + 
  scale_x_continuous("Year", limits = c(1980,2024), breaks = seq(1940,2020,20)) +
  scale_y_continuous("Frequency", limits = c(0,2100), breaks = seq(0,2000,500))
gg_hist_year_biome6

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
  rel_heights = c(1, 0.4),
  # labels = c("",  "l"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome12_both

### Histogram of data over years -----------------------------------------------
gg_hist_year_biome12 <- ggplot(data_unm_biome, aes(x=year)) + 
  geom_histogram(color="#FFDB6D", fill="#FFDB6D") + theme_bw() + 
  labs(title = "Mediterranean Forests") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10)) + 
  scale_x_continuous("Year", limits = c(1960,2024), breaks = seq(1940,2020,20)) +
  scale_y_continuous("Frequency", limits = c(0,4500), breaks = seq(0,5000,1000))
gg_hist_year_biome12

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
  height = 10
)

ggsave(
  filename = here::here("manuscript/figures/fig1_lqmm.png"),
  plot = fig1_lqmm,
  width = 11,
  height = 10
)

## SI Figure histogram over years ----------------------------------------------
fig_hist_year <- cowplot::plot_grid(
  gg_hist_year_biome1,
  gg_hist_year_biome2,
  gg_hist_year_biome4,
  gg_hist_year_biome5,
  gg_hist_year_biome6,
  gg_hist_year_biome12,
  ncol = 3
)
fig_hist_year

fig_hist_year <- gg_hist_year_biome1 +
gg_hist_year_biome2 +
gg_hist_year_biome4 +
gg_hist_year_biome5 +
gg_hist_year_biome6 +
gg_hist_year_biome12 +
  plot_layout(ncol = 3) + 
  plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 10)) 

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
