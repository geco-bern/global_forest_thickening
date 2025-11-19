# Fixed effects plot -----------------------------------------------------------
# This script analyses the main environmental drivers affecting the STL changes

## Load packages ---------------------------------------------------------------
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
library(ingestr)
library(DescTools)
library(corrplot)
library(ggokabeito)

remotes::install_github("https://github.com/valentinitnelav/plotbiomes")
library(plotbiomes)

# load data
data_fil_biomes <- readRDS(here::here("data/data_fil_biomes.rds"))
data_fil_biomes <- readRDS(here::here("data/75perc/data_fil_biomes.rds"))

# correlation among variables
M <- as.matrix(data_fil_biomes[, c(27, 29, 30, 33)] %>% distinct())
corrplot(cor(M, use = "pairwise.complete.obs"), method = "number")

# Data distribution ------------------------------------------------------------
settings_worldclim <- list(varnam = c("tavg", "prec"))

df_plots <- data_fil_biomes |>
  select(sitename = plotID, lon, lat) |>
  distinct()

df_worldclim <- ingest(
  df_plots,
  source    = "worldclim",
  settings  = settings_worldclim,
  dir       = "~/data/archive/worldclim_fick_2017/data/"
)

# mean over months
df_worldclim_mean <- df_worldclim |>
  mutate(
    mat = purrr::map_dbl(data, ~ mean(.$tavg)),
    map = purrr::map_dbl(data, ~ sum(.$prec))
  ) |>
  select(-data)

ggplot() +
  # add biome polygons
  geom_polygon(
    data = Whittaker_biomes,
    aes(
      x = temp_c,
      y = precp_cm,
      fill = biome
    ),
    # adjust polygon borders
    colour = "gray98",
    linewidth = 0.5
  ) +
  theme_bw() +

  # fill the polygons with predefined colors
  scale_fill_manual(
    name = "Whittaker biomes",
    breaks = names(Ricklefs_colors),
    labels = names(Ricklefs_colors),
    values = Ricklefs_colors
  ) +

  # add the temperature - precipitation data points
  geom_point(
    data = df_worldclim_mean,
    aes(
      x = mat,
      y = map / 10
    ),
    alpha = 0.2
  )

# plotbiomes::whittaker_base_plot() +
ggplot() +
  geom_hex(data = df_worldclim_mean, aes(x = mat, y = map / 10), bins = 50) +
  theme_classic()

# LMM with lmer() --------------------------------------------------------------
## Fit model -------------------------------------------------------------------
# with all environmental factors and their interaction with time as predictors
mod_lmer_env <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    scale(year) * scale(PBR) +
    (1 | dataset / plotID) + (1 | species),
  data = data_fil_biomes,
  na.action = "na.exclude"
)

## Visualise fixed effects -----------------------------------------------------
out <- summary(mod_lmer_env)

df_coef <- round(out$coefficients[, c(1, 2, 5)], 4) |>
  as.data.frame() |>
  rename(
    pval = `Pr(>|t|)`,
    std = `Std. Error`,
    est = Estimate
  ) |>
  rownames_to_column(var = "var") |>
  mutate(
    pvalue = ifelse(pval > 0.1, "", pval),
    pvalue = ifelse(pval < 0.05, "*", pvalue),
    pvalue = ifelse(pval < 0.01, "**", pvalue),
    pvalue = ifelse(pval < 0.001, "***", pvalue)
  )

df_coef_plot <- df_coef |>
  mutate(
    eff = ifelse(row_number() %in% 1:7, "Main effect", "Interaction terms"),
    eff = as_factor(eff)
  ) |>
  mutate(
    varnew = ifelse(var == "scale(year)", "year", var),
    varnew = ifelse(var == "scale(ai)", "MI", varnew),
    varnew = ifelse(var == "scale(ndep)", "Ndep", varnew),
    varnew = ifelse(var == "scale(ORGC)", "ORGC", varnew),
    varnew = ifelse(var == "scale(PBR)", "PBR", varnew),
    varnew = ifelse(var == "scale(year):scale(ai)", "MI", varnew),
    varnew = ifelse(var == "scale(year):scale(ndep)", "Ndep", varnew),
    varnew = ifelse(var == "scale(year):scale(ORGC)", "ORGC", varnew),
    varnew = ifelse(var == "scale(year):scale(PBR)", "PBR", varnew),
    varnew = as_factor(varnew),
    varnew = fct_relevel(varnew, c("PBR", "Ndep", "ORGC", "MI"))
  ) |>
  filter(varnew == "MI" | varnew == "Ndep" | varnew == "PBR" | varnew == "ORGC")

## Save model object and coefficients table ------------------------------------
saveRDS(mod_lmer_env, file = here::here("data/mod_lmer_env.rds"))
saveRDS(df_coef_plot, file = here::here("data/df_coef_plot.rds"))

## Plot effects ----------------------------------------------------------------
fig2a <- df_coef_plot |>
  slice(1:4) |> # only main effects
  ggplot() +
  geom_point(
    aes(
      x = varnew,
      y = est
    ),
    size = 3
  ) +
  geom_errorbar(
    aes(
      x = varnew,
      ymin = est - 1.96 * std,
      ymax = est + 1.96 * std
    ),
    width = 0
  ) +
  labs(x = "", y = "Coefficient (scaled)", title = "Main effects") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_classic() +
  scale_x_discrete(
    labels = c(
      "MI" = "Moisture\n Index",
      "ORGC" = "Organic\n carbon",
      "Ndep" = "Nitrogen\n deposition",
      "PBR" = "Phosporus\n availability"
    )
  ) +
  coord_flip()

fig2b <- df_coef_plot |>
  slice(5:8) |> # only main effects
  ggplot() +
  geom_point(
    aes(
      x = varnew,
      y = est
    ),
    size = 3
  ) +
  geom_errorbar(
    aes(
      x = varnew,
      ymin = est - 1.96 * std,
      ymax = est + 1.96 * std
    ),
    width = 0
  ) +
  labs(x = "", y = "Coefficient (scaled)", title = "Interactions with year") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_classic() +
  scale_x_discrete(
    labels = NULL
  ) +
  coord_flip()

cowplot::plot_grid(
  fig2a,
  fig2b,
  labels = letters[1:2],
  rel_widths = c(1, 0.83)
)

ggsave(
  here::here("manuscript/figures/fig2.pdf"),
  width = 6,
  height = 3
)

# # out-of-the-box method
# plot_model(
#   mod_lmer_env,
#   type = "est",
#   terms = c("scale(ai)", "scale(ndep)", "scale(ORGC)", "scale(PBR)"),
#   title = "Scaled Predictors",
#   se = TRUE,
#   show.values = TRUE,
#   axis.lim = c(-0.01, 0.01)
# )

## Various diagnostics ---------------------------------------------------------
summary(mod_lmer_env)
r.squaredGLMM(mod_lmer_env)
AIC(mod_lmer_env)

plot(allEffects(mod_lmer_env))

ggplot() +
  geom_point(
    data = data_fil_biomes |>
      filter(country == "Switzerland"),
    aes(x = year, y = ndep, color = dataset),
    alpha = 0.5,
    size = 1.5,
    inherit.aes = FALSE
  )

plot_model(mod_lmer_env, type = "pred", show.data = TRUE, dot.size = 1.0, terms = c("logQMD", "ai"))
plot_model(mod_lmer_env, type = "pred", show.data = TRUE, dot.size = 1.0, terms = c("logQMD", "ndep"))

plot_model(mod_lmer_env, type = "pred", show.data = TRUE, dot.size = 1.0, terms = c("year", "ai"))
plot_model(mod_lmer_env, type = "pred", show.data = TRUE, dot.size = 1.0, terms = c("year", "ndep"))

plot_model(mod_lmer_env)

# LMM with lqmm() --------------------------------------------------------------

data_unm <- readRDS(here::here("data/data_unm.rds"))

data_unm <- data_unm |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD),
    ai_sc = scale(ai),
    ndep_sc = scale(ndep),
    orgc_sc = scale(ORGC),
    pbr_sc = scale(PBR)
  )

### Identify disturbed plots ---------------------------------------------------
data_unm <- data_unm |>
  identify_disturbed_plots()

breaks <- get_breaks(data_unm$year)

df_disturbed <- data_unm |>
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

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed_biome12 <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",
    title = bquote(bold("d") ~ ~"Mediterranean Forests")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

# Additional filter: remove plots with no change in ln(N)
data_unm <- data_unm |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# remove disturbed plots
data_unm_including_disturbed <- data_unm
data_unm <- data_unm |>
  filter(ndisturbed == 0)

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm,
  type = "normal",
  control = lqmmControl(startQR = TRUE)
)

# with all environmental factors and their interaction with time as predictors
set.seed(123)
mod_lqmm_env <- lqmm(
  logDensity ~ logQMD_sc +
    year_sc * ai_sc +
    year_sc * ndep_sc +
    year_sc * orgc_sc +
    year_sc * pbr_sc,
  random = ~1,
  group = plotID,
  tau = 0.70,
  data = data_unm |>
    drop_na(),
  type = "normal",
  control = list(
    # LP_max_iter = 2000,
    # LP_tol_ll = 1e-3,
    startQR = TRUE
  )
)

set.seed(123)
mod_lqmm_env <- lqmm(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    scale(year) * scale(PBR),
  random = ~1,
  group = plotID,
  tau = 0.70,
  data = data_unm |>
    drop_na(),
  type = "normal",
  control = list(
    # LP_max_iter = 2000,
    # LP_tol_ll = 1e-3,
    startQR = TRUE
  )
)

out <- summary(mod_lqmm_env)

plot_lqmm_bybiome(
  data_unm,
  mod_lqmm_env
)

## Fit model -------------------------------------------------------------------

# # with all environmental factors and their interaction with time as predictors
mod_lqmm_env <- lqmm(
  logDensity ~ logQMD +
    year * ai +
    year * ndep +
    year * ORGC +
    year * PBR,
  random = ~1,
  group = plotID,
  tau = 0.80,
  data = data_unm_biome |>
    select(logQMD, logDensity, year, ai, ndep, ORGC, PBR, plotID) |>
    drop_na(),
  type = "normal",
  control = list(
    LP_max_iter = 2000,
    LP_tol_ll = 1e-4
  )
)
#

## Visualise fixed effects -----------------------------------------------------
# out <- coef(mod_lqmm_env) |>
#   as.data.frame() |>
#   rownames_to_column(var = "var") |>
#   select(var, value = `0.9`)
#
# # this takes too long
# out <- summary(mod_lqmm_env)
#
#
# pval(mod_lqmm_env)
#
# df_coef <- round(out$coefficients[,c(1,2,5)],4) |>
#   as.data.frame() |>
#   rename(
#     pval = `Pr(>|t|)`,
#     std = `Std. Error`,
#     est = Estimate
#   ) |>
#   rownames_to_column(var = "var") |>
#   mutate(pvalue = ifelse(pval>0.1,"", pval),
#          pvalue = ifelse(pval<0.05,"*", pvalue),
#          pvalue = ifelse(pval<0.01,"**", pvalue),
#          pvalue = ifelse(pval<0.001,"***",pvalue))
#
# df_coef_plot <- df_coef |>
#   mutate(eff = ifelse(row_number() %in% 1:7, "Main effect", "Interaction terms"),
#          eff = as_factor(eff)) |>
#   mutate(varnew = ifelse(var == "scale(year)", "year", var),
#          varnew = ifelse(var == "scale(ai)", "MI", varnew),
#          varnew = ifelse(var == "scale(ndep)", "Ndep", varnew),
#          varnew = ifelse(var == "scale(ORGC)", "ORGC", varnew),
#          varnew = ifelse(var == "scale(PBR)", "PBR", varnew),
#          varnew = ifelse(var == "scale(year):scale(ai)", "MI", varnew),
#          varnew = ifelse(var == "scale(year):scale(ndep)", "Ndep", varnew),
#          varnew = ifelse(var == "scale(year):scale(ORGC)", "ORGC", varnew),
#          varnew = ifelse(var == "scale(year):scale(PBR)", "PBR", varnew),
#          varnew = as_factor(varnew),
#          varnew = fct_relevel(varnew,c("PBR","Ndep","ORGC","MI"))) |>
#   filter(varnew == "MI"|varnew == "Ndep"|varnew ==  "PBR"|varnew == "ORGC")
#
# fig2a <- df_coef_plot |>
#   slice(1:4) |>  # only main effects
#   ggplot() +
#   geom_point(
#     aes(
#       x = varnew,
#       y = est
#     ),
#     size = 3
#   ) +
#   geom_errorbar(
#     aes(
#       x = varnew,
#       ymin = est - 1.96 * std,
#       ymax = est + 1.96 * std
#     ),
#     width = 0
#   ) +
#   labs(x = "", y = "Coefficient (scaled)", title = "Main effects") +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic()  +
#   scale_x_discrete(
#     labels = c(
#       "MI" = "Moisture\n Index",
#       "ORGC" = "Organic\n carbon",
#       "Ndep" = "Nitrogen\n deposition",
#       "PBR" = "Phosporus\n availability")
#   ) +
#   coord_flip()
#
# fig2b <- df_coef_plot |>
#   slice(5:8) |>  # only main effects
#   ggplot() +
#   geom_point(
#     aes(
#       x = varnew,
#       y = est
#     ),
#     size = 3
#   ) +
#   geom_errorbar(
#     aes(
#       x = varnew,
#       ymin = est - 1.96 * std,
#       ymax = est + 1.96 * std
#     ),
#     width = 0
#   ) +
#   labs(x = "", y = "Coefficient (scaled)", title = "Interactions with year") +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic()  +
#   scale_x_discrete(
#     labels = NULL
#   ) +
#   coord_flip()
#
# cowplot::plot_grid(
#   fig2a,
#   fig2b,
#   labels = letters[1:2],
#   rel_widths = c(1, 0.83)
# )
#
# ggsave(
#   here::here("manuscript/figures/fig2.pdf"),
#   width = 6,
#   height = 3
# )
