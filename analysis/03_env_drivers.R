# Fixed effects plot -----------------------------------------------------------
# This script analyses the main environmental drivers affecting the STL changes

## Load packages ---------------------------------------------------------------
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
library(here)
library(broom)
library(kableExtra)
library(modelsummary)

# load data
data_fil_biomes <- read_rds(here("data/inputs/data_fil75_biomes.rds"))

# LMM with lmer() --------------------------------------------------------------
## Fit model -------------------------------------------------------------------
# with all environmental factors and their interaction with time as predictors
mod_lmer_env_complete <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    scale(year) * scale(PBR) +
    scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = data_fil_biomes,
  na.action = "na.exclude"
)

write_rds(mod_lmer_env_complete, file = here("data/mod_lmer_env_complete.rds"))

# fit model without PBR - this turns out as the best one
mod_lmer_env_nopbr <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    # scale(year) * scale(PBR) + # excluded because neither fixed nor interactive effect is significant
    scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = data_fil_biomes,
  na.action = "na.exclude"
)

write_rds(mod_lmer_env_nopbr, file = here("data/mod_lmer_env_nopbr.rds"))

# fit model without PBR and ORGC
mod_lmer_env_nopbr_noorgc <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    # scale(year) * scale(ORGC) +
    # scale(year) * scale(PBR) + # excluded because neither fixed nor interactive effect is significant
    scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = data_fil_biomes,
  na.action = "na.exclude"
)

write_rds(mod_lmer_env_nopbr_noorgc, file = here("data/mod_lmer_env_nopbr_noorgc.rds"))

# fit model without PBR and CNrt
mod_lmer_env_nopbr_nocn <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    # scale(year) * scale(PBR) + # excluded because neither fixed nor interactive effect is significant
    # scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = data_fil_biomes,
  na.action = "na.exclude"
)

write_rds(mod_lmer_env_nopbr_nocn, file = here("data/mod_lmer_env_nopbr_nocn.rds"))


## Write summary of model fits to latex table
options("modelsummary_format_numeric_latex" = "plain")

modelsummary(
  list(
    "Complete" = mod_lmer_env_complete,
    "No PBR" = mod_lmer_env_nopbr,
    "No PBR, ORGC" = mod_lmer_env_nopbr_noorgc,
    "No PBR, C:N" = mod_lmer_env_nopbr_nocn
  ),
  output = here("manuscript/tables/mods_env.tex"),
  # estimate  = "p.value",
  estimate = "{estimate}{stars}",
  # estimate  = "{estimate} [{conf.low}, {conf.high}], {p.value}",
  title = "Regression Results",
  statistic = "std.error",
  coef_omit = "Intercept"
)

# select best model for further analyses looking at AIC and BIC in latex table
mod_lmer_env <- mod_lmer_env_nopbr

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
    eff = ifelse(row_number() %in% 1:8, "Main effect", "Interaction terms"),
    eff = as_factor(eff)
  ) |>
  mutate(
    varnew = ifelse(var == "scale(year)", "year", var),
    varnew = ifelse(var == "scale(tavg)", "MAT", varnew),
    varnew = ifelse(var == "scale(ai)", "MI", varnew),
    varnew = ifelse(var == "scale(ndep)", "Ndep", varnew),
    varnew = ifelse(var == "scale(ORGC)", "ORGC", varnew),
    varnew = ifelse(var == "scale(PBR)", "PBR", varnew),
    varnew = ifelse(var == "scale(CNrt)", "C:N", varnew),
    varnew = ifelse(var == "scale(year):scale(tavg)", "MAT", varnew),
    varnew = ifelse(var == "scale(year):scale(ai)", "MI", varnew),
    varnew = ifelse(var == "scale(year):scale(ndep)", "Ndep", varnew),
    varnew = ifelse(var == "scale(year):scale(ORGC)", "ORGC", varnew),
    varnew = ifelse(var == "scale(year):scale(PBR)", "PBR", varnew),
    varnew = ifelse(var == "scale(year):scale(CNrt)", "C:N", varnew),
    varnew = as_factor(varnew),
    varnew = fct_relevel(varnew, c("C:N", "Ndep", "ORGC", "MI", "MAT"))
  ) |>
  filter(varnew == "MAT" | varnew == "MI" | varnew == "Ndep" |
    varnew == "PBR" | varnew == "ORGC" | varnew == "C:N")

## Save model object and coefficients table ------------------------------------
saveRDS(mod_lmer_env, file = here("data/mod_lmer_env.rds"))
saveRDS(df_coef_plot, file = here("data/df_coef_plot.rds"))

# # Figure 2 (OLD VERSION)
# fig2 <- ggplot(df_coef_plot) +
#   geom_bar(aes(y = varnew, weight = est, fill = eff), position = position_stack(reverse = TRUE), width = .5) +
#   geom_text(aes(y = varnew, x = est, label = pvalue),
#     color = "white", size = 5,
#     position = position_stack(vjust = 0.5), vjust = 0.75
#   ) +
#   xlab("Coefficients") +
#   ylab("Environmental drivers") +
#   labs(fill = "Predictors") +
#   geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = .8) +
#   theme_classic() +
#   theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     axis.text = element_text(size = 10), axis.title = element_text(size = 10),
#     axis.text.y = element_text(hjust = 0.5),
#     legend.text = element_text(size = 9.5), legend.title = element_text(size = 9.5),
#     plot.title = element_text(size = 10),
#     legend.key = element_rect(fill = NA, color = NA),
#     legend.position = c(0.85, 0.15),
#     legend.direction = "vertical",
#     legend.box = "horizontal",
#     legend.margin = margin(2, 2, 2, 2),
#     legend.key.size = unit(.6, "cm"),
#     legend.box.margin = margin(1, 1, 1, 1)
#   ) +
#   scale_x_continuous(limits = c(-0.04, 0.055), breaks = seq(-0.05, 0.05, 0.025)) +
#   scale_y_discrete(labels = c(
#     "MI" = "Moisture \nIndex",
#     "Tavg" = "Mean \nTemperature",
#     "ORGC" = "Organic \ncarbon",
#     "Ndep" = "Nitrogen \ndeposition",
#     "PBR" = "Phosporus \navailability"
#   )) +
#   scale_fill_okabe_ito()
# fig2
# ggsave(paste0(here(), "/manuscript/figures/fig2_color.png"), width = 8, height = 5, dpi = 300)

## Plot effects ----------------------------------------------------------------
fig2a <- df_coef_plot |>
  slice(1:5) |> # only main effects
  ggplot() +
  geom_point(
    aes(
      x = varnew,
      y = est
    ),
    size = 2
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
      "MAT" = "Mean\n Temperature",
      "MI" = "Moisture\n Index",
      "ORGC" = "Organic\n carbon",
      "Ndep" = "Nitrogen\n deposition",
      "PBR" = "Phosporus\n availability"
    )
  ) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
  ) +
  coord_flip()

fig2b <- df_coef_plot |>
  slice(6:10) |> # only main effects
  ggplot() +
  geom_point(
    aes(
      x = varnew,
      y = est
    ),
    size = 2
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
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
  ) +
  coord_flip()

cowplot::plot_grid(
  fig2a,
  fig2b,
  labels = letters[1:2],
  rel_widths = c(1, 0.83)
)

ggsave(
  here("manuscript/figures/fig2_tavg.png"),
  width = 6,
  height = 3,
  dpi = 300
  )

ggsave(
  here("manuscript/figures/fig2_tavg.pdf"),
  width = 6,
  height = 3
)

# LMM with lqmm() --------------------------------------------------------------

# data_unm <- readRDS(here("data/inputs/data_unm.rds")) |>
#   tidyr::drop_na(
#     year, logQMD, tavg, ai, ndep, ORGC, PBR, CNrt, lai
#   ) |>
#   mutate(
#     year_sc = scale(year),
#     logQMD_sc = scale(logQMD),
#     tavg_sc = scale(tavg),
#     ai_sc = scale(ai),
#     ndep_sc = scale(ndep),
#     ORGC_sc = scale(ORGC),
#     PBR_sc = scale(PBR),
#     CNrt_sc = scale(CNrt),
#     lai_sc = scale(lai)
#   )
#
# visdat::vis_miss(data_unm, warn_large_data = FALSE)
#
# mod_lqmm_env <- lqmm(
#   logDensity ~ logQMD_sc +
#     year_sc * tavg_sc +
#     year_sc * ai_sc +
#     year_sc * ndep_sc +
#     year_sc * ORGC_sc +
#     year_sc * PBR_sc +
#     year_sc * CNrt_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.9,
#   data = data_unm,
#   type = "normal",
#   control = lqmmControl(
#     LP_max_iter = 5000, # inner loop iterations
#     LP_tol_ll   = 1e-05, # inner loop tolerance
#     startQR     = TRUE
#   )
# )


data_unm <- readRDS(here("data/inputs/data_unm.rds"))

data_unm <- data_unm |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD),
    tavg_sc = scale(tavg),
    ai_sc = scale(ai),
    ndep_sc = scale(ndep),
    orgc_sc = scale(ORGC),
    pbr_sc = scale(PBR)
  )

### Identify disturbed plots ---------------------------------------------------
source(here("R/identify_disturbed_plots.R"))
source(here("R/get_breaks.R"))

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
  filter(nplots > 30) |>
  mutate(fdisturbed = ndisturbed / nplots) |>
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed <- df_disturbed |>
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year"
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

gg_fdisturbed

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
    year_sc * tavg_sc +
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
    LP_max_iter = 5000,
    LP_tol_ll   = 1e-05,
    startQR     = TRUE
  )
)

summary(mod_lqmm_env)

out <- summary(mod_lqmm_env)

plot_lqmm_bybiome(
  data_unm,
  mod_lqmm_env
)
